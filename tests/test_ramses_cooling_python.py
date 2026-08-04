#!/usr/bin/env python3
"""
Verify the RAMSESCooling Python class against a RAMSES-exported cooling table
read via the ozymandias Fortran wrapper.

The comparison is done on the full table grid (n1=161 nH bins × n2=101 T/mu bins):
  - primordial H+He cooling coefficient vs mytable.cool
  - photo-heating coefficient         vs mytable.heat
  - metal cooling coefficient (Z=1)   vs mytable.metal
  - net cooling (without Compton)     vs cool + Z*metal - heat (from table)

Compton cooling/heating (mytable.cool_com, heat_com) is NOT implemented in
RAMSESCooling and is excluded from the comparison; instead it is shown
separately on the diagnostic plots.

Usage:
    python test_ramses_cooling_python.py /path/to/cooling_NNNNN.out [z_redshift]

If z_redshift is omitted the script reads aexp from info_NNNNN.txt in the same
directory and computes z = 1/aexp - 1.

Two PNG files are saved next to the cooling file:
    cooling_NNNNN_python_components.png   — component-by-component comparison
    cooling_NNNNN_python_net.png          — net cooling maps + equilibrium curves
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# RAMSESCooling lives in ozy/cooling.py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ozy.cooling import RAMSESCooling  # noqa: E402 (path setup must precede)


# ─────────────────────────────────────────────────────────────────────────────
# helpers
# ─────────────────────────────────────────────────────────────────────────────
def _infer_z(cool_file):
    """Parse aexp from info_NNNNN.txt in the same directory; return z = 1/aexp - 1."""
    dirname = os.path.dirname(os.path.abspath(cool_file))
    for fname in sorted(os.listdir(dirname)):
        if fname.startswith('info_') and fname.endswith('.txt'):
            with open(os.path.join(dirname, fname)) as f:
                for line in f:
                    if 'aexp' in line and '=' in line:
                        aexp = float(line.split('=')[1].strip())
                        return 1.0 / aexp - 1.0
    raise ValueError(f"No info_*.txt with aexp found in {dirname!r}")


def _rel_err(py_val, ref_val, floor=1e-300):
    return (py_val - ref_val) / np.maximum(np.abs(ref_val), floor)


def _print_stats(label, py_val, ref_val, mask=None):
    if mask is None:
        mask = ref_val > 0
    err = np.abs(_rel_err(py_val[mask], ref_val[mask]))
    print(f"  {label:<30s}  median |err| = {np.median(err):.3e}   "
          f"max |err| = {np.max(err):.3e}   "
          f"frac > 10% = {(err > 0.1).mean():.3f}")


# ─────────────────────────────────────────────────────────────────────────────
# main comparison
# ─────────────────────────────────────────────────────────────────────────────
def run_comparison(cool_file, z):
    import sys, os as _os
    _amr_dir = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                             'ozy', 'amr')
    if _amr_dir not in sys.path:
        sys.path.insert(0, _amr_dir)
    from amr2_pkg import cooling_module

    # disable Fortran self-shielding so comparison uses boost = 1 on both sides
    try:
        cooling_module.self_shielding = False
        print("Fortran self_shielding set to False")
    except AttributeError:
        print("[warn] cooling_module.self_shielding not settable; "
              "default may introduce self-shielding differences")

    mytable = cooling_module.cooling_table()
    cooling_module.retrieve_table(cool_file, mytable)

    n1, n2   = mytable.n1, mytable.n2
    nH_ax    = 10.**np.asarray(mytable.nh)   # (n1,)  — nH [cm^-3]
    T2_ax    = 10.**np.asarray(mytable.t2)   # (n2,)  — T/mu [K]
    mu_2d    = np.asarray(mytable.mu)        # (n1, n2)

    print(f"\nTable: n1={n1} (nH), n2={n2} (T/mu)  |  z = {z:.4f}")
    print(f"  nH   : {nH_ax[0]:.2e} – {nH_ax[-1]:.2e} cm^-3")
    print(f"  T/mu : {T2_ax[0]:.2e} – {T2_ax[-1]:.2e} K")

    # reference arrays (all stored as log10; convert to linear)
    cool_ref  = 10.**np.asarray(mytable.cool)       # (n1, n2)  [erg s^-1 cm^3]
    heat_ref  = 10.**np.asarray(mytable.heat)
    metal_ref = 10.**np.asarray(mytable.metal)      # at Z_sun = 1
    ccomp_ref = 10.**np.asarray(mytable.cool_com)   # [erg s^-1 cm^3] × nH (Compton)
    hcomp_ref = 10.**np.asarray(mytable.heat_com)

    # 2D grid — T_phys uses mu from the table (T/mu is the table axis variable)
    nH_2d     = nH_ax[:, None] * np.ones((n1, n2))   # (n1, n2)
    T_phys_2d = T2_ax[None, :] * mu_2d               # (n1, n2)  physical temperature

    # RAMSESCooling.get_cooling_heating expects 1D arrays (its boost broadcast is
    # `shape (3,1) × shape (N,)`) — flatten, call, then reshape back to (n1, n2)
    nH_flat    = nH_2d.ravel()
    T_flat     = T_phys_2d.ravel()
    boost_flat = np.ones(n1 * n2)

    def _call_py(Z):
        Lr, Hv = py.get_cooling_heating(T_flat, nH_flat, z, Z, boost_flat)
        return Lr.reshape(n1, n2), Hv.reshape(n1, n2)

    py = RAMSESCooling()

    # ── primordial cooling and photoheating (Z=0 → metal=0) ────────────────
    L_rad_Z0, H_uvb = _call_py(0.0)
    prim_py = L_rad_Z0 / nH_2d**2   # (n1, n2) [erg s^-1 cm^3]
    heat_py = H_uvb    / nH_2d**2

    # ── metal cooling coefficient at Z=1 ────────────────────────────────────
    L_rad_Z1, _ = _call_py(1.0)
    metal_py = L_rad_Z1 / nH_2d**2 - prim_py   # (n1, n2) [erg s^-1 cm^3] at Z=1

    print("\nComponent statistics (absolute relative error):")
    _print_stats("Primordial cool",  prim_py,  cool_ref)
    _print_stats("Photoheating",     heat_py,  heat_ref)
    _print_stats("Metal cool (Z=1)", metal_py, metal_ref)

    # ── verify f_courty implementation: ratio should be ~1 after the fix ──────
    # metal_py now includes (cc07 + fine-structure) × f_courty, so the ratio
    # metal_ref / metal_py should be close to 1 across the grid.
    safe_metal_py = np.maximum(metal_py, 1e-300)
    f_courty_implied = metal_ref / safe_metal_py   # (n1, n2) — should be ~1 now
    f_courty_vs_nH = np.median(f_courty_implied, axis=1)  # (n1,)
    print("\nResidual metal ratio metal_tbl / metal_py  (median over T axis; should be ~1):")
    for i in [0, 10, 30, 60, 90, 120, 140, 160]:
        print(f"  nH = {nH_ax[i]:.2e} cm^-3  →  ratio ≈ {f_courty_vs_nH[i]:.4e}")

    # ── net cooling comparison at three metallicities ───────────────────────
    Z_vals  = [0.001, 0.1, 1.0]
    net_ref = {}   # table: cool + Z*metal - heat  (no Compton)
    net_py  = {}

    print("\nNet cooling statistics (no Compton):")
    for Z in Z_vals:
        net_ref[Z] = cool_ref + Z * metal_ref - heat_ref
        L_rad_Z, _ = _call_py(Z)
        net_py[Z]  = L_rad_Z / nH_2d**2 - heat_py
        _print_stats(f"Net cool Z={Z:.3f}", net_py[Z], net_ref[Z],
                     mask=(np.abs(net_ref[Z]) > 1e-60))

    # ── Figure 1: component comparisons ─────────────────────────────────────
    fig1, axes = plt.subplots(3, 3, figsize=(16, 12))
    fig1.suptitle(f'RAMSESCooling component verification   z = {z:.3f}   '
                  f'(no self-shielding)', fontsize=13)

    rows = [
        ('Primordial H+He cool',  prim_py,  cool_ref,
         r'$\Lambda_{\rm prim}$ [erg s$^{-1}$ cm$^3$]'),
        ('Photo-heating',          heat_py,  heat_ref,
         r'$\Gamma_{\rm UV}$ [erg s$^{-1}$ cm$^3$]'),
        ('Metal cool  $Z=1\,Z_\odot$', metal_py, metal_ref,
         r'$\Lambda_{\rm metal}$ [erg s$^{-1}$ cm$^3$]'),
    ]

    for r, (title, py_v, ref_v, clabel) in enumerate(rows):
        ax_ref, ax_py, ax_err = axes[r]

        pos_mask = ref_v > 0
        vmin = np.percentile(ref_v[pos_mask], 2)
        vmax = ref_v[pos_mask].max()
        kw   = dict(shading='auto', cmap='viridis')

        im0 = ax_ref.pcolormesh(nH_ax, T2_ax, ref_v.T,
                                norm=LogNorm(vmin=vmin, vmax=vmax), **kw)
        im1 = ax_py.pcolormesh(nH_ax, T2_ax,
                               np.maximum(py_v, vmin * 1e-3).T,
                               norm=LogNorm(vmin=vmin, vmax=vmax), **kw)

        err2d = _rel_err(py_v, ref_v)
        lim   = max(np.percentile(np.abs(err2d), 98), 1e-4)
        im2   = ax_err.pcolormesh(nH_ax, T2_ax, err2d.T,
                                  vmin=-lim, vmax=lim,
                                  shading='auto', cmap='RdBu_r')

        for ax, im, cb_lbl, ttl in [
            (ax_ref, im0, clabel,          f'{title}  —  Table'),
            (ax_py,  im1, clabel,          f'{title}  —  Python'),
            (ax_err, im2, 'Relative error', f'{title}  —  Relative error'),
        ]:
            ax.set_xscale('log'); ax.set_yscale('log')
            ax.set_xlabel(r'$n_H$ [cm$^{-3}$]', fontsize=9)
            ax.set_ylabel(r'$T/\mu$ [K]', fontsize=9)
            ax.set_title(ttl, fontsize=9)
            plt.colorbar(im, ax=ax, label=cb_lbl, fraction=0.046, pad=0.04)

    fig1.tight_layout()
    out1 = cool_file.replace('.out', '_python_components.png')
    fig1.savefig(out1, dpi=150, bbox_inches='tight')
    print(f"\nSaved {out1}")

    # ── Figure 2: net cooling maps + Compton + equilibrium curves ───────────
    fig2, axes2 = plt.subplots(2, 3, figsize=(16, 9))
    fig2.suptitle(f'Net cooling: Table vs Python   z = {z:.3f}   '
                  f'(no Compton in Python)', fontsize=13)

    compton_net = (ccomp_ref - hcomp_ref) / nH_2d   # [erg s^-1 cm^3]

    for col, Z in enumerate(Z_vals):
        ax_top = axes2[0, col]
        ax_bot = axes2[1, col]

        ref_net = net_ref[Z]                          # table without Compton
        ref_full = ref_net + compton_net              # table including Compton

        # top panel: table net cooling + equilibrium curves
        absf = np.abs(ref_full)
        im = ax_top.pcolormesh(nH_ax, T2_ax, absf.T,
                               norm=LogNorm(vmin=1e-40, vmax=1e-20),
                               shading='auto', cmap='plasma')
        plt.colorbar(im, ax=ax_top,
                     label=r'$|\Lambda_{\rm net}|$ [erg s$^{-1}$ cm$^3$]',
                     fraction=0.046, pad=0.04)

        # equilibrium curves (where sign of net cooling changes)
        T_eq_ref = T2_ax[np.argmin(np.abs(ref_full), axis=1)]
        T_eq_py  = T2_ax[np.argmin(np.abs(net_py[Z]), axis=1)]
        ax_top.plot(nH_ax, T_eq_ref, 'w-',  lw=2.0, label='Table eq. (with Compton)')
        ax_top.plot(nH_ax, T_eq_py,  'r--', lw=2.0, label='Python eq. (no Compton)')
        ax_top.set_xscale('log'); ax_top.set_yscale('log')
        ax_top.set_xlabel(r'$n_H$ [cm$^{-3}$]', fontsize=9)
        ax_top.set_ylabel(r'$T/\mu$ [K]', fontsize=9)
        ax_top.set_title(f'Table $|\\Lambda_{{\\rm net}}|$   $Z = {Z}\\,Z_\\odot$', fontsize=9)
        if col == 0:
            ax_top.legend(fontsize=7)

        # bottom panel: relative error of Python vs Table (both without Compton)
        err = _rel_err(net_py[Z], ref_net, floor=np.abs(ref_net).max() * 1e-12)
        lim = max(np.percentile(np.abs(err), 98), 0.01)
        im2 = ax_bot.pcolormesh(nH_ax, T2_ax, err.T,
                                vmin=-lim, vmax=lim,
                                shading='auto', cmap='RdBu_r')
        plt.colorbar(im2, ax=ax_bot, label='(Python - Table) / |Table|',
                     fraction=0.046, pad=0.04)
        ax_bot.set_xscale('log'); ax_bot.set_yscale('log')
        ax_bot.set_xlabel(r'$n_H$ [cm$^{-3}$]', fontsize=9)
        ax_bot.set_ylabel(r'$T/\mu$ [K]', fontsize=9)
        ax_bot.set_title(f'Relative error  $Z = {Z}\\,Z_\\odot$', fontsize=9)

    fig2.tight_layout()
    out2 = cool_file.replace('.out', '_python_net.png')
    fig2.savefig(out2, dpi=150, bbox_inches='tight')
    print(f"Saved {out2}")

    return {
        'prim_py': prim_py,   'cool_ref': cool_ref,
        'heat_py': heat_py,   'heat_ref': heat_ref,
        'metal_py': metal_py, 'metal_ref': metal_ref,
        'net_py': net_py,     'net_ref': net_ref,
        'nH_ax': nH_ax,       'T2_ax': T2_ax,
        'compton_net': compton_net,
    }


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    cool_file = sys.argv[1]
    z = float(sys.argv[2]) if len(sys.argv) > 2 else _infer_z(cool_file)
    print(f"Cooling file : {cool_file}")
    print(f"Redshift z   : {z:.4f}")
    run_comparison(cool_file, z)

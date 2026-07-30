# ozy/__init__.py — Cleaned version

import importlib
import sys
import warnings
import os
import types

# --- compiled extension modules (full package paths) ---
_EXT_MODULES = {
    "_amr2_pkg": "ozy.amr._amr2_pkg",
    "_part2_pkg": "ozy.part._part2_pkg",
    "_vis_pkg": "ozy.visualisation._vis_pkg",
    # Additional compiled wrapper modules used elsewhere in the package
    "_tutils_pkg": "ozy.TreeMaker._tutils_pkg",
    "_hutils_pkg": "ozy.HaloFinder._hutils_pkg",
    # Some wrappers are built without the '_pkg' suffix (f2py/f90wrap variations)
    "_tutils": "ozy.TreeMaker._tutils",
    "_hutils": "ozy.HaloFinder._hutils",
}

# Import compiled extensions and alias top-level names to prevent duplicates
for topname, fullpath in _EXT_MODULES.items():
    try:
        mod = importlib.import_module(fullpath)
        sys.modules.setdefault(topname, mod)
        base_name = topname.replace("_pkg", "")
        if base_name not in sys.modules:
            sys.modules[base_name] = mod
    except Exception:
        warnings.warn(
            f"Could not import compiled extension '{topname}' from '{fullpath}'. "
            "Some wrapper functions may be unavailable."
        )

# If the generated wrapper package uses absolute imports like
# `import amr2_pkg.filtering` (common with f90wrap/f2py output),
# create a lightweight top-level package placeholder with the correct
# __path__ before importing the package-relative wrapper module. This
# lets absolute imports inside the wrapper resolve to the package's
# directory while keeping the canonical module name `ozy.<subpkg>.<pkg>`.
base_dir = os.path.dirname(__file__)
# Map wrapper package names to their expected subdirectory under `ozy`.
_WRAPPER_SUBDIR = {
    'amr2_pkg': 'amr',
    'part2_pkg': 'part',
    'tutils_pkg': 'TreeMaker',
    'hutils_pkg': 'HaloFinder',
    'vis_pkg': 'visualisation',
}

for wrapper_name in _WRAPPER_SUBDIR:
    candidate = os.path.join(base_dir, _WRAPPER_SUBDIR[wrapper_name], wrapper_name)
    if os.path.isdir(candidate) and wrapper_name not in sys.modules:
        placeholder = types.ModuleType(wrapper_name)
        placeholder.__path__ = [candidate]
        sys.modules[wrapper_name] = placeholder

# Import the wrapper packages under the top-level names (e.g. 'amr2_pkg')
# and then alias them to their package-relative names to avoid executing
# the same registration code twice under different module identities.
for wrapper_name, pkg_rel in (
    ('amr2_pkg', 'ozy.amr.amr2_pkg'),
    ('part2_pkg', 'ozy.part.part2_pkg'),
    ('tutils_pkg', 'ozy.TreeMaker.tutils_pkg'),
    ('hutils_pkg', 'ozy.HaloFinder.hutils_pkg'),
    ('vis_pkg', 'ozy.visualisation.vis_pkg'),
):
    try:
        # import as top-level package name so internal absolute imports
        # like `import amr2_pkg.filtering` resolve to the same module object
        top_mod = importlib.import_module(wrapper_name)
        sys.modules.setdefault(wrapper_name, top_mod)
        # alias into package-relative name
        sys.modules.setdefault(pkg_rel, top_mod)
        # also create a shorter historical alias if useful
        short = wrapper_name.replace('_pkg', '')
        sys.modules.setdefault(short, top_mod)
    except Exception:
        # not fatal; wrapper package might not be present/built
        pass

# --- Public API: package-relative imports ---
from .loader import load
from .main import Snapshot, CosmoSnapshot
from .driver import drive
from . import variables_settings
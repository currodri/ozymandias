#
# Makefile for HaloMaker
#
# 1. general options 
# --- BIG_RUN option for a normal (as opposed to resimulation) simulation whatever the i/o format.
# --- RENORM option if you are analyzing the DM part of a hydro simulation and want to renormalize 
# 	the mass of your dark matter particle so as to include the baryons (by pretending they are distributed like the DM):
# 	this might be useful when one wants to compare a pure DM simulation to its DM + hydro counterpart
# --- BN98 : if set, use a redshift-dependent value for rho_threshold to identify halos in adaptaHOP, 
#            this value is proportional to the critical overdensity given by the fit of Bryan&Norman(1998)  
#            if set, the value of rho_threshold given in the input file parameters is not used
# --- SIMPLE_OUTPUTS : output simple stuff for debugging and checks.
# --- CONTAM : compute and output contamination flag for each (sub-)halo -> tells whether a halo contains a low-res particle. This is for DM zooms. 
#
# 2/ RAMSES related OPTIONS:
#    ---------------------
# --- METALS flag if RAMSES has computed metallicity from stellar feedback and you are reading version 3.0 i/o format
# --- STARS : if set, read star particles instead of DM particles
# ------ REALTIME : is set (with STARS), this will convert conformal time of birth of each particle into an age in yr. 
# ------ JOAO : output galaxy catalog in special format usable to construct SED's etc. 
# 
# 3/ GADGET related OPTION:
#    ---------------------
# --- Test_Gd flag for testing Gadget reading routines only, especially a resimulation or hydro. That's a debugging option use to use if HaloMaker cannot read a gadget snapshot.

COMPILE = "ifx"
# options for HaloFinder 
HOPTIONS = -DBIG_RUN -DBN98
ZOPTIONS = -DCONTAM # no DBIG_RUN for zooms
# options for GalFinder
GOPTIONS = -DMETALS -DSIMPLE_OUTPUTS -DSTARS -DREALTIME -DJOAO


ifeq ($(COMPILE),"ifort")
# Intel Fortran Classic (ifort)
 F90C    = ifort -qopenmp
 DEBUG   = -g -fpe0 -traceback -C -ftrapuv #-warn all -stand f90 
 OPTIM   = -O2 -fp-model precise -qopenmp #-xavx #-pg #-fast 
 FFLAGS  = ${OPTIM} -free -cpp #-convert big_endian 
else ifeq ($(COMPILE),"ifx")
# Intel Fortran LLVM-based (ifx)
  F90C    = ifx -qopenmp
  DEBUG   = -g -fpe0 -traceback -C -ftrapuv
  OPTIM   = -O2 -fp-model precise -qopenmp
  FFLAGS  = $(OPTIM) -free -cpp
endif

ifeq ($(COMPILE),"gfortran")
 F90C    = gfortran  -ffree-line-length-none -fopenmp ###-ffree-form 
 DEBUG   = -g -ffpe-trap=invalid -fbacktrace -fbounds-check -frange-check -Wall
 OPTIM   = -O2
 FFLAGS  = ${OPTIM} -cpp
endif

.SUFFIXES: .s .z .o .f90 

#
# Rules:
#

%.o: %.f90
	$(F90C) $(FFLAGS) $(HOPTIONS) -c $*.f90

%.z: %.f90
	$(F90C) $(FFLAGS) $(ZOPTIONS) -c $*.f90

%.s: %.f90
	$(F90C) $(FFLAGS) $(GOPTIONS) -c $*.f90

default: HaloFinder HaloFinder_zoom GalFinder clean

clean:	
	rm -f *.o
	rm -f *.mod

HaloFinder: halo_defs.o utils.o subbox.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o  HaloMaker.o
	${F90C} ${FFLAGS} ${HOPTIONS}  halo_defs.o utils.o subbox.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o HaloMaker.o -o HaloFinder

HaloFinder_zoom: halo_defs.z utils.z subbox.z num_rec.z input_output.z fof_mod.z compute_neiKDtree_mod.z compute_halo_props.z  HaloMaker.z
	${F90C} ${FFLAGS} ${ZOPTIONS}  halo_defs.o utils.o subbox.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o HaloMaker.o -o HaloFinder_zoom

GalFinder: conformal_time.s halo_defs.s utils.s subbox.s num_rec.s input_output.s fof_mod.s compute_neiKDtree_mod.s compute_halo_props.s  HaloMaker.s
	${F90C} ${FFLAGS} ${GOPTIONS}  conformal_time.o halo_defs.o utils.o subbox.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o HaloMaker.o -o GalFinder


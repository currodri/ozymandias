#=======================================================================
#                   define the compiler names
#=======================================================================

CC       = gcc
#CC       = ifort
F90      = gfortran
#F90      = ifx
#F90      =  /opt/intel/composer_xe_2015.3.187/bin/intel64/ifort
PYTHON   = python

#=======================================================================
#                     additional flags
#=======================================================================

ifeq ($(F90),gfortran)
	FPP      = gfortran -E
	FPP_F90FLAGS = -x f95-cpp-input -fPIC -g
	F90FLAGS = -Ofast -fPIC -cpp -ffree-line-length-256 -DIMASS -DNOIFPORT
    FCOMP    = gfortran
    LIBS     =
endif

ifeq ($(F90),ifx)

	FPP      = gfortran -E # gfortran f90wrap temp files only. not compilation
	FPP_F90FLAGS = -x f95-cpp-input -fPIC
	F90FLAGS = -O3 -cpp -fpscomp logicals -fPIC -DIMASS # use 1 and 0 for True and False
	#F90FLAGS = -g -traceback -check bounds -cpp -fpscomp logicals -DIMASS -fPIC # use 1 and 0 for True and False
    FCOMP    = intelem # for f2py
    LIBS =
endif

CFLAGS  = -fPIC #     ==> universal for ifort, gfortran, pgi

#=======================================================================
#=======================================================================

UNAME = $(shell uname)

ifeq (${UNAME}, Darwin)
  LIBTOOL = libtool -static -o
else
  LIBTOOL = ar src
endif

# ======================================================================
# PROJECT CONFIG, do not put spaced behind the variables
# ======================================================================
# Python module name
PYTHON_MODN = part2
# mapping between Fortran and C types
KIND_MAP = kind_map

#=======================================================================
#
#=======================================================================

VPATH = ../fortran_commons


#=======================================================================
#       List all source files required for the project
#=======================================================================

# names (without suffix), f90 sources
LIBSRC_SOURCES = local utils dictionary_module linalg_module coordinates_module cooling_module \
 				read_amr_module variables statistics filtering\
				cosmology_module profiles_module integrator_module \
				export_module

# file names
LIBSRC_FILES = $(addsuffix .f90,${LIBSRC_SOURCES})

# object files
LIBSRC_OBJECTS = $(addsuffix .o,${LIBSRC_SOURCES})

# only used when cleaning up
LIBSRC_FPP_FILES = $(addsuffix .fpp,${LIBSRC_SOURCES})

#=======================================================================
#       List all source files that require a Python interface
#=======================================================================

# names (without suffix), f90 sources
LIBSRC_WRAP_SOURCES = dictionary_module linalg_module coordinates_module cooling_module \
 				read_amr_module variables statistics filtering\
				cosmology_module profiles_module integrator_module \
				export_module

# file names
LIBSRC_WRAP_FILES = $(addsuffix .f90,${LIBSRC_WRAP_SOURCES})

# object files
LIBSRC_WRAP_OBJECTS = $(addsuffix .o,${LIBSRC_WRAP_SOURCES})

# fpp files
LIBSRC_WRAP_FPP_FILES = $(addsuffix .fpp,${LIBSRC_WRAP_SOURCES})

#=======================================================================
#                 Relevant suffixes
#=======================================================================

.SUFFIXES: .f90 .fpp

#=======================================================================
#
#=======================================================================

.PHONY: all clean

%.o: %.mod

all: _${PYTHON_MODN}.so _${PYTHON_MODN}_pkg.so


clean:
	-rm -f ${LIBSRC_OBJECTS} ${LIBSRC_FPP_FILES} libsrc.a _${PYTHON_MODN}*.so \
	_${PYTHON_MODN}_pkg.so *.mod *.fpp f90wrap*.f90 f90wrap*.o *.o ${PYTHON_MODN}.py
	-rm -rf ${PYTHON_MODN}_pkg/
	-rm -rf src.*/ .f2py_f2cmap .libs/ __pycache__/


.f90.o:
	${F90} ${F90FLAGS} -c $< -o $@


.c.o:
	${CC} ${CFLAGS} -c $< -o $@


.f90.fpp:
	${FPP} ${FPP_F90FLAGS} $<  -o $@


libsrc.a: ${LIBSRC_OBJECTS}
	${LIBTOOL} $@ $?


_${PYTHON_MODN}.so: libsrc.a ${LIBSRC_FPP_FILES}
	f90wrap -m ${PYTHON_MODN} ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP}
	f2py-f90wrap --fcompiler=$(FCOMP) --build-dir . -c -m _${PYTHON_MODN} -L. -lsrc f90wrap*.f90

_${PYTHON_MODN}_pkg.so: libsrc.a ${LIBSRC_FPP_FILES}
	f90wrap -m ${PYTHON_MODN}_pkg ${LIBSRC_WRAP_FPP_FILES} -k ${KIND_MAP} -P
	f2py-f90wrap --fcompiler=$(FCOMP) --build-dir . -c -m _${PYTHON_MODN}_pkg -L. -lsrc f90wrap*.f90


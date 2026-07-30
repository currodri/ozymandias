#########################################################
# This is one of my generic Makefiles			
# Author: Sergio Martin-Alvarez  			
#########################################################
COMPILE="gfortran"

INNAME=HalosExtractor.f90
OUTNAME=HalosExtractor.out
#########################################################
# option when compiling with ifort
ifeq ($(COMPILE),"gfortran")
 F90C = gfortran
 OPTIM   = -O3
 GENFLAGS  = ${OPTIM} -cpp -mcmodel=medium
endif
#########################################################
HAlosExtractor:
	${F90C} $(INNAME) ${GENFLAGS} ${LIBS} -o $(OUTNAME)	
	cp $(OUTNAME) $(BINFOLDER)$(OUTNAME)
clean:
	rm *~
	rm *.o
#########################################################

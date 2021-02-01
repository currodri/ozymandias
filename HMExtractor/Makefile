#########################################################
# This is one of my generic Makefiles			
# Author: Sergio Martin-Alvarez  			
#########################################################
COMPILE="ifort"

INNAME=HalosExtractor.f90
OUTNAME=HalosExtractor.out
#########################################################
# option when compiling with ifort
ifeq ($(COMPILE),"ifort")
 F90C = ifort
 OPTIM   = -O2
 GENFLAGS  = ${OPTIM} -fpp -mcmodel=medium
endif
#########################################################
HAlosExtractor:
	${F90C} $(INNAME) ${GENFLAGS} ${LIBS} -o $(OUTNAME)	
	cp $(OUTNAME) $(BINFOLDER)$(OUTNAME)
clean:
	rm *~
	rm *.o
#########################################################

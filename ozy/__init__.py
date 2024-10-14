import os
OZYPATH = os.path.realpath(__file__).split('/__init__.py')[0]
AMRPATH = os.path.join(OZYPATH,'amr')
PARTPATH = os.path.join(OZYPATH,'part')
VISPATH = os.path.join(OZYPATH,'visualisation')
import sys
sys.path.append(OZYPATH)
sys.path.append(AMRPATH)
sys.path.append(PARTPATH)
sys.path.append(VISPATH)

from ozy.loader import load
from ozy.main import OZY
from ozy.driver import drive

def print_art():
    from art import text2art

    from ozy.__version__ import VERSION
    copywrite = '   (C) 2021 F. Rodriguez Montero'
    version   = '   Version %s' % VERSION

    art =  text2art("Ozymandias","ogre")
    print('\n%s\n%s\n%s\n' % (art, copywrite, version))

def check_mpi_and_print():
    """Check if running under MPI and print art only if rank 0."""
    try:
        # Attempt to import MPI and check rank
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        if rank == 0:
            print_art()

    except ImportError:
        # If mpi4py is not available, just print normally
        print_art()

# Automatically call the function to check MPI and print if necessary
check_mpi_and_print()
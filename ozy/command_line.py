import argparse
import os

import h5py
import numpy as np


def run():
    """Main script for package entry point.
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input file or input directory')
    parser.add_argument('-o', '--output', type=str, help='Output file name')
    parser.add_argument('-lvar', '--linking_var', type=str, help='Object variable used for the linking between halos and galaxies')
    parser.add_argument('-nproc', type=int, help='Set number of processors for fof6d and group property calculation', default=1)
    args = parser.parse_args()
    var_dict = vars(args)
    
    if os.path.isdir(args.input):
        run_multiple_ozy(args.input, var_dict)
        return
    if not os.path.isfile(args.input):
        raise IOError('%s not a valid file!' % args.input)
    if os.path.isdir(args.input):
        run_multiple_ozy(args.input, var_dict)
        return
    ozy_file = False
    # Check if input file is a an ozymandias catalogue.
    try:
        hd = h5py.File(args.input, 'r')
        if 'ozy' in hd.attrs.keys():
            ozy_file = True
        hd.close()
    except:
        pass
    
    if ozy_file:
        open_ozy_file(args.input)
    else:
        run_ozy(args.input, var_dict)
        
def open_ozy_file(infile):
    import IPython

    from .loader import load
    
    obj = load(infile)
    
    IPython.embed(header="OZY file loaded into the 'obj' variable.")
    
def run_ozy(infile, args):
    
    if args['output'] is not None:
        if args['output'].endswith('.hdf5'):
            outfile = args['output']
        else:
            outfile = '%s.hdf5' % args['output']
    elif infile.endswith('.txt'):
        outfile = 'ozy_%s.hdf5' % (infile.split('/')[-1][-9:-4]) # Just to get the index of the RAMSES snap.
    else:
        # The case we only receive the output like 'output_00047'.
        if 'output_' in infile:
            outfile = infile.replace('output_', 'ozy_') + '.hdf5'
        # The case we just receive the index of the output, like '00047'.
        else:
            outfile = 'ozy_%s.hdf5' % infile
            

    if not 'output_' in infile and len(infile) == 5:
        infile = 'output_' + infile
    else:
        raise ImportError('Snapshot format not supported. Please check!')
            
    from .main import OZY
    obj = OZY(infile)
    
    obj.build_HaloMaker(**args)
    obj.save(outfile)
    
def run_multiple_ozy(dir, args):
    import glob

    # Look for output folders.
    infiles = glob.glob('output_*')
    if len(infiles) == 0:
        raise IOError('Could not find output folders in this directory: %s!' % dir)
    
    for f in infiles:
        try:
            run_ozy(f, args)
        except:
            print('Failed on %s' % f)
            pass

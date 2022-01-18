#!/bin/bash

##################################################################################
# SAMPLE SCRIPT FOR RUNNING VELOCIRAPTOR AND TREEFROG IN COSMO RAMSES SIMS       #
#                                                                                #
# This script allows for an easy handling of batch submission of VELOCIraptor    #
# and the following TreeFrog analysis.                                           #
# WARNING: This has only being prepared for Glamdring!                           #
# By: F. Rodriguez Montero 08/01/2022                                            #
##################################################################################

########################################################################################
# Code help display
H=$1
INP=$3
if [ "$H" == "-h" ] || [ "$H" == "-help" ] || [ "$H" == "--help" ] || [ -z $INP ] 2>/dev/null; then
    echo "-------------------------------------------------------------"
    echo "stfrun.sh - Manual:"
    echo
    echo "Purpose:"
    echo "This script automatically executes a modified version of VELOCIraptor"
    echo "and TreeFrog on the particle data of a RAMSES simulation."
    echo
    echo "Usage:"
    echo "stfrun.sh ptype isnap fsnap"    
    echo "Inputs:"
    echo "ptype=$1  # ptype corresponds to the execution type"
    echo "mode=$2  # mode corresponds to the execution mode"
    echo "isnap=$3 # isnap corresponds to the minimum ID to be calculated"
    echo "fsnap=$4 # fsnap corresponds to the maximum ID to be calculated"
    echo
    echo "Available execution types:"
    echo "dm: Using only dark matter particles"
    echo "stars: Using only star particles"
    echo "Available execution modes:"
    echo "all: Run both VELOCIraptor and TreeFrog"
    echo "vel: Only run VELOCIraptor"
    echo "tree: Only run TreeFrog"
    echo 
    exit
fi
########################################################################################


echo "###########################################################################"
echo
echo "              STF FOR RAMSES COSMOLOGICAL HYDRO SIMULATION"
echo "                                                                      v0.3 "
echo "###########################################################################"

# Simulation model or name to be used
ptype=$1

# Running mode
mode=$2

# Initial and final snapshot indexes
isnap=$3
fsnap=$4
nsnaps=`echo $isnap" "$fsnap|awk '{print $2-$1+1}'`

########################################################################################
# Review input arguments 
if ! [ "$isnap" -eq "$isnap" ] 2>/dev/null; then
  isnap=0
fi
if ! [ "$fsnap" -eq "$fsnap" ] 2>/dev/null; then
  fsnap=0
fi
while ! [ "$isnap" -gt "0" ] || ! [ "$fsnap" -ge "$isnap" ]; do
    if ! [ "$isnap" -gt "0" ]; then
        echo "isnap should be higher than 0, introduce a value for isnap:"
        read isnap
    fi
    if ! [ "$fsnap" -ge "$isnap" ]; then
        echo "fsnap should be higher than fsnap, introduce a value for fsnap:"
        read fsnap
    fi
done
while [ "$ptype" != "dm" ] && [ "$ptype" != "stars" ]; do
        echo "ptype not recognised, reintroduce ptype (dm or stars):"
        read ptype
done
while [ "$mode" != "all" ] && [ "$mode" != "vel" ] && [ "$mode" != "tree" ]; do
        echo "mode not recognised, reintroduce mode (all, vel or tree):"
        read mode
done
############################################################

# Number of input files
nfiles=1

# Input directory
indir=..

# Output dir
outdir=.

# Base configuration file for VELOCIraptor to be used
if [ "$ptype" == "dm" ]; then
    paramfile=${outdir}/NH_halo.cfg
    echo "Using DM particles"
else
    paramfile=${outdir}/NH_galaxy.cfg
    echo "Using star particles"
fi

# VELOCIraptor executable
stfexe=${outdir}/stf

# TreeFrog executable
treefrogexe=${outdir}/treefrog

if [ "$mode" == "all" ] || [ "$mode" == "vel" ]; then
    echo "Running VELOCIraptor in the snapshot range $isnap to $fsnap (nsnaps=$nsnaps)"

    for ((j=$isnap; j<=$fsnap; j++))
    do
        jj=`printf "%05d" $j`
        cp $paramfile $ptype.output_$jj.param;
        sed -i 's%Output=OUTNAME%Output='"$outdir"'/'"$ptype"'.c'"$i"'.output_'"$jj"'%g' $outdir/$ptype.output_$jj.param;
        sed -i 's%Snapshot_value=SNVALUE%Snapshot_value='"$j"'%g' $outdir/$ptype.output_$jj.param;
        ifile=`printf "%s/output_%05d" $indir $j`
        $stfexe -i $ifile -s $nfiles -C $outdir/$ptype.output_$jj.param -I 4 -t $jj -o $outdir/$ptype.output_$jj > $outdir/$ptype.output_$jj.log;
    done
fi

# TreeFrog commands
if [ "$mode" == "all" ] || [ "$mode" == "tree" ]; then
    echo "Running TreeFrog in the snapshot range $isnap to $fsnap (nsnaps=$nsnaps)"

    # Number of input VELOCIraptor files (set by number of mpi threads) per snapshot
    numfiles=1

    # Base configuration file for VELOCIraptor to be used
    if [ "$ptype" == "dm" ]; then
        paramfile=${outdir}/NH_halo_tree.cfg
        logname=$outdir/halo_tree.log
    else
        paramfile=${outdir}/NH_galaxy_tree.cfg
        logname=$outdir/galaxy_tree.log
    fi

    rm $outdir/halolist.txt
    for entry in $ptype.output_*.properties
    do
        tmp=${entry#*$ptype.output_}
        jj=${tmp%.*}
        echo $outdir/$ptype.output_$jj >> $outdir/halolist.txt
    done

    $treefrogexe -i $outdir/halolist.txt -s $nsnaps -N $numfiles -o $outdir/$ptype -C $paramfile > $logname
fi

########################################################################################
echo
echo "                          Finished, now closing script"
echo "############################################################################## "
echo "############################################################################## "
echo "############################################################################## "
echo
########################################################################################
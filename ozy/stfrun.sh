#!/bin/bash

##################################################################################
# SAMPLE SCRIPT FOR RUNNING VELOCIRAPTOR AND TREEFROG IN COSMO RAMSES SIMS       #
#                                                                                #
# This script allows for an easy handling of batch submission of VELOCIraptor    #
# and the following TreeFrog analysis.                                           #
# WARNING: This has only being prepared for Glamdring!                           #
# By: F. Rodriguez Montero 08/01/2022                                            #
##################################################################################

echo "###########################################################################"
echo
echo "              STF FOR RAMSES COSMOLOGICAL HYDRO SIMULATION"
echo "                                                                      v0.1 "
echo "###########################################################################"

# Simulation model or name to be used
simname=$1

# Initial and final snapshot indexes
isnap=$2
fsnap=$3
nsnaps=`echo $isnap" "$fsnap|awk '{print $2-$1+1}'`

# Number of input files
nfiles=1

# Input directory
indir=..

# Output dir
outdir=.

# Base configuration file for VELOCIraptor to be used
paramfile=${outdir}/NH_galaxy.cfg

# VELOCIraptor executable
stfexe=${outdir}/stf

# TreeFrog executable
treefrogexe=${outdir}/treefrog

echo "Running VELOCIraptor in the snapshot range $isnap to $fsnap (nsnaps=$nsnaps)"

for ((j=$isnap; j<=$fsnap; j++))
do
    jj=`printf "%05d" $j`
    cp $paramfile $simname.output_$jj.param;
    sed -i 's%Output=OUTNAME%Output='"$outdir"'/'"$simname"'.c'"$i"'.output_'"$jj"'%g' $outdir/$simname.output_$jj.param;
    sed -i 's%Snapshot_value=SNVALUE%Snapshot_value='"$j"'%g' $outdir/$simname.output_$jj.param;
    ifile=`printf "%s/output_%05d" $indir $j`
    $stfexe -i $ifile -s $nfiles -C $outdir/$simname.output_$jj.param -I 4 -t $jj -o $outdir/$simname.output_$jj > $outdir/$simname.output_$jj.log;
done

# TreeFrog commands

echo "Running TreeFrog in the snapshot range $isnap to $fsnap (nsnaps=$nsnaps)"

# Number of input VELOCIraptor files (set by number of mpi threads) per snapshot
numfiles=1

# Base configuration file for VELOCIraptor to be used
paramfile=${outdir}/NH_galaxy_tree.cfg

rm $outdir/halolist.txt
for ((j=$isnap; j<=$fsnap; j++))
do
    jj=`printf "%05d" $j`
    echo $outdir/$simname.output_$jj >> $outdir/halolist.txt
done
$treefrogexe -i $outdir/halolist.txt -s $nsnaps -N $numfiles -o $outdir/$simname.tree -C $paramfile > $outdir/tree.log
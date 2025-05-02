#!/bin/bash -l

# Batch script to run a serial job under SGE.

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=24:0:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=32G

# Request 15 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=10G

# Set the name of the job.
#$ -N FCLS-pure

# Set the working directory to somewhere in your scratch space.  
#  This is a necessary step as compute nodes cannot write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID.
#$ -wd /home/ucapunh/Scratch/FCLS

#$ -t 1-10
#$ -V

export JULIA_PKG_USE_CLI_GIT=true

module load julia/1.11.1

XI=$1
D=$2

echo "Args: $XI $D"

proj=$HOME/Scratch/FCLS

# Your work should be done in $TMPDIR
cd $TMPDIR

# Run the application and put the output into a file called date.txt
julia --project=$proj --startup-file=no -e "using Pkg; Pkg.instantiate()"
julia --project=$proj --startup-file=no -t auto $proj/main_pure.jl $SGE_TASK_ID $XI $D

# Preferably, tar-up (archive) all output files onto the shared scratch area
outpath=${proj}/data/raw/purified/xi=${XI}_D=${D}
mkdir -p $outpath
tar -zcvf ${outpath}/${JOB_ID}.${SGE_TASK_ID}.tar.gz $TMPDIR

# Make sure you have given enough time for the copy to complete!

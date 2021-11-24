#!/bin/bash

#SBATCH --job-name="SnapHiC_example"
#SBATCH -N 2
#SBATCH -o output.log
#SBATCH -e output.errorlog
#SBATCH -p skylake
#SBATCH --ntasks-per-node=15
#SBATCH --mem=128g
#SBATCH -t 12:00:00

### SET THE DIRECTIVES ABOVE IF YOU ARE WORKING IN HPC WITH A JOB SCHEDULER
### AND LOAD THE REQUIRED MODULES IF NEEDED. REQUIRES PYTHON3.6+ WITH 
### MODULES TO BE INSTALLED USING "pip install -r requirements.txt"
module purge
#module load intel/19.4
#module load intelmpi/2018.1.163  #loads mpi
module load openmpi_3.1.4/intel_18.2 #loads mpi

unset PYTHONPATH
export PATH=/nas/longleaf/home/wjzhong/.conda/envs/hic/bin:$PATH #adds python with required packages to path
export PATH="/etc/alternatives":$PATH #adds java to path

export OPENBLAS_MAIN_FREE=1

############################################################################
###                            User Variables                            ###
############################################################################
snapHiC_G_dir="/proj/yunligrp/users/Wujuan/HiC/SnapHiC-G"	#where the snapHiC is located on your system
parallelism="parallel" 					#options are "parallel" "threaded" "singleproc"
number_of_processors=30					#required only if parallelism is set to parallel or threaded
indir="/proj/yunligrp/users/weifangl/SnapHiC/ODC" #directory containing input files (e.g. *.pairs files)
suffix="_indexed_contacts.txt.gz" 			#all input files should have the same suffix. it can be an empty string "", or ".txt"
outdir="/21dayscratch/scr/w/j/wjzhong/HiC/results/odc_output"						#directory where output files will be stored
chrs="2 4" 						#2 integer numbers, the column number of chromosomes in the input files. (e.g. "3 5") starting from 1
pos="3 5" 							#2 integer numbers, the column number of mapped position of read-pairs. (e.g.  "4 6") starting from 1
chrlen="/proj/yunligrp/users/Wujuan/HiC/SnapHiC-G/ext/hg19.chrom.sizes" 		#path to the chrom.sizes file"
genome="hg19"  						#genomeID that will be used for genereation of ".hic" file 
filter_file="/proj/yunligrp/users/Wujuan/HiC/SnapHiC-G/ext/hg19_filter_regions.txt" 	#regions to be filtered, for example due to low mappability
steps="hic interaction postprocess" 			#steps to run the pipeline. Recommended (1) "bin rwr" at first, (2) then  "hic interaction postprocess"
prefix="odc_100"
tss_file="/proj/yunligrp/users/Wujuan/HiC/SnapHiC-G/ext/hg19.refGene.transcript.TSS.061421.txt"  
############################################################################


if [[ "$parallelism" == "parallel" ]]; then
	mpirun -np $number_of_processors python $snapHiC_G_dir/snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --filter-file $filter_file --steps $steps --prefix $prefix --parallel --tss-file $tss_file 
elif [[ "$parallelism" == "threaded" ]]; then
	python $snapHiC_G_dir/snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --filter-file $filter_file --steps $steps --prefix $prefix --threaded -n $number_of_processors --tss-file $tss_file 
else
	python $snapHiC_G_dir/snap.py -i $indir -s $suffix -o $outdir -c $chrs -p $pos -l $chrlen -g $genome --filter-file $filter_file --steps $steps --prefix $prefix --tss-file $tss_file 
fi

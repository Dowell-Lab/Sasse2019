#!/bin/bash
#SBATCH --job-name=DAStk				# Job name
#SBATCH --mail-type=FAIL			 	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=margaret.gruca@colorado.edu	 	# Where to send mail
#SBATCH --nodes=1					# Run on a single node
#SBATCH --ntasks=8   					# Number of CPU (processer cores i.e. tasks)
#SBATCH --mem=20gb					# Memory limit
#SBATCH -p short					# Job queue
#SBATCH --time=2:00:00				  	# Time limit hrs:min:sec
#SBATCH --output=/scratch/Users/magr0763/Gerber/ChIP-seq/e_and_o/%x.%j.out
#SBATCH --error=/scratch/Users/magr0763/Gerber/ChIP-seq/e_and_o/%x.%j.err



##################### LOAD MODULES ######################################

module load python/3.6.3
module load python/3.6.3/numpy/1.14.1
module load python/3.6.3/scipy/0.17.1
#module load python/2.7.14/pandas/0.18.1 ## Not in Fiji for Python3 -- set in local path

########################################################################

###################### SET PATHS #######################################

export PATH=~/.local/bin/:$PATH
BEDS='/scratch/Users/magr0763/Gerber/ChIP-seq/analysis/macs2/'
MOTIF='/scratch/Users/magr0763/Gerber/motifs_HOCOMOCOv11_nGRE_p1e-5'
PROJECT=/scratch/Users/magr0763/Gerber/ChIP-seq/analysis/DAStk/ChIP
ROOTNAME=$1
COND1=$2
COND2=$3
DAStk=/Users/magr0763/.local/lib/python3.6/site-packages/DAStk/

################## PRINT JOB INFO ######################################

printf "\nPath: $PATH"
printf "\nSample ID: $ROOTNAME"
printf "\nMotif Dir: $MOTIF"
printf "\nBeds: $BEDS"
printf "\nDAStk: $DAStk"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

#######################################################################

cd $PROJECT

echo process_atac
date
date

process_atac \
 --prefix Beas2B \
 --threads 8 \
 --atac-peaks ${BEDS}/${ROOTNAME}_summits.bed\
 --motif-path $MOTIF \

echo process_atac done
date
date


#echo differential MD score
#date
#date
#
#
#differential_md_score \
# --prefix $ROOTNAME \
# -1 $COND1 \
# -2 $COND2 \
# --p-value 0.00001 \
# -b
#
#echo MD_score computed
#date
#date

#!/bin/bash
#SBATCH --job-name=bbsplit  # Job name
#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=margaret.gruca@colorado.edu # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=64 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=12:00:00 # Time limit hrs:min:sec
#SBATCH -p short
#SBATCH --mem=200gb # Memory limit
#SBATCH --output=/scratch/Users/magr0763/Gerber/ChIP-seq2/e_and_o/%x.%j.out
#SBATCH --error=/scratch/Users/magr0763/Gerber/ChIP-seq2/e_and_o/%x.%j.err

ref=/scratch/Users/magr0763/genomes/mapping/
project=/scratch/Users/magr0763/Gerber/ChIP-seq2/
fwd1=Ctrl_INPUT1_A_R1
rev1=Ctrl_INPUT1_A_R2
fwd2=Ctrl_INPUT1_B_R1
rev2=Ctrl_INPUT1_B_R2
rootname1=Ctrl_INPUT1_A
rootname2=Ctrl_INPUT1_B
rootname=Ctrl_INPUT1

module load bbmap/38.05

bbsplit.sh -Xmx200g \
in1=${project}/trimmed/${fwd1}.trim.fastq.gz \
in2=${project}/trimmed/${rev1}.trim.fastq.gz \
ref=${ref}/hg38/hg38.fa,${ref}/p_aeruginosa_F5677.fa \
basename=/scratch/Users/magr0763/Gerber/ChIP2/bbsplit/${rootname}_out_%.fastq.gz \
outu1=${project}/bbsplit/${fwd1}.clean.trim.fastq.gz \
outu2=${project}/bbsplit/${rev1}.clean.trim.fastq.gz \
scafstats=${project}/bbsplit/${rootname1}.scafstats.txt \
refstats=${project}/bbsplit/${rootname1}.refstats.txt  \

#!/bin/bash
#SBATCH --job-name=multicov # Job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=margaret.gruca@colorado.edu # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=06:00:00 # Time limit hrs:min:sec
#SBATCH -p short
#SBATCH --mem=10gb # Memory limit
#SBATCH --output=/scratch/Users/magr0763/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/magr0763/e_and_o/%x_%j.err

### Retrieve/use all modules loaded ###
module load bedtools/2.25.0

bamdir1=/scratch/Users/magr0763/Gerber/nf_GRO-seq/mapped/bams/
bam1=${bamdir1}B2B10GROseq-E.sorted.bam
bam2=${bamdir1}B2B10GROseq-E2.sorted.bam
bam3=${bamdir1}B2B10GROseq-D.sorted.bam
bam4=${bamdir1}B2B10GROseq-D2.sorted.bam
bam5=${bamdir1}B2B10GROseq-TE.sorted.bam
bam6=${bamdir1}B2B10GROseq-TE2.sorted.bam
bam7=${bamdir1}B2B10GROseq-TD.sorted.bam
bam8=${bamdir1}B2B10GROseq-TD2.sorted.bam
bam9=${bamdir1}B2B30GROseq-E.sorted.bam
bam10=${bamdir1}B2B30GROseq-E2.sorted.bam
bam11=${bamdir1}B2B30GROseq-D.sorted.bam
bam12=${bamdir1}B2B30GROseq-D2.sorted.bam
bam13=${bamdir1}B2B30GROseq-TE.sorted.bam
bam14=${bamdir1}B2B30GROseq-TE2.sorted.bam
bam15=${bamdir1}B2B30GROseq-TD.sorted.bam
bam16=${bamdir1}B2B30GROseq-TD2.sorted.bam
BED=/scratch/Users/magr0763/Gerber/nf_GRO-seq/analysis/tfit_merged/Sasse2019.tfit_merged_bidirs_noTSS.bed
out=/scratch/Users/magr0763/Gerber/nf_GRO-seq/analysis/tfit_merged/Sasse2019.counts_tfit_merged_bidirs_noTSS.bed

multiBamCov \
 -bams ${bam1} ${bam2} ${bam3} ${bam4} ${bam5} ${bam6} ${bam7} ${bam8} ${bam9} ${bam10} ${bam11} ${bam12} ${bam13} ${bam14} ${bam15} ${bam16} \
 -bed $BED \
 > ${out}

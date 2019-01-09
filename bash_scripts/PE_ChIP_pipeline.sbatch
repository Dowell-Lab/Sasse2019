#!/bin/bash
#SBATCH --job-name=PE-ChIP # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=margaret.gruca@colorado.edu # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=32 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH -p long
#SBATCH --mem=200gb # Memory limit
#SBATCH --output=/scratch/Users/magr0763/Gerber/ChIP-seq2/e_and_o/%x_%j.out
#SBATCH --error=/scratch/Users/magr0763/Gerber/ChIP-seq2/e_and_o/%x_%j.err

#This script will take you from a raw, unedited fastq file to tdf

ref=/scratch/Users/magr0763/genomes/mapping/hg38
project=/scratch/Users/magr0763/Gerber/ChIP-seq2/
fwd1=GR-IA1_shGR_V2_A_R1
rev1=GR-IA1_shGR_V2_A_R2
fwd2=GR-IA1_shGR_V2_B_R1
rev2=GR-IA1_shGR_V2_B_R2
rootname1=GR-IA1_shGR_V2_A
rootname2=GR-IA1_shGR_V2_B
rootname=GR-IA1_shGR_V2
#cntrl_1=
#cntrl_2=


pwd; hostname; date
date


#Quality check of fastqc file before trimming, includes reverse complements

module load fastqc/0.11.5
fastqc ${project}/fastq/${fwd1}.fastq.gz -o ${project}/fastq/FastQC/
fastqc ${project}/fastq/${rev1}.fastq.gz -o ${project}/fastq/FastQC/
fastqc ${project}/fastq/${fwd2}.fastq.gz -o ${project}/fastq/FastQC/
fastqc ${project}/fastq/${rev2}.fastq.gz -o ${project}/fastq/FastQC/

echo pre-trim qual
date
date

#Trimming fastqc file reads --> can use trimmomatic (below) or trimgalore

module load java/1.8.0_101
module load bbmap/38.05 
bbduk.sh -Xmx20g \
overwrite= t \
in=${project}/fastq/${fwd1}.fastq.gz \
in2=${project}/fastq/${rev1}.fastq.gz \
out=${project}/trimmed/${fwd1}.trim.fastq.gz \
out2=${project}/trimmed/${rev1}.trim.fastq.gz \
ref=/scratch/Users/magr0763/bbmap/resources/adapters.fa \
ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
maq=10 minlen=20 \
stats=${project}/trimmed/${rootname1}.trimstats.txt \
refstats=${project}/trimmed/${rootname1}.refstats.txt \
ehist=${project}/trimmed/${rootname1}.ehist.txt \
tpe tbo \
#maxgc=1 \
#literal=CCCGTGTTGAGTCAAATTAAGCCGCAGGCTCCACTCCTGGTGGTGCCCTT \

module load java/1.8.0_101
module load bbmap/38.05
bbduk.sh -Xmx20g \
overwrite= t \
in=${project}/fastq/${fwd2}.fastq.gz \
in2=${project}/fastq/${rev2}.fastq.gz \
out=${project}/trimmed/${fwd2}.trim.fastq.gz \
out2=${project}/trimmed/${rev2}.trim.fastq.gz \
ref=/scratch/Users/magr0763/bbmap/resources/adapters.fa \
ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
maq=10 minlen=20 \
stats=${project}/trimmed/${rootname2}.trimstats.txt \
refstats=${project}/trimmed/${rootname2}.refstats.txt \
ehist=${project}/trimmed/${rootname2}.ehist.txt \
tpe tbo \
#maxgc=1 \
#literal=CCCGTGTTGAGTCAAATTAAGCCGCAGGCTCCACTCCTGGTGGTGCCCTT \

echo bbduk trim done
date
date


#Quality check trimmed reads --> if first QC bad, check these before moving forward and adjust parameters if needed

module load fastqc/0.11.5
fastqc ${project}/trimmed/${fwd1}.trim.fastq.gz -o ${project}/trimmed/FastQC/
fastqc ${project}/trimmed/${rev1}.trim.fastq.gz -o ${project}/trimmed/FastQC/
fastqc ${project}/trimmed/${fwd2}.trim.fastq.gz -o ${project}/trimmed/FastQC/
fastqc ${project}/trimmed/${rev2}.trim.fastq.gz -o ${project}/trimmed/FastQC/

echo post-trim qual
date
date

### Splitting reads between two reference genomes to decontaminate

module load bbmap/38.05

bbsplit.sh -Xmx200g \
in1=${project}/trimmed/${fwd1}.trim.fastq.gz \
in2=${project}/trimmed/${rev1}.trim.fastq.gz \
ref=/scratch/Users/magr0763/genomes/mapping/hg38/hg38.fa,/scratch/Users/magr0763/genomes/mapping/p_aeruginosa_F5677.fa \ \
basename=/scratch/Users/magr0763/Gerber/ChIP-seq2/bbsplit/${rootname1}_out_%_#.fastq.gz \
outu1=${project}/bbsplit/${fwd1}.clean.trim.fastq.gz \
outu2=${project}/bbsplit/${rev1}.clean.trim.fastq.gz \
scafstats=${project}/bbsplit/${rootname1}.scafstats.txt \
refstats=${project}/bbsplit/${rootname1}.refstats.txt  \

bbsplit.sh -Xmx200g \
in1=${project}/trimmed/${fwd2}.trim.fastq.gz \
in2=${project}/trimmed/${rev2}.trim.fastq.gz \
ref=/scratch/Users/magr0763/genomes/mapping/hg38/hg38.fa,/scratch/Users/magr0763/genomes/mapping/p_aeruginosa_F5677.fa \
basename=/scratch/Users/magr0763/Gerber/ChIP-seq2/bbsplit/${rootname2}_out_%_#.fastq.gz \
outu1=${project}/bbsplit/${fwd2}.clean.trim.fastq.gz \
outu2=${project}/bbsplit/${rev2}.clean.trim.fastq.gz \
scafstats=${project}/bbsplit/${rootname2}.scafstats.txt \
refstats=${project}/bbsplit/${rootname2}.refstats.txt  \

echo split
date
date

#Map trimmed fastq files to reference genome --> sam files include all reads

module load samtools/1.3.1
module load hisat2/2.1.0

hisat2 -p 32 \
 --very-sensitive \
 --no-spliced-alignment \
 -x ${ref}/HISAT2_indexes/genome \
 -1 ${project}/bbsplit/${rootname1}_out_hg38_1.fastq.gz \
 -2 ${project}/bbsplit/${rootname1}_out_hg38_2.fastq.gz \
 > ${project}/mapped/sams/${rootname1}.trim.sam \
 2>${project}/mapped/sams/${rootname1}.trim.stderr 

hisat2 -p 32 \
 --very-sensitive \
 --no-spliced-alignment \
 -x ${ref}/HISAT2_indexes/genome \
 -1 ${project}/bbsplit/${rootname2}_out_hg38_1.fastq.gz \
 -2 ${project}/bbsplit/${rootname2}_out_hg38_2.fastq.gz \
 > ${project}/mapped/sams/${rootname2}.trim.sam \
 2>${project}/mapped/sams/${rootname2}.trim.stderr

cat ${project}/mapped/sams/${rootname1}.trim.sam <(grep -v '^@' ${project}/mapped/sams/${rootname2}.trim.sam) > ${project}/mapped/sams/${rootname}.cat.trim.sam

echo mapped, sam
date
date

#Word count sam files, convert sam files to bam files (compressed, binary sam files)

wc -l ${project}/mapped/sams/${rootname}.cat.trim.sam > ${project}/mapped/sams/${rootname}.cat.trim.sam.wc

samtools view -S -b -o ${project}/mapped/bams/${rootname}.cat.trim.bam ${project}/mapped/sams/${rootname}.cat.trim.sam 2> ${project}/mapped/bams/${rootname}.cat.trim.bam.err

samtools flagstat ${project}/mapped/bams/${rootname}.cat.trim.bam > ${project}/mapped/bams/${rootname}.cat.trim.bam.flagstat 2> ${project}/mapped/bams/${rootname}.cat.trim.bam.flagstat.err

echo bam
date
date

#Sort bam files

samtools sort -m 200G ${project}/mapped/bams/${rootname}.cat.trim.bam > ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam

samtools flagstat ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam > ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam.flagstat 2> ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam.flagstat.err

echo sorted.bam
date
date

echo index bam
date
date

#Index sorted bam files for use in mutlicov OR if you decide to sort out a specific region (e.g. use samtools view to select only chromosome 1)

samtools index ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam.bai

echo indexed bam
date
date


module load preseq

#c_curve : can use inputs from file types [bed/bam] and will plot the estimated complexity of a sample. Specify -B for sorted bam files, -P for paired end read files

preseq c_curve -B \
 -o ${project}/qc/preseq/${rootname}.c_curve.txt \
 ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam

echo c_curve
echo complexity
date

#lc_extrap : can use inputs from from tyles [bed/bam] and will estimate future yields for complexity if the sample is sequenced at higher read depths. Specify -B for sorted bam files, -P for paired end reads

preseq lc_extrap -B \
 -o ${project}/qc/preseq/${rootname}.lc_extrap.txt \
 ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam

echo lc_extrap
echo future yield
date

#analyze read distributions using RSeQC -- will give you number of reads over different regions of genome dependent on the input annotation file (typically will be over genes, but you can imagine the applicability to eRNA calls eventually for quick comparison of relative eRNA coverage between samples)

module load python/2.7.14/rseqc

read_distribution.py  -i ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam -r /scratch/Users/magr0763/genomes/GeneAnnotations/hg38_refseq.bed > ${project}/qc/rseqc/${rootname}.read_dist.txt

echo rseqc
date
date

#BedGraph generator -- generating positive and negative strand mapping; awk coverts values to negative values for negative strand mapping

module load bedtools/2.25.0
genomeCoverageBed -bg -strand + -ibam ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam -g hg38 > ${project}/mapped/bedgraph/${rootname}.cat.trim.unsorted.BedGraph

sortBed -i ${project}/mapped/bedgraph/${rootname}.cat.trim.unsorted.BedGraph > ${project}/mapped/bedgraph/${rootname}.cat.trim.BedGraph

echo BedGraph
date
date

#Read count correcting (rcc) -- eliminating sequencing bias 

module load python/2.7.14
python /Users/magr0763/mappingandvis/readcountcorrectBG.py ${project}/mapped/bedgraph/${rootname}.cat.trim.BedGraph ${project}/mapped/bams/sorted/${rootname}.cat.trim.sorted.bam.flagstat ${project}/mapped/bedgraph/${rootname}.cat.trim.rcc.BedGraph

echo readcountcorrectedbedgraph
date
date

#Generate tdfs (compressed bedgraphs) viewable in IGV

/opt/igvtools/2.3.75/igvtools toTDF ${project}/mapped/bedgraph/${rootname}.cat.trim.rcc.BedGraph ${project}/mapped/tdfs/${rootname}.cat.trim.rcc.tdf /Users/magr0763/igv/genomes/sizes/hg38.chrom.sizes
echo tdf
date
date


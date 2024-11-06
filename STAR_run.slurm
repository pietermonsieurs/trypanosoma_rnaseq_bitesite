#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_tryps

module load STAR/2.7.10b-GCC-11.3.0

## parameter settings
threads=28
genome_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/ref_genome/STAR/
output_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/results/star/
gtf_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/ref_genome/Mmusculus_Tbrucei_combined.gtf

## get the fastq_file_R2
echo $fastq_file_R1
fastq_file_R2=${fastq_file_R1/_R1.fastq.gz/_R2.fastq.gz}
file_prefix=${fastq_file_R1%_R1.fastq.gz}
file_prefix=${file_prefix##*/}
out_prefix=${output_dir}/${file_prefix}_

echo $fastq_file_R2
echo $out_prefix

## run STAR
STAR \
    --readFilesCommand zcat \
    --runThreadN ${threads} \
    --genomeDir ${genome_dir} \
    --readFilesIn ${fastq_file_R1} ${fastq_file_R2} \
    --quantMode GeneCounts \
    --outFileNamePrefix ${out_prefix} \
    --sjdbGTFfile ${gtf_file} \
    --outSAMtype BAM SortedByCoordinate


## run for all R1 fastq files
# cd /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/fastq
# for fastq_file_R1 in ${PWD}/*_R1.fastq.gz; do sbatch --export=fastq_file_R1=${fastq_file_R1} /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/bin/STAR_run.slurm; done
#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --time=2:00:00
#SBATCH -A ap_itg_tryps

module load BWA

# parameter settings
fastq_screen_bin=/user/antwerpen/205/vsc20587/data/software/FastQ-Screen-0.14.1/fastq_screen
threads=2
out_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/results/fastq_screen/

## test fastq file to optimize workflow
# fastq_file=/user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq/data/HG7FFDSXY_104099-001-001R_GATTGCTC-GTGAAGTG_L001_R1.fastq.gz

$fastq_screen_bin \
    --conf /data/antwerpen/205/vsc20587/software/FastQ-Screen-0.14.1/fastq_screen.conf \
    --aligner bwa \
    --bwa "-k 50" \
    --threads $threads \
    --outdir $out_dir \
    ${fastq_file}

## run for all fastq_files
# cd /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/fastq/
# for fastq_file in *R1.fastq.gz; do sbatch --export=fastq_file=${fastq_file} /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/bin/novogene_fastq_screen.slurm; done


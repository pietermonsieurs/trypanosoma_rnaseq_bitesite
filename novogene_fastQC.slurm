#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=4
#SBATCH --time=0:50:00
#SBATCH --job-name=FastQC
#SBATCH -A ap_itg_tryps

module load BWA
module load BioTools
module load Java


# parameter settings
out_dir=/user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/results/fastqc/
fastqc_bin=/user/antwerpen/205/vsc20587/data/software/FastQC/fastqc
threads=4

# run FastQC based on the worker moduel 
$fastqc_bin --outdir ${out_dir} --threads $threads $fastq_file




## run atools
# cd /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/fastq/
# for fastq_file in *_R1.fastq.gz; do sbatch --export=fastq_file=${PWD}/${fastq_file} /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/bin/novogene_fastQC.slurm; done


## run multiQC
# cd /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/results/fastqc
# module load Python/3.8.3-intel-2020a
# export PYTHONPATH="/user/antwerpen/205/vsc20587/data/software/python_lib/lib/python3.8/site-packages/:${PYTHONPATH}"
# /user/antwerpen/205/vsc20587/data/software/python_lib/bin/multiqc ./



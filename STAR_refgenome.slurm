#!/bin/bash -l

#SBATCH --ntasks=1 --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH -A ap_itg_tryps

module load STAR/2.7.10b-GCC-11.3.0

## run the genome indexing tool of STAR
STAR \
    --runThreadN 28 \
    --runMode genomeGenerate \
    --genomeDir /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/ref_genome/STAR/ \
    --genomeFastaFiles /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/ref_genome/Mmusculus_Tbrucei_combined.fa \
    --sjdbGTFfile /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/ref_genome/Mmusculus_Tbrucei_combined.gtf \
    --sjdbOverhang 100


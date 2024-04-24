## check and load modules
module spider STAR
module load STAR/2.7.10b-GCC-11.3.0

## get the reference genomes
wget https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gff3/mus_musculus/Mus_musculus.GRCm39.111.gff3.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.gtf.gz
wget https://tritrypdb.org/common/downloads/release-67/TbruceiTREU927/fasta/data/TriTrypDB-67_TbruceiTREU927_Genome.fasta
wget https://tritrypdb.org/common/downloads/release-67/TbruceiTREU927/gff/data/TriTrypDB-67_TbruceiTREU927.gff
/user/antwerpen/205/vsc20587/data/software/gffread/gffread TriTrypDB-67_TbruceiTREU927.gff -T -o TriTrypDB-67_TbruceiTREU927.gtf


## concatenate the fasta and gff file
gzip TriTrypDB-67_TbruceiTREU927_Genome.fasta
cat Mus_musculus.GRCm39.dna.primary_assembly.fa.gz TriTrypDB-67_TbruceiTREU927_Genome.fasta.gz > Mmusculus_Tbrucei_combined.fa.gz

gzip TriTrypDB-67_TbruceiTREU927.gtf
cat Mus_musculus.GRCm39.111.gtf.gz TriTrypDB-67_TbruceiTREU927.gtf.gz >  Mmusculus_Tbrucei_combined.gtf.gz

gunzip Mmusculus_Tbrucei_combined.fa.gz
gunzip Mmusculus_Tbrucei_combined.gtf.gz

## run the genome indexing tool of STAR
STAR \
    --runThreadN 2 \
    --runMode genomeGenerate \
    --genomeDir /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/ref_genome/STAR/ \
    --genomeFastaFiles /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/ref_genome/Mmusculus_Tbrucei_combined.fa \
    --sjdbGTFfile /user/antwerpen/205/vsc20587/scratch/trypanosoma_rnaseq_bitesite/data/ref_genome/Mmusculus_Tbrucei_combined.gtf \
    --sjdbOverhang 100
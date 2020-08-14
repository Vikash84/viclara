#!/bin/bash

#SBATCH --partition batch
#SBATCH -w compute04
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=nf-viclara
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=J.Juma@cgiar.org

# clean environment
module purge

#
source ~/miniconda3/etc/profile.d/conda.sh
conda activate

# load module
module load nextflow/20.04.1
export NXF_OPTS='-Xms512M -Xmx4G'
#module load fastqc/0.11.7
#module load trimmomatic/0.38
#module load bbmap/38.67
#module load bwa/0.7.15
#module load bowtie2/2.2.5
#module load multiqc/1.4
#module load kraken/2.0.8-beta
# /home/jjuma/Work/Sam_Oyola/mNGS_06022020/NCDV_1
# /home/jjuma/PhD_RVF2019/vir-meta/testdataset

nextflow viclara.nf run --reads "/home/jjuma/Work/Sam_Oyola/mNGS_06022020/NCDV_1/*_R{1,2}_001.fastq.gz" \
  -profile conda \
  --reference /var/scratch/jjuma/viral_metagenomics/references_db/Galgal6/GCF_000002315.6_GRCg6a_genomic.fna \
  --removePhiX \
  --hg /var/scratch/jjuma/viral_metagenomics/references_db/GRch38/GCF_000001405.39_GRCh38.p13_genomic.fna \
  --aligner bwa \
  --adapters /export/apps/trimmomatic/0.38/adapters/NexteraPE-PE.fa \
  --classify \
  --krn2_db /var/scratch/jjuma/kraken2_db/viral \
  --krn2_task download-taxonomy \
  --krn2_library viral \
  --outdir /var/scratch/jjuma/viclara-nf_output \
  --work-dir /var/scratch/jjuma/ \
  -resume

#!/bin/bash
#SBATCH --job-name=nextflow_genomeQC   # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lsh2503120@lshtm.ac.uk # Where to send mail
#SBATCH --ntasks=1                    # Run on a single core
#SBATCH --mem=30gb                     # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
pwd; hostname; date

echo "Running nextflow genome QC pipleine"

source ~/.bashrc

conda deactivate 
conda activate nextflow 

nextflow run nf-core/genomeqc \
 --input ./samplesheet.csv \
  --outdir ./results \
   -profile conda
    -resume 

date
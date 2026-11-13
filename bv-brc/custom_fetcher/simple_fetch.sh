#!/bin/bash
#SBATCH --job-name=retrieve_pneumococcal_genomes    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lsh2503120@lshtm.ac.uk # Where to send mail
#SBATCH --ntasks=1                    # Run on a single core
#SBATCH --mem=15gb                     # Job memory request
#SBATCH --time=12:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
pwd; hostname; date

echo "Running pneumococcal genome retrieval"


python3 ./fetch_genomes.py download ./accession_genomes.acc ~/pneumococcal-genomes


date


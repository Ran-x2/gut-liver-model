#!/bin/bash
#SBATCH --job-name=rxr456_RNA_job
#SBATCH --output=D3INTL_job_output.txt
#SBATCH --error=D3INTL_job_error.txt
#SBATCH --time=48:00:00
#SBATCH --ntasks=48
#SBATCH --mem=128G
#SBATCH --mail-user=rxr456@case.edu
#SBATCH --mail-type=ALL 

module load Java
module load R
module load STAR
module load GCC
export PATH=$PATH:/home/rxr456/RSEM/RSEM-1.3.3 
echo 'export PATH=/home/rxr456/perl5/bin:$PATH' >> ~/.bashrc
source ~/.bashrc  
echo 'export PERL5LIB=/home/rxr456/perl5/lib/perl5:$PERL5LIB' >> ~/.bashrc
source ~/.bashrc
source /home/rxr456/miniconda3/etc/profile.d/conda.sh #make sure you have the conda and python correctly installed
STORAGE_DIR="/mnt/vstor/SOM_PATH_DKB50/members/rxr456"
TRIM_DIR="$STORAGE_DIR"/trim

# Loop through the input files in $TRIM_DIR
for SAMPLE_DIR in TRIM_DIR; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR" "_trimmed_1.fastq")
    R1="$TRIM_DIR/${SAMPLE_NAME}_trimmed_1.fastq"
    R2="$TRIM_DIR/${SAMPLE_NAME}_trimmed_2.fastq"
    
    # Check if the sample has already been processed
        echo "Running STAR and RSEM on $SAMPLE_NAME..."
        echo "$R1"
        echo "$R2"
        # Run STAR and RSEM

        # how to find STAR path: where STAR. and substitue the code
        rsem-calculate-expression --star --star-path /usr/local/easybuild_allnodes/software/STAR/2.7.9a-GCC-11.2.0/bin \
        --paired-end --estimate-rspd --append-names --output-genome-bam -p 48 "$R1" "$R2" \
        /home/rxr456/RSEM_ref/RSEM_ref "$STORAGE_DIR/RSEM_quant/$SAMPLE_NAME"        
done



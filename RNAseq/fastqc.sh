module load FastQC

for SAMPLE_DIR in "$RAW_DIR"/*.fq.gz; do
    fastqc -o $STORAGE_DIR/fastqc/ -t 48 $SAMPLE_DIR
    done
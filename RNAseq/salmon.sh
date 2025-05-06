RAW_DATA_DIR="/mnt/h/RNAseq/raw"
INDEX_DIR="/mnt/h/RNAseq/index/transcripts_index"
OUTPUT_DIR="/mnt/h/RNAseq/aligned"

mkdir -p "$OUTPUT_DIR"
for SAMPLE_DIR in "$RAW_DATA_DIR"/D*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    R1=$(ls "$SAMPLE_DIR"/*_1.fq.gz)
    R2=$(ls "$SAMPLE_DIR"/*_2.fq.gz)
    SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/$SAMPLE_NAME"

    echo "Running Salmon on $SAMPLE_NAME..."
    echo "$R1"
    echo "$R2"
    salmon quant -i "$INDEX_DIR" -l A -1 "$R1" -2 "$R2" -p 8 --validateMappings -o "$SAMPLE_OUTPUT_DIR"
done

echo "Salmon quantification complete for all samples!"
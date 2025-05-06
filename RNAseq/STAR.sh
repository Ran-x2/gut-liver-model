#use your RSEM ref is fine
STAR --genomeDir /home/rxr456/STAR_ref \ 
--readFilesIn "$TRIM_DIR/D4_INT_HEP_C_2_trimmed_1.fastq" "$TRIM_DIR/D4_INT_HEP_C_2_trimmed_2.fastq"  \
--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
--twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM \
--runThreadN 48 --outFileNamePrefix "$STORAGE_DIR/STAR/D4_INT_HEP_5OP_1"
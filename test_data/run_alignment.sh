
LAST_BIN=/vol/projects/pbohn/tools/last-1450/bin
BAM2FQ=/vol/projects/pbohn/tools/pb_functions/bam2fq/target/release/bam2fq
TEST_DIR=/vol/projects/pbohn/RNA_probing/run121_panel_RNA4/notebooks/data_processing/scripts/bamButler/test_data
LAST_REF_PREFIX=/vol/projects/pbohn/RNA_probing/run121_panel_RNA4/references/individual/HCV_IRES_LAST/HCV_IRES
REF=/vol/projects/pbohn/RNA_probing/run121_panel_RNA4/references/individual/HCV_IRES.fa

$BAM2FQ --output-prefix $TEST_DIR/unaligned --output-stats-prefix to_be_deleted_ $TEST_DIR/unaligned.bam
rm to_be_deleted_*
$LAST_BIN/lastal -Qkeep -m20 $LAST_REF_PREFIX $TEST_DIR/unaligned.fastq.gz | $LAST_BIN/last-split -m1 > $TEST_DIR/aligned.maf
$LAST_BIN/maf-convert sam $TEST_DIR/aligned.maf > $TEST_DIR/aligned.sam
rm $TEST_DIR/aligned.maf
samtools view -S -bt $REF.fai $TEST_DIR/aligned.sam > $TEST_DIR/aligned.bam
rm $TEST_DIR/aligned.sam
samtools calmd --output-fmt BAM $TEST_DIR/aligned.bam $REF > $TEST_DIR/aligned_MD.bam
rm $TEST_DIR/aligned.bam
samtools sort -O bam $TEST_DIR/aligned_MD.bam > $TEST_DIR/aligned_MD_sorted.bam
rm $TEST_DIR/aligned_MD.bam
samtools index $TEST_DIR/aligned_MD_sorted.bam
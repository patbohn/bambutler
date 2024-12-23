// tests/integration_test.rs
use bambutler::{create_read_index, process_bam_file};
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use tempfile::TempDir;

mod common;
use common::{count_cigar_ops, test_data_dir};

#[test]
fn test_create_read_index() -> Result<()> {
    let unaligned_path = test_data_dir().join("unaligned.bam");
    let index = create_read_index(&unaligned_path)?;
    
    // Create a Vec<u8> for the test read name
    let test_read_name = b"5bddecba-5f37-4b05-b3f3-170e77949d6f".to_vec();
    assert!(index.contains_key(&test_read_name), "Index should contain test read");
    
    if let Some(read) = index.get(&test_read_name) {
        assert!(!read.sequence.is_empty(), "Read should have sequence");
        assert!(!read.qualities.is_empty(), "Read should have quality scores");
    }

    Ok(())
}

#[test]
fn test_process_bam_file() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let test_dir = test_data_dir();
    
    let aligned_path = test_dir.join("aligned_MD_sorted.bam");
    let unaligned_path = test_dir.join("unaligned.bam");
    
    let unaligned_index = create_read_index(&unaligned_path)?;
    
    let stats = process_bam_file(
        &aligned_path,
        &unaligned_index,
        &temp_dir.path().to_path_buf(),
    )?;
    
    let output_path = temp_dir.path().join(
        aligned_path.file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .replace(".bam", "_converted.bam")
    );
    
    let mut output_bam = bam::Reader::from_path(&output_path)?;
    let mut record = bam::Record::new();
    
    while let Some(result) = output_bam.read(&mut record) {
        result?;
        assert_eq!(
            count_cigar_ops(&record, 5),
            0,
            "Output BAM should not contain hard clips"
        );
        
        assert!(!record.seq().as_bytes().is_empty(), "Record should have sequence");
        assert!(!record.qual().is_empty(), "Record should have quality scores");
        
        let qname = record.qname().to_vec();
        if let Some(unaligned) = unaligned_index.get(&qname) {
            for (tag, _) in &unaligned.tags {
                assert!(
                    record.aux(tag.as_slice()).is_ok(),
                    "Missing tag {} from unaligned read",
                    String::from_utf8_lossy(tag)
                );
            }
        }
    }

    assert!(stats.reads_processed > 0, "Should have processed some reads");
    assert!(
        stats.reads_modified <= stats.reads_processed,
        "Modified reads should not exceed total reads"
    );

    Ok(())
}
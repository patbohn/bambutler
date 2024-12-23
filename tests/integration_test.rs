// tests/integration_test.rs
use bambutler::{create_read_index, process_bam_file};
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use tempfile::TempDir;
use rust_htslib::bam::record::Aux;

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

    // Test mandatory tags (mv as B:c array and ts as integer)
    // plus an optional tag (pi)
    let transfer_tags = vec![
        "mv".to_string(),  // movement tag (B:c type)
        "ts".to_string(),  // timestamp (integer)
        "pi".to_string()   // optional tag
    ];
    
    let stats = process_bam_file(
        &aligned_path,
        &unaligned_index,
        &temp_dir.path().to_path_buf(),
        &transfer_tags
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
                       // Check mandatory tags with correct types
           if let Ok(mv_tag) = record.aux(b"mv") {
            match mv_tag {
                Aux::ArrayU8(_) => (), // B:c type
                _ => panic!("mv tag should be of type B:c (ArrayU8)")
            }
        }
        
        if let Ok(ts_tag) = record.aux(b"ts") {
            match ts_tag {
                Aux::I32(_) | Aux::U32(_) | Aux::I16(_) | Aux::U16(_) | 
                Aux::I8(_) | Aux::U8(_) => (), // Integer types
                _ => panic!("ts tag should be an integer type")
            }
        }
        
        // For pi tag, just check if it exists when it should
        if unaligned.tags.iter().any(|(t, _)| t == b"pi") {
            assert!(record.aux(b"pi").is_ok(), "pi tag should be transferred when present");
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
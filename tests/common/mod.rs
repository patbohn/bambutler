// tests/common/mod.rs
use rust_htslib::bam::record::Record;
use std::path::PathBuf;

pub fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_data")
}

pub fn count_cigar_ops(record: &Record, op_type: u32) -> usize {
    record.cigar()
        .iter()
        .filter(|op| op.char() as u32 == op_type)
        .count()
}

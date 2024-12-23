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

pub fn compare_tags(record1: &Record, record2: &Record) -> bool {
    let tags1: Vec<_> = record1.aux_iter()
        .filter_map(Result::ok)
        .collect();
    let tags2: Vec<_> = record2.aux_iter()
        .filter_map(Result::ok)
        .collect();
    
    tags1.len() == tags2.len() && tags1.iter().all(|(tag1, val1)| {
        record2.aux(tag1)
            .map(|val2| val1 == &val2)
            .unwrap_or(false)
    })
}
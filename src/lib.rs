use rust_htslib::{bam, bam::record::Record, bam::Read};
use rustc_hash::FxHashMap;
use anyhow::{Context, Result};
use std::path::PathBuf;
use log::{info, warn};
use rayon::prelude::*;

// Move your structs and type definitions here
#[derive(Clone)]
pub struct UnalignedRead {
    pub sequence: Vec<u8>,
    pub qualities: Vec<u8>,
    pub tags: Vec<(String, String)>,
}

// Custom enum to store tag values
#[derive(Clone)]
pub enum TagValue {
    Char(i8),
    Int(i32),
    Float(f32),
    String(String),
    Array(Vec<u8>),
}

impl From<Aux<'_>> for TagValue {
    fn from(aux: Aux) -> Self {
        match aux {
            Aux::Char(c) => TagValue::Char(c),
            Aux::Int8(i) => TagValue::Int(i32::from(i)),
            Aux::UInt8(i) => TagValue::Int(i32::from(i)),
            Aux::Int16(i) => TagValue::Int(i32::from(i)),
            Aux::UInt16(i) => TagValue::Int(i32::from(i)),
            Aux::Int32(i) => TagValue::Int(i),
            Aux::UInt32(i) => TagValue::Int(i32::try_from(i).unwrap_or(0)),
            Aux::Float(f) => TagValue::Float(f),
            Aux::String(s) => TagValue::String(String::from_utf8_lossy(s).into_owned()),
            Aux::Array(array) => TagValue::Array(array.to_owned()),
            _ => TagValue::String("".to_string()),
        }
    }
}


#[derive(Default)]
pub struct Stats {
    pub reads_processed: usize,
    pub reads_modified: usize,
    pub reads_missing: usize,
}

impl Stats {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Create an index of reads from the unaligned BAM file
pub fn create_read_index(path: &PathBuf) -> Result<FxHashMap<Vec<u8>, UnalignedRead>> {
    info!("Creating index from unaligned BAM file...");
    let mut bam = bam::Reader::from_path(path)?;
    let mut index = FxHashMap::default();
    let mut buffer = Record::new();

    while let Some(result) = bam.read(&mut buffer) {
        result?;
        let name = buffer.qname().to_vec();
        
        // Extract tags
        let tags = buffer
            .aux_iter()
            .filter_map(|result| result.ok())
            .map(|(tag, value)| (
                String::from_utf8_lossy(tag).to_string(),
                TagValue::from(value)
            ))
            .collect();

        index.insert(name, UnalignedRead {
            sequence: buffer.seq().as_bytes(),
            qualities: buffer.qual().to_vec(),
            tags,
        });
    }

    info!("Indexed {} reads", index.len());
    Ok(index)
}

/// Convert CIGAR string from hard clips to soft clips
pub fn convert_cigar(cigar: &[u8]) -> Vec<u8> {
    cigar.iter()
        .map(|&op| if op == 5 { 4 } else { op })
        .collect()
}

/// Process a single BAM file
pub fn process_bam_file(
    input_path: &PathBuf,
    unaligned_index: &FxHashMap<Vec<u8>, UnalignedRead>,
    output_dir: &PathBuf,
) -> Result<Stats> {
    let mut stats = Stats::new();
    let mut input = bam::Reader::from_path(input_path)?;
    
    // Create output path
    let output_name = input_path
        .file_name()
        .context("Invalid input filename")?
        .to_str()
        .context("Invalid UTF-8 in filename")?
        .replace(".bam", "_converted.bam");
    let output_path = output_dir.join(output_name);
    
    let header = bam::Header::from_template(input.header());
    let mut output = bam::Writer::from_path(&output_path, &header, bam::Format::Bam)?;
    
    let mut buffer = Record::new();
    
    while let Some(result) = input.read(&mut buffer) {
        result?;
        stats.reads_processed += 1;

        if stats.reads_processed % 100_000 == 0 {
            info!("Processed {} reads...", stats.reads_processed);
        }

        let name = buffer.qname().to_vec();
        let has_hard_clips = buffer
            .cigar()
            .iter()
            .any(|op| op.0 == 5);

        match unaligned_index.get(&name) {
            Some(unaligned) => {
                let mut new_record = Record::new();
                
                // Copy basic fields
                new_record.set_qname(&name);
                new_record.set_pos(buffer.pos());
                new_record.set_mapq(buffer.mapq());
                new_record.set_flags(buffer.flags());
                
                // Convert CIGAR if needed
                if has_hard_clips {
                    new_record.set_cigar(&convert_cigar(buffer.cigar_bytes()));
                    stats.reads_modified += 1;
                } else {
                    new_record.set_cigar(buffer.cigar());
                }

                // Set sequence and qualities from unaligned read
                new_record.set_seq(&unaligned.sequence);
                new_record.set_qual(&unaligned.qualities);

                // Transfer original tags
                for result in buffer.aux_iter() {
                    if let Ok((tag, value)) = result {
                        match TagValue::from(value) {
                            TagValue::Char(c) => new_record.push_aux(tag, c)?,
                            TagValue::Int(i) => new_record.push_aux(tag, i)?,
                            TagValue::Float(f) => new_record.push_aux(tag, f)?,
                            TagValue::String(s) => new_record.push_aux(tag, s.as_bytes())?,
                            TagValue::Array(arr) => new_record.push_aux(tag, &arr)?,
                        }
                    }
                }

                // Add new tags from unaligned read
                for (tag, value) in &unaligned.tags {
                    if !buffer.aux(tag.as_bytes()).is_ok() {
                        match value {
                            TagValue::Char(c) => new_record.push_aux(tag.as_bytes(), *c)?,
                            TagValue::Int(i) => new_record.push_aux(tag.as_bytes(), *i)?,
                            TagValue::Float(f) => new_record.push_aux(tag.as_bytes(), *f)?,
                            TagValue::String(s) => new_record.push_aux(tag.as_bytes(), s.as_bytes())?,
                            TagValue::Array(arr) => new_record.push_aux(tag.as_bytes(), arr)?,
                        }
                    }
                }

                output.write(&new_record)?;
            }
            None => {
                stats.reads_missing += 1;
                output.write(&buffer)?;
            }
        }
    }

    info!(
        "File stats: processed={}, modified={}, missing={}",
        stats.reads_processed, stats.reads_modified, stats.reads_missing
    );
    
    Ok(stats)
}
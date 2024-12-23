use rust_htslib::{bam, bam::record::{Record, Aux, Cigar, CigarString, CigarStringView}, bam::Read};
use rustc_hash::FxHashMap;
use anyhow::{Context, Result, anyhow};
use std::path::PathBuf;
use log::info;
use String;


#[derive(Clone)]
pub enum TagValue {
    Integer(i32),
    Float(f32),
    String(String),
    ByteArray(Vec<u8>), // For B:c tags
}

impl TagValue {
    pub fn from_aux(aux: &Aux) -> Result<Self> {
        match aux {
            Aux::I8(v) => Ok(TagValue::Integer(*v as i32)),
            Aux::U8(v) => Ok(TagValue::Integer(*v as i32)),
            Aux::I16(v) => Ok(TagValue::Integer(*v as i32)),
            Aux::U16(v) => Ok(TagValue::Integer(*v as i32)),
            Aux::I32(v) => Ok(TagValue::Integer(*v)),
            Aux::U32(v) => Ok(TagValue::Integer(*v as i32)),
            Aux::Float(v) => Ok(TagValue::Float(*v)),
            Aux::String(v) => Ok(TagValue::String(v.to_string())),
            Aux::ArrayU8(v) => {
                // Collect iterator directly into Vec<u8>
                Ok(TagValue::ByteArray(v.iter().collect()))
            },
            _ => Err(anyhow!("Unsupported tag type")),
        }
    }

    pub fn to_aux(&self) -> Aux<'_> {
        match self {
            TagValue::Integer(v) => Aux::I32(*v),
            TagValue::Float(v) => Aux::Float(*v),
            TagValue::String(v) => Aux::String(v),
            TagValue::ByteArray(v) => {
                // Convert Vec<u8> to AuxArray using into()
                Aux::ArrayU8(v.as_slice().into())
            },
        }
    }
}

#[derive(Clone)]
pub struct UnalignedRead {
    pub sequence: Vec<u8>,
    pub qualities: Vec<u8>,
    pub tags: Vec<(Vec<u8>, TagValue)>,
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
        
        // Convert tags to owned values
        let tags: Vec<_> = buffer
            .aux_iter()
            .filter_map(|r| r.ok())
            .filter_map(|(tag, aux)| {
                TagValue::from_aux(&aux)
                    .ok()
                    .map(|value| (tag.to_vec(), value))
            })
            .collect();

        index.insert(name, UnalignedRead {
            sequence: buffer.seq().as_bytes().to_vec(),
            qualities: buffer.qual().to_vec(),
            tags,
        });
    }

    info!("Indexed {} reads", index.len());
    Ok(index)
}

/// Convert CIGAR string from hard clips to soft clips
fn convert_cigar(cigar: CigarStringView) -> CigarString {
    let cigar_vec: Vec<Cigar> = cigar.iter().map(|&op| {
            match op {
                Cigar::HardClip(op_len) => Cigar::SoftClip(op_len), // Convert HardClip to SoftClip
                _ => op,
        }
    }).collect();
    CigarString(cigar_vec)
}

/// Process a single BAM file
pub fn process_bam_file(
    input_path: &PathBuf,
    unaligned_index: &FxHashMap<Vec<u8>, UnalignedRead>,
    output_dir: &PathBuf,
    transfer_tags: &[String],
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

        match unaligned_index.get(&name) {
            Some(unaligned) => {
                // Copy basic fields
                let mut new_record = buffer.clone();
                

                // Convert CIGAR if needed and prepare record data
                let cigar = convert_cigar(buffer.cigar());
                new_record.set(
                    &name,
                    Some(&bam::record::CigarString(cigar.to_vec())),
                    &unaligned.sequence,
                    &unaligned.qualities
                );

                // Add new tags from unaligned read
                for tag_name in transfer_tags {
                    let tag_bytes = tag_name.as_bytes();
                    if tag_bytes.len() == 2 {
                        // Only transfer if aligned read doesn't have the tag
                        if !buffer.aux(tag_bytes).is_ok() {
                            // Find matching tag in unaligned read
                            if let Some((_, value)) = unaligned.tags.iter()
                                .find(|(t, _)| t == tag_bytes) {
                                new_record.push_aux(tag_bytes, value.to_aux())?;
                            }
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
        "File stats: processed={}, missing={}",
        stats.reads_processed, stats.reads_missing
    );
    
    Ok(stats)
}
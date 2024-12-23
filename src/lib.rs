use rust_htslib::{bam, bam::record::{Record, Aux}, bam::Read};
use rustc_hash::FxHashMap;
use anyhow::{Context, Result};
use std::path::PathBuf;
use log::info;

#[derive(Clone)]
pub struct UnalignedRead {
    pub sequence: Vec<u8>,
    pub qualities: Vec<u8>,
    pub tags: Vec<(Vec<u8>, Aux)>,  // Store the Aux directly
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
        
        // Store tags directly
        let tags: Vec<_> = buffer
            .aux_iter()
            .filter_map(Result::ok)
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
fn convert_cigar(cigar: &[u32]) -> Vec<u32> {
    cigar.iter()
        .map(|&op| {
            let op_type = op >> 4;
            let op_len = op & 0xf;
            // If hard clip (5), convert to soft clip (4)
            if op_type == 5 {
                (4 << 4) | op_len
            } else {
                op
            }
        })
        .collect()
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
        let has_hard_clips = buffer
            .raw_cigar()
            .iter()
            .any(|&op| (op >> 4) == 5);  // 5 is BAM_CHARD_CLIP

        match unaligned_index.get(&name) {
            Some(unaligned) => {
                let mut new_record = Record::new();
                
                // Copy basic fields
                new_record.set_qname(&name);
                new_record.set_pos(buffer.pos());
                new_record.set_mapq(buffer.mapq());
                new_record.set_flags(buffer.flags());
                
                // Convert CIGAR if needed and prepare record data
                let cigar = if has_hard_clips {
                    stats.reads_modified += 1;
                    convert_cigar(buffer.raw_cigar())
                } else {
                    buffer.raw_cigar().to_vec()
                };

                // Create the full BAM record data
                let mut data = Vec::new();
                data.extend_from_slice(&name);
                data.push(0u8); // null terminator for name
                data.extend(cigar.iter().flat_map(|x| x.to_le_bytes()));
                data.extend(&unaligned.sequence);
                data.extend(&unaligned.qualities);

                // Set the complete record data
                new_record.set_data(&data);

                // Copy original tags
                for result in buffer.aux_iter() {
                    if let Ok((tag, value)) = result {
                        new_record.push_aux(tag, value)?;
                    }
                }

                // Add new tags from unaligned read
                for tag_name in transfer_tags {
                    let tag_bytes = tag_name.as_bytes();
                    if tag_bytes.len() == 2 {
                        // Only transfer if aligned read doesn't have the tag
                        if !buffer.aux(tag_bytes).is_ok() {
                            // Find matching tag in unaligned read
                            if let Some((_, value)) = unaligned.tags.iter()
                                .find(|(t, _)| t == tag_bytes) {
                                new_record.push_aux(tag_bytes, value.clone())?;
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
        "File stats: processed={}, modified={}, missing={}",
        stats.reads_processed, stats.reads_modified, stats.reads_missing
    );
    
    Ok(stats)
}
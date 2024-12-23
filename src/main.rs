use structopt::StructOpt;
use anyhow::Result;
use std::path::PathBuf;
use log::info;
use rayon::prelude::*;

// Import our library functionality
use bambutler::{create_read_index, process_bam_file};

#[derive(Debug, StructOpt)]
#[structopt(name = "bambutler", about = "Convert hard clips to soft clips and transfer tags")]
struct Opts {
    /// Input aligned BAM files
    #[structopt(long, parse(from_os_str))]
    aligned_bams: Vec<PathBuf>,

    /// Input unaligned BAM file
    #[structopt(long, parse(from_os_str))]
    unaligned_bam: PathBuf,

    /// Output directory
    #[structopt(long, parse(from_os_str))]
    output_dir: PathBuf,

    /// Tags to transfer (comma-separated list, e.g. "mv,ts,ns,pi")
    #[structopt(long, default_value = "", use_delimiter = true, value_delimiter = ',')]
    transfer_tags: Vec<String>,
}

fn main() -> Result<()> {
    env_logger::init();
    let opts = Opts::from_args();

    // Create output directory if it doesn't exist
    std::fs::create_dir_all(&opts.output_dir)?;

    // Create index from unaligned BAM
    let unaligned_index = create_read_index(&opts.unaligned_bam)?;

    // Process BAM files in parallel
    let results: Result<Vec<_>> = opts.aligned_bams
        .par_iter()
        .map(|path| process_bam_file(path, &unaligned_index, &opts.output_dir, &opts.transfer_tags))
        .collect();

    // Aggregate statistics
    let stats = results?;
    let total_processed: usize = stats.iter().map(|s| s.reads_processed).sum();
    let total_modified: usize = stats.iter().map(|s| s.reads_modified).sum();
    let total_missing: usize = stats.iter().map(|s| s.reads_missing).sum();

    info!("Complete! Total stats:");
    info!("  Processed: {}", total_processed);
    info!("  Modified: {}", total_modified);
    info!("  Missing: {}", total_missing);

    Ok(())
}
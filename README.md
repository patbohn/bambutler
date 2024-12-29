# Bam hard to soft clip Converter and tag transfer

This tool was built to facilitate usage of LAST alignments (which outputs hard-clipped alignments without additional tags) with remora resquiggling (requiring soft-clipped alignments). Surprisingly I couldn't find a tool that was able to do this transfer, so out of necessity I decided to create one myself. The complexity here lies in the fact that there were over 100 LAST alignment files, whose reads were randomly distributed within the 5 M read unaligned BAM output from dorado. I did not find a way around loading in the read data to be transfered into memory first, however even for these 5 million reads, this worked decently well with this implementation (likely because it stores them in their native data types). 

Build it as every other Rust tool (`cargo build --release`), and use it as described in its help text (`target/release/bambutler -h`)

```bash
bambutler 0.1.0
Convert hard clips to soft clips and transfer tags

USAGE:
    bambutler [OPTIONS] --output-dir <output-dir> --unaligned-bam <unaligned-bam>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --aligned-bams <aligned-bams>...      Input aligned BAM files
        --output-dir <output-dir>             Output directory
        --transfer-tags <transfer-tags>...    Tags to transfer (comma-separated list, e.g. "mv,ts,ns,pi") [default: ]
        --unaligned-bam <unaligned-bam>       Input unaligned BAM file
```
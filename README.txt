This is a simple program developed to filter records from a .bam file which correspond to alignments in specific regions of interest on a contig of interest.
The contigs/regions of interest are specified using a .bed file.
The program also recomputes MAPQ scores taking into account only the filtered records.
It is written using Rust. 

To build the program read_filter:

1. Clone this repository or otherwise download it.
2. Install Cargo. The best way to do this is to install Rust: 
`curl https://sh.rustup.rs -sSf | sh`
3. Navigate to the the directory on your local machine where you have cloned the repository
4. Run the command:
`cargo build --release`
5. This will create a new directory called `target`. Navigate to `target/release` and over there will be a program named read_filter
6. Run this program either directly from this direcory (or move it beforehand to your preferred location):
`read_filter [OPTIONS] --bam_file <FILE> --bed_file <FILE> --output <FILE>`
7. Further help (including a list of options available) can be found by running:
`read_filter --help`.

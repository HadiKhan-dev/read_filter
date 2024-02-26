#![allow(warnings)]

pub mod filtering;

use std::str;
use std::collections::HashMap;
use clap::{arg, command, value_parser, ArgAction, Command};


fn main() {

    let matches = command!()
    .version("1.0.0")
    .author("Hadi Khan")
    .about("A simple program to filter reads in a BAM file and update the MAPQ score")
    .arg(
        arg!(
            -f --bam_file <FILE> "BAM file to filter, REQUIRED"
        )
        .required(true))
        //.value_parser(value_parser!(PathBuf)))
    .arg(
        arg!(
            -b --bed_file <FILE> "BED file containing regions of interest, REQUIRED"
        )
        .required(true))
        //.value_parser(value_parser!(PathBuf)))
    .arg(
        arg!(
            -o --output <FILE> "Path for the output file to write, REQUIRED"
        )
        .required(true))
        //.value_parser(value_parser!(PathBuf)))
    .arg(
        arg!(
            -t --temperature <VALUE> "Temperature for recalculating MAPQ score"
        )
        .value_parser(value_parser!(f64))
        .default_value("10")
        
    )
    // .arg(
    //     arg!(
    //         -q --queue_limit <VALUE> "Maximum size of queue of batched reads waiting to be processed"
    //     )
    //     .value_parser(value_parser!(u32))
    //     .default_value("1000000")
    // )
    .arg(
        arg!(
            --htslib_threads <VALUE> "Number of threads to use for the htslib based reader/writer"
        )
        .value_parser(value_parser!(u32))
        .default_value("1")
    )
    .arg(
        arg!(
            --processing_threads <VALUE> "Number of threads to use for processing name batched reads"
        )
        .value_parser(value_parser!(u32))
        .default_value("1")
    )
    .arg(
        arg!(
            --drop_secondary "Flag to drop the secondary reads in the BAM file produced by the program"
        ))
    .arg(
        arg!(
            --drop_sequence "Flag to drop any sequence data present in the original BAM file"
        ))
    .get_matches();

    let bam_file = matches.get_one::<String>("bam_file").unwrap();
    let bed_file = matches.get_one::<String>("bed_file").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();

    let secondary_drop = matches.get_flag("drop_secondary");
    let seq_drop = matches.get_flag("drop_sequence");

    let full_command_vec: Vec<String> = std::env::args().into_iter().collect();
    let joined = full_command_vec.join(" ");
    let full_command = joined.as_str();


    let mut keep_secondary = true;
    let mut keep_sequence = true;

    if secondary_drop {
        keep_secondary = false;
    }

    if seq_drop {
        keep_sequence = false;
    }

    let temperature = *matches.get_one::<f64>("temperature").unwrap();

    let htslib_threads = *matches.get_one::<u32>("htslib_threads").unwrap();

    let processing_threads = *matches.get_one::<u32>("processing_threads").unwrap();

    //let queue_limit = *matches.get_one::<u32>("queue_limit").unwrap();

    filtering::filter_records(bam_file,bed_file,
        output_file,processing_threads,htslib_threads,
        temperature,keep_secondary,keep_sequence,full_command);

}
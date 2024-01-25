#![allow(warnings)]

use std::collections::HashMap;
use std::time::Instant;
use std::fs::File;
use std::path::PathBuf;
use std::io::{BufReader,BufRead};
use std::str;
use threadpool::ThreadPool;
use rust_htslib::tpool::ThreadPool as HtslibThreadPool;
use rust_htslib::{bam, bam::Read,bam::Record};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::Writer;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::header::HeaderRecord;

pub fn check_keep(record_ref: &Record,
    keep_locs: &HashMap<i32,Vec<[i32;2]>>, keep_secondary: bool,
    ) -> bool {

   if record_ref.is_supplementary() {
       return false;
   }

   if record_ref.is_unmapped() {
       return false;
   }

   if !keep_secondary {
       if record_ref.is_secondary() {
           return false;
       }
   }



   let record_id = record_ref.tid();
   let record_start = record_ref.reference_start() as i32;

   let regions_of_interest = keep_locs.get(&record_id);

   match regions_of_interest {
       Some(reg_int) => {

           for block in reg_int {
               if (block[0] <= record_start) && (record_start <= block[1]) {
                   return true;
               }
           }
           return false;
       }
       None => {return false}
   }

}

pub fn get_new_mapq(alignment_scores: &Vec<i64>,temp: f64,max_score: i64) -> Vec<u8> {

    let num_matches = alignment_scores.len();

    let mut new_qs: Vec<i64> = vec![-1;num_matches];
    let mut expd: Vec<f64> = vec![-1.0;num_matches];
    let mut probs: Vec<f64> = Vec::new();
    let mut prob_wrong: Vec<f64> = Vec::new();
    let mut phred_scaled: Vec<f64> = Vec::new();
    let mut scores: Vec<u8> = Vec::new();
    let mut full_sum: f64 = 0.0;

    for i in 0..num_matches {
        if alignment_scores[i] == -1 {
            continue;
        }
        expd[i] = ((alignment_scores[i] as f64)/(temp*(max_score as f64))).exp();
        full_sum += expd[i];
    }

    for i in 0..num_matches {
        if alignment_scores[i] == -1 {
            probs.push(-1.0);
            prob_wrong.push(-1.0);
            phred_scaled.push(-1.0);
            scores.push(255);
            continue;
        }
        probs.push(expd[i]/full_sum);
        prob_wrong.push(1.0-expd[i]/full_sum);
        phred_scaled.push(-10.0*((1.0-expd[i]/full_sum).log10()));

        let mut scaled_score = 1.0*(-10.0*((1.0-expd[i]/full_sum).log10()));

        if scaled_score > 60.0 {
            scaled_score = 60.0;
        }
        let scaled_u8 = scaled_score as u8;
        scores.push(scaled_u8);
    }

    //println!("Probs: {:?}",probs);
    //println!("Phred: {:?}",phred_scaled);
    //println!("Scores: {:?}",scores);
    return scores;

}

pub fn write_records(bam_writer: &mut Writer, records: Vec<Record>) -> () {
    for rec in records {
        bam_writer.write(&rec);
    }
    //bam_writer.;

}

pub fn do_analysis(same_records: &Vec<Record>,
    keep_locs: &HashMap<i32,Vec<[i32;2]>>, temperature: f64,
    keep_secondary: bool, keep_sequence: bool, 
    ) -> Vec<Record> {


    let mut new_records: Vec<Record> = Vec::new();

    let mut keep_flags: Vec<u8> = vec![0;same_records.len()];
    let mut alignment_score: Vec<i64> = vec![0;same_records.len()];

    for i in 0..same_records.len() {
        let rec = &same_records[i];
        if check_keep(rec,keep_locs,keep_secondary) {
            keep_flags[i] = 1;
        }
    }


    for i in 0..same_records.len() {
        if keep_flags[i] == 0 {
            alignment_score[i] = -1;
        } else {
            let rec = &same_records[i];

                if let Aux::U8(as_value) = rec.aux(b"AS").unwrap() {
                    alignment_score[i] = as_value as i64;
                } else if let Aux::U16(as_value) = rec.aux(b"AS").unwrap() {
                    alignment_score[i] = as_value as i64;
                } else if let Aux::U32(as_value) = rec.aux(b"AS").unwrap() {
                    alignment_score[i] = as_value as i64;
                }
        }
    }

    let read_length = same_records[0].seq_len_from_cigar(true) as i64;

    let new_qualities: Vec<u8> = get_new_mapq(&alignment_score,
        temperature,read_length);

    for i in 0..same_records.len() {
        if keep_flags[i] == 1{
            let mut cloned_record = same_records[i].clone();
            cloned_record.set_mapq(new_qualities[i]);
            let record_name = same_records[i].qname().clone();

            if !keep_sequence {
                //let mut new_rec = cloned_record.clone();
                let cig_string = cloned_record.cigar().clone().take();
                let cigar_seq = Some(&cig_string);

                cloned_record.set(record_name,cigar_seq,&[],&[]);
                //cloned_record = new_rec;
            }
            new_records.push(cloned_record);
        }
    }

    return new_records;
}

pub fn split_directions(records: Vec<Record>) -> Vec<Vec<Record>> {

    let mut forward_records = vec![];
    let mut backward_records = vec![];

    for record in records {
        if record.is_first_in_template() {
            //println!("{:?}",record.flags());
            //println!("FIRST");
            forward_records.push(record);
        } else if record.is_last_in_template() {
            //println!("{:?}",record.flags());
            //println!("SECOND");
            backward_records.push(record);
        } else {
            //println!("WEIRD");
        }
    }
    return vec![forward_records,backward_records];
}

pub fn filter_records(bam_filename: &str, bed_filename: &str, output_filename: &str,
    num_filtering_threads: u32,num_bam_threads: u32, temperature: f64,
    keep_secondary: bool, keep_sequence: bool,
full_command: &str) -> () {

    //Set up basic parameters like name of the program etc.
    let PROGRAM_NAME = "read_filter";
    let PROGRAM_NAME_LENGTH = PROGRAM_NAME.len();

    //Set up the threadpool for multiprocessing
    let htslib_pool = HtslibThreadPool::new(num_bam_threads).unwrap();

    //let filtering_pool: ThreadPool = ThreadPool::new(num_filtering_threads as usize);


    //Set up the sender and receiver channels
    //let (tx, rx) = channel();

    let start_time = Instant::now();
    
    //Read the files from disc
    let mut bed_file = File::open(bed_filename).unwrap();
    let bed_reader = BufReader::new(bed_file);

    let mut bam_reader = bam::Reader::from_path(bam_filename).unwrap();
    bam_reader.set_thread_pool(&htslib_pool);
    let header = bam::Header::from_template(bam_reader.header());

    //Create a new header which we will add info about the program we are running to
    let mut new_header = header.clone();

    let header_bts = new_header.to_bytes();
    let header_lines: Vec<&str> = str::from_utf8(&header_bts).unwrap().split("\n").collect();

    //Set up the structure of the new record we are going to add to the header
    let mut cur_count: i32 = 0;
    let mut last_id: &str = "";
    for line in header_lines {
        if &line[..3] != "@PG" {
            continue
        }
        let line_elements: Vec<&str> = line.split("\t").collect();
        
        for tag in line_elements {
            if &tag[..3] == "ID:" {
                last_id = &tag[3..];
                if last_id.len() < PROGRAM_NAME_LENGTH {
                    continue
                }
                if &last_id[..PROGRAM_NAME_LENGTH] == PROGRAM_NAME {
                    cur_count += 1;
                }
            }
        }
    }
    let mut id_string = PROGRAM_NAME.clone().to_string();
    if cur_count > 0 {
        id_string.push_str(".");
        id_string.push_str(&cur_count.to_string());
    }

    //Create the new record for the header
    let mut program_info = HeaderRecord::new(b"PG");
    program_info.push_tag(b"ID",id_string);
    program_info.push_tag(b"PN",PROGRAM_NAME);
    program_info.push_tag(b"PP",last_id);
    program_info.push_tag(b"VN","1.0.0");
    program_info.push_tag(b"CL",full_command);

    //Add the new record to our header
    new_header.push_record(&program_info);

    //Create a writer to write the new BAM file to, using our new header
    let mut output: Writer = bam::Writer::from_path(output_filename, &new_header, bam::Format::Bam).unwrap();
    output.set_thread_pool(&htslib_pool);
    //Set up hashmaps to keep track of which portions of contigs to keep and the reference ids of these contigs
    let mut reference_names: HashMap<String,i32> = HashMap::new();
    let mut keep_locations: HashMap<i32,Vec<[i32;2]>> = HashMap::new();

    //Populate the reference id hashmap
    let mut seq_counter: i32 = 0;
    for (key, records) in header.to_hashmap() {
        if key == "SQ" {
            for record in records {
                reference_names.insert(record["SN"].clone(),seq_counter);
                seq_counter += 1;
            }
        }
    }

    //Populate the contig portions hashmap
    for line in bed_reader.lines() {
        let readline = line.unwrap();
        let spl: Vec<&str> = readline.split("\t").collect();

        let haplotype = spl[0].to_string();
        let haplotype_id = reference_names.get(&haplotype).unwrap();
        let pos_start = spl[1].parse::<i32>().unwrap();
        let pos_end = spl[2].parse::<i32>().unwrap();

        if !keep_locations.contains_key(&haplotype_id) {
            keep_locations.insert(*haplotype_id,Vec::new());
        }

        keep_locations.get_mut(&haplotype_id).unwrap().push([pos_start,pos_end]);
        
    }


    //The main logic of the function:
    let mut last_record = Record::new();
    let mut cur_record = Record::new();
    let mut ct = 0;
    let mut pt_ct = 0;
    let mut same_vec:Vec<Record> = Vec::new();

    //Iterate over our records
    while let Some(r) = bam_reader.read(&mut cur_record) {
        
        //If we are at the first record then start populating our vecotr of records with the same name (i.e. for the same read)
        if ct == 0 {
            same_vec.push(cur_record.clone());
            last_record = cur_record.clone();
            ct += 1;
            continue;
        }

        let last_name = last_record.qname();
        let this_name = cur_record.qname();

        //If the current record is for the same read as the last one, then add it to the list of records for this read
        if last_name == this_name {
            same_vec.push(cur_record.clone());
            last_record = cur_record.clone();
            ct += 1;
        } else { //Otherwise we have all records for the last read

            //Analyse these records and filter them/update mapq scores etc.
            let new_recs = do_analysis(&same_vec,
                &keep_locations,temperature,keep_secondary,keep_sequence);

            let separated_recs = split_directions(new_recs);
            
            //Write the list of filtered records to our output file
            for group in separated_recs {
                write_records(&mut output,group);
            }
            
            //start a new list of records for the current read
            same_vec = vec![cur_record.clone()];

            //update the last record before looping over to the next record
            last_record = cur_record.clone();
            ct += 1;
        }

        // pt_ct += 1;
        // if pt_ct == 100000 {
        //     println!("{}",ct);
        //     pt_ct = 0;
        // }
        
        // if ct > 10000000 {
        //     break
        // }
    }

    //Once we have read all the record, repeat the process for the records corresponding to the very last read
    let new_recs = do_analysis(&same_vec,
        &keep_locations,temperature,keep_secondary,keep_sequence);
    let separated_recs = split_directions(new_recs);

    for group in separated_recs {
        write_records(&mut output,group);
    }

    let elapsed = start_time.elapsed();
    println!("Elapsed time: {:.2?}",elapsed);


}




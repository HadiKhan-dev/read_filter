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
use rust_htslib::bam::record::{Aux, Cigar, CigarString, Seq};
use rust_htslib::bam::header::HeaderRecord;

pub fn check_keep(record_ref: &Record,
    keep_locs: &HashMap<i32,Vec<[i32;2]>>,
    ) -> bool {

   if record_ref.is_supplementary() {
       return false;
   }

   if record_ref.is_unmapped() {
       return false;
   }

   let record_id = record_ref.tid();
   let record_start = record_ref.reference_start() as i32;
   let record_end = record_ref.reference_end() as i32;

   let regions_of_interest = keep_locs.get(&record_id);

   match regions_of_interest {
       Some(reg_int) => {

           for block in reg_int {
               if (block[0] <= record_start) && (record_start <= block[1]) && 
               (block[0] <= record_end) && (record_end <= block[1]) {
                   return true;
               }
           }
           return false;
       }
       None => {return false}
   }

}

pub fn get_new_mapq(alignment_scores: &Vec<i64>,temp: f64) -> Vec<u8> {

    let num_matches = alignment_scores.len();

    let mut new_qs: Vec<i64> = vec![-1;num_matches];
    let mut betas: Vec<f64> = vec![-1.0;num_matches];
    let mut probs: Vec<f64> = Vec::new();
    let mut prob_wrong: Vec<f64> = Vec::new();
    let mut phred_scaled: Vec<f64> = Vec::new();
    let mut scores: Vec<u8> = Vec::new();
    let mut log_full_sum: f64 = 0.0;


    for i in 0..num_matches {
        if alignment_scores[i] == -1 {
            continue;
        }
        betas[i] = (alignment_scores[i] as f64)/temp;

        //Update our running total of the log of the sum of e^(AS_i/T) so far
        if betas[i]-log_full_sum <= 20.0 {
            log_full_sum += (betas[i]-log_full_sum).exp().ln_1p();
        } else {
            log_full_sum = betas[i]+(log_full_sum-betas[i]).exp().ln_1p();
        }
    }

    for i in 0..num_matches {
        if alignment_scores[i] == -1 {
            probs.push(-1.0);
            prob_wrong.push(-1.0);
            phred_scaled.push(-1.0);
            scores.push(255);
            continue;
        }
        //Calculate probability of read actually coming from alignment location, probability of it not coming from there and the PHRED scaled score
        probs.push((betas[i]-log_full_sum).exp());
        prob_wrong.push(1.0-(betas[i]-log_full_sum).exp());
        phred_scaled.push(-10.0*((1.0-(betas[i]-log_full_sum).exp()).log10()));

        let mut scaled_score = 1.0*(-10.0*((1.0-(betas[i]-log_full_sum).exp()).log10()));

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

pub fn write_records(bam_writer: &mut Writer, records: Vec<Record>, keep_secondary: bool) -> () {
    for rec in records {
        if !keep_secondary {
            if rec.is_secondary() {
            continue
            }
        }
        bam_writer.write(&rec);
    }
    //bam_writer.;

}

pub fn do_analysis(same_records: &Vec<Record>,
    keep_locs: &HashMap<i32,Vec<[i32;2]>>, temperature: f64,
    keep_sequence: bool, 
    ) -> Vec<Record> {

    //If we are passed in an empty set of records to analyse, jusr return an empty vector
    if same_records.len() == 0 {
        return vec![];
    }

    let mut new_records: Vec<Record> = Vec::new();

    let mut keep_flags: Vec<u8> = vec![0;same_records.len()];
    let mut alignment_score: Vec<i64> = vec![0;same_records.len()];

    //Flag up those records that we wish to keep
    for i in 0..same_records.len() {
        let rec = &same_records[i];
        if check_keep(rec,keep_locs) {
            keep_flags[i] = 1;
        }
    }

    //Get the alignment scores for all our records
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
        temperature);
    
    let mut best_match_index = -1;
    let mut best_match_score = -1;

    //Find the new best march
    for i in 0..new_qualities.len() {
        if ((new_qualities[i] as isize) > best_match_score) & (new_qualities[i] != 255){
            best_match_index = i as isize;
            best_match_score = (new_qualities[i] as isize);
        }
    }

    //Code for transferring over the sequence data
    let mut found_seq_string = String::new();
    let mut found_sequence = "";

    //Find the sequence data for the read by reading the primary alignment
    for i in 0..same_records.len() {

        if same_records[i].seq().len() > found_seq_string.len() {
            found_seq_string = String::from_utf8(same_records[i].seq().as_bytes()).expect("");
            found_sequence = found_seq_string.as_str();
        }
    }

    for i in 0..same_records.len() {
        if keep_flags[i] == 1{
            let mut cloned_record: Record = same_records[i].clone();
            let record_name = same_records[i].qname().clone();

            let cig_string = cloned_record.cigar().clone().take();
            let mut cigar_seq = Some(&cig_string);

            //Remove sequence data
            cloned_record.set(record_name,cigar_seq,&[],&[]);

            cloned_record.set_mapq(new_qualities[i]);

            if (i as isize) != best_match_index {
                cloned_record.set_secondary();
            } else {

                //Keep track of whether the new record we are going to set as the primary was previously non-primary
                let mut changed_primary = false;

                if cloned_record.is_secondary() {
                    changed_primary = true;
                }

                cloned_record.unset_secondary();
                if keep_sequence {

                    let mut new_cigar_vec : Vec<Cigar> = Vec::new();


                    for item in cigar_seq.unwrap() {
                        let mut to_add: Cigar;
                        match item {
                            Cigar::HardClip(y) => {
                                to_add = Cigar::SoftClip(*y);
                            }
                            _ => {
                                to_add = item.clone();
                            }
                        }
                        new_cigar_vec.push(to_add);
                    }

                    //Create a reference to create the new Cigar String object
                    let new_cigar_vec_ref = &CigarString(new_cigar_vec);
                    let new_cigar_seq = Some(new_cigar_vec_ref);


                    cloned_record.set(record_name,new_cigar_seq,
                        found_sequence.as_bytes(),
                        &vec![ 255 as u8; found_sequence.len()]);
                }
            }
            new_records.push(cloned_record);
        }
    }

    return new_records;
}

pub fn split_directions(records: Vec<Record>) -> Vec<Vec<Record>> {
    /*
    Split the records for the forward read from the records for the reverse read
     */
    let mut forward_records = vec![];
    let mut backward_records = vec![];

    for record in records {

        if record.is_first_in_template() {
            forward_records.push(record);
        } else if record.is_last_in_template() {
            backward_records.push(record);
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

            //Break up the records by whether they are for the foward read in the pair or the reverse read
            let separated_recs = split_directions(same_vec);

            for group in separated_recs {
                //Analyse these records and filter them/update mapq scores etc.
                let new_recs = do_analysis(&group,
                    &keep_locations,temperature,keep_sequence);
                //Write the list of filtered records to our output file
                write_records(&mut output,new_recs,keep_secondary);
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
        
        if ct > 10000000 {
            break
        }
    }

    //Once we have read all the records, repeat the process for the records corresponding to the very last read
    let separated_recs = split_directions(same_vec);

    for group in separated_recs {
        //Analyse these records and filter them/update mapq scores etc.
        let new_recs = do_analysis(&group,
            &keep_locations,temperature,keep_sequence);
        //Write the list of filtered records to our output file
        write_records(&mut output,new_recs,keep_secondary);
    }
    let elapsed = start_time.elapsed();
    println!("Elapsed time: {:.2?}",elapsed);


}




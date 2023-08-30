use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::{collections::HashMap, env};

use anyhow::{anyhow, Result};
use clap::{command, error::ErrorKind, Parser};
use regex::Regex;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct GetOptions {
    input_sample_id: String,
    file_suffix: String,
    human_reference_genome: String,
    te_reference_genome: String,
    input_directory: String,
    output_directory: String,
    number_of_reads: Option<u32>,
    data_type: String,
    sequencing_type: String,
    length_insert_size: Option<f32>,
    l_std_insert_size: Option<f32>,
    read_len: Option<u32>,
    threads: Option<u32>,
    split: Option<u32>,
    multiple_bam: bool,
    bwa_mem: bool,
    genotype: bool,
}
impl GetOptions {
    fn normalize(&mut self) {
        if self.file_suffix.is_empty() {
            self.file_suffix = ".bam".to_string();
        }
        if !self.read_len.is_some() {
            self.read_len = Some(100);
        }
        if !self.split.is_some() {
            self.split = Some(20);
        }
        if !self.threads.is_some() {
            self.threads = Some(1);
        }

        if self.input_directory.is_empty() {
            self.input_directory = env::current_dir()
                .unwrap()
                .into_os_string()
                .into_string()
                .unwrap(); //perl add /
        }
        if self.output_directory.is_empty() {
            self.output_directory = env::current_dir()
                .unwrap()
                .into_os_string()
                .into_string()
                .unwrap(); //perl add /
        }
        if !self.number_of_reads.is_some() {
            self.number_of_reads = Some(3);
        }
        if self.data_type.is_empty() {
            self.data_type = "WGS".to_string();
        }
        if self.sequencing_type.is_empty() {
            self.sequencing_type = String::from("paired-end");
        }
        if self.sequencing_type.eq_ignore_ascii_case("single-end") {
            self.length_insert_size = Some(500f32);
        }
        // if self.length_insert_size.is_some() && self.l_std_insert_size.is_some(){
        //   self.std
        // }

        // GetOptions { input_sample_id: self.input_sample_id, file_suffix: self.file_suffix, human_reference_genome: self.human_reference_genome, te_reference_genome: self.te_reference_genome, input_directory: self.input_directory, output_directory: self.output_directory, number_of_reads: self.number_of_reads, data_type: self.data_type, sequencing_type: self.sequencing_type, length_insert_size: self.length_insert_size, l_std_insert_size: self.l_std_insert_size, read_len: self.read_len, threads: self.threads, split: self.split, multiple_bam: self.multiple_bam, bwa_mem: self.bwa_mem, genotype: self.genotype }
    }
}
fn main() -> Result<()> {
    let mut std_insert_size: f32 = 0.0;
    // my $bowtie2_d="";
    // my $tophat_d="";
    // my $bwa_d="";
    // my $samtools_d="";
    // my $SE_MEI_d="";
    let bowtie2_d = "";
    let tophat_d = "";
    let bwa_d = "";
    let samtools_d = "";
    let se_mei_d = "";

    let mut directory = String::new();
    let mut args = GetOptions::parse();
    args.normalize();
    if args.length_insert_size.is_some() && args.l_std_insert_size.is_some() {
        std_insert_size = args.l_std_insert_size.unwrap();
    }

    println!();
    // Step 1

    println!(
        "########  #######   ##       ##                        ##     ##\n\
         ########  ##    ##  ##       ##                        ##     ##\n\
         ##        ##    ##   ##     ##                         ##     ##\n\
         ##        ##   ##    ##     ##                         ##     ##\n\
         ########  ######      ##   ##      #####     ####      ##     ##\n\
         ########  #####       ##   ##     ##   ##   ##  ##     ##     ##\n\
         ##        ## ##        ## ##     ##        ##    ##    ##     ##\n\
         ##        ##  ##       ## ##     ##        ##    ##    ##     ##\n\
         ########  ##   ##       ###       ##   ##   ##   ##    ## ##  ## ##\n\
         ########  ##    ##      ###        #####     #### ##    ##     ##"
    );

    println!("\n\n# ERVcaller\n");
    println!(
        "# Please contact Xun Chen Ph.D. for questions and help:\n# Email: xunchen85@gmail.com or Xun.Chen@uvm.edu\n\n"
    );

    println!("\nStep 1: Loading...\n=====================================");

    if args.input_sample_id.is_empty() {
        eprintln!(
            "{} Error \n# No samples are provided ",
            ErrorKind::InvalidValue
        );
        std::process::exit(1);
    }

    let argss: Vec<String> = env::args().collect();

    if let Some(program_name) = argss.get(0) {
        directory.push_str(program_name);
        println!("{}", directory);
    }

    let mut components: Vec<&str> = directory.split('/').collect();

    if components.len() > 1 {
        components.pop(); // Remove the last component
        directory = components.join("/");
        directory.push('/');
    } else {
        directory = String::new();
    }

    let mut tsd_min_len = 100;
    let mut alignment_score = 30;
    let mut human_genome = &args.human_reference_genome;
    let mut human_genome_tophat = &args.human_reference_genome;

    let mut genome: HashMap<String, String> = HashMap::new();
    let mut order = 0;
    let min_insertsize = 0;
    let thread_1 = args.threads.unwrap() - 1;
    let double_length_insertsize = 0.0f64; // Define the correct type
    let cmd = String::new(); // You can use String for dynamic strings
    let mut bp1_tmp3: Vec<String> = Vec::new(); // Define the correct type
    let header3 = "Sample_ID Is_Split_mode ... Genotype\n".to_string();

    /////////////////////////////////
    /// ////////////////////////////
    if !PathBuf::from(&args.output_directory).is_dir() {
        fs::create_dir(&args.output_directory).expect("Failed to create directory");
    }

    //line 213
    if let Err(err) = env::set_current_dir(&args.output_directory) {
        eprintln!("Failed to change directory: {}", err);
    }

    //////// 2.1 Check input file
    let temp_directory = format!("{}_temp", &args.input_sample_id);

    if !Path::new(&temp_directory).exists() {
        if let Err(err) = fs::create_dir(&temp_directory) {
            eprintln!("Failed to create directory: {}", err);
        }
    }
    println!("\nStep 2: Detecting TE insertions...\n=====================================\n");
    let input_file_1 = format!(
        "{}{}_1.{}",
        &args.input_directory, &args.input_sample_id, &args.file_suffix
    );
    let input_file_2 = format!(
        "{}{}_2.{}",
        &args.input_directory, &args.input_sample_id, &args.file_suffix
    );
    let sequencing_type = "paired-end"; // Replace with the sequencing type

    let file_suffix_pattern = Regex::new(r"fq|fastq").unwrap(); // Regex pattern for "fq" or "fastq"
    let bam_suffix_pattern = Regex::new(r"bam|sam").unwrap(); // Regex pattern for "fq" or "fastq"

    let input_file = format!(
        "{}{}.{}",
        &args.input_directory, &args.input_sample_id, &args.file_suffix
    );

    if Path::new(&input_file_1).exists()
        && Path::new(&input_file_2).exists()
        && sequencing_type == "paired-end"
        && file_suffix_pattern.is_match(&args.file_suffix)
    {
        println!("~~~~~ paired-end reads in fastq format were loaded");
    } else if Path::new(&input_file).exists()
        && sequencing_type == "single-end"
        && file_suffix_pattern.is_match(&args.file_suffix)
    {
        println!("~~~~~ single-end read in fastq format was loaded");
    } else if Path::new(&input_file).exists()
        && sequencing_type == "paired-end"
        && bam_suffix_pattern.is_match(&args.file_suffix)
    {
        println!("~~~~~ paired-end reads in bam format were loaded\n");
        if args.genotype {
            let bai_file = format!(
                "{}{}.{}.bai",
                &args.input_directory, &args.input_sample_id, &args.file_suffix
            );

            if Path::new(&bai_file).exists() {
                println!("~~~~~ the input bam file was indexed");
            } else {
                println!("~~~~~ the input bam file was not indexed, please index the bam file using samtools for performing the validation or genotyping function");
                std::process::exit(1);
            }
        }
    } else if Path::new(&input_file).exists()
        && sequencing_type == "single-end"
        && bam_suffix_pattern.is_match(&args.file_suffix)
    {
        println!("~~~~~ single-end reads in bam format were loaded");
        if args.genotype {
            let bai_file = format!(
                "{}{}.{}.bai",
                &args.input_directory, &args.input_sample_id, &args.file_suffix
            );

            if Path::new(&bai_file).exists() {
                println!("~~~~~ the input bam file was indexed");
            } else {
                println!("~~~~~ the input bam file was not indexed, please index the bam file using samtools for performing the validation or genotyping function");
                std::process::exit(1);
            }
        }
    } else if Path::new(&input_file).exists() && args.multiple_bam {
        println!("~~~~~ a list of multiple BAM files were loaded");
    } else {
        eprintln!("# Error: could not find the input data under the provided sampleID");
        println!(
            "Input: {}{}{}",
            &args.input_directory, &args.input_sample_id, &args.file_suffix
        );
        std::process::exit(1);
    }

    ////// Step 2.2 Extract supporting reads
    if args.file_suffix.eq_ignore_ascii_case("bam") || args.multiple_bam {
        let order = 1;
        convert_bamtofastq(&args.input_sample_id);
        if !args.bwa_mem {
            align_to_hg(&format!("{}_h1", &args.input_sample_id), ".1fq");
            let order = 2;
            convert_bamtofastq(&format!("{}_h1", &args.input_sample_id));
            if std::path::Path::new(&format!("{}_h1_sm.bam", &args.input_sample_id)).exists() {
                std::fs::rename(
                    &format!("{}_h1_sm.bam", &args.input_sample_id),
                    &format!("{}_sm.bam", &args.input_sample_id),
                )
                .expect("Failed to rename file");
            }
            if std::path::Path::new(&format!("{}_h1_su.bam", &args.input_sample_id)).exists() {
                std::fs::rename(
                    &format!("{}_h1_su.bam", &args.input_sample_id),
                    &format!("{}_su.bam", &args.input_sample_id),
                )
                .expect("Failed to rename file");
            }
        }

        //
        if args.split.is_some() || &args.sequencing_type == "single-end" {
            Command::new("gunzip")
                .arg("-c")
                .arg(format!(
                    "{}_soft.fastq.gz >{}_1sf.fastq",
                    &args.input_sample_id, &args.input_sample_id
                ))
                .output()
                .expect("Failed to execute command");

            // Capture the output of the first gunzip and write it to output_1sf

            if !args.bwa_mem {
                Command::new("gunzip")
                    .arg("-c")
                    .arg(format!(
                        "{}_h1_soft.fastq.gz >>{}_1sf.fastq",
                        &args.input_sample_id, &args.input_sample_id
                    ))
                    .output()
                    .expect("Failed to execute command");

                // let mut output_1sf_append = std::fs::OpenOptions::new()
                //     .append(true)
                //     .open(format!("{}_1sf.fastq", input_sample_id))
                //     .expect("Failed to open file for appending");

                // Capture the output of the second gunzip and append it to output_1sf_append
            }
        }

        if &args.sequencing_type == "paired-end" {
            println!("paird-end sequence type");

            if !args.bwa_mem {
                // system("mv ${input_sampleID}_h1_h1_1.1fq ${input_sampleID}_1.1fq");
                // system("mv ${input_sampleID}_h1_h1_2.1fq ${input_sampleID}_2.1fq");
                Command::new("mv")
                    .arg(format!("{}_h1_h1_1.1fq", &args.input_sample_id))
                    .arg(format!("{}_1.1fq", &args.input_sample_id))
                    .output()
                    .expect("Failed to execute command");

                Command::new("mv")
                    .arg(format!("{}_h1_h1_2.1fq", &args.input_sample_id))
                    .arg(format!("{}_2.1fq", &args.input_sample_id))
                    .output()
                    .expect("Failed to execute command");
            } else {
                Command::new("mv")
                    .args(&[
                        format!("{}_h1_1.1fq", &args.input_sample_id),
                        format!("{}_1.1fq", &args.input_sample_id),
                    ])
                    .status()
                    .expect("Failed to execute command");
                Command::new("mv")
                    .args(&[
                        format!("{}_h1_2.1fq", &args.input_sample_id),
                        format!("{}_2.1fq", &args.input_sample_id),
                    ])
                    .status()
                    .expect("Failed to execute command");
            }
        } else {
            println!("Not paird-end sequence type");
            if !args.bwa_mem {
                Command::new("mv")
                    .args(&[
                        format!("{}_h1_h1.1fq", &args.input_sample_id),
                        format!("{}.1fq", &args.input_sample_id),
                    ])
                    .status()
                    .expect("Failed to execute command");
            } else {
                //bwa_MEM
                Command::new("mv")
                    .args(&[
                        format!("{}_h1.1fq", &args.input_sample_id),
                        format!("{}.1fq", &args.input_sample_id),
                    ])
                    .status()
                    .expect("Failed to execute command");
            }
        }
    } else {
        order = 0;
        align_to_hg(&args.input_sample_id, &args.file_suffix);
        convert_bamtofastq(&args.input_sample_id);
        if args.split.is_some() || &args.sequencing_type == "single-end" {
            Command::new("gunzip")
                .arg("-c")
                .arg(format!(
                    "{}_soft.fastq.gz >{}_1sf.fastq",
                    &args.input_sample_id, &args.input_sample_id
                ))
                .output()
                .expect("Failed to execute command");
        }
        if &args.sequencing_type == "paired-end" {
            Command::new("mv")
                .args(&[
                    format!("{}_h1_1.1fq", &args.input_sample_id),
                    format!("{}_1.1fq", &args.input_sample_id),
                ])
                .status()
                .expect("Failed to execute command");

            Command::new("mv")
                .args(&[
                    format!("{}_h1_2.1fq", &args.input_sample_id),
                    format!("{}_2.1fq", &args.input_sample_id),
                ])
                .status()
                .expect("Failed to execute command");
        }
    }

    ////// Filter split reads
    let sf1_file =
        File::open(format!("{}_1sf.fastq", &args.input_sample_id)).expect("Failed to open file");
    let sf2_file = File::create(format!("{}_1sf.fastq2", &args.input_sample_id))
        .expect("Failed to create file for appending");

    let mut reader = BufReader::new(sf1_file);
    let mut writer = BufWriter::new(sf2_file);

    let mut tmp1 = String::new();
    while let Ok(bytes_read) = reader.read_line(&mut tmp1) {
        if bytes_read == 0 {
            break;
        }

        if tmp1.starts_with("@soft") {
            let tmp1_parts: Vec<&str> = tmp1.split('|').collect();
            if let Some(tmp1_third_part) = tmp1_parts.get(2) {
                let tmp1_third_part_value =
                    tmp1_third_part.trim().parse::<i32>().unwrap_or_default();
                if tmp1_third_part_value % 4 >= 2 {
                    // Handle single-end case if needed
                    writer.write_all(tmp1.as_bytes())?;
                    reader.read_line(&mut tmp1)?;
                    writer.write_all(tmp1.as_bytes())?;
                    reader.read_line(&mut tmp1)?;
                    writer.write_all(tmp1.as_bytes())?;
                    reader.read_line(&mut tmp1)?;
                    writer.write_all(tmp1.as_bytes())?;
                } else {
                    reader.read_line(&mut tmp1)?;
                    reader.read_line(&mut tmp1)?;
                    reader.read_line(&mut tmp1)?;
                }
            }
        }

        tmp1.clear();
    }

    Command::new("mv")
        .args(&[
            format!("{}_1sf.fastq2", &args.input_sample_id),
            format!("{}_1sf.fastq", &args.input_sample_id),
        ])
        .status()?;

    // 2.3 Chimeric reads amd Split reads
    println!("\nChimeric and split reads...\n=====================================\n");
    if &args.sequencing_type == "paired-end" {
        Command::new(format!("{}bwa",&bwa_d))
            .arg(format!("mem -t {} -k 19 -r 1.5 -c 100000 -m 50 -T 30 -h 10000 -a -Y -M $TE_reference_genomes {}_1.1fq {}_2.1fq >{}_vsu.sam", &args.threads.unwrap(), &args.input_sample_id, &args.input_sample_id, &args.input_sample_id))
            .output()?;
    } else {
        Command::new("touch")
            .arg(format!("{}_vsu.sam", &args.input_sample_id))
            .status()?;
        //File::create(format!("{}_vsu.sam", &args.input_sample_id))?;
    }
    Command::new("touch")
            .arg(format!("{}_all_breakpoint", &args.input_sample_id))
            .status()?;

    









    Ok(())
}

fn align_to_hg(format: &str, arg: &str) {
    todo!()
}

fn convert_bamtofastq(input_sample_id: &str) {
    todo!()
}

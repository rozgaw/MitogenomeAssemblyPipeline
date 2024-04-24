#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 10:12:04 2024

@author: haleyatkins
"""

## fish mitogenome assembly 
import os
import argparse
import glob

#for SRA files 
def convert_sra_to_fastq(main):
    untrimmed_fastq_dir = os.path.join(main, "untrimmed_fastqs")
    os.makedirs(untrimmed_fastq_dir, exist_ok=True)
    os.chdir(untrimmed_fastq_dir)

    fastq_dump = "fastq-dump -I --split-files"
    for sra in glob.glob(os.path.join(main, 'SRAfiles', '*')):
        fastqdump_cmd = f"{fastq_dump} {sra}"
        os.system(fastqdump_cmd)

def get_files(main):
    os.makedirs(os.path.join(main, 'SRAfiles'), exist_ok=True)
    os.chdir(os.path.join(main, 'SRAfiles'))
    file_url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR18689888/SRR18689888"
    os.system(f"wget {file_url}")

def fastQC(main):
    untrimmed_fastq_dir = os.path.join(main, "untrimmed_fastqs")
    qc_dir = os.path.join(main, "fastqc")
    os.makedirs(qc_dir, exist_ok=True)
    os.chdir(qc_dir)

    # Check if there are any .fastq files present in untrimmed_fastqs
    fastq_files = glob.glob(os.path.join(untrimmed_fastq_dir, "*.fastq"))
    if not fastq_files:
        print("No .fastq files found in the untrimmed_fastqs directory. Skipping FastQC.")
        return

    fastqc_command = f"fastqc {os.path.join(untrimmed_fastq_dir, '*.fastq')}"
    os.system(fastqc_command)

#fastp trimming 
def fastp(main):
    
    untrimmed_fastq_dir = os.path.join(main, "untrimmed_fastqs")
    trimmed_fastq_dir = os.path.join(main, "trimmed_fastqs")
    os.makedirs(trimmed_fastq_dir, exist_ok=True)

    #list all untrimmed fastq files
    untrimmed_fastq_files = glob.glob(os.path.join(untrimmed_fastq_dir, "*.fastq"))

    #iterate over untrimmed fastq files and generate corresponding trimmed filenames
    for untrimmed_fastq in untrimmed_fastq_files:
        # Paired-end reads
        trimmed_fastq_base = os.path.basename(untrimmed_fastq).replace("_1.fastq", "_1_trimmed.fastq").replace("_2.fastq", "_2_trimmed.fastq")
        trimmed_fastq_1 = os.path.join(trimmed_fastq_dir, trimmed_fastq_base.replace("_1_trimmed.fastq", "trimmed_1.fastq"))
        trimmed_fastq_2 = os.path.join(trimmed_fastq_dir, trimmed_fastq_base.replace("_2_trimmed.fastq", "trimmed_2.fastq"))
        #command for paired-end reads with adaptor trimming only, no quality trimming 
        fastp_command = f"fastp -i {untrimmed_fastq} -I {untrimmed_fastq.replace('_1.fastq', '_2.fastq')} -o {trimmed_fastq_1} -O {trimmed_fastq_2} -Q"
        # Execute fastp command
        os.system(fastp_command)
        # Execute fastp command
        os.system(fastp_command)
 
def GetOrganelle(main):
    os.chdir(os.path.join(main, "trimmed_fastqs"))
    getorg_dir = os.path.join(main, "GetOrganelle_assembly")
    os.makedirs(getorg_dir, exist_ok=True)
    os.chdir(getorg_dir)
    go_command = "get_organelle_from_reads.py -s cytb.fasta -1 *_1_trimmed.fastq -2 *_2_trimmed.fastq -R 10 -k 21,45,65,85,105 -F animal_mt -o mito_fish_assembly"
    os.system(go_command)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", help="specify the main directory", default=os.getcwd())
    args = parser.parse_args()

    # using the cwd as the parent directory of the pipeline output folder
    main_directory = args.directory
    os.makedirs(main_directory, exist_ok=True)
    get_files(main_directory)
    convert_sra_to_fastq(main_directory)
    fastQC(main_directory)
    fastp(main_directory)
    GetOrganelle(main_directory)

if __name__ == "__main__":
    main()
    



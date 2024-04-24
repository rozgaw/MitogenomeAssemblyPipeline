#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 20:57:19 2024

@author: haleyatkins
"""

import os
import argparse
import glob

#convert sra to fastq files 
def convert_sra_to_fastq(main):
    #create directory for the untrimmed fastq files 
    untrimmed_fastq_dir = os.path.join(main, "untrimmed_fastqs")
    os.makedirs(untrimmed_fastq_dir, exist_ok=True)
    #change into the untrimmed fastq directory just made
    os.chdir(untrimmed_fastq_dir)
    #for SRA files only! converts SRA into fastq files 
    fastq_dump = "fastq-dump -I --split-files"
    
    #create the path for the SRA files
    sra_file = os.path.join(main, 'SRAfiles', 'your_sra_file.sra')
    #calls the fastq dump with specified SRA URL
    fastqdump_cmd = f"{fastq_dump} {sra_file}"
    #executes fastq command
    os.system(fastqdump_cmd)


def get_files(main, file1=None, file2=None):
    #create the directory to store the SRA files
    sra_dir = os.path.join(main, 'SRAfiles')
    os.makedirs(sra_dir, exist_ok=True)
    #change into the SRA directory
    os.chdir(sra_dir)
    #if both FASTQ files are given copy them into SRA directory
    if file1 and file2:
        os.system(f"cp {file1} {file2} .")
    else:
        #prompt the user for the SRA URL or the paths to the FASTQ files
        file_input = input("Enter the SRA URL or provide the paths to the two FASTQ files: ")
        if file_input.startswith("https://"):
            #if SRA URL is given download the SRA file with wget
            os.system(f"wget {file_input}")
        else:
            #if the input is not an URL assume it's the paths to the fastq files
            #splits input into a list of strings for fastq1 and fastq2
            file1, file2 = file_input.split()
            #executes the command to copy the files to the SRA directory
            os.system(f"cp {file1} {file2} .")
            
#fastQC for quality check
def fastQC(main):
    #make untrimmed fastqs directory for the fastqs and fastQC outputs
    untrimmed_fastq_dir = os.path.join(main, "untrimmed_fastqs")
    qc_dir = os.path.join(main, "fastqc")
    os.makedirs(qc_dir, exist_ok=True)
    #change into the fastqc directory 
    os.chdir(qc_dir)

    #fastq_files = glob.glob(os.path.join(untrimmed_fastq_dir, "*.fastq"))
    #if not fastq_files:
        #print("No .fastq files found in the untrimmed_fastqs directory. Skipping FastQC.")
        #return
    #fastqc command for the fastqs in untrimmed directory 
    fastqc_command = f"fastqc {os.path.join(untrimmed_fastq_dir, '*.fastq')}"
    # executes the fastqc command 
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
 
#getorganelle assembly 
def GetOrganelle(main):
    #change into the trimmed fastq directory 
    os.chdir(os.path.join(main, "trimmed_fastqs"))
    #create the assembly directory
    getorg_dir = os.path.join(main, "GetOrganelle_assembly")
    os.makedirs(getorg_dir, exist_ok=True)
    #go into the getorganelle assembly directory 
    os.chdir(getorg_dir)
    #getorganelle command with the cytb.fasta file 
    go_command = "get_organelle_from_reads.py -s cytb.fasta -1 *_1_trimmed.fastq -2 *_2_trimmed.fastq -R 10 -k 21,45,65,85,105 -F animal_mt -o mito_fish_assembly"
    #executes the getorganelle command 
    os.system(go_command)

def main():
    parser = argparse.ArgumentParser()
    #allows the user to specify the main directory
    parser.add_argument("--directory", help="specify the main directory", default=os.getcwd())
    #argument to give path for fastq1
    parser.add_argument("--file1", help="First FASTQ file")
    #argument to give path for fastq2
    parser.add_argument("--file2", help="Second FASTQ file")
    args = parser.parse_args()

    #using cwd as the parent directory 
    main_directory = args.directory
    
    #check if both fastq files are provided directly
    if args.file1 and args.file2:
        file1 = args.file1
        file2 = args.file2
    else:
        #set file paths to None if not given
        file1 = None
        file2 = None

    #create main directory
    os.makedirs(main_directory, exist_ok=True)
    #call get_files function to handle SRA URL or provided FASTQ files
    get_files(main_directory, file1, file2)
    
    #call the normal pipeline functions
    fastQC(main_directory)
    fastp(main_directory)
    GetOrganelle(main_directory)

if __name__ == "__main__":
    main()
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

def convert_sra_to_fastq(main):
    os.makedirs(os.path.join(main, "untrimmed_fastqs"), exist_ok=True)
    os.chdir(os.path.join(main, "untrimmed_fastqs"))

    fastq_dump = "fastq-dump -I --split-files"
    for sra in glob.glob(os.path.join(main, 'SRAfiles', '*')):
        fastqdump_cmd = f"{fastq_dump} {sra}"
        os.system(fastqdump_cmd)

def get_files(main):
    os.makedirs(os.path.join(main, 'SRAfiles'), exist_ok=True)
    os.chdir(os.path.join(main, 'SRAfiles'))
    fastq_files = ["https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR18689888/SRR18689888"]
    for file_url in fastq_files:
        os.system(f"wget {file_url}")

def fastQC(main):
    os.chdir(os.path.join(main, "untrimmed_fastqs"))
    fastqc_command = "fastqc *fastq"

def GetOrganelle(main):
    os.chdir(os.path.join(main, "untrimmed_fastqs"))
    go_command = "get_organelle_from_reads.py -1 *1.fastq -2 *2.fastq -R 10 -k 21,45,65,85,105 -F animal_mt -o animal_mt_out"

def main():
    parser = argparse.ArgumentParser(description="Automated pipeline for processing SRA files.")
    parser.add_argument("directory", type=str, help="Directory where the files are located.")
    args = parser.parse_args()

    main_directory = args.directory
    get_files(main_directory)
    convert_sra_to_fastq(main_directory)
    fastQC(main_directory)
    GetOrganelle(main_directory)

if __name__ == "__main__":
    main()
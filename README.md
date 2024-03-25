# MitogenomeAssemblyPipeline

1. Download SRA FASTQ files
     * Make sure you have the SRA toolkit installed
     * In the terminal run: "prefetch [SRA ID]
          * In our case, we ran: "prefetch SRR18689888"
     * Next to split the paired end reads run: "fastq-dump -I --split-files SRR18689888"
          * This should leave you with two fastq files. One "[SRR ID]_1.fastq" for the forward read and the other "[SRR ID]_2.fastq" for the reverse read.
          * In our case the files are "SRR18689888_1.fastq" and "SRR18689888_2.fastq"

3. Evaluate read quality using __FastQC__
     * Input: Raw sequencing data: FASTQ files
     * Output: Summary graphs and tables to help assess read quality
4. Assemble the mitochondrial genome using __SMART2__
     * Inputs:
          * Sample Name
          * Seed gene file, FASTA format
          * Illumina DNA-Seq reads, paired-end FASTQ format
     * Outputs:
          * Consensus sequences FASTA file (assembled mitogenome)
          * Annotation
          * Compressed tar archive: images and plots of assemblies graphs, clusters, and annotation
          * Report file
5. Run alignment of raw FASTQ files with the assembled mitogenome from SMART2 (FASTA) using __Bowtie2__
     * Raw illumina reads first trimmed using __fastp__ to remove sequencing adaptors
     * Inputs:
          * Mitogenome FASTA as reference genome
          * Raw FASTQ files to be aligned
     * Outputs: SAM file of alignments
     * Reads filtered to include only reads that aligned to the mitogenome - only those reads  will be used for the following steps

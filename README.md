# MitogenomeAssemblyPipeline

1. Download SRA FASTQ files
     * Make sure you have the SRA toolkit installed
     * In the terminal run: "prefetch [SRA ID]"
          * In our case, we ran: "prefetch SRR18689888"
     * Next to split the paired end reads run: "fastq-dump -I --split-files SRR18689888"
          * This should leave you with two fastq files. One "[SRA ID]_1.fastq" for the forward read and the other "[SRA ID]_2.fastq" for the reverse read.
          * In our case the files are "SRR18689888_1.fastq" and "SRR18689888_2.fastq"

2. Evaluate read quality using __FastQC__
     * Input: Raw sequencing data: FASTQ files
     * Command: fastqc SRR18689888_1.fastq SRR18689888_2.fastq
     * Output: Summary graphs and tables to help assess read quality
     * To look at the visualizations in FastQC:
     	* open SRR18689888_1_fastqc.html
     	* open SRR18689888_2_fastqc.html
     	* opens the html files in your browser to view the visualizations
3. Assemble the mitochondrial genome using __SMART2__
    * Inputs:
          * Sample Name: "Name the assembly however you want" 
          * Seed gene file, FASTA format: "sequence.fasta"
		* This is highly conserved sequence from across different fish species
		* In our case, we're using a cytochrome b partial cds sequence of the species we're assembling
          * Illumina DNA-Seq reads, paired-end FASTQ format: "SRR18689999_1.fastq" and "SRR18689888_2.fastq"
     * Outputs:
          * Consensus sequences FASTA file (assembled mitogenome): "scaffold_seqs.fasta"
          * Annotation file: Result file 
          * Compressed tar archive: images and plots of assemblies graphs, clusters, and annotation
          * Report file
          * Download the mitogenome assembly-Results folder.zip to access output files
4. Run alignment of raw FASTQ files with the assembled mitogenome from SMART2 (FASTA) using __Bowtie2__
     * Raw illumina reads first trimmed using __fastp__ to remove sequencing adaptors
     * Command: "fastp -i SRR18689888_1.fastq -I SRR18689888_2.fastq -o out.SRR18689888_1.fastq -O out.SRR18689888_2.fastq"
     * Create the Bowtie2 index using the assembled mitogenome from SMART2:
          * Command: "bowtie2-build scaffold_seqs.fasta mitogenome_index"
     * Align the Trimmed FASTQ Files to the Mitogenome:
          * Command: "bowtie2 -x mitogenome_index -1 out.SRR18689888_1.fastq -2 out.SRR18689888_2.fastq -S alignment.sam --al-conc SRR18689888_mapped_%.fq"
     * Inputs:
          * Mitogenome FASTA (scaffold_seqs.fasta) as reference genome
          * Trimmed FASTQ files to be aligned
     * Outputs: SAM file of alignments (Can convert to BAM file for efficiency - samtools view -b -F 4 alignment.sam > aligned_reads.bam)
     	* Output file: alignment.sam (OR aligned_reads.bam)
     * Output file will contain sequencing reads that are aligned to the mitogenome - only those reads  will be used for the following steps


NEW PLAN: TRYING OUT GETORGANELLE
      * first must download GetOrganelle
           * Use the following page for download instructions: https://github.com/Kinggerm/GetOrganelle/wiki/Installation#installation
	   

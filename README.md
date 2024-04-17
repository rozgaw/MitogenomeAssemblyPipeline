# Fish Mitogenome Assembly Pipeline

## Dependencies & Packages 
* SRA toolkit <https://github.com/ncbi/sra-tools> 




1. Download SRA FASTQ files
     * Make sure you have the SRA toolkit installed
     * In the terminal run: "prefetch [SRA ID]"
          * In our case, we ran: "prefetch SRR18689888"
     * Next to split the paired end reads run: "fastq-dump -I --split-files SRR18689888"
          * This should leave you with two fastq files. One "[SRA ID]_1.fastq" for the forward read and the other "[SRA ID]_2.fastq" for the reverse read.
          * In our case the files are "SRR18689888_1.fastq" and "SRR18689888_2.fastq"

2. Evaluate read quality using __FastQC__
     * Input: Raw sequencing data: FASTQ files
     * Command (in terminal): fastqc SRR18689888_1.fastq SRR18689888_2.fastq
     * Output: Summary graphs and tables to help assess read quality
     * To look at the visualizations in FastQC:
     	* open SRR18689888_1_fastqc.html
     	* open SRR18689888_2_fastqc.html
     	* opens the html files in your browser to view the visualizations

3. Trim sequences using __fastp__
   * Input: Raw sequencing data: FASTQ files
   * Command: "fastp -i SRR18689888_1.fastq -I SRR18689888_2.fastq -o out.SRR18689888_1.fastq -O out.SRR18689888_2.fastq"
   * Command to trim adaptors only (No quality trimming): "fastp -i SRR18689888_1.fastq -I SRR18689888_2.fastq -o out.SRR18689888_1.fastq -O out.SRR18689888_2.fastq -Q"
   * Output: trimmed fastq files out.file_1.fastq out.file_2.fastq

4. Assemble the mitogenome using __GetOrganelle__
   * Use the following page for download instructions: https://github.com/Kinggerm/GetOrganelle/wiki/Installation#installation
   * Input: adaptor trimmed fastq files (out.SRR18689888_1.fastq, out.SRR18689888_2.fastq)
   * Command: "get_organelle_from_reads.py -s cytb.fasta -1 out.SRR18689888_1.fastq -2 out.SRR18689888_2.fastq -R 10 -k 21,45,65,85,105 -F animal_mt -o fish_assembly"
   * Output: 
   	    * *.path_sequence.fasta, each fasta file is an assembled genome
   	    * *.selected_graph.gfa, the organelle-only assembly graph 
   	    * get_org.log.txt, the log file 
   	    * extended_K*.assembly_graph.fastg, the raw assembly graph
   	    * extended_K*.assembly_graph.fastg.extend_embplant_pt-embplant_mt.fastg, a simplified assembly graph
   	    * extended_K*.assembly_graph.fastg.extend_embplant_pt-embplant_mt.csv, a tab-format contig label file for bandage visualization

   * Most important file is the *.fasta file and all of teh other files can deleted/ignored if the full genome is complete (you can find this information in the log file too)

5. Annotate the assembled genome using __MitoFish__ (__MitoAnnotator__)
   * Input: Mitogenome file in FASTA format
     	* Needs to be less than 100 Kb
     	* Note if DNA is circular (complete)
     	* Note if visualization is desired (slower run time)
   * Output: Download the files from MitoAnnotator, it only keeps them on the server for 10 days 
   	    * Annotation.pdf (visualization)
   	    * genes.fasta
   	    * raw.fasta
   	    * NCBI.txt
   	    * Annotation.txt
  
*** OPTIONAL: Necessary for us to show that GetOrganelle is a comparable tool for mitgenome assembly but not necessary for running the pipeline itself *** 

6. Compare GetOrganelle assembly to SMART2 assembly using __BLAST__
   * Map the assembly output to the paper's output from SMART2 using BLAST to see how they compare
   * The published assembly can be accessed using GenBank accession no. ON310810

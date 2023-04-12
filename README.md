# Bacterial-WGS
## Analysis of ONT and Illumina WGS data
## A. Illumina Pipeline
1. Data acqusition in fastq format and unzipping
2. Quality control (QC) of raw data: This step involves checking the quality of the raw sequencing data using quality control software, such as FastQC/MultiQC, to assess the overall quality, sequence length distribution, base quality, and sequence duplication.

3. Read trimming and filtering: In this step, low-quality bases and reads are removed to improve the overall quality of the sequencing data. Popular tools for trimming and filtering include Fastp, Trimmomatic, Sickle and BBDuk. Here we will use  Fastp.

4. Genome assembly: This step involves the construction of contigs from the filtered reads using assembly software, such as SPAdes, ABySS, RadTag, BWA, Bowtie or Velvet. This step generates the draft genome assembly of the bacterium. Here we will use SPAdes. de novo assembly

5. Quality assessment of genome assembly: This step involves assessing the quality of the genome assembly using tools like QUAST, which calculates assembly metrics such as contig N50, genome size, and number of contigs. Tools applied QUAST, SAM & Quali Maps.

6. Polishing refers to the process of refining and improving the quality of a genome sequence assembly. It involves correcting errors in the initial assembly, closing gaps between contigs, and improving the accuracy of base calls. Tools applied  Pilon, Quiver, and Racon.

7. Pathogen identification: This step involves identifying the bacterial species or strains present in the sequenced sample using tools such as Kraken, MetaPhlAn, or Centrifuge.

8. Multilocus sequence typing (MLST) is a typing technique of multiple loci, using DNA sequences of internal fragments of multiple housekeeping genes to characterize isolates of microbial species/to classify and identify different strains of bacteria.

9. SNPs/InDels Identification involves variants calling ie identifying single nucleotide polymorphisms (SNPs) and small insertions or deletions (indels) in the genome assembly. Tools like FreeBayes, GATK, BCF, rtg can be used for this purpose. Viewing can be done by IGV, UCSC etc

10. Genome annotation: In this step, the genome assembly is annotated with gene prediction software such as Prokka or RAST, which predict genes and annotate them with functional information ie checks non-coding RNAs, CRISPRs, pathogenic and susceptibility genes/virulence. Blasting is also done here.

11. Comparative genomics: Comparative genomics involves comparing the newly sequenced genome with existing genomes to identify unique features, virulence factors, and antibiotic resistance genes. Tools like Roary or Prokka can be used for this purpose.
        i. Evolutionary analysis use phylogenetic tree/mapping can be constructed using the identified SNPs to analyze the 
          evolutionary relationship between different bacterial strains. ANI-Dendogram which uses dREP gives genetic 
          relatedness between two bacterial genomes. It calculates the average percentage of nucleotide sequence identity 
          between the two genomes, taking into account both conserved and divergent regions.
        ii. Pangenome analysis use Roary.
        iii. Genome visualization in ring structures use BRIG.
        iv. AMR genes use RGI, CARD, PCA, ResFinder
        v. Functional annotation COG, GO, KEGG
        vi. DIvergence Time estimation
12. Metadata - linking genotypic data with phenotypic traits, Plasmids, Mutations etc

# Steps
1. Log into hpc via ssh & interactive computing
```
interactive -w compute05 -c 3
``` 
2. Creating directory
```
cd /var/scratch/
mkdir -p $USER/bacteria-wgs/crpa
cd $USER/bacteria-wgs/crpa
```
The mkdir -p command is used to create a directory and its parent directories (if they do not exist) in a single command.

The -p option allows the mkdir command to create the parent directories if they do not exist.

3. Data retieval: 
Use scp username@remote:/path/to/data /path/to/local/directory
```
scp -r woguta@hpc.ilri.cgiar.org:/path/to/data  .

```
The -r option is used to copy the directory recursively. The . at the end of the command specifies the current directory as the destination.
space fullstop " ." specifies your current folder.

4. Viewing  the files: 
 cd into current directory and view all the retrieved data by run
 ```
less filename.fastq
head -n 10 filename.fastq
tail -n 10 filename.fastq
```
5. Unzip .gz files: 
Use the gunzip or gzip command depending on your system
```
gunzip -k *.gz
```
it keeps the .gz files too, to remove/delete .gz files run within your "current directory!"
```
rm -f *.gz
```
alternatively, run this code once without k as not to keep .gz files
```
gunzip *.gz
```
View the files as in 3 above

6 Load modules, cd database
```
module load fastqc/0.11.9
module load fastp/0.22.0
module load krona/2.8.1
module load centrifuge/1.0.4
module load kraken/2.1.2
module load spades/3.15
module load quast/5.0.2
module load samtools/1.15.1
module load BUSCO/5.2.2
module load bowtie2/2.5.0
module load bedtools/2.29.0
module load bamtools/2.5.1
module load ivar/1.3.1
module load snpeff/4.1g
module load bcftools/1.13
module load nextclade/2.11.0
module load R/4.2
module load prokka/1.11
module load blast/2.12.0+
```
7. Run fastqc on one sample
```
fastqc
        -t 4
        -o ./results/fastqc/
        -f fastq ./raw_data/Fastq/AS-26335-C1-C_S4_L001_R1_001.fastq
                ./raw_data/Fastq/AS-26335-C1-C_S4_L001_R2_001.fastq
```
Here, the -t option specifies the number of CPU threads to use, -o specifies the output directory, -f specifies the file format (optional for FASTQ files), and the two positional arguments specify the paths to the two input FASTQ files.

8. Run fastp

a) for one file
```
fastp --in1 ./raw_data/Fastq/AS-26335-C1-C_S4_L001_R1_001.fastq.gz \
	--in2 ./raw_data/Fastq/AS-26335-C1-C_S4_L001_R2_001.fastq.gz \
	--out1 ./results/fastp/AS-26335-C1-C_S4_L001_R1_001.trim.fastq.gz \
	--out2 ./results/fastp/AS-26335-C1-C_S4_L001_R2_001.trim.fastq.gz \
	--json results/fastp/AS-26335-C1-C_S4_L001.fastp.json \
	--html results/fastp/AS-26335-C1-C_S4_L001.fastp.html \
	--failed_out ./results/fastp/AS-26335-C1-C_S4_L001_fail.fastq.gz \
	--thread 4 \
	-5 -3 -r \
	--detect_adapter_for_pe \
	--qualified_quality_phred 20 \
	--cut_mean_quality 20\
	--length_required 15 \
	--dedup \
	|& tee ./results/fastp/AS-26335-C1-C_S4_L001.fastp.log
        
```
--json: output a json file that summarizes the processing statistics of the reads
--html: output an html file that visualizes the processing statistics of the reads
--failed_out: output the reads that failed to pass the processing steps to a separate file
--thread: number of threads to use, here 3
-5 -3 -r: options for adapter trimming and quality filtering
--detect_adapter_for_pe: automatically detect and trim adapters for paired-end reads
--qualified_quality_phred: quality threshold for base calling
--cut_mean_quality: quality threshold for sliding window trimming
--length_required: minimum length required for reads to be kept after processing
--dedup: deduplicate identical reads
|& tee: save the command output to a log file as well as display it in the terminal

b) Rerun the fastqc on the trimmed one file
```
fastqc
        -t 4
        -o ./results/fastqc/
        -f fastq ./raw_data/Fastq/AS-26335-C1-C_S4_L001_R1_001.trim.fastq
                ./raw_data/Fastq/AS-26335-C1-C_S4_L001_R2_001.trim.fastq
```

b) Fastqc loop for all the the files

```
nano
```
Create shebang scripts
```
#!/bin/bash

#bash generate-fastqc-reports.sh ./results/fastqc/ ./raw_data/Fastq/

#Make directory to store the results
mkdir -p ./results/fastqc/

#load fastqc module
module load fastqc/0.11.4

#fastqc reports directory
REPORT_DIR=$1

#fastq files directory
FASTQ_DIR=$2

#run fastqc tool
for file in $FASTQ_DIR/*.fastq; do
        fastqc ${file} -o ${REPORT_DIR} -f fastq
        done
````
To use this script, you can save it to a file (e.g., run_fastqc.sh), make it executable (chmod +x run_fastqc1.sh), and then run it with the directory paths as arguments
```
bash run_fastqc1.sh
```

Running fastp for all the files/trimming all the files
i) First Pathway
```
#!/bin/bash

#Set the input and output directories
INPUT_DIR=./raw_data/Fastq/
OUTPUT_DIR=./results/fastp/

#load modules
module load fastp/0.22.0

# make output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

#Loop through all files in the input directory
for R1 in $INPUT_DIR/*R1_001.fastq.gz
do
    R2=${R1/R1_001.fastq.gz/R2_001.fastq.gz}
    NAME=$(basename ${R1} _R1_001.fastq.gz)
    fastp --in1 ${R1} \
          --in2 ${R2} \
          --out1 ${OUTPUT_DIR}/${NAME}.R1.trim.fastq.gz \
          --out2 ${OUTPUT_DIR}/${NAME}.R2.trim.fastq.gz \
          --json ${OUTPUT_DIR}/${NAME}.fastp.json \
          --html ${OUTPUT_DIR}/${NAME}.fastp.html \
          --failed_out ${OUTPUT_DIR}/${NAME}.failed.fastq.gz \
          --thread 4 \
          -5 -3 -r \
          --detect_adapter_for_pe \
          --qualified_quality_phred 20 \
          --cut_mean_quality 20 \
          --length_required 15 \
          --dedup \
          |& tee ${OUTPUT_DIR}/${NAME}.fastp.log
done
```
Save and run
```
 sbatch -w compute05 run_fastp1.sh
```
iii) Second pathway sbatch - sunmitiing jobs to the cluster as you do other things (the best option)
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J fastp
#SBATCH -n 4

#load modules
module load fastp/0.22.0

#making directories
INPUT_DIR=./raw_data/Fastq/
OUTPUT_DIR=./results/fastp/

# make output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

for R1 in $INPUT_DIR/*R1_001.fastq.gz
do
    R2=${R1/R1_001.fastq.gz/R2_001.fastq.gz}
    NAME=$(basename ${R1} _R1_001.fastq.gz)
    fastp --in1 ${R1} \
          --in2 ${R2} \
          --out1 ${OUTPUT_DIR}/${NAME}.R1.trim.fastq.gz \
          --out2 ${OUTPUT_DIR}/${NAME}.R2.trim.fastq.gz \
          --json ${OUTPUT_DIR}/${NAME}.fastp.json \
          --html ${OUTPUT_DIR}/${NAME}.fastp.html \
          --failed_out ${OUTPUT_DIR}/${NAME}.failed.fastq.gz \
          --thread 4 \
          -5 -3 -r \
          --detect_adapter_for_pe \
          --qualified_quality_phred 20 \
          --cut_mean_quality 20 \
          --length_required 15 \
          --dedup \
          |& tee ${OUTPUT_DIR}/${NAME}.fastp.log
done
```
run the job as saved

```
sbatch -w compute05 run_fastp2.sh
```
Checking the submitted job
```
ls -lth
```
Check the nano slurm number and squeue

```
nano slurm number

```
Then 
```
squeue
```
These will loop through all the R1 fastq files in the INPUT_DIR directory and run fastp on each pair of R1 and R2 files, outputting the trimmed files and the reports to the OUTPUT_DIR directory. The ${R1/R1_001.fastq.gz/R2_001.fastq.gz} line replaces the _R1_001.fastq.gz suffix in the filename with _R2_001.fastq.gz to get the R2 file. The ${basename ${R1} _R1_001.fastq.gz} line extracts the basename of the R1 file without the _R1_001.fastq.gz suffix to use as the name of the output files.

Do fastqc for all the trimmed files
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J fastqc
#SBATCH -n 4

# load modules
module load fastqc/0.11.9

# set input and output directories
INPUT_DIR=./results/fastp
OUTPUT_DIR=./results/fastqc

# make output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# loop over all the trimmed fastq files in the input directory
for file in $INPUT_DIR/*.trim.fastq.gz; do

    # get the filename without the extension
    filename=$(basename "$file" .trim.fastq.gz)

    # run fastqc on the file and save the output to the output directory
    fastqc -t 4 -o $OUTPUT_DIR $file

done
```
run the job as saved

```
sbatch -w compute05 run_fastqc_trim.sh
```
9. Genome Assembly using SPAdes for the trimmed files

i) for one sample
```
spades.py -k 27 \
-1 ./results/fastp/AS-27566-C1_S5_L001.R1.trim.fastq.gz \
-2 ./results/fastp/AS-27566-C1_S5_L001.R2.trim.fastq.gz \
-o ./results/spades \
-t 4 \
-m 100
```
k-mer size: 27 (-k 27)
Input read files: AS-27566-C1_S5_L001.R1.trim.fastq.gz and AS-27566-C1_S5_L001.R2.trim.fastq.gz (-1 and -2)
Output directory: ./results/spades (-o)
Number of threads: 4 (-t 4)
Memory limit: 100 GB (-m 100)

For all samples via loop
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J SPAdes
#SBATCH -n 4

# Load modules
module load spades/3.15

# Define input and output directories
INPUT_DIR=./results/fastp
OUTPUT_DIR=./results/spades

#make output directory
mkdir -p "${OUTPUT_DIR}"

# Iterate over all files in input directory
for file in ${INPUT_DIR}/*_R1.trim.fastq.gz
do
  # Extract sample name from file name
  SAMPLE=$(basename "${file}" _R1.trim.fastq.gz)

  # Run spades.py
  spades.py -k 27 \
            -1 ${INPUT_DIR}/${SAMPLE}_R1.trim.fastq.gz \
            -2 ${INPUT_DIR}/${SAMPLE}_R2.trim.fastq.gz \
            -o ${OUTPUT_DIR}/${SAMPLE} \
            -t 4 \
            -m 100 \
            |& tee ${OUTPUT_DIR}/${SAMPLE}.spades.log || echo "${SAMPLE} failed"
done
```
run the loop saved as run_spades.sh
```
 sbatch -w compute05 run_spades.sh
 ```
10. View the fasta files

i) sequence/contigs
```
less -S ./results/spades/contigs.fasta
```
ii) First 10 files/head
```
grep '>' ./results/spades/contigs.fasta | head
```
iii) Check thru
```
cat ./results/spades/contigs.fasta
```
11.Genome Assessment [input file contigs.fasta)

I) Genome contiguity

Checks length/cutoff for the longest contigs that contain 50% of the total genome length measured as contig N50. i.e  involves evaluating the accuracy and completeness of the genome assembly using metrics such as N50 length, scaffold and contig numbers, and genome size. Tool used QUAST.

a) For one sample
```
quast.py \
./results/spades/contigs.fasta \
-t 4 \
-o ./results/quast
```
"quast.py": This is the command to run the QUAST software.

"/results/spades/contigs.fasta": This is the path and filename of the input genome assembly file in FASTA format.

"-t 4": This specifies that the software should use four threads to run the analysis, which can speed up the process.

"-o /results/quast": This specifies the output directory where the results of the analysis will be stored.

Inspect the quast report

i) Download to your local wkd from hpc
```
cp -i AS-27566-C1-C_S23_L001/*.html ~/
```
ii) To local computer homepage/wkd
```
scp woguta@hpc.ilri.cgiar.org:~/AS-27566-C1-C_S23_L001.html .
```
b) Loop for all the samples
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J QUAST
#SBATCH -n 4

# Set the path to the directory containing the input files
input_dir=./results/spades

# Set the path to the directory where the output will be saved
output_dir=./results/quast

# Make the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop through all FASTA files in the input directory
for file in ${input_dir}/*.fasta; do
    filename=$(basename "$file")
    output_path="${output_dir}/${filename%.*}"
    
    # Run the quast.py command on the current file
    quast.py "$file" -t 4 -o "$output_path"
done
```
Savs and run
```
 sbatch -w compute05 run_quast.sh
 ```
In this loop, the for statement iterates over all files in the input directory that end with the extension .fasta. For each file, the loop extracts the filename using the basename command, and creates a corresponding output path by replacing the file extension with the output directory using ${filename%.*}. The quast.py command is then run on the current file, using 4 threads and the output path created earlier.

#!/usr/bin/bash -l: This is a shebang line that specifies the shell to use for the script.

#SBATCH -p batch: This is an SLURM directive that specifies the partition or queue to run the job in. In this case, the job will be submitted to the "batch" partition.

#SBATCH -J Quast: This sets the name of the job to "Quast".

#SBATCH -n 4: This specifies the number of CPU cores to allocate for the job. In this case, the job will use 4 cores.

input_dir=./results/spades: This sets the path to the directory containing the input files.

output_dir=./results/quast: This sets the path to the directory where the output will be saved.

for file in ${input_dir}/*.fasta; do: This begins a for loop that iterates over all files in the input directory that end with the extension .fasta.

filename=$(basename "$file"): This extracts the filename from the full file path.

output_path="${output_dir}/${filename%.*}": This creates the output path by replacing the file extension with the output directory.

quast.py "$file" -t 4 -o "$output_path": This runs the quast.py command on the current file, using 4 threads and the output path created earlier.

 II) Genome completeness
 
Assesses the presence or absence of highly conserved genes (orthologs) in an assembly/ ensures that all regions of the genome have been sequenced and assembled. It's performed using BUSCO (Benchmarking Universal Single-Copy Orthologs). Ideally, the sequenced genome should contain most of these highly conserved genes. they're lacking or less then the genome is not complete.

I) For one sample 

i) Using genome search (used here)
```
busco \
-i ./results/spades/contigs.fasta \
-m genome \
-o AS-27566-C1_S5_L001_busco \
-l bacteria \
-c 4 \
-f
```
ii) for unknown spp using busco bacterial database 
```
busco \
-i ./results/spades/contigs.fasta \
-m bacteria \
-o AS-27566-C1-C_S23_L001_busco \
-l bacteria_odb10 \
-c 4 \
-f
```
This command will run BUSCO on the contigs.fasta file using the "bacteria_odb10" database, outputting the results to a directory called AS-27566-C1-C_S23_L001_busco, using 4 CPUs, and overwriting any previous results in that directory (-f flag).

i) View short summary 
```
less -S AS-27566-C1_S5_L001_busco/short_summary.specific.bacteria_odb10.AS-27566-C1_S5_L001_busco.txt
```

ii) View full summary 
```
 less -S AS-27566-C1_S5_L001_busco/run_bacteria_odb10/full_table.tsv
```
iii) List and view a amino acid of protein sequence
```
 ls -lht AS-27566-C1_S5_L001_busco/run_bacteria_odb10/busco_sequences/
```
b) Loop for all files?
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J BUSCO
#SBATCH -n 4

# Load modules
module load busco/5.2.2

# Define input and output directories
INPUT_DIR=./results/spades
OUTPUT_DIR=./results/busco

#make output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Run BUSCO on all files in input directory
for file in ${INPUT_DIR}/*.fasta
do
  # Extract sample name from file name
  SAMPLE=$(basename "${file}" .fasta)

  # Run BUSCO
  busco \
    -i ${file} \
    -m genome \
    -o ${OUTPUT_DIR}/${SAMPLE}_busco \
    -l bacteria \
    -c 4 \
    -f \
    |& tee ${OUTPUT_DIR}/${SAMPLE}_busco.log || echo "${SAMPLE} failed"
done
```
Save & run_busco.sh
```
 sbatch -w compute05 run_busco.sh
 ```
The input_dir and output_dir variables are set to the input and output directories, respectively, and the database variable is set to the name of the BUSCO lineage-specific database to use. The script loops through all FASTA files in the input_dir directory, sets the output path for each file, and then runs the busco command on the current file using the appropriate parameters. The output is saved to a directory named after the input file, with the ".busco" suffix added.
 
III) Genome annotation

In genome annotation, the goal is to identify and label the features of on a genome sequence.

First, unload the modules to avoid conflict between the loaded dependencies
```
module purge
```
```
module load prokka/1.11
```
```
cd ./results/prokka
```
i) Run for one sample
```
prokka ./results/spades/contigs.fasta \
--outdir ./results/prokka \
--cpus 4 \
--mincontiglen 200 \
--centre C \
--locustag L \
--compliant \
--force
```

In this command, ./results/spades/contigs.fasta is the path to the input genome assembly file, --outdir specifies the output directory, --cpus specifies the number of CPUs to use, --mincontiglen specifies the minimum length of contigs to keep, --centre sets the sequencing centre abbreviation in the GenBank output, --locustag sets the locus tag prefix, --compliant enforces compliance with NCBI submission guidelines, and --force overwrites any existing output files.

Protein abundance
```
grep -o "product=.*" L_*.gff | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt
```
This command searches for the string "product=" in all files in the current directory that start with "L_" and have the extension ".gff". Then it uses the sed command to remove the string "product=" from each line of output. The resulting lines are sorted, counted, sorted in reverse numerical order, and finally written to a file called "protein_abundances.txt".

To break it down further, here's what each part of the command does:
grep -o "product=.*" L_*.gff: searches for the string "product=" in all files that match the pattern "L_*.gff", and outputs only the matching parts of the lines.
sed 's/product=//g': removes the string "product=" from each line of output.
sort: sorts the lines of output alphabetically.
uniq -c: counts the number of occurrences of each unique line of output.
sort -nr: sorts the lines of output by the count, in reverse numerical order (i.e., highest count first).
> protein_abundances.txt: writes the output to a file called "protein_abundances.txt".

View protein abundances
```
less -S protein_abundances.txt
```
ii) fo all files.
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J Prokka
#SBATCH -n 4

# Set the path to the directory containing the input files
input_dir=./results/spades

# Set the path to the directory where the output will be saved
output_dir=./results/prokka

# Make the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop through all FASTA files in the input directory
for file in ${input_dir}/*.fasta; do
    filename=$(basename "$file")
    output_path="${output_dir}/${filename%.*}"
    
    # Run the Prokka command on the current file
    prokka "$file" --outdir "$output_path" --cpus 4 --mincontiglen 200 --centre C --locustag L --compliant --force
done
```
Save and run
```
 sbatch -w compute05 run_prokka.sh
 ```

12. Species Identification
```
module load blast/2.12.0+
```
```
cd ./results/blast
```
i) For one sample
```
blastn \
-task megablast \
-query ./results/spades/contigs.fasta \
-db /export/data/bio/ncbi/blast/db/v5/nt \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-culling_limit 5 \
-num_threads 4 \
-evalue 1e-25 \
-out ./results/blast/contigs.fasta.vs.nt.cul5.1e25.megablast.out
```

This is a command to run blastn with the megablast algorithm on the contigs.fasta file generated by the SPAdes assembler. It searches against the non-redundant nucleotide database (nt) from NCBI using a strict output format (outfmt) that includes the query sequence ID, taxonomy IDs, bit score, standard deviation, scientific names, kingdom, and sequence title. It sets a culling limit of 5 to remove redundant hits and uses 4 threads for the search. The e-value threshold is set to 1e-25, and the output is saved to the contigs.fasta.vs.nt.cul5.1e25.megablast.out file in the blast directory.

-task megablast: uses the megablast algorithm, which is optimized for highly similar sequences.
-query ./results/spades/contigs.fasta: specifies the input FASTA file containing the assembled contigs to be searched.
-db /export/data/bio/ncbi/blast/db/v5/nt: specifies the path to the NCBI nt database to be searched.
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle': specifies the output format as a custom tabular format with columns for query sequence ID, taxonomic IDs of hits, bitscore, standard deviation, scientific names, kingdom names, and hit descriptions.
-culling_limit 5: filters out hits that have more than 5 other hits with better scores.
-num_threads 4: specifies the number of threads to be used for the search.
-evalue 1e-25: sets the e-value threshold for reporting significant hits to 1e-25.
-out ./results/blast/contigs.fasta.vs.nt.cul5.1e25.megablast.out: specifies the output file for the BLASTN results.

ii) Loop
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J blastn
#SBATCH -n 4

# Set the input directory
input_dir="./results/spades"

# Set the output directory
output_dir="./results/blast"

# Set the BLAST database path
db_path="/export/data/bio/ncbi/blast/db/v5/nt"

# Make the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop over all files in the input directory
for file in "${input_dir}"/*.fasta; do
    # Extract the filename without the extension
    filename=$(basename "${file%.*}")
    
    # Perform the BLASTN search
    blastn -task megablast \
        -query "${file}" \
        -db "${db_path}" \
        -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
        -culling_limit 5 \
        -num_threads 4 \
        -evalue 1e-25 \
        -out "${output_dir}/${filename}.vs.nt.cul5.1e25.megablast.out"
done
```
Save and run_blastn_search.sh
```
run_blastn.sh
```

View the blast
```
less -S ./results/blast/contigs.fasta.vs.nt.cul5.1e25.megablast.out
```
13. AMR Identification

Use Resistance Gene Identifier (RGI) which applies Comprehensive Antibiotic Resistance Database (CARD) as a reference to predict antibiotic resistome(s) from protein or nucleotide data based on homology and SNP models.
```
module purge

module load rgi/6.0.2
```
```
cd ./results/rgi

ln -s /var/scratch/global/bacteria-wgs/databases/localDB .
```
This creates symbolic link to the RGI database located in /var/scratch/global/bacteria-wgs/databases/localDB in the ./results/rgi directory using the ln -s command, to be used as a reference database for RGI for comparison against the contigs or assembled genomes fromthe bacterial WGS data to identify antibiotic resistance genes.

i) Run AMR for one sample
```
# Perform RGI analysis

rgi main --input_sequence ./results/spades/contigs.fasta \
--output_file ./results/rgi/AS-27566-C1-C_S23_L001_rgi \
--local \
-a BLAST \
-g PRODIGAL \
--clean \
--low_quality \
--num_threads 4 \
--split_prodigal_jobs
	
	
# samples and AMR genes organized alphabetically:
rgi heatmap --input ./results/rgi \
--output ./results/rgi/AS-27566-C1-C_S23_L001_rgi_alphabetic.png
```
The first command is using the rgi main command to run RGI analysis on the assembled contigs file (./results/spades/contigs.fasta). The results will be saved in a file called AS-27566-C1-C_S23_L001_rgi in the ./results/rgi directory. The --local option indicates that the local database will be used, the -a BLAST option specifies the use of BLAST algorithm for homology searches, and the -g PRODIGAL option specifies the use of PRODIGAL for gene prediction. The --clean and --low_quality options are used for quality control, and the --split_prodigal_jobs option is used to split the gene prediction jobs into multiple threads to speed up the process.

The second command is using the rgi heatmap command to generate a heatmap of the RGI results. The --input option specifies the directory where the RGI results are saved, and the --output option specifies the file name and location of the heatmap image. The heatmap will be organized alphabetically by sample and AMR gene. The resulting image will be saved as a PNG file in the ./results/rgi directory with the name AS-27566-C1-C_S23_L001_rgi_alphabetic.png.

Loop for all samples
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J rgi
#SBATCH -n 4

# Load modules
module load rgi/5.1.0

# Define input and output directories
INPUT_DIR=./results/spades
OUTPUT_DIR=./results/rgi

# Create output directory if it does not exist
mkdir -p "${OUTPUT_DIR}"

# Iterate over all files in input directory
for file in ${INPUT_DIR}/*.fasta
do
  # Extract sample name from file name
  SAMPLE=$(basename "${file}" .fasta)
  
  # Perform RGI analysis
  rgi main --input_sequence ${file} \
  --output_file ${OUTPUT_DIR}/${SAMPLE}_rgi \
  --local \
  -a BLAST \
  -g PRODIGAL \
  --clean \
  --low_quality \
  --num_threads 4 \
  --split_prodigal_jobs

  # Generate heatmap for AMR genes
  rgi heatmap --input ${OUTPUT_DIR}/${SAMPLE}_rgi \
  --output ${OUTPUT_DIR}/${SAMPLE}_rgi_heatmap.png
done
```
This script will loop over all files in the input directory with the extension .fasta, extract the sample name from the file name, perform RGI analysis on each file, and generate a heatmap for the AMR genes. The output files will be saved in the output directory with the sample name and appropriate extensions.

Summarize the results using rgi tab
```
rgi tab --input ./results/rgi/AS-27566-C1-C_S23_L001_rgi.txt --output ./results/rgi/AS-27566-C1-C_S23_L001_rgi_summary.tsv
```
This creates a tab-delimited file named AS-27566-C1-C_S23_L001_rgi_summary.tsv in the ./results/rgi/ directory. The file is viewd in a text editor or spreadsheet software like Microsoft Excel or Google Sheets. The file contains information about the predicted AMR genes in the sample, including their names, types, and percent identities.

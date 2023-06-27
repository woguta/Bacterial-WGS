# Analysis of llumina WGS data
## Illumina Pipeline
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

## NB

## To change directories and files timestamps
```
# Display the initial timestamps
echo "Initial timestamps:"
ls -lR

# Touch all directories
find . -type d -exec touch {} +

# Touch all files
find . -type f -exec touch {} +

# Display the updated timestamps
echo "Updated timestamps:"
ls -lR
```
# Steps
1. Log into hpc via ssh & interactive computing
```
interactive -w compute06 -c 3
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
 sbatch -w compute06 run_fastp1.sh
```
iii) Second pathway sbatch - sunmitiing jobs to the cluster as you do other things (the best option)
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J fastp
#SBATCH -n 4

#module purge
module purge

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
sbatch -w compute06 run_fastp2.sh
```
Checking the submitted job
```
ls -lth
```
Check the nano slurm number and squeue

```
nano slurm(number)

```
or 
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
sbatch -w compute06 run_fastqc_trim.sh
```
9. De-novo Genome Assembly using SPAdes for the trimmed files

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

# module purge
module purge

# Load modules
module load spades/3.15

# Define input and output directories
INPUT_DIR=./results/fastp
OUTPUT_DIR=./results/spades

#make output directory
mkdir -p "${OUTPUT_DIR}"

# Iterate over all files in input directory
for file in ${INPUT_DIR}/*.R1.trim.fastq.gz
do
  # Extract sample name from file name
  SAMPLE=$(basename "${file}" .R1.trim.fastq.gz)

# Check if the corresponding read 2 file exists
  file2=${INPUT_DIR}/${SAMPLE}.R2.trim.fastq.gz
  if [ ! -f ${file2} ]; then
    echo "Error: ${file2} not found!"
    continue
  fi
  
  # Make output directory for this sample
  mkdir -p "${OUTPUT_DIR}/${SAMPLE}"
  
 # Run spades.py
  spades.py -k 27 \
            -1 ${INPUT_DIR}/${SAMPLE}.R1.trim.fastq.gz \
            -2 ${INPUT_DIR}/${SAMPLE}.R2.trim.fastq.gz \
            -o ${OUTPUT_DIR}/${SAMPLE} \
            -t 4 \
            -m 100 \
            |& tee ${OUTPUT_DIR}/${SAMPLE}.spades.log || echo "${SAMPLE} failed"
done
```
run the loop saved as run_spades.sh
```
 sbatch -w compute06 -c 4 run_spades_allsamples.sh
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

iv) For plasmid assembly
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J PlasmidSPAdes
#SBATCH -n 8

# module purge
module purge

# Load modules
module load spades/3.15

# Define input and output directories
input_dir=./results/fastp
output_dir=./results/plasmid_assembly

# Make output directory
mkdir -p "${output_dir}"

# Iterate over all files in the input directory
for file in ${input_dir}/*.R1.trim.fastq.gz; do
  # Extract sample name from file name
  sample=$(basename "${file}" .R1.trim.fastq.gz)

  # Check if the corresponding read 2 file exists
  file2=${input_dir}/${sample}.R2.trim.fastq.gz
  if [ ! -f ${file2} ]; then
    echo "Error: ${file2} not found!"
    continue
  
  fi
  
  # Make output directory for this sample
  mkdir -p "${output_dir}/${sample}"

  # Run SPAdes for plasmid assembly
  spades.py --plasmid \
            -1 ${input_dir}/${sample}.R1.trim.fastq.gz \
            -2 ${input_dir}/${sample}.R2.trim.fastq.gz \
            -o ${output_dir}/${sample} \
            -t 8 \
            -m 100 \
            --debug
done
```
vii) Identify the plasmids/annotate

When plasmid contigs.fasta are available:
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J PlasmidID
#SBATCH -n 16

# Enable debugging and abort on error
 set -e

# Load modules
module purge
module load plasmidid/1.6.3

#Define input, output directories
input_dir="./results/fastp"
output_db="./plasmidIDdb"
dateNow=`date +"%Y-%m-%d"`
#dateNow="2023-06-21"

# Make output directory
mkdir -p "${output_db}"

echo "Downloading PlasmidIDdb"

#Download plasmidID database and define db path
download_plasmid_database.py -o "${output_db}"
#plasmid_db_path="./plasmidIDdb/2023-06-21_plasmids.fasta"
plasmid_db_path="./plasmidIDdb/${dateNow}_plasmids.fasta"

echo "Running PlasmidID"

for R1 in "${input_dir}"/*.R1.trim.fastq.gz; do
  if [ -f "${R1}" ]; then
    # Extract the sample name
    sample=$(basename "${R1}" .R1.trim.fastq.gz)

    # Define R2 file path
    R2="${input_dir}/${sample}.R2.trim.fastq.gz"

    # Define the contigs path
    contigs_dir="./results/plasmid_assembly/${sample}"
    contigs_fasta="${contigs_dir}/contigs.fasta"

    # Check if contigs directory exists
    if [ ! -d "${contigs_dir}" ]; then
      echo "ERROR: Contigs directory does not exist for sample ${sample}. Skipping sample ${sample}."
      continue
    fi

    # Check if contigs.fasta is readable and non-empty
    if [ ! -r "${contigs_fasta}" ] || [ ! -s "${contigs_fasta}" ]; then
      echo "ERROR: '${contigs_fasta}' is not a readable non-empty FASTA file for sample ${sample}. Skipping sample ${sample}."
      continue
    fi

    # Display contigs file path
    echo "Contigs for ${sample}: ${contigs_fasta}"

    # Run PlasmidID
    echo "Running PlasmidID for sample ${sample}"
    echo "Plasmid database path: ${plasmid_db_path}"
    plasmidID \
      -1 "${R1}" \
      -2 "${R2}" \
      -d "${plasmid_db_path}" \
      -c "${contigs_fasta}" \
      --no-trim \
      -g PSEUDO \
      -s "${sample}" \
      -T 16
  else
    echo -e "\tERROR!!! File not found: ${input_dir}/*R1.trim.fastq.gz"
    exit
  fi
done
```

When doing plasmid contigs.fasta assembly (Recommended!)

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J plasmidID
#SBATCH -n 16

# Exit with error reporting
set -e

# Load modules
module purge
module load plasmidid/1.6.3

# Define/make directories and database path
input_dir="./results/fastp"
plasmid_db_path="./plasmidIDdb/2023-06-21_plasmids.fasta"

echo "Running PlasmidID"

for R1 in "${input_dir}"/*.R1.trim.fastq.gz; do
  if [ -f "${R1}" ]; then
    # Extract the sample name
    sample=$(basename "${R1}" .R1.trim.fastq.gz)

    # Define R2 file path
    R2="${input_dir}/${sample}.R2.trim.fastq.gz"

    # Create the output directory for plasmid assembly
    assembly_output_dir="${output_dir}/${sample}"
    mkdir -p "${assembly_output_dir}"

    # Run PlasmidID for plasmid assembly
    echo "Running PlasmidID for sample ${sample}"
    echo "Plasmid database path: ${plasmid_db_path}"
    echo "Output directory for plasmid assembly: ${assembly_output_dir}"
    plasmidID \
      -1 "${R1}" \
      -2 "${R2}" \
      -d "${plasmid_db_path}" \
      -s "${sample}" \
      -g PSEUDO \
      --no-trim \
      -T 16
  else
    echo -e "\tERROR!!! File not found: ${input_dir}/*R1.trim.fastq.gz"
    exit
  fi
done
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
Save and run
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

In genome annotation, the goal is to identify and label the features of/on a genome sequence.

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
 sbatch -w compute05 c4 run_prokka.sh
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

#module purge to avoid clashing
module purge

#module load
module load blast/2.12.0+

# Set the input directory
input_dir="./results/spades"

# Set the output directory
output_dir="./results/blast"

# Set the BLAST database path
db_path="/export/data/bio/ncbi/blast/db/v5/nt"

# Make the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop over all files in the input directory
for file in "${input_dir}"/contigs.fasta; do
    # Extract the filename without the extension
    filename=$(basename "${file%.*}")
    
# Make output directory for this sample
 mkdir -p "${OUTPUT_DIR}/${filename}"
  
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
sbatch run_blastn.sh
```
iii) Doing all the *.contigs.fasta & *.scaffold.fasta
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J blastn
#SBATCH -n 4

#module purge to avoid clashing
module purge

#module load
module load blast/2.12.0+

# Set the input directory
input_dir="./results/spades"

# Set the output directory
output_dir="./results/blast"

# Set the BLAST database path
db_path="/export/data/bio/ncbi/blast/db/v5/nt"

# Make the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop over all sample directories in the input directory
for sample_dir in "${input_dir}"/*/; do
    # Extract the sample name from the directory path
    sample=$(basename "${sample_dir}")
    
    # Make output directory for this sample
    mkdir -p "${output_dir}/${sample}"
  
    # Perform the BLASTN search against scaffold.fasta
    blastn -task megablast \
           -query "${sample_dir}/scaffolds.fasta" \
           -db "${db_path}" \
           -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
           -culling_limit 5 \
           -num_threads 4 \
           -evalue 1e-25 \
           -out "${output_dir}/${sample}/scaffolds.vs.nt.cul5.1e25.megablast.out"
    
    # Perform the BLASTN search against contigs.fasta
    blastn -task megablast \
           -query "${sample_dir}/contigs.fasta" \
           -db "${db_path}" \
           -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
           -culling_limit 5 \
           -num_threads 4 \
           -evalue 1e-25 \
           -out "${output_dir}/${sample}/contigs.vs.nt.cul5.1e25.megablast.out"
done
```
This assumes that the input_dir/
├── sample1/
│   ├── scaffolds.fasta
│   └── contigs.fasta
├── sample2/
│   ├── scaffolds.fasta
│   └── contigs.fasta
├── sample3/
│   ├── scaffolds.fasta
│   └── contigs.fasta
└── ...

View the blast
```
less -S ./results/blast/${sample}/contigs.fasta.vs.nt.cul5.1e25.megablast.out
```
```
less -S ./results/blast/${sample}/scaffolds.fasta.vs.nt.cul5.1e25.megablast.out
```
iv) Extracting blast data
```
#load neccessary modules
module purge
module load R/4.2

#open r console/studio in the terminal
R

#remove all list saved in R in the terminal
rm(list = ls())

# Set the input directory where the sample directories are located
input_dir <- "./results/blast"

# Get a list of all sample directories in the input directory
sample_dirs <- list.dirs(input_dir, recursive = FALSE)

# Create a vector of available blast result files
blast_data_types <- c("contigs", "scaffolds")

for (blast_data_type in blast_data_types) {
  # Create an empty data frame to store the aggregated data
  agg_df <- data.frame()

  # Loop over each sample directory
  for (sample_dir in sample_dirs) {
    # Get a list of all blast result files in the current sample directory
    file_list <- list.files(path = sample_dir, pattern = paste0(blast_data_type, ".vs.nt.cul5.1e25.megablast.out"), full.names = TRUE)

    # Loop over each blast result file
    for (file in file_list) {
      # Read the blast data from the file
      blast_data <- read.delim(file, sep = "\t", header = FALSE)

      # Skip sample if no blast data was found
      if (nrow(blast_data) == 0) {
        next
      }

      # Extract the sample name from the file path
      sample_name <- tools::file_path_sans_ext(basename(sample_dir))

      # Add the sample name as a column in the data frame
      blast_data$sample_name <- sample_name

      # Append the data to the aggregated data frame
      agg_df <- rbind(agg_df, blast_data)
    }
  }

  # Set the output file name
  output_file <- paste0("combined_", blast_data_type, "_blast_data.csv")

  # Write the aggregated data to a CSV file
  write.csv(agg_df, file = output_file, row.names = FALSE)
}}
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
Summarize the results using rgi Tab
```
rgi tab --input ./results/rgi/AS-27566-C1-C_S23_L001_rgi.txt --output ./results/rgi/AS-27566-C1-C_S23_L001_rgi_summary.tsv
```
This creates a tab-delimited file named AS-27566-C1-C_S23_L001_rgi_summary.tsv in the ./results/rgi/ directory. The file is viewd in a text editor or spreadsheet software like Microsoft Excel or Google Sheets. The file contains information about the predicted AMR genes in the sample, including their names, types, and percent identities.

ii) Loop for all samples
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J RGI
#SBATCH -n 4

#module purge to avoid conflicts
module purge

#load required modules
module load rgi/6.0.2

#soft landing to the required RGI database within the pwd
ln -s /var/scratch/gkibet/card/localDB .

# Define input and output directories
input_dir="./results/spades"
output_dir="./results/rgi"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Loop over all sample directories in the input directory
for sample_dir in "${input_dir}"/*/; do
    # Extract the sample name from the directory path
    sample=$(basename "${sample_dir}")

    # Make output directory for this sample if doesn't exit
    mkdir -p "${output_dir}/${sample}"

    # Perform RGI analysis on contigs
    rgi main \
        -i "${sample_dir}/contigs.fasta" \
        -o "${output_dir}/${sample}/contigs.rgi.tsv" \
        -t contig \
        -- ./database/localDB \
        -a BLAST \
        -g PRODIGAL \
        --low_quality \
        --num_threads 4 \
        --split_prodigal_jobs \
	--clean \
	--debug        

    # Generate heatmap for AMR genes for contigs
    rgi heatmap \
        --input "${output_dir}/${sample}/contigs.rgi.tsv" \
        --output "${output_dir}/${sample}_heatmap.png" \
        --debug
done
```
This  script performs RGI (Resistance Gene Identifier) analysis on contigs and scaffolds generated from genome assembly for multiple samples. The script first purges modules and then loads the necessary module RGI, and sets up the RGI database by creating a symbolic link to the CARD (Comprehensive Antibiotic Resistance Database) localDB.

The script then defines input and output directories for the analysis and creates output directories for each sample. The RGI analysis is performed on both contigs and scaffolds, with the results output to separate files for each sample. A heatmap is also generated for each sample showing the presence of AMR (Antimicrobial Resistance) genes.

After processing all the samples, the script combines the results from all contigs.rgi.csv files into a single file called combined_contigs_rgi_csv, and all scaffolds.rgi.csv files into a single file called combined_scaffolds_rgi_csv. Finally, the script merges the combined_contigs_rgi_csv and combined_scaffolds_rgi_csv files into a single file called combined_contigs_scaffolds_rgi_csv.

iii) Run abricate amr genes analysis for all samples

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J abricate
#SBATCH -n 4

#module purge to avoid clashing
module purge

# Load the required modules
module load abricate/1.0.1

# Set the input directory
input_dir="./results/spades"

# Set the output directory
output_dir="./results/abricate"

# Set the ABRICATE databases to use
databases=("ncbi" "card" "argannot" "resfinder" "megares" "plasmidfinder" "vfdb")

# Set the ABRICATE and RGI minimum identity and coverage thresholds
min_identity="80"
min_coverage="60"

# Make the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Set the summary file name
summary_file="${output_dir}/abricate_rgi_summary.csv"

# Create the header for the summary file
echo "sample,db,query_id,accession,coverage,identity,gene,start,end,hit_length,amr_gene_family,drug_class,mechanism,group,reference" > "${summary_file}"

# Loop over all sample directories in the input directory
for sample_dir in "${input_dir}"/*/; do
    # Extract the sample name from the directory path
    sample=$(basename "${sample_dir}")
    
    # Make output directory for this sample
    mkdir -p "${output_dir}/${sample}"
  
    # Run ABRICATE on the contigs files
    for db in "${databases[@]}"; do
        abricate --db "${db}" \
                 --minid "${min_identity}" \
                 --mincov "${min_coverage}" \
                 "${sample_dir}/contigs.fasta" \
                 > "${output_dir}/${sample}/contigs_${db}.abricate.tsv"
done
```

##Extracting AMR data

i) From RGI results
```
#!usr/bin/bashrc -l

#load neccessary modules
module purge
module load R/4.2

#open r console/studio in the terminal
R

# Set the input directory where the sample directories are located
input_dir <- "./results/rgi"

# Get a list of all sample directories in the input directory
sample_dirs <- list.dirs(input_dir, recursive = FALSE)

# create a vector of available for rgi results
rgi_data_types <- c("json", "txt")

for (rgi_data_type in rgi_data_types) {

# Create an empty data frame to store the aggregated data
agg_df <- data.frame()

# Loop over each sample directory
for (sample_dir in sample_dirs) {
  # Get a list of all rgi result files in the current sample directory
  file_list <- list.files(path = sample_dir, pattern = paste0("contigs_",rgi_data_type,".rgi.tsv.json|txt"), full.names = TRUE)

  # Loop over each rgi result file
  for (file in file_list) {
    # Read the rgi data from the file
    rgi_data <- read.delim(file, sep = "\t")
	
	# skip sample if no rgi data was found
	if (nrow(rgi_data) == 0) {
	next
	}

    # Extract the sample name from the file path
    sample_name <- tools::file_path_sans_ext(basename(sample_dir))

    # Add the sample name as a column in the data frame
    rgi_data$sample_name <- sample_name

    # Append the data to the aggregated data frame
    agg_df <- rbind(agg_df, rgi_data)
  }
}

# Set the output file name
output_file <- paste0("combined_",rgi_data_type,"_rgi_data.csv")

# Write the aggregated data to a CSV file
write.csv(agg_df, file = output_file, row.names = FALSE)
}
done
```

ii) From ABRicate results
```
ABRicate
#!usr/bin/bashrc -l

#load neccessary modules
module purge
module load R/4.2

#open R console/studio in the terminal
R

#Path to the Rscript
#export PATH=".:$PATH"

# Set the input directory where the sample directories are located
input_dir <- "./results/abricate"

# Get a list of all sample directories in the input directory
sample_dirs <- list.dirs(input_dir, recursive = FALSE)

# create a vector of available abricate results
abricate_data_types <- c("ncbi", "argannot", "card", "plasmidfinder", "resfinder", "vfdb", "megares")

for (abricate_data_type in abricate_data_types) {

# Create an empty data frame to store the aggregated data
agg_df <- data.frame()

# Loop over each sample directory
for (sample_dir in sample_dirs) {
  
  # Get a list of all abricate result files in the current sample directory
  file_list <- list.files(path = sample_dir, pattern = paste0("contigs_",abricate_data_type,".abricate.tsv"), full.names = TRUE)

  # Loop over each abricate result file
   for (file in file_list) {
    
    # Read the abricate data from the file
    abricate_data <- read.delim(file, sep = "\t")

    # skip sample if no abricate data was found
        if (nrow(abricate_data) == 0) {
        next
        }

    # Extract the sample name from the file path
    sample_name <- tools::file_path_sans_ext(basename(sample_dir))

    # Add the sample name as a column in the data frame
    abricate_data$sample_name <- sample_name
 
   # Append the data to the aggregated data frame
    agg_df <- rbind(agg_df, abricate_data)
  }
}

# Set the output file name
output_file <- paste0("combined_",abricate_data_type,"_abricate_data.csv")

# Write the aggregated data to a CSV file
write.csv(agg_df, file = output_file, row.names = FALSE)
}
done
```

##Transforming/extracting AMR counts per sample
```
#load required libraries
library(tidyverse)

# read in the data
rgi_data <- read.csv("combined_json_rgi_data.csv")

# select the columns of interest
rgi_data_subset <- rgi_data %>% select(Best_Hit_ARO, sample_name)

# covert data to long format and get counts
ARO_counts <- rgi_data_subset %>%
  group_by(sample_name)%>%
  count(Best_Hit_ARO)

# convert the data to wide format
ARO_counts_wide <- ARO_counts %>% 
  pivot_wider(names_from = Best_Hit_ARO, values_from = n)

# convert NAs to zeros
ARO_counts_wide[is.na(ARO_counts_wide)] <- 0

# write the data to csv file
write.csv(ARO_counts_wide, "ARO_counts.csv", row.names = FALSE)
```
14. Identification of virulence factors

The ability of a bacterium to colonise a host and contribute to its pathogenecity, or ability to cause disease, are known as virulence factors. In order to gather data regarding the virulence factors of bacterial pathogens, we will make use of the virulence factor database (VFDB), which is an integrated and comprehensive resource. We will scan through our contigs against the VFDB quickly to find any sequences that include known virulence factors.

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J blast_vf
#SBATCH -n 4

#module purge
module purge

# Load modules
module load blast/2.11.0

# Define input and output directories
INPUT_DIR=./results/spades
OUTPUT_DIR=./results/blast

# Create output directory if it does not exist
mkdir -p "${OUTPUT_DIR}"

# Run BLAST
blastn \
    -query ${INPUT_DIR}/contigs.fasta \
    -db /var/scratch/global/bacteria-wgs/databases/VFDB/vfdb_seta_nt \
    -out ${OUTPUT_DIR}/virulence_factors_blast_VFDB.out \
    -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -num_threads 4
```
-query: specifies the input query file in fasta format.
-db: specifies the path to the VFDB blast database to be used for the search.
-out: specifies the output file name and location for the blast results.
-outfmt: specifies the format of the output file. In this case, it is set to format 6, which provides a tab-delimited output containing the query ID, subject taxonomy ID, bit score, standard deviation of the bit score, subject scientific name, subject kingdom, and subject title.
-num_threads: specifies the number of threads to be used for the search. In this case, it is set to 4.

Loop for all files
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J blast_vf
#SBATCH -n 4

#module purge
module purge

# Load modules
module load blast/2.11.0

# Define input and output directories
INPUT_DIR=./results/spades
OUTPUT_DIR=./results/blast

# Create output directory if it does not exist
mkdir -p "${OUTPUT_DIR}"

# Loop over all files in the input directory
for FILE in ${INPUT_DIR}/*.fasta
do
  # Get the sample name from the file name
  SAMPLE_NAME=$(basename ${FILE} .fasta)

  # Run BLAST
  blastn \
    -query ${FILE} \
    -db /var/scratch/${USER}/bacteria-wgs/databases/VFDB/vfdb_seta_nt \
    -out ${OUTPUT_DIR}/${SAMPLE_NAME}_virulence_factors_blast_VFDB.out \
    -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -num_threads 4
done
```
This loop uses a wildcard (*) to match all files with a .fasta extension in the INPUT_DIR. The loop then extracts the sample name from the file name (assuming that the file name ends with .fasta) and uses it to name the output file in the OUTPUT_DIR. The loop then runs the same blastn command for each sample, outputting the results to a separate file for each sample.

15. Genome reference assembly

a) Load our reference genome from NCBI, take fasta and gff files
	
- On a web browser, open the link NCBI.

- Type 'Pseudomonas aeruginosa' on the search box and select 'Genome' database.

- Right click on the genome FASTA and and gff, then select 'copy links'.

- Use wget to fetch the files as follows:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz .
```
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.gff.gz .
```
b) Gunzip and rename: pseudomonas_aeruginosa_pao1_substrain_genome.fasta and pseudomonas_aeruginosa_pao1_substrain.gff

c) Indexing the ref_genome using samtools faidx produces a .fai file consisting of five tab-separated columns: chrname, seqlength, first-base offset, seqlinewidth without \n (newline character) and seqlinewidth with\n. This is essential for samtools' operations and bcftools.

```
module load samtools/1.15.1
```
```
samtools faidx pseudomonas_aeruginosa_pao1_substrain_genome.fasta
```
d) In order to allow easy access of genome regions during read mapping we will index the reference genome using bowtie2, bowtie2-build outputs six files 1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2 all constituting the index and are all needed to align reads to the reference which will no longer be needed by bowtie after this indexing. The output is in binary format and can not be visualized like text.

```
# Load modules
module purge
module load bowtie2/2.5.0

# Define the input file and output prefix
input_file="./pseudomonas_aeruginosa_pao1_substrain_genome.fasta"
output_dir="./../results/bowtie/PAO1"

#Create output directory if it doesn't exist
mkdir -p "${output_dir}"

# Run Bowtie 2 indexing
bowtie2-build \
    --threads 4 \
    "${input_file}" \
    "${output_dir}"/PAO1
```
e) Run bowtie for reference genome assembly
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J bowtie2
#SBATCH -n 10

# Load modules
module purge
module load bowtie2/2.5.0
module load samtools/1.15.1

# Define the input_dir where the samples are located
input_dir="./results/fastp/"
output_dir="./results/bowtie/"

# Make output_dir if it doesn't exist
mkdir -p "${output_dir}"

# Loop over all sample files in input_dir
for sample_file in "${input_dir}"/*.R1.trim.fastq; do
    # Extract the sample name from the directory path
    sample=$(basename "${sample_file%.*}")

    # Define input file paths for the current sample
    index_prefix="./results/bowtie/PAO1/PAO1"
    fastq_1="${sample_file}"
    fastq_2="${sample_file/R1.trim.fastq/R2.trim.fastq}"
    output_prefix="./results/bowtie/${sample}"

    echo "Processing sample: ${sample}"
    echo "Fastq 1: ${fastq_1}"
    echo "Fastq 2: ${fastq_2}"
    echo "Index Prefix: ${index_prefix}"
    echo "Output prefix: ${output_prefix}"

    # Run Bowtie2 alignment and save BAM output
    bowtie2 -x ${index_prefix} \
            -1 ${fastq_1} \
            -2 ${fastq_2} \
            --threads 10 \
            --un-conc-gz ${output_prefix}.unmapped.fastq.gz \
            --local \
            --very-sensitive-local \
            2> ${output_prefix}.bowtie2.log \
            | samtools view -@ 10 -F4 -bhS -o ${output_prefix}.bam -
done
```
f) Sort the bam files output from bowtie results. In order to perform any analysis and visualize the alingment, the alignment should be sorted in the order which they occur in the reference genome based on their alignment coordinates. Indexing the sorted bam files allows viiewing using various tools like IGV to extract alignments or information like mutation. Soting gives sotred.bam file while indexing generates an index file with .bam.bai extention.

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J samtools_sort
#SBATCH -n 5

# Load modules
module purge
module load samtools/1.15.1

# Define input and output directories
input_dir="./results/bowtie/"
output_dir="./results/bowtie/sorted/"

# Make output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over all sample files in the input directory
for bam_file in "$input_dir"/*.R1.trim.fastq.bam; do
    # Extract the sample name from the file path
    sample=$(basename "${bam_file%.*}")

    # Define the output file path for the sorted BAM
    sorted_bam="$output_dir/${sample}.sorted.bam"

    echo "Sorting BAM for sample: ${sample}"

    # Sort the BAM file using samtools
    samtools sort -@ 4 \
        -o "$sorted_bam" \
        -T "$output_dir/${sample}" \
        "$bam_file"

    echo "BAM sorting completed for sample: ${sample}"

    #Index the sorted BAM file
    samtools index -@4 "$sorted_bam"
    
    echo "BAM indexing completed for sample: ${sample}"
done
```
g) Assessment/Viewing of the sequences using samtools view

-View the whole alignment 'BAM' file.

```
samtools view ./results/bowtie/sorted/${sample}.sorted.bam | less -S
```

-Count the number of reads in the alignment.
```
samtools view ./results/bowtie/sorted/${sample}.sorted.bam | wc -l
```

-Count the number of reads mapping for a gene e.g elastase lasB NP_252413.1

```
samtools view ./results/bowtie/sorted/${sample}.sorted.bam | grep "NP_252413.1" | wc -l
```
h) Computing genome coverage using bedtools
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J bedtools_coverage
#SBATCH -n 5

# Load modules
module purge
module load bedtools/2.29.0

# Define input and output directories
input_dir="./results/bowtie/sorted"
output_dir="./results/bowtie/coverage/"

# Make output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over all sorted BAM files in the directory
for bam_file in "$input_dir"/*.sorted.bam; do
    # Extract the sample name from the file name
    sample=$(basename "${bam_file%.*}")

    # Define the output file path
    output_file="$output_dir/${sample}.coverage"

    echo "Processing sample: ${sample}"
    echo "Output file: ${output_file}"

    # Calculate genome coverage using bedtools
    bedtools genomecov -d -ibam "$bam_file" > "$output_file"

    echo "Genome coverage calculated for ${sample}"
done
```

The output of this command is a three column file with chromosome, Position and depth coverage. View the ouput with the command below:

```
less ./results/bowtie/coverage/${sample}.coverage
```
h) Plot the genome coverage in R

```
# Load modules
module purge
module load R/4.3
R
library(ggplot2)

# Define input and output directories
input_dir <- "./results/bowtie/coverage"
output_dir <- "./results/bowtie/coverage/plots"

# Make output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Get a list of sorted coverage files in the directory
coverage_files <- list.files(path = input_dir, pattern = "\\.sorted\\.coverage$", full.names = TRUE)

# Loop over all sorted coverage files
for (coverage_file in coverage_files) {
  # Extract the sample name from the file name
  sample <- basename(tools::file_path_sans_ext(coverage_file))
  
  # Define the output file path
  output_file <- file.path(output_dir, paste0(sample, "_genome_coverage.png"))
  
  cat("Processing sample:", sample, "\n")
  cat("Output file:", output_file, "\n")
  cat("Generating genome coverage plot for", sample, "\n")
  
  # Read the coverage data into a data frame
  coverage_data <- read.table(coverage_file, header = FALSE, col.names = c("Chromosome", "Position", "Coverage"))
  
  # Plot the genome coverage
  ggplot2::ggplot(coverage_data, ggplot2::aes(x = Position, y = Coverage)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Position", y = "Coverage") +
    ggplot2::ggtitle("Genome Coverage") +
    ggplot2::theme_minimal()
  
  # Save the plot as a PNG file
  ggplot2::ggsave(output_file, width = 10, height = 6, dpi = 300)
  
  cat("Genome coverage plot generated for", sample, "\n")
}
```
14. Variants callings

i) Using Ivar and Samtools mpileup

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J ivar_variants_calling
#SBATCH -n 5

#load necessary modules
module purge
module load samtools/1.15.1
module load ivar/1.3.1

# Specify the directory containing the sorted BAM files
input_dir="./results/bowtie/sorted_indexed"
output_dir="./results/ivar_variants"

# Make output directory if not available
mkdir -p "${output_dir}"

# Iterate over each sorted BAM file in the directory
for bam_file in "$input_dir"/*.sorted.bam; do
    # Extract the file name without the extension
    filename=$(basename "$bam_file" .sorted.bam)

    # Run samtools mpileup and ivar variants for each BAM file
    samtools mpileup \
        --ignore-overlaps \
        --count-orphans \
        --no-BAQ \
        --max-depth 0 \
        --min-BQ 0 \
        --reference ./genome_ref/pseudomonas_aeruginosa_pao1_substrain_genome.fasta \
        "$bam_file" \
        | ivar variants \
            -t 0.25 \
            -q 20 \
            -m 10 \
            -g ./genome_ref/pseudomonas_aeruginosa_pao1_substrain.gff \
            -r ./genome_ref/pseudomonas_aeruginosa_pao1_substrain_genome.fasta \
            -p "${output_dir}/${filename}.variants"

    echo "Processed ${bam_file}"
done
```
ii) Comprehensive variant calling using GATK and Vep
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J GATKVep_VC
#SBATCH -n 10

# Load necessary modules
module purge
module load samtools/1.15.1
module load java/17
module load picard/2.27.5
module load gatk/4.4.0.0
module load vep/109.3

# Specify the directory containing the sorted BAM files
input_dir="./results/bowtie/sorted_indexed"
output_dir="./results/gatkvep_variants"

# Make output directory if not available
mkdir -p "${output_dir}"

# Specify the path to the reference genome and GFF files
reference="./genome_ref/pseudomonas_aeruginosa_pao1_substrain_genome.fasta"
gff="./genome_ref/pseudomonas_aeruginosa_pao1_substrain_genome.gff"

# Create the sequence dictionary file for the reference genome
gatk CreateSequenceDictionary -R "$reference" -O "${reference%.*}.dict"

# Specify the path to the sequence dictionary file
dict_file="${reference%.*}.dict"

# Create an empty sample map file
sample_map_file="${output_dir}/sample_map.txt"
> "$sample_map_file"

# Iterate over each sorted BAM file in the input directory
for bam_file in "$input_dir"/*.fastq.sorted.bam; do
    # Extract the file name without extension
    sample_name=$(basename "$bam_file" .sorted.bam)

    # Get the desired sample name based on the filename
    desired_sample_name="${sample_name}"

    # Update the desired sample name in the BAM file
    samtools view -H "$bam_file" | sed -e "s/SM:${sample_name}/SM:${desired_sample_name}/" | samtools view -bS - > "${output_dir}/${desired_sample_name}.sorted.bam"
 
   # Check if the BAM file was successfully updated
    if [ $? -ne 0 ]; then
        echo "Error: Failed to update sample name in BAM file: ${bam_file}"
        exit 1
    fi

    #Index the sorted BAM file
    samtools index -@4 "${output_dir}/${desired_sample_name}.sorted.bam"

    echo "Sample Name: ${desired_sample_name}.sorted.bam"
    echo "Indexed Sample Name: ${output_dir}/${desired_sample_name}.bam.fai"

    # Add the desired sample name to the sample map file
    echo "${desired_sample_name} ${output_dir}/${desired_sample_name}.sorted.bam" >> "$sample_map_file"
    echo "${desired_sample_name} ${output_dir}/${desired_sample_name}.sorted.bam.fai" >> "$sample_map_file"

    # Define input files
    input_file="${output_dir}/${desired_sample_name}.sorted.bam"

    # Define output file
    output_file="${output_dir}/${desired_sample_name}.vcf.gz"

     echo $reference
     echo $input_file
     echo $output_file
     echo $dict_file

    # Step 1: Variant Calling using GATK HaplotypeCaller
    gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R "$reference" \
        -I "$input_file" \
        -O "$output_file" \
        --sequence-dictionary "$dict_file"

    echo "Variant Calling completed for ${desired_sample_name}"
done

# Step 2: Joint Genotyping using GATK GenomicsDBImport and GenotypeGVCFs
gatk --java-options "-Xmx4g" GenomicsDBImport \
    --genomicsdb-workspace-path "${output_dir}/genomicsdb" \
    --sample-name-map "$sample_map_file"

gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R "$reference" \
    -V gendb://"${output_dir}/genomicsdb" \
    -O "${output_dir}/joint_genotyped.vcf.gz"

# Step 3: Variant Filtering using GATK VariantFiltration
gatk --java-options "-Xmx4g" VariantFiltration \
    -R "$reference" \
    -V "${output_dir}/joint_genotyped.vcf.gz" \
    --filter-expression "QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "PASS_FILTER" \
    -O "${output_dir}/filtered_variants.vcf.gz"

# Step 4: Variant Annotation using VEP (Variant Effect Predictor)
vep -i "${output_dir}/filtered_variants.vcf.gz" \
    -o "${output_dir}/annotated_variants.vcf" \
    --fasta "$reference" \
    --cache \
    --dir_cache /export/apps/vep/109.3/cache \
    --species "Pseudomonas aeruginosa" \
    --force_overwrite \
    --merged \
    --no_stats \
    --no_intergenic \
    --flag_pick_allele_gene \
    --symbol \
    --canonical \
    --biotype \
    --pubmed \
    --vcf_info_field "ANN" \
    --gff "$gff"

echo "Joint Genotyping and annotation completed"
echo "Done"
```
iii) a. Make snpEff database for P.aeruginosa PAO1
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J snpEff_db
#SBATCH -n 16

#Load modules
module purge
module load java/17

# Exit with error debug
set -e

mkdir -p ./database/snpEff/data/PAO1/
cp -r ./genome_ref/pseudomonas_aeruginosa_PAO1.gff ./database/snpEff/data/PAO1/genes.gff
cp -r ./genome_ref/pseudomonas_aeruginosa_PAO1.fasta ./database/snpEff/data/PAO1/sequences.fa
echo -e "# P. aeruginosa bacterial genome, version PseudomonasPAO1\nPAO1.genome: PAO1" > ./database/snpEff/data/PAO1/snpEff.config

#Build the database
java -Xmx4g -jar /export/apps/snpeff/5.1g/snpEff.jar build \
        -config ./database/snpEff/data/PAO1/snpEff.config \
        -dataDir ./../ \
        -gff3 \
        -v PAO1
```
iii) b. Comprehensive indexing of ref_genome, aligning reads, variants calling and annotation done at once using bwa, samtools, bcftools, snpeff &snsift
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J BwaBcf_VC
#SBATCH -n 16

# Exit immediately if a command exits with a non-zero status
set -e

# Load necessary modules
module purge
module load bwa/0.7.17
module load samtools/1.15.1
module load bcftools/1.15.1
module load snpeff/5.1
module load java/17

# Specify the directory containing the trimmed fastq files
input_dir="./results/fastp"
output_dir="./results/bwamem"
output_dir2="./results/bcf_variants"

# Define the path to the reference genome FASTA file
reference_file="./genome_ref/pseudomonas_aeruginosa_PAO1.fasta"

# Make output directories if not available
mkdir -p "${output_dir}"
mkdir -p "${output_dir2}"

# Define the output prefix for the index files and its directory
reference_prefix="${output_dir}/PAO1/PAO1"
mkdir -p "${reference_prefix}"

# Step 1: Index the reference genome
echo "Indexing the reference genome..."
bwa index -p "${reference_prefix}" -a bwtsw "${reference_file}"

# Step 2: Perform alignment and generate sorted BAM files for each sample
echo "Performing alignment..."
for sample_file in "${input_dir}"/*.R1.trim.fastq.gz; do
    sample_name=$(basename "${sample_file}" .R1.trim.fastq.gz)
    echo "Processing sample: ${sample_name}"
    bwa mem -t 16 "${reference_prefix}" "${sample_file}" "${input_dir}/${sample_name}.R2.trim.fastq.gz" |
    samtools view -bS - |
    samtools sort -@ 16 -o "${output_dir}/${sample_name}.sorted.bam" -
    samtools index "${output_dir}/${sample_name}.sorted.bam"
done

# Some stats
echo "Generating alignment statistics..."
for sorted_bam in "${output_dir}"/*.sorted.bam; do
    sample_name=$(basename "${sorted_bam}" .sorted.bam)
    samtools flagstat "${sorted_bam}"
done

# Step 3: Perform variant calling on each sample
echo "Performing variant calling..."
for sorted_bam in "${output_dir}"/*.sorted.bam; do
    sample_name=$(basename "${sorted_bam}" .sorted.bam)
    echo "Processing sample: ${sample_name}"
    bcftools mpileup -Ou -f "${reference_file}" "${sorted_bam}" |
    bcftools call --threads 16 --ploidy 1 -Ou -mv -o "${output_dir2}/${sample_name}.vcf"
done

# Step 4: Filter and report the SNVs variants in variant calling format (VCF)
echo "Filtering variants..."
for vcf_file in "${output_dir2}"/*.vcf; do
    sample_name=$(basename "${vcf_file}" .vcf)
    echo "Processing sample: ${sample_name}"
    bcftools filter --threads 16 -i 'DP>=10' -Ov "${vcf_file}" > "${output_dir2}/${sample_name}.filtered.vcf"
done

# Step 5: Index the filtered VCF files
echo "Indexing filtered VCF files..."
for filtered_vcf in "${output_dir2}"/*.filtered.vcf; do
    echo "Indexing: ${filtered_vcf}"
    bcftools index --threads 16 "${filtered_vcf}"
done

# Step 6: Variant Annotation with GFF file
# Define input files and directories
vcf_dir="${output_dir2}"
annotated_dir="./results/annotated_variants"
reference_file="./genome_ref/pseudomonas_aeruginosa_PAO1.fasta"
gff_file="./genome_ref/pseudomonas_aeruginosa_PAO1.gff"
snpeff_jar="/export/apps/snpeff/5.1g/snpEff.jar"

# Create output directory for annotated variants
mkdir -p "$annotated_dir"

# Run SnpEff for variant annotation
echo "Performing variant annotation..."
for vcf_file in "$vcf_dir"/*.filtered.vcf; do
    sample_name=$(basename "${vcf_file}" .filtered.vcf)
    echo "Processing sample: $sample_name"
    bcftools view --threads 16 "$vcf_file" |
    java -Xmx4g -jar "$snpeff_jar" \
           -config ./database/snpEff/data/PAO1/snpEff.config \
	   -dataDir ./../ \
	   -v PAO1 ${vcf_file} > "${annotated_dir}/${sample_name}.snpeff.vcf"
done

echo "Variant annotation complete."

# Step 7: Rename summary.html and genes.txt and zip vcf files
mv ./snpEff_summary.html "${annotated_dir}/${sample_name}.snpeff.summary.html"
mv ./snpEff_genes.txt "${annotated_dir}/${sample_name}.snpeff.genes.txt"

# Compress vcf
bgzip -c "${annotated_dir}/${sample_name}.snpeff.vcf" > "${annotated_dir}/${sample_name}.snpeff.vcf.gz"

# Create tabix index - Samtools
tabix -p vcf -f "${annotated_dir}/${sample_name}.snpeff.vcf.gz"

# Generate VCF files
bcftools stats "${annotated_dir}/${sample_name}.snpeff.vcf.gz" > "${annotated_dir}/${sample_name}.snpeff.stats.txt"

echo "Variant summar renaming, compressing complete."

# Step 8: Variant Extraction with SnpSift
# Define input files and directories
snpeff_dir="${output_dir2}"
extracted_dir="./results/extracted_variants"

# Create output directory for extracted variants
mkdir -p "$extracted_dir"

# Run SnpSift for variant extraction
echo "Performing variant extraction..."
for snpeff_file in "$snpeff_dir"/*.snpeff.vcf.gz; do
    sample_name=$(basename "${snpeff_file}" .snpeff.vcf.gz)
    echo "Processing sample: $sample_name"
    bcftools view -t 16 "$snpeff_file" |
    java -Xmx4g -jar "/export/apps/snpeff/5.1g/snpSift.jar" \
        extractFields \
        -s "," \
        -e "." \
        /dev/stdin \
        "ANN[*].GENE" "ANN[*].GENEID" \
        "ANN[*].IMPACT" "ANN[*].EFFECT" \
        "ANN[*].FEATURE" "ANN[*].FEATUREID" \
        "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \
        "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \
        "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \
        "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \
        "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \
        "LOF[*].GENE" "LOF[*].GENEID" "LOF[*].NUMTR" "LOF[*].PERC" \
        "NMD[*].GENE" "NMD[*].GENEID" "NMD[*].NUMTR" "NMD[*].PERC" \
        > "${extracted_dir}/${sample_name}.snpsift.txt"
done
```

iv) Extracting the annotated_variants from snpsift.txt

```
#!usr/bin/bashrc -l

#load neccessary modules
module purge
module load R/4.3

#open r console/studio in the terminal
R

# Set the input directory where the sample directories are located
input_dir <- "./extracted_variants"

# Get a list of all sample directories in the input directory
sample_dirs <- list.dirs(input_dir, recursive = FALSE)

# create a vector of available for rgi results
annotated_data_type <- c("txt")

for annotated_data_type {

# Create an empty data frame to store the aggregated data
agg_df <- data.frame()

# Loop over each sample directory
for (sample_dir in sample_dirs) {
  # Get a list of all annotated result files in the current sample directory
  file_list <- list.files(path = sample_dir, pattern = paste0("snpsift.txt"), full.names = TRUE)

  # Loop over each annotated result file
  for (file in file_list) {
    # Read the annotated data from the file
    annotated_data <- read.delim(file, sep = "\t")
	
	# skip sample if no annotated data was found
	if (nrow(annotated_data) == 0) {
	next
	}

    # Extract the sample name from the file path
    sample_name <- tools::file_path_sans_ext(basename(sample_dir))

    # Add the sample name as a column in the data frame
    annoatated_data$sample_name <- sample_name

    # Append the data to the aggregated data frame
    agg_df <- rbind(agg_df, annotated_data)
  }
}

# Set the output file name
output_file <- paste0("combined_annotated_data.csv")

# Write the aggregated data to a CSV file
write.csv(agg_df, file = output_file, row.names = FALSE)
}
done

echo "Variant extraction complete."
```

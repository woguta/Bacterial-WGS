# The ONT WGS pipeline

Is a bioinformatics pipeline designed for analyzing bacterial genomes using Oxford Nanopore Technologies (ONT) sequencing data. The pipeline includes several steps, which are outlined below:

1.  Basecalling: Raw ONT sequencing data (in the form of electrical signals) are converted into nucleotide sequences using basecalling algorithms such as Guppy.

2.  Quality control: The quality of the basecalled sequences is assessed using tools such as FastQC or NanoPlot. Sequences with poor quality are filtered out. Here NanoPlot is used

3.  Assembly: The filtered sequences are assembled into contigs using a variety of tools such as Flye, Canu, and Shasta.

4.  Polishing: The assembled contigs are polished to improve their accuracy and reduce errors using tools such as Racon and Medaka.

5.  Scaffold: The polished contigs are scaffolded to link the contigs and produce a more complete genome sequence using tools such as SSPACE-LongRead.

6.  Gap Filling: Gaps between contigs are filled using tools such as Pilon.

7.  Annotation: The assembled and polished genome sequence is annotated to identify genes, functional elements, and other features using tools such as Prokka and RAST.

8.  Comparative analysis: The annotated genomes are compared to reference genomes and other bacterial genomes to identify similarities, differences, and evolutionary relationships using tools such as Mauve and Roary.

9.  Visualization: The results of the analysis are visualized using tools such as Artemis and Circos.

# Steps

1.  Log into hpc via ssh & interactive computing
```
interactive -w compute05 -c 3
```

2.  Creating directory
```
cd /var/scratch/
mkdir -p $USER/bacteria-wgs/flair
cd $USER/bacteria-wgs/flair
```

The mkdir -p command is used to create a directory and its parent directories (if they do not exist) in a single command.

The -p option allows the mkdir command to create the parent directories if they do not exist.

3.  .  Data retrieval

i) Use scp username@remote:/path/to/data /path/to/local/directory

```
scp -r woguta@hpc.ilri.cgiar.org:/path/to/data  .
```

The -r option is used to copy the directory recursively. The . at the end of the command specifies the current directory as the destination. space fullstop " ." specifies your current folder.

ii) Base calling using gyppyplex5+
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J guppy6+
#SBATCH -n 8

# Load the module
module load guppy/6.0.1

# Define input and output directories
input_dir="./raw_data/fast5"
output_dir="./raw_data/fastq"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Loop over all sample directories in the input directory
for sample_dir in "${input_dir}"/*/; do
    # Extract the sample name from the directory path
    sample=$(basename "${sample_dir}")

    # Define input fast5 file path for this sample
    input_fast5="${sample_dir}"

    # Run guppy_basecaller
    guppy_basecaller \
        --device "cuda:0" \
        --compress_fastq \
        --input_path "${input_fast5}" \
        --save_path "${output_dir}/${sample}" \
        --config dna_r9.4.1_450bps_hac.cfg \
        --num_callers 8 \
        --cpu_threads_per_caller 8 \
        --gpu_runners_per_device 8 \
        --chunk_size 1000 \
        --records_per_fastq 1000 \
        --recursive \
        --disable_pings \
        --num_read_chunks 8 \
        --qscore_filtering \
        --min_qscore 7 \
        --trim_barcodes
done
```
4. Load modules
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
5. Run fastqc on

i.      Two samples
```
fastqc \
        -t 4 \
        --extract \
        -o ./results/fastqc \
        -f fastq ./raw_data/Fastq/barcode02.all.fastq \
                 ./raw_data/Fastq/barcode48.all.fastq
```

ii.    Create the loop for all samples
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J fastqc
#SBATCH -n 4

# Set the input and output directories
INPUT_DIR=./raw_data/Fastq/
OUTPUT_DIR=./results/fastqc/

# Make directory to store the results
mkdir -p "$OUTPUT_DIR"

# Load FastQC module
module load fastqc/0.11.9

# Run FastQC on all ONT sequencing files in the directory "raw_data/merged_fastq_pass/"
for file in $INPUT_DIR/*.fastq; do
    fastqc -t 4 --extract -o $OUTPUT_DIR -f fastq $file
done
```
Save as run_fastqc and execute
```
sbatch -w compute05 run_fastqc
```
Copy to tour working directory and then to home directory
```
cp ./results/fastqc/barcode02.all_fastqc.html ~/
```
```
scp woguta@hpc.ilri.cgiar.org:*.all_fastqc.html .
```
```
explorer.exe .
```
6. Quality Control using NanoPlot
Install Nanoplpot
```
pip install NanoPlot
```
Upgrade
```
pip install NanoPlot --upgrade
```
or install the latest nanoplot version at once
```
pip3 install nanoplot
```
to update pip3
```
pip3 install --upgrade pip
```
check for the path and storage of the nanoplot
```
pip show nanoplot
```
If errors are observed/can't download, then install psyam first then nanoplot after
```
conda install -c bioconda pysam
```
```
conda install -c bioconda nanoplot
```
check the nanoplot 
```
conda list nanoplot
```
Run NnoPlot
i) Run two samples
```
# Add NanoPlot path to PATH variable
export PATH="$PATH:~/nanoplot"
NanoPlot \
    -t 4 \
    --fastq \
    ./raw_data/Fastq/barcode02.all.fastq \
    ./raw_data/Fastq/barcode48.all.fastq \
    -o ./results/nanoplot/ \
    --plots dot
```
ii) Loop for all samples
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J nanoplot
#SBATCH -n 4

# Set the input and output directories
INPUT_DIR=./raw_data/Fastq/
OUTPUT_DIR=./results/nanoplot/

# Make directory to store the results
mkdir -p "$OUTPUT_DIR"

# Add nanoplot path to PATH variable
export PATH="$PATH:~/nanoplot"

# Loop through all fastq files in input directory
for file in $INPUT_DIR/*.fastq; do
    # Get filename without extension
    filename=$(basename "$file" .fastq)
    # Run nanoplot with 4 threads and output to results directory
    NanoPlot --threads 4 --outdir $OUTPUT_DIR $file
    # Move generated plots to a subdirectory
    mkdir -p $OUTPUT_DIR/plots
    mv $OUTPUT_DIR/$filename* $OUTPUT_DIR/plots/
done
```
8. Adapters and Barcodes Trimmimg

ONT fastq files based called by guppyplex has been trimmed

Trimming is done to improve its quality. Some popular tools for trimming ONT data are:

Porechop: Porechop is a tool designed to trim adapters from ONT reads. It can also perform size selection and quality filtering.

NanoFilt: NanoFilt is a tool for filtering and trimming ONT reads based on quality scores, read length, and other parameters.

Filtlong: Filtlong is a tool for filtering and trimming long reads, including ONT reads, based on quality scores, length, and identity to a reference.

LongQC: LongQC is a tool that provides quality control metrics for long-read sequencing data, including ONT reads. It can also perform filtering and trimming based on read length and quality scores.

PycoQC: PycoQC is a tool that provides quality control metrics for ONT sequencing data. It can also trim reads based on quality scores.

Nanopack: Nanopack is a collection of tools for analyzing ONT sequencing data, including tools for quality control, filtering, and trimming.

Run porechop for all samples
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J porechop
#SBATCH -n 8

#load necessary modules
module purge
module load porechop/0.3.2pre

# Define input and output directories
input_dir="./raw_data/merged_fastq_pass"
output_dir="./results/porechop"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Loop over all fastq.gz files in the input directory
for fastq_gz_file in "${input_dir}"/*.all.fastq.gz;
do
    # Extract the filename without the extension
    filename=$(basename "${fastq_gz_file}")
    filename="${filename%.*}"

# Run Porechop to trim adapters and barcodes
    porechop -i "${fastq_gz_file}" \
             -o "${output_dir}/${filename}.trimmed.fastq.gz"
done
```
NB Renaming barcode*all.fastq to barcode.fastq
```
i=1
for barcode in *.all.fastq; do
    newname="barcode_${i}.fastq"
    mv "$barcode" "$newname"
    i=$((i+1))
done
```
or
```
for file in barcode*.all.fastq.trimmed.fastq.gz; do
  mv "$file" "${file%%.all.fastq.trimmed.fastq.gz}.trimmed.fastq.gz"
done
```
This script will rename all files in the current directory that end with .all.fastq to the format barcode_x.fastq, where x is a sequential number starting from 1.
9. Run Nanoplot again on the trimmed.fastq.gz files

10. Denovo genome assembly for all the samples, here we're using flye
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J flye
#SBATCH -n 8

#load required modules
module purge
module load flye/2.9

# Define input and output directories
input_dir="./results/porechop"
output_dir="./results/flye"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Loop over all sample directories in the input directory
for sample_dir in "${input_dir}"/*.fastq;
do
    # Extract the sample name from the directory path
    sample=$(basename "${sample_dir}" .fastq)

    # Make output directory for this sample if it does not exist
    mkdir -p "${output_dir}/${sample}"

    # Define input fastq file path for this sample
    input_fastq="${sample_dir}"

   #set genome size range
   if [ "$species" == "staphylococcus" ]; then
       genome_size="3m"
   elif [ "$species" == "enterobacter" ]; then
       genome_size="6m"
   elif [ "$species" == "e.coli" ]; then
    genome_size="6m"
   elif [ "$species" == "pseudomonas" ]; then
       genome_size="7m"
   elif [ "$species" == "default" ]; then
       genome_size="5m"
  fi

    # Run Flye with recommended parameters for bacterial denovo genome assembly
    flye --nano-corr "${input_fastq}" \
         --out-dir "${output_dir}/${sample}" \
         --genome-size "${genome_size}" \
         --threads 8 \
         --scaffold \
         --iteration 3 \
         --iterations 5 \
         --meta \
         --debug
done
```
flye: This is the command to run Flye, the assembler for long-read sequencing data.

--nano-raw "${input_fastq}": This option specifies the input FASTQ file to be assembled by Flye. The variable ${input_fastq} contains the path to the FASTQ file for the current barcode being processed.

--nano-raw "${input_fastq}": nano-corr signifies corected data here were guppyplex and porechop were used. In summary, nano-raw data is the raw signal data generated by the nanopore sequencer, while nano-corr data is the basecalled data that has been corrected for errors. Nano-raw data is typically used as input for basecalling, while nano-corr data is typically used as input for downstream analysis, such as genome assembly.

--out-dir "${output_dir}/${barcode}": This option specifies the output directory for Flye to write the assembly results to. The variable ${output_dir} contains the top-level directory where all the barcode-specific output directories will be created, and ${barcode} is the current barcode being processed.

--genome-size "3m-8m": This option specifies the expected size range of the genome being assembled. In this case, the expected size is between 3 and 8 million base pairs.

--threads 4: This option specifies the number of CPU threads to use during the assembly process. In this case, Flye will use 4 threads to perform the assembly.

--plasmids: This option enables Flye to identify and assemble plasmid sequences in the input data.

--scaffold: This option enables Flye to scaffold the contigs into larger sequences using paired-end reads.

--polish-target all: This option enables Flye to polish the assembled genome using all available data.

--iteration 3: This option specifies the iteration number for polishing the genome. In this case, Flye will perform polishing during the third iteration of the assembly process.

--iterations 4: This option specifies the maximum number of assembly iterations Flye will perform. In this case, Flye will perform up to 4 iterations to refine the assembly.

--debug: This option enables Flye to output debugging information to the console giving warning and troubleshooting issues with the assembly or fine-tuning the assembly parameters.

For PLasmid assembly use flye/2.8.1

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J plasmidflye
#SBATCH -n 5

#load required modules
module purge
module load flye/2.8.1

# Define input and output directories
input_dir="./results/porechop"
output_dir="./results/plasmidflye"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Loop over all sample directories in the input directory
for sample_dir in "${input_dir}"/*.fastq;
do
    # Extract the sample name from the directory path
    sample=$(basename "${sample_dir}" .fastq)

    # Make output directory for this sample if it does not exist
    mkdir -p "${output_dir}/${sample}"

    # Define input fastq file path for this sample
    input_fastq="${sample_dir}"

#set genome size range
   if [ "$species" == "staphylococcus" ]; then
       genome_size="3m"
   elif [ "$species" == "enterobacter" ]; then
       genome_size="6m"
   elif [ "$species" == "e.coli" ]; then
      genome_size="6m"
   elif [ "$species" == "pseudomonas" ]; then
      genome_size="7m"
   elif [ "$species" == "default" ]; then
      genome_size="5m"
  fi

  # Run Flye with recommended parameters for bacterial denovo genome assembly
    flye --plasmid \
         --nano-corr "${input_fastq}" \
         --out-dir "${output_dir}/${sample}" \
         --genome-size "${genome_size}" \
         --threads 5 \
         --iteration 3 \
         --iterations 4 \
         --meta \
         --debug
done
```
PLasmid assembly and Annotation using assembly.fasta or contigs.fasta if plasmid assembly uisng flye palsmid failed (Recommended)/Used here
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J Flair_plasmidID
#SBATCH -n 3

# Exit with error reporting
set -e

# Load modules
module purge
module load plasmidid/1.6.5

# Define directories and database path
flye_dir="./results/flye"
output_db="./plasmidIDdb"
dateNow=$(date +"%Y-%m-%d")
#dateNow="2023-06-29"

# Make the output directory if it doesn't exist
mkdir -p "${output_db}"

# Download PlasmidID database and set database path
download_plasmid_database.py -o "${output_db}"
plasmid_db_path="${output_db}/${dateNow}_plasmids.fasta"
#plasmid_db_path="./plasmidIDdb/2023-06-29_plasmids.fasta"

echo "Running PlasmidID"
# Loop through all assembly.fasta files in the input directory
for sample_dir in "${flye_dir}"/*; do
    if [ -d "$sample_dir" ]; then
        sample=$(basename "$sample_dir")

        # Define contigs path
        assembly_fasta="${sample_dir}/assembly.fasta"

        if [ -f "$assembly_fasta" ]; then
            output_path="${flye_dir}/${sample}"

            # Echo assembly file for the current sample
            echo "Assembly fasta for ${sample}: ${assembly_fasta}"

            # Run PlasmidID
            echo "Running PlasmidID for sample ${sample}"
            echo "Plasmid database path: ${plasmid_db_path}"
            plasmidID \
                -c "$assembly_fasta" \
                -d "${plasmid_db_path}" \
                -s "${sample}" \
                --no-trim \
                -T 3
        else
            echo "ERROR!!! contigs_file not found for sample: ${sample}"
        fi
    fi
done
```

If the script exits with error 1, make a new directory in which samples not meeting the threshold/giving errors are removed. Then run this script
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J Flair_plasmidID
#SBATCH -n 5

# Exit with error reporting
#set -e

# Load modules
module purge
module load plasmidid/1.6.5

# Define directories and database path
flye_dir="./results/flye1"
output_db="./plasmidIDdb"
#dateNow=$(date +"%Y-%m-%d")
#dateNow="2023-06-29"

# Make the output directory if it doesn't exist
mkdir -p "${output_db}"

# Download PlasmidID database and set database path
#download_plasmid_database.py -o "${output_db}"
#plasmid_db_path="${output_db}/${dateNow}_plasmids.fasta"
plasmid_db_path="./plasmidIDdb/2023-06-29_plasmids.fasta"

run_plasmidID() {
   trap ' echo Error $? occurred  on $LINENO && exit 1 ' ERR
   plasmidID \
        -c "$assembly_fasta" \
        -d "${plasmid_db_path}" \
        -s "${sample}" \
        --no-trim \
        -T 5
}

echo "Running PlasmidID"
# Loop through all assembly.fasta files in the input directory
for sample_dir in "${flye_dir}"/*; do
    if [ -d "$sample_dir" ]; then
        sample=$(basename "$sample_dir")

        # Define contigs path
        assembly_fasta="${sample_dir}/assembly.fasta"

        if [ -f "$assembly_fasta" ]; then
            output_path="${flye_dir}/${sample}"

            # Echo assembly file for the current sample
            echo "Assembly fasta for ${sample}: ${assembly_fasta}"

            # Check if plasmidID results already exist for ${sample}
            if [ -f "./NO_GROUP/${sample}/${sample}_final_results.tab" ]
            then
               echo "Results already exist for sample ${sample}. Skipping."
            continue
            else
               echo "Proceeding with analysis of ${sample}..."
            fi

             # Run PlasmidID
            echo "Running PlasmidID for sample ${sample}"
            echo "Plasmid database path: ${plasmid_db_path}"
            run_plasmidID
            if [ $? -eq 0 ]; then
                echo "PlasmidID analysis of ${sample} successful..."
            else
                echo -e "\tPlasmidID analysis of ${sample} was NOT successful...\n\tS$      continue
            fi
            else
                echo -e "ERROR!!! contigs_file not found for sample: ${sample}"
            exit
        fi
    fi
done
```
11. Species Identification using blast

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J blastn
#SBATCH -n 8

#module load necessary modules
module purge
module load blast/2.12.0+

# Set the input directory
input_dir="./results/flye"

# Set the output directory
output_dir="./results/blast"

# Set the BLAST database path
db_path="/export/data/bio/ncbi/blast/db/v5/nt"

# Make the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop over all samples in the input directory
for sample_dir in "${input_dir}"/*/; do
    # Extract the sample name from the directory path
    sample=$(basename "${sample_dir}")

# Define input fasta path for this sample
    input_fasta="${sample_dir}"/assembly.fasta

# Perform the BLASTN search
 blastn -task megablast \
        -query "${input_fasta}" \
        -db "${db_path}" \
        -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
        -culling_limit 5 \
        -num_threads 8 \
        -evalue 1e-25 \
        -out "${output_dir}/${sample}.assembly.vs.nt.cul5.1e25.megablast.out"
done
```
NB Cancelling a slurm job using scancel <job_ID>

```
scancel 3114
```
NB: Deleting all files within a given away from it, but be careful oh!
```
rm -rf ./path/to/the/file/folder/*
```
Extracting blast results
```
#remove all list saved in R in the terminal
rm(list = ls())

# Set the input directory where the sample directories are located
input_dir <- "./results/blast"

# Get a list of all file directories in the input directory
file_dirs <- list.dirs(input_dir, recursive = FALSE)

# Create a vector of available blast result files
blast_data_files <- list.files(path = input_dir, pattern = ".*\\.vs\\.nt\\.cul5\\.1e25\\.megablast\\.out$", full.names = TRUE)

# Create an empty data frame to store the aggregated data
agg_df <- data.frame()

# Loop over each blast result file
for (blast_file in blast_data_files) {
  # Read the blast data from the file
  blast_data <- read.delim(blast_file, sep = "\t", header = FALSE)

  # Skip file if no blast data was found
  if (nrow(blast_data) == 0) {
    next
  }

  # Extract the sample name from the file path
  sample_name <- tools::file_path_sans_ext(basename(blast_file))

  # Add the sample name as a column in the data frame
  blast_data$sample_name <- sample_name

  # Append the data to the aggregated data frame
  agg_df <- rbind(agg_df, blast_data)
}

# Set the output file name
output_file <- "./results/combined_blast_data.csv"

# Write the aggregated data to a CSV file
write.csv(agg_df, file = output_file, row.names = FALSE)
```
12. AMR Identification
i) Using RGI
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J rgi
#SBATCH -n 8

# Stop the script if any command fails
set -e

#load necessary modules
module purge
module load rgi/6.0.2

#Set the input and output directories
input_dir="./results/flye"
output_dir="./results/rgi"

#Make the output directory if it doesn't exist
mkdir -p "${output_dir}"

#Loop over all samples in the input directory
for sample_dir in "${input_dir}"/*/;
do
    # Extract the sample name from the directory path without extension
    sample=$(basename "${sample_dir}")

    #Make output directory for this sample if it doesn't exist
    mkdir -p "${output_dir}/${sample}"

    #Define input fasta path for this sample
    input_fasta="${sample_dir}/assembly.fasta"

    #Perform rgi analysis using contigs in assembly.fasta
    rgi main \
        -i "${input_fasta}" \
        -o "${output_dir}/${sample}/${sample}" \
        -t contig \
        --local \
        -a BLAST \
        -g PRODIGAL \
        --low_quality \
        --num_threads 8 \
        --split_prodigal_jobs \
        --clean \
        --debug 
done
```
OR

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J rgi
#SBATCH -n 8

#load necessary modules
module purge
module load rgi/6.0.2

#Set the input directory
input_dir="./results/flye"

#Set the output directory
output_dir="./results/rgi2"

#Make the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop over all sample directories in the input directory
for sample_dir in "${input_dir}"/*/;
do
    # Extract the sample name from the directory path
    sample=$(basename "${sample_dir}")

    # Make output directory for this sample if it doesn't exist
    mkdir -p "${output_dir}/${sample}"

    #Perform rgi analysis for contigs
    rgi main \
        -i "${sample_dir}/assembly.fasta" \
        -o "${output_dir}/${sample}/assembly.rgi.tsv" \
        -t contig \
        --local \
        -a BLAST \
        -g PRODIGAL \
        --low_quality \
        --num_threads 4 \
        --split_prodigal_jobs \
        --clean \
        --debug
   done
done
```
ii) Using Abricate
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J abricate
#SBATCH -n 8

# Load the required modules
module purge
module load abricate/1.0.1

# Set the input directory
input_dir="./results/flye"

# Set the output directory
output_dir="./results/abricate"

# Set the ABRICATE databases to use
databases=("ncbi" "card" "argannot" "resfinder" "megares" "plasmidfinder" "ecoli_vf" "ecoh" "vfdb")

# Set the ABRICATE and RGI minimum identity and coverage thresholds
min_identity="80"
min_coverage="60"

# Set the RGI database to use
rgi_database="CARD"

# Make the output directory if it doesn't exist
mkdir -p "${output_dir}"

# Loop over all sample directories in the input directory
for sample_dir in "${input_dir}"/*/; do
    # Extract the sample name from the directory path
    sample=$(basename "${sample_dir}")

    # Make output directory for this sample
    mkdir -p "${output_dir}/${sample}"

    # Run ABRICATE on the contigs file
    for db in "${databases[@]}"; do
        abricate --db "${db}" \
                 --minid "${min_identity}" \
                 --mincov "${min_coverage}" \
                 "${sample_dir}/assembly.fasta" \
                 > "${output_dir}/${sample}/contigs_${db}.abricate.tsv"
    done
done
```
Extracting AMR data Results

i) From RGI results
```
#!usr/bin/bashrc -l

#load neccessary modules
module purge
module load R/4.2

#open r console/studio in the terminal
R

#Path to the Rscript
#export PATH=".:$PATH"

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
   file_list <- list.files(path = sample_dir, pattern = paste0("contigs_",rgi_data_type,".json|txt"), full.names = TRUE)

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
#!usr/bin/bash -l

#load neccessary module
module load R/4.2
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
13. Variants Calling eg SNVs

i) Build snpEff config database for each bacterial species. wget the genome.fna  and gff files 

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J snpEff_db
#SBATCH -n 3

#Load modules
module purge
module load snpeff/4.1g
module load java/17

# Exit with error debug
set -e

#Define input directories
input_dir="./genome_ref"

#snpEff config for P. aeruginosa
echo "Processing snpEff config for P. aeruginosa"
mkdir -p ./database/snpEff/data/PAO1/
cp -rf $input_dir/P_aeruginosa/P_aureginosa_PA01_genomic.gff ./database/snpEff/data/PAO1/genes.gff
cp -rf $input_dir/P_aeruginosa/P_aureginosa_PA01_genomic.fna ./database/snpEff/data/PAO1/sequences.fa
echo -e "# P. aeruginosa bacterial genome, version PseudomonasPAO1\nPAO1.genome: PAO1" > ./database/snpEff/data/PAO1/snpEff.config

#Build the database
java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar build \
        -config ./database/snpEff/data/PAO1/snpEff.config \
        -dataDir ./../ \
        -gff3 \
        -v PAO1

#snpEff config S.aureus
echo "Processing snpEff config for S. aureus"
mkdir -p ./database/snpEff/data/NCTC8325/
cp -rf $input_dir/S_aureus/S_aureus_subsp_aureus_NCTC8325.gff ./database/snpEff/data/NCTC8325/genes.gff
cp -rf $input_dir/S_aureus/S_aureus_subsp_aureus_NCTC8325.fasta ./database/snpEff/data/NCTC8325/sequences.fa
echo -e "# S. aureus bacterial genome, version StaphylococcusNCTC8325\nNCTC8325.genome: NCTC8325" > ./database/snpEff/data/NCTC8325/snpEff.config

#Build the database
java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar build \
        -config ./database/snpEff/data/NCTC8325/snpEff.config \
        -dataDir ./../ \
        -gff3 \
        -v NCTC8325

#snpEff config S. sciuri
echo "Processing snpEff config for S. sciuri"
mkdir -p ./database/snpEff/data/ASM220916v2/
cp -rf $input_dir/S_sciuri/Staphylococcus_sciuri_ASM220916v2.gff ./database/snpEff/data/ASM220916v2/genes.gff
cp -rf $input_dir/S_sciuri/Staphylococcus_sciuri_ASM220916v2.fna ./database/snpEff/data/ASM220916v2/sequences.fa
echo -e "# S. sciuri bacterial genome, version StaphylococcusASM220916v2\nASM220916v2.genome: ASM220916v2" > ./database/snpEff/data/ASM220916v2/snpEff.config

#Build the database
java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar build \
        -config ./database/snpEff/data/ASM220916v2/snpEff.config \
        -dataDir ./../ \
        -gff3 \
        -v ASM220916v2

#snpEff config S. warneri
echo "Processing snpEff config for S. warneri"
mkdir -p ./database/snpEff/data/22.1/
cp -rf ./genome_ref/S_warneri/Staphylococcus_warneri_22.1.gff ./database/snpEff/data/22.1/genes.gff
cp -rf ./genome_ref/S_warneri/Staphylococcus_warneri_22.1.fna ./database/snpEff/data/22.1/sequences.fa
echo -e "# S. warneri bacterial genome, version Staphylococcus22.1\n22.1.genome: 22.1" > ./database/snpEff/data/22.1/snpEff.config

#Build the database
java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar build \
        -config ./database/snpEff/data/22.1/snpEff.config \
        -dataDir ./../ \
        -gff3 \
        -v 22.1

#snpEff config E. coli
echo "Processing snpEff config for E. coli"
mkdir -p ./database/snpEff/data/MG1655/
cp -rf $input_dir/E_coli/E_coli_str_K-12_substr_MG1655_genomic.gff ./database/snpEff/data/MG1655/genes.gff
cp -rf $input_dir/E_coli/E_coli_str_K-12_substr_MG1655_genomic.fasta ./database/snpEff/data/MG1655/sequences.fa
echo -e "# E. coli bacterial genome, version EscherichiaMG1655\nMG1655.genome: MG1655" > ./database/snpEff/data/MG1655/snpEff.config

#Build the database
java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar build \
        -config ./database/snpEff/data/MG1655/snpEff.config \
        -dataDir ./../ \
        -gff3 \
        -v MG1655

#snpEff config E. cloacae
echo "Processing snpEff config for E. cloacae"
mkdir -p ./database/snpEff/data/RS35/
cp -rf $input_dir/E_cloacae/Enterobacter_cloacae_strain_RS35.gff ./database/snpEff/data/RS35/genes.gff
cp -rf $input_dir/E_cloacae/Enterobacter_cloacae_strain_RS35.fna ./database/snpEff/data/RS35/sequences.fa
echo -e "# E. cloacae bacterial genome, version EnterobacterRS35\nRS35.genome: RS35" > ./database/snpEff/data/RS35/snpEff.config

#Build the database
java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar build \
        -config ./database/snpEff/data/RS35/snpEff.config \
        -dataDir ./../ \
        -gff3 \
        -v RS35
```
ii) Run snpEff for each spp

a) E .coli samples

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J Ecoli_snpEff
#SBATCH -n 8

# Exit immediately if a command exits with a non-zero status
set -e

# Load necessary modules
module purge
module load bwa/0.7.17
module load samtools/1.15.1
module load bcftools/1.15.1
module load snpeff/4.1g
module load java/17

# Specify the directories
input_dir="./results/porechop"
output_dir="./results/bwamem_ecoli"
output_dir2="./results/bcf_ecoli_variants"

# Define the path to the reference genome FASTA file
reference_file="./database/snpEff/data/MG1655/sequences.fa"

# Make output directories if not available
mkdir -p "${output_dir}"
mkdir -p "${output_dir2}"

# Define the output prefix for the index files and its directory
reference_prefix="${output_dir}/MG1655/MG1655"
mkdir -p "${reference_prefix}"

# Step 1: Index the reference genome
echo "Indexing the reference genome..."
bwa index -p "${reference_prefix}" -a bwtsw "${reference_file}"

echo "Processing E. coli samples"
#Define E. coli samples
ecoli_samples=(
  "${input_dir}/barcode_18.trimmed.fastq.gz"
  "${input_dir}/barcode_19.trimmed.fastq.gz"
  "${input_dir}/barcode_20.trimmed.fastq.gz"
  "${input_dir}/barcode_25.trimmed.fastq.gz"
  "${input_dir}/barcode_26.trimmed.fastq.gz"
  "${input_dir}/barcode_27.trimmed.fastq.gz"
  "${input_dir}/barcode_29.trimmed.fastq.gz"
  "${input_dir}/barcode_30.trimmed.fastq.gz"
  "${input_dir}/barcode_31.trimmed.fastq.gz"
  "${input_dir}/barcode_32.trimmed.fastq.gz"
  "${input_dir}/barcode_33.trimmed.fastq.gz"
  "${input_dir}/barcode_34.trimmed.fastq.gz"
  "${input_dir}/barcode_35.trimmed.fastq.gz"
  "${input_dir}/barcode_36.trimmed.fastq.gz"
  "${input_dir}/barcode_37.trimmed.fastq.gz"
  "${input_dir}/barcode_38.trimmed.fastq.gz"
  "${input_dir}/barcode_39.trimmed.fastq.gz"
  "${input_dir}/barcode_40.trimmed.fastq.gz"
  "${input_dir}/barcode_41.trimmed.fastq.gz"
  "${input_dir}/barcode_42.trimmed.fastq.gz"
  "${input_dir}/barcode_43.trimmed.fastq.gz"
  "${input_dir}/barcode_34.trimmed.fastq.gz"
)

echo "E. coli samples: ${ecoli_samples[@]}"
for sample in "${ecoli_samples[@]}"; do
  # Extract the sample name
  sample_name=$(basename "${sample}" .trimmed.fastq.gz)

  # Define input file path
  input_file="${sample}"

# Step 2: Perform alignment and generate sorted BAM files for each sample name
echo "Performing alignment..."
echo "Processing sample: ${sample_name}"
echo "Input file for ${sample_name}: ${input_file}"
#Perform alignment
  bwa mem -t 16 "${reference_prefix}" "${input_file}" |
  samtools view -bS - |
  samtools sort -@ 8 -o "${output_dir}/${sample_name}.sorted.bam" -
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
    bcftools call --threads 8 --ploidy 1 -Ou -mv -o "${output_dir2}/${sample_name}.vcf"
done

# Step 4: Filter and report the SNVs variants in variant calling format (VCF)
echo "Filtering variants..."
for vcf_file in "${output_dir2}"/*.vcf; do
    sample_name=$(basename "${vcf_file}" .vcf)
    echo "Processing ${sample_name}: ${vcf_file}"
    bcftools filter --threads 8 -i 'DP>=10' -Ov "${vcf_file}" > "${output_dir2}/${sample_name}.filtered.vcf"
done

#Step 5: Compress the filtered VCF files #bgzip -c file.vcf > file.vcf.gz
echo "Compressing filtered VCF files..."
for filtered_vcf in "${output_dir2}"/*.filtered.vcf; do
    echo "Compressing: ${filtered_vcf}"
    bgzip "${filtered_vcf}"
done

# Step 6: Index the filtered.vcf.gz VCF files
echo "Indexing filtered VCF files..."
for filtered_vcf_gz in "${output_dir2}"/*.filtered.vcf.gz; do
    echo "Indexing: ${filtered_vcf_gz}"
    bcftools index -c --threads 8 "${filtered_vcf_gz}"
done

# Step 7: Variant Annotation with GFF file
# Define input files and directories
vcf_dir="${output_dir2}"
annotated_dir="./results/annotated_ecoli_variants"
reference_file="./database/snpEff/data/MG1655/sequences.fa"
gff_file="../database/snpEff/data/MG1655/genes.gff"
snpeff_jar="/export/apps/snpeff/4.1g/snpEff.jar"

# Create output directory for annotated variants
mkdir -p "$annotated_dir"

# Run SnpEff for variant annotation #barcode_33.filtered.vcf.gz
echo "Performing variant annotation..."
for vcf_file in "$vcf_dir"/*.filtered.vcf.gz; do
    sample_name=$(basename "${vcf_file}" .filtered.vcf.gz)
    echo "Processing sample: $sample_name"
    bcftools view --threads 8 "$vcf_file" |
    java -Xmx4g -jar "$snpeff_jar" \
           -config ./database/snpEff/data/MG1655/snpEff.config \
           -dataDir ./../ \
           -v MG1655 ${vcf_file} > "${annotated_dir}/${sample_name}.snpEff.vcf"

echo "Variant annotation complete."

# Step 8: Rename summary.html and genes.txt and zip vcf files
mv ./snpEff_summary.html "${annotated_dir}/${sample_name}.snpEff.summary.html"
mv ./snpEff_genes.txt "${annotated_dir}/${sample_name}.snpEff.genes.txt"

# Compress vcf
bgzip -c "${annotated_dir}/${sample_name}.snpEff.vcf" > "${annotated_dir}/${sample_name}.snpEff.vcf.gz"

# Create tabix index - Samtools
tabix -p vcf -f "${annotated_dir}/${sample_name}.snpEff.vcf.gz"

# Generate VCF files
bcftools stats "${annotated_dir}/${sample_name}.snpEff.vcf.gz" > "${annotated_dir}/${sample_name}.snpEff.stats.txt"

echo "Variant summary renaming, compressing complete."

done

# Step 9: Variant Extraction with SnpSift
# Define input files and directories
snpeff_dir="${annotated_dir}"
extracted_dir="./results/extracted_variants"

# Create output directory for extracted variants
mkdir -p "$extracted_dir"

# Run SnpSift for variant extraction #barcode_1.filtered.snpEff.genes.txt
echo "Performing variant extraction..."
for snpeff_file in "$snpeff_dir"/*.snpEff.vcf.gz; do
    sample_name=$(basename "${snpeff_file}" .snpEff.vcf.gz)
    echo "Processing ${sample_name}: ${snpeff_file}"
    bcftools view --threads 5 "$snpeff_file" |
    java -Xmx4g -jar "/export/apps/snpeff/4.1g/SnpSift.jar" \
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

echo "SnpSift variant extraction complete"
```
b) Other samples here given by S.warneri

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J Flair_Pseudomonas_snpEff
#SBATCH -n 5

# Exit immediately if a command exits with a non-zero status
set -e

# Load necessary modules
module purge
module load bwa/0.7.17
module load samtools/1.15.1
module load bcftools/1.15.1
module load snpeff/4.1g
module load java/17

# Specify the directories
input_dir="./results/porechop/"
output_dir="./results/bwamem_warneri"
output_dir2="./results/bcf_warneri_variants"

# Define the path to the reference genome FASTA file
reference_file="./database/snpEff/data/22.1/sequences.fa"

# Make output directories if not available
mkdir -p "${output_dir}"
mkdir -p "${output_dir2}"

# Define the output prefix for the index files and its directory
reference_prefix="${output_dir}/22.1/22.1"
mkdir -p "${reference_prefix}"

# Step 1: Index the reference genome
echo "Indexing the reference genome..."
bwa index -p "${reference_prefix}" -a bwtsw "${reference_file}"

echo "Processing Pseudo samples"
for sample in "${input_dir}"/barcode_16.trimmed.fastq.gz; do
  #Extract the sample name
  sample_name=$(basename "${sample}" .trimmed.fastq.gz)

  #Define input file path
  input_file="${sample}"

  #Step 2: Perform alignment and generate sorted BAM files for each sample name
  echo "Performing alignment..."
  echo "Input file for ${sample_name}: ${input_file}"
  bwa mem -t 5 "${reference_prefix}" "${input_file}" |
  samtools view -bS - |
  samtools sort -@ 5 -o "${output_dir}/${sample_name}.sorted.bam" -
  samtools index "${output_dir}/${sample_name}.sorted.bam"
done

#Some stats
echo "Generating alignment statistics..."
for sorted_bam in "${output_dir}"/*.sorted.bam; do
  sample_name=$(basename "${sorted_bam}" .sorted.bam)
  samtools flagstat "${sorted_bam}"

  #Step 3: Perform variant calling on each sample
  echo "Proceeding with analysis of ${sample_name}..."
  echo "Variant calling for ${sample_name}: ${sorted_bam}"
  bcftools mpileup -Ou -f "${reference_file}" "${sorted_bam}" |
  bcftools call --threads 5 --ploidy 1 -Ou -mv -o "${output_dir2}/${sample_name}.vcf"
done

# Step 4: Filter and report the SNVs variants in variant calling format (VCF)
echo "Filtering variants..."
for vcf_file in "${output_dir2}"/*.vcf; do
    sample_name=$(basename "${vcf_file}" .vcf)
    echo "Processing sample ${sample_name}: ${vcf_file}"
    bcftools filter --threads 5 -i 'DP>=10' -Ov "${vcf_file}" > "${output_dir2}/${sample_name}.filtered.vcf"
done

#Step 5: Compress the filtered VCF files#bgzip -c file.vcf > file.vcf.gz
echo "Compressing filtered VCF files..."
for filtered_vcf in "${output_dir2}"/*.filtered.vcf; do
    echo "Compressing: ${filtered_vcf}"
    bgzip "${filtered_vcf}"
done

# Step 6: Index the filtered.vcf.gz VCF files
echo "Indexing filtered VCF files..."
for filtered_vcf_gz in "${output_dir2}"/*.filtered.vcf.gz; do
    echo "Indexing: ${filtered_vcf_gz}"
    bcftools index -c --threads 5 "${filtered_vcf_gz}"
done

#Step 7: Variant Annotation with GFF file
#Define input files and directories
vcf_dir="${output_dir2}"
annotated_dir="./results/annotated_pseudo_variants"
reference_file="./database/snpEff/data/22.1/sequences.fa"
gff_file="../database/snpEff/data/22.1/genes.gff"
snpeff_jar="/export/apps/snpeff/4.1g/snpEff.jar"

#Create output directory for annotated variants
mkdir -p "$annotated_dir"

#echo "Variant annotation for VCF files..."
for filtered_vcf_gz in "${output_dir2}"/*.filtered.vcf.gz; do
  sample_name=$(basename "${filtered_vcf_gz}" .filtered.vcf.gz)
  echo "Performing variant SnpEff annotation ${sample_name}: ${filtered_vcf_gz}"
  bcftools view --threads 5 "$filtered_vcf_gz" |
  java -Xmx4g -jar "$snpeff_jar" \
        -config ./database/snpEff/data/22.1/snpEff.config \
        -dataDir ./../ \
        -v 22.1 ${filtered_vcf_gz} > "${annotated_dir}/${sample_name}.snpEff.vcf"

#  echo "Variant annotation complete."

#Step 8: Rename summary.html and genes.txt and zip vcf files
mv ./snpEff_summary.html "${annotated_dir}/${sample_name}.snpEff.summary.html"
mv ./snpEff_genes.txt "${annotated_dir}/${sample_name}.snpEff.genes.txt"

# Compress vcf
bgzip -c "${annotated_dir}/${sample_name}.snpEff.vcf" > "${annotated_dir}/${sample_name}.snpEff.vcf.gz"

# Create tabix index - Samtools
tabix -p vcf -f "${annotated_dir}/${sample_name}.snpEff.vcf.gz"

# Generate VCF files
bcftools stats "${annotated_dir}/${sample_name}.snpEff.vcf.gz" > "${annotated_dir}/${sample_name}.snpEff.stats.txt"

echo "Variant summary, renaming and compressing complete for: ${staph_samples[@]}"
done

# Step 9: Variant Extraction with SnpSift
# Define input files and directories
snpeff_dir="${annotated_dir}"
extracted_dir="./results/extracted_variants"

# Create output directory for extracted variants
mkdir -p "$extracted_dir"

# Run SnpSift for variant extraction #barcode_1.filtered.snpEff.genes.txt
echo "Performing variant extraction..."
for snpeff_file in "$snpeff_dir"/*.snpEff.vcf.gz; do
    sample_name=$(basename "${snpeff_file}" .snpEff.vcf.gz)
    echo "Processing ${sample_name}: ${snpeff_file}"
    bcftools view --threads 5 "$snpeff_file" |
    java -Xmx4g -jar "/export/apps/snpeff/4.1g/SnpSift.jar" \
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

echo "SnpSift variant extraction complete"
```
c) Extract the snpSift text variant results by R

```
#!/usr/bin/Rscript

# Load necessary modules
module purge
module load R/4.3

# Open R console/studio in the terminal
R

# Remove all objects saved in R workspace
rm(list = ls())

# Set the input directory where the sample directories are located
input_dir <- "./results/extracted_variants"


# Get a list of all sample files in the input directory
file_list <- list.files(input_dir, full.names=TRUE)

# Create an empty data frame to store the aggregated data
agg_df <- data.frame()

  # Loop over each annotated result file
  for (file in file_list) {
    # Read the annotated data from the file
    annotated_data <- read.delim(file, sep = "\t", header = TRUE)

    # Skip sample if no annotated data was found
    if (nrow(annotated_data) == 0) {
      next
    }

    # Extract the sample name from the file path
    sample_name <- tools::file_path_sans_ext(basename(file))

    # Add the sample name as a column in the data frame
    annotated_data$sample_name <- sample_name

    # Append the data to the aggregated data frame
    agg_df <- rbind(agg_df, annotated_data)
  }

# Set the output file name
output_file <- paste0("combined_annotated_data.csv")

# Check if the aggregated data frame is empty
if (nrow(agg_df) > 0) {
  # Write the aggregated data to a CSV file
  write.csv(agg_df, file = output_file, row.names = FALSE)
  cat("Variant extraction complete. Data saved to:", output_file, "\n")
} else {
  cat("No data extracted. Output file not created.\n")
}
```

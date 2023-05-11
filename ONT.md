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
        -out "${output_dir}/${sample}contigs.vs.nt.cul5.1e25.megablast.out"
done
```
NB Cancelling a slurm job using scancel <job_ID>

```
scancel 3114
```

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

Use scp username@remote:/path/to/data /path/to/local/directory

```
scp -r woguta@hpc.ilri.cgiar.org:/path/to/data  .
```

The -r option is used to copy the directory recursively. The . at the end of the command specifies the current directory as the destination. space fullstop " ." specifies your current folder.

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
8. Trimmimg
Trimming is done to improve its quality. Some popular tools for trimming ONT data are:

Porechop: Porechop is a tool designed to trim adapters from ONT reads. It can also perform size selection and quality filtering.

NanoFilt: NanoFilt is a tool for filtering and trimming ONT reads based on quality scores, read length, and other parameters.

Filtlong: Filtlong is a tool for filtering and trimming long reads, including ONT reads, based on quality scores, length, and identity to a reference.

LongQC: LongQC is a tool that provides quality control metrics for long-read sequencing data, including ONT reads. It can also perform filtering and trimming based on read length and quality scores.

PycoQC: PycoQC is a tool that provides quality control metrics for ONT sequencing data. It can also trim reads based on quality scores.

Nanopack: Nanopack is a collection of tools for analyzing ONT sequencing data, including tools for quality control, filtering, and trimming.

i) Run two samples
```
# Add nanoplot path to PATH variable
export PATH="$PATH:~/nanoplot"
nanoplot \
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
    nanoplot --threads 4 --outdir $OUTPUT_DIR $file
    # Move generated plots to a subdirectory
    mkdir -p $OUTPUT_DIR/plots
    mv $OUTPUT_DIR/$filename* $OUTPUT_DIR/plots/
done
```
ONT fastq files based called by guppyplex has been trimmed
NB Renaming barcode*all.fastq to barcode.fastq
```
i=1
for barcode in *.all.fastq; do
    newname="barcode_${i}.fastq"
    mv "$barcode" "$newname"
    i=$((i+1))
done
```
This script will rename all files in the current directory that end with .all.fastq to the format barcode_x.fastq, where x is a sequential number starting from 1.

Denovo genome assembly for all the samples, here we're using flye
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J flye
#SBATCH -n 4

#module purge to avoid conflicts
module purge

#load required modules
module load flye/2.9

# Define input and output directories
input_dir="./raw_data/Fastq"
output_dir="./results/flye"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Loop over all barcode directories in the input directory
for barcode_file in "${input_dir}"/barcode_*.fastq;
do
    # Extract the barcode name from the directory path
    barcode=$(basename "${barcode_file}" .fastq)

    # Set read permission for the input file
    chmod +r "${barcode_file}"

    # Define input FASTQ file path for this barcode
    input_fastq="${barcode_file}"

    # Make output directory for this barcode if it does not exist
    mkdir -p "${output_dir}/${barcode}"

    # Run Flye with recommended parameters for bacterial denovo genome assembly
    flye --nano-raw "${input_fastq}" \
         --out-dir "${output_dir}/${barcode}" \
         --genome-size "3m-8m" \
         --threads 4 \
         --plasmids \
         --scaffold \
         --polish-target all \
         --iteration 3 \
         --iterations 4 \
         --debug
done
```
flye: This is the command to run Flye, the assembler for long-read sequencing data.

--nano-raw "${input_fastq}": This option specifies the input FASTQ file to be assembled by Flye. The variable ${input_fastq} contains the path to the FASTQ file for the current barcode being processed.

--out-dir "${output_dir}/${barcode}": This option specifies the output directory for Flye to write the assembly results to. The variable ${output_dir} contains the top-level directory where all the barcode-specific output directories will be created, and ${barcode} is the current barcode being processed.

--genome-size "3m-8m": This option specifies the expected size range of the genome being assembled. In this case, the expected size is between 3 and 8 million base pairs.

--threads 4: This option specifies the number of CPU threads to use during the assembly process. In this case, Flye will use 4 threads to perform the assembly.

--plasmids: This option enables Flye to identify and assemble plasmid sequences in the input data.

--scaffold: This option enables Flye to scaffold the contigs into larger sequences using paired-end reads.

--polish-target all: This option enables Flye to polish the assembled genome using all available data.

--iteration 3: This option specifies the iteration number for polishing the genome. In this case, Flye will perform polishing during the third iteration of the assembly process.

--iterations 4: This option specifies the maximum number of assembly iterations Flye will perform. In this case, Flye will perform up to 4 iterations to refine the assembly.

--debug: This option enables Flye to output debugging information to the console giving warning and troubleshooting issues with the assembly or fine-tuning the assembly parameters.

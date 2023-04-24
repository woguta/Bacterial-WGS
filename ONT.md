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
        -f fastq ./raw_data/merged_fastq_pass/barcode02.all.fastq \
                 ./raw_data/merged_fastq_pass/barcode48.all.fastq
```

ii.    Create the loop for all samples
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J fastqc
#SBATCH -n 4

# Set the input and output directories
INPUT_DIR=./raw_data/merged_fastq_pass/
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

If errors are observed/can't download, then install psyam first then nanoplot after
```
conda install -c bioconda pysam
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
nanoplot \
    -t 4 \
    --fastq \
    ./raw_data/merged_fastq_pass/barcode02.all.fastq \
    ./raw_data/merged_fastq_pass/barcode48.all.fastq \
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
INPUT_DIR=./raw_data/merged_fastq_pass/
OUTPUT_DIR=./results/nanoplot/

# Make directory to store the results
mkdir -p "$OUTPUT_DIR"

# Run nanoplot on all fastq files in the input directory
for file in $INPUT_DIR/*.fastq; do
    nanoplot \
        -t 4 \
        --fastq \
        $file \
        -o $OUTPUT_DIR \
        --plots dot
done
```

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
it keeps the .gz files too, to remove/delete .gz files run within your current directory 
```
rm -f *.gz
```
alternatively, run this code once without k
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

b) Rerun the fastqc on the trimmed one file
```
fastqc
        -t 4
        -o ./results/fastqc/
        -f fastq ./raw_data/Fastq/AS-26335-C1-C_S4_L001_R1_001.trim.fastq
                ./raw_data/Fastq/AS-26335-C1-C_S4_L001_R2_001.trim.fastq
```

b) Fastqc loop shell in bash for all the the files

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
To use this script, you can save it to a file (e.g., run_fastqc.sh), make it executable (chmod +x run_fastqc.sh), and then run it with the directory paths as arguments
```
bash run_fastqc.sh ./results/fastqc/ ./raw_data/Fastq/
```

#Running fastp for all the files/trimming all the files
i) First Pathway
```
#!/bin/bash

#Set the input and output directories
INPUT_DIR=./raw_data/Fastq/
OUTPUT_DIR=./results/fastp/

#load modules
module load fastp/0.22.0

#Loop through all files in the input directory
for FILENAME in ${INPUT_DIR}*.fastq.gz; do

  #Extract the sample name from the file name
  SAMPLE_NAME=$(basename ${FILENAME} .fastq.gz)

  #Run fastp on the sample
  fastp \
    --in1 ${INPUT_DIR}${SAMPLE_NAME}.fastq.gz \
    --in2 ${INPUT_DIR}${SAMPLE_NAME}_2.fastq.gz \
    --out1 ${OUTPUT_DIR}${SAMPLE_NAME}.trim.fastq.gz \
    --out2 ${OUTPUT_DIR}${SAMPLE_NAME}_2.trim.fastq.gz \
    --json ${OUTPUT_DIR}${SAMPLE_NAME}.fastp.json \
    --html ${OUTPUT_DIR}${SAMPLE_NAME}.fastp.html \
    --failed_out ${OUTPUT_DIR}${SAMPLE_NAME}_failed.fastq.gz \
    --thread 4 \
    -5 -3 -r \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --cut_mean_quality 20 \
    --length_required 15 \
    --dedup \
    |& tee ${OUTPUT_DIR}${SAMPLE_NAME}.fastp.log

done
```
Save the scrip as run_fastp.sh

run the script
```
bash run_fastp.sh
```
ii) Second Pathway
```
#making directories
INPUT_DIR=./raw_data/Fastq/
OUTPUT_DIR=./results/fastp/

#load modules
module load fastp/0.22.0

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
run_fastp2.sh
```
iii) Third pathway sbatch - sunmitiing jobs to the cluster as you do other things (the best option)
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
mkdir -p $OUTPUT_DIR

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

sequence/contigs
```
ls -lht ./results/spades/contigs.fasta
```
First 10 files/head
```
grep '>' contigs.fasta | head
```
11.Genome Assessment 

i) Genome contiguity

Checks length/cutoff for the longest contigs that contain 50% of the total genome length measured as contig N50. i.e  involves evaluating the accuracy and completeness of the genome assembly using metrics such as N50 length, scaffold and contig numbers, and genome size. Tool used QUAST.

a) For one sample
```
quast.py \
./results/spades/contigs.fasta \
-t 4 \
-o ./results/quast
```

Inspect the quast report

i) to your local wkd from hpc
```
cp -i AS-27566-C1-C_S23_L001/*.html ~/
```
ii) To local computer homepage/wkd
```
scp woguta@hpc.ilri.cgiar.org:~/AS-27566-C1-C_S23_L001.html .
```
b) For all the samples
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J QUAST
#SBATCH -n 4

# Set the path to the directory containing the input files
input_dir=./results/spades

# Set the path to the directory where the output will be saved
output_dir=./results/quast

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
 
 ii) Genome completeness
 
Assesses the presence or absence of highly conserved genes (orthologs) in an assembly/ ensures that all regions of the genome have been sequenced and assembled. It's performed using BUSCO (Benchmarking Universal Single-Copy Orthologs). Ideally, the sequenced genome should contain most of these highly conserved genes. they're lacking or less then the genome is not complete.

a) For one sample 

i) Using genome ref
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

b) Loop for all files?
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J BUSCO
#SBATCH -n 4

# Set the path to the directory containing the input files
input_dir=./results/spades

# Set the path to the directory where the output will be saved
output_dir=./results/busco

# Set the path to the BUSCO lineage-specific database
database=bacteria_odb10

# Loop through all FASTA files in the input directory
for file in ${input_dir}/*.fasta; do
    filename=$(basename "$file")
    output_path="${output_dir}/${filename%.*}"
    
    # Run the BUSCO command on the current file
    busco -i "$file" -m bacteria -o "$output_path" -l "$database" -c 4 -f
done
```
Save us run_busco.sh
```
 sbatch -w compute05 run_busco.sh
 ```

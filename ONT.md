# The ONT WGS (Whole Genome Sequencing) bacterial pipeline

Is a bioinformatics pipeline designed for analyzing bacterial genomes using Oxford Nanopore Technologies (ONT) sequencing data. The pipeline includes several steps, which are outlined below:

1.  Basecalling: Raw ONT sequencing data (in the form of electrical signals) are converted into nucleotide sequences using basecalling algorithms such as Guppy.

2.  Quality control: The quality of the basecalled sequences is assessed using tools such as FastQC and NanoPlot. Sequences with poor quality are filtered out.

3.  Assembly: The filtered sequences are assembled into contigs using a variety of tools such as Flye, Canu, and Shasta.

4.  Polishing: The assembled contigs are polished to improve their accuracy and reduce errors using tools such as Racon and Medaka.

5.  Scaffold: The polished contigs are scaffolded to link the contigs and produce a more complete genome sequence using tools such as SSPACE-LongRead.

6.  Gap Filling: Gaps between contigs are filled using tools such as Pilon.

7.  Annotation: The assembled and polished genome sequence is annotated to identify genes, functional elements, and other features using tools such as Prokka and RAST.

8.  Comparative analysis: The annotated genomes are compared to reference genomes and other bacterial genomes to identify similarities, differences, and evolutionary relationships using tools such as Mauve and Roary.

9.  Visualization: The results of the analysis are visualized using tools such as Artemis and Circos.

# Steps

# Log into hpc via ssh & interactive computing
```
interactive -w compute05 -c 3
```

# 1.  Creating directory
```
cd /var/scratch/
mkdir -p $USER/bacteria-wgs/flair
cd $USER/bacteria-wgs/flair
```

The mkdir -p command is used to create a directory and its parent directories (if they do not exist) in a single command.

The -p option allows the mkdir command to create the parent directories if they do not exist.

# 2.  Data retieval

Use scp username@remote:/path/to/data /path/to/local/directory

```
scp -r woguta@hpc.ilri.cgiar.org:/path/to/data  .
```

The -r option is used to copy the directory recursively. The . at the end of the command specifies the current directory as the destination. space fullstop " ." specifies your current folder.


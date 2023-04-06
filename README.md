# Bacterial-WGS
## Analysis of ONT and Illumina WGS data
## A. Illumina
1. Creating directory
```
cd /var/scratch/
mkdir -p $USER/bacteria-wgs/crpa
cd $USER/bacteria-wgs/crpa

```
2. Data retieval
use scp username@remote:/path/to/data /path/to/local/directory
```
scp -r woguta@hpc.ilri.cgiar.org:/export/data/ilri/miseq/users/Sam/BacteriaWGS_WO/Fastq/ .

```
The -r option is used to copy the directory recursively. The . at the end of the command specifies the current directory as the destination.
3. Viewing  the files
 cd into current directory and view all the retrieved data by run
 ```
run less filename.fastq
run head -n 10 filename.fastq
run tail -n 10 filename.fastq
```
4. Unzip .gz files
Use the gunzip or gzip command depending on your system
```
gunzip -k *.gz
```
it keeps the .gz files too, to remove/delete .gz files run within your current directory 
```
rm -f *.gz
```
alternative run this code once without k
```
gunzip *.gz
```
View the files as in 3 above


# Steps

## Log into hpc via ssh & interactive computing
```
interactive -w compute05 -c 3
```

## 1.  Creating directory
```
cd /var/scratch/
mkdir -p $USER/bacteria-wgs/flair
cd $USER/bacteria-wgs/flair
```

The mkdir -p command is used to create a directory and its parent directories (if they do not exist) in a single command.

The -p option allows the mkdir command to create the parent directories if they do not exist.

## 2.  Data retieval

Use scp username@remote:/path/to/data /path/to/local/directory
```
scp -r woguta@hpc.ilri.cgiar.org:/path/to/data  .
```

The -r option is used to copy the directory recursively. The . at the end of the command specifies the current directory as the destination. space fullstop " ." specifies your current folder.


# fqvendorfail
Drop casava vendor failed fastq reads

# Install

1. Clone the repo
2. `make`

This tool needs you to have zlib installed and uses the provided version of
Heng Li's seqtk [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk).

# Usage

```
$ ./fqvendorfail

Usage: fqvendorfail -o <output prefix> <R1.fq> [<R2.fq>]
    -o STRING output prefix for filtered fastq files.
```

Provide an output prefix (e.g., `/path/my_file`) and 1 fastq file if single-end
or two if paired-end. The tool will add the `_R1.fq.gz` and `_R2.fq.gz` to the 
output prefix you provide.

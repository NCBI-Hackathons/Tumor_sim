# Tumor_sim
Simulation of Tumor Genomes -- Initiated at the 2017 NYGC-NCBI Hackathon

Goal: Generate a simulated tumor genome, based on a simulated normal genome BAM file created from a user-provided genome file as reference.

The software package is written in Python and contained in the [lib folder](https://github.com/NCBI-Hackathons/Tumor_sim/tree/master/lib). 

Dependencies: found in requirements.txt
***
## Installation
```
wget ...
```

## Usage
```

```

## Test
To download hg38: wget  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

To download hg19, use the following commands:

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -xzf chromFa.tar.gz
cat chr*.fa > hg19.fa
gzip hg19.fa

A subsampled version of hg38 is also provided in the [data folder](https://github.com/NCBI-Hackathons/Tumor_sim/tree/master/data) of this repository

### Resource and references:

[![Travis](https://api.travis-ci.org/NCBI-Hackathons/Tumor_sim.svg?branch=master)](https://travis-ci.org/NCBI-Hackathons/Tumor_sim)

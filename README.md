# Tumor_sim
A python package that simulates the structural variations of cancer genomes. -- Initiated at the 2017 NYGC-NCBI Hackathon 

### Background 
There is no "ground truth" for detecting structural variations in cancer genomics. Complex rearrangements occur per individual cancer samples, with no distinction in the sequence of mutations (SNPs, copy number variations and indels to CNVs, transpositions, and duplication) that lead to any specific structural arrangement. Simulating the possible SNVs and CNVs in a unique variation of a reference genome followed by simulation tumor model of cancer-related variations will be beneficial in
better understanding the pathways of cancer progression after a mutation event, as found in cell division, aneuploidy, and other relative phenomena.

### Goal
Generate a simulated tumor genome, based on a user-provided genome file as reference. 
***
## Installation
*Dependencies:* found in requirements.txt, please be sure to have them installed on your system. The software package is written in Python and contained in the [lib folder](https://github.com/NCBI-Hackathons/Tumor_sim/tree/master/lib). 
```
> pip install -r requirements.txt
```

## Usage
The user can provide their reference genome in FASTA format as input file to `simulate_endToEnd.py`, or use our default example fasta.

```
> cd lib
> python simulate_endToEnd.py 
> or (optional)
> python simulate_endToEnd.py [-usage] <path/to/input_file>
```

To view the help options, type -h:
```
> python simuate_endToEnd.py -h
Usage: simulate_endToEnd.py [-h] [--input_fasta INPUT_FASTA]
                            [--output_tumor_fasta OUTPUT_TUMOR_FASTA]
                            [--output_normal_fasta OUTPUT_NORMAL_FASTA]

Simulate cancer genomic structural variations

optional arguments:
  -h, --help            show this help message and exit
  --input_fasta INPUT_FASTA
                        file path for the input (default genome) fasta
  --output_tumor_fasta OUTPUT_TUMOR_FASTA
                        file path for the output tumor (cancer genome) fasta
  --output_normal_fasta OUTPUT_NORMAL_FASTA
                        file path for the output normal (SNV-added) fasta
```

### How it works
simulate_endToEnd.py calls upon mutation_orchestra.py to generate the random mutations for the normal unique genome case.
mutation_orchestra.py uses the lower level mutation_creator.py for the simpler mutations and logs the distinction of indels, translocations, duplications and inversions.

## Test

To run unit tests, run `nosetests` from the top-lelel directory of the project. If you encounter errors, such as test_main failing because it cannot find the data in the data folder, it's most likely because you are not running tests from the top-level directory.
A subsampled version of hg38 is also provided in the [data folder](https://github.com/NCBI-Hackathons/Tumor_sim/tree/master/data) of this repository.
To download the reference FASTA hg38 or hg19, use the following commands:

For reference fasta `hg38.fa` (approximately 3.0 GB):

```
> wget  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
> tar -xzf hg38.chromFa.tar.gz
> cd chroms
> cat chr*.fa > hg38.fa
> rm chr*.fa
```
For reference fasta `hg19.fa` (approximately 3.1 GB):

```
> wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
> tar -xzf chromFa.tar.gz
> cat chr*.fa > hg19.fa
> rm chr*.fa
```

[![Travis](https://api.travis-ci.org/NCBI-Hackathons/Tumor_sim.svg?branch=master)](https://travis-ci.org/NCBI-Hackathons/Tumor_sim)

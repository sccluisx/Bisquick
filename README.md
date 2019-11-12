# Bisquick
## Kmer based algorithm to calculate methylation rates from dna sequence reads


### Introduction
Bisquick is an algorithm to calculate methylation rates on a given genome and a reads file

### Instalation
#### Requirements
* Cmake
* Modern C++ compiler
* [Seqan Library](https://seqan.readthedocs.io/en/seqan-v2.4.0/ "Seqan") v. 2.4 
* [Boost Library](https://www.boost.org/)

Bisquick is designed for linux64 platform. 


1. Clone this repository
2. Create a Release directory `mkdir Release`
3. Move to the Release directory `cd Release`
4. Execute Cmake `cmake ..`
5. Compile the code `make`

### Usage

#### SYNOPSIS
    bisquick [OPTIONS]

#### DESCRIPTION
    Bisulfite analysis as fast as baking a pie!

#### OPTIONS
    -h, --help
          Display the help message.
    --version
          Display version information.
    -k, --kmersize INTEGER
          Indicate the length of the kmer. Default: 7.
    -gp, --genomePath STRING
          Directory path of the genome directory.
    -r, --readsfile STRING
          File of the reads.
    -o, --output STRING
          Output directory.

#### VERSION
    Last update: November 2019
    bisquick version: 0.1
    SeqAn version: 2.4.0
    
### Usage Example
For this example we will use the data in the testdata folder

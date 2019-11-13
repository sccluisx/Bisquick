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

For this example we will use the data in the [testdata](https://github.com/sccluisx/Bisquick/tree/master/testdata testdata) folder which contains the dna sequence reference of the human chromosome 21.
For this example we will simulate Bisulfite reads using  the tool [Sherman](https://github.com/FelixKrueger/Sherman), download it from the repo and enter the next command

`./Sherman -l 300 -n 1000 -g testdata/chrm21/ref -o testdata/chrm21/bisulfite_reads -CG 40 -CH 0 -e 0 `

this command will genearte 1000 reads  of 300 nucleotides long each one, and will introduce a Bisulfite conversion rate for cytosines in CG-context of 40% 

Now we can use bisquick
`bisquick -k 50 -gp testdata/chrm21/ref -r testdata/chrm21/bisulfite_reads/simulated.fastq -o test > test_sumary.txt`

In this case we are running bisquick with a 50-mer as parameter and bisquick will generate two outputs 

* test_sumary.txt: A summary report of the program execution. Most important will report the calculated methylation rate between other things like, execution time of building the index, execution time of processing the reads, size of the compressed kmer hashtable, etc...

* a BED file containing the methylation rate of every CpG position of the genome reference following the format

`<chromosome>  <start position> <end position> <methylation percentage>`

This BEDfile can be used to produce a bedgraph in the genome browser of the UCSC and visualize it in a [custom track](https://genome.ucsc.edu/cgi-bin/hgCustom "genome browser")

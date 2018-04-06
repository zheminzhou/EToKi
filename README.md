# EToKi (Enterobase Tool Kit)
all methods related to Enterobase

## External dependencies:
#### phylo module
* raxml

#### prepare module
* bbmap

#### MLSType module
* blastn
* usearch

#### assemble module
* samtools
* bwa
* bowtie2
* kraken
* gatk
* megahit
* pilon
* spades

#### ortho module
* mcl
* fasttree
* vsearch

## Python version: 2.7.9

## PIP Packages:
* ete3==3.0.0b17
* hashlib==20081119
* numpy==1.11.3
* pandas==0.18.0
* psutil==5.1.0

## Installation: 
Assume all the external dependencies are in system PATH.
```
git clone https://github.com/zheminzhou/EToKi.git
cd EToKi
pip install -r requirements.txt
python EToKi.py configure
```
Specify the links to external commands if they are not in the system PATH. Use
```
python EToKi.py configure -h
```
For additional information.

## Usage:   EToKi.py <command> [options]

```
Commands:
  configure        Configure external dependencies
  prepare          Preprocessing for short reads
  assemble         de novo / reference-guided asembly for either metagenomic or non-metagenomic reads
  ortho            Pan-genome prediction
  MLSTdb           Create database for MLST typing
  MLSType          MLST nomenclature
  phylo            Infer phylogeny and ancestral states from genomic alignments or SNP matrix
  RecHMM           Identify Recombination sketches.


Use EToKi.py <command> -h for help in each command.
```

## Examples: 

### 1. phylogeny + ancestral reconstruction + recombination detection
```
cd examples
python ../EToKi.py phylo -t all -p phylo_out -m phylo_rec.fasta
python ../EToKi.py RecHMM -d phylo_out.mutations.gz -p rec_out
```

Outputs are:

* examples/phylo_out.matrix.gz - a table of mutations in the alignment. 
* examples/phylo_out.labeled.nwk - phylogeny with labeled internal nodes. 
* examples/phylo_out.ancestral_states.gz - Ancestral states of internal nodes. 
* examples/phylo_out.mutations.gz - Occurences of mutations on different branches. 
* examples/rec_out.best.model.report - Estimated parameters for the recombinations. 
* examples/rec_out.recombination.region - Identified recombination regions. 

#### NOTE: New RecHMM identifies three categories of recombinations:
1. External: Recombination with an external source. High SNP densities and low homoplasies. 
2. Internal: Recombination from an internal source. Normal SNP densities and high homoplasies. 
3. Mixed:    Repetitive imports from external sources to different branches of the tree. High SNP densities and high homoplasies. 

You can use legacy RecHMM by setting --task 0, which considers only external source. 
### 2. short reads preprocess + assembly + polish + consensus quality + evaluation
```
cd examples
python ../EToKi.py prepare --pe A_R1.fastq.gz,A_R2.fastq.gz -p prep_out
python ../EToKi.py assemble --pe prep_out_L1_R1.fastq.gz,prep_out_L1_R2.fastq.gz --se prep_out_L1_R3.fastq.gz -p asm_out
```

Outputs are:

* examples/prep_out_L1_R?.fastq.gz - Short reads after quality filtering
* asm_out.result.fastq - Assembly with consensus quality information
* asm_out.result.fasta - Assembly in fasta format

#### NOTE: example is copied from sample data in SPAdes 3.10


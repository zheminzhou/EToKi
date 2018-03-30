# EToKi (Enterobase Tool Kit)
all methods related to Enterobase

External dependencies:
* raxml
* bbmap
* usearch
* samtools
* kraken
* gatk
* megahit
* mcl
* pilon
* fasttree
* bwa
* bowite2
* spades
* blastn
* vsearch


Python version: 2.7.9

Packages:
ete3==3.0.0b17
hashlib==20081119
numpy==1.11.3
pandas==0.18.0
psutil==5.1.0

Installation: 
> git clone 
> EToKi.py EnConf

Usage:   EToKi.py <command> [options]

Commands:
  EnConf            configure external dependencies
  EnPrep            Preprocessing for short reads
  EnBler            de novo / reference-guided asembly for either metagenomic or non-metagenomic reads
  EnOrth            Pan-genome prediction
  EnServ            Create database for MLST typing
  EnSign            MLST nomenclature
  EnPhyl            Infer phylogeny and ancestral states from genomic alignments or SNP matrix
  RecHMM            Identify Recombination sketches.

Use EToKi.py <command> -h for help in different commands.



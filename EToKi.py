import argparse, numpy as np, sys

def etoki() :
    try :
        if sys.argv[1] == 'RecHMM' :
            sys.argv[1] = 'EnRaRe'
       
        exec 'from {0} import {0}'.format(sys.argv[1])
    except :
        print '''
Program: EToKi (Enterobase Tool Kit)

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

Use EToKi.py <command> -h to get help for each command.
'''
        sys.exit(0)
    eval(sys.argv[1])(sys.argv[2:])


if __name__ == '__main__' :
    etoki()

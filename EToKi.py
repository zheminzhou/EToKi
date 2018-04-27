import numpy as np, sys

def etoki() :
    try :
        exec 'from modules.{0} import {0}'.format(sys.argv[1])
        sys.argv[0] = ' '.join(sys.argv[:2])
    except Exception as e :
        print '''
Program: EToKi (Enterobase Tool Kit)

Usage:   EToKi.py <command> [options]

Commands:
  configure        Configure external dependencies
  prepare          Preprocessing for short reads
  assemble         de novo / reference-guided asembly for either metagenomic or non-metagenomic reads
  ortho            Pan-genome prediction
  MLSTdb           Create database for MLST typing
  MLSType          MLST nomenclature
  phylo            Infer phylogeny and ancestral states from genomic alignments or SNP matrix
  RecHMM           Identify Recombination sketches
  RecFilter        Remove Recombination sketches
  BrRefine         Correct tree using RecHMM outpus

Use EToKi.py <command> -h to get help for each command.
'''
        if len(sys.argv) > 1 :
            import traceback
            traceback.print_exception(*sys.exc_info())
        sys.exit(0)

    eval(sys.argv[1])(sys.argv[2:])


if __name__ == '__main__' :
    etoki()

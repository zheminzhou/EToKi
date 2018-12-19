import sys


def etoki():
    try:
        exec('from modules.{0} import {0}'.format(sys.argv[1]))
        sys.argv[0] = ' '.join(sys.argv[:2])
    except :
        print('''
Program: EToKi (Enterobase Tool Kit)

Usage:   EToKi.py <command> [options]

Commands:
  0)configure        Configure external dependencies
  1)prepare          Preprocessing for short reads
  10)assemble        de novo / reference-guided asembly for either metagenomic or non-metagenomic reads
  11)ortho           Pan-genome prediction (prepare for wgMLST scheme)
  12)MLSTdb          Create database for MLST typing
  13)MLSType         MLST nomenclature
  14)cgMLST          evaluate wgMLST genes using a reference set of genomes (working)
  15)hierCC          hierarchical cgMLST clustering
  20)align           align genomes onto a reference
  21)toVCF           combine multiple alignments into matrix (working)
  22)phylo           Infer phylogeny and ancestral states from genomic alignments or SNP matrix
  23)RecHMM          Identify Recombination sketches

Use EToKi.py <command> -h to get help for each command.
''')
        if len(sys.argv) > 1:
            import traceback
            traceback.print_exception(*sys.exc_info())
        sys.exit(0)

    eval(sys.argv[1])(sys.argv[2:])


if __name__ == '__main__' :
    etoki()

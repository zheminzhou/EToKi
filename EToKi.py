import sys, os
from modules.configure import logger
#refMasker, profile2matrix
#            ['metaCaller','secondary SNP caller for metagenomic samples with intra-species diversity'], 

commands = [['configure', 'configure external dependencies if you do not use implemented version.'],    #
            ['prepare',   'trim, collapse, downsize and rename in short reads.'],                       #
            ['assemble',  'de novo or reference-guided assembly for genomic or metagenomic reads'],     #
            ['ortho',     'Pan-genome prediction (also used to find genes for wgMLST schemes)'],        #
            ['MLSTdb',    'select exemplar genes from pan genes as references for MLST typing'],        #
            ['MLSType',   'MLST nomenclature using a local set of references'],                         #
            ['cgMLST',    'pick core genes from wgMLST genes'], 
            ['hierCC',    'generate hierarchical clusters from cgMLST profiles'],                       #
            ['evalHCC',   'evaluate HierCC to find the most stable clusters'], 
            ['align',     'align multiple queried genomes to a single reference'], 
            ['phylo',     'infer phylogeny and ancestral states from genomic alignments or SNP matrix'], 
            ['RecHMM',    'identify recombination sketches from a SNP matrix'], 
            ['EBEis',     'ab initio serotype prediction for Escherichia'],                             #
            ['uberBlast', 'A merged BLAST-like results from Blastn, uBlastp, minimap2 and mmseqs'],     #
            ['clust',     'cluster of short sequences using mmseqs linclust']]                          #


def etoki():
    try:
        if len(sys.argv) <= 1 :
            raise ValueError
        try :
            exec('from modules.{0} import {0}'.format(sys.argv[1]))
        except ImportError as e :
            logger(str(e))
            raise ValueError
        else:
            sys.argv[0] = ' '.join(sys.argv[:2])
            
    except ValueError as e :
        sys.stdout.write('''
Program: EToKi (Enterobase Tool Kit)

Usage:   EToKi.py <command> [options]

Commands:
'''
              + '\n'.join(['    {0} {1}'.format(cmd[0].ljust(12), cmd[1]) for cmd in commands]) + 
'''
Use EToKi.py <command> -h to get help for each command.
''')
    else :
        eval(sys.argv[1])(sys.argv[2:])


if __name__ == '__main__' :
    etoki()

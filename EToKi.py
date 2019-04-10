import sys, os
from modules.configure import logger
#refMasker, profile2matrix
#            ['metaCaller','secondary SNP caller for metagenomic samples with intra-species diversity'], 

commands = [['configure', 'install and/or configure 3rd party programs'],    #
            ['prepare',   'trim, collapse, downsize and rename the short reads.'],                       #
            ['assemble',  'de novo or reference-guided assembly for genomic or metagenomic reads'],     #
            ['ortho',     'pan-genome (and wgMLST scheme) prediction'],        #
            ['MLSTdb',    'Set up exemplar alleles and database for MLST schemes'],        #
            ['MLSType',   'MLST nomenclature using a local set of references'],                         #
#            ['cgMLST',    'pick core genes from wgMLST genes'], 
            ['hierCC',    'generate hierarchical clusters from cgMLST profiles'],                       #
#            ['evalHCC',   'evaluate HierCC to find the most stable clusters'], 
            ['align',     'align multiple queried genomes to a single reference'], 
            ['phylo',     'infer phylogeny and ancestral states from genomic alignments'], 
            ['RecHMM',    'identify recombination sketches from a SNP matrix'], 
            ['RecFilter', 'Remove recombination sketches from a SNP matrix'], 
            ['EBEis',     'ab initio serotype prediction for Escherichia coli and Shigella spp.'],                             #
            ['isCRISPOL', 'ab initio prediction of CRISPOL array for Salmonella enterica serovar Typhimurium'],                #
            ['uberBlast', 'Use Blastn, uBlastp, minimap2 and/or mmseqs to identify similar sequences'],     #
            ['clust',     'linear-time clustering of short sequences using mmseqs linclust']]                          #


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

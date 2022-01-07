#! /usr/bin/python3
import sys, os
import argparse
from modules.configure import logger

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


commands = [['configure', 'install and/or configure 3rd party programs'],    #
            ['prepare',   'trim, collapse, downsize and rename the short reads.'],                       #
            ['assemble',  'de novo or reference-guided assembly for genomic or metagenomic reads'],     #
            ['MLSTdb',    'Set up exemplar alleles and database for MLST schemes'],        #
            ['MLSType',   'MLST nomenclature using a local set of references'],                         #
            ['MLSTsum',   'Summarise MLSType results and assign new allele designations'],
            ['cgMLST',    'Select a list of genes for the cgMLST scheme'],
            ['align',     'align multiple queried genomes to a single reference'], 
            ['phylo',     'infer phylogeny and ancestral states from genomic alignments'], 
            ['EBEis',     'in silico serotype prediction for Escherichia coli and Shigella spp.'],                             #
            ['uberBlast', 'Use Blastn, uBlastp, minimap2 and/or mmseqs to identify similar sequences'],     #
            ['clust',     'linear-time clustering of short sequences using mmseqs linclust'],                          #
            ['isCRISPOL', 'in silico prediction of CRISPOL array for Salmonella enterica serovar Typhimurium'],                #
            #            ['hierCC',    'generate hierarchical clusters from cgMLST profiles'],                       #
#            ['RecHMM',    'identify recombination sketches from a SNP matrix'], 
#            ['RecFilter', 'Remove recombination sketches from a SNP matrix'], 
]

def etoki():
    parser = MyParser('EToKi')
    subparser = parser.add_subparsers(title='sub-commands', dest='cmd')
    for cmd, help in commands :
        subparser.add_parser(cmd, help=help)
    arg, others = parser.parse_known_args(sys.argv[1:2])
    if arg.cmd is None :
        parser.print_help()
        sys.exit(2)
    try :
        sys.argv[0] += ' ' + arg.cmd
        exec('from modules.{0} import {0}'.format(arg.cmd))
        eval(arg.cmd)(sys.argv[2:])
    except ValueError as e :
        parser.print_help()


if __name__ == '__main__' :
    etoki()

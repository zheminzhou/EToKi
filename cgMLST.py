import sys, subprocess, json, pandas as pd, numpy as np, shutil, os
from configure import transeq, logger, readFasta
try :
    import ujson as json
except :
    import json

def cgMLST(allele_profile, allele_file) :
    def get_allele_info(alleles) :
        allele_aa = transeq(alleles)
        allele_stat = {}
        for n, s in alleles.iteritems() :
            locus, allele_id = n.rsplit('_', 1)
            if locus not in allele_stat :
                allele_stat[locus] = {}

            if len(s) % 3 > 0 :
                pseudo = 2      # frameshift
            else :
                aa = allele_aa.get(n+'_1', 'A')
                if aa[:-1].find('X') >= 0 :
                    pseudo = 3  # premature
                elif s[:3] not in ('ATG', 'GTG', 'TTG') :
                    pseudo = 4  # no start
                elif aa[-1] != 'X' :
                    pseudo = 5  # no stop
                else :
                    pseudo = 6  # intact
            allele_stat[locus][allele_id] = [len(s), pseudo]
        return allele_stat


    matrix = pd.read_csv(allele_profile, sep='\t', header=None, dtype=str).as_matrix()
    loci = np.array([not m.startswith('#') for m in matrix[0]])
    data = matrix[1:, loci]
    data[ np.in1d(data, ['-', 'n', 'N']).reshape(data.shape) ] = '0'

    data = data.astype(int)
    data[data < 0] = 0
    loci = matrix[0][loci]
    genomes = matrix[1:, 0]

    allele_stat = get_allele_info(readFasta(allele_file))

    genome_stat = { genome:[0 for l in loci] for genome in genomes }
    locus_stat  = [ [locus, len(allele_stat[locus]), 
                     np.mean([ v[0] for v in allele_stat[locus].values() ]), 
                     np.min([ v[0] for v in allele_stat[locus].values() ]), 
                     np.max([ v[0] for v in allele_stat[locus].values() ])] for locus in loci ]
    for g, d in zip(genomes, data) :
        for i, dd in enumerate(d) :
            genome_stat[g][i] = dd*10 + allele_stat.get(loci[i], {}).get(str(dd), [0, 0])[-1]
    return genome_stat, locus_stat


if __name__ == '__main__' :
    genome_stat, locus_stat = cgMLST(sys.argv[1], sys.argv[2])
    json.dump({genome:genome_stat, locus:locus_stat}, open(sys.argv[3], 'wb') )

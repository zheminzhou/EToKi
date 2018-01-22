import sys, subprocess, json, pandas as pd, numpy as np, shutil, os
from EnConf import transeq, logger, readFasta
try :
    import ujson as json
except :
    import json

def EnCore(allele_profile, allele_file) :
    def get_allele_info(alleles) :
        allele_aa = transeq(alleles)
        allele_stat = {}
        for n, s in alleles.iteritems() :
            locus, allele_id = n.rsplit('_', 1)
            if locus not in allele_stat :
                allele_stat[locus] = {}

            if len(s) % 3 > 0 or allele_aa.get(n+'_1', 'M')[:-1].find('X') >= 0 :
                pseudo = 1
            else :
                pseudo = 0
            allele_stat[locus][allele_id] = [len(s), pseudo]
        return allele_stat


    matrix = pd.read_csv(allele_profile, sep='\t', header=None, dtype=str).as_matrix()
    loci = np.array([not m.startswith('#') for m in matrix[0]])
    data = matrix[1:, loci]
    data[ np.in1d(data, ['-', 'n', 'N']).reshape(data.shape) ] = '0'

    data = data.astype(int)
    loci = matrix[0][loci]
    genomes = matrix[1:, 0]

    allele_stat = get_allele_info(readFasta(allele_file))

    genome_stat = { genome:{ locus:[0, 0, 0] for locus in loci } for genome in genomes }
    for g, d in zip(genomes, data) :
        for l, dd in zip(loci, d) :
            genome_stat[g][l][2] = dd
            if str(dd) in allele_stat.get(l, {}) :
                genome_stat[g][l][0] = 2 if allele_stat[l][str(dd)][-1] else 3
                genome_stat[g][l][1] = allele_stat[l][str(dd)][0]
            else :
                genome_stat[g][l][0] = 1 if dd > 0 else 0
    return genome_stat


if __name__ == '__main__' :
    genome_stat = EnCore(sys.argv[1], sys.argv[2])
    with open('static/data_source.js', 'wb') as fout :
        fout.write( 'var matrix=' + json.dumps(genome_stat).replace('},', '},\n') + '\n')
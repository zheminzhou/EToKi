import sys, subprocess, json, pandas as pd, numpy as np
from pyLib.bio_parser import transeq

def get_allele_seq(scheme, locus, unique_allele) :
    seqs = {}
    a= unique_allele.tolist() + [1]
    for id in xrange(0, len(a), 150) :
        cmd = 'curl http://xxxxxx/{0}/alleles?locus={1}&allele_id=in%20({2})&fieldnames=value,value_id'.format(scheme, locus, ','.join([str(x) for x in a[id:(id+150)]]))
        res = json.loads(subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0])
        for r in res :
            seqs['{0}_{1}'.format(locus, r['allele_id'])] = r['seq']
    return seqs

def get_pseudo_allele( allele_seq) :
    pseudo = {n:0 for n in allele_seq}
    for n, s in allele_seq.iteritems() :
        if len(s) % 3 != 0 :
            pseudo[n] = 1
        else :
            aa = transeq(s)
            if 'X' in aa[:-1] or aa[-1] != 'X':
                pseudo[n] = 1
    return pseudo

if __name__ == '__main__' :
    scheme, allele_table = sys.argv[1:3]
    matrix = pd.read_csv(allele_table, sep='\t', header=None, dtype=str).as_matrix()
    data = matrix[1:, 2:].astype(int)
    print '#Locus\tLength\tN_allele\tPresence\tPseudogene'
    for locus, alleles in zip(matrix[0][2:], data.T) :
        unique_allele = np.unique(alleles)
        unique_allele = unique_allele[unique_allele>0]
        allele_seq = get_allele_seq(scheme, locus, unique_allele)
        length = len(allele_seq.get('{0}_1'.format(locus), 0))
        presence = np.sum(alleles > 0 )
        N_allele = unique_allele.size
        p = get_pseudo_allele(allele_seq)
        pseudo = np.sum([p.get('{0}_{1}'.format(locus, a), 0)  for a in alleles ])
        print '{0}\t{1}\t{2}\t{3}\t{4}'.format( locus, length, N_allele, presence, pseudo )

import sys, csv, numpy as np, os
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from numba import njit, jit

try:
    from configure import transeq, uopen
except :
    from .configure import transeq, uopen

try :
    import ujson as json
except :
    import json

def readFasta(fasta, filter=None) :
    sequence = {}
    with uopen(fasta) as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                if not filter or name in filter :
                    sequence[name] = []
            elif len(line) > 0 and not line.startswith('#') and name in sequence :
                sequence[name].extend(line.strip().split())
    for s in sequence :
        sequence[s] = (''.join(sequence[s])).upper()
    return sequence

@njit
def seq_status(seq, aa_seq) :
    if len(seq) % 3 > 0 :
        return 2
    stop = len(seq)
    for i in range( int(len(seq)/3)-1, -1, -1 ) :
        aa_seq[i] = seq[3*i]*676 + seq[3*i+1]*26 + seq[3*i+2]
        if aa_seq[i] in (12844, 12850, 13000) :
            stop = i
    if stop < int(len(seq)/3)-1 :
        return 3
    if aa_seq[0] not in (500, 4556, 13344) :
        return 4
    if stop == len(seq) :
        return 5
    else :
        return 6


def cgMLST(allele_profile, allele_file) :
    def get_allele_info(alleles) :
        if os.path.isfile(allele_file + '.stat') :
            return json.load(open(allele_file + '.stat'))
        allele_stat = {}
        for n, s in alleles.items() :
            locus, allele_id = n.rsplit('_', 1)
            if locus not in allele_stat :
                allele_stat[locus] = {}

            aa_seq = np.zeros(int(len(s)/3), dtype=int)
            seq = np.array(list(s)).astype(bytes).view(np.uint8)-65
            pseudo = seq_status(seq, aa_seq)

            allele_stat[locus][allele_id] = int(allele_id)*1000000 + len(s)*10 + pseudo
        json.dump(allele_stat, open(allele_file + '.stat', 'w'))
        return allele_stat

    with uopen(allele_profile) as fin :
        data = np.array([d for d in csv.reader(fin, delimiter='\t')])

    loc_col = np.array([False] + [not m.startswith('#') for m in data[0, 1:]])
    loci = data[0, loc_col]
    genomes = data[1:, 0]

    data = data[1:, loc_col]
    data[ np.in1d(data, ['-', 'n', 'N']).reshape(data.shape) ] = '0'
    data = data.astype(int)
    data[data < 0] = 0

    allele_names = {'{0}_{1}'.format(locus, allele) for locus, alleles in zip(loci, data.T) for allele in alleles if allele > 0}
    alleles = readFasta(allele_file, allele_names)

    allele_stat = get_allele_info(alleles)
    for g, d in zip(loci, data.T) :
        if g in allele_stat :
            alleles = np.zeros(max(list(map(int, list(allele_stat[g].keys())))+[np.max(d)])+1, int)
            for a, s in allele_stat[g].items() :
                alleles[int(a)] = s
            d[:] = alleles[d]
        else :
            d[:] = 0
    loci = loci[np.sum(data>0, 0)>0]
    data = data[:, np.sum(data>0, 0)>0]
    print('Start with {1} genes in {0} genomes'.format(*data.shape))
    iterations = [{'genePresence':0.5, 'intactCDS':0.5, 'genomeProp':0.4}, 
                  {'genePresence':0.8, 'intactCDS':0.8, 'genomeProp':0.7}, 
                  {'genePresence':0.95, 'intactCDS':0.95, 'oddsRatio':3.}]
    colPresence = np.zeros(data.shape[1], dtype=int)
    for ite, cuts in enumerate(iterations) :
        print('====== Iteration {0} ======'.format(ite))
        
        if 'genePresence' in cuts :
            print('Remove genes that present in < {0} of genomes'.format(cuts['genePresence']))
        if 'intactCDS' in cuts :
            print('Remove genes that are intact in < {0} of genomes.'.format(cuts['intactCDS']))

        pP = np.sum(data>0, 0).astype(float)/data.shape[0]
        pI = np.sum(data % 10 >= 4, 0).astype(float)/np.sum(data>0, 0)
        p = (pP >= cuts.get('genePresence', 0.)) & (pI >= cuts.get('intactCDS', 0.)) & (colPresence == ite)

        if 'oddsRatio' in cuts :
            print('Remove genes that are significantly variable (> {0} sigma) in a Gaussian process regression.'.format(cuts['oddsRatio']))
            x = np.sum((data%1000000/10).astype(int), 0).astype(float)/np.sum(data>0, 0)
            y = np.apply_along_axis(lambda d:np.unique(d[d>0]).size, 0, data)*100./np.sum(data>0, 0)
            x0, y0 = x[p], y[p]
            x1, y1 = x[colPresence == ite], y[colPresence == ite]
            
            kernel = 100.*RBF(length_scale=10.0, length_scale_bounds=(1e-3, 1e3)) + 1.0*WhiteKernel(1e-1, noise_level_bounds=(1e-5, 1e2))
            gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
            gp.fit(x0[:, np.newaxis], y0)
            y_pred, sigma = gp.predict(x[:, np.newaxis], return_std=True)
            oddsRatio = (y - y_pred)/sigma
            p &= (oddsRatio <= cuts['oddsRatio'])

        colPresence[ p ] = ite+1
        print('Remain {0} genes.'.format(np.sum(p)))

        if 'genomeProp' in cuts :
            print('Remove genomes that contain < {0} of genes.'.format(cuts['genomeProp']))
            d2 = data[:, p]
            rowPresence = np.sum(d2>0, 1) >= cuts['genomeProp'] * d2.shape[1]
            removedGenomes = zip(genomes[~rowPresence], np.sum(d2>0, 1)[~rowPresence])
            print('\n'.join(['!!! Removed genomes: {0}:{1}'.format(*r) for r in removedGenomes]))
            genomes = genomes[rowPresence]
            data = data[rowPresence]
            l = np.sum(data>0, 0) > 0
            loci, colPresence, data = loci[l], colPresence[l], data[:, l]
            print('Remain {0} genomes.'.format(data.shape[0]))

    print ('\n====== Results ======')
    print ('#\tGene\t%Presence\t%Intact\tAve.Len\tpVar\tpExp\t# Sigma\tCI(99.7%) Low\tCI(99.7%) High')
    for p, locus, ppP, ppI, xx, yy, py, dy in zip(colPresence, loci, pP, pI, x, y, y_pred, sigma) :
        print ('#{0}\t{1}\t{2:.6f}\t{3:.6f}\t{4:.0f}\t{5:.6f}\t{6:.6f}\t{7:.6f}\t{8:.6f}\t{9:.6f}'.format('cgMLST' if p > ite else 'Filter_'+str(p), locus, ppP*100, ppI*100, xx, yy/100., py/100., (yy-py)/dy, (py-3*dy)/100., (py+3*dy)/100.))
        


if __name__ == '__main__' :
    cgMLST(sys.argv[1], sys.argv[2])

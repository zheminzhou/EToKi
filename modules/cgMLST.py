import sys, csv, numpy as np, os
import tempfile
from multiprocessing import Pool
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from numba import njit, jit

try:
    from configure import transeq, uopen, asc2int
except :
    from .configure import transeq, uopen, asc2int

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
    for n, s in sequence.items() :
        sequence[n] = (''.join(s)).upper()
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

def get_allele_info(allele_npz) :
    alleles = json.load(open(allele_npz, 'rt'))
    #alleles = np.load(allele_npz)['alleles']
    output = []
    for n, s in alleles :
        locus, allele_id = n.rsplit('_', 1)
        if len(s) % 3 > 0:
            pseudo = 2
        else:
            aa_seq = np.zeros(int(len(s) / 3), dtype=int)
            seq = np.array(list(s)).astype(bytes).view(np.uint8) - 65
            pseudo = seq_status(seq, aa_seq)
        output.append([locus, allele_id, int(allele_id)*1000000 + len(s)*10 + pseudo])
    np.savez_compressed(allele_npz, alleles=np.array(output, dtype=object), allow_pickle=True)
    return allele_npz

def cgMLST(args) :
    pool = Pool(8)
    params = getParams(args)
    profile_file, allele_files, prefix = params['profile'], params['alleles'], params['output']

    profile = pd.read_csv(profile_file, sep='\t', header=0, na_filter=None, dtype=str)
    profile = profile.set_index(profile.columns[0])
    loc_col = ~profile.columns.str.startswith('#')
    profile = profile.loc[:, loc_col]
    profile.values[profile.isin(['n', 'N']) | np.vectorize(lambda s:s.startswith('-'))(profile.values) ] = '0'

    allele_names = {'{0}_{1}'.format(k, v) for i, p in profile.iterrows() for k, v in p.items() if not v.startswith('-')}
    alleles = {}
    for allele_file in allele_files :
        alleles.update(readFasta(allele_file, allele_names))

    allele_stat = {}
    with tempfile.TemporaryDirectory(dir='.', prefix='CG_') as tmpdir:
        allele_list = list(alleles.items())
        fnames = []
        for i in np.arange(8):
            fnames.append(os.path.join(tmpdir, '{0}.npz'.format(i)))
            json.dump(allele_list[i::8], open(fnames[i], 'w'))
            #np.savez_compressed(fnames[i], alleles=allele_list[i::8])
        for fname in pool.imap_unordered(get_allele_info, fnames):
            res = np.load(fname, allow_pickle=True)['alleles']
            for locus, allele_id, stat in res :
                if locus not in allele_stat :
                    allele_stat[locus] = {}
                allele_stat[locus][allele_id] = stat

    genomes = profile.index
    data = profile.values
    loci = profile.columns
    for g, d in zip(loci, data.T) :
        if g in allele_stat :
            alleles = {dd:0 for dd in np.unique(d)}
            for a in alleles :
                if a in allele_stat[g] :
                    alleles[a] = allele_stat[g][a]
            d[:] = [ alleles[dd] for dd in d ]
        else :
            d[:] = 0
    data = data.astype(int)
    loci = loci[np.sum(data>0, 0)>0]
    data = data[:, np.sum(data>0, 0)>0]
    print('Start with {1} genes in {0} genomes'.format(*data.shape))
    iterations = [{'genePresence':0.5, 'intactCDS':0.5, 'genomeProp':0.4}, 
                  {'genePresence':0.8, 'intactCDS':0.8, 'genomeProp':0.6},
                  {'genePresence':params['genepresence'], 'intactCDS':params['intactcds'], 'oddsRatio':params['oddratio'], 'geneLength':params['genelength']}]
    colPresence = np.zeros(data.shape[1], dtype=int)
    x = np.sum((data % 1000000 / 10).astype(int), 0).astype(float) / np.sum(data > 0, 0)


    for ite, cuts in enumerate(iterations) :
        print('====== Iteration {0} ======'.format(ite))
        
        if 'genePresence' in cuts :
            print('Remove genes that present in < {0} of genomes'.format(cuts['genePresence']))
        if 'intactCDS' in cuts :
            print('Remove genes that are intact in < {0} of genomes.'.format(cuts['intactCDS']))

        pP = np.sum(data>0, 0).astype(float)/data.shape[0]
        pI = np.sum(data % 10 >= 4, 0).astype(float)/np.sum(data>0, 0)
        p = (pP >= cuts.get('genePresence', 0.)) & (pI >= cuts.get('intactCDS', 0.)) & (colPresence == ite)

        x = np.sum((data % 1000000 / 10).astype(int), 0).astype(float) / np.sum(data > 0, 0)
        if cuts.get('geneLength', 0) > 0 :
            print('Remove genes that have an average allelic length < {0} bps.'.format(cuts['geneLength']))
            p &= (x >= cuts['geneLength'])

        if 'oddsRatio' in cuts :
            print('Remove genes that are significantly variable (> {0} sigma) in a Gaussian process regression. This can take a long time.'.format(cuts['oddsRatio']))
            y = np.apply_along_axis(lambda d: np.unique(d[d > 0]).size, 0, data) * 100. / np.sum(data > 0, 0)
            x0, y0 = x[p], y[p]
            #x1, y1 = x[colPresence == ite], y[colPresence == ite]
            
            kernel = 100.*RBF(length_scale=10.0, length_scale_bounds=(1e-3, 1e3)) + 1.0*WhiteKernel(1e-1, noise_level_bounds=(1e-5, 1e2))
            gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
            gp.fit(x0[:, np.newaxis], y0)
            y_pred, sigma = gp.predict(x[:, np.newaxis], return_std=True)
            oddsRatio = (y - y_pred)/sigma
            if cuts['oddsRatio'] > 0 :
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
    with open('{0}.cgMLST'.format(prefix), 'wt') as fout :
        fout.write('#\tGene\t%Presence\t%Intact\tAve.Len\tpVar\tpExp\t# Sigma\tCI(99.7%) Low\tCI(99.7%) High\n')
        for p, locus, ppP, ppI, xx, yy, py, dy in sorted(zip(colPresence, loci, pP, pI, x, y, y_pred, sigma), key=lambda x:(-x[0], x[1])) :
            fout.write('#{0}\t{1}\t{2:.6f}\t{3:.6f}\t{4:.0f}\t{5:.6f}\t{6:.6f}\t{7:.6f}\t{8:.6f}\t{9:.6f}\n'.format('cgMLST' if p > ite else 'Filter_'+str(p+1), locus, ppP*100, ppI*100, xx, yy/100., py/100., (yy-py)/dy, (py-cuts['oddsRatio']*dy)/100., (py+cuts['oddsRatio']*dy)/100.))
        
def getParams(args) :
    import argparse
    parser = argparse.ArgumentParser(description='cgMLST. Find highly conserved genes as cgMLST candidates. ')
    parser.add_argument('-o', '--output',     help='[Required] Output - prefix for the outputs. ', required=True)
    parser.add_argument('-p', '--profile', help='[Required] Input - summarised profiles by MLSTsum. Can be specified multiple times. ', required=True)
    parser.add_argument('--genepresence', help='[Default: 0.95] Proportion of genome presence for a cgMLST locus. ', default=0.95, type=float)
    parser.add_argument('--intactcds', help='[Default: 0.94] Proportion of intact CDS for a cgMLST locus. ', default=0.94, type=float)
    parser.add_argument('--genelength', help='[Default: 0] minimum length of a cgMLST locus. ', default=0, type=int)
    parser.add_argument('--oddratio', help='[Default: 3.] a control of sequence variability for a cgMLST locus. ', default=3., type=float)
    parser.add_argument('alleles', nargs='*',help='Input - allele sequences')

    return parser.parse_args(args).__dict__


if __name__ == '__main__' :
    #pool = Pool(8)
    cgMLST(sys.argv[1:])

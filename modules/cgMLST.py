import sys, subprocess, json, pandas as pd, numpy as np, shutil, os
try:
    from configure import transeq, logger, readFasta, uopen
except :
    from .configure import transeq, logger, readFasta, uopen
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel

try :
    import ujson as json
except :
    import json
    
    

def cgMLST(allele_profile, allele_file) :
    def get_allele_info(allele_file) :
        if os.path.isfile(allele_file + '.stat') :
            return json.load(open(allele_file + '.stat'))
        alleles = readFasta(allele_file)
        allele_aa = transeq(alleles)
        allele_stat = {}
        for n, s in alleles.items() :
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
            allele_stat[locus][allele_id] = int(allele_id)*1000000 + len(s)*10 + pseudo
        json.dump(allele_stat, open(allele_file + '.stat', 'w'))
        return allele_stat


    matrix = pd.read_csv(allele_profile, sep='\t', header=None, dtype=str).values
    loci = np.array([not m.startswith('#') for m in matrix[0]])
    data = matrix[1:, loci]
    data[ np.in1d(data, ['-', 'n', 'N']).reshape(data.shape) ] = '0'

    data = data.astype(int)
    data[data < 0] = 0
    loci = matrix[0, loci]
    genomes = matrix[1:, 0]

    allele_stat = get_allele_info(allele_file)

    for g, d in zip(loci, data.T) :
        if g in allele_stat :
            alleles = np.zeros(max(map(int, list(allele_stat[g].keys())))+1, int)
            for a, s in allele_stat[g].items() :
                alleles[int(a)] = s
            d[:] = alleles[d]
        else :
            d[:] = 0
    loci = loci[np.sum(data>0, 0)>0]
    data = data[:, np.sum(data>0, 0)>0]
    print('Start with {1} genes in {0} genomes'.format(*data.shape))
    iterations = [{'genePresence':0.6, 'intactCDS':0.6, 'genomeProp':0.5}, 
                  {'genePresence':0.8, 'intactCDS':0.8, 'genomeProp':0.7}, 
                  {'genePresence':0.98, 'intactCDS':0.94, 'oddsRatio':3.}]
    for ite, cuts in enumerate(iterations) :
        print('====== Iteration {0} ======'.format(ite))
        
        colPresence = np.ones(data.shape[1], dtype=bool)
        if 'genePresence' in cuts :
            print('Remove genes that present in < {0} of genomes'.format(cuts['genePresence']))
            pP = np.sum(data>0, 0).astype(float)/data.shape[0]
            colPresence &= ( pP >= cuts['genePresence'] )
            
        if 'intactCDS' in cuts :
            print('Remove genes that are intact in < {0} of genomes.'.format(cuts['intactCDS']))
            pI = np.sum(data % 10 >= 4, 0).astype(float)/np.sum(data>0, 0)
            colPresence &= (pI >= cuts['intactCDS'] )
                
        if 'oddsRatio' in cuts :
            print('Remove genes that are significantly variable (> {0} sigma) in a Gaussian process regression.'.format(cuts['oddsRatio']))
            x = np.sum((data%1000000/10).astype(int), 0).astype(float)/np.sum(data>0, 0)
            y = np.apply_along_axis(lambda d:np.unique(d).size, 0, data)*100./np.sum(data>0, 0)
            x0, y0 = x[colPresence], y[colPresence]
            
            kernel = 100.*RBF(length_scale=1.0, length_scale_bounds=(1e-3, 1e3)) + 1.0*WhiteKernel(1e-1, noise_level_bounds=(1e-5, 1e2))
            gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
            gp.fit(x0[:, np.newaxis], y0)
            y_pred, sigma = gp.predict(x[:, np.newaxis], return_std=True)
            oddsRatio = (y - y_pred)/sigma
            colPresence &= (oddsRatio <= cuts['oddsRatio'])

        print('Remain {0} genes.'.format(np.sum(colPresence)))
        if ite >= len(iterations) - 1 :
            bins = np.arange(0.5,1.01, 0.005)
            counts = np.bincount(np.digitize(pP, bins-0.000001), minlength=bins.size)
            print ('======')
            for b, c in zip(bins[counts>0], counts[counts>0]) :
                print ('#GenePresence\t{0:.1f}%\t{1}'.format(b*100-0.5, c))
            print ('======')
            counts = np.bincount(np.digitize(pI, bins-0.000001), minlength=bins.size)
            for b, c in zip(bins[counts>0], counts[counts>0]) :
                print ('#Pseudogenes\t{0:.1f}%\t{1}'.format(abs(100.5-b*100), c))
            print ('======')
            bins = np.arange(-20., max(oddsRatio)+0.2, 0.2)
            counts = np.bincount(np.digitize(oddsRatio, bins), minlength=bins.size)
            for b, c in zip(bins[counts>0], counts[counts>0]) :
                print ('#oddsRatio\t{0:.1f}\t{1}'.format(b, c))


            print ('\n====== Details ======')
            print ('#\tGene\t%Presence\t%Intact\tAve.Len\t%Variants\tExpected\tOdds\tCI(99.7%) Low\tCI(99.7%) High')
            for p, locus, ppP, ppI, xx, yy, py, dy in zip(colPresence, loci, pP, pI, x, y, y_pred, sigma) :
                print ('#{0}\t{1}\t{2:.1f}\t{3:.1f}\t{4:.0f}\t{5:.2f}\t{6:.2f}\t{7:.2f}\t{8:.2f}\t{9:.2f}'.format('cgMLST' if p else '', locus, ppP*100, ppI*100, xx, yy, py, (yy-py)/dy, py-3*dy, py+3*dy))

        loci = loci[colPresence]
        data = data[:, colPresence]

        if 'genomeProp' in cuts :
            print('Remove genomes that contain < {0} of genes.'.format(cuts['genomeProp']))
            rowPresence = np.sum(data>0, 1) >= cuts['genomeProp'] * data.shape[1]
            removedGenomes = zip(genomes[~rowPresence], np.sum(data>0, 1)[~rowPresence])
            print('\n'.join(['!!! Removed genomes: {0}:{1}'.format(*r) for r in removedGenomes]))
            genomes = genomes[rowPresence]
            data = data[rowPresence]
            print('Remain {0} genomes.'.format(data.shape[0]))
        
    #return data


if __name__ == '__main__' :
    cgMLST(sys.argv[1], sys.argv[2])

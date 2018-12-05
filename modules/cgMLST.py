import sys, subprocess, json, pandas as pd, numpy as np, shutil, os
from configure import transeq, logger, readFasta, uopen
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel

try :
    import ujson as json
except :
    import json
    
    

def cgMLST(allele_profile, allele_file) :
    def get_allele_info(allele_file) :
        alleles = readFasta(allele_file)
        if os.path.isfile(allele_file + '.faa.json') :
            allele_aa = json.load(open(allele_file + '.faa.json'))
        else :
            allele_aa = transeq(alleles)
            json.dump(allele_aa, open(allele_file + '.faa.json', 'w'))
        allele_stat = {}
        for n, s in alleles.items() :
            locus, allele_id = n.rsplit('_', 1)
            allele_id = int(allele_id)
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
            allele_stat[locus][allele_id] = allele_id*1000000 + len(s)*10 + pseudo
        return allele_stat


    matrix = pd.read_csv(allele_profile, sep='\t', header=None, dtype=str).values
    loci = np.array([not m.startswith('#') for m in matrix[0]])
    data = matrix[1:, loci]
    data[ np.in1d(data, ['-', 'n', 'N']).reshape(data.shape) ] = '0'

    data = data.astype(int)
    data[data < 0] = 0
    loci = matrix[0][loci]
    genomes = matrix[1:, 0]

    allele_stat = get_allele_info(allele_file)

    for g, d in zip(loci, data.T) :
        if g in allele_stat :
            alleles = np.zeros(max(list(allele_stat[g].keys()))+1, int)
            for a, s in allele_stat[g].items() :
                alleles[a] = s
            d[:] = alleles[d]
        else :
            d[:] = 0
    print('Start with {1} genes in {0} genomes'.format(*data.shape))
    for ite, cuts in enumerate([{'genePresence':0.5, 'genomeProp':0.5}, 
                                {'genePresence':0.6, 'genomeProp':0.6}, 
                                {'genePresence':0.7, 'genomeProp':0.7}, 
                                {'genePresence':0.8, 'genomeProp':0.8}, 
                                {'genePresence':0.98, 'intactCDS':0.94, 'oddRatio':3}]) :
        if 'genePresence' in cuts :
            print('Remove genes that present in < {0} of genomes'.format(cuts['genePresence']))
            colPresence = np.sum(data>0, 0) >= cuts['genePresence'] * data.shape[0]
            loci = loci[colPresence]
            data = data[:, colPresence]
            print('Remain {0} genes.'.format(data.shape[1]))
        if 'intactCDS' in cuts :
            print('Remove genes that are intact in < {0} of genomes.'.format(cuts['intactCDS']))
            colPresence = np.sum(data % 10 >= 4, 0) >= cuts['genePresence'] * np.sum(data>0, 0)
            loci = loci[colPresence]
            data = data[:, colPresence]
            print('Remain {0} genes.'.format(data.shape[1]))
        if 'oddRatio' in cuts :
            print('Remove genes that are significantly variable in a piecewise regression.'.format(cuts['oddRatio']))
            x = np.sum((data%1000000/10).astype(int), 0)/np.sum(data>0, 0)
            y = np.apply_along_axis(lambda d:np.unique(d).size, 0, data)*100./np.sum(data>0, 0)

            kernel = 1.0 * RBF(length_scale=1) + 1.0 * WhiteKernel()
            gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
            gp.fit(x[:, np.newaxis], y)
            y_pred, sigma = gp.predict(x.reshape(-1, 1), return_std=True)
            print (data)
        if 'genomeProp' in cuts :
            print('Remove genomes that contain < {0} of genes.'.format(cuts['genomeProp']))
            rowPresence = np.sum(data>0, 1) >= cuts['genomeProp'] * data.shape[1]
            removedGenomes = genomes[~rowPresence]
            if removedGenomes.size :
                print('!!! Removed genomes: {0}'.format(','.join(removedGenomes)))
            genomes = genomes[rowPresence]
            data = data[rowPresence]
            print('Remain {0} genomes.'.format(data.shape[0]))
        print('====== SUMMARY ======')
        #print (data)
        
    return data


if __name__ == '__main__' :
    genome_stat, locus_stat = cgMLST(sys.argv[1], sys.argv[2])
    json.dump({genome:genome_stat, locus:locus_stat}, open(sys.argv[3], 'wb') )

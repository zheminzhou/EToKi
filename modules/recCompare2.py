from ete3 import Tree
import numpy as np, pandas as pd, sys, os, re
try:
    from configure import uopen
except :
    from .configure import uopen


def readRec(fname) :
    rec = {}
    with uopen(fname) as fin :
        for line in fin :
            part = line.strip().split('\t')
            key = tuple([part[0], part[1]])
            if key not in rec :
                rec[key] = []
            rec[key].append([int(part[2]), int(part[3])])
    return rec
if __name__ == '__main__' :
    snpFile, recFile1, recFile2 = sys.argv[1:]
    
    rec1 = readRec(recFile1)
    rec2 = readRec(recFile2)

    # load SNP matrix
    with uopen(snpFile) as fin :
        for line in fin :
            if line.startswith('##') : 
                continue
            else :
                header = line.strip().split('\t')
                header = { h:i for i, h in enumerate(header) }
                mat = pd.read_csv(fin, sep='\t', header=None, dtype=str).values.T
                break
    
    headers = list(header.keys())
    results = np.zeros([int(len(headers)*(len(headers)-1)/2), 4], dtype=int)
    
    k = -1
    for i1, t1 in enumerate(headers) :
        for i2 in np.arange(i1+1, len(headers)) :
            k += 1
            if k % 100 == 0 :
                sys.stderr.write('# {0}\n'.format(k))
            t2 = headers[i2]
            key = tuple(sorted([t1, t2]))
            snps = mat[1, (mat[header[t1]] != mat[header[t2]]) & (mat[header[t1]] != '-') & (mat[header[t1]] != '-')].astype(int)
            
            tags = np.zeros(snps.shape, dtype=int)
            for rId, r in enumerate([rec1.get(key, []), rec2.get(key, [])]) :
                for rr in r :
                    tags[(snps >= rr[0]) & (snps <= rr[1])] += 2**rId

            tags = np.bincount(tags)
            results[k, :tags.size] = tags
    print ('#TargetFile\t#Accuracies\t#Precisions\t#Sensitivities')
    print ('{0}\t{1}\t{2}\t{3}'.format( recFile2, np.mean(np.sum(results.T[3]+results.T[0])/np.sum(results)), np.mean(np.sum(results.T[3])/np.sum(results.T[2:])), np.mean(np.sum(results.T[3])/np.sum(results.T[1]+results.T[3])) ))
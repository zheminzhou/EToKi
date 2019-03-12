# MGplacer
import sys, os, numpy as np, pandas as pd, subprocess, re
import pickle
from scipy.stats import poisson
from ete3 import Tree
try :
    from configure import uopen
except :
    from .configure import uopen
try :
    xrange(3)
except :
    xrange = range

def prepareTree(treeFile) :
    tre = Tree(treeFile, format=1)
    for node in tre.traverse('postorder') :
        if not node.is_leaf() and node.name != '' and node.up :
            child = node.remove_child(node.children[-1])
            node.add_child(name=node.name, dist=0.0, support=0.0)
            node.add_child(child)
            node.name = ''
    return tre

def readMutations(mutationFile) :
    accLength = 0
    seqRange = {}
    missingBlocks = []
    with uopen(mutationFile) as fin :
        for line in fin :
            if line.startswith('##') :
                if line.startswith('## Sequence_length:') :
                    part = line.strip().split()
                    seqRange[part[2]] = [accLength, accLength+int(part[3])]
                    accLength += int(part[3])
                elif line.startswith('## Missing_region:') :
                    part = line.strip().split()
                    acc = seqRange[part[2]][0]
                    s, e = int(part[3]) + acc, int(part[4]) + acc
                    if [s-1, e] == seqRange[part[2]] :
                        diff = e - s + 1
                        seqRange.pop(part[2])
                        for n in seqRange :
                            if seqRange[n][0] >= e : seqRange[n][0] -= diff
                            if seqRange[n][1] >= e : seqRange[n][1] -= diff
                        for blk in missingBlocks :
                            if blk[0] >= e : blk[0] -= diff
                            if blk[1] >= e : blk[1] -= diff
                    else :
                        if e - s + 1 >= 500 :
                            missingBlocks.append([s, e])
            elif line.startswith('#') :
                mat = pd.read_csv(fin, sep='\t', header=None, dtype=str).values
                break
    mat = mat[np.vectorize(lambda x:len(x))(mat.T[4]) == 4]
    sites = np.vstack([np.vectorize(lambda x:seqRange.get(x, [-1])[0])(mat.T[1]), mat.T[2].astype(int)])
    mat = mat[sites[0]>=0]
    sites = np.sum(sites[:, sites[0]>=0], 0)
    mutBlocks = list(zip(mat.T[0], sites))
    return seqRange, missingBlocks, mutBlocks

def readRecRegions(recFile, seqRange) :
    recBlocks = []
    with uopen(recFile) as fin :
        for line in fin :
            part = line.strip().split()
            if part[0] == 'Importation' :
                if part[2] not in seqRange : continue
                acc = seqRange[part[2]][0]
                recBlocks.append([part[1], acc+int(part[3]), acc+int(part[4])])
            elif len(part) == 3 :
                try :
                    recBlocks.append([part[0], int(part[1]), int(part[2])])
                except :
                    pass
    return recBlocks

def preparePhandango(prefix, tre, totalLength, missingBlocks, mutBlocks, recBlocks) :
    tre.write(outfile=prefix + '.internal.tre', format=1)
    leaves = {leaf:1 for leaf in tre.get_leaf_names()}
    headLeaf = tre.get_leaf_names()[0]
    #with open(prefix + '.missing.gff', 'w') as fout :
        #for block in missingBlocks :
            #fout.write('Missing\tEMBL\tCDS\t{0}\t{1}\t.\t+\t.\tID="missing_region"\n'.format(*block))

    
    with open(prefix + '.rec_tabular.txt', 'w') as fout :
        fout.write('LIST OF FOREIGN GENOMIC SEGMENTS:\nStart\tEnd\tOrigin\tHomeCluster\tBAPSIndex\tStrainName\n')
        for block in missingBlocks :
            fout.write('{1}\t{2}\t2\t2\tMissings\t{0}\n'.format(headLeaf, *block))
        #for block in mutBlocks :
            #fout.write('{1}\t{1}\t2\t2\tMutation\t{0}\n'.format(*block))
        for block in recBlocks :
            fout.write('{1}\t{2}\t1\t1\tImports\t{0}\n'.format(*block))

if __name__ == '__main__' :
    prefix, treeFile, mutationFile, recFile = sys.argv[1:]
    
    tre = prepareTree(treeFile)
    
    seqRange, missingBlocks, mutBlocks = readMutations(mutationFile)
    
    recBlocks = readRecRegions(recFile, seqRange)
    
    preparePhandango(prefix, tre, max([v[1] for v in seqRange.values()]), missingBlocks, mutBlocks, recBlocks)
    
    
# hierBG.py
# hierarchical Burst-likc Grouping of MLST allelic profiles
#
# Author: Zhemin Zhou
# Lisence: GPLv3
#
# New assignment: hierBG.py <allelic_profile> <memmap_dump> > <text_output>
# Incremental assignment: hierBG.py <allelic_profile> <memmap_dump> <old_dump> > <text_output>
# Input format:
# ST_id gene1 gene2
# 1 1 1
# 2 1 2
# ...

import pandas as pd, numpy as np, sys, time, os, numba as nb
from multiprocessing import Pool
try:
    xrange(3)
except :
    xrange = range

def get_distance(idx) :
    global mat
    n_loci = mat.shape[1] - 1
    if idx == 0 or idx >= mat.shape[0] :
        return np.zeros(shape=[0, 4], dtype=int)
    profile = mat[idx]
    d1 = (n_loci * np.sum((profile[1:] != mat[:idx, 1:]) & (profile[1:] > 0), 1)/np.sum(profile[1:] > 0).astype(float)+0.5).astype(int)+1
    d2 = n_loci - np.sum((profile[1:] == mat[:idx, 1:]) & (profile[1:] > 0), 1)+1
    d1[d1>d2] = d2[d1>d2]
    dists = np.vstack([[idx for i in np.arange(idx)], np.arange(idx), d1, d2]).astype(int).T
    return dists[np.argsort(dists.T[2], kind='mergesort')]

@nb.jit(nopython=True)
def assignment(dists, res) :
    for id in xrange(len(dists)) :
        idx, ref, d1, d2 = dists[id]
        for d in xrange(d1, n_loci+1) :
            if res[idx, d] != res[ref, d] :
                if d >= res[idx, 0] :
                    if d >= d2 :
                        if res[idx, d] < res[ref, d] :
                            grps = [res[idx, d], res[ref, d]]
                        else :
                            grps = [res[ref, d], res[idx, d]]
                        #grps = sorted([res[idx, d], res[ref, d]])
                        res[:idx, d][res[:idx, d] == grps[1]] = grps[0]
                        res[idx, d] = grps[0]
                else :
                    if res[idx, d] < res[ref, d] :
                        res[:idx, d][res[:idx, d] == res[ref, d]] = res[idx, d]
                    else :
                        res[idx, d:] = res[ref, d:]
                        break
            else :
                break
        if res[idx, 0] > d1 :
            res[idx, 0] = d1
    return

fixedGroup = True
if __name__ == '__main__' :
    ot = time.time()
    profile_file, cluster_file, old_cluster = (sys.argv[1:] + [''])[:3]
    mat = pd.read_csv(profile_file, sep='\t', header=None, dtype=str).values
    allele_columns = np.array([i == 0 or (not h.startswith('#')) for i, h in enumerate(mat[0])])
    mat = mat[1:, allele_columns].astype(int)
    n_loci = mat.shape[1] - 1
    
    sys.stderr.write('{0}: Loaded in allelic profiles with dimension: {1} and {2}. The first column is assumed to be type id. \n'.format(time.time() - ot, *mat.shape))
#    absence = np.sum(mat <= 0, 1)
#    mat = mat[np.argsort(absence, kind='mergesort')]

    res = np.memmap(cluster_file, shape=mat.shape, dtype=int, mode='w+')
    res.T[0] = n_loci+1
    res[0] = mat[0, 0]
    if os.path.isfile(old_cluster) :
        cls = np.memmap(old_cluster, dtype=int, mode='r')
        n_row = cls.size / mat.shape[1]
        cls.shape = (n_row, mat.shape[1])
        typed = {c[0]:id for id, c in enumerate(cls) if c[0] > 0}
        if len(typed) > 0 :
            sys.stderr.write('{0}: Loaded in {1} old hBG assignments. \n'.format(time.time()-ot, len(typed)))
            mat_idx = np.array([t in typed for t in mat.T[0]])
            mat[:] = np.vstack([mat[mat_idx], mat[(mat_idx) == False] ])
    else :
        typed = {}
    
    sys.stderr.write('{0}: Start hBG assignments\n'.format(time.time()-ot))
    pool = Pool(10)
    for index in xrange(0, mat.shape[0], 100) :
        to_run = []
        for idx in np.arange(index, index+100) :
            if idx < mat.shape[0] :
                if mat[idx, 0] in typed :
                    res[idx, :] = cls[typed[mat[idx, 0]], :]
                else :
                    res[idx, 1:] = mat[idx, 0]
                    to_run.append(idx)
        if len(to_run) == 0 :
            continue
        if not fixedGroup :
            dists = np.vstack( pool.map(get_distance, to_run) )
        else :
            dists = np.vstack( [d[0] for d in pool.map(get_distance, to_run)] )
        assignment(dists, res)
        sys.stderr.write('{0}: Assigned {1} of {2} types into hBGs.\n'.format(time.time() - ot, index, mat.shape[0]))
    res.T[0] = mat.T[0]
    sys.stdout.write('#ST_id\t{0}\n'.format('\t'.join(['g'+str(id) for id in np.arange(n_loci)])))

    for r in res[np.argsort(res.T[0])] :
        sys.stdout.write('\t'.join([str(rr) for rr in r])+'\n')
    del res

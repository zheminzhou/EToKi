# hierCC.py
# hierarchical Clustering Complex of MLST allelic profiles
#
# Author: Zhemin Zhou
# Lisence: GPLv3
#
# New assignment: hierCC.py <allelic_profile> <memmap_dump> > <text_output>
# Incremental assignment: hierBG.py <allelic_profile> <memmap_dump> <old_dump> > <text_output>
# Input format:
# ST_id gene1 gene2
# 1 1 1
# 2 1 2
# ...

import os, sys, time, argparse
import pandas as pd, numpy as np, numba as nb
from multiprocessing import Pool
from configure import uopen

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

def get_args(args) :
    parser = argparse.ArgumentParser(description='''hierCC takes allelic profile (as in https://pubmlst.org/data/) and
work out specialised single linkage clustering result of all the profiles in the list.''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', '--profile', help='[INPUT; REQUIRED] name of the profile file. Can be GZIPed.', required=True)
    parser.add_argument('-o', '--output', help='[OUTPUT; REQUIRED] Prefix for the output files. These include a NUMPY and TEXT verions of the same clustering result', required=True)
    parser.add_argument('-i', '--incremental', help='[INPUT; optional] The NUMPY version of an old clustering result', default='')
    parser.add_argument('-d', '--delta', help='[optional] comma delimited list of threshold (delta). All values are included by default.', default=None)
    parser.add_argument('--immutable', help='[optional] Use a immutable clustering system. The designations of old profiles are immutable. Faster but leads to non-optimal assignment.', default=False, action='store_true')

    return parser.parse_args(args)    

def hierCC(args) :
    params = get_args(args)
    ot = time.time()
    profile_file, cluster_file, old_cluster = params.profile, params.output+'.npy', params.incremental
    
    global mat, n_loci
    mat = pd.read_csv(profile_file, sep='\t', header=None, dtype=str).values
    allele_columns = np.array([i == 0 or (not h.startswith('#')) for i, h in enumerate(mat[0])])
    mat = mat[1:, allele_columns].astype(int)
    n_loci = mat.shape[1] - 1
    
    sys.stderr.write('{0}: Loaded in allelic profiles with dimension: {1} and {2}. The first column is assumed to be type id. \n'.format(time.time() - ot, *mat.shape))
    if not params.immutable :
        absence = np.sum(mat <= 0, 1)
        mat = mat[np.argsort(absence, kind='mergesort')]

    res = np.memmap(cluster_file, shape=mat.shape, dtype=int, mode='w+')
    res.T[0] = n_loci+1
    res[0] = mat[0, 0]
    if os.path.isfile(old_cluster) :
        cls = np.memmap(old_cluster, dtype=int, mode='r')
        n_row = cls.size / mat.shape[1]
        cls.shape = (n_row, mat.shape[1])
        typed = {c[0]:id for id, c in enumerate(cls) if c[0] > 0}
        if len(typed) > 0 :
            sys.stderr.write('{0}: Loaded in {1} old hierCC assignments. \n'.format(time.time()-ot, len(typed)))
            mat_idx = np.array([t in typed for t in mat.T[0]])
            mat[:] = np.vstack([mat[mat_idx], mat[(mat_idx) == False] ])
    else :
        typed = {}
    
    sys.stderr.write('{0}: Start hierCC assignments\n'.format(time.time()-ot))
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
        if not params.immutable :
            dists = np.vstack( pool.map(get_distance, to_run) )
        else :
            dists = np.vstack( [d[0] for d in pool.map(get_distance, to_run) if len(d)] )
        assignment(dists, res)
        sys.stderr.write('{0}: Assigned {1} of {2} types into hierCC.\n'.format(time.time() - ot, index, mat.shape[0]))
    res.T[0] = mat.T[0]
    

    if not params.delta :
        with uopen(params.output+'.hierCC.gz', 'w') as fout :
            fout.write('#ST_id\t{0}\n'.format('\t'.join(['d'+str(id) for id in np.arange(n_loci)])))
            for r in res[np.argsort(res.T[0])] :
                fout.write('\t'.join([str(rr) for rr in r])+'\n')
    else :
        deltas = map(int, params.delta.split(','))
        with uopen(params.output+'.hierCC.gz', 'w') as fout :
            fout.write('#ST_id\t{0}\n'.format('\t'.join(['d'+str(id) for id in deltas])))
            for r in res[np.argsort(res.T[0])] :
                fout.write('\t'.join([str(r[id+1]) for id in deltas])+'\n')
    del res
    sys.stderr.write('NUMPY clustering result (for incremental hierCC): {0}.npy\n'.format(params.output))
    sys.stderr.write('TEXT clustering result (for visual inspection): {0}.hierCC.gz\n'.format(params.output))


mat, n_loci = None, None
if __name__ == '__main__' :
    hierCC(sys.argv[1:])
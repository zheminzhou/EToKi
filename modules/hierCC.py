# hierCC.py
# hierarchical Clustering Complex of MLST allelic profiles
#
# Author: Zhemin Zhou
# Lisence: GPLv3
#
# New assignment: hierCC.py -p <allelic_profile> -o <output_prefix>
# Incremental assignment: hierCC.py -p <allelic_profile> -o <output_prefix> -i <old_cluster_npy>
# Input format:
# ST_id gene1 gene2
# 1 1 1
# 2 1 2
# ...

import os, sys, time, argparse
import pandas as pd, numpy as np, numba as nb
from multiprocessing import Pool
try:
    from configure import uopen, logger, xrange
except :
    from .configure import uopen, logger, xrange


@nb.jit(nopython=True)
def assignment2(dists, res) :
    for id in xrange(len(dists)) :
        idx, ref, jd = dists[id]
        res[idx, jd:] = res[ref, jd:]

def get_distance2(idx) :
    global mat
    n_loci = mat.shape[1] - 1
    if idx == 0 or idx >= mat.shape[0] :
        return np.zeros(shape=[0, 3], dtype=int)
    profile = mat[idx]
    s = np.sum((profile[1:] == mat[:idx, 1:]) & (profile[1:] > 0), 1)
    ql = np.sum(profile[1:] > 0)
    rl = np.sum(mat[:idx, 1:] > 0, 1)
    rll = n_loci - np.max([(n_loci - ql)*3, int((n_loci - ql)+n_loci*0.03+0.5)])
    rl[rl < rll] = rll
    rl[rl > ql ] = ql
    d = ((rl - s).astype(float)*n_loci/rl+0.5).astype(int)+1
    dists = np.vstack([np.repeat(idx, idx), np.arange(idx), d]).astype(int).T
    return dists[np.argmin(dists.T[2])]

#@nb.jit(nopython=True)
def assignment3(dists, res, encode, presence) :
    n_loci = mat.shape[1] - 1
    for id in xrange(len(dists)) :
        idx, ref, _, s, ql = dists[id]
        gl = presence[encode[res[ref, 1:]], np.arange(1, res.shape[1])]
        gl[gl > ql] = ql
        d = (n_loci*(gl-s).astype(float)/gl+0.5).astype(int)
        jd = np.argmax( (d - np.arange(n_loci)) <= 0 ) + 1
        res[idx, jd:] = res[ref, jd:]
        gl = presence[encode[res[idx, jd:]], np.arange(jd, presence.shape[1])]
        presence[encode[res[idx, jd:]][gl < ql], np.arange(jd, presence.shape[1])[gl < ql]] = ql
        
        

def get_distance3(idx) :
    global mat
    n_loci = mat.shape[1] - 1
    if idx == 0 or idx >= mat.shape[0] :
        return np.zeros(shape=[0, 3], dtype=int)
    profile = mat[idx]
    s = np.sum((profile[1:] == mat[:idx, 1:]) & (profile[1:] > 0), 1)
    ql = np.sum(profile[1:] > 0)
    rl = np.sum(mat[:idx, 1:] > 0, 1)
    rll = n_loci - np.max([(n_loci - ql)*3, int((n_loci - ql)+n_loci*0.03+0.5)])
    rl[rl < rll] = rll
    rl[rl > ql ] = ql
    d = ((rl - s).astype(float)*n_loci/rl+0.5).astype(int)
    dists = np.vstack([np.repeat(idx, idx), np.arange(idx), d, s, np.repeat(ql, idx)]).astype(int).T
    return dists[np.argmin(dists.T[2])]


def get_distance(idx) :
    global mat
    n_loci = mat.shape[1] - 1
    if idx == 0 or idx >= mat.shape[0] :
        return np.zeros(shape=[0, 4], dtype=int)
    profile = mat[idx]
    ql = np.max(1.0, np.sum(profile[1:] > 0).astype(float))
    d1 = (n_loci * np.sum((profile[1:] != mat[:idx, 1:]) & (profile[1:] > 0), 1)/ql+0.5).astype(int)+1
    d2 = n_loci - np.sum((profile[1:] == mat[:idx, 1:]) & (profile[1:] > 0), 1)+1
    d1[d1>d2] = d2[d1>d2]
    dists = np.vstack([np.repeat(idx, idx), np.arange(idx), d1, d2]).astype(int).T
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
    profile_file, cluster_file, old_cluster = params.profile, params.output+'.npz', params.incremental
    
    global mat, n_loci
    mat = pd.read_csv(profile_file, sep='\t', header=None, dtype=str).values
    allele_columns = np.array([i == 0 or (not h.startswith('#')) for i, h in enumerate(mat[0])])
    mat = mat[1:, allele_columns].astype(int)
    n_loci = mat.shape[1] - 1
    
    logger('{0}: Loaded in allelic profiles with dimension: {1} and {2}. The first column is assumed to be type id.'.format(time.time() - ot, *mat.shape))
    if not params.immutable :
        absence = np.sum(mat <= 0, 1)
        mat = mat[np.argsort(absence, kind='mergesort')]

    if params.immutable :
        presence = np.repeat(np.sum(mat[:, 1:]>0, 1), mat.shape[1]).reshape(mat.shape)
        presence.T[0] = mat.T[0]
        encode = np.zeros(np.max(mat.T[0])+1, dtype=int)
        encode[mat.T[0]] = np.arange(mat.shape[0])
    else :
        presence = None
    if os.path.isfile(old_cluster) :
        od = np.load(old_cluster)
        cls = od['hierCC']
        if params.immutable :
            if 'presence' in od :
                old_presence = od['presence'] 
                if old_presence is not None and old_presence.size > 0 :
                    for c in old_presence :
                        if c[0] > 0 :
                            presence[encode[c[0]], :] = c
                else :
                    old_presence = None
            else :
                old_presence = None

        typed = {c[0]:id for id, c in enumerate(cls) if c[0] > 0}
        if len(typed) > 0 :
            logger('{0}: Loaded in {1} old hierCC assignments.'.format(time.time()-ot, len(typed)))
            mat_idx = np.array([t in typed for t in mat.T[0]])
            mat[:] = np.vstack([mat[mat_idx], mat[(mat_idx) == False] ])
    else :
        typed = {}

    logger('{0}: Start hierCC assignments'.format(time.time()-ot))
    pool = Pool(10)
    
    res = np.repeat(mat.T[0], mat.shape[1]).reshape(mat.shape)
    res[1:, 0] = n_loci+1
    for index in xrange(0, mat.shape[0], 100) :
        to_run = []
        for idx in np.arange(index, index+100) :
            if idx < mat.shape[0] :
                if mat[idx, 0] in typed :
                    res[idx, :] = cls[typed[mat[idx, 0]], :]
                    if params.immutable and old_presence is None :
                        gl = presence[encode[res[idx, 1:]], np.arange(1, presence.shape[1])]
                        ql = np.sum(mat[idx, 1:] > 0)
                        presence[encode[res[idx, 1:]][gl < ql], np.arange(1, presence.shape[1])[gl < ql]] = ql
                else :
                    to_run.append(idx)
        if len(to_run) == 0 :
            continue
        if not params.immutable :
            dists = np.vstack( pool.map(get_distance, to_run) )
            assignment(dists, res)
        else :
            dists = np.vstack( pool.map(get_distance3, to_run) )
            assignment3(dists, res, encode, presence)

        logger('{0}: Assigned {1} of {2} types into hierCC.'.format(time.time() - ot, index, mat.shape[0]))
    res.T[0] = mat.T[0]
    np.savez_compressed(cluster_file, hierCC=res, presence=presence)

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
    logger('NUMPY clustering result (for incremental hierCC): {0}.npz'.format(params.output))
    logger('TEXT  clustering result (for visual inspection): {0}.hierCC.gz'.format(params.output))


mat, n_loci = None, None
if __name__ == '__main__' :
    hierCC(sys.argv[1:])

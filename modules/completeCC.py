# hierCC.py
# comlete linkage Clustering of MLST allelic profiles
#
# Author: Zhemin Zhou
# Lisence: GPLv3
#
# New assignment: completeCC.py -p <allelic_profile> -o <output_prefix>
# Incremental assignment: completeCC.py -p <allelic_profile> -o <output_prefix> -i <old_npz> --immutable
# Input format:
# ST_id gene1 gene2
# 1 1 1
# 2 1 2
# ...

import os, sys, time, argparse
from scipy.cluster.hierarchy import linkage
import pandas as pd, numpy as np, numba as nb
from multiprocessing import Pool
import scipy.spatial.distance as ssd

try:
    from configure import uopen, logger, xrange
except:
    from .configure import uopen, logger, xrange

def get_distances(prefix, profile, pool) :
    logger('Calculating pairwise distance ...')
    np.save(prefix+'.profile.npy', profile)
    indices = np.array([[profile.shape[0]*(v/10.)**2+0.5, profile.shape[0]*((v+1)/10.)**2+0.5] for v in np.arange(10, dtype=float)], dtype=int)
    subfiles = pool.map(parallel_distance, [[prefix+'.profile.npy', prefix+'.dist.{0}.npy', idx] for idx in indices])
    prof_dist = np.hstack([ np.load(subfile) for subfile in subfiles ])
    prof_dist += prof_dist.T
    for subfile in subfiles :
        os.unlink(subfile)
    return prof_dist

def parallel_distance(callup) :
    prof_file, sub_prefix, index_range = callup
    profiles = np.load(prof_file)

    res = profile_distance(profiles, index_range)
    subfile = sub_prefix.format(index_range[0])
    np.save(subfile, res)
    return subfile

def profile_distance(mat, index_range=None) :
    if index_range is None :
        index_range = [0, mat.shape[0]]

    distances = np.zeros(shape=[mat.shape[0], index_range[1] - index_range[0]])
    for i2, idx in enumerate(np.arange(*index_range)) :
        profile = mat[idx]
        s = np.sum((profile[1:] == mat[idx+1:, 1:]) & (profile[1:] > 0), 1)
        ql = np.sum(profile[1:] > 0)
        rl = np.sum(mat[idx+1:, 1:] > 0, 1)
        rll = n_loci - np.max([(n_loci - ql) * 3, int((n_loci - ql) + n_loci * 0.03 + 0.5)])
        rl[rl < rll] = rll
        rl[rl > ql] = ql
        rl[rl < 0.5] = 0.5
        diffs = ((rl - s).astype(float) * n_loci / rl + 0.5).astype(int)
        distances[idx+1:, i2] = diffs
    return distances

def get_args(args):
    parser = argparse.ArgumentParser(description='''completeCC takes allelic profile (as in https://pubmlst.org/data/) and
work out specialised complete linkage clustering result of all the profiles in the list.''',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', '--profile', help='[INPUT; REQUIRED] name of the profile file. Can be GZIPed.',
                        required=True)
    parser.add_argument('-o', '--output',
                        help='[OUTPUT; REQUIRED] Prefix for the output files. These include a NUMPY and TEXT verions of the same clustering result',
                        required=True)
    parser.add_argument('-i', '--incremental', help='[INPUT; optional] The NUMPY version of an old clustering result',
                        default='')
    parser.add_argument('-d', '--delta',
                        help='[optional] comma delimited list of threshold (delta). All values are included by default.',
                        default=None)
    parser.add_argument('--immutable',
                        help='[optional] Use a immutable clustering system. The designations of old profiles are immutable. Faster but leads to non-optimal assignment.',
                        default=False, action='store_true')

    return parser.parse_args(args)


def hierCC(args):
    params = get_args(args)
    ot = time.time()
    profile_file, cluster_file, old_cluster = params.profile, params.output + '.completeCC.npz', params.incremental

    global mat, n_loci
    mat = pd.read_csv(profile_file, sep='\t', header=None, dtype=str).values
    allele_columns = np.array([i == 0 or (not h.startswith('#')) for i, h in enumerate(mat[0])])
    mat = mat[1:, allele_columns].astype(int)
    n_loci = mat.shape[1] - 1

    logger(
        '{0}: Loaded in allelic profiles with dimension: {1} and {2}. The first column is assumed to be type id.'.format(
            time.time() - ot, *mat.shape))

    if os.path.isfile(old_cluster):
        od = np.load(old_cluster, allow_pickle=True)
        cls = od['completeCC']

        typed = {c[0]: id for id, c in enumerate(cls) if c[0] > 0}
        if len(typed) > 0:
            logger('{0}: Loaded in {1} old completeCC assignments.'.format(time.time() - ot, len(typed)))
            mat_idx = np.array([t in typed for t in mat.T[0]])
            mat[:] = np.vstack([mat[mat_idx], mat[(mat_idx) == False]])
    else:
        typed = {}

    logger('{0}: Get pairwise distances'.format(time.time() - ot))
    pool = Pool(10)
    dist = get_distances(params.output, mat, pool)
    logger('{0}: Get complex linkage clustering'.format(time.time() - ot))
    cls = linkage(ssd.squareform(dist), method='complete')
    logger('{0}: Start completeCC assignments'.format(time.time() - ot))
    nodes = np.arange(dist.shape[0]*2-1, dtype=int)
    descendents = [ [i] for i in np.arange(dist.shape[0]) ] + [None for i in np.arange(dist.shape[0]-1)]
    res = np.repeat(mat.T[0], mat.shape[1]).reshape(mat.shape)
    for idx, c in enumerate(cls.astype(int)) :
        n_id = idx + dist.shape[0]
        d = sorted([int(c[0]), int(c[1])], key=lambda x:descendents[x][0])
        min_id = descendents[d[0]][0]
        descendents[n_id] = descendents[d[0]] + descendents[d[1]]
        for tgt in descendents[d[1]] :
            res[tgt, c[2]+1:] = res[min_id, c[2]+1:]
    np.savez_compressed(cluster_file, completeCC=res)

    if not params.delta:
        with uopen(params.output + '.completeCC.gz', 'w') as fout:
            fout.write('#ST_id\t{0}\n'.format('\t'.join(['d' + str(id) for id in np.arange(n_loci)])))
            for r in res[np.argsort(res.T[0])]:
                fout.write('\t'.join([str(rr) for rr in r]) + '\n')
    else:
        deltas = map(int, params.delta.split(','))
        with uopen(params.output + '.completeCC.gz', 'w') as fout:
            fout.write('#ST_id\t{0}\n'.format('\t'.join(['d' + str(id) for id in deltas])))
            for r in res[np.argsort(res.T[0])]:
                fout.write('\t'.join([str(r[id + 1]) for id in deltas]) + '\n')
    del res
    logger('NUMPY clustering result (for incremental completeCC): {0}.completeCC.npz'.format(params.output))
    logger('TEXT  clustering result (for visual inspection): {0}.completeCC.gz'.format(params.output))


mat, n_loci = None, None
if __name__ == '__main__':
    hierCC(sys.argv[1:])

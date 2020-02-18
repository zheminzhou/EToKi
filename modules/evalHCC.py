import sys, os, pandas as pd, numpy as np, numba as nb
import argparse
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import silhouette_score, normalized_mutual_info_score

try:
    from configure import uopen, logger
except :
    from .configure import uopen, logger
    

def shannon_index(cluster) :
    logger('Calculating Shannon Index...')
    h = np.zeros(cluster.shape[1])
    for cid, c in enumerate(cluster.T) :
        cls_cnt = np.unique(c, return_counts=True)[1]/float(c.size)
        h[cid] = -np.sum(cls_cnt*np.log(cls_cnt))/np.log(c.size)
    return h

def get_similarity2(data) :
    method, cc1, cc2 = data
    if np.unique(cc1).size == 1 and  np.unique(cc1).size == 1 :
        return 1.
    return eval(method)(cc1, cc2)

def get_similarity(method, cluster, stepwise) :
    logger('Calculating similarities...')
    similarity = np.ones([cluster.shape[1], cluster.shape[1]], dtype=np.float64)
    for i1, cc1 in enumerate(cluster.T) :
        if i1 % 10 == 0 :
            logger('    similarities between level {0} and greater levels'.format(i1 * stepwise))
        similarity[i1, i1+1:] = pool.map(get_similarity2, [ [method, cc1, cc2] for cc2 in cluster.T[i1+1:] ])
        similarity[i1+1:, i1] = similarity[i1, i1+1:]
    similarity[similarity>1.] = 1.
    similarity[similarity<0.001] = 0.001
    return similarity

def get_silhouette(prefix, profile, cluster, stepwise) :
    logger('Calculating pairwise distance ...')
    np.save(prefix+'.profile.npy', profile)
    indices = np.array([[profile.shape[0]*(v/10.)**2+0.5, profile.shape[0]*((v+1)/10.)**2+0.5] for v in np.arange(10, dtype=float)], dtype=int)
    #subfiles = list(map(parallel_distance, [['evalHCC.profile.npy', 'evalHCC.dist.{0}.npy', idx] for idx in indices]))
    subfiles = pool.map(parallel_distance, [[prefix+'.profile.npy', prefix+'.dist.{0}.npy', idx] for idx in indices])
    prof_dist = np.hstack([ np.load(subfile) for subfile in subfiles ])
    prof_dist += prof_dist.T
    np.save(prefix+'.dist.npy', prof_dist)
    for subfile in subfiles :
        os.unlink(subfile)

    logger('Calculating Silhouette score ...')
    silhouette = np.array(pool.map(get_silhouette2, [ [prefix+'.dist.npy', tag] for tag in cluster.T ]))
    for subfile in [prefix+'.profile.npy', prefix+'.dist.npy'] :
        os.unlink(subfile)
    return silhouette

def get_silhouette2(data) :
    dist_file, tag = data
    s = np.unique(tag).size
    if 2 <= s < tag.shape[0] :
        return silhouette_score(np.load(dist_file), tag, metric = 'precomputed')
    else :
        return 0.


def parallel_distance(callup) :
    prof_file, sub_prefix, index_range = callup
    profiles = np.load(prof_file)

    res = profile_distance(profiles, index_range)
    subfile = sub_prefix.format(index_range[0])
    np.save(subfile, res)
    return subfile

def profile_distance(profiles, index_range=None) :
    if index_range is None :
        index_range = [0, profiles.shape[0]]

    presences = (profiles > 0)
    distances = np.zeros(shape=[profiles.shape[0], index_range[1] - index_range[0]])
    for i2, id in enumerate(np.arange(*index_range)) :
        profile, presence = profiles[id], presences[id]
        comparable = (presences[id+1:] * presence)
        diffs = (np.sum((profiles[id+1:] != profile) & comparable, axis=1).astype(float) + .5) / (np.sum(comparable, axis=1) + 1.)
        distances[id+1:, i2] = diffs
        #distances[id, :i2] = diffs[index_range[0]:index_range[0]+id]
    distances = -np.log(1-distances)
    return distances

def get_args(args) :
    parser = argparse.ArgumentParser(description='''evalHCC evaluates HierCC results using varied statistic summaries.''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', '--profile', help='[INPUT; REQUIRED] name of the profile file. Can be GZIPed.', required=True)
    parser.add_argument('-c', '--cluster', help='[INPUT; REQUIRED] name of the hierCC file. Can be GZIPed.', required=True)
    parser.add_argument('-o', '--output', help='[OUTPUT; REQUIRED] Prefix for the output files.', required=True)
    parser.add_argument('-s', '--stepwise', help='[DEFAULT: 10] Evaluate every <stepwise> levels.', default=10, type=int)

    return parser.parse_args(args)    


def evalHCC(args) :
    args = get_args(args)

#def evalHCC(profile, cluster, stepwise, ave_gene_length=1000.) :
    with uopen(args.profile) as fin :
        logger('Loading profiles ...')                
        profile_header = fin.readline().strip().split('\t')
        ST_col = np.where([p.find('#ST')>=0 for p in profile_header])[0].tolist()
        if len(ST_col) <= 0 :
            ST_col = [0]
        cols = ST_col + np.where([not h.startswith('#') for h in profile_header])[0].tolist()
        profile = pd.read_csv(fin, sep='\t', header=None, index_col=0, usecols=cols)
        profile_names = profile.index.values
        profile = profile.values
    
    with uopen(args.cluster) as fin :
        logger('Loading hierCC ...')                        
        cluster_header = fin.readline().strip().split('\t')
        cols = [0] + np.where([not h.startswith('#') for h in cluster_header])[0].tolist()
        cluster = pd.read_csv(fin, sep='\t', header=None, index_col=0, usecols=cols)
        cluster_names = cluster.index.values
        cluster = cluster.values
        s = np.arange(0, cluster.shape[1], args.stepwise)
        cluster = cluster[:, s]

    presence = np.in1d(cluster_names, profile_names)
    cluster, cluster_names = cluster[presence], cluster_names[presence]
    order = {n:id for id, n in enumerate(cluster_names)}
    profile_order = np.array([ [id, order[n]] for id, n in enumerate(profile_names) if n in order ])
    profile_order = profile_order[np.argsort(profile_order.T[1]), 0]
    profile_names = profile_names[profile_order]
    profile = profile[profile_order]
    
    with open(args.output, 'w') as fout :
        shannon = shannon_index(cluster)
        similarity = get_similarity('adjusted_rand_score', cluster, args.stepwise)

        levels = [ id*args.stepwise for id in range(len(shannon)) ]
        for level, s in zip(levels, shannon) :
            fout.write('#Shannon\t{0}\t{1:.3f}\n'.format(level, np.abs(s)))
        fout.write('\n\n#ARI\t{0}\n'.format('\t'.join([str(lvl) for lvl in levels])))
        for level, sim in zip(levels, similarity) :
            fout.write('{0}\t{1}\n'.format(level, '\t'.join([ '{0:.3f}'.format(s) for s in sim ])))
        fout.write('\n\n')
        silhouette = get_silhouette(args.output, profile, cluster, args.stepwise)
        for level, s in zip(levels, silhouette) :
            fout.write('#Silhouette\t{0}\t{1:.3f}\n'.format(level, np.abs(s)))

    #np.savez_compressed('evalHCC.npz', shannon=shannon, similarity=similarity, silhouette=silhouette)
    #logger('Done. Results saved in evalHCC.npz')


import multiprocessing
pool = multiprocessing.Pool(10)
if __name__ == '__main__' :
    evalHCC(sys.argv[1:])

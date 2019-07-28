import sys, os, pandas as pd, numpy as np, numba as nb
from sklearn.metrics.cluster import adjusted_rand_score
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

def get_silhouette(profile, cluster, stepwise, ave_gene_length=1.) :
    logger('Calculating pairwise distance ...')
    np.save('evalHCC.profile.npy', profile)
    indices = np.array([[profile.shape[0]*(v/10.)**2+0.5, profile.shape[0]*((v+1)/10.)**2+0.5] for v in np.arange(10, dtype=float)], dtype=int)
    #subfiles = list(map(parallel_distance, [['evalHCC.profile.npy', 'evalHCC.dist.{0}.npy', ave_gene_length, idx] for idx in indices]))
    subfiles = pool.map(parallel_distance, [['evalHCC.profile.npy', 'evalHCC.dist.{0}.npy', ave_gene_length, idx] for idx in indices])
    prof_dist = np.hstack([ np.load(subfile) for subfile in subfiles ])
    prof_dist += prof_dist.T
    np.save('evalHCC.dist.npy', prof_dist)
    logger('Calculating Silhouette score ...')
    silhouette = np.array(pool.map(get_silhouette2, [ ['evalHCC.dist.npy', tag] for tag in cluster.T ]))
    for subfile in subfiles + ['evalHCC.profile.npy', 'evalHCC.dist.npy'] :
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
    prof_file, sub_prefix, ave_gene_length, index_range = callup
    profiles = np.load(prof_file)

    res = profile_distance(profiles, ave_gene_length, index_range)
    subfile = sub_prefix.format(index_range[0])
    np.save(subfile, res)
    return subfile

def profile_distance(profiles, ave_gene_length, index_range=None) :
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
    distances = (1-(1-distances)**(1./ave_gene_length))*ave_gene_length
    return distances


def evaluate(profile, cluster, stepwise, ave_gene_length=1000.) :
    with uopen(profile) as fin :
        logger('Loading profiles ...')                
        profile_header = fin.readline().strip().split('\t')
        ST_col = np.where([p.find('#ST')>=0 for p in profile_header])[0].tolist()
        if len(ST_col) <= 0 :
            ST_col = [0]
        cols = ST_col + np.where([not h.startswith('#') for h in profile_header])[0].tolist()
        profile = pd.read_csv(fin, sep='\t', header=None, index_col=0, usecols=cols)
        profile_names = profile.index.values
        profile = profile.values
    
    with uopen(cluster) as fin :
        logger('Loading hierCC ...')                        
        cluster_header = fin.readline().strip().split('\t')
        cols = [0] + np.where([not h.startswith('#') for h in cluster_header])[0].tolist()
        cluster = pd.read_csv(fin, sep='\t', header=None, index_col=0, usecols=cols)
        cluster_names = cluster.index.values
        cluster = cluster.values
        s = np.arange(0, cluster.shape[1], stepwise)
        cluster = cluster[:, s]

    presence = np.in1d(cluster_names, profile_names)
    cluster, cluster_names = cluster[presence], cluster_names[presence]
    order = {n:id for id, n in enumerate(cluster_names)}
    profile_order = np.array([ [id, order[n]] for id, n in enumerate(profile_names) if n in order ])
    profile_order = profile_order[np.argsort(profile_order.T[1]), 0]
    profile_names = profile_names[profile_order]
    profile = profile[profile_order]
    
    shannon = shannon_index(cluster)

    similarity = get_similarity('adjusted_rand_score', cluster, stepwise)

    silhouette = get_silhouette(profile, cluster, stepwise, ave_gene_length)

    np.savez_compressed('evalHCC.npz', shannon=shannon, similarity=similarity, silhouette=silhouette)
    logger('Done. Results saved in evalHCC.npz')


import multiprocessing
pool = multiprocessing.Pool(10)
if __name__ == '__main__' :
    profile, cluster, stepwise = sys.argv[1:]
    evaluate(profile, cluster, int(stepwise))
import pandas as pd, numpy as np, sys, os, gzip
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import v_measure_score, silhouette_score, f1_score
import matplotlib.pyplot as plt
try:
    from configure import uopen
except :
    from .configure import uopen
    

def shannon_index(cls) :
    h = np.zeros(cls.shape[0])
    for cid, c in enumerate(cls) :
        cls_cnt = np.unique(c, return_counts=True)[1]/float(c.size)
        h[cid] = -np.sum(cls_cnt*np.log(cls_cnt))/np.log(c.size)
    return h

def cluster_distance_AR(cluster) :
    dist = np.ones([cluster.shape[1], cluster.shape[1]], dtype=np.float64)
    for i1, cc1 in enumerate(cluster.T) :
        for i2 in np.arange(i1+1, cluster.shape[1]) :
            cc2 = cluster.T[i2]
            dist[i2, i1] = dist[i1, i2] = adjusted_rand_score(cc1, cc2)
    return dist

def cluster_distance_V(cluster) :
    dist = np.zeros([cluster.shape[1], cluster.shape[1]], dtype=np.float64)
    for i1, cc1 in enumerate(cluster.T) :
        for i2 in np.arange(i1+1, cluster.shape[1]) :
            cc2 = cluster.T[i2]
            dist[i2, i1] = dist[i1, i2] = 1-f1_score(cc1, cc2, average='micro')
    return dist

def profile_distance(profile, genes=None) :
    if genes is None :
        genes = np.ones(profile.shape[1], dtype=bool)
    dist = np.zeros([ profile.shape[0], profile.shape[0] ], dtype=np.float64)
    pp = (profile[:, genes] > 0)
    for id, p in enumerate(profile[:, genes]) :
        ppp = pp[id]
        d = (p != profile[:id][:, genes]) * (pp[:id] * ppp)
        dist[:id, id] = dist[id, :id] = (1 - (1 - (np.sum(d, 1).astype(float)+.5)/(np.sum(pp[:id] * ppp, 1)+1.))**0.001)*1000.
        #dist[:id, id] = dist[id, :id] = np.sum(d, 1).astype(float)/np.sum(pp[:id] * ppp, 1)
    return dist

options = 'shannon,f1,silhouette'
if __name__ == '__main__' :
    options, profile, cluster = sys.argv[1:4]
    options = options.split(',')
    stepwise = int(sys.argv[4])
    with uopen(profile) as fin :
        profile_header = fin.readline().strip().split('\t')
        cols = [0] + np.where([not h.startswith('#') for h in profile_header])[0].tolist()
        profile = pd.read_csv(fin, sep='\t', header=None, index_col=0, usecols=cols)
        profile_names = profile.index.values
        profile = profile.values#.astype(int)
    
    with uopen(cluster) as fin :
        cluster_header = fin.readline().strip().split('\t')
        cols = [0] + np.where([not h.startswith('#') for h in cluster_header])[0].tolist()
        cluster = pd.read_csv(fin, sep='\t', header=None, index_col=0, usecols=cols)
        cluster_names = cluster.index.values
        cluster = cluster.values#.astype(int)
        stepwise = np.arange(0, cluster.shape[1], stepwise)
        cluster = cluster[:, stepwise]

    presence = np.in1d(cluster_names, profile_names)
    cluster, cluster_names = cluster[presence], cluster_names[presence]
    order = {n:id for id, n in enumerate(cluster_names)}
    profile_order = np.array([ [id, order[n]] for id, n in enumerate(profile_names) if n in order ])
    profile_order = profile_order[np.argsort(profile_order.T[1]), 0]
    profile_names = profile_names[profile_order]
    profile = profile[profile_order]
    
    np.save('source_profile', profile)
    np.save('source_cluster', cluster)
    np.save('names', profile_names)
    prof_dist = None
    if 'shannon' in options :
        plt.figure()
        h = shannon_index(cluster.T)
        np.save('shannon', h)
        plt.scatter(stepwise, h, lw=2)
    if 'f1' in options :
        plt.figure()
        dist = cluster_distance_V(cluster)
        np.save('f1_dist', dist)
        dist[dist>1.] = 1.
        dist[dist<0.001] = 0.001
        dist2 = np.log10(dist)
        heatmap = plt.imshow(dist2, cmap='hot')
        plt.colorbar(heatmap)
    if 'silhouette' in options :
        prof_dist = profile_distance(profile)
        silhouette = []
        for s, tag in zip(stepwise, cluster.T) :
            if np.unique(tag).size < profile.shape[0] and np.unique(tag).size > 1 :
                silhouette.append([s, silhouette_score(prof_dist, tag, metric='precomputed')])
        np.save('silhouette', silhouette)
        plt.figure()
        silhouette = np.array(silhouette)
        plt.scatter(silhouette.T[0], silhouette.T[1], lw=2)
    plt.show()
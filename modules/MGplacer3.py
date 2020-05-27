# MGplacer
import sys, os, numpy as np, pandas as pd, subprocess, re
from scipy.stats import poisson
from ete3 import Tree
import theano, theano.tensor as t
import click, pymc3 as pm

try:
    from configure import uopen, asc2int, xrange
except:
    from .configure import uopen, asc2int, xrange


# read ancestral state
def readAncestral(fname, nSample):
    sites = []
    data = []
    conv = np.bincount([65, 67, 67, 71, 71, 71, 84, 84, 84, 84]) - 1
    sys.stderr.write('Start reading Matrix: \n')
    with uopen(fname) as fin:
        for line in fin:
            if line.startswith('##'):
                pass
            elif line.startswith('#'):
                header = line.strip().split('\t')
                nodes = np.array([not h.startswith('#') for h in header], dtype=bool)
                nodeNames = np.array(header)[nodes]
                break
        for i, mat in enumerate(pd.read_csv(fin, sep='\t', header=None, dtype=str, chunksize=20000)):
            sys.stderr.write('Reading Matrix - Line: {0}      \r'.format(i * 20000))
            mat = mat.values
            mat = mat[np.in1d(mat[:, np.where(nodes)[0][0]], ['A', 'C', 'G', 'T'])]
            data.append(conv[mat[:, nodes].astype(bytes).view(np.uint8)])
            for id, (cont, site) in enumerate(mat[:, :2]):
                sites.append([(cont, int(site)), ] + [1.] * nSample)
    sys.stderr.write('Read Matrix DONE. Total SNP sites: {0}                \n'.format(len(sites)))
    return np.array(sites), nodeNames, np.ascontiguousarray(np.vstack(data).T, dtype=np.uint8)


# read bam file
def parseBAMs(fnames, sites):
    samtools = 'samtools'
    conv = np.bincount([65, 67, 67, 71, 71, 71, 84, 84, 84, 84]) - 1
    knownSamples = np.zeros([sites.shape[0], len(fnames), 4])

    from tempfile import NamedTemporaryFile
    tmpFile = NamedTemporaryFile(delete=False, dir='.')
    siteMap = {}
    for i, s in enumerate(sites):
        tmpFile.write('{0}\t{1}\n'.format(*s[0]).encode('utf-8'))
        siteMap[s[0]] = i
    tmpFile.close()

    sys.stderr.write('Reading BAM file(s). \n'.format(len(sites)))
    for fnId, fname in enumerate(fnames):
        p = subprocess.Popen('{0} mpileup -l{2} {1}'.format(samtools, fname, tmpFile.name).split(),
                             universal_newlines=True, stdout=subprocess.PIPE)
        for line in p.stdout:
            cont, site, _, _, bases = line.strip('\n').split('\t')[:5]
            j = siteMap.get((cont, int(site)), -1)
            if j > -1:
                bases = re.sub(r'\^.|\$', '', bases).upper()
                c = re.findall(r'(\w*)[+-](\d+)(\w+)', bases)
                if len(c) > 0:
                    bases = ''.join([c[0][0]] + [cc[2][int(cc[1]):] for cc in c])
                bases = conv[np.array(list(bases), dtype=bytes).view(np.int8)]
                bases = bases[bases >= 0]
                baseCounts = np.bincount(bases, minlength=4)
                knownSamples[j, fnId, :] = baseCounts
    os.unlink(tmpFile.name)

    sampleDepths = np.sum(knownSamples, 2)
    sampleMean = np.max([np.min([np.mean(sampleDepths, 0), np.median(sampleDepths, 0)], 0),
                         sampleDepths.shape[0] / np.sum(1 / (sampleDepths + 0.5), 0) - 0.5], 0)
    siteWeight = np.min([(1 - poisson.cdf(sampleMean / 2., sampleDepths)) / (
                1 - poisson.cdf(sampleMean / 2., np.max([sampleMean / 2., [1] * sampleMean.size], 0))), \
                         poisson.cdf(sampleMean * 2., sampleDepths) / \
                         poisson.cdf(sampleMean * 2., np.max([sampleMean * 2.,[1] * sampleMean.size],0))], 0)
    siteWeight[siteWeight > 1] = 1
    sites[:, 1:] = siteWeight * sampleDepths

    sampleDepths[sampleDepths <= 0] = 1
    for i in np.arange(knownSamples.shape[2]):
        knownSamples[:, :, i] *= siteWeight
    return knownSamples


def readTree(treefile, knownNodes):
    nodeMap = {n: i for i, n in enumerate(knownNodes)}

    branches = []
    tre = Tree(treefile, format=1)
    cRoot = min([c for c in tre.children if c.name in nodeMap], key=lambda c: nodeMap[c.name])
    for node in tre.traverse('postorder'):
        if node.name in nodeMap and node != cRoot:
            n, d = node, 0.
            while n.up and n.up.name not in nodeMap:
                d += n.dist
                n = n.up
            if n.up:
                branches.append([nodeMap[n.up.name], nodeMap[node.name], d + node.dist, 0])
            else:
                branches.append([nodeMap[cRoot.name], nodeMap[node.name], d + cRoot.dist, 0])
    branches = np.array(branches)
    branches.T[3] = branches.T[2] / np.sum(branches.T[2])
    return branches

max_dist = 1.

@click.command()
@click.argument('ancestralfile')
@click.argument('bamfile')
@click.argument('treefile')
@click.argument('maxgenotype', type=int, default=3)
def main(ancestralfile, bamfile, treefile, maxgenotype=3):
    bamFiles = bamfile.split(',')
    nSample = len(bamFiles)
    maxGenotype = maxgenotype

    #if not os.path.isfile('test1'):
    knownSites, knownNodes, knownMatrix = readAncestral(ancestralfile, nSample)
    branches = readTree(treefile, knownNodes)
    knownSamples = parseBAMs(bamFiles, knownSites).transpose([0, 2, 1])

    exists = (np.sum(knownSites.T[1:], 0) >= 0.01) & np.any(knownSamples > 0, (1,2))
    knownSites = knownSites[exists]
    knownMatrix = knownMatrix[:, exists]
    knownSamples = knownSamples[exists]
    for br in branches :
        diff_sites = knownMatrix[br[0].astype(int)] != knownMatrix[br[1].astype(int)]
        snv = knownMatrix[br[1].astype(int), diff_sites]
        presence = np.sum(knownSamples[diff_sites, snv], 1)
        if np.sum(presence > 0)*100 <= presence.shape[0] :
            knownMatrix[br[1].astype(int)] = 100
            branches[branches.T[0] == br[1], 0] = br[0]
            br[:2] = -1
    branches = branches[branches.T[0] > -1]
    exists = [len(set(np.unique(mat)) - {100}) > 1 for mat in knownMatrix.T]
    knownSites = knownSites[exists]
    knownMatrix = knownMatrix[:, exists]
    knownSamples = knownSamples[exists]

    snps = np.zeros(knownMatrix.shape[1])
    for br in branches:
        snps += knownMatrix[br[0].astype(int)] != knownMatrix[br[1].astype(int)]
    max_dist = np.sum(snps <= np.median(snps)) / snps.shape[0] #np.sum(snps)
    #max_dist = 1.

    m, x, y, w = [], [], [], []
    for mi, (site, mat, sample) in enumerate(zip(knownSites, knownMatrix.T, knownSamples)):
        c = [i for i in np.unique(mat) if i != 100]
        for i in c :
            n = sample[i]
            m.append(mi)
            x.append(mat == i)
            y.append(n/len(c))
            w.append(site[1:]/len(c))
    x = np.array(x, dtype=float)
    m = np.array(m, dtype=int)
    w = np.sum(np.array(w, dtype=float), 1)
    y = np.round(np.sum(np.array(y, dtype=float), 1)/w, 2)
    tmp, idx = np.unique(np.hstack([x, y[:, np.newaxis]]), axis=0, return_inverse=True)
    x0, y0, w0, m0 = tmp[:, :-1], tmp.T[-1], np.zeros(tmp.shape[0], dtype=float), [[] for i in np.arange(tmp.shape[0])]
    for i in np.arange(tmp.shape[0]) :
        m0[i] = m[idx == i]
        w0[i] = np.sum(w[idx == i])
    x, y, w, m = x0, y0, w0, m0
    sys.stderr.write(
        'Sites that are not covered by any read are ignored. {1} distinct categories of SNVs in {0} sites remain.\n'.format(np.sum(exists),
                                                                                                     x.shape[0]))

    x2 = t._shared(x)
    branches2 = t._shared(branches)

    for nGenotype in np.arange(1, maxGenotype + 1):
        sys.stderr.write(
            '\n----------\nRunning MCMC with assumption of {0} genotype(s) present in the sample.\n'.format(nGenotype))
        with pm.Model() as model:
            brs = pm.Flat('brs', shape=nGenotype if nGenotype > 1 else ())
            props2 = pm.Dirichlet('props2', a=1./np.arange(1, nGenotype+1, dtype=float)) \
                if nGenotype > 1 else \
                pm.Exponential('props2', lam=1, shape=())
            props = pm.Deterministic('props', props2*pm.math.sum(props2-0.05) + 0.05)
            genotypes = getBranchGenotype(x2, branches2, brs)
            sigma = pm.Gamma('sigma', alpha=2, beta=0.5) \
                if nGenotype == 0 else \
                pm.Gamma('sigma', alpha=0.5, beta=2)

            mu = pm.math.sum(genotypes * props, 1)
            dist = pm.math.abs_(mu - y)

            restricted_y = pm.math.clip(dist, 0, max_dist)
            lk = pm.Deterministic('lk', pm.math.sum(w * pm.Normal.dist(mu=0., sigma=sigma).logp(restricted_y)))

            rec_y = pm.math.minimum(dist, 1-dist)
            hybrid_score = pm.Deterministic('hybrid_score', pm.math.sqrt(pm.math.sum(pm.math.sqr(rec_y))/len(y)) )

            pm.Potential('likelihood', lk)

            step_br = TreeWalker(brs, branches)
            step_others = pm.step_methods.Metropolis(vars=[sigma, props])
            trace = pm.sample(progressbar=True, draws=5000, tune=15000, step=[step_br, step_others], chains=8, cores=8,
                              compute_convergence_checks=False)

        trace_logp = np.array([ np.mean([ [t['likelihood'], t['hybrid_score']] for t in strace ], 0) for strace in trace._straces.values() ])
        sys.stderr.write('Done.\n----------\n'.format(nGenotype))
        # select traces
        trace_id = np.argmax(trace_logp.T[0])
        logp = trace_logp[trace_id]

        sys.stdout.write(
            '----------\nNo. Genotypes:\t{0}\tlogp:\t{1}\thybrid_score:\t{2}\n'.format(nGenotype, logp[0], logp[1]))
        sigma = trace.get_values('sigma', chains=trace_id)
        sigma = np.sort(sigma)
        sys.stdout.write('Sigma\tMean:\t{0:.6E}\tCI95%:\t[ {1:.6E} - {2:.6E} ]\n'.format(np.mean(sigma),
                                                                                         sigma[int(sigma.size * 0.025)],
                                                                                         sigma[
                                                                                             int(sigma.size * 0.975)]))

        if nGenotype < 2:
            br_locs = trace.get_values('brs', chains=trace_id)[:, np.newaxis]
            props = trace.get_values('props', chains=trace_id)[:, np.newaxis]
        else:
            br_locs = trace.get_values('brs', chains=trace_id)
            props = trace.get_values('props', chains=trace_id)

        props /= np.sum(props, 1)[:, np.newaxis]
        for id, (br_loc, prop) in enumerate(zip(br_locs.T, props.T)):
            prop = np.sort(prop)
            if nGenotype > 0:
                brs, locs = br_loc.astype(int), br_loc % 1
                brNames, brCounts = np.unique(brs, return_counts=True)
                brCounts = brCounts.astype(float) / np.sum(brCounts)
                idx = np.argsort(-brCounts)
                brNames, brCounts = brNames[idx], brCounts[idx]
                sys.stdout.write(
                    '\tGenotype {0}:\tMean proportion:\t{1:.4f}\tCI95%:\t[ {2:.4f} - {3:.4f} ]\n'.format(id + 1,
                                                                                                         np.mean(prop), \
                                                                                                         prop[int(
                                                                                                             prop.size * 0.025)],
                                                                                                         prop[int(
                                                                                                             prop.size * 0.975)]))
                for br, cnt in zip(brNames, brCounts):
                    if cnt >= 0.01 or cnt >= 0.3 * brCounts[0]:
                        lc = np.sort(locs[brs == br])
                        sys.stdout.write(
                            '\t\t{0:.2f} %\t{1} - {2}\tLocation:\t{3:.4f}\tCI95%:\t[ {4:.4f} - {5:.4f} ]\n'.format(
                                cnt * 100, \
                                knownNodes[int(branches[br, 0])], \
                                knownNodes[int(branches[br, 1])], \
                                np.mean(lc), lc[int(lc.size * 0.025)], \
                                lc[int(lc.size * 0.975)]))
    sys.stderr.write('All DONE\n')

# no in use
def waic(trace, start=0):
    def traveseTraces(trace):
        observations = []
        for chain in trace._straces.values():
            observations.append([])
            for tr in chain:
                observations[-1].append(tr['likelihood'])
        return observations

    observations = traveseTraces(trace)
    nObs = observations[0][0].size
    res = []
    for observations2 in observations:
        lppd, var1 = np.zeros(nObs), np.zeros(nObs)
        for i in np.arange(start, nObs, 1000):
            obs = np.array([observation[i:i + 1000] for observation in observations2])
            max_obs = np.max(obs, 0)
            lppd[i:i + 1000] = np.log(np.mean(np.exp(obs - max_obs), 0)) + max_obs
            var1[i:i + 1000] = np.var(obs, 0)

        logp = [np.sum(observation) for observation in observations2]
        max_logp = np.max(logp)
        w = np.exp((logp - max_logp) / np.log(nObs))
        wbic = np.sum(w * logp) / np.sum(w)
        res.append([np.sum(lppd), np.sum(lppd) - np.sum(var1), wbic])
        x = np.mean(observations2, 0)
        # print(np.unique(x.astype(int), return_counts=True))
    return res


@theano.compile.ops.as_op(itypes=[t.dmatrix, t.dmatrix, t.dvector], otypes=[t.dmatrix])
def getBranchGenotype2(x, branches, brs):
    br2 = brs.astype(int)
    locs = brs - br2
    return x[:, branches[br2, 0].astype(int)] * (1 - locs) + x[:, branches[br2, 1].astype(int)] * locs


def getBranchGenotype(x, branches, brs):
    return getBranchGenotype2(x, branches, brs) \
        if brs.type == t.dvector \
        else getBranchGenotype1(x, branches, brs)


@theano.compile.ops.as_op(itypes=[t.dmatrix, t.dmatrix, t.dscalar], otypes=[t.dmatrix])
def getBranchGenotype1(x, branches, brs):
    br2 = brs.astype(int)
    locs = brs - br2
    return (x[:, branches[br2, 0].astype(int)] * (1 - locs) + x[:, branches[br2, 1].astype(int)] * locs)[:, np.newaxis]

class TreeWalker(pm.step_methods.Metropolis):
    def __init__(self, var, tree):
        super().__init__([var], tune=False, proposal_dist=pm.step_methods.metropolis.UniformProposal)
        self.metrop_select = pm.step_methods.arraystep.metrop_select

        self.var = var.name
        self.scaling = 100.
        self.counter = 0

        n2br = {int(nodes[1]): br for br, nodes in enumerate(tree)}
        self.tree = {int(br): [[], []] for br, nodes in enumerate(tree)}
        for br, nodes in enumerate(tree):
            neighbors = np.where((tree.T[0] == nodes[0]) & (tree.T[1] != nodes[1]))[0]
            if neighbors.size > 0:
                self.tree[br][0].extend([[0, n] for n in neighbors])
            if int(nodes[0]) in n2br:
                b2 = n2br[int(nodes[0])]
                self.tree[b2][1].append([0, br])
                self.tree[br][0].append([1, b2])

    def astep(self, br0):
        brs = np.copy(br0)
        i = self.counter % brs.size
        br, loc = int(brs[i]), brs[i] % 1
        others = np.concatenate([brs[:i], brs[i+1:]])

        moves = np.random.rand(int(self.scaling + 0.5)) * 2. - 1.
        for m in moves:
            loc += m
            while loc > 0.999999999:
                n = self.tree[br][1]
                if len(n):
                    x, br = n[np.random.choice(len(n))]
                    loc -= 0.999999999
                else:
                    loc = -1e-9
            while loc < 0:
                n = self.tree[br][0]
                if len(n):
                    x, br = n[np.random.choice(len(n))]
                    loc = 1 + loc if x else -loc
                else:
                    loc = 0
        brs[i] = br + loc
        incompatibles = set([b2[1] for b in self.tree[br] for b2 in b ] + [br])
        accept1 = set(others.astype(int).tolist()) & incompatibles
        accept = self.delta_logp(brs, br0) if len(accept1) < 1 else -999.
        br_new, accepted = self.metrop_select(accept, brs, br0)
        self.accepted += accepted
        if self.counter == 100:
            if self.accepted / float(self.counter) > 0.2:
                self.scaling = np.min([self.scaling * 1.5, 1000.])
            elif self.accepted / float(self.counter) < 0.05:
                self.scaling = np.max([self.scaling / 1.5, 1.])
            self.counter = self.accepted = 0
        self.counter += 1

        stats = {
            'tune': self.tune,
            'accept': np.exp(accept),
        }
        return br_new, [stats]


if __name__ == '__main__':
    main()

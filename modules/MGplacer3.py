# MGplacer
import sys, os, numpy as np, pandas as pd, subprocess, re, numba as nb
from scipy.stats import poisson
from ete3 import Tree
import theano, theano.tensor as t
import click, pymc3 as pm
import logging

theano.gof.compilelock.set_lock_status(False)

ch = logging.StreamHandler(sys.stdout)
ch.setFormatter(logging.Formatter('%(asctime)s %(message)s'))
ch.setLevel(logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(ch)

try:
    from configure import uopen, asc2int, xrange, externals
except:
    from .configure import uopen, asc2int, xrange, externals

conv = np.bincount([65, 67, 67, 71, 71, 71, 84, 84, 84, 84], minlength=256) - 1

# read ancestral state
def readAncestral(fname):
    '''outputs: site list, branch list, 2D mutations'''
    sites, data = [], []
    logger.warning('Start reading Matrix: ')
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
            logger.warning('Reading Matrix - Line: {0}      \r'.format(i * 20000))
            mat = mat.values
            mat = mat[np.in1d(mat[:, np.where(nodes)[0][0]], ['A', 'C', 'G', 'T'])]
            data.append(conv[np.vectorize(ord)(mat[:, nodes])])
            for id, (cont, site) in enumerate(mat[:, :2]):
                sites.append([(cont, int(site)), 1., 0.])
    logger.warning('Read Matrix DONE. Total SNP sites: {0}                '.format(len(sites)))
    return np.array(sites, dtype=object), nodeNames, np.ascontiguousarray(np.vstack(data).T, dtype=np.int8)


# read bam file
def parseBAMs(bamfile, sites):
    bamFiles = bamfile.split(',')
    ''''''
    samtools = externals['samtools']
    knownSamples = np.zeros([sites.shape[0], 4])

    from tempfile import NamedTemporaryFile
    tmpFile = NamedTemporaryFile(delete=False, dir='.')
    siteMap = {}
    for i, s in enumerate(sites):
        tmpFile.write('{0}\t{1}\n'.format(*s[0]).encode('utf-8'))
        siteMap[s[0]] = i
    tmpFile.close()

    logger.warning('Reading BAM file(s). '.format(len(sites)))
    for fnId, fname in enumerate(bamFiles):
        p = subprocess.Popen('{0} mpileup -AB -q0 -Q0 -l{2} {1}'.format(samtools, fname, tmpFile.name).split(),
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
                knownSamples[j, :] += baseCounts
    os.unlink(tmpFile.name)

    # outliers are > 1.5 IQR and >3X depth
    sampleDepths = np.sum(knownSamples, 1)
    iqr = max(np.quantile(sampleDepths, 0.75, 0) - np.quantile(sampleDepths, 0.25, 0), 2.01)

    sampleDepthRange = np.hstack([np.quantile(sampleDepths, 0.25, 0) - 1.5 * iqr, np.quantile(sampleDepths, 0.75, 0) + 1.5 * iqr])
    validValues = ((sampleDepths >= sampleDepthRange.T[0]) & (sampleDepths <= sampleDepthRange.T[1])).astype(float)
    sampleMean = np.max([1.00001, np.sum(sampleDepths*validValues, 0)/np.sum(validValues, 0)])

    siteWeight = ((sampleDepths >= 0.5 * sampleMean) & (sampleDepths <= 2. * sampleMean)).astype(float)
    siteWeight[sampleDepths > sampleMean*2.] = poisson.pmf(sampleDepths[sampleDepths > sampleMean*2.], sampleMean*2.) / \
        poisson.pmf(np.max([np.round(sampleMean*2.), 1]), sampleMean*2.)
    siteWeight[sampleDepths < sampleMean/2.] = poisson.pmf(sampleDepths[sampleDepths < sampleMean/2.], sampleMean/2.) / \
        poisson.pmf(np.max([np.round(sampleMean/2.), 1]), sampleMean/2.)
    siteWeight[sampleDepths <= 0] = 0.
    sites[:, 1] = siteWeight

    return knownSamples


def readTree(treefile, knownNodes):
    '''output: branches[start_node, end_node, length, proportion]'''
    nodeMap, in_use = {n: i for i, n in enumerate(knownNodes)}, {n:False for n in knownNodes}
    branches = []
    tre = Tree(treefile, format=1)
    cRoot = min([c for c in tre.children if c.name in nodeMap], key=lambda c: nodeMap[c.name])
    for node in tre.traverse('postorder'):
        if node.name in nodeMap and node != cRoot:
            in_use[node.name] = True
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

def eval(genotypes, offsets, props) :
    indices = np.arange(genotypes.size).reshape(genotypes.shape)
    indices = np.array([i.ravel() for i in np.meshgrid(*indices)]).T
    combinations = np.take(genotypes, indices)
    freqs = np.take(np.array([1-offsets, offsets]).T, indices).prod(1)
    combinations, reps = np.unique(combinations, axis=0, return_inverse=True)
    freqs = np.bincount(reps, freqs)
    combinations = np.array([np.bincount(c, props, 4) for c in combinations])
    return combinations, freqs

def remove_unrelavant_branches(branches, knownSites, knownMatrix, knownSamples, to_del) :
    # remove branches that are not changing
    for td, br in zip(to_del, branches) :
        diff_sites = knownMatrix[br[0].astype(int)] != knownMatrix[br[1].astype(int)]
        snv = knownMatrix[br[1].astype(int), diff_sites]
        presence = knownSamples[diff_sites, snv] > 0
        if td or np.count_nonzero(presence) == 0 or np.sum(knownSites[diff_sites, 1]*presence)*100 <= np.sum(knownSites[diff_sites, 1]) :
            knownMatrix[br[1].astype(int)] = 100
            branches[branches.T[0] == br[1], 0] = br[0]
            br[:2] = -1
    branches = branches[branches.T[0] > -1]
    exists = [len(set(np.unique(mat)) - {100}) > 1 for mat in knownMatrix.T]
    knownSites = knownSites[exists]
    knownMatrix = knownMatrix[:, exists].astype(np.int8)
    knownSamples = knownSamples[exists]
    return branches, knownSites, knownMatrix, knownSamples

@click.command()
@click.option('-a', '--ancestralfile', help='ancestral file generated by EToKi.py phylo --task ancestral')
@click.option('-b', '--bamfile', help='read mapping results in BAM format')
@click.option('-t', '--treefile', help='the labelled tree that is generated together with the ancestral file')
@click.option('-n', '--maxgenotype', help='Assumed maximum number of genotypes in the short reads. Default: 2', type=int, default=2)
@click.option('-m', '--mingenotype', help='Assumed minimum number of genotypes in the short reads. Default: 1', type=int, default=1)
@click.option('-r', '--homoplasy', help='Assumed homoplastic ratio in the dataset. Default: 0.', type=float, default=0.)
@click.option('-x', '--rmbranch', help='Remove some branches from the dataset. Default: None', default='')
def main(ancestralfile, bamfile, treefile, mingenotype, maxgenotype, homoplasy, rmbranch):
    global rec_level
    rec_level = homoplasy
    maxGenotype = maxgenotype
    # prepare data: noticing sites, tree, and sample traces
    knownSites, knownNodes, knownMatrix = readAncestral(ancestralfile)
    branches = readTree(treefile, knownNodes)
    knownSamples = parseBAMs(bamfile, knownSites)

    # remove unimportant sites
    exists = knownSites.T[1] >= 0.001
    knownSites = knownSites[exists]
    knownMatrix = knownMatrix[:, exists]
    knownSamples = knownSamples[exists]

    # update weights with read depths
    knownSites.T[2] = np.sum(knownSamples, 1)
    knownSites.T[1] *= knownSites.T[2]
    knownSamples /= np.sum(knownSamples, 1)[:, np.newaxis]

    # arrange matrix such that same patterns are together
    for m, s in zip(knownMatrix.T, knownSamples) :
        tr = np.argsort(-np.bincount(m, minlength=4))
        s[:] = s[tr]
        m[:] = np.argsort(tr)[m]
    idx = np.lexsort( np.vstack([knownMatrix, knownSamples.T]) )
    knownSites = knownSites[idx]
    knownMatrix = knownMatrix[:, idx]
    knownSamples = knownSamples[idx]

    # remove unrelavant branches
    to_del = np.in1d(branches.T[1], np.where(np.in1d(knownNodes, rmbranch.split(',')))[0])
    branches, knownSites, knownMatrix, knownSamples = remove_unrelavant_branches(branches, knownSites, knownMatrix, knownSamples, to_del)
    weights = knownSites.T[1:].astype(float)

    # identify uniq patterns
    matrix_pattern, sample_pattern = 0, 0
    patterns = np.zeros(knownSites.shape[0], dtype=np.int8)
    for i in np.arange(1, patterns.size) :
        if np.all(knownMatrix.T[i] == knownMatrix.T[i-1]) :
            if np.all(knownSamples[i] == knownSamples[i-1]) :
                patterns[i] = 2
            else :
                patterns[i] = 1
                sample_pattern += 1
        else :
            matrix_pattern += 1
    logger.warning('Found {0} SNP patterns and {1} read patterns from {2} SNP sites'.format(
        matrix_pattern, sample_pattern+matrix_pattern, knownSites.shape[0]))

    brMatrix = np.array([knownMatrix[int(br[0])] + knownMatrix[int(br[1])]*4 for br in branches ], dtype=np.int8)
    patterns2 = t._shared(patterns)
    weights2 = t._shared(weights)
    brMatrix2 = t._shared(brMatrix)
    knownSamples2 = t._shared(knownSamples)
    inferredHeterogeneity = np.zeros(2, dtype=float)
    inferredHeterogeneity2 = t._shared(inferredHeterogeneity)

    init_vals = {'brs':np.array(np.random.rand(mingenotype)*branches.shape[0]),
                 'props_raw':np.ones(mingenotype)/mingenotype,
                 'sigma':np.array(1.)}
    for nGenotype in np.arange(mingenotype, maxGenotype + 1):
        logger.warning(
            '----------  Running MCMC with assumption of {0} genotype(s) present in the sample.'.format(nGenotype))

        with pm.Model() as model:
            brs = pm.Flat('brs', shape=nGenotype if nGenotype > 1 else (),
                          testval=init_vals['brs']
                          if nGenotype > 1 else init_vals['brs'][0] )
            props_raw = pm.Dirichlet('props_raw', a=init_vals['props_raw'], testval=init_vals['props_raw'] ) \
                if nGenotype > 1 \
                else pm.DiscreteUniform('props_raw', upper=1, lower=1, testval=init_vals['props_raw'][0] )

            props = pm.Deterministic('props', props_raw*(1.-0.0*nGenotype) + 0.0)
            sigma = pm.Gamma('sigma', alpha=1., beta=1., testval=init_vals['sigma'] )

            lk = pm.Deterministic('lk', getGenotypesAndLK(inferredHeterogeneity2, brMatrix2, \
                                                          weights2, knownSamples2, \
                                                          sigma, brs, props, patterns2))
            pm.Deterministic('hetero_dist', inferredHeterogeneity2)
            pm.Potential('likelihood', lk)

            step_br = TreeWalker(brs, branches)
            step_others = pm.step_methods.Metropolis(vars=[sigma, props_raw])
            trace = pm.sample(progressbar=True, draws=5000, tune=5000*nGenotype,
                              step=[step_br, step_others], chains=8, cores=8, return_inferencedata=False,
                              compute_convergence_checks=False)

        trace_logp = np.array([ np.median([ t['likelihood'] for t in strace ], 0) for strace in trace._straces.values() ])
        logger.warning('Done.'.format(nGenotype))
        # select traces
        trace_id = np.argmax(trace_logp)

        if nGenotype <= 1 :
            ls = {'brs':np.array([trace.get_values('brs', chains=trace_id)[-1]]),
                   'props_raw':np.array([trace.get_values('props_raw', chains=trace_id)[-1]]),
                   'sigma':trace.get_values('sigma', chains=trace_id)[-1]}
        else :
            ls = {'brs':trace.get_values('brs', chains=trace_id)[-1],
                   'props_raw':trace.get_values('props_raw', chains=trace_id)[-1],
                   'sigma':trace.get_values('sigma', chains=trace_id)[-1]}
        init_vals = {'brs':np.concatenate([ ls['brs'], np.random.rand(1)*branches.shape[0] ]),
                       'props_raw':np.concatenate([ls['props_raw'], [0.5]])/(0.5+np.sum(ls['props_raw'])),
                       'sigma':ls['sigma']}

        logp = trace_logp[trace_id]
        lk = trace.get_values('likelihood', chains=trace_id)
        lk = np.sort(lk)
        ad = trace.get_values('hetero_dist', chains=trace_id).T[0]
        ad = np.sort(ad)
        hd = trace.get_values('hetero_dist', chains=trace_id).T[1]
        hd = np.sort(hd)
        logger.info('----------')
        logger.info(
            'No. Genotypes:\t{0}\tlogp:\t{1:.2f} [ {2:.2f} - {3:.2f} ]\thybrid_score:\t{4:.8f}\t{5:.8f}'.format(
                nGenotype, logp, lk[int(lk.size * 0.025)], lk[int(lk.size * 0.975)],
                np.median(ad), np.median(hd) ))

        sigma = trace.get_values('sigma', chains=trace_id)
        sigma = np.sort(sigma)
        logger.info(
            'Sigma\tMean:\t{0:.4E}\tCI95%:\t[ {1:.4E} - {2:.4E} ]'.format(
                np.median(sigma), sigma[int(sigma.size * 0.025)], sigma[int(sigma.size * 0.975)]))

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
                logger.info(
                    '\tGenotype {0}:\tMean proportion:\t{1:.4f}\tCI95%:\t[ {2:.4f} - {3:.4f} ]'.format(
                        id + 1, np.median(prop), prop[int(prop.size * 0.025)], prop[int(prop.size * 0.975)]))

                for br, cnt in zip(brNames, brCounts):
                    if cnt >= 0.01 or cnt >= 0.3 * brCounts[0]:
                        lc = np.sort(locs[brs == br])
                        logger.info(
                            '\t\t{0:.2f} %\t{1} - {2}\tLocation:\t{3:.4f}\tCI95%:\t[ {4:.4f} - {5:.4f} ]'.format(
                                cnt * 100,
                                knownNodes[int(branches[br, 0])], knownNodes[int(branches[br, 1])],
                                np.median(lc), lc[int(lc.size * 0.025)], lc[int(lc.size * 0.975)]))
    logger.info('All DONE')

@nb.njit(fastmath=True, cache=True)
def greedyType(genotype, branchEnd, loc, gN, cache) :
    genotype[:, :] = 0.
    for gId in range(gN) :
        p = 1.
        xId = 0
        for bId in range(branchEnd.shape[0]) :
            lId = (gId >> bId) & 1
            be = (branchEnd[bId] >> (lId << 1)) & 3
            cache[bId] = be
            xId = (xId << 2) + be

            if lId == 1 :
                p *= loc[bId]
            else :
                p *= 1. - loc[bId]
        genotype[xId, -3] += p
        genotype[xId, :branchEnd.shape[0]] = cache[:branchEnd.shape[0]]
    combs = np.sort(branchEnd)
    n_comb = 1
    for i in np.arange(1, combs.shape[0]):
        if combs[i] != combs[i - 1]:
            n_comb += 1
    return n_comb

rec_level = 0.
@nb.njit(fastmath=True, cache=True)
def gibbsType(encodeType, sample, props, sigma, w, r, stage, cache, n_comb) :
    ret = -1
    if stage <= 1 :
        p2 = np.log(1./n_comb)
        for gId in range(encodeType.shape[0]) :
            if encodeType[gId, -3] > 0 :
                cache[:] = 0.
                for bId in range(props.size) :
                    cache[int(encodeType[gId, bId])] += props[bId]
                encodeType[gId, -5] = np.sum(np.abs(cache-sample)) #np.sqrt(np.sum(np.square(cache-sample)))
                encodeType[gId, -4] = np.sum(np.abs(np.sort(-cache)-np.sort(-sample))) #np.sqrt(np.sum(np.square(np.sort(-cache)-np.sort(-sample))))
                p = p2 + np.log(1./(sigma*2.506628275)) \
                    - .5*((1-rec_level)*(encodeType[gId, -5]**2) + rec_level*(encodeType[gId, -4]**2))/(sigma**2.)
                p = w[0] * p
                encodeType[gId, -2] = p + np.log(encodeType[gId, -3])
            else :
                encodeType[gId, -2] = -2147483647.
        encodeType[:, -1] = np.exp(encodeType[:, -2] - np.max(encodeType[:, -2]))
        encodeType[:, -1] = encodeType[:, -1]/np.sum(encodeType[:, -1])
    for gId in range(encodeType.shape[0]) :
        if encodeType[gId, -1] > 0 :
            if r <= encodeType[gId, -1] :
                ret = gId
                break
            else :
                r -= encodeType[gId, -1]
    return encodeType[ret, -2], np.square(encodeType[ret, -5])*w[0], (np.square(encodeType[ret, -4])*w[0] if w[1] > 1 else 0.)

@nb.njit(fastmath=True, cache=True)
def getTypes(inferredHeterogeneity, branchEnds, knownSamples, weights, loc, props, sigma, rands, pattern) :
    cache = np.zeros((4,), dtype=np.float64)
    encodeType = np.zeros((4 ** branchEnds.shape[1], branchEnds.shape[1] + 5),
                          dtype=np.float64)  # 2 bits for each genotype. A:00 C:01 G:10 T:11
    inferredHeterogeneity[:] = 0.
    gN = 2**branchEnds.shape[1] # two (ancestral or derived) possibles for each genotype
    lglk = 0.
    for i in range(branchEnds.shape[0]) :
        p = pattern[i]
        if i > 0 and p < 2 :
            if p == 0 and np.all(branchEnds[i-1] == branchEnds[i]) :
                p = 1
            if p == 1 and np.all(knownSamples[i-1] == knownSamples[i]) :
                p = 2

        # list all possible nucleotide combinations
        if p == 0 :
            n_comb = greedyType(encodeType, branchEnds[i], loc, gN, cache)
        lglk2, ih1, ih2 = gibbsType(encodeType, knownSamples[i], props, sigma, weights.T[i], rands[i], p, cache, n_comb)
        lglk += lglk2
        inferredHeterogeneity[0] += ih1
        inferredHeterogeneity[1] += ih2
    inferredHeterogeneity[:] = np.sqrt(inferredHeterogeneity/np.sum(weights[0]))
    return lglk

@theano.compile.ops.as_op(itypes=[t.dvector, t.bmatrix, t.dmatrix, t.dmatrix, t.dscalar, t.dvector, t.dvector, t.bvector], otypes=[t.dscalar])
def getGenotypesAndLK2(inferredHeterogeneity, brMatrix, weights, knownSamples, sigma, brs, props, pattern) :
    br = brs.astype(int)
    loc = brs - br
    branchEnds = brMatrix[br].T
    lk = getTypes(inferredHeterogeneity, branchEnds, knownSamples,
             weights, loc, props, sigma,
             np.random.rand(branchEnds.shape[0]), pattern)
    return np.array(lk)

@theano.compile.ops.as_op(itypes=[t.dvector, t.bmatrix, t.dmatrix, t.dmatrix, t.dscalar, t.dscalar, t.dscalar, t.bvector], otypes=[t.dscalar])
def getGenotypesAndLK1(inferredHeterogeneity, brMatrix, weights, knownSamples, sigma, brs, props, pattern) :
    br = brs.reshape(1).astype(int)
    loc = brs.reshape(1) - br
    props = np.array([props])
    branchEnds = brMatrix[br].T
    lk = getTypes(inferredHeterogeneity, branchEnds, knownSamples,
             weights, loc, props, sigma,
             np.random.rand(branchEnds.shape[0]), pattern)
    return np.array(lk)

def getGenotypesAndLK(inferredHeterogeneity, brMatrix, weights, knownSamples,
                      sigma, brs, props, patterns) :
    return getGenotypesAndLK2(inferredHeterogeneity, brMatrix, weights, knownSamples,
                              sigma, brs, props, patterns) \
        if brs.type == t.dvector \
        else getGenotypesAndLK1(inferredHeterogeneity, brMatrix, weights, knownSamples,
                                sigma, brs, props, patterns)

@theano.compile.ops.as_op(itypes=[t.bvector, t.bvector, t.dmatrix], otypes=[t.dscalar])
def getHeteroDist(inferredHeterogeneity, knownHeterogeneity, weight) :
    hd = np.sum(np.abs(knownHeterogeneity - inferredHeterogeneity) * weight.T[0])/np.sum(weight.T[0])
    return np.array(hd)

class TreeWalker(pm.step_methods.Metropolis):
    def __init__(self, var, tree):
        super().__init__([var], tune=False, proposal_dist=pm.step_methods.metropolis.UniformProposal)
        self.metrop_select = pm.step_methods.arraystep.metrop_select

        self.var = var.name
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
        if self.counter == 0 :
            self.scaling = 100.
        brs = np.copy(br0)
        i = self.counter % brs.size
        br, loc = int(brs[i]), brs[i] % 1
        others = np.concatenate([brs[:i], brs[i+1:]])

        moves = np.random.rand(int(self.scaling + 0.5)) * 2. - 1.
        for m in moves:
            loc += m
            while loc > 0.99999999:
                n = self.tree[br][1]
                if len(n):
                    x, br = n[np.random.choice(len(n))]
                loc -= 0.99999999
                #else:
                #    loc = -1e-8
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
        self.accepted += len(accept1) > 0 or accepted
        self.counter += 1
        if self.counter % 100 == 0:
            if self.accepted / 100. > 0.2:
                self.scaling = np.min([self.scaling * 1.5, 1000.])
            elif self.accepted / 100. <= 0.05:
                self.scaling = np.max([self.scaling / 1.5, 10.])
            self.accepted = 0

        stats = {
            'tune': self.tune,
            'accept': np.exp(accept),
        }
        return br_new, [stats]


if __name__ == '__main__':
    main()

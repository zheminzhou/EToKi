# MGplacer
import sys, os, numpy as np, pandas as pd, subprocess, re, pickle
from numba import jit, prange
from scipy.stats import poisson
from ete3 import Tree
if sys.version_info[0] >= 3:
    xrange = range
try :
    from configure import uopen
except :
    from .configure import uopen

# read ancestral state
def readAncestral(fname, nSample) :
    sites = []
    with uopen(fname) as fin :
        for line in fin :
            if line.startswith('##') :
                pass
            elif line.startswith('#') :
                header = line.strip().split('\t')
                nodes = np.array([not h.startswith('#') for h in header], dtype=bool)
                nodeNames = np.array(header)[nodes]
                mat = pd.read_csv(fin, sep='\t', header=None, dtype=str).values
                break
    mat = mat[np.in1d(mat.T[nodes][0], ['A', 'C', 'G', 'T'])]
    
    conv = np.bincount([ 67,  71, 71,  84, 84, 84])
    data = conv[mat[:, nodes].astype(bytes).view(np.int8)]
    for id, (cont, site) in enumerate(mat[:, :2]) :
        sites.append([(cont, int(site)), ] + [1.]*nSample)
    return np.array(sites), nodeNames, np.ascontiguousarray(data.T, dtype=np.int8)


# read bam file
def parseBAMs(fnames, sites) :
    samtools = 'samtools'
    conv = np.bincount([ 65, 67, 67, 71, 71, 71, 84, 84, 84, 84]) - 1
    knownSamples = np.zeros([sites.shape[0], len(fnames), 4])

    from tempfile import NamedTemporaryFile
    tmpFile = NamedTemporaryFile(delete=False, dir='.')
    siteMap = {}
    for i, s in enumerate(sites) :
        tmpFile.write('{0}\t{1}\n'.format(*s[0]).encode('utf-8'))
        siteMap[s[0]] = i
    tmpFile.close()

    for fnId, fname in enumerate(fnames) :
        p = subprocess.Popen('{0} mpileup -l{2} {1}'.format(samtools, fname, tmpFile.name).split(), universal_newlines=True, stdout=subprocess.PIPE)
        for line in p.stdout :
            cont, site, _, _, bases = line.strip('\n').split('\t')[:5]
            j = siteMap.get((cont, int(site)), -1)
            if j > -1 :
                bases = re.sub(r'\^.|\$', '', bases).upper()
                c = re.findall(r'(\w*)[+-](\d+)(\w+)', bases)
                if len(c) > 0 :
                    bases = ''.join([c[0][0]]+[cc[2][int(cc[1]):] for cc in c])
                bases = conv[np.array(list(bases), dtype=bytes).view(np.int8)]
                bases = bases[bases >= 0]
                baseCounts = np.bincount(bases, minlength=4)
                knownSamples[j, fnId, :] = baseCounts
    os.unlink(tmpFile.name)
    
    sampleDepths = np.sum(knownSamples[:, :, :], 2)
    sampleMean = np.mean(sampleDepths, 0)
    siteWeight = np.min([(1-poisson.cdf(sampleMean/np.sqrt(2), sampleDepths))/(1-poisson.cdf(sampleMean/np.sqrt(2), sampleMean/np.sqrt(2))), \
                         poisson.cdf(sampleMean*np.sqrt(2), sampleDepths)/poisson.cdf(sampleMean*np.sqrt(2), sampleMean*np.sqrt(2))], 0)
    siteWeight[siteWeight > 1] = 1
    for sw, sd, sm in zip(siteWeight.T, sampleDepths.T, sampleMean.T) :
        sw[sd > sm*np.sqrt(2)] /= sd[sd > sm*np.sqrt(2)]/(sm*np.sqrt(2))
    knownSites[:, 1:] = siteWeight
    return knownSamples

def readTree(fname, knownNodes) :
    nodeMap = { n:i for i,n in enumerate(knownNodes) }
    
    branches = []
    tre = Tree(treeFile, format=1)
    cRoot = min([c for c in tre.children if c.name in nodeMap], key=lambda c:nodeMap[c.name])
    for node in tre.traverse('postorder') :
        if node.name in nodeMap and node != cRoot :
            n, d = node, 0.
            while n.up and n.up.name not in nodeMap :
                d += n.dist
                n = n.up
            if n.up :
                branches.append([nodeMap[n.up.name], nodeMap[node.name], d + node.dist, 0])
            else :
                branches.append([nodeMap[cRoot.name], nodeMap[node.name], d + cRoot.dist, 0])
    branches = np.array(branches)
    branches.T[3] = branches.T[2] / np.sum(branches.T[2])
    return branches

def getLikelihood(initState, genoState, weightedSamples, pGenotype, recRate, misRate) :
    nSite, nType, nSample, nGenotype = initState.shape
    genoState[:] = .25*recRate
    genoState[initState == 1] = 1-.75*recRate
    genoState *= (pGenotype.T*(1-misRate)).T
    return getLikelihood2(initState, genoState, weightedSamples, misRate ) # - pGenotype.size - misRate.size - np.sum(np.sum(np.sum(initState, axis=(2,3))>0, 1) * knownSites.T)

@jit(nopython=True)
def greedyPath( positives, allPaths ) :
    (nSite, nPos), nPaths, nGenotypes = positives.shape, range(allPaths.shape[1]), range( allPaths.shape[2] )
    for iSite in range(nSite) :
        pos, path = positives[iSite], allPaths[iSite]
        for i in nPaths :
            p = path[i]
            for j in nGenotypes :
                k = (nPos**j)
                p[j] = pos[int(i /k) % nPos]
    return 

@jit(nopython=True)
def likelihood(allPaths, states, observations, rate, defaultRate) :
    nSite, nType, nSample, nGenotype = states.shape
    nPaths = range(allPaths.shape[1])
    lk = 0
    for iSite in range(nSite) :
        paths = allPaths[iSite]
        state = states[iSite]
        obs = observations[iSite]
        
        siteLK = -99999
        
        for pi in nPaths :
            path = paths[pi]
            rate[:] = defaultRate
            for i in range(nGenotype) :
                p = path[i]
                rate[p] += state[p, :, i]
            rate = np.log(rate)
            pLK = np.sum((rate - np.sum(rate, 0)/4)*obs)
            if pLK > siteLK :
                siteLK = pLK
        lk += siteLK
    return lk

def getLikelihood2(initState, genoState, weightedSamples, misRate) :
    nSite, nType, nSample, nGenotype = initState.shape
    altCount = np.sum(np.sum(weightedSamples, 2) > 0, 1)
    
    lk = 0
    defaultRate = np.zeros([nType, nSample], dtype=float)
    defaultRate[:] = misRate
    
    for altNum in xrange(1, np.max(altCount)+1) :
        altSites = (altCount == altNum)
        observations, states = weightedSamples[altSites], genoState[altSites]
        positives = np.where(np.sum(observations, 2) > 0)[1].reshape([observations.shape[0], altNum])
        allPaths = np.empty([observations.shape[0], altNum**nGenotype, nGenotype], dtype=int)
        greedyPath(positives, allPaths)
        lk += likelihood(allPaths, states, observations, np.zeros([nType, nSample], dtype=float), defaultRate)
    return lk

def updateGenotypes(genotypes, initState, branches, knownMatrix, nChoice=1) :
    initState2 = np.copy(initState)
    genotypes2 = np.copy(genotypes)
    nGenotype = genotypes.shape[0]
    nSite = np.arange(initState.shape[0])
    ids = np.random.choice(nGenotype, replace=False, size=nChoice)
    acc = np.cumsum(branches.T[3])
    brIds = [ np.argmax(r < acc) for r in np.random.rand(nChoice) ]
    genotypes2[ids] = brIds
    for id, brId in zip(ids, brIds) :
        initState2[:, :, :, id] = 0
        initState2[nSite, knownMatrix[int(branches[brId][0])], :, id] = 1
        initState2[nSite, knownMatrix[int(branches[brId][1])], :, id] = 1
    return genotypes2, initState2
    

if __name__ == '__main__' :
    ancestralFile, bamFile, treeFile = sys.argv[1:]
    
    bamFiles = bamFile.split(',')
    # number of possible genotypes
    nGenotype = 4
    nType = 4
    nSample = len(bamFiles)
    
    if not os.path.isfile('test1') :
        knownSites, knownNodes, knownMatrix = readAncestral(ancestralFile, nSample)
        knownSamples = parseBAMs(bamFiles, knownSites).transpose([0,2,1])
        pickle.dump([knownSites, knownNodes, knownSamples, knownMatrix], open('test', 'wb'))
    knownSites, knownNodes, knownSamples, knownMatrix = pickle.load(open('test', 'rb'))

    nSite = knownSites.shape[0]

    # read tree
    branches = readTree(treeFile, knownNodes)
    nBranch = branches.shape[0]
    
    # global parameters
    globalParams = dict(
        recRate = 0.1, 
        misRate = np.array([0.01] * nSample), 
        pGenotype = np.random.rand(nSample*nGenotype).reshape([nSample, nGenotype]), 
        genotypes = np.zeros(nGenotype, dtype=int), 
    )
    globalParams['pGenotype'] /= np.sum(globalParams['pGenotype'], 1).reshape(globalParams['pGenotype'].shape[0], 1)
    
    # states
    initState = np.zeros([nSite, nType, nSample, nGenotype], dtype=int) 
    genoState = np.zeros([nSite, nType, nSample, nGenotype])
    
    weightedSamples = (knownSamples.T * knownSites[:, 1:].astype(float).T).T
    
    globalParams['genotypes'], initState = updateGenotypes(globalParams['genotypes'], initState, branches, knownMatrix, nChoice=nGenotype)
    import time
    
    lk = getLikelihood(initState, genoState, weightedSamples, globalParams['pGenotype'], globalParams['recRate'], globalParams['misRate'])
    
    print(-1, lk, globalParams)
    nChain = 100000
    for ite in xrange(nChain) :
        # update initState
        for i in range(3) :
            newGenotypes, newInitState = updateGenotypes(globalParams['genotypes'], initState, branches, knownMatrix, nChoice=1)
            if np.sum(newGenotypes) != np.sum(globalParams['genotypes']) :
                newLK = getLikelihood(newInitState, genoState, weightedSamples, globalParams['pGenotype'], globalParams['recRate'], globalParams['misRate'])
                if newLK >= lk or np.random.ranf() < np.exp(newLK-lk) :
                    globalParams['genotypes'][:] = newGenotypes
                    initState = newInitState
        # update parameters
        newRecRate = min(max((np.random.ranf()-0.5)*2*0.02 + globalParams['recRate'], 0.01), 0.99)
        if newRecRate != globalParams['recRate'] :  
            newLK = getLikelihood(initState, genoState, weightedSamples, globalParams['pGenotype'], newRecRate, globalParams['misRate'])
            if newLK >= lk or np.random.ranf() < np.exp(newLK-lk) :
                globalParams['recRate'] = newRecRate
                lk = newLK
        newMisRate = np.copy(globalParams['misRate'])
        mrId = np.random.choice(newMisRate.shape[0]) if newMisRate.shape[0] > 0 else 0
        rr = min(max((np.random.ranf()-0.5)*2*0.02 + newMisRate[mrId], 0.01), 0.99)
        if rr != newMisRate[mrId] :
            newMisRate[mrId] = rr
            newLK = getLikelihood(initState, genoState, weightedSamples, globalParams['pGenotype'], globalParams['recRate'], newMisRate)
            if newLK >= lk or np.random.ranf() < np.exp(newLK-lk) :
                globalParams['misRate'][:] = newMisRate
                lk = newLK
        if globalParams['pGenotype'].shape[1] > 1 :
            newPGenotype = np.copy(globalParams['pGenotype'])
            pgId = np.random.choice(newPGenotype.shape[0]) if newPGenotype.shape[0] > 0 else 0
            sId, tId = np.random.choice(newPGenotype.shape[1], 2, replace=False)
            pgTr = np.min([newPGenotype[pgId][sId]-0.02, 0.98-newPGenotype[pgId][tId], np.random.ranf()*0.02])
            if pgTr > 0 :
                newPGenotype[pgId][sId] -= pgTr
                newPGenotype[pgId][tId] += pgTr
                newLK = getLikelihood(initState, genoState, weightedSamples, newPGenotype, globalParams['recRate'], globalParams['misRate'])
                if newLK >= lk or np.random.ranf() < np.exp(newLK-lk) :
                    globalParams['pGenotype'][:] = newPGenotype
                    lk = newLK
        print(ite, lk, globalParams, [ knownNodes[int(branches[gt][1])] for gt in globalParams['genotypes'] ])
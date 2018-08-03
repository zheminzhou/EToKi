# MGplacer
import sys, os, numpy as np, pandas as pd, subprocess, re, pickle
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


def getLikelihood(initState, genoState, knownSites, knownSamples, pGenotype, recRate, misRate) :
    nSite, nGenotype, _ = initState.shape
    genoState[:] = np.where(initState == 1, 1-0.75*recRate, 0.25*recRate)
        
    for weights, observations, state in zip(knownSites[:, 1:], knownSamples, genoState) :
        primaryState = (pGenotype*(1-misRate.reshape(misRate.size, 1))).reshape(list(pGenotype.shape) + [1])*state
        secondaryState = (0.25*misRate).reshape(misRate.size, 1)*pGenotype
        
        alters = np.where(np.sum(observations, 0) > 0)[0]
        t1 = observations[:, alters]
        t2 = np.sum(observations, 1).reshape(2, 1) - t1
        alterCounts = np.array(list(zip(t1, t2))).transpose(0, 2, 1)
        
        alterChoice = np.zeros([alters.size**nGenotype, nGenotype], dtype=np.int8)
        for id in np.arange(alterChoice.shape[1]) :
            alterChoice.T[id] = np.where(alters >= 0)[0][((np.arange(alterChoice.shape[0]) / (alters.size**id)).astype(int) % alters.size)]
        
        alterCounts = 1
        tmp = genoState*(1 - 0.25*m)*p.reshape([3,1,1,1])
        genoState[:] = np.where(initState==1, 1-3/4*recRate, recRate/4)
    
    for s, observeration, weights in zip(genoState, knownSamples, knownSites[:, 1:]) :
        for obs, w, p, m in zip(observeration, weights, pGenotype, misRate) :
            
            print(s, o)
    return

def updateGenotypes(genotypes, initState, branches, knownMatrix, nChoice=1) :
    nGenotype = genotypes.shape[0]
    nSite = np.arange(initState.shape[0])
    ids = np.random.choice(nGenotype, replace=False, size=nChoice)
    acc = np.cumsum(branches.T[3])
    brIds = [ np.argmax(r < acc) for r in np.random.rand(nChoice) ]
    genotypes[ids] = brIds
    for id, brId in zip(ids, brIds) :
        initState[nSite, id, knownMatrix[int(branches[brId][0])]] = 1
        initState[nSite, id, knownMatrix[int(branches[brId][1])]] = 1
    return
    

if __name__ == '__main__' :
    ancestralFile, bamFile, treeFile = sys.argv[1:]
    
    bamFiles = bamFile.split(',')
    # number of possible genotypes
    nGenotype = 3
    nType = 4
    nSample = len(bamFiles)
    
    if not os.path.isfile('test') :
        knownSites, knownNodes, knownMatrix = readAncestral(ancestralFile, nSample)
        knownSamples = parseBAMs(bamFiles, knownSites)
        pickle.dump([knownSites, knownNodes, knownSamples, knownMatrix], open('test', 'wb'))
    knownSites, knownNodes, knownSamples, knownMatrix = pickle.load(open('test', 'rb'))

    nSite = knownSites.shape[0]

    # read tree
    branches = readTree(treeFile, knownNodes)
    nBranch = branches.shape[0]
    
    # global parameters
    globalParams = dict(
        recRate = 0.1, 
        misRate = np.random.rand(nSample) * 0.02, 
        pGenotype = np.random.rand(nSample*nGenotype).reshape([nSample, nGenotype]), 
        genotypes = np.zeros(nGenotype, dtype=int), 
    )
    globalParams['pGenotype'] /= np.sum(globalParams['pGenotype'], 1).reshape(globalParams['pGenotype'].shape[0], 1)
    
    # states
    initState = np.zeros([nSite, nGenotype, nType], dtype=np.int8) 
    genoState = np.zeros([nSite, nGenotype, nType])
    emitState = np.zeros([nSite, nType]) # emission likelihood combination of all genotypes
    
    updateGenotypes(globalParams['genotypes'], initState, branches, knownMatrix, nChoice=nGenotype)
    getLikelihood(initState, genoState, knownSites, knownSamples, globalParams['pGenotype'], globalParams['recRate'], globalParams['misRate'])
    
    nChain = 100000
    for ite in xrange(nChain) :
        # update initState
        
        # update parameters
        updateParams()


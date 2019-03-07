import os, sys, re, subprocess, numpy as np, pandas as pd
from sklearn.mixture import GaussianMixture
try :
    from  configure import uopen, externals, readFasta, xrange
except :
    from .configure import uopen, externals, readFasta, xrange

# load current matrix
# infer core genomic region
def loadMatrix(fname) :
    sequences, snps = {}, {}
    with uopen(fname) as fin :
        for line in fin :
            if line.startswith('##') :
                part = line.strip().split()
                if line.startswith('## Sequence_length') :
                    sequences[part[2]] = np.zeros(int(part[3]), dtype=np.int8)
                elif line.startswith('## Missing_region') :
                    sequences[part[2]][int(part[3])-1:int(part[4])] = -1
            else :
                headers = line.strip().split('\t')
                break
        matrix = pd.read_csv(fin, header=None, sep='\t').values
        encode = {'A':1, 'C':2, 'G':4, 'T':8}
        matrix.T[4] = list(map(lambda d:encode.get(d[0], -9999) + encode.get(d[-1], -9999), matrix.T[4].tolist()))
        matrix = matrix[matrix.T[4] > 0]
        matrix.T[2] -= 1
        for m in matrix :
            sequences[m[1]][m[2]] |= m[4]
            if (m[1], m[2]) not in snps :
                snps[(m[1], m[2])] = {}
            if m[4] not in snps[(m[1], m[2])] :
                snps[(m[1], m[2])][m[4]] = []
            snps[(m[1], m[2])][m[4]].append(m[0])
    return sequences, snps

# load bam file
def loadBam(prefix, reference, bams, sequences, snps) :
    sites = []
    p = subprocess.Popen('samtools mpileup -ABQ0 {0}'.format(' '.join(bams)).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    for line in p.stdout :
        part = line.strip().split('\t')
        s = int(part[1])-1
        if s % 100000 == 0 :
            sys.stdout.write('# {0}\n'.format(s))        
        if sequences[part[0]][s] > 0 or s % 5 == 0:
            bases = ''.join(part[4::3])
            bases = re.sub('[\*\$\+-]', '', re.sub(r'\^.', '', bases.upper()))
            bases = re.split('(\d+)', bases)
            for i in range(1, len(bases), 2) :
                bases[i+1] = bases[i+1][int(bases[i]):]
            types, cnts = np.unique(list(''.join(bases[::2])), return_counts=True)
            if np.sum(cnts) >= 3 :
                if types.size > 1 :
                    cnts.sort()
                    sites.append([cnts[-1], np.sum(cnts[:-1])])
                else :
                    sites.append([cnts[0], 0])
    sites = np.array(sites)
    ave_depth = np.max([np.median(np.sum(sites, 1)), 2.])
    sys.stdout.write('{3}: Average read depth: {0}; Sites between {1} and {2} will be used for hybrid estimation.\n'.format(ave_depth, ave_depth/2., ave_depth*3., prefix))
    sites = sites[(ave_depth/2. <= np.sum(sites, 1)) & (np.sum(sites, 1) <= ave_depth*3)]
    
    m=GaussianMixture(n_components=1, covariance_type='tied')
    m.fit(sites)
    best_model = [m.bic(sites), m]
    for n_components in xrange(2, 6) :
        sys.stdout.write('# Testing {0} components.\n'.format(n_components))
        m=GaussianMixture(n_components=n_components, covariance_type='tied')
        for i in xrange(20) :
            m.fit(sites)
            bic = m.bic(sites)
            if bic < best_model[0] :
                best_model = [bic, m]
                m=GaussianMixture(n_components=n_components, covariance_type='tied')
    m = best_model[1]
    mId = np.argmax(m.means_.T[1] / np.sum(m.means_, 1))
    sys.stdout.write('{3}: Find {0} GMM components. The most divergent group is {1} and counts for {2} of total sites.\n'.format(m.n_components, m.means_[mId].tolist(), m.weights_[mId], prefix))
    mDiv = m.means_[mId][0]/np.sum(m.means_[mId])
    mDiv = 10*np.log10([[mDiv, 1-mDiv], [1-mDiv, mDiv]])
    
    seq = {n:list(s) for n,s in readFasta(reference).items()}
    qual = {n:[0] * len(s) for n, s in seq.items()}
    
    lowQ, lowC, highQ = 0, 0, 0
    p = subprocess.Popen('samtools mpileup -ABQ0 {0}'.format(' '.join(bams)).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    for line in p.stdout :
        part = line.strip().split('\t')
        s = int(part[1])-1
        if s % 100000 == 0 :
            sys.stdout.write('# {0}\n'.format(s))        
        bases = ''.join(part[4::3])
        bases = re.sub('[\*\$\+-]', '', re.sub(r'\^.', '', bases.upper()))
        bases = re.split('(\d+)', bases)
        for i in range(1, len(bases), 2) :
            bases[i+1] = bases[i+1][int(bases[i]):]
        types, cnts = np.unique(list(''.join(bases[::2])), return_counts=True)
        if types.size > 0 :
            depth = np.sum(cnts)
            if cnts.size == 1 :
                g, mId = [cnts[0], 0], 0
            elif cnts.size > 1 :
                mId = np.argmax(cnts)
                g = [cnts[mId], depth-cnts[mId]]
            seq[part[0]][s] = types[mId]
            if depth >= 3 and depth/3. <= ave_depth <= depth*3. :
                q = min(40, max(1, int(round(np.sum(g * mDiv[0]) - np.sum(g * mDiv[1]), 0))))
                qual[part[0]][s] = q
                if q < 10 :
                    lowQ += 1
                else :
                    highQ += 1
            else :
                lowC += 1
    qual = {n:''.join([ chr(ss+33) for ss in s ]) for n, s in qual.items() }
    with open(prefix+'.fastq', 'w') as fout :
        for n, s in seq.items() :
            fout.write('@{0}\n{1}\n+\n{2}\n'.format(n, ''.join(s), qual[n]))
    sys.stdout.write('{0}: {1} good sites; {2} low covered sites; {3} low quality sites;\n'.format(prefix, highQ, lowC, lowQ))        
    
    return

def MGplacer0(args) :
    sequences, snps = loadMatrix(args[1])
    loadBam(args[0], args[2], args[3:], sequences, snps)

if __name__ == '__main__' :
    MGplacer0(sys.argv[1:])
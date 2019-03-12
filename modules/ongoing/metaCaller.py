import os, sys, re, subprocess, numpy as np, pandas as pd
from scipy.stats import binom
try :
    from  configure import uopen, externals, readFasta, xrange, logger
except :
    from .configure import uopen, externals, readFasta, xrange, logger

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
    sequence = readFasta(reference)
    sequence = {n:[s, [0]*len(s)] for n, s in sequence.items()}

    sites = {}
    for bam in bams :
        if bam is not None :
            depth = subprocess.Popen('{samtools} depth -q 0 -Q 0 {bam}'.format(bam=bam, **externals).split(), stdout=subprocess.PIPE, universal_newlines=True)
            try:
                d = pd.read_csv(depth.stdout, sep='\t').values
                sites.update({ cName:1 for cName in np.unique(d.T[0]) })
            except :
                pass

    sequence = {n:s for n, s in sequence.items() if n in sites}
    with open('{0}.mapping.reference.fasta'.format(prefix), 'w') as fout :
        for n, s in sorted(sequence.items()) :
            fout.write('>{0}\n{1}\n'.format(n, '\n'.join([ s[0][site:(site+100)] for site in xrange(0, len(s[0]), 100)])))

    bam_opt = ' '.join(['--bam {0}'.format(b) for b in bams if b is not None])
    pilon_cmd = '{pilon} --fix snps,indels,gaps --vcf --output {prefix}.mapping --genome {prefix}.mapping.reference.fasta {bam_opt}'.format(prefix=prefix, bam_opt=bam_opt, **externals)
    subprocess.Popen( pilon_cmd.split(), stdout=subprocess.PIPE, universal_newlines=True ).communicate()

    uncertains = []
    with open('{0}.mapping.vcf'.format(prefix)) as fin :
        for line in fin :
            if line.startswith('#') : continue
            part = line.strip().split('\t')
            if sequences[part[0]][int(part[1])-1] >= 0 :
                if len(part[3]) == 1 and len(part[4]) == 1 :
                    pp = part[7].split(';')
                    dp = float(pp[0][3:])
                    if dp >= 3 :
                        qd = int(pp[4][3:])
                        if part[-1] == '0/1' or qd < 10 :
                            bcs = sorted([float(bc) for bc in pp[5][3:].split(',')])
                            uncertains.append([bcs[-1], np.sum(bcs[:-1])])
    uncertains = np.array(uncertains)
    p = np.sum(uncertains.T[0])/np.sum(uncertains)
    qPerRead = 10*(np.log10(p) - np.log10(1-p))
    for n in sequence :
        sequence[n][0] = list(sequence[n][0])
        
    highQ, lowQ, lowC = 0, 0, 0
    with open('{0}.mapping.vcf'.format(prefix)) as fin :
        for line in fin :
            if line.startswith('#') : continue
            part = line.strip().split('\t')
            if len(part[3]) == 1 and len(part[4]) == 1 :
                s = int(part[1]) - 1
                pp = part[7].split(';')
                dp = float(pp[0][3:])
                qd = int(pp[4][3:])
                if part[-1] == '0/1' or qd < 10 :
                    bcs = np.array([int(bc) for bc in pp[5][3:].split(',')])
                    if np.sum(bcs) > 0 :
                        sequence[part[0]][0][s] = ['A', 'C', 'G', 'T'][np.argmax(bcs)]
                    else :
                        sequence[part[0]][0][s] = part[3]
                    if dp < 3 :
                        lowC += 1
                    else :
                        bcs.sort()
                        bcs = [bcs[-1], np.sum(bcs[:-1])]
                        q1 = binom.cdf(bcs[0], bcs[0]+bcs[1], p)
                        q2 = qPerRead * (bcs[0] - bcs[1]) if q1 >= 0.05 else 1
                        if q2 >= 10 :
                            highQ += 1
                        else :
                            lowQ += 1
                        sequence[part[0]][1][s] = min(40, max(1, int(q2)))
                else :
                    if dp < 3 :
                        lowC += 1
                    else :
                        if qd >= 10 :
                            highQ += 1
                        else:
                            lowQ += 1
                        sequence[part[0]][1][s] = qd
                    if part[-1] == '1/1' :
                        sequence[part[0]][0][s] = part[4]

    logger('{0}: Expected mix-up: {1} {2} ; Got highQ {3} ; lowQ {4} ; lowC {5}'.format(prefix, uncertains.shape[0], p, highQ, lowQ, lowC))
    with open('{0}.metaCaller.fastq'.format(prefix), 'w') as fout :
        p = prefix.rsplit('/', 1)[-1]
        for n, (s, q) in sequence.items() :
            fout.write( '@{0}\n{1}\n+\n{2}\n'.format( p+'_'+n, ''.join(s), ''.join([chr(qq+33) for qq in q]) ) )
    os.unlink( '{0}.mapping.vcf'.format(prefix) )
    os.unlink( '{0}.mapping.fasta'.format(prefix) )
    os.unlink( '{0}.mapping.reference.fasta'.format(prefix) )
    return '{0}.metaCaller.fastq'.format(prefix)

def metaCaller(args) :
    (prefix, mutations, reference), bams = args[:3], args[3:]
    sequences, snps = loadMatrix(mutations)
    loadBam(prefix, reference, bams, sequences, snps)

if __name__ == '__main__' :
    metaCaller(sys.argv[1:])
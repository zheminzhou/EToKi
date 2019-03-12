import os, sys, numpy as np
try:
    from .configure import uopen
except :
    from configure import uopen

conservedSites = np.zeros(255, dtype=int)
def parseSNP(fout, sequences, core) :
    contigName = sequences[0][1]
    seqNames = np.array([seqId for seqId, name, seq in sequences])
    bases = np.array([list(seq) for seqId, name, seq in sequences ], dtype=bytes, order='F').T#.view(dtype=np.uint8)
    bases[~np.in1d(bases, [b'A', b'C', b'G', b'T']).reshape(bases.shape)] = b'-'
    nSeq = bases.shape[1]
    missings = [[contigName, -1, -1]]
    for site, base in enumerate(bases) :
        types, counts = np.unique(base, return_counts=True)
        if counts[types == b'-'] > core * nSeq :
            if site > missings[-1][1] :
                missings.append([contigName, site+1, site+1])
            else :
                missings[-1][2] = site + 1
        else :
            types = types[types != b'-'].view(np.uint8)
            if len(types) == 1 :
                conservedSites[types[0]] += 1
            else :
                fout.write('{0}\t{1}\t{2}\n'.format( contigName, site+1, b'\t'.join(base).decode('utf-8') ))
    return [contigName, bases.shape[0]], missings[1:]

def readXFasta(fasta_file, core=0.8) :
    seqs = []
    nameMap = {}
    contLens, missingLens = [], []
    with uopen(fasta_file) as fin, open(fasta_file+'.tmp', 'w') as fout :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                seqName = name.split(':', 1)[0]
                contName = name.split(':', 1)[-1]
                if seqName in nameMap :
                    seqId = nameMap[seqName]
                    seqs[seqId] = [seqName, contName, []]
                else :
                    seqId = nameMap[seqName] = len(seqs)
                    seqs.append([seqName, contName, []])
            elif line.startswith('=') :
                for seq in seqs :
                    seq[2] = ''.join(seq[2]).upper()
                contName = [seq[1] for seq in seqs if seq[2] != ''][0]
                contLen = [len(seq[2]) for seq in seqs if seq[2] != ''][0]
                for seq in seqs :
                    if seq[2] == '' :
                        seq[2] = '-' * contLen
                cLen, mLen = parseSNP(fout, seqs, core)
                contLens.append(cLen)
                missingLens.extend(mLen)
                seqs = [ [n, '', []] for n in nameMap ]
            else :
                seqs[seqId][2].extend(line.strip().split())
        if len(seqs[0][2]) > 0 :
            cLen, mLen = parseSNP(fout, seqs, core)
            contLens.append(cLen)
            missingLens.extend(mLen)
    sys.stdout.write('## Constant_bases: {0} {1} {2} {3}\n'.format(conservedSites[ord('A')], conservedSites[ord('C')], conservedSites[ord('G')], conservedSites[ord('T')]))
    for cLen in contLens :
        sys.stdout.write('## Sequence_length: {0} {1}\n'.format(*cLen))
    for mLen in missingLens :
        sys.stdout.write('## Missing_region: {0} {1} {2}\n'.format(*mLen))
    sys.stdout.write('#Seq\t#Site\t{0}\n'.format('\t'.join([n for n, i in sorted(nameMap.items(), key=lambda n:n[1])])))
    with uopen(fasta_file+'.tmp') as fin :
        for line in fin :
            sys.stdout.write(line)
    os.unlink(fasta_file+'.tmp')
    return 


if __name__ == '__main__' :
    readXFasta(sys.argv[1])
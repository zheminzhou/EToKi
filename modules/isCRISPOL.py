import sys, os, subprocess, re, pandas as pd, numpy as np
try :
    from configure import externals
except :
    from .configure import externals
    
crispolDB = os.path.join(os.path.dirname(__file__), 'CRISPOL.db')

def refSync(regions, maxGap=600, maxGapDiff = 600) :
    regions.sort(key=lambda r:(r[0], r[6], r[1], r[8]))
    syncRegions = []
    lenRegions = len(regions)
    xId = 0
    for i1, r in enumerate(regions) :
        if (r[7] - r[6] + 1) >= r[12] * 0.6 :
            syncRegions.append(r + [xId, r[11]/r[12]] + ([r[8], r[9]] if r[8] > 0 else [-r[9], -r[8]]))
            xId += 1
        for i2 in range(i1+1, lenRegions) :
            r2 = regions[i2]
            g1 = r2[6] - r[7] - 1
            if r[0] != r2[0] or g1 > maxGap :
                break
            e1 = min(r2[6]-r[6], r2[7]-r[7])
            if (r2[7] - r[6]+1) < r[12]*0.6 or e1 <= 30 :
                continue
            c1 = 1. if g1 <= 0 else (r[7]-r[6]+1 + r2[7]-r2[6]+1)/(r2[7] - r[6]+1)
            if c1 <= 0.6 :
                continue
            if r[1] == r2[1] and (r[8] > 0) == (r2[8] > 0) :
                g2 = r2[8] - r[9] - 1
                e2 = min(r2[8]-r[8], r2[9]-r[9])
                c2 = 1. if g2 <= 0 else (r[9]-r[8]+1 + r2[9]-r2[8]+1)/(r2[9] - r[8]+1)
                if g2 <= maxGap and e2 >= 30 and c2 >= 0.6 and abs(g1 - g2) <= maxGapDiff :
                    g = min(g2, g1)
                    if g >= 0 :
                        score = r[11] + r2[11]
                    else :
                        score = r[11] + r2[11] - (-g)*min(r[11]/(r[7]-r[6]+1), r2[11]/(r2[7]-r2[6]+1))
                    syncRegions.extend([ r + [xId, score/r[12]] + ([r[8], r2[9]] if r[8] > 0 else [-r2[9], -r[8]]), \
                                         r2 + [xId, score/r[12]] + ([r[8], r2[9]] if r[8] > 0 else [-r2[9], -r[8]]) ])
                    xId += 1
                    continue
            g2 = ( r[13]-r[9] - 1 if r[8] > 0 else -1 - r[9] ) + ( r2[8]-1 if r2[8] > 0 else r[13]+r2[8]-1 )
            if g2 < 100 :
                if g1 >= 0 :
                    score = r[11] + r2[11]
                else :
                    score = r[11] + r2[11] - (-g1)*min(r[11]/(r[7]-r[6]+1), r2[11]/(r2[7]-r2[6]+1))
                syncRegions.extend([ r + [xId, score/r[12]] + ([r[8], r[9]] if r[8] > 0 else [-r[9], -r[8]]), \
                                     r2 + [xId, score/r[12]] + ([r2[8], r2[9]] if r2[8] > 0 else [-r2[9], -r2[8]]) ])
                xId += 1
    return syncRegions

def qrySync(regions) :
    deleted = {}
    regions.sort(key=lambda r:(r[1], r[16]))
    lenRegions = len(regions)
    for i1, r1 in enumerate(regions) :
        if r1[14] in deleted : continue
        toDel = {}
        for i2 in range(i1+1, lenRegions) :
            r2 = regions[i2]
            if r1[14] == r2[14] or r2[14] in deleted :
                continue
            if r1[1] != r2[1] or r2[16] > r1[17] :
                break
            d = r1[17] - r2[16] + 1
            if d >= (r1[17] - r1[16] + 1) * 0.5 or d >= (r2[17] - r2[16] + 1) * 0.5 :
                if r1[15] >= r2[15] :
                    toDel[r2[14]] = 1
                else :
                    deleted[r1[14]] = 1
                    break
        if r1[14] not in deleted :
            deleted.update(toDel)
    return sorted([r for r in regions if r[14] not in deleted], key=lambda r:(r[14], r[6]))

def getCRISPOL(refs, regions) :
    spacers = []
    with open(refs) as fin :
        for line in fin :
            if line.startswith('>') :
                name = line.strip().split()
                name[0] = name[0][1:]
                name[1] = int(name[1])
                if name[1] > 0 :
                    spacers.append(name)
    spacers = { s[0]:[0, int(s[1])] for s in spacers }
    spacers['CORE'][0] = 2
    
    regions.sort(key=lambda x:(x[1], x[8]))
    lenRegion = len(regions)
    for rId, region in enumerate(regions) :
        if region[0] not in spacers :
            continue
        if rId > 0 :
            r2 = regions[rId-1]
            if r2[1] == region[1] and (r2[8] > 0) == (region[8] > 0) and abs(r2[9] - region[8]) < 10 :
                spacers[region[0]][0] = 1
        if spacers[region[0]][0] == 0 and rId < lenRegion - 1 :
            r2 = regions[rId+1]
            if r2[1] == region[1] and (r2[8] > 0) == (region[8] > 0) and abs(r2[8] - region[9]) < 10 :
                spacers[region[0]][0] = 1
        if spacers[region[0]][0] and region[0].find('var') > 0 :
            v = region[0].split('var')[0]
            spacers[v][0] = 1
    spacers = sorted([ [so, sn, st] for sn, (st, so) in spacers.items() ])
        
    return spacers
    

def blast2region(qry, method='blastn', minIdentity=92, minCover=19) :
    if method in ('blastn', 'tblastn') :
        subprocess.Popen('{makeblastdb} -dbtype nucl -in {0}'.format(qry, **externals).split(), stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    else :
        subprocess.Popen('{makeblastdb} -dbtype prot -in {0}'.format(qry, **externals).split(), stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    blast = subprocess.Popen([externals[method], '-task', 'blastn', '-db', qry, '-query', crispolDB, '-outfmt', "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen", \
                              '-evalue', '0.1'], stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines=True)
    
    outputs = pd.read_csv(blast.stdout, sep='\t', header=None).values
    outputs = outputs[(outputs.T[2]>=minIdentity) & (outputs.T[7]-outputs.T[6]+1 >= minCover)]
    outputs[outputs.T[8]>outputs.T[9], 8:10] = -outputs[outputs.T[8]>outputs.T[9], 8:10]
    if method == 'tblastn' :
        outputs.T[6] = o(utputs.T[6]-1)*3+1
        outputs.T[7] = outputs.T[7]*3
        outputs.T[12] *= 3
    elif method == 'blastx' :
        outputs[outputs.T[8] > 0, 8] = (outputs.T[8]-1)*3+1
        outputs[outputs.T[8] > 0, 9] = outputs.T[9]*3
        outputs[outputs.T[8] <= 0, 8] = outputs.T[8]*3
        outputs[outputs.T[8] <= 0, 9] = (outputs.T[9]+1)*3-1
        outputs.T[13] *= 3
    os.unlink(qry+'.nhr')
    os.unlink(qry+'.nin')
    os.unlink(qry+'.nsq')
    outputs = refSync(outputs.tolist())
    outputs = qrySync(outputs)
    return outputs

def isCRISPOL(args) :
    import argparse
    parser = argparse.ArgumentParser(description='''
EToKi.py isCRISPOL
''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('assemblies', metavar='N', help='FASTA files containing assemblies of S. enterica Typhimurium.', nargs='*')
    params = parser.parse_args(args)
    
    for qry in params.assemblies :
        region = blast2region(qry)
        crispol = getCRISPOL(crispolDB, region)
        print('{0}\t{1}\t{2}'.format(qry, ','.join([c[1] for c in crispol]), ','.join([str(c[2]) for c in crispol]) ))

if __name__ == '__main__' :
    isCRISPOL(sys.argv[1:])

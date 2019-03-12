import os, sys, tempfile, time, shutil, numpy as np, pandas as pd, re
from numba import jit
from subprocess import Popen, PIPE
from multiprocessing.pool import ThreadPool
from operator import itemgetter
try:
    from .configure import externals, logger, xrange, readFastq, transeq, blosum62, rc, asc2int
except :
    from configure import externals, logger, xrange, readFastq, transeq, blosum62, rc, asc2int

makeblastdb = externals['makeblastdb']
blastn = externals['blastn']
usearch = externals['usearch']
mmseqs = externals['mmseqs']
minimap2 = externals['minimap2']

@jit(nopython=True)
def tab2overlaps(tabs, ovl_l, ovl_p, nTab, overlaps) :
    ovlId = 0
    for i1 in xrange(overlaps[-1, 0], nTab) :
        t1 = tabs[i1]
        ovl_l2 = min(ovl_l, ovl_p*(t1[3]-t1[2]+1))
        if i1 > overlaps[-1, 0] :
            i2r = xrange(i1+1, nTab)
        else :
            i2r = xrange(overlaps[-1, 1], nTab)
        for i2 in i2r :
            t2 = tabs[i2]
            if t1[0] != t2[0] or t2[2] > t1[3] : break
            ovl = min(t1[3], t2[3]) - t2[2] + 1
            if ovl >= ovl_l2 or ovl >= ovl_p*(t2[3]-t2[2]+1) :
                overlaps[ovlId, :] = [t1[1], t2[1], ovl]
                ovlId += 1
                if ovlId == 1000000 :
                    overlaps[-1, :2] = [i1, i2]
                    break
        if ovlId == 1000000 :
            break
    if ovlId < 1000000 :
        overlaps[-1, :] = -1
    return overlaps


def _linearMerge(data) :
    matches, params = data
    grpCol = pd.Series(data= [[]] * matches.shape[0])
    matches = np.hstack([matches, grpCol.values[:, np.newaxis]])
    gapDist, lenDiff = params[1:]
    gene, geneLen = matches[0][0], matches[0][12]
    tailing = 20
    
    def resolve_edges(edges) :
        grps = []
        for id, m1 in edges[0] :
            for jd, m2 in edges[1] :
                if (m1[1] == m2[1] and max(abs(m1[8]), abs(m1[9])) > min(abs(m2[8]), abs(m2[9])) ) or \
                   abs(m1[2]-m2[2]) > 0.3 or m1[6] >= m2[6] or m1[7] >= m2[7] or m2[6]-m1[7]-1 >=gapDist:
                    continue
                rLen = m2[7] - m1[6] + 1
                g1 = -m1[9]-1 if m1[9] < 0 else m1[13] - m1[9]
                g2 =  m1[8]-1 if m1[8] > 0 else m1[13] + m1[8]
                qLen = m1[9]-m1[8]+1 + m2[9]-m2[8]+1 + g1 + g2
                if g1+g2 >= gapDist or min(rLen, qLen)*lenDiff < max(rLen, qLen) :
                    continue
                overlap = sorted([m1[7] - m2[6] + 1, -g1-g2], reverse=True)

                rLen1, rLen2 = m1[7] - m1[6] + 1, m2[7] - m2[6] + 1
                if overlap[0] > 0 :
                    score = m1[11] + m2[11] - overlap[0] * min( float(m1[11])/rLen1, float(m2[11])/rLen2 )
                    ident = (m1[2]*rLen1 + m2[2]*rLen2 - overlap[0] * min(m1[2], m2[2]))/(rLen1 + rLen2 - overlap[0])
                else :
                    score = m1[11] + m2[11]
                    ident = (m1[2]*rLen1 + m2[2]*rLen2)/(rLen1 + rLen2)
                if overlap[1] < 0 :
                    score +=  overlap[1]/3.
                if score > m1[11] and score > m2[11] :
                    grps.append( [ score, ident, rLen, 1, id, jd ] )
        return grps
    
    groups = []
    prev, edges = matches[0][1], [[], []]
    nSave = len(matches)
    
    for id, m1 in enumerate(matches) :
        rLen1 = m1[7] - m1[6] + 1
        groups.append([ m1[11], m1[2], rLen1, 0, id ])
        if m1[6] > tailing and ((m1[8] > 0 and m1[8] - 1 <= gapDist) or (m1[8] < 0 and m1[13] + m1[8] < gapDist)) :   # any hit within the last 150 bps to either end of a scaffold is a potential fragmented gene
            edges[1].append([id, m1])
        if m1[7] <= m1[12] - tailing :
            if (m1[8] > 0 and m1[13]-m1[9] <= gapDist) or (m1[8] < 0 and -1-m1[9] < gapDist) :
                edges[0].append([id, m1])
            for jd in xrange(id+1, nSave) :
                m2 = matches[jd]
                if m1[1] != m2[2] or (m1[8] < 0 and m2[8] > 0) or m2[8] - m1[9] -1 >= gapDist :    # maximum 300bps between two continuous hits in the same scaffold
                    break
                rLen, qLen = m2[7]-m1[6]+1, m2[9]-m1[8]+1
                if abs(m1[2]-m2[2]) > 0.3 or m1[9] >= m2[9] or m1[6] >= m2[6] or m1[7] >= m2[7] or m2[6] - m1[7] -1 >= gapDist \
                   or min(rLen, qLen)*lenDiff < max(rLen, qLen) :
                    continue
                rLen2 = m2[7] - m2[6] + 1
                overlap = sorted([m1[7]-m2[6]+1, m1[9]-m2[8]+1], reverse=True)
                if overlap[0] > 0 :
                    score = m1[11] + m2[11] - overlap[0] * min( float(m1[11])/rLen1, float(m2[11])/rLen2 )
                    ident = (m1[2]*rLen1 + m2[2]*rLen2 - overlap[0]*min(m1[2], m2[2]))/(rLen1 + rLen2 - overlap[0])
                else :
                    score = m1[11] + m2[11]
                    ident = (m1[2]*rLen1 + m2[2]*rLen2)/(rLen1 + rLen2)
                if overlap[1] < 0 :
                    score +=  overlap[1]/3.
                if score > m1[11] and score > m2[11] :
                    groups.append( [ score, ident, rLen, 0, id, jd ] )
    if len(edges[0]) and len(edges[1]) :
        groups.extend(resolve_edges(edges))
    if len(groups) > len(matches) :
        groups.sort(reverse=True)
        usedMatches, usedGroups = {}, []
        for grp in groups :
            if (grp[4], 4) in usedMatches or (grp[-1], 5) in usedMatches :
                continue
            if grp[3] > 0 :
                if (grp[4], 5) in usedMatches or (grp[-1], 4) in usedMatches :
                    continue
            if grp[4] != grp[-1] :
                lMat, rMat = matches[grp[4]], matches[grp[-1]]
                il, im = sorted([grp[4], grp[-1]])
                skp = 0
                for i in xrange(il+1, im) :
                    if matches[i][1] in {lMat[1], rMat[1]} :
                        if (i, 4) in usedMatches or (i, 5) in usedMatches :
                            skp = 1
                            break
                if skp :
                    continue
                for i in xrange(il+1, im) :
                    if matches[i][1] in {lMat[1], rMat[1]} :
                        usedMatches[(i, 4)] = usedMatches[(i, 5)] = 0
            usedGroups.append(grp)
            usedMatches[(grp[4], 4)] = usedMatches[(grp[-1], 5)] = 1
            if grp[3] > 0 :
                usedMatches[(grp[4], 5)] = usedMatches[(grp[-1], 4)] = 1

        usedGroups.sort(key=itemgetter(4), reverse=True)
        for gId in xrange(len(usedGroups)-1) :
            g1, g2 = usedGroups[gId:gId+2]
            if g1[4] == g2[-1] :
                m = matches[g1[4]]
                score = g1[0] + g2[0] - m[11]
                length = g1[2] + g2[2] - (m[7]-m[6]+1)
                iden = (g1[1]*g1[2] + g2[1]*g2[2] - min(g1[1],g2[1])*(m[7]-m[6]+1))/length
                usedGroups[gId+1] = [score, iden, length, 0, g2[4]] + g1[4:]
                g1[1] = -1
    else :
        usedGroups = groups
        usedMatches = {(k, k): 1 for k in np.arange(matches.shape[0])}
    for g in usedGroups :
        if g[1] >= 0 :
            ids = [matches[i][15] for i in g[4:]]
            for i in g[4:] :
                matches[i, -1] = g[:3] + ids
    ids = { k[0] for k, v in usedMatches.items() if v == 1 }
    matches = matches[np.array(list(ids))]
    return matches


def cigar2score(data) :
    cigar, rSeq, qSeq, frame, mode, gapOpen, gapExtend = data
    frame = (frame-1) % 3
    gap, rBlk, qBlk = [], [], []
    qId, rId = 0, 0
    for n, t in cigar :
        if t == 'M' :
            rBlk.append(rSeq[rId:rId+n])
            qBlk.append(qSeq[qId:qId+n])
            rId, qId = rId + n, qId + n
        else :
            if t == 'D' :
                gap.append(n) 
                rId += n
            elif t == 'I' :
                gap.append(n) 
                if mode > 1 :
                    qBlk.append(qSeq[qId:qId+n])
                    rBlk.append([-1]*n)
                qId += n
    nGap, bGap = len(gap), np.sum(gap)
    qAln = np.concatenate(qBlk)
    rAln = np.concatenate(rBlk)
    if mode == 1 :
        nMatch = np.sum(qAln == rAln)
        nMismatch = qAln.size - nMatch
        return float(nMatch)/(nMatch + nMismatch+bGap), nMatch*2 - nMismatch*2 - nGap*(gapOpen-gapExtend) - bGap*gapExtend
    else :
        qAln, rAln = qAln[frame:], rAln[frame:]
        if qAln.size % 3 :
            qAln, rAln = qAln[:-(qAln.size % 3)], rAln[:-(qAln.size % 3)]
        qAln, rAln = qAln.reshape(-1, 3), rAln.reshape(-1, 3)
        if mode == 3 :
            match = (qAln == rAln)
            nMatch = np.sum(np.sum(match, 0) * (9./7., 9./7., 3./7.))
            nMismatch = np.sum(rAln >= 0) - nMatch
            return float(nMatch)/(nMatch + nMismatch+bGap), nMatch*2 - nMismatch*2 - nGap*(gapOpen-gapExtend) - bGap*gapExtend
        else :
            s = (~np.any(rAln < 0, 1), )
            qAln, rAln = qAln[s], rAln[s]
            qCodon = np.sum(qAln * (25,5,1), 1)
            rCodon = np.sum(rAln * (25,5,1), 1)
            qAA, rAA = gtable[qCodon], gtable[rCodon]
            nMatch = np.sum(qAA == rAA)*3.
            nTotal = qAA.size * 3. + bGap
            score = np.sum(blosum62[(qAA << 5) + rAA])
            return nMatch/nTotal, score - nGap*(gapOpen-gapExtend) - bGap*gapExtend
nucEncoder = np.zeros(255, dtype=int)
nucEncoder[:] = 2
nucEncoder[(np.array(['A', 'C', 'G', 'T']).view(asc2int),)] = (0, 1, 3, 4)
gtable = np.array(list('KNXKNTTXTTXXXXXRSXRSIIXMIQHXQHPPXPPXXXXXRRXRRLLXLLXXXXXXXXXXXXXXXXXXXXXXXXXEDXEDAAXAAXXXXXGGXGGVVXVVXYXXYSSXSSXXXXXXCXWCLFXLF')).view(asc2int).astype(int)-65

def poolBlast(params) :
    def parseBlast(fin, min_id, min_cov) :
        blastab = pd.read_csv(fin, sep='\t',header=None)
        blastab[2] /= 100.
        blastab = blastab[(blastab[2] >= min_id) & (blastab[7]-blastab[6]+1 >= min_cov) ]
        blastab[14] = np.array(list(map(getCIGAR, zip(blastab[15], blastab[14]))))
        blastab = blastab.drop(columns=[15])
        return blastab
    
    blastn, refDb, qry, min_id, min_cov = params
    blast_cmd = '{blastn} -db {refDb} -query {qry} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen qseq sseq" -num_alignments 1000 -task blastn -evalue 1e-2 -dbsize 5000000 -reward 2 -penalty -3 -gapopen 6 -gapextend 2'.format(
        blastn=blastn, refDb=refDb, qry=qry)
    blastab = parseBlast(Popen(blast_cmd, stdout=PIPE, shell=True, universal_newlines=True).stdout, min_id, min_cov)
    return blastab



def getCIGAR(data) :
    ref, qry = data
    if qry.find('-') < 0 and ref.find('-') < 0 :
        cigar = [[len(qry), 'M']] 
    else :
        tag = np.array(['M', 'I', 'D'])
        cigar = np.concatenate([[-1], (np.array(list(qry)) == '-')*2 + (np.array(list(ref)) == '-'), [-1]])
        pos = np.where(np.diff(cigar) != 0)[0]
        cigar = [ list(v) for v in zip(np.diff(pos), tag[cigar[pos[:-1]+1]]) ]
    return cigar


class RunBlast(object) :
    def __init__(self) :
        self.qrySeq = self.refSeq = None
    def run(self, ref, qry, methods, min_id, min_cov, n_thread=8, re_score=0, filter=[False, 0.9, 0.], linear_merge=[False, 300.,1.2], return_overlap=[True, 300, 0.6], fix_end=[6., 6.]) :
        tools = dict(blastn=self.runBlast, ublast=self.runUBlast, ublastself=self.runUblastSELF, minimap=self.runMinimap, minimapasm=self.runMinimapASM, mmseq=self.runMMseq)
        self.min_id = min_id
        self.min_cov = min_cov
        self.n_thread = n_thread
        self.pool = ThreadPool(n_thread)
        blastab = []
        self.dirPath = tempfile.mkdtemp(prefix='NS_', dir='.')
        try :
            for method in methods :
                if method.lower() in tools :
                    blastab.append(tools[method.lower()](ref, qry))
        finally :
            shutil.rmtree(self.dirPath)
            if blastab :
                blastab = pd.concat(blastab, ignore_index=True)
                blastab = np.hstack([blastab.values, np.arange(blastab.shape[0], dtype=int)[:, np.newaxis]])
        
        if re_score :
            blastab=self.reScore(ref, qry, blastab, re_score)
        if filter[0] :
            blastab=self.ovlFilter(blastab, filter)
        if linear_merge[0] :
            blastab=self.linearMerge(blastab, linear_merge)
        self.fixEnd(blastab, *fix_end)
        if return_overlap[0] :
            overlap = self.returnOverlap(blastab, return_overlap)
            blastab = pd.DataFrame(blastab).sort_values([0,1,11]).values
            return blastab, overlap
        else :
            blastab = pd.DataFrame(blastab).sort_values([0,1,11]).values
            return blastab
            
    def returnOverlap(self, blastab, param) :
        logger('Calculate overlaps.')
        
        ovl_l, ovl_p = param[1:]
        contigs = { tab[1]:id for id, tab in enumerate(blastab) }
        tabs = [ [contigs[tab[1]], tab[15]] + sorted( [tab[8], tab[9]] ) for tab in blastab ]
        tabs = np.array(sorted(tabs, key=itemgetter(0, 2, 3)), dtype=int)
        overlaps = np.empty(shape=[1000001, 3], dtype=int)
        overlaps[-1, :] = [0, 1, -1]
        res = []
        while overlaps[-1, 0] >= 0 :
            logger('Searching {0} / {1} tabs'.format(overlaps[-1, 0], len(tabs)))
            overlaps[:-1, :] = -1
            overlaps = tab2overlaps(tabs, ovl_l, ovl_p, len(tabs), overlaps)
            res.append(overlaps[overlaps.T[2] > 0][:])
        res = np.vstack(res)
        logger('Identified {0} overlaps.'.format(len(res)))
        return res
    
    def reScore(self, ref, qry, blastab, mode, perBatch=10000) :
        if not self.qrySeq :
            self.qrySeq, self.qryQual = readFastq(qry)
        if not self.refSeq :
            self.refSeq, self.refQual = readFastq(ref)
        for k, v in self.qrySeq.items() :
            self.qrySeq[k] = nucEncoder[np.array(list(v)).view(asc2int)]
        for k, v in self.refSeq.items() :
            self.refSeq[k] = nucEncoder[np.array(list(v)).view(asc2int)]
        
        nTab = len(blastab)
        for bId in xrange(0, blastab.shape[0], perBatch) :
            logger('Update scores: {0} / {1}'.format(bId, nTab))
            tabs = blastab[bId:bId+perBatch]
            #scores = np.array([ cigar2score([t[14], self.refSeq[str(t[1])][t[8]-1:t[9]] if t[8] < t[9] else 4 - self.refSeq[str(t[1])][t[9]-1:t[8]][::-1], self.qrySeq[str(t[0])][t[6]-1:t[7]], t[6], mode, 6, 1]) for t in tabs ])
            scores = np.array(list(map(cigar2score, ( [t[14], self.refSeq[str(t[1])][t[8]-1:t[9]] if t[8] < t[9] else 4 - self.refSeq[str(t[1])][t[9]-1:t[8]][::-1], self.qrySeq[str(t[0])][t[6]-1:t[7]], t[6], mode, 6, 1] for t in tabs ))))
            tabs.T[2], tabs.T[11] = scores.T
        return blastab

    def ovlFilter(self, blastab, params) :
        coverage, delta = params[1:]
        logger('Run filtering. Start with {0} hits.'.format(len(blastab)))
        blastab[blastab.T[8] > blastab.T[9], 8:10] *= -1

        blastab = pd.DataFrame(blastab).sort_values(by=[1,0,8,6]).values
        for i, t1 in enumerate(blastab) :
            if t1[2] < 0 : continue
            toDel = []
            for j in xrange(i+1, blastab.shape[0]) :
                t2 = blastab[j]
                if t2[2] < 0 : continue
                if np.any(t1[:2] != t2[:2]) or t1[9] < t2[8] :
                    break
                c = min(t1[9], t2[9]) - t2[8] + 1
                if (c >= coverage*(t1[9]-t1[8]+1) and t2[11] - t1[11] >= delta) :
                    t1[2] = -1.
                    break
                elif (c >= coverage*(t2[9]-t2[8]+1) and t1[11] - t2[11] >= delta) :
                    toDel.append(j)
                elif c >= (t1[9]-t1[8]+1) and c < coverage*(t2[9]-t2[8]+1) :
                    c2 = min(t1[7], t2[7]) - max(t2[6], t1[6]) + 1
                    if c2 >= (t1[7]-t1[6]+1) and c2 < coverage*(t2[7]-t2[6]+1) :
                        t1[2] == -1
                        break
                elif c >= (t2[9]-t2[8]+1) and c < coverage*(t1[9]-t1[8]+1) :
                    c2 = min(t1[7], t2[7]) - max(t2[6], t1[6]) + 1
                    if c2 >= (t2[7]-t2[6]+1) and c2 < coverage*(t1[7]-t1[6]+1) :
                        toDel.append(j)
            if t1[2] >= 0 :
                for j in toDel :
                    blastab[j][2] = -1.
        blastab = blastab[blastab.T[2] >= 0]
        blastab[blastab.T[8] < 0, 8:10] *= -1
        logger('Done filtering. End with {0} hits.'.format(blastab.shape[0]))
        return blastab
    def linearMerge(self, blastab, params) :
        logger('Start merging neighboring regions.')
        blastab[blastab.T[8] > blastab.T[9], 8:10] *= -1
        blastab = pd.DataFrame(blastab).sort_values([0,1,8,6]).values
        blastab = np.vstack(list(map(_linearMerge, [[matches, params] for matches in np.split(blastab, np.where(np.diff(np.unique(blastab.T[0], return_inverse=True)[1]))[0]+1 )])))
        blastab[blastab.T[8] < 0, 8:10] *= -1
        logger('Finish merging neighboring regions.')
        return blastab

    def fixEnd(self, blastab, se, ee) :
        for p in blastab :
            e1, e2 = p[6] - 1, p[12] - p[7]
            cigar = p[14]
            if p[9] > p[8] :
                if 0 < e1 <= se :
                    d = min(p[6]-1, p[8]-1)
                    p[6], p[8], cigar[0][0] = p[6]-d, p[8]-d, cigar[0][0]+d
                if 0 < e2 <= ee :
                    d = min(p[12]-p[7], p[13]-p[9])
                    p[7], p[9], cigar[-1][0] = p[7]+d, p[9]+d, cigar[-1][0]+d
            else :
                if 0 < e1 <= se :
                    d = min(p[6]-1, p[13]-p[8])
                    p[6], p[8], cigar[0][0] = p[6]-d, p[8]+d, cigar[0][0]+d
                if 0 < e2 <= ee :
                    d = min(p[12]-p[7], p[9]-1)
                    p[7], p[9], cigar[-1][0] = p[7]+d, p[9]-d, cigar[-1][0]+d
            p[14] = ''.join( '{0}{1}'.format(n, t) for n, t in cigar )

    def runBlast(self, ref, qry) :
        logger('Run BLASTn starts')
        if not self.qrySeq :
            self.qrySeq, self.qryQual = readFastq(qry)
        if not self.refSeq :
            self.refSeq, self.refQual = readFastq(ref)
        refDb = refNA = os.path.join(self.dirPath, 'refNA')
        if self.refQual is not None :
            with open(refNA, 'w') as fout :
                for n,s in self.refSeq.items() :
                    fout.write('>{0}\n{1}\n'.format(n, s))
        else :
            refNA = ref
        Popen('{makeblastdb} -dbtype nucl -in {refNA} -out {refDb}'.format(makeblastdb=makeblastdb, refNA=refNA, refDb = refDb).split(), stderr=PIPE, stdout=PIPE, universal_newlines=True).communicate()
        qrys = [ os.path.join(self.dirPath, 'qryNA.{0}'.format(id)) for id in range(self.n_thread)]
        qrySeq = sorted(list(self.qrySeq.items()), key=lambda s:-len(s[1]))
        for id, q in enumerate(qrys) :
            with open(q, 'w') as fout :
                for n, s in qrySeq[id::self.n_thread] :
                    fout.write('>{0}\n{1}\n'.format(n, s))
        blastab = pd.concat(self.pool.map(poolBlast, [ [blastn, refDb, q, self.min_id, self.min_cov] for q in qrys ]))
        logger('Run BLASTn finishes. Got {0} alignments'.format(blastab.shape[0]))
        return blastab

    def runMinimap(self, ref, qry) :
        logger('Run Minimap starts')
        p = Popen('{0} -ct{3} -k13 -w5 -A2 -B3 -O8,16 -E2,1 -r50 -p.001 -N500 -f2000,10000 --end-bonus 5 -n1 -m19 -s40 -g200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(\
            minimap2, ref, qry, self.n_thread).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True)
        blastab = self.__parseMinimap(p.stdout, self.min_id, self.min_cov)
        logger('Run Minimap finishes. Got {0} alignments'.format(blastab.shape[0]))
        return blastab

    def __parseMinimap(self, fin, min_id, min_cov) :
        blastab = []
        for line in fin :
            part = line.strip().split('\t')
            if float(part[9])/float(part[10]) < min_id or float(part[10]) < min_cov :
                continue
            direction = part[4]
            ref_sites = [int(part[2])+1, int(part[3])]
            qry_sites = [int(part[7])+1, int(part[8])] if direction == '+' else [int(part[8]), (int(part[7])+1)]
            if direction == '+' :
                cigar = [[int(n), t] for n, t in re.findall(r'(\d+)([A-Z])', part[-1][5:])]
            else :
                cigar = [[int(n), t] for n, t in reversed(re.findall(r'(\d+)([A-Z])', part[-1][5:]))]
            tab = [part[0], part[5], float(part[9])/float(part[10]), int(part[10]), 0, 0] + ref_sites + qry_sites + [0.0, float(part[13][5:]), int(part[1]), int(part[6]), cigar]
            blastab.append(tab)
        return pd.DataFrame(blastab)

    def runMinimapASM(self, ref, qry) :
        logger('Run MinimapASM starts')        
        p = Popen('{0} -ct{3} --frag=yes -A2 -B8 -O20,40 -E3,2 -r20 -g200 -p.000001 -N5000 -f1000,5000 -n2 -m30 -s30 -z200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(
            minimap2, ref, qry, self.n_thread).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True)
        blastab = self.__parseMinimap(p.stdout, self.min_id, self.min_cov)
        logger('Run MinimapASM finishes. Got {0} alignments'.format(blastab.shape[0]))
        return blastab

    def runUblastSELF(self, ref, qry) :
        return self.runUBlast(ref, qry, nhits=100, frames='F')
    def runUBlast(self, ref, qry, nhits=6, frames='7') :
        logger('Run uBLAST starts')        
        def parseUBlast(fin, refseq, qryseq, min_id, min_cov) :
            blastab = pd.read_csv(fin, sep='\t',header=None)
            blastab[2] /= 100.
            blastab = blastab[blastab[2] >= min_id]
            blastab[3], blastab[4] = blastab[3]*3, blastab[4]*3
            
            qf, rf = blastab[0].str.rsplit(':', 1, expand=True), blastab[1].str.rsplit(':', 1, expand=True)
            if np.all(qf[0].str.isdigit()) :
                qf[0] = qf[0].astype(int)
            if np.all(rf[0].str.isdigit()) :
                rf[0] = rf[0].astype(int)
            blastab[0], qf = qf[0], qf[1].astype(int)
            blastab[1], rf = rf[0], rf[1].astype(int)
            blastab[6], blastab[7] = blastab[6]*3+qf-3, blastab[7]*3+qf-1
            
            blastab[12], blastab[13] = blastab[0].apply(lambda x:len(qryseq[str(x)])), blastab[1].apply(lambda x:len(refseq[str(x)]))
            blastab[14] = [ [[3*vv[0], vv[1]] for vv in v ] for v in map(getCIGAR, zip(blastab[15], blastab[14]))]
            
            rf3 = (rf <= 3)
            blastab.loc[rf3, 8], blastab.loc[rf3, 9] = blastab.loc[rf3, 8]*3+rf[rf3]-3, blastab.loc[rf3, 9]*3+rf[rf3]-1
            blastab.loc[~rf3, 8], blastab.loc[~rf3, 9] = blastab.loc[~rf3, 13]-(blastab.loc[~rf3, 8]*3+rf[~rf3]-3-3)+1, blastab.loc[~rf3, 13]-(blastab.loc[~rf3, 9]*3+rf[~rf3]-3-1)+1
            d = np.max([blastab[7] - blastab[12], blastab[9] - blastab[13], 1-blastab[9], np.zeros(blastab.shape[0], dtype=int)], axis=0)
            blastab[7] -= d

            def ending(x, y) :
                x[-1][0] -= y
            np.vectorize(ending)(blastab[14], d)
            d[~rf3] *= -1
            blastab[9] -= d
            blastab = blastab[blastab[7]-blastab[6]+1 >= min_cov].drop(columns=[15,16])
            return blastab
        
        refAA = os.path.join(self.dirPath, 'refAA')
        qryAA = os.path.join(self.dirPath, 'qryAA')
        aaMatch = os.path.join(self.dirPath, 'aaMatch')
        
        if not self.qrySeq :
            self.qrySeq, self.qryQual = readFastq(qry)
        if not self.refSeq :
            self.refSeq, self.refQual = readFastq(ref)

        qryAASeq = transeq(self.qrySeq, frame='F')
        with open(qryAA, 'w') as fout :
            for n, ss in sorted(qryAASeq.items()) :
                _, id, s = min([ (len(s[:-1].split('X')), id, s) for id, s in enumerate(ss) ])
                fout.write('>{0}:{1}\n{2}\n'.format(n, id+1, s))
        
        refAASeq = transeq(self.refSeq, frames)
        toWrite = []
        for n, ss in sorted(refAASeq.items()) :
            for id, s in enumerate(ss) :
                toWrite.append('>{0}:{1}\n{2}\n'.format(n, id+1, s))
        
        blastab = []
        for id in xrange(4) :
            with open(refAA, 'w') as fout :
                for line in toWrite[id::4] :
                    fout.write(line)
        
            ublast_cmd = '{usearch} -threads {n_thread} -db {refAA} -ublast {qryAA} -mid {min_id} -evalue 1 -accel 0.9 -maxhits {nhits} -userout {aaMatch} -ka_dbsize 5000000 -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+raw+ql+tl+qrow+trow+qstrand'.format(
                usearch=usearch, refAA=refAA, qryAA=qryAA, aaMatch=aaMatch, n_thread=self.n_thread, min_id=self.min_id, nhits=nhits)
            p = Popen(ublast_cmd.split(), stderr=PIPE, stdout=PIPE, universal_newlines=True).communicate()
            if os.path.getsize(aaMatch) > 0 :
                blastab.append(parseUBlast(open(aaMatch), self.refSeq, self.qrySeq, self.min_id, self.min_cov))
        blastab = pd.concat(blastab)
        logger('Run uBLAST finishes. Got {0} alignments'.format(blastab.shape[0]))
        return blastab
    
    def runMMseq(self, ref, qry) :
        logger('Run MMSeqs starts')
        def parseMMSeq(fin, refseq, qryseq, min_id, min_cov) :
            blastab = pd.read_csv(fin, sep='\t', header=None)
            blastab = blastab[blastab[2] >= min_id]
            qlen = blastab[0].apply(lambda r:len(qryseq[r]))
            rlen = blastab[1].apply(lambda r:len(refseq[r]))
            cigar = blastab[14].apply( lambda x: [ [int(n)*3, t] for n, t in re.findall(r'(\d+)([A-Z])', x)] )
            ref_sites = pd.concat([ 3*(blastab[6]-1)+1, 3*blastab[7] ], keys=[0,1], axis=1)
            d = ref_sites[1] - qlen
            d[d<0] = 0
            def ending(x, y) :
                x[-1][0] -= y
            np.vectorize(ending)(cigar, d)
            ref_sites[1] -= d

            direction = (blastab[8] < blastab[9])
            qry_sites = pd.concat([ blastab[8], blastab[9]-d ], axis=1)
            qry_sites[~direction] = pd.concat([ blastab[8]-d, blastab[9] ], axis=1)[~direction]
            
            blastab = pd.DataFrame(np.hstack([ blastab[[0, 1, 2]], np.apply_along_axis(lambda x:x[1]-x[0]+1, 1, ref_sites.values)[:, np.newaxis], pd.DataFrame(np.zeros([blastab.shape[0], 2], dtype=int)), ref_sites, qry_sites, blastab[[10, 11]], qlen[:, np.newaxis], rlen[:, np.newaxis], cigar[:, np.newaxis] ]))
            return blastab[blastab[3] >= min_cov]

        tmpDir = os.path.join(self.dirPath, 'tmp')
        refNA = os.path.join(self.dirPath, 'refNA')
        qryNA = os.path.join(self.dirPath, 'qryNA')

        refCDS = os.path.join(self.dirPath, 'refCDS')
        qryAA = os.path.join(self.dirPath, 'qryAA')
        aaMatch = os.path.join(self.dirPath, 'aaMatch2')

        Popen('{0} createdb {1} {2} --dont-split-seq-by-len'.format(mmseqs, ref, refNA).split(), stdout=PIPE).communicate()
        Popen('{0} createdb {1} {2} --dont-split-seq-by-len'.format(mmseqs, qry, qryNA).split(), stdout=PIPE).communicate()
        Popen('{0} translatenucs {1} {2}'.format(mmseqs, qryNA, qryAA).split(), stdout=PIPE).communicate()
        for ite in range(9) :
            if os.path.isdir(tmpDir) :
                shutil.rmtree(tmpDir)
            p = Popen('{0} search {1} {2} {3} {4} -a --alt-ali 30 -s 6 --translation-table 11 --threads {5} --min-seq-id 0.5 -e 10'.format(\
                mmseqs, qryAA, refNA, aaMatch, tmpDir, self.n_thread).split(), stdout=PIPE)
            p.communicate()
            if p.returncode == 0 :
                break
            if ite > 2 :
                Popen('{0} extractorfs {2} {3}'.format(mmseqs, qryAA, refNA, refCDS).split(), stdout=PIPE).communicate()
                p = Popen('{0} search {1} {2} {3} {4} -a --alt-ali 30 -s 6 --translation-table 11 --threads {5} --min-seq-id 0.5 -e 10'.format(\
                    mmseqs, qryAA, refCDS, aaMatch, tmpDir, self.n_thread).split(), stdout=PIPE)
                p.communicate()
                if p.returncode == 0 :
                    break
            time.sleep(1)
        Popen('{0} convertalis {1} {2} {3} {3}.tab --threads {4} --format-output'.format(\
            mmseqs, qryAA, refNA, aaMatch, self.n_thread).split() + ['query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,raw,qlen,tlen,cigar'], stdout=PIPE).communicate()
        
        if not self.qrySeq :
            self.qrySeq, self.qryQual = readFastq(qry)
        if not self.refSeq :
            self.refSeq, self.refQual = readFastq(ref)
        blastab = parseMMSeq(open(aaMatch+'.tab'), self.refSeq, self.qrySeq, self.min_id, self.min_cov)
        logger('Run MMSeqs finishes. Got {0} alignments'.format(blastab.shape[0]))
        return blastab


def uberBlast(args) :
    import argparse
    parser = argparse.ArgumentParser(description='Five different alignment methods. ')
    parser.add_argument('-r', '--reference',     help='[INPUT; REQUIRED] filename for the reference. This is normally a genomic assembly. ', required=True)
    parser.add_argument('-q', '--query',  help='[INPUT; REQUIRED] filename for the query. This can be short-reads or genes or genomic assemblies. ', required=True)
    parser.add_argument('-o', '--output', help='[OUTPUT; Default: None] save result to a file or to screen (stdout). Default do nothing. ', default=None)
    parser.add_argument('--blastn',       help='Run BLASTn. Slowest. Good for identities between [70, 100]', action='store_true', default=False)
    parser.add_argument('--ublast',       help='Run uBLAST on tBLASTn mode. Fast. Good for identities between [30-100]', action='store_true', default=False)
    parser.add_argument('--ublastSELF',   help='Run uBLAST on tBLASTn mode. Fast. Good for identities between [30-100]', action='store_true', default=False)
    parser.add_argument('--minimap',      help='Run minimap. Fast. Good for identities between [90-100]', action='store_true', default=False)
    parser.add_argument('--minimapASM',   help='Run minimap on assemblies. Fast. Good for identities between [90-100]', action='store_true', default=False)
    parser.add_argument('--mmseq',        help='Run mmseq2 on tBLASTn mode. Fast. Good for identities between [60-100]', action='store_true', default=False)
    
    parser.add_argument('--min_id', help='[DEFAULT: 0.3] Minimum identity before reScore for an alignment to be kept', type=float, default=0.3)
    parser.add_argument('--min_cov', help='[DEFAULT: 40] Minimum length for an alignment to be kept', type=float, default=40.)

    parser.add_argument('-s', '--re_score', help='[DEFAULT: 0] Re-interpret alignment scores and identities. 0: No rescore; 1: Rescore with nucleotides; 2: Rescore with amino acid; 3: Rescore with codons', type=int, default=0)
    parser.add_argument('-f', '--filter', help='[DEFAULT: False] Remove secondary alignments if they overlap with any other regions', default=False, action='store_true')
    parser.add_argument('--filter_cov', help='[DEFAULT: 0.9] ', default=0.9, type=float)
    parser.add_argument('--filter_score', help='[DEFAULT: 0] ', default=0., type=float)
    parser.add_argument('-m', '--linear_merge', help='[DEFAULT: False] Merge consective alignments', default=False, action='store_true')
    parser.add_argument('--merge_gap', help='[DEFAULT: 300] ', default=300., type=float)
    parser.add_argument('--merge_diff', help='[DEFAULT: 1.2] ', default=1.2, type=float)
    parser.add_argument('-O', '--return_overlap', help='[DEFAULT: False] Report overlapped alignments', default=False, action='store_true')
    parser.add_argument('--overlap_length', help='[DEFAULT: 300] Minimum overlap to report', default=300, type=float)
    parser.add_argument('--overlap_proportion', help='[DEFAULT: 0.6] Minimum overlap proportion to report', default=0.6, type=float)
    parser.add_argument('-e', '--fix_end', help='[FORMAT: L,R; DEFAULT: 0,0] Extend alignment to the edges if the un-aligned regions are <= [L,R] basepairs.', default='0,0')
    parser.add_argument('-t', '--n_thread', help='[DEFAULT: 8] Number of threads to use. ', type=int, default=8)
    
    args = parser.parse_args(args)
    methods = []
    for method in ('blastn', 'ublast', 'ublastSELF', 'minimap', 'minimapASM', 'mmseq') :
        if args.__dict__[method] :
            methods.append(method)
    for opt in ('fix_end',) :
        args.__dict__[opt] = args.__dict__[opt].split(',')
        args.__dict__[opt][-2:] = list(map(float, args.__dict__[opt][-2:]))
    data = RunBlast().run(args.reference, args.query, methods, args.min_id, args.min_cov, args.n_thread, args.re_score, \
                             [args.filter, args.filter_cov, args.filter_score], \
                             [args.linear_merge, args.merge_gap, args.merge_diff], \
                             [args.return_overlap, args.overlap_length, args.overlap_proportion], \
                             args.fix_end)
    if args.output :
        fout = sys.stdout if args.output.upper() == 'STDOUT' else open(args.output, 'w')
        for t in data :
            fout.write ('\t'.join([str(tt) for tt in t ]) + '\n')
        fout.close()
    return data

if __name__ == '__main__' :
    blastab = uberBlast(sys.argv[1:])

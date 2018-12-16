import os, sys, tempfile, time, shutil, numpy as np, pandas as pd, re
from subprocess import Popen, PIPE
from multiprocessing import Pool
from operator import itemgetter
try:
    from .configure import externals, logger, xrange, readFastq, transeq, blosum62, rc, asc2int
except :
    from configure import externals, logger, xrange, readFastq, transeq, blosum62, rc, asc2int

formatdb = externals['formatdb']
blastn = externals['blastn']
ublast = externals['ublast']
mmseqs = externals['mmseqs']
minimap2 = externals['minimap2']


def cigar2score(data) :
    cigar, rSeq, qSeq, frame, mode, gapOpen, gapExtend = data
    frame = (frame-1) % 3
    qSeq, rSeq = np.array(list(qSeq)).view(asc2int), np.array(list(rSeq)).view(asc2int)
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
                    rBlk.append([0]*n)
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
            nMismatch = np.sum(rAln > 0) - nMatch
            return float(nMatch)/(nMatch + nMismatch+bGap), nMatch*2 - nMismatch*2 - nGap*(gapOpen-gapExtend) - bGap*gapExtend
        else :
            s = (~np.any(rAln == 0, 1), )
            qAln, rAln = qAln[s], rAln[s]
            qCodon = np.sum(nucEncoder[qAln] << (4,2,0), 1)
            rCodon = np.sum(nucEncoder[rAln] << (4,2,0), 1)
            qCodon[qCodon < 0] = 50
            rCodon[rCodon < 0] = 50
            qAA, rAA = gtable[qCodon], gtable[rCodon]
            nMatch = np.sum(qAA == rAA)*3.
            nTotal = qAA.size * 3. + bGap
            score = np.sum(blosum62[(qAA << 5) + rAA])
            return nMatch/nTotal, score - nGap*(gapOpen-gapExtend) - bGap*gapExtend

def poolBlast(params) :
    blastn, refDb, qry, min_id, min_cov = params
    blast_cmd = '{blastn} -db {refDb} -query {qry} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen qseq sseq" -task blastn -evalue 1e-3 -dbsize 5000000 -reward 2 -penalty -3 -gapopen 6 -gapextend 2'.format(
        blastn=blastn, refDb=refDb, qry=qry)
    blastab = parseBlast(Popen(blast_cmd, stdout=PIPE, shell=True, universal_newlines=True).stdout, min_id, min_cov)
    return blastab


nucEncoder = np.empty(255, dtype=int)
nucEncoder.fill(-100)
nucEncoder[(np.array(['A', 'C', 'G', 'T']).view(asc2int),)] = (0, 1, 2, 3)
gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF -')).view(asc2int).astype(int)-65

def parseBlast(fin, min_id, min_cov) :
    blastab = []
    for line in fin :
        part = line.strip().split('\t')
        part[3:10] = list(map(int, part[3:10]))
        part[2] = float(part[2])/100.
        if part[2] <= min_id or part[7]-part[6]+1 <= min_cov :
            continue
        part[12:14] = list(map(int, part[12:14]))
        part[11] = float(part[11])
        part[14] = getCIGAR(part[15], part[14])
        blastab.append(part[:15] + [0])
    return blastab

def getCIGAR(ref, qry) :
    if qry.find('-') < 0 and ref.find('-') < 0 :
        cigar = [[len(qry), 'M']] 
    else :
        tag = np.array(['M', 'I', 'D'])
        cigar = np.concatenate([[-1], (np.array(list(qry)) == '-')*2 + (np.array(list(ref)) == '-'), [-1]])
        pos = np.where(np.diff(cigar) != 0)[0]
        cigar = [ list(v) for v in zip(np.diff(pos), tag[cigar[pos[:-1]+1]])]
    return cigar


class RunBlast(object) :
    def __init__(self) :
        self.qrySeq = self.refSeq = None
    def run(self, ref, qry, methods, min_id, min_cov, n_thread=8, re_score=0, filter=[False, 0.9, 0.], linear_merge=[False, 300.,1.2], fixEnd=[6., 6.]) :
        tools = dict(blastn=self.runBlast, ublast=self.runUBlast, minimap=self.runMinimap, minimapasm=self.runMinimapASM, mmseq=self.runMMseq)
        self.min_id = min_id
        self.min_cov = min_cov
        self.n_thread = n_thread
        blasttab = []
        self.dirPath = tempfile.mkdtemp(prefix='NS_', dir='.')
        try :
            for method in methods :
                if method.lower() in tools :
                    blasttab.extend(tools[method.lower()](ref, qry))
        finally :
            shutil.rmtree(self.dirPath)
        if re_score :
            blasttab=self.reScore(blasttab, re_score)
        if filter[0] :
            blasttab=self.ovlFilter(blasttab, filter)
        if linear_merge[0] :
            self.linearMerge(blasttab, linear_merge)
        if np.max(fix_end) :
            self.fixEnd(blasttab)
        return blasttab
    
    def reScore(self, blasttab, mode, perBatch=5000) :
        pool = Pool(self.n_thread)
        if not self.qrySeq :
            self.qrySeq, self.qryQual = readFastq(qry)
        if not self.refSeq :
            self.refSeq, self.refQual = readFastq(ref)
        nTab = len(blasttab)
        for bId in xrange(0, len(blasttab), perBatch) :
            logger('Update scores: {0} / {1}'.format(bId, nTab))
            tabs = blasttab[bId:bId+perBatch]
            scores = pool.map(cigar2score, [ [t[14], self.refSeq[t[1]][t[8]-1:t[9]] if t[8] < t[9] else rc(self.refSeq[t[1]][t[9]-1:t[8]]), self.qrySeq[t[0]][t[6]-1:t[7]], t[6], mode, 6, 1] for t in tabs ])
            for t, s in zip(tabs, scores) :
                t[2], t[11] = s[0], s[1]
        return blasttab

    def ovlFilter(self, blastab, params) :
        coverage, delta = params[1:]
        logger('Run filtering. Start with {0} hits.'.format(len(blastab)))
        for t in blastab :
            if t[8] > t[9] :
                t[8:10] = -t[8], -t[9]
        blastab.sort(key=itemgetter(1,0,8,6,15))
        nTab = len(blastab)
        for i, t1 in enumerate(blastab) :
            if t1[2] < 0 : continue
            toDel = []
            for j in xrange(i+1, nTab) :
                t2 = blastab[j]
                if t2[2] < 0 : continue
                if t1[:2] != t2[:2] or t1[9] < t2[8] :
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
        blastab = [b for b in blastab if b[2] >= 0]
        for t in blastab :
            if t[8] < 0 :
                t[8:10] = -t[8], -t[9]
        logger('Done filtering. End with {0} hits.'.format(len(blastab)))
        return blastab
    def linearMerge(self, blastab, params) :
        for t in blastab :
            if t[8] > t[9] :
                t[8:10] = -t[8], -t[9]
        genes = {}
        for t in blastab :
            if t[0] not in genes :
                genes[t[0]] = [t]
            else :
                genes[t[0]].append(t)
        return dict(list(map(self.__linearMerge, [[matches, params] for matches in genes.values()])))
    
    def __linearMerge(self, data) :
        matches, params = data
        gapDist, lenDiff = params[1:]
        gene, geneLen = matches[0][0], matches[0][12]
        #min_match_len = [min(params['clust_match_prop']*ref_len, max(params['match_prop'] * ref_len, params['match_len']), max(params['match_prop2'] * ref_len, params['match_len2'])), 
                         #min(params['clust_match_prop']*ref_len, max(params['match_prop2'] * ref_len, params['match_len2'])) ]
        
        def resolve_edges(edges) :
            grps = []
            for m1 in edges[0] :
                for m2 in edges[1] :
                    new_len = m2[5] - m1[4] + 1
                    if new_len >= min_match_len[1] :
                        if m1[2] == m2[2] :
                            s1, e1 = (m1[7], m1[8]) if m1[7] > 0 else (-m1[8], -m1[7])
                            s2, e2 = (m2[7], m2[8]) if m2[7] > 0 else (-m2[8], -m2[7])
                            overlap = min(e1, e2) - max(s1, s2) + 1
                            if overlap >= 0.5*params['synteny_ovl_len'] or overlap >= 0.5*params['synteny_ovl_prop']*min(e1-s1+1, e2-s2+1) :
                                continue
                        overlap = m1[5]-m2[4]+1
                        if overlap >= params['synteny_ovl_len'] or overlap <= -params['synteny_gap'] or overlap >= params['synteny_ovl_prop']*min(m1[5]-m1[4]+1, m2[5]-m2[4]+1) or -overlap >= min(m1[5]-m1[4]+1, m2[5]-m2[4]+1) :
                            continue
                        if overlap > 0 :
                            score = m1[8] + m2[8] - overlap * min( float(m1[8])/(m1[5]-m1[4]+1), float(m2[8])/(m2[5]-m2[4]+1) )
                            ident = (m1[3]*(m1[5]-m1[4]+1) + m2[3]*(m2[5]-m2[4]+1) - overlap*min(m1[3], m2[3]))/((m1[5]-m1[4]+1)+(m2[5]-m2[4]+1) - overlap)
                        else :
                            score = m1[8] + m2[8]
                            ident = (m1[3]*(m1[5]-m1[4]+1) + m2[3]*(m2[5]-m2[4]+1))/((m1[5]-m1[4]+1)+(m2[5]-m2[4]+1))
                        grps.append( [new_len, score, ident ] + ([m1[11], m2[11]] if m1[3] < m2[3] else [m2[11], m1[11]]) )
                        m1[12], m2[12] = 1, 1
            return grps
    
        matches.sort(key=itemgetter(1,8, 6))
        groups = {}
        prev, edges = matches[0][1], [[], []]
        nSave = len(matches)
        
        for id, m1 in enumerate(matches) :
            if m1[6] > 20 and ((m1[8] > 0 and m1[8] - 1 <= gapDist) or (m1[8] < 0 and m1[13] + m1[8] <= gapDist)) :   # any hit within the last 150 bps to either end of a scaffold is a potential fragmented gene
                edges[1].append(m1)
            if m1[5] < m1[12]-20 :
                if (m1[8] > 0 and m1[13]-m1[9] <= gapDist) or (m1[8] < 0 and -1-m1[9] <= gapDist) :
                    edges[0].append(m1)
                for jd in xrange(id+1, nSave) :
                    m2 = matches[jd]
                    if m1[:2] != m2[:2] or (m1[8] < 0 and m2[8] > 0) or m2[8] - m1[9] -1 >= gapDist :    # maximum 300bps between two continuous hits in the same scaffold
                        break
                    if m1[9] >= m2[9] or m1[6] > m2[6] or m1[7] > m2[7] or m2[6] - m1[7] -1 >= gapDist \
                       or min(m2[7]-m1[6]-1, m2[9]-m1[8]-1)*lenDiff < max(m2[7]-m1[6]-1, m2[9]-m1[8]-1) :
                        continue
                    new_len = m2[7] - m1[6] + 1

                    overlap = max(m1[7]-m2[6]+1, m1[9]-m2[8]+1)   ## to improve
                    if overlap > 0 :
                        score = m1[8] + m2[8] - overlap * min( float(m1[8])/match_len, float(m2[8])/match_len2 )
                        ident = (m1[3]*match_len + m2[3]*match_len2 - overlap*min(m1[3], m2[3]))/(match_len + match_len2 - overlap)
                    else :
                        score = m1[8] + m2[8]
                        ident = (m1[3]*match_len + m2[3]*match_len2)/(match_len + match_len2)
                    groups[m1[1]].append( [ new_len, score, ident ] + ([m1[11], m2[11]] if m1[3] < m2[3] else [m2[11], m1[11]]) )
                    m1[12], m2[12] = 1, 1
        
        #groups[prev].extend(resolve_edges(edges))
    
        max_score = max([g[:3] for genome, group in groups.items() for g in group] + [[0, 1, 100.0]], key=itemgetter(1))
        max_score = [float(max_score[0])/max_score[1], max_score[2]/100.0]
        for genome, group in groups.items() :
            for grp in group :
                grp[:3] = [1, max_score[0] * grp[1] * ((float(grp[0])/ref_len)**2), min(grp[2]/max_score[1], 100.0)]
            group.sort(key=lambda x:x[1], reverse=True)
            group[:] = [grp for grp in group if grp[1]*3. >= group[0][1] or grp[2] >= group[0][2]-10]
        match_inuse = {mm: 1 for mat in groups.values() for m in mat for mm in m[3:]}
        for mat in save :
            if mat[11] not in match_inuse :
                mat[12:] = [0]

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
        Popen('{formatdb} -dbtype nucl -in {refNA} -title {refDb}'.format(formatdb=formatdb, refNA=refNA, refDb = refDb).split(), stderr=PIPE, stdout=PIPE, universal_newlines=True).communicate()
        pool = Pool(self.n_thread)
        qrys = [ os.path.join(self.dirPath, 'qryNA.{0}'.format(id)) for id in range(self.n_thread)]
        qrySeq = sorted(list(self.qrySeq.items()), key=lambda s:-len(s[1]))
        for id, q in enumerate(qrys) :
            with open(q, 'w') as fout :
                for n, s in qrySeq[id::self.n_thread] :
                    fout.write('>{0}\n{1}\n'.format(n, s))
        blastab = []
        for btabs in pool.map(poolBlast, [ [blastn, refDb, q, self.min_id, self.min_cov] for q in qrys ]) :
            blastab.extend(btabs)
        logger('Run BLASTn finishes. Got {0} alignments'.format(len(blastab)))
        return blastab

    def runMinimap(self, ref, qry) :
        logger('Run Minimap starts')
        p = Popen('{0} -ct{3} -k13 -w5 -A2 -B3 -O8,16 -E2,1 -r50 -p.001 -N500 -f2000,10000 --end-bonus 5 -n1 -m19 -s40 -g200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(\
            minimap2, ref, qry, self.n_thread).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True)
        blastab = self.__parseMinimap(p.stdout, self.min_id, self.min_cov)
        logger('Run Minimap finishes. Got {0} alignments'.format(len(blastab)))
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
            tab = [part[0], part[5], float(part[9])/float(part[10]), int(part[10]), 0, 0] + ref_sites + qry_sites + [0.0, float(part[13][5:]), int(part[1]), int(part[6]), cigar, 1]
            blastab.append(tab)
        return blastab

    def runMinimapASM(self, ref, qry) :
        logger('Run MinimapASM starts')        
        p = Popen('{0} -ct{3} --frag=yes -A2 -B8 -O20,40 -E3,2 -r20 -g200 -p.000001 -N5000 -f1000,5000 -n2 -m30 -s30 -Z200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(
            minimap2, ref, qry, self.n_thread).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True)
        blastab = self.__parseMinimap(p.stdout, self.min_id, self.min_cov)
        logger('Run MinimapASM finishes. Got {0} alignments'.format(len(blastab)))
        return blastab


    def runUBlast(self, ref, qry) :
        logger('Run uBLAST starts')        
        def parseUBlast(fin, refseq, qryseq, min_id, min_cov) :
            blastab = []
            for line in fin :
                part = line.strip().split('\t')
                part[3:10] = list(map(int, part[3:10]))
                part[2] = float(part[2])/100.
                if part[2] < min_id : continue
                part[3:5] = part[3]*3, part[4]*3
                
                part[0], rf = part[0].rsplit(':', 1)
                part[1], qf = part[1].rsplit(':', 1)
                part[11], part[12], part[13] = float(part[11]), len(qryseq[part[0]]), len(refseq[part[1]])
                rf, qf = int(rf), int(qf)
                part[6], part[7] = part[6]*3+rf-3, part[7]*3+rf-1
                part[14] = [[3*v[0], v[1]] for v in getCIGAR(part[15], part[14])]
    
                if qf <= 3 :
                    part[8], part[9] = part[8]*3+qf-3, part[9]*3+qf-1
                    d = max(part[7] - part[12], part[9]-part[13])
                    if d > 0 :
                        part[7], part[9], part[14][-1][0] = part[7]-d, part[9]-d, part[14][-1][0]-d
                else :
                    part[8], part[9] = part[13]-(part[8]*3+qf-3-3)+1, part[13]-(part[9]*3+qf-3-1)+1
                    d = max(part[7] - part[12], 1 - part[9])
                    if d > 0 :
                        part[7], part[9], part[14][-1][0] = part[7]-d, part[9]+d, part[14][-1][0]-d
                if part[7] - part[6] + 1 < min_cov : continue
                blastab.append(part[:15] + [2])
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
            for n, ss in qryAASeq.items() :
                _, id, s = min([ (len(s[:-1].split('X')), id, s) for id, s in enumerate(ss) ])
                fout.write('>{0}:{1}\n{2}\n'.format(n, id+1, s))
        refAASeq = transeq(self.refSeq)
        with open(refAA, 'w') as fout :
            for n, ss in refAASeq.items() :
                for id, s in enumerate(ss) :
                    fout.write('>{0}:{1}\n{2}\n'.format(n, id+1, s))
        ublast_cmd = '{ublast} -threads {n_thread} -db {refAA} -ublast {qryAA} -evalue 1e-3 -accel 0.9 -maxhits 30 -userout {aaMatch} -ka_dbsize 5000000 -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+raw+ql+tl+qrow+trow+qstrand'.format(
            ublast=ublast, refAA=refAA, qryAA=qryAA, aaMatch=aaMatch, n_thread=self.n_thread)
        p = Popen(ublast_cmd.split(), stderr=PIPE, stdout=PIPE, universal_newlines=True).communicate()
        blastab = parseUBlast(open(aaMatch), self.refSeq, self.qrySeq, self.min_id, self.min_cov)
        logger('Run uBLAST finishes. Got {0} alignments'.format(len(blastab)))
        return blastab
    
    def runMMseq(self, ref, qry) :
        logger('Run MMSeqs starts')
        def parseMMSeq(fin, refseq, qryseq, min_id, min_cov) :
            blastab = []
            for line in fin :
                part = line.strip().split('\t')
                if float(part[2]) < min_id : continue
                rlen = len(refseq[part[1]])
                direction = '+' if int(part[8]) < int(part[9]) else '-'
                cigar = [ [int(n)*3, t] for n, t in re.findall(r'(\d+)([A-Z])', part[14])]
                
                ref_sites = [3*(int(part[6])-1)+1, 3*(int(part[7]))]
                d = max(ref_sites[1] - len(qryseq[part[0]]), 0)
                cigar[-1][0] -= d
                ref_sites[1] -= d
                if ref_sites[1]-ref_sites[0]+1 < min_cov : continue
                qry_sites = [int(part[8]), int(part[9])-1-d] if direction == '+' else [(rlen-int(part[9])+1), (rlen-int(part[8])+2+d)]
                tab = [part[0], part[1], float(part[2]), 3*(int(part[3])+2 if direction == '+' else abs(int(part[3]))+4), 0, 0, ] + ref_sites + \
                    qry_sites + [float(part[10]), float(part[11]), len(qryseq[part[0]]), rlen, cigar, 3]
                blastab.append(tab)
            return blastab
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
            mmseqs, qryAA, refNA, aaMatch, self.n_thread).split() + ['query target pident alnlen mismatch gapopen qstart qend tstart tend evalue raw qlen tlen cigar'], stdout=PIPE).communicate()
        
        if not self.qrySeq :
            self.qrySeq, self.qryQual = readFastq(qry)
        if not self.refSeq :
            self.refSeq, self.refQual = readFastq(ref)
        blastab = parseMMSeq(open(aaMatch+'.tab'), self.refSeq, self.qrySeq, self.min_id, self.min_cov)
        logger('Run MMSeqs finishes. Got {0} alignments'.format(len(blastab)))        
        return blastab


def uberBlast(args) :
    import argparse
    parser = argparse.ArgumentParser(description='Five different alignment methods. ')
    parser.add_argument('-r', '--reference',     help='[INPUT; REQUIRED] filename for the reference. This is normally a genomic assembly. ', required=True)
    parser.add_argument('-q', '--query',  help='[INPUT; REQUIRED] filename for the query. This can be short-reads or genes or genomic assemblies. ', required=True)
    parser.add_argument('--blastn',       help='Run BLASTn. Slowest. Good for identities between [70, 100]', action='store_true', default=False)
    parser.add_argument('--ublast',       help='Run uBLAST on tBLASTn mode. Fast. Good for identities between [30-100]', action='store_true', default=False)
    parser.add_argument('--minimap',      help='Run minimap. Fast. Good for identities between [90-100]', action='store_true', default=False)
    parser.add_argument('--minimapASM',   help='Run minimap on assemblies. Fast. Good for identities between [90-100]', action='store_true', default=False)
    parser.add_argument('--mmseq',        help='Run mmseq2 on tBLASTn mode. Fast. Good for identities between [60-100]', action='store_true', default=False)
    
    parser.add_argument('--min_id', help='[DEFAULT: 0.3] Minimum identity before reScore for an alignment to be kept', type=float, default=0.3)
    parser.add_argument('--min_cov', help='[DEFAULT: 40] Minimum length for an alignment to be kept', type=float, default=40.)

    parser.add_argument('-s', '--re_score', help='[DEFAULT: 0] Re-interpret alignment scores and identities. 0: No rescore; 1: Rescore with nucleotides; 2: Rescore with amino acid; 3: Rescore with codons', type=int, default=0)
    parser.add_argument('-f', '--filter', help='[DEFAULT: False] Remove secondary alignments if they overlap with any other regions', default=True, action='store_true')
    parser.add_argument('--filter_cov', help='[DEFAULT: 90] ', default=0.9, type=float)
    parser.add_argument('--filter_score', help='[DEFAULT: 0] ', default=0., type=float)
    parser.add_argument('-m', '--linear_merge', help='[DEFAULT: False] Merge consective alignments', default=False, action='store_true')
    parser.add_argument('--merge_gap', help='[DEFAULT: 300] ', default=300., type=float)
    parser.add_argument('--merge_diff', help='[DEFAULT: 1.2] ', default=1.2, type=float)
    parser.add_argument('-e', '--fix_end', help='[FORMAT: L,R; DEFAULT: 0,0] Extend alignment to the edges if the un-aligned regions are <= [L,R] basepairs.', default='0,0')
    parser.add_argument('-t', '--n_thread', help='[DEFAULT: 8] Number of threads to use. ', type=int, default=8)
    
    args = parser.parse_args(args)
    methods = []
    for method in ('blastn', 'ublast', 'minimap', 'minimapASM', 'mmseq') :
        if args.__dict__[method] :
            methods.append(method)
    for opt in ('fix_end',) :
        args.__dict__[opt] = args.__dict__[opt].split(',')
        args.__dict__[opt][-2:] = list(map(float, args.__dict__[opt][-2:]))
    blastab = RunBlast().run(args.reference, args.query, methods, args.min_id, args.min_cov, args.n_thread, args.re_score, [args.filter, args.filter_cov, args.filter_score], [args.linear_merge, args.merge_gap, args.merge_diff], args.fix_end)

if __name__ == '__main__' :
    uberBlast(sys.argv[1:])
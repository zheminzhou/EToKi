import os, sys, numpy as np, tempfile, shutil, re
from subprocess import Popen, PIPE
from operator import itemgetter
try:
    from configure import externals, rc, uopen, xrange, get_md5
except :
    from .configure import externals, rc, uopen, xrange, get_md5

usearch = externals['usearch']
makeblastdb = externals['makeblastdb']
blastn = externals['blastn']


def transeq(seq, frames=[1,2,3,4,5,6]) :
    gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF '))
    if isinstance(seq, dict) :
        names, seqs = [], []
        for n,s in seq.items() :
            names.append(n)
            seqs.append( 'Y' * ((3-len(s)%3)%3) + re.sub(r'[^ACGT]', 'X', s.upper()) + 'Z' * ((3-len(s)%3)%3) )
    else :
        names, seqs = [n for n,s in seq], [re.sub(r'[^ACGT]', 'X', s.upper()) + 'X' * ((3-len(s)%3)%3) for n, s in seq]
    trans_seq = [[n, []] for n in names]
    mseqs = np.unique(list('   '.join(seqs)) + [' ', 'A', 'C', 'G', 'T', 'X', 'Y', 'Z'], return_inverse=True)[1][:-8]
    seqs = mseqs[mseqs != 6]
    rseqs = np.flip(5-mseqs[mseqs != 7], 0)
    seqs[seqs == 0], seqs[seqs >= 5] = -100000, -100
    rseqs[rseqs <= 0], rseqs[rseqs == 5] = -100, -100000
    seqs -= 1; rseqs -= 1
    for frame in frames :
        if frame <= 3 :
            codons = np.concatenate([seqs[frame-1:], np.array([-100]*(frame-1), dtype=int)]).reshape(-1, 3)
        else :
            codons = np.concatenate([rseqs[frame-4:], np.array([-100]*(frame-4), dtype=int)]).reshape(-1, 3)
        codons[:, 0] <<= 4
        codons[:, 1] <<= 2
        codons = np.sum(codons, 1)
        codons[codons < -1000000] = 64
        codons[codons < 0] = 50
        tseq = ''.join(gtable[codons]).split(' ')
        if frame <= 3 :
            for ts, tt in zip(trans_seq, tseq) :
                ts[1].append(tt)
        else :
            for ts, tt in zip(trans_seq, tseq[::-1]) :
                ts[1].append(tt)
    return dict(trans_seq) if isinstance(seq, dict) else trans_seq

class dualBlast(object) :
    def readFasta(self, fasta) :
        sequence = {}
        with open(fasta) as fin :
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    sequence[name] = []
                elif len(line) > 0 :
                    sequence[name].extend(line.strip().split())
        for s in sequence :
            sequence[s] = (''.join(sequence[s])).upper()
        return sequence
    
    def readFastq(self, fastq) :
        sequence, qual = {}, {}
        with open(fastq) as fin :
            line = fin.readline()
            if line.startswith('>') :
                sequence = self.readFasta(fastq)
                return sequence, None
        with open(fastq) as fin :
            for lineId, line in enumerate(fin) :
                if lineId % 4 == 0 :
                    name = line[1:].strip().split()[0]
                    sequence[name] = []
                    qual[name] = []
                elif lineId % 4 == 1 :
                    sequence[name].extend(line.strip().split())
                elif lineId % 4 == 3 :
                    qual[name].extend(line.strip().split())
        for s in sequence :
            sequence[s] = (''.join(sequence[s])).upper()
            qual[s] = ''.join(qual[s])
        return sequence, qual
    
    def getCIGAR(self, ref, qry) :
        if ref.find('-') < 0 and qry.find('-') < 0 :
            cigar = [[len(ref), 'M']] 
        else :
            tag = np.array(['M', 'D', 'I'])
            cigar = np.concatenate([[-1], (np.array(list(ref)) == '-')*2 + (np.array(list(qry)) == '-'), [-1]])
            pos = np.where(np.diff(cigar) != 0)[0]
            cigar = [ list(v) for v in zip(np.diff(pos), tag[cigar[pos[:-1]+1]])]
        return cigar
    
    def parseBlast(self, fin, min_iden, min_len) :
        blastab = []
        for line in fin :
            part = line.strip().split('\t')
            part[3:10] = list(map(int, part[3:10]))
            part[2] = float(part[2])/100.
            if part[2] < min_iden or part[3] < min_len :
                continue
            part[12:14] = list(map(int, part[12:14]))
            part[11] = float(part[11])
            part[14] = self.getCIGAR(part[14], part[15])
            blastab.append(part[:15])
        return blastab
    def parseUBlast(self, fin, qryseq, refseq, min_iden, min_len) :
        blastab = []
        for line in fin :
            part = line.strip().split('\t')
            part[3:10] = list(map(int, part[3:10]))
            part[2] = float(part[2])/100.
            part[3:5] = part[3]*3, part[4]*3
            if part[2] < min_iden or part[3] < min_len :
                continue            
            
            part[0], rf = part[0].rsplit(':', 1)
            part[1], qf = part[1].rsplit(':', 1)
            part[11], part[12], part[13] = float(part[11]), len(refseq[part[0]]), len(qryseq[part[1]])
            rf, qf = int(rf), int(qf)
            part[6], part[7] = part[6]*3+rf-3, part[7]*3+rf-1
            part[14] = [[3*v[0], v[1]] for v in self.getCIGAR(part[14], part[15])]

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
            
            blastab.append(part[:15])
        return blastab
    
    def run_comparison(self, dirPath, qry, ref, min_iden, min_len, n_thread=6) :
        qryNA = os.path.join(dirPath, 'qryNA')
        qryAA = os.path.join(dirPath, 'qryAA')
        refAA = os.path.join(dirPath, 'refAA')
        naMatch = os.path.join(dirPath, 'naMatch')
        aaMatch = os.path.join(dirPath, 'aaMatch')
        
        refSeq = self.readFasta(ref)
        qrySeq, qual = self.readFastq(qry)
        with open(qryNA, 'w') as fout :
            for n,s in qrySeq.items() :
                fout.write('>{0}\n{1}\n'.format(n, s))
        
        Popen('{makeblastdb} -dbtype nucl -in {qry}'.format(makeblastdb=makeblastdb, qry=qryNA).split(), stderr=PIPE, stdout=PIPE).communicate()

        refs = [ [os.path.join(dirPath, 'ref.{0}'.format(id)), os.path.join(dirPath, 'ref.{0}.out'.format(id)), []] for id in range(n_thread)]
        
        refSeq2 = sorted(list(refSeq.items()), key=lambda s:-len(s[1]))
        for id, (r, o, p) in enumerate(refs) :
            with open(r, 'w') as fout :
                for n, s in refSeq2[id::n_thread] :
                    fout.write('>{0}\n{1}\n'.format(n, s))
            blast_cmd = '{blastn} -db {qry} -query {ref} -out {out} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen qseq sseq" -task blastn -evalue 1e-3 -dbsize 5000000 -reward 2 -penalty -2 -gapopen 6 -gapextend 2'.format(
                blastn=blastn, qry=qryNA, ref=r, out=o)
            p.append(Popen(blast_cmd, stdout=PIPE, shell=True))
        
        refAASeq = transeq(refSeq, frames=[1,2,3])
        with open(refAA, 'w') as fout :
            for n, ss in refAASeq.items() :
                for id, s in enumerate(ss) :
                    if len(s[:-1].split('X')) < 2 :
                        fout.write('>{0}:{1}\n{2}\n'.format(n, id+1, s))
        qryAASeq = transeq(qrySeq)
        with open(qryAA, 'w') as fout :
            for n, ss in qryAASeq.items() :
                for id, s in enumerate(ss) :
                    fout.write('>{0}:{1}\n{2}\n'.format(n, id+1, s))
        with open(naMatch, 'w') as fout :
            for r, o, p in refs :
                p[0].communicate()
                fout.write(open(o).read())
                os.unlink(o)
        ublast_cmd = '{usearch} -threads {n_thread} -db {qryAA} -ublast {refAA} -evalue 1e-3 -accel 0.9 -maxhits 6 -userout {aaMatch} -ka_dbsize 5000000 -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+raw+ql+tl+qrow+trow+qstrand'.format(
            usearch=usearch, qryAA=qryAA, refAA=refAA, aaMatch=aaMatch, n_thread=n_thread)
        pp = Popen(ublast_cmd.split(), stderr=PIPE, stdout=PIPE)
        pp.communicate()
        
        blastab = self.parseBlast(open(naMatch), min_iden, min_len)
        blastab.extend(self.parseUBlast(open(aaMatch), qrySeq, refSeq, min_iden, min_len))
        self.fixEnd(blastab, 6, 9)
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

class blastParser(object) :
    def linear_merge(self, blasttab, min_iden, min_frag_prop, min_frag_len, max_dist=300, diag_diff=1.5, max_diff=200, **params) :
        for part in blasttab :
            if part[8] > part[9] :
                part[8], part[9] = -part[8], -part[9]

        blasttab.sort(key=itemgetter(0,1,6,8,-11))
        nB = len(blasttab)
        for id, p1 in enumerate(blasttab) :
            if p1[0] == '' : continue
            for jd in xrange(id+1, nB) :
                p2 = blasttab[jd]
                if p2[0] == '' : continue
                if (p1[0], p1[1]) != (p2[0], p2[1]) or p2[6] - p1[6] > 4 :
                    break
                d = abs(p1[6]-p2[6]) + abs(p1[7]-p2[7]) + abs(p1[8]-p2[8]) + abs(p1[9]-p2[9])
                if d <= 5 :
                    if p1[2] >= p2[2] :
                        p2[0] = ''
                    else :
                        p1[0] = ''
                        break
        blasttab = [p for p in blasttab if p[0] != '']
        
        nB = len(blasttab)
        syntenies = []
        for id, p1 in enumerate(blasttab) :
            for jd in xrange(id+1, nB)  :
                p2 = blasttab[jd]

                if p1[0] != p2[0] or p1[1] != p2[1] or p2[6] - p1[7] > max_dist :
                    break
                elif p1[8] < 0 < p2[8] or p2[8] < p1[8] + 15 or p2[9] < p1[9] + 15 or p2[7]<p1[7]+15 or p2[6] < p1[6]+15 or p2[8] - p1[9] > max_dist :
                    continue
                m, n = p2[7] - p1[6] + 1, p2[9] - p1[8] + 1
                if m < min_frag_len or m < min_frag_prop*p1[12] or max(m, n) - min(m, n) > max_diff or max(m,n) > diag_diff * min(m, n) :
                    continue
                o_len = 0 if p2[6] > p1[7] else p1[7] - p2[6] + 1
                p1_len = p1[7] - p1[6] + 1 - o_len
                p2_len = p2[7] - p2[6] + 1 - o_len
                
                iden = (p1[2]*p1_len + p2[2]*p2_len +max(p1[2], p2[2])*o_len)/(p1_len+p2_len+o_len)
                if iden < min_iden :
                    continue

                p1s, p2s = p1[11]/(p1[7]-p1[6]+1), p2[11]/(p2[7]-p2[6]+1)
                dist = max(p2[6]-p1[7]-1, p2[8]-p1[9]-1, 0)
                score = p1s*p1_len +p2s*p2_len+max(p1s, p2s)*o_len - dist
                if score > 0 :
                    syntenies.append([id, jd, iden, score])
        syn_score = {}
        for id , syn in enumerate(syntenies) :
            if syn[0] not in syn_score and syn[1] not in syn_score :
                p1, p2 = blasttab[syn[0]][:], blasttab[syn[1]]
                c1, c2 = p1[-1], p2[-1]                
                r_dist, q_dist = p2[6]-p1[7]-1, p2[8]-p1[9]-1
                if min(r_dist, q_dist) < 0 :
                    p1s, p2s = p1[11]/(p1[7]-p1[6]+1), p2[11]/(p2[7]-p2[6]+1)
                    if p1s <= p2s :
                        cc = [ [int(n), t] for n,t in re.findall(r'(\d+)([A-Z])', c1)]
                        i = -1
                    else :
                        cc = [ [int(n), t] for n,t in re.findall(r'(\d+)([A-Z])', c2)]
                        i = 0

                    while min(r_dist, q_dist) < 0 :
                        if cc[i][1] == 'M' :
                            d = min(cc[i][0], max(-r_dist, -q_dist, 0))
                            r_dist, q_dist = d + r_dist, d + q_dist
                        elif cc[i][1] == 'D' :
                            d = min(cc[i][0], max(-r_dist, -q_dist, 0))
                            r_dist += d
                        elif cc[i][1] == 'I' :
                            d = min(cc[i][0], max(-r_dist, -q_dist, 0))
                            q_dist += d
                        else :
                            raise 'unknown'
                        if d >= cc[i][0] :
                            cc = cc[1:] if i == 0 else cc[:-1]
                        else :
                            cc[i][0] -= d
                    if i == -1 :
                        c1 = ''.join([ '{0}{1}'.format(*c) for c in cc ])
                    else :
                        c2 = ''.join([ '{0}{1}'.format(*c) for c in cc ])
                gap=[]
                if r_dist > 0 :
                    gap.append('{0}D'.format(r_dist))
                if q_dist > 0 :
                    gap.append('{0}I'.format(q_dist))
                if len(gap) == 0 :
                    cc1 = re.findall(r'(^.*?)(\d+)([A-Z]$)', c1)[0]
                    cc2 = re.findall(r'(^\d+)([A-Z])(.*$)', c2)[0]
                    if cc1[2] == cc2[1] :
                        c1 = '{0}{1}{2}'.format(cc1[0], int(cc1[1])+int(cc2[0]), cc1[2])
                        c2 = cc2[2]
                p1[7], p1[9], p1[2], p1[11], p1[14] = p2[7], p2[9], syn[2], syn[3], ''.join([c1] + gap + [c2])
                blasttab.append(p1)
            syn_score[syn[0]] = id
            syn_score[syn[1]] = id
            
        for part in blasttab :
            if part[8] < 0 :
                part[8], part[9] = -part[8], -part[9]
            x = ['{0}D'.format(part[6]-1), part[-1]] if part[6] > 1 else [part[-1]]
            if part[7] < part[12] :
                x.append('{0}D'.format(part[12]-part[7]))
            if len(x) : part[-1] = ''.join(x)
        return blasttab
    def parse_blast(self, hits, parameters) :
        for part in hits:
            allele = part[0].rsplit('_', 1)
            r_start, r_end = 1, part[12]
            direct = 1 if part[8] < part[9] else -1
            tailing = min(r_start - part[6], part[7] - r_end)

            q_start, q_end = part[8], part[9]

            ## ignore the hit if it is too small
            if direct*(q_end - q_start) + 1 < parameters['min_frag_prop'] * (r_end - r_start + 1) or direct*(q_end - q_start)+1 < parameters['min_frag_len'] :
                part[0] = ''
            else :
                part.extend([r_start, r_end, q_start, q_end, float(part[11]), allele[0], allele[1], tailing])
                #              15       16      17      18        19            20          21        22
        hits = sorted([hit for hit in hits if hit[0] != ''], key = itemgetter(22, 19), reverse=True)
        alleles = {'__non_specific__':[]}

        for part in hits:
            ## check whether the allele has overlapped with other better hits
            if part[19] <= 0 : continue
            allele_full = '{0}_{1}'.format(*part[20:])
            r_start, r_end, q_start, q_end = part[15:19]
            overlap, to_move = -1, []

            if part[20] not in alleles :
                alleles[part[20]] = []
            else :
                s, e = sorted([q_start, q_end])
                for id, region in enumerate( alleles[part[20]] ):
                    if region[0] == '' : continue
                    if region[2] == part[1] and (min(e, region[4]) - max(s, region[3]) + 1) >= parameters['merging_prop']*(e-s+1) :
                        if region[1] < part[2] and (min(e, region[4]) - max(s, region[3]) + 1) >= parameters['merging_prop']*(region[4]-region[3]+1) :
                            if region[1] < parameters['min_iden'] or (part[22] >= 0 and part[2] - region[1] >= parameters['merging_error']):
                                to_move.append(id)
                                continue
                            else :
                                region[1] = part[2]
                        overlap = id

                    elif region[1] - part[2] >= parameters['merging_error'] :
                        overlap = 999999999
                        break
                    elif part[2] - region[1] >= parameters['merging_error'] :
                        to_move.append(id)
                if overlap < 999999999 :
                    for id in reversed(to_move) :
                        alleles['__non_specific__'].append( alleles[part[20]].pop(id)[:-1] + [''] )



            if overlap < 0 :
                # insert a new region if there is no overlap
                if q_start < q_end :
                    alleles[part[20]].append([part[20], float(part[2]), part[1], q_start, q_end, '+', r_start - part[6], part[7]-r_end, '', allele_full+':'+part[14]])
                else :
                    alleles[part[20]].append([part[20], float(part[2]), part[1], q_end, q_start, '-', r_start - part[6], part[7]-r_end, '', allele_full+':'+part[14]])
            elif overlap == 999999999 :
                if q_start < q_end :
                    alleles['__non_specific__'].append([part[20], float(part[2]), part[1], q_start, q_end, '+', r_start - part[6], part[7]-r_end, '', ''])#part[0]])
                else :
                    alleles['__non_specific__'].append([part[20], float(part[2]), part[1], q_end, q_start, '-', r_start - part[6], part[7]-r_end, '', ''])#part[0]])
        return alleles

    def inter_loci_overlap(self, alleles, parameters) :
        regions = [reg for region in alleles.values() for reg in region]
        # sort with contig name and start points
        regions.sort(key=itemgetter(2,3))

        for id, regi in enumerate(regions) :
            if regi[0] == '' : continue
            todel, deleted = [], 0
            for jd in xrange(id+1, len(regions)) :
                regj = regions[jd]
                if regj[0] == '' or regi[0] == regj[0]: continue
                if regi[2] != regj[2] or regj[3] > regi[4] :
                    break
                overlap = min(regi[4], regj[4]) - regj[3] + 1
                if (regi[-1] != '' and float(overlap) >= parameters['merging_prop'] * (regi[4]-regi[3]+1)) or \
                   (regj[-1] != '' and float(overlap) >= parameters['merging_prop'] * (regj[4]-regj[3]+1)) :
                    delta = regi[1] - regj[1]
                    if delta > 0.05 :
                        todel.append(jd)
                    elif delta <= -0.05 :
                        deleted = 1
                        break
            if deleted == 0 :
                for jd in todel:
                    regions[jd][0] = ''
            else :
                regi[0] = ''

        return [{'locus':reg[0], 'identity':reg[1], 'CIGAR':reg[9], 'coordinates':[reg[2], int(reg[3]), int(reg[4]), reg[5]], 'flanking':reg[6:8], 'status':reg[8], 'accepted':(0 if reg[8] == '' else 128)} for reg in regions if reg[0] != '' and reg[-1] != '']

    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    def get_seq(self, seq, header, start, end, direction) :
        fasta = seq[header]
        if direction == '+' :
            seq = fasta[(start-1):end].upper()
        else :
            seq = ''.join([self.complement.get(x.upper(), 'N') for x in reversed(fasta[(start-1):end])])
        return seq
    def get_qual(self, qual, header, start, end, direction='+', force=False) :
        if qual == None:
            return 20 if force else 0
        try:
            s = max(1, min(start, end))
            e = min(len(qual[header]), max(start, end))
            return ord(min(qual[header][(s-1):e])) - 33
        except :
            return 0

    def lookForORF(self, seq, rec) :
        coordinates, edges = rec['coordinates'], rec['flanking'][:]
        seq = self.get_seq(seq, *coordinates)
        #if (len(seq) - sum(edges)) % 3 == 0 :
        startCodon, stopCodon = 0, 0
        for s in xrange(edges[0]%3, len(seq), 3) :
            c = seq[s:s+3]
            if c in {'ATG', 'TTG', 'GTG'} :
                startCodon, edges[0] = 1, 0
                e0 = s+3
                new_s = coordinates[1] + s if coordinates[3] == '+' else coordinates[2] - s
                break
        if not startCodon :
            return 6
        for e in xrange(e0, len(seq), 3) :
            c = seq[e:e+3]
            if c in {'TAG', 'TAA', 'TGA'} :
                stopCodon, edges[1] = 1, 0
                new_e = coordinates[1] + e+2 if coordinates[3] == '+' else coordinates[2] - e-2
                break
        if startCodon and stopCodon and abs(e-s) + 1 >= 0.6 * (abs(coordinates[2]-coordinates[1])+1 - sum(rec['flanking'])) :
            c2 = (new_s, new_e) if coordinates[3] == '+' else (new_e, new_s)
            if c2[0] != coordinates[1] or c2[1] != coordinates[2] :
                rec['CIGAR'] = rec['CIGAR'].rsplit(':', 1)[0] + ':EXEMPT'
            coordinates[1:3] = c2
            rec['flanking'] = edges
            return 0
        return 6


    def form_alleles(self, regions, qrySeq, qryQual, genome_id, accepted, argument) :
        alleles = {}
        regions.sort(key=lambda x:x['identity'], reverse=True)
        regions.sort(key=lambda x:min(x['flanking'] + [0]), reverse=True)
        for region in regions:
            if sum(region['flanking']) >= -30 and (argument.get('ORF', False) or argument.get('CDS', False)) and region['CIGAR'] != 'intergenic' :
                flag = self.lookForORF(qrySeq, region)
                region['accepted'] = region['accepted'] | flag
            region['seq'] = self.get_seq(qrySeq, *region['coordinates'])
            region['id'] = ''
            region['value_md5'] = get_md5(region['seq'])
            if min(region['flanking']) >= 0 and len(re.findall(r'[^ACGT]', region['seq'])) == 0 :  ## add proportional check
                region['accepted'] = region['accepted'] | 1
            else :
                region['status'] += '{Fragmented}'
                region['accepted'] = region['accepted'] | 64
                region['allele_id'] = -1

            if region['locus'] in alleles :
                if region['accepted'] & 64 > 0 :
                    if alleles[ region['locus'] ]['accepted'] & 64 > 0 :
                        if 'secondary' not in alleles[ region['locus'] ] :
                            alleles[ region['locus'] ]['secondary'] = []
                        alleles[ region['locus'] ]['secondary'].append( dict(coordinates =region['coordinates'], seq=region['seq'], identity=region['identity']) )
                elif alleles[ region['locus'] ] ['accepted'] & 32 == 0 :
                    alleles[ region['locus'] ] ['status'] += '{Duplicated}'
                    alleles[ region['locus'] ] ['seq'] = 'DUPLICATED'
                    alleles[ region['locus'] ] ['value_md5'] = get_md5('DUPLICATED')
                    alleles[ region['locus'] ] ['accepted'] = (alleles[ region['locus'] ] ['accepted'] | 32) & (~1)
                    alleles[ region['locus'] ] ['allele_id'] = -1
                    if 'secondary' not in alleles[ region['locus'] ] :
                        alleles[ region['locus'] ]['secondary'] = []
                    alleles[ region['locus'] ]['secondary'].append( dict(coordinates =region['coordinates'], seq=region['seq'], identity=region['identity']) )
            else :
                if accepted == 0 or self.get_qual(qryQual, *region['coordinates']) < 10:
                    region['accepted'] = region['accepted'] | 2
                region['reference'] = 'MLSType:'+genome_id
                alleles[region['locus']] = region
            if region['accepted'] & 2 > 0 :
                region['accepted'] = region['accepted'] & (~1)
        for locus, allele in alleles.items() :
            if allele['accepted'] & 65 == 64 :
                allele_len = allele['coordinates'][2] - allele['coordinates'][1] + 1
                for ale in allele.get('secondary', {}) :
                    allele_len += ale['coordinates'][2] - ale['coordinates'][1] + 1
                if allele_len < argument['min_frag_prop'] :
                    alleles.pop(locus)
            if 'identity' in allele and allele['identity'] < argument['min_iden'] :
                allele['allele_id'] = -1
                allele['accepted'] = (allele['accepted'] & (~1)) | 256
                allele['status'] += '{Low identities:'+str(allele['identity'])+'}'
                if allele['accepted'] & 224 > 0 :
                    alleles.pop(locus, None)
        return alleles
    def intergenic(self, regions, lenRange) :
        inter_blocks = []
        regions.sort(key=lambda r:r['coordinates'])
        prev = regions[0]
        for region in regions[1:] :
            if prev['coordinates'][0] == region['coordinates'][0] :
                d = region['coordinates'][1] - prev['coordinates'][2] - 1
                if max(lenRange[0], 0) <= d <= lenRange[1] :
                    e1 = prev['flanking'][1] if prev['coordinates'][3] == '+' else prev['flanking'][0]
                    e2 = region['flanking'][0] if region['coordinates'][3] == '+' else region['flanking'][1]
                    if e1 == e2 and e1 == 0 :
                        block = {'status':prev['status']+region['status'], 'flanking':[0,0], 'accepted':0, 'identity':min(prev['identity'], region['identity']), 'CIGAR':'intergenic'}
                        if prev['coordinates'][3] == region['coordinates'][3] :
                            if region['coordinates'][3] == '+' :
                                block['locus'] = '{0}:R:{1}:L'.format(prev['locus'], region['locus'])
                                block['coordinates'] = [prev['coordinates'][0], prev['coordinates'][2]+1, region['coordinates'][1]-1, '+']
                            else :
                                block['locus'] = '{0}:R:{1}:L'.format(region['locus'], prev['locus'])
                                block['coordinates'] = [prev['coordinates'][0], prev['coordinates'][2]+1, region['coordinates'][1]-1, '-']
                        else :
                            d = 'R' if prev['coordinates'][3] == '+' else 'L'
                            
                            if prev['locus'] < region['locus'] :
                                block['locus'] = '{0}:{2}:{1}:{2}'.format(prev['locus'], region['locus'], d)
                                block['coordinates'] = [prev['coordinates'][0], prev['coordinates'][2]+1, region['coordinates'][1]-1, '+']
                            else :
                                block['locus'] = '{0}:{2}:{1}:{2}'.format(region['locus'], prev['locus'], d)
                                block['coordinates'] = [prev['coordinates'][0], prev['coordinates'][2]+1, region['coordinates'][1]-1, '-']
                        inter_blocks.append(block)
            prev = region
        return sorted(regions+inter_blocks, key=lambda r:r['coordinates'])
            

class seqOperation(object) :
    def readSequence(self, fname) :
        seq = {}
        with open(fname) as fin :
            header = fin.read(1)
            fin.seek(0, 0)
            if header == '@' :
                for id, line in enumerate(fin) :
                    if id % 4 == 0 :
                        name = line[1:].strip().split()[0]
                        seq[name]= [None, None]
                    elif id % 4 == 1 :
                        seq[name][0] = line.strip()
                    elif id % 4 == 3 :
                        seq[name][1] = [ord(q)-33 for q in line.strip()]
            else :
                for line in fin :
                    if line.startswith('>') :
                        name = line[1:].strip().split()[0]
                        seq[name] = [[], None]
                    else :
                        seq[name][0].extend( line.strip().split() )
                for n, s in seq.iteritems() :
                    s[0] = ''.join(s[0])
                    s[1] = [0 for ss in s[0]]
        return seq

    def write_refsets(self, reference) :
        ref_aa = '{0}.refset.aa'.format(parameters['unique_key'])
        refseq = self.readSequence(reference)
        refamino = transeq({n:s[0] for n,s in refseq.iteritems() }, 1)
        with open(ref_aa, 'w') as fout :
            for n, s in refamino.iteritems() :
                if s[:-1].find('X') == -1 :
                    fout.write('>{0}\n{1}\n'.format(n, s))
        return ref_aa

    def write_query(self, query) :
        fna, faa = '{0}.query.na'.format(parameters['unique_key']), '{0}.query.aa'.format(parameters['unique_key'])
        qryseq = self.readSequence(query)
        qryamino = transeq({n:s[0] for n,s in qryseq.iteritems()}, frame=7)
        with open(fna, 'w') as fout:
            for n, s in qryseq.iteritems() :
                fout.write('>{0}\n{1}\n'.format(n, s[0]))
        with open(faa, 'w') as fout :
            for n, s in qryamino.iteritems() :
                fout.write('>{0}\n{1}\n'.format(n, s))
        return qryseq, fna, faa

def nomenclature(genome, refAllele, parameters) :
    # write query
    dirPath = tempfile.mkdtemp(prefix='NS_', dir='.')
    try :
        qry = os.path.join(dirPath, 'query')
        ref = os.path.join(dirPath, 'reference')
        with open(qry, 'w') as fout :
            fout.write(genome)
        with open(ref, 'w') as fout :
            fout.write(refAllele)
        blasttab = dualBlast().run_comparison(dirPath, qry, ref, parameters['min_iden']-0.1, parameters['min_frag_len']-10)
    
        # filter
        blasttab_parser = blastParser()
        blasttab = blasttab_parser.linear_merge(blasttab, **parameters)
        loci = blasttab_parser.parse_blast(blasttab, parameters)
        regions = blasttab_parser.inter_loci_overlap(loci, parameters)
        #regions = blasttab_parser.intergenic(regions, parameters.get('intergenic',[30,600]))
    
        # submission
        qrySeq, qryQual = dualBlast().readFastq(qry)
        alleles = blasttab_parser.form_alleles(regions, qrySeq, qryQual, parameters['unique_key'], not parameters['query_only'], parameters)
        for field, allele in alleles.items() :
            allele['id'] = allele['value_md5'] if allele['accepted'] < 8 else '-'+allele['value_md5']
        if parameters.get('database', None) is not None :
            import pandas as pd
            database = pd.read_csv(parameters['database'], header=None, index_col=0)
            for field, allele in alleles.items() :
                try :
                    allele_info = database.loc[allele['value_md5']].values
                    if len(allele_info.shape) > 1 :
                        allele['id'] = str(allele_info[allele_info.T[0] == field][0, 1])
                    elif len(allele_info.shape) == 1 and allele_info[0] == field :
                        allele['id'] = str(allele_info[1])
                except :
                    pass
    finally:
        shutil.rmtree(dirPath)
    allele_seq = '\n'.join([ '>{0} value_md5={2} id={7} CIGAR={6} accepted={3} reference={4} identity={5} coordinates={8}\n{1}'.format(locus, allele['seq'], allele['value_md5'], allele['accepted'], allele['reference'], allele['identity'], allele['CIGAR'], allele['id'], '{0}:{1}..{2}:{3}'.format(*allele['coordinates']), allele['status'], allele['identity']) for locus, allele in sorted(alleles.items()) ])
    return allele_seq


parameters = {}
def MLSType(args) :
    parameters = getParams(args)
    alleles = nomenclature(open(parameters['genome']).read(), open(parameters['refAllele']).read(), parameters)
    if parameters['output'] is not None :
        fout = open(parameters['output'], 'w') if parameters['output'].upper() != 'STDOUT' else sys.stdout
        fout.write(alleles + '\n')
        fout.close()
    return alleles


def getParams(args) :
    import argparse
    parser = argparse.ArgumentParser(description='MLSType. Find and designate MLST alleles from a queried assembly. ')
    parser.add_argument('-i', '--genome',     help='[REQUIRED] Input - filename for genomic assembly. ', required=True)
    parser.add_argument('-r', '--refAllele',  help='[REQUIRED] Input - fasta file for reference alleles. ', required=True)
    parser.add_argument('-k', '--unique_key', help='[REQUIRED] An unique identifier for the assembly. ', required=True)
    parser.add_argument('-d', '--database',   help='[OPTIONAL] Input - lookup table of existing alleles. ', default=None)
    parser.add_argument('-o', '--output',     help='[DEFAULT: No output] Output - filename for the generated alleles. Specify to STDOUT for screen output. ', default=None)
    parser.add_argument('-q', '--query_only', help='[DEFAULT: False] Do not submit new allele, only query. ', action='store_true', default=False)
    parser.add_argument('-f', '--force', help='[DEFAULT: False] Force to accept low quality alleles. ', action='store_true', default=False)
    parser.add_argument('-m', '--min_iden',   help='[DEFAULT: 0.65 ] Minimum identities between refAllele and genome. ', type=float, default=0.65)
    parser.add_argument('-p', '--min_frag_prop', help='[DEFAULT: 0.6 ] Minimum covereage of a fragment. ', type=float, default=0.6)
    parser.add_argument('-l', '--min_frag_len',  help='[DEFAULT: 50 ] Minimum length of a fragment. ', type=float, default=50)
    
    parser.add_argument('-x', '--intergenic',  help='[DEFAULT: 50,500 ] Call alleles in intergenic region if to closely located loci differ by a distance between two numbers. ', default='50,500')

    parser.add_argument('--merging_prop',  help='[DEFAULT: 0.5 ] Two hits are conflicted if they cover by this proportion. ', type=float, default=0.5)
    parser.add_argument('--merging_error',  help='[DEFAULT: 0.05 ] Remove a secondary hit if its similarity is lower than another overlapped region by this value. ', type=float, default=0.05)
    
    parser.add_argument('--max_dist',  help='[DEFAULT: 300 ] Synteny block: Ignore if two alignments seperate by at least this value. ', type=float, default=300)
    parser.add_argument('--diag_diff',  help='[DEFAULT: 1.2 ] Synteny block: Ignore if the lengths of the resulted block differ by X fold between qry and ref. ', type=float, default=1.2)
    parser.add_argument('--max_diff',  help='[DEFAULT: 200 ] Synteny block: Ignore if the lengths of the resulted block differ by this value between qry and ref. ', type=float, default=200)
    args = parser.parse_args(args).__dict__
    args['intergenic'] = [float(d) for d in args['intergenic'].split(',')]
    
    return args
    

if __name__ == '__main__' :
    alleles = MLSType(sys.argv[1:])
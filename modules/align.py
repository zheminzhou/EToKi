# align multiple genomes onto a single reference, using minimap2
# remove short repetitive regions
# call SNPs and short indels
import os, sys, numpy as np, argparse, subprocess, re, gzip
from multiprocessing import Pool
try :
    xrange(1)
except :
    xrange = range

def parseArgs(argv) :
    parser = argparse.ArgumentParser(description='''Align multiple genomes onto a single reference. ''')
    parser.add_argument('-r', '--reference', help='[REQUIRED; INPUT] reference genomes to be aligned against. Use <Tag>:<Filename> format to assign a tag to the reference.', required=True)
    parser.add_argument('-p', '--prefix', help='[OUTPUT] prefix for all outputs.', default='Enlign')
    parser.add_argument('-a', '--alignment', help='[OUTPUT] Generate core genomic alignments in FASTA format', default=False, action='store_true')
    parser.add_argument('-m', '--matrix', help='[OUTPUT] Do not generate core SNP matrix', default=True, action='store_false')
    parser.add_argument('-c', '--core', help='[PARAM] percentage of presences for core genome. [DEFAULT: 0.95]', type=float, default=0.95)
    parser.add_argument('-n', '--n_proc', help='[PARAM] number of processes to use. [DEFAULT: 5]', default=5, type=int)
    parser.add_argument('queries', metavar='queries', nargs='+', help='queried genomes. Use <Tag>:<Filename> format to feed in a tag for each genome. Otherwise filenames will be used as tags for genomes. ')
    args = parser.parse_args(argv)
    args.reference = args.reference.split(':', 1) if args.reference.find(':')>0 else [os.path.basename(args.reference), args.reference]
    args.queries = sorted([ [qt, qf] for qt, qf in [ qry.split(':', 1) if qry.find(':')>0 else [os.path.basename(qry), qry] for qry in args.queries ] if qt != args.reference[0] ])
    return args

def alignAgainst(data) :
    prefix, minimap2, db, (rtag, reference), (tag, query) = data
    try :
        qrySeq, qryQual = readFastq(query)
    except :
        return [tag, query]
    refSeq, refQual = readFastq(reference)
    proc = subprocess.Popen('{0} -c -t1 --frag=yes -A2 -B8 -O20,40 -E3,2 -r20 -g200 -p.000001 -N5000 -f1000,5000 -n2 -m30 -s30 -z200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(
                                minimap2, db, query).split(), stdout=subprocess.PIPE, universal_newlines=True)
    alignments = []
    for lineId, line in enumerate(proc.stdout) :
        part = line.strip().split('\t')
        part[1:4] = [int(p) for p in part[1:4]]
        part[6:11] = [int(p) for p in part[6:11]]
        part[11] = float(part[13][5:])
        part[12], part[13] = lineId, part[11]/part[10]
        part[14:17] = [[], [], []]
        alignments.append(part)
    proc.wait()
    
    deleteChain = {}
    nItem = len(alignments)
    
    alignments.sort(key=lambda x:x[:4])
    for i1, p1 in enumerate(alignments) :
        for i2 in xrange(i1+1, nItem) :
            p2 = alignments[i2]
            if p1[0] != p2[0] : break

            s, e = max(p1[2], p2[2]), min(p1[3], p2[3])
            if s > e+10 :
                break
            if (e-s) >= 0.9 * (p1[3]-p1[2]) and p2[13] - 0.1 >= p1[13] :
                deleteChain[p1[12]] = deleteChain.get(p1[12], set([])) | set([p2[12]])
            if (e-s) >= 0.9 * (p2[3]-p2[2]) and p1[13] - 0.1 >= p2[13] :
                deleteChain[p2[12]] = deleteChain.get(p2[12], set([])) | set([p1[12]])
    alignments.sort(key=lambda x:x[5:9])
    for i1, p1 in enumerate(alignments) :
        for i2 in xrange(i1+1, nItem) :
            p2 = alignments[i2]
            if p1[5] != p2[5] : break

            s, e = max(p1[7], p2[7]), min(p1[8], p2[8])
            if s > e+10 :
                break
            
            if (e-s) >= 0.9 * (p1[8]-p1[7]) and p2[13] - 0.05 >= p1[13] :
                deleteChain[p1[12]] = deleteChain.get(p1[12], set([])) | set([p2[12]])
            if (e-s) >= 0.9 * (p2[8]-p2[7]) and p1[13] - 0.05 >= p2[13] :
                deleteChain[p2[12]] = deleteChain.get(p2[12], set([])) | set([p1[12]])

    deleted = {}
    for p in sorted(alignments, key=lambda x:x[11], reverse=True) :
        id = p[12]
        if id in deleteChain :
            for jd in deleteChain[id] :
                if jd not in deleted :
                    deleted[id] = 1
                    break
    alignments = [p for p in alignments if p[12] not in deleted]
    
    # repeats in qry
    nItem = len(alignments)
    alignments.sort(key=lambda x:x[:4])
    for i1, p1 in enumerate(alignments) :
        for i2 in xrange(i1+1, nItem) :
            p2 = alignments[i2]
            if p1[0] != p2[0] : break
            s, e = max(p1[2], p2[2]), min(p1[3], p2[3])
            if e > s :
                p1[16].append([s, e])
                p2[16].append([s, e])
            else :
                break
    # repeats in ref
    alignments.sort(key=lambda x:x[5:9])
    for i1, p1 in enumerate(alignments) :
        for i2 in xrange(i1+1, nItem) :
            p2 = alignments[i2]
            if p1[5] != p2[5] : break
            s, e = max(p1[7], p2[7]), min(p1[8], p2[8])
            if e > s :
                p1[15].append([s, e])
                p2[15].append([s, e])
            else :
                break
    
    maskedRegion = {}
    refRepeat = []
    for p in alignments :
        # prepare a unique set of repeat region
        qryRepeat = []
        if len(p[16]) > 0 :
            qryRepeat.append(p[16][0])
            for pp in p[16][1:] :
                if pp[0] > qryRepeat[-1][1]+20 :
                    qryRepeat.append(pp)
                elif pp[1] > qryRepeat[-1][1]:
                    qryRepeat[-1][1] = pp[1]
        ref = [refSeq[p[5]], refQual[p[5]]]
        qry = [qrySeq[p[0]], qryQual[p[0]]]
        cigar = p[-1][5:]
        d = 1 if p[4] == '+' else -1
        if d < 0 :
            qryRepeat = [[q[1], q[0], -1, -1] for q in qryRepeat]
        else :
            qryRepeat = [[q[0], q[1], -1, -1] for q in reversed(qryRepeat)]

        mut = []
        alnSite = [p[7], p[2] if d > 0 else p[3]-1]
        for cl, ct in re.findall(r'(\d+)([MID])', cigar) :
            cl = int(cl)
            if ct == 'M' :
                # extract aligned sequences
                r = ref[0][alnSite[0]:alnSite[0]+cl]
                r1 = ref[1][alnSite[0]:alnSite[0]+cl]
                q = qry[0][alnSite[1]:alnSite[1]+cl] if d > 0 else rc(qry[0][(alnSite[1]-cl+1):(alnSite[1]+1)])
                q1 = qry[1][alnSite[1]:alnSite[1]+cl] if d > 0 else ''.join(reversed(qry[0][(alnSite[1]-cl+1):(alnSite[1]+1)]))

                e =[alnSite[0]+cl, alnSite[1]+cl*d]
                for qid in xrange(len(qryRepeat)-1, -1, -1) :
                    qr = qryRepeat[qid]
                    if d*qr[0] <= d*e[1] :
                        if qr[2] == -1 :
                            qr[2] = alnSite[0] + d*(qr[0] - alnSite[1])
                        if d*qr[1] <= d*e[1] :
                            qr[3] = alnSite[0] + d*(qr[1] - alnSite[1])
                            p[15].append(qr[2:])
                            del qryRepeat[qid]
                    else :
                        break
                for id, (rr, rr1, qq, qq1) in enumerate(np.array([list(r), list(r1), list(q), list(q1)]).T) :
                    if ord(rr1) < 43 or ord(qq1) < 43 :
                        maskedRegion[(p[5], alnSite[0]+id)] = 0
                    if rr != qq and rr != 'N' and qq != 'N' :
                        mut.append([alnSite[0]+id, alnSite[1]+id*d, rr, qq, p[4]])
                alnSite = e
            elif ct == 'I' :
                q = qry[0][alnSite[1]:alnSite[1]+cl] if d < 0 else rc(qry[0][(alnSite[1]-cl+1):(alnSite[1]+1)] )
                q1 = qry[1][alnSite[1]:alnSite[1]+cl] if d > 0 else ''.join(reversed(qry[0][(alnSite[1]-cl+1):(alnSite[1]+1)] ))
                
                e = alnSite[1] + cl*d
                for qid in xrange(len(qryRepeat)-1, -1, -1) :
                    qr = qryRepeat[qid]
                    if d*qr[0] <= d*e :
                        if qr[2] == -1 :
                            qr[2] = alnSite[0]
                        if d*qr[1] <= d*e :
                            qr[3] = alnSite[0]
                            p[15].append(qr[2:])
                            del qryRepeat[qid]
                    else :
                        break
                
                if ord(min(list(q1))) >= 43 :
                    mut.append([alnSite[0], min(alnSite[1], e), '.', '+' + q, p[4]])
                for site in xrange(alnSite[0], alnSite[0]+2) :
                    maskedRegion[(p[5], site)] = 0
                alnSite[1] = e
            elif ct == 'D' :
                r = ref[0][alnSite[0]:alnSite[0]+cl]
                r1 = ref[1][alnSite[0]:alnSite[0]+cl]
                if ord(min(list(r1))) >= 43 :
                    mut.append([alnSite[0], int(alnSite[1]+0.5*d), '.', '-' + r, p[4]])
                for site in xrange(alnSite[0], alnSite[0]+2) :
                    maskedRegion[(p[5], site)] = 0
                alnSite[0]+=cl
        p[14] = mut
        refRepeat.extend([ [p[5], pp[0], pp[1]] for pp in p[15] ])

    repeats = []
    if len(refRepeat) :
        refRepeat.sort()
        repeats = [refRepeat[0]]
        for p in refRepeat[1:] :
            if p[0] != repeats[-1][0] or p[1] - 20 > repeats[-1][2] :
                repeats.append(p)
            elif p[2] > repeats[-1][2] :
                repeats[-1][2] = p[2]

    for p in repeats :
        for site in xrange(p[1], p[2]) :
            maskedRegion[(p[0], site)] = 1

    repeats = []
    for cont, site in sorted(maskedRegion) :
        if len(repeats) == 0 or repeats[-1][0] != cont or repeats[-1][2]+1 < site :
            repeats.append([cont, site, site])
        else :
            repeats[-1][2] = site
  
    mutations = []
    alignments = [aln for aln in alignments if aln[9] >= 100]
    for aln in alignments :
        for m in aln[14] :
            if len(m[3]) == 1 :
                if (aln[5], m[0]) not in maskedRegion :
                    mutations.append([aln[5], aln[0]] + m)
            elif maskedRegion.get((aln[5], m[0]), 0) != 1 :
                if m[3].startswith('-') and maskedRegion.get((aln[5], m[0]+len(m[3])-2), 0) > 0 :
                    continue
                mutations.append([aln[5], aln[0]] + m)
    with uopen(prefix + '.gff.gz', 'w') as fout :
        fout.write('##gff-version 3\n')
        fout.write('## Reference: {0}\n'.format(reference))
        fout.write('## Query: {0}\n'.format(query))
        fout.write('## Tag: {0}\n'.format(tag))
        for aln in alignments :
            if aln[5] == aln[0] and aln[2] == aln[7] and aln[3] == aln[8] :
                fout.write('{0}\trefMapper\tmisc_feature\t{1}\t{2}\t{3}\t{4}\t.\t/inference="Self%20Alignments"\n'.format(
                    aln[5], aln[7]+1, aln[8], aln[9], aln[4], aln[0], aln[2]+1, aln[3], 
                ))
            else :
                fout.write('{0}\trefMapper\tmisc_feature\t{1}\t{2}\t{3}\t{4}\t.\t/inference="Aligned%20with%20{5}:{6}-{7}"\n'.format(
                    aln[5], aln[7]+1, aln[8], aln[9], aln[4], aln[0], aln[2]+1, aln[3], 
                ))
                
        for p in repeats :
            fout.write('{0}\trefMapper\tunsure\t{1}\t{2}\t.\t+\t.\t/inference="Uncertain%20base%20calling%20or%20ambigious%20alignment"\n'.format(
                p[0], p[1]+1, p[2]+1, 
            ))
        for mut in mutations :
            e1 = mut[2] if not mut[5].startswith('-') else mut[2] + len(mut[5]) - 2
            e2 = mut[3] if not mut[5].startswith('+') else mut[3] + len(mut[5]) - 2
            if len(mut[5]) > 26 :
                mut[5] = '{0}[{1}bps]'.format(mut[5][0], len(mut[5])-1)

            fout.write('{0}\trefMapper\tvariation\t{1}\t{2}\t.\t+\t.\t/replace="{7}";/compare="{3}:{4}-{5}:{8}";/origin="{6}"\n'.format(
                mut[0], mut[2]+1, e1+1, mut[1], mut[3]+1, e2+1, mut[4], mut[5], mut[6]
            ))

    return [tag, prefix + '.gff.gz']

def readMap(data) :
    mTag, mFile = data
    presences, absences, mutations = [], [], []
    
    aligns = ['', 0, 0]
    miss = -1
    with uopen(mFile) as fin :
        for line in fin :
            if line.startswith('##') :
                if line.startswith('## Reference: ') :
                    ref = line.split(' ')[-1]
                elif line.startswith('## Query: ') :
                    qry = line.split(' ')[-1]
                    if ref == qry :
                        miss = -99999999
            else :
                break
    print(mTag, mFile)
    with uopen(mFile) as fin :
        for line in fin :
            if line.startswith('#') : continue
            part = line.strip().split('\t')
            part[3:5] = [int(part[3]), int(part[4])]
            if part[2] == 'misc_feature' :
                if len(presences) == 0 or presences[-1][0] != part[0] or presences[-1][2] < part[3] :
                    presences.append([part[0], part[3], part[4]])
                elif presences[-1][2] < part[4] :
                    presences[-1][2] = part[4]
            elif part[2] == 'unsure' :
                absences.append([part[0], part[3], part[4], miss])
            elif part[2] == 'variation' :
                alt = re.findall(r'replace="([^"]+)"', part[8])
                ori = re.findall(r'origin="([^"]+)"', part[8])
                if len(alt) and len(ori) :
                    mutations.append([mTag, part[0], part[3], ori[0], alt[0]])
    return presences, absences, mutations

def getMatrix(prefix, reference, alignments, core, matrixOut, alignmentOut) :
    refSeq, refQual = readFastq(reference)
    coreSites = { n:np.zeros(len(refSeq[n]), dtype=int) for n in refSeq }
    matSites = { n:np.zeros(len(refSeq[n]), dtype=int) for n in refSeq }
    alnId = { aln[0]:id for id, aln in enumerate(alignments) }
    res = pool.map(readMap, alignments)
    
    matrix = {}
    for presences, absences, mutations in res :
        for mut in mutations :
            j = alnId[mut[0]]
            site = tuple(mut[1:3])
            if site not in matrix :
                matrix[site] = [[], []]
                matSites[mut[1]][mut[2]-1] = mut[2]
            if len(mut[4]) == 1 :
                if len(matrix[site][0]) == 0 :
                    matrix[site][0] = ['-' for id in alnId]
                matrix[site][0][j] = mut[4]
            else :
                if len(matrix[site][1]) == 0 :
                    matrix[site][1] = ['-' for id in alnId]
                matrix[site][1][j] = mut[4]
    for (mTag, mFile), (presences, absences, mutations) in zip(alignments, res) :
        j = alnId[mTag]
        for n, s, e in presences :
            coreSites[n][s-1:e] +=1
            mutations = matSites[n][s-1:e]
            for kk in mutations[mutations > 0] :
                k = (n, kk)
                if len(matrix[k][0]) and matrix[k][0][j] == '-' :
                    matrix[k][0][j] = '.'
                if len(matrix[k][1]) and matrix[k][1][j] == '-' :
                    matrix[k][1][j] = '.'
        for n, s, e, m in absences :
            coreSites[n][s-1:e] -=1
            mutations = matSites[n][s-1:e]
            for kk in mutations[mutations > 0] :
                k = (n, kk)
                if len(matrix[k][0]) and matrix[k][0][j] == '.' :
                    matrix[k][0][j] = '-'
                if len(matrix[k][1]) and matrix[k][1][j] == '.' :
                    matrix[k][1][j] = '-'
    pres = np.unique(np.concatenate(list(coreSites.values())), return_counts=True)
    pres = [pres[0][pres[0] > 0], pres[1][pres[0] > 0]]
    coreNum = len(alignments) * core
    for p, n in zip(*pres) :
        sys.stderr.write('#{2} {0} {1}\n'.format(p, n, '' if p > coreNum else '#'))

    missings = []
    coreBases = {'A':0, 'C':0, 'G':0, 'T':0}
    for n in sorted(coreSites) :
        sites = coreSites[n]
        for site, num in enumerate(sites) :
            cSite = (n, site+1)
            if num < coreNum and cSite in matrix and len(matrix[cSite][1]) > 0 :
                num = np.sum(matrix[cSite][1] != '-')
                matrix[cSite][0] = []
            if num < coreNum :
                matrix.pop(cSite, None)
                if len(missings) == 0 or missings[-1][0] != n or missings[-1][2] + 1 < cSite[1] :
                    missings.append([n, cSite[1], cSite[1]])
                else :
                    missings[-1][2] = cSite[1]
            else :
                b = refSeq[n][cSite[1]-1]
                if cSite in matrix and len(matrix[cSite][0]) :
                    matrix[cSite][0] = [ (b if s == '.' else s) for s in matrix[cSite][0]]
                else :
                    coreBases[b] = coreBases.get(b, 0) + 1
                    
    outputs = {}
    if matrixOut :
        outputs['matrix'] = prefix + '.matrix.gz'
        with uopen(prefix + '.matrix.gz', 'w') as fout :
            fout.write('## Constant_bases: {A} {C} {G} {T}\n'.format(**coreBases))
            for n in refSeq :
                fout.write('## Sequence_length: {0} {1}\n'.format(n, len(refSeq[n])))
            for region in missings :
                fout.write('## Missing_region: {0} {1} {2}\n'.format(*region))
            fout.write('\t'.join(['#Seq', '#Site'] + [ mTag for mTag, mFile in alignments ]) + '\n')
            for site in sorted(matrix) :
                bases = matrix[site]
                if len(bases[0]) :
                    fout.write('{0}\t{1}\t{2}\n'.format(site[0], site[1], '\t'.join(bases[0])))
                if len(bases[1]) :
                    fout.write('{0}\t{1}\t{2}\n'.format(site[0], site[1], '\t'.join(bases[1])))
    if alignmentOut :
        outputs['alignment'] = prefix + '.fasta.gz'
        sequences = []
        for (mTag, mFile), (presences, absences, mutations) in zip(alignments, res) :
            j = alnId[mTag]
            seq = { n:['-']*len(s) for n, s in refSeq.items() } if j > 0 else { n:list(s) for n, s in refSeq.items() }
            if j :
                for n, s, e in presences :
                    seq[n][s-1:e] = refSeq[n][s-1:e]
                for n, s, e, c in absences :
                    seq[n][s-1:e] = '-' * (e-s+1)
            for site in matrix :
                bases = matrix[site]
                if len(bases[0]) :
                    seq[site[0]][site[1]-1] = bases[0][j]
            sequences.append(seq)
        with uopen(prefix + '.fasta.gz', 'w') as fout :
            for id, n in enumerate(sorted(refSeq)) :
                if id :
                    fout.write('=\n')
                for (mTag, mFile), seq in zip(alignments, sequences) :
                    fout.write('>{0}:{1}\n{2}\n'.format(mTag, n, ''.join(seq[n])))
    return outputs

def runAlignment(prefix, reference, queries, core, minimap2) :
    #alignments = list(map(alignAgainst, [[prefix +'.' + query[0].rsplit('.', 1)[0] + '.' + str(id), minimap2, prefix + '.mmi', reference, query] for id, query in enumerate(queries)]))
    alignments = pool.map(alignAgainst, [[prefix +'.' + query[0].rsplit('.', 1)[0] + '.' + str(id), minimap2, prefix + '.mmi', reference, query] for id, query in enumerate(queries)])

    try :
        os.unlink(reference + '.mmi')
    except :
        pass
    return alignments

def prepReference(prefix, reference, minimap2, pilercr, trf, **args) :
    # prepare minimap2 reference
    if reference :
        subprocess.Popen('{0} -k15 -w5 -d {2}.mmi {1}'.format(minimap2, reference, prefix).split(), stderr=subprocess.PIPE).communicate()
    sys.stdout.write(prefix)
    return

def align(argv) :
    args = parseArgs(argv)
    
    global pool
    pool = Pool(args.n_proc)
    #print(args.reference)
    refMask = prepReference(args.prefix, args.reference[1], **externals)
    alignments = runAlignment(args.prefix, args.reference, [args.reference] + args.queries, args.core, externals['minimap2'])
    outputs = {'mappings': dict(alignments)}
    if args.matrix or args.alignment :
        outputs.update(getMatrix(args.prefix, args.reference[1], alignments, args.core, args.matrix, args.alignment))
    import json
    sys.stdout.write(json.dumps(outputs, indent=2, sort_keys=True))
    return outputs

pool = None
if __name__ == '__main__' :
    from configure import externals, readFastq, rc, uopen
    from uberBlast import uberBlast
    align(sys.argv[1:])
else :
    from .configure import externals, readFastq, rc, uopen
    from .uberBlast import uberBlast
    

# align multiple genomes onto a single reference, using minimap2
# remove short repetitive regions
# call SNPs and short indels
import os, sys, numpy as np, argparse, subprocess, re, gzip, _collections
from multiprocessing import Pool
try :
    from .configure import readFastq, readFasta, xrange
except :
    from configure import readFastq, readFasta, xrange

def parseArgs(argv) :
    parser = argparse.ArgumentParser(description='''Align multiple genomes onto a single reference. ''')
    parser.add_argument('-r', '--reference', help='[REQUIRED; INPUT] reference genomes to be aligned against. Use <Tag>:<Filename> format to assign a tag to the reference.', required=True)
    parser.add_argument('-g', '--gff', help='[INPUT] gff for reference genomes. Used only to annotate SNPs', default=None)
    parser.add_argument('-G', '--gcode', help='[INPUT] translation table for coding sequences. default:1', default=1)
    parser.add_argument('--qry_list', help='[INPUT] name of a list file that contains genome queries (one per line).', default=None)
    parser.add_argument('-p', '--prefix', help='[OUTPUT] prefix for all outputs.', default='Enlign')
    parser.add_argument('-a', '--alignment', help='[OUTPUT] Generate core genomic alignments in FASTA format', default=False, action='store_true')
    parser.add_argument('-m', '--matrix', help='[OUTPUT] Do not generate core SNP matrix', default=True, action='store_false')
    parser.add_argument('-l', '--last', help='Activate to use LAST as aligner. [DEFAULT: minimap2]', default=False, action='store_true')
    parser.add_argument('-c', '--core', help='[PARAM] percentage of presences for core genome. [DEFAULT: 0.95]', type=float, default=0.95)
    parser.add_argument('-n', '--n_proc', help='[PARAM] number of processes to use. [DEFAULT: 5]', default=5, type=int)
    parser.add_argument('-q', '--lowq', help='[OPTIONAL; INPUT] Genome of low quality. [DEFAULT: ]', default=[], action='append')
    parser.add_argument('-d', '--divergent', help='[OPTIONAL] Alignment for divergent sequences. [DEFAULT: False]', default=False, action='store_true')
    parser.add_argument('queries', metavar='queries', nargs='*', help='queried genomes. Use <Tag>:<Filename> format to feed in a tag for each genome. Otherwise filenames will be used as tags for genomes. ')
    args = parser.parse_args(argv)
    args.reference = args.reference.split(':', 1) if args.reference.find(':')>0 else [os.path.basename(args.reference), args.reference]
    args.lowq = sorted([ [qt, qf] for qt, qf in [ qry.split(':', 1) if qry.find(':')>0 else [os.path.basename(qry), qry] for qry in args.lowq ] if qt != args.reference[0] ])
    args.aligner = externals['minimap2'] if not args.last else [externals['lastdb'], externals['lastal']]
    if args.divergent :
        global divergent
        divergent = True
    return args

class last_package(object) :
    @staticmethod
    def run_lastal( refdb, query, output, lastal ) :
        cmd = '{0} -j4 -r1 -q2 -a7 -b1 {1} {2}'.format( lastal, refdb, query )
        lastal_run = subprocess.Popen( cmd.split(), stdout=subprocess.PIPE, universal_newlines=True )
        with open(output, 'w') as fout:
            fout.write(lastal_run.communicate()[0])
        if lastal_run.returncode != 0 :
            fastq = readFastq(query)
            with open(output+'.qry', 'w') as fout :
                for n, (s, q) in fastq.items() :
                    fout.write('@{0}\n{1}\n+\n{2}\n'.format(n, s, re.sub(r'[!"#$%&\']', '(', q)))
            cmd = '{0} -Q1 -j4 -r1 -q2 -a7 -b1 {1} {2}'.format( lastal, refdb, output + '.qry' )
            lastal_run = subprocess.Popen( cmd.split(), stdout=subprocess.PIPE )
            with open(output, 'w') as fout:
                fout.write(lastal_run.communicate()[0])
            os.unlink(output + '.qry')
        return output
    @staticmethod
    def call_mutation( comparison, rep_mask = 3 ) :
        seq = [comparison[6], comparison[12]]
        qual = comparison[13]
        direct = [1 if comparison[4] == '+' else -1, 1 if comparison[10] == '+' else -1]
        coord = [comparison[2] *direct[0] -1, comparison[8]*direct[1] -1]
        low_complexity = [[coord[0]-999, coord[0]]]
        low_qual = []
        mutations = []
        for id, s1 in enumerate(seq[0]) :
            s2 = seq[1][id]
            if s1 != '-' :
                coord[0] += 1
            if s2 != '-' :
                coord[1] += 1
            if qual[id] == 1 :
                if len(low_qual) == 0 or low_qual[-1][3] + 3 < coord[1] :
                    low_qual.append([coord[0], coord[0], coord[1], coord[1], 0])
                else :
                    low_qual[-1][1] = coord[0]
                    low_qual[-1][3] = coord[1]
            if s1.islower() or s2.islower() :
                if len(low_complexity) == 0 or low_complexity[-1][1] < coord[0]-1 :
                    low_complexity.append([coord[0], coord[0]])
                else :
                    low_complexity[-1][1] = coord[0]
            if s1.upper() != s2.upper() : #and s1.upper() != 'N' and s2.upper() != 'N' :
                ms1 = coord[0] 
                ms2 = coord[1] 
                if len(mutations) < 1 or ((ms1 != mutations[-1][0] or mutations[-1][4][-1] != '-') and (ms2 != mutations[-1][2] or mutations[-1][5][-1] != '-')) :
                    mutations.append([ms1, ms1, ms2, ms2, s1.upper(), s2.upper()])
                    if id == 0 :
                        mutations[-1].append( max(qual[id:(id+2)]) )
                    else :
                        mutations[-1].append( max(qual[(id-1):(id+2)]) )
                else :
                    mutations[-1][1] = ms1
                    mutations[-1][3] = ms2
                    mutations[-1][4] += s1.upper()
                    mutations[-1][5] += s2.upper()
                    if id < len(seq[0]) -1 and mutations[-1][6] == 0 :
                        mutations[-1][6] = qual[id+1]
        low_complexity = [ reg for reg in low_complexity + [[coord[0]+1, coord[0]+1000]] if reg[1]-reg[0] +1 >= 50 ]
        for lq in low_qual :
            if lq[0] > comparison[2] :
                lq[0] -= 1
            if lq[1] < comparison[3] :
                lq[1] += 1
        comparison[6] = low_qual + low_complexity[1:-1]
        comparison[12] = ''
        comparison = comparison[:13]
        if len(mutations) > 0 :
            for mut in mutations :
                if mut[6] == 0 :
                    for reg in low_complexity :
                        if mut[1] < reg[0]-rep_mask :
                            break 
                        elif mut[0] <= reg[1]+rep_mask :
                            mut[6] = 1
            comparison.extend(mutations)
        ref_ins = [len(mut[5]) for mut in mutations if mut[5][0] == '-']
        qry_ins = [len(mut[5]) for mut in mutations if mut[4][0] == '-']
        gap_open = len(ref_ins) + len(qry_ins)
        mut_num = len(mutations) - gap_open
        comparison[0] = (comparison[3]-comparison[2]+1) - 3 * mut_num - 7*gap_open - sum(qry_ins) - 2*sum(ref_ins) 
    
        return comparison

    # 0-6    score, ref_cont, ref_start, ref_end, ref_direct, ref_len, ref_seq
    # 7-12  qry_cont, qry_start, qry_end, qry_direct, qry_len, qry_seq
    # 13+   (mutations ...) 
    @staticmethod
    def sub_comparison(comparison, ref_coords=None, qry_coords=None) :
        if (ref_coords is None) == (qry_coords is None) :
            raise Exception('Exact only one set of coords, either qry or ref')
        direct = [1 if comparison[4] == '+' else -1, 1 if comparison[10] == '+' else -1]
        mutations = []
        if ref_coords is not None :
            if direct[0] > 0 :
                rc = [max(ref_coords[0], comparison[2]), min(ref_coords[1], comparison[3])]
            else :
                rc = [max(-ref_coords[1], -comparison[2]), min(-ref_coords[0], -comparison[3])]
            if rc[0] > rc[1] :
                return []
            qc = [0, 0]
            mutations = [ mut[:] for mut in comparison[13:] if mut[1]>= rc[0] and mut[0] <= rc[1] ]
            if len(mutations) > 0 :
                if mutations[0][0] <= rc[0] and mutations[0][5][0] == '-' :
                    qc[0] = mutations[0][2] + 1
                    mutations[0][4] = mutations[0][4][(rc[0]-mutations[0][0]):]
                    mutations[0][5] = mutations[0][5][(rc[0]-mutations[0][0]):]
                    mutations[0][0] = rc[0]
                elif mutations[0][4][0] == '-' :
                    qc[0] = mutations[0][2] - (mutations[0][0]-rc[0]+1)
                elif mutations[0][5][0] == '-' :
                    qc[0] = mutations[0][2] - (mutations[0][0]-rc[0]-1)
                else :
                    qc[0] = mutations[0][2] - (mutations[0][0]-rc[0])
                if mutations[-1][0] == rc[1] and mutations[-1][4][0] == '-' :
                    qc[1] = mutations[-1][2]-1
                    mutations.pop(-1)
                else :
                    qc[1] = mutations[-1][3] + (rc[1] - mutations[-1][1])
                    if mutations[-1][1] >= rc[1] and mutations[-1][5][0] == '-' :
                        mutations[-1][1] = rc[1]
                        mutations[-1][4] = mutations[-1][4][:(rc[1]-mutations[-1][0]+1)]
                        mutations[-1][5] = mutations[-1][5][:(rc[1]-mutations[-1][0]+1)]
            else :
                pre_mut = [ mut for mut in [[comparison[2]*direct[0], comparison[2], comparison[8]*direct[1], comparison[8]]] + comparison[13:] if mut[0] <= rc[0] ][-1]
                delta = -pre_mut[0] + pre_mut[2]
                qc = [rc[0]+delta, rc[1] + delta ]
        else :
            if direct[1] > 0 :
                qc = [max(qry_coords[0], comparison[8]), min(qry_coords[1], comparison[9])]
            else :
                qc = [max(-qry_coords[1], -comparison[8]), min(-qry_coords[0], -comparison[9])]
            if qc[0] > qc[1] :
                return []
            rc = [0, 0]
            mutations = [ mut[:] for mut in comparison[13:] if mut[3]>= qc[0] and mut[2] <= qc[1] ]
            if len(mutations) > 0 :
                if mutations[0][2] <= qc[0] and mutations[0][4][0] == '-' :
                    rc[0] = mutations[0][0] + 1
                    mutations[0][4] = mutations[0][4][(qc[0]-mutations[0][2]):]
                    mutations[0][5] = mutations[0][5][(qc[0]-mutations[0][2]):]
                    mutations[0][2] = qc[0]
                elif mutations[0][4][0] == '-' :
                    rc[0] = mutations[0][0] - (mutations[0][2]-qc[0]-1)
                elif mutations[0][5][0] == '-' :
                    rc[0] = mutations[0][0] - (mutations[0][2]-qc[0]+1)
                else :
                    rc[0] = mutations[0][0] - (mutations[0][2]-qc[0])
                if mutations[-1][2] >= qc[1] and mutations[-1][5][0] == '-' :
                    rc[1] = mutations[-1][0]-1
                    mutations.pop(-1)
                else :
                    rc[1] = mutations[-1][1] + (qc[1]-mutations[-1][3])
                    if mutations[-1][3] >= qc[1] and mutations[-1][4][0] == '-' :
                        mutations[-1][3] = qc[1]
                        mutations[-1][4] = mutations[-1][4][:(qc[1]-mutations[-1][2]+1)]
                        mutations[-1][5] = mutations[-1][5][:(qc[1]-mutations[-1][2]+1)] 
            else :
                pre_mut = [ mut for mut in [[comparison[2]*direct[0], comparison[2], comparison[8]*direct[1], comparison[8]]] + comparison[13:] if mut[2] <= qc[0] ][-1]
                delta = -pre_mut[0] + pre_mut[2]
                rc = [qc[0]-delta, qc[1] - delta ]
    
        ref_ins = [len(mut[5]) for mut in mutations if mut[5][0] == '-']
        qry_ins = [len(mut[5]) for mut in mutations if mut[4][0] == '-']
        gap_open = len(ref_ins) + len(qry_ins)
        mut_num = len(mutations) - gap_open
        score = (rc[1]-rc[0]+1) - 3 * mut_num - 7*gap_open - sum(qry_ins) - 2*sum(ref_ins) 
        return [score, comparison[1], abs(rc[0]), abs(rc[1]), comparison[4], comparison[5], comparison[6], comparison[7], abs(qc[0]), abs(qc[1]), comparison[10], comparison[11], comparison[12]] + mutations

    @staticmethod
    def make_alignment( filename ) :
        comparisons = []
        with open(filename, 'r') as fin:
            for line in fin:
                if line[0] == 'a' :
                    comparison = [ int(line.split(' ', 2)[1][6:]) ]
                elif line[0] == 's' :
                    part = line.strip().split()[1:]
                    part[1:5] = [int(part[1]), int(part[2]), part[3], int(part[4])]
                    if part[3] == '+' :
                        part[1:3] = [part[1]+1, part[1]+part[2]]
                    else :
                        part[1:3] = [part[4]-part[1], part[4]-part[1]-part[2]+1]
                    comparison.extend(part)
                    if len(comparison) >= 13 :
                        comparison.append([int((m in 'nN') or (n in 'Nn')) for m, n in zip(comparison[6], comparison[12])])
                elif line[0] in 'pq' :
                    part = line.strip().split()
                    comparison[13] = [ max(comparison[13][id], int(b in '!"#$%&\'()*+,-./')) for id, b in enumerate(part[-1])]
                elif len(line.strip()) == 0 :
                    if comparison[0] >= 200 : 
                        comparisons.append( last_package.call_mutation(comparison) )
        
        # remove significant low identity regions in query
        comparisons.sort(key=lambda x: min(x[8:10]) )
        comparisons.sort(key=lambda x: x[7] )
    
        low_q = []
        for id, regi in enumerate(comparisons) :
            if len(regi) == 0 : continue
            for jd in xrange(id+1, len(comparisons)) :
                regj = comparisons[jd]
                if len(regj) == 0 : continue
                if regi[7] != regj[7] : break
                si, ei = sorted(regi[8:10])
                sj, ej = sorted(regj[8:10])
                s = max(si, sj)
                e = min(ei, ej)
                if e >= s :
                    overlap_i = last_package.sub_comparison(regi, qry_coords=[s, e])
                    overlap_j = last_package.sub_comparison(regj, qry_coords=[s, e])
                    
                    if overlap_i[0] < 0.95 * overlap_j[0] and ( regi[0] < regj[0] or ei < ej ) : 
                        if s - si >= 30 :
                            comparisons[id] = last_package.sub_comparison(regi, qry_coords=[si, s-1])
                            #if overlap_i[3] < overlap_i[2] and overlap_i[4] == '+' :
                                #print comparisons[id]
                            overlap_i[12] = 'E'
                            low_q.append(overlap_i)
                            if overlap_i[3] >= overlap_i[2]  :
                                regi = comparisons[id]
                            if len(regi) == 0: break
                        else :
                            comparisons[id][12] = 'E'
                            break
                    elif overlap_i[0] * 0.95 > overlap_j[0] :
                        if ej - e >= 30 :
                            comparisons[jd] = last_package.sub_comparison(regj, qry_coords=[e+1, ej])
                            overlap_j[12] = 'E'
                            if overlap_j[3] >= overlap_j[2]  :
                                low_q.append(overlap_j)
                        else :
                            comparisons[jd][12] = 'E'
                    elif s == si and e == ei and regj[0] > regi[0]*3 and overlap_i[0] <= overlap_j[0] :
                        comparisons[id][12] = 'E'
                        break
                    elif s == sj and e == ej and regi[0] > regj[0]*3 and overlap_i[0] >= overlap_j[0] :
                        comparisons[jd][12] = 'E'
                    else :
                        comparisons[id][12] = 'D'
                        comparisons[jd][12] = 'D'
                else :
                    break
    
        # remove significant low identity regions in reference
        comparisons = sorted([x for x in comparisons if len(x) > 0] + low_q, key=lambda x: x[2] )
        comparisons.sort(key=lambda x: x[1] )
    
        for id, regi in enumerate(comparisons) :
            if len(regi) == 0 : continue
            for jd in xrange(id+1, len(comparisons)) :
                regj = comparisons[jd]
                if len(regj) == 0 : continue
                if regi[1] != regj[1] : break
                si, ei = regi[2:4]
                sj, ej = regj[2:4]
                s = max(si, sj)
                e = min(ei, ej)
                if e >= s :                
                    overlap_i = last_package.sub_comparison(regi, ref_coords=[s, e])
                    overlap_j = last_package.sub_comparison(regj, ref_coords=[s, e])
                    
                    if overlap_i[0] < 0.95 * overlap_j[0] and ( regi[0] < regj[0] or ei < ej ) : 
                        if s - si >= 30 :
                            comparisons[id] = last_package.sub_comparison(regi, ref_coords=[si, s-1])
                            regi = comparisons[id]
                            if len(regi) == 0: break
                        else :
                            comparisons[id] = []
                            break
                    elif overlap_i[0] * 0.95 > overlap_j[0] :
                        if ej - e >= 30 :
                            comparisons[jd] = last_package.sub_comparison(regj, ref_coords=[e+1, ej])
                        else :
                            comparisons[jd] = []
                    elif overlap_i[0] == overlap_j[0] and len(overlap_i) == len(overlap_j) :
                        if si == sj and ei == ej:
                            diff = 0
                            for i, i_snp in enumerate(overlap_i[13:]) :
                                j_snp = overlap_j[13+i]
                                if i_snp[0] != j_snp[0] or i_snp[5] != j_snp[5] :
                                    diff = 1
                                    break
                            if diff == 0 :
                                if comparisons[id][12] in 'DE':
                                    comparisons[id] = [] 
                                    break
                                else: 
                                    comparisons[jd] = []
                        elif si <= sj and ei >= ej :
                            diff = 0
                            for i, i_snp in enumerate(overlap_i[13:]) :
                                j_snp = overlap_j[13+i]
                                if i_snp[0] != j_snp[0] or i_snp[5] != j_snp[5] :
                                    diff = 1
                                    break
                            if diff == 0 :
                                comparisons[jd] = []
                        elif si >= sj and ei <= ej:
                            diff = 0
                            for i, i_snp in enumerate(overlap_i[13:]) :
                                j_snp = overlap_j[13+i]
                                if i_snp[0] != j_snp[0] or i_snp[5] != j_snp[5] :
                                    diff = 1
                                    break
                            if diff == 0 :
                                comparisons[id] = []
                                break
                else :
                    break
            if len(comparisons[id]) > 0 :
                regi = comparisons[id]
                regi[6] = [lq for lq in regi[6] if lq[1] >= regi[2] and lq[0] <= regi[3]]
        
        # mark repetitive regions in query
        repeats = []
        mutations = {}
        comparisons = sorted([x for x in comparisons if len(x) > 0 and x[12] != 'E'], key=lambda x: min(x[8:10]) )
        comparisons.sort(key=lambda x: x[7] )
        
        for id, regi in enumerate(comparisons) :
            for jd in xrange(id+1, len(comparisons)) :
                regj = comparisons[jd]
                if regi[7] != regj[7] : break
                si, ei = sorted(regi[8:10])
                sj, ej = sorted(regj[8:10])
                s = max(si, sj)
                e = min(ei, ej)
                if e >= s :
                    for mut in regi[13:] :
                        if abs(min(mut[2:4])) <= e and abs(max(mut[2:4])) >= s :
                            mut[6] = 1
                    for mut in regj[13:] :
                        if abs(min(mut[2:4])) <= e and abs(max(mut[2:4])) >= s :
                            mut[6] = 1
                    overlap_i = last_package.sub_comparison(regi, qry_coords=[s, e])
                    overlap_j = last_package.sub_comparison(regj, qry_coords=[s, e])
                    regi[6] = [lq for lq in regi[6] if lq[0] <overlap_i[2] or lq[1] > overlap_i[3]]
                    regj[6] = [lq for lq in regj[6] if lq[0] <overlap_j[2] or lq[1] > overlap_j[3]]
                    repeats.append(overlap_i[1:4] + [0])
                    repeats.append(overlap_j[1:4] + [0])
            
        # identify repetitive regions in the reference
        comparisons.sort(key=lambda x: x[2] )
        comparisons.sort(key=lambda x: x[1] )
    
        for id, regi in enumerate(comparisons) :
            if len(regi) == 0 : continue
            for jd in xrange(id+1, len(comparisons)) :
                regj = comparisons[jd]
                if regi[1] != regj[1] : break
                si, ei = sorted(regi[2:4])
                sj, ej = sorted(regj[2:4])
                s = max(si, sj)
                e = min(ei, ej)
                if e >= s :
                    overlap_i = last_package.sub_comparison(regi, ref_coords=[s, e])
                    overlap_j = last_package.sub_comparison(regj, ref_coords=[s, e])
                    if len(overlap_i) == len(overlap_j) :
                        diff = 0
                        for i, i_snp in enumerate(overlap_i[13:]) :
                            j_snp = overlap_j[13+i]
                            if i_snp[0] != j_snp[0] or i_snp[5] != j_snp[5] :
                                diff = 1
                                break
                        if diff == 1 :
                            for mut in regi[13:] :
                                if abs(mut[0]) <= e and abs(mut[1]) >= s :
                                    mut[6] = 1
                            for mut in regj[13:] :
                                if abs(mut[0]) <= e and abs(mut[1]) >= s :
                                    mut[6] = 1
                            regi[6] = [lq for lq in regi[6] if lq[0] <s or lq[1] > e]
                            regj[6] = [lq for lq in regj[6] if lq[0] <s or lq[1] > e]
                            repeats.append([regi[1], s, e, 0])
                    else :
                        for mut in regi[13:] :
                            if abs(mut[0]) <= e and abs(mut[1]) >= s :
                                mut[6] = 1
                        for mut in regj[13:] :
                            if abs(mut[0]) <= e and abs(mut[1]) >= s :
                                mut[6] = 1
                        regi[6] = [lq for lq in regi[6] if lq[0] <s or lq[1] > e]
                        regj[6] = [lq for lq in regj[6] if lq[0] <s or lq[1] > e]
                        repeats.append([regi[1], s, e, 0])
            for mut in regi[13:] :
                if mut[6] == 0 :
                    if regi[1] not in mutations : 
                        mutations[ regi[1] ] = {}
                    if mut[0] not in mutations[ regi[1] ]:
                        mutations[ regi[1] ] [ mut[0] ] = {}
                    if mut[5] not in mutations[ regi[1] ] [ mut[0] ] :
                        mutations[ regi[1] ] [ mut[0] ] [ mut[5] ] = [regi[7], regi[10]] + mut
                    else : 
                        mutations[ regi[1] ] [ mut[0] ] [ mut[5] ].extend([regi[7], regi[10]] + mut)
            repeats.extend([[regi[1]]+ lq[:2] + [1] for lq in regi[6]])
        
        repeats.sort(key=lambda x:x[1])
        repeats.sort(key=lambda x:x[0])
        repetitive_regions = []
        for rep in repeats:
            if len(repetitive_regions) == 0 or repetitive_regions[-1][0] != rep[0] or repetitive_regions[-1][2]+1 < rep[1] :
                repetitive_regions.append(rep)
            elif rep[2] > repetitive_regions[-1][2] :
                repetitive_regions[-1][2] = rep[2]
                if repetitive_regions[-1][3] > 0 :
                    repetitive_regions[-1][3] = rep[3]
        nocall = {}
        for r in repetitive_regions + [c[1:] for c in comparisons if float(c[0])/(abs(c[3]-c[2])+1) < 0.7 or c[0] < 200] :
            for s in xrange(r[1], r[2]+1) :
                nocall[(r[0], s)] = 1
        mutations = { contig:{ site:alters for site, alters in variation.items() if (contig, site) not in nocall } for contig, variation in mutations.items() }
        comparisons = [c for c in comparisons if float(c[0])/(abs(c[3]-c[2])+1) >= 0.7 and c[0] >= 200]
        return comparisons, repetitive_regions, mutations

    @staticmethod
    def write_down(filename, regions, repeats, mutations, reference, query, tag) :
        with uopen(filename, 'w') as fout:
            fout.write('##gff-version 3\n')
            fout.write('## Reference: {0}\n'.format(reference))
            fout.write('## Query: {0}\n'.format(query))
            fout.write('## Tag: {0}\n'.format(tag))

            fout.write('\n'.join(['{0}\trefMapper\tmisc_feature\t{1}\t{2}\t{3}\t{4}\t.\t{5}'.format(r[1], r[2], r[3], r[0], r[10], '/inference="Aligned%20with%20{0}:{1}-{2}"'.format(*r[7:10])) for r in regions]) + '\n')
            fout.write('\n'.join(['{0}\trefMapper\tunsure\t{1}\t{2}\t.\t+\t.\t/inference="{3}"'.format(r[0], r[1], r[2], 'Repetitive%20region' if r[3] == 0 else 'Uncertain%20base%20calling%20or%20ambigious%20alignment') for r in repeats]) + '\n')
            for contig, variation in sorted(mutations.items()):
                for site, alters in sorted(variation.items()) :
                    for alter, source in alters.items() :
                        if source[6][0] == '-' :
                            difference = '+{0}'.format(source[7])
                            origin = '.'
                        elif source[7][0] == '-' :
                            difference = '-{0}'.format(source[6])
                            origin = source[6]
                        else :
                            difference = source[7]
                            origin = source[6]
                        compare = ''
                        for id in xrange(0, len(source), 9) :
                            compare += '{0}:{1}-{2}:{3};'.format(source[id+0], abs(source[id+4]), abs(source[id+5]), source[id+1])
                        fout.write('{0}\trefMapper\tvariation\t{1}\t{2}\t.\t+\t.\t{3}\n'.format(contig, source[2], source[3], '/replace="{0}";/compare="{1}";/origin="{2}"'.format(difference, compare[:-1], origin)))


def lastAgainst(qry_tag, qry_file, refdb, prefix, ref_file, lastal) :
    if not os.path.isfile('{0}.{1}.lastal'.format(prefix, qry_tag)) :
        output = last_package.run_lastal(refdb, qry_file, '{0}.{1}.lastal'.format(prefix, qry_tag), lastal )
    else :
        output = '{0}.{1}.lastal'.format(prefix, qry_tag)
    regions, repeats, mutations = last_package.make_alignment( output )
    os.unlink(output)
    outfile = '{0}.gff.gz'.format(prefix)
    last_package.write_down(outfile, regions, repeats, mutations, ref_file, qry_file, qry_tag)
    return [qry_tag, outfile]


def alignAgainst(data) :
    prefix, aligner, db, (rtag, reference), (tag, query) = data
    if isinstance(aligner, list) :
        return lastAgainst(tag, query, db, prefix, reference, aligner[1])
    try :
        qrySeq, qryQual = readFastq(query)
    except :
        return [tag, query]
    if os.path.isfile( '{0}.gff.gz'.format(prefix) ) :
        return [tag, prefix + '.gff.gz']
    refSeq, refQual = readFastq(reference)
    
    if not divergent :
        proc = subprocess.Popen('{0} -k13 -w5 -c -t1 --frag=yes -A1 -B14 -O24,60 -E2,1 -r100 -g1000 -P -N5000 -f1000,5000 -n2 -m50 -s200 -z200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(
                                    aligner, db, query).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    else :
#        print('diver')
        proc = subprocess.Popen('{0} -k13 -w5 -c -t1 --frag=yes -A1 -B4 -O8,16 -E2,1 -r20k,40k --rmq -g10k -P -N5000 -f1000,5000 -n2 -m50 -s100 -z200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(
                                    aligner, db, query).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

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
            qryRepeat = [[q[1]-1, q[0]-1, -1, -1] for q in qryRepeat]
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
                q1 = qry[1][alnSite[1]:alnSite[1]+cl] if d > 0 else ''.join(reversed(qry[1][(alnSite[1]-cl+1):(alnSite[1]+1)]))

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
                q1 = qry[1][alnSite[1]:alnSite[1]+cl] if d > 0 else ''.join(reversed(qry[1][(alnSite[1]-cl+1):(alnSite[1]+1)] ))
                
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
                for site in xrange(alnSite[0], alnSite[0]+len(r)) :
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
    #print(mTag, mFile)
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
    return mTag, presences, sorted(absences), mutations

def readGFF(fname) :
    cds = []
    names = {}
    fnames = fname.split(',')
    fname = fnames[0]

    for fn in fnames :
        with uopen(fn) as fin :
            for line in fin :
                if line.startswith('#') :
                    continue
                else :
                    part = line.strip().split('\t')
                    if len(part) > 2 :
                        name = re.findall(r'locus_tag=([^;]+)', part[8])
                        if len(name) == 0 :
                            parent = re.findall(r'Parent=([^;]+)', part[8])
                            if len(parent) and parent[0] in names :
                                name = names[parent[0]]
                        # if len(name) == 0 :
                        #     name = re.findall(r'Name=([^;]+)', part[8])
                        if len(name) == 0 :
                            name = re.findall(r'ID=([^;]+)', part[8])
                        if part[2] == 'CDS' :
                            assert len(name) > 0, 'Error: CDS has no name. {0}'.format(line)
                            gname = name[0]
                            # if gname not in cds :
                            cds.append([gname, part[0], int(part[3]), int(part[4]), part[6], 0, [], int(part[3]), int(part[4]), fname])
                            # cds[gname] = [gname, part[0], int(part[3]), int(part[4]), part[6], 0, [], int(part[3]), int(part[4])]
                            # elif part[0] == cds[gname][1] and fname == cds[gname][0] :
                            #     cds[gname].extend([int(part[3]), int(part[4])])
                            #     cds[gname][3] = max(cds[gname][3], int(part[4]))
                        else :
                            ids = re.findall(r'ID=([^;]+)', part[8])
                            if len(ids) :
                                names[ids[0]] = name
    return sorted(cds, key=lambda c:c[1:3])



def getMatrix(prefix, reference, alignments, lowq_aligns, core, matrixOut, alignmentOut, gff, gcode) :
    genes = readGFF(gff) if gff else {}
    refSeq, refQual = readFastq(reference[1])
    coreSites = { n:np.zeros(len(refSeq[n]), dtype=int) for n in refSeq }
    matSites = { n:np.zeros(len(refSeq[n]), dtype=int) for n in refSeq }
    alnId = { aln[0]:id for id, aln in enumerate(alignments+lowq_aligns) }

    matrix = {}
    presence_dict = {}
    for gid, aln in enumerate((alignments, lowq_aligns)) :
        for mTag, presences, absences, mutations in pool.imap_unordered(readMap, aln) :
            j = alnId[mTag]
            presence_dict[(gid, j)] = [presences, absences]

            for mut in mutations :
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
    
    for (gid, j), (presences, absences) in presence_dict.items() :
        for n, s, e in presences :
            if gid == 0 :
                coreSites[n][s-1:e] +=1
            mutations = matSites[n][s-1:e]
            for kk in mutations[mutations > 0] :
                k = (n, kk)
                if len(matrix[k][0]) and matrix[k][0][j] == '-' :
                    matrix[k][0][j] = '.'
                if len(matrix[k][1]) and matrix[k][1][j] == '-' :
                    matrix[k][1][j] = '.'
        for n, s, e, m in absences :
            if gid == 0 :
                coreSites[n][s-1:e] -= (1 if mTag not in reference else len(alignments)+1)
            mutations = matSites[n][s-1:e]
            for kk in mutations[mutations > 0] :
                k = (n, kk)
                if len(matrix[k][0]) and matrix[k][0][j] == '.' :
                    matrix[k][0][j] = '-'
                if len(matrix[k][1]) and matrix[k][1][j] == '.' :
                    matrix[k][1][j] = '-'
    
    coreNum = max(len(alignments) * core, 1)
    for n in sorted(coreSites) :
        sites = coreSites[n]
        for site, num in enumerate(sites) :
            cSite = (n, site+1)
            if 1 <= num < coreNum and cSite in matrix and len(matrix[cSite][1]) > 0 :
                sites[site] = np.sum(matrix[cSite][1] != '-')
                matrix[cSite][0] = []
    
    pres = np.concatenate(list(coreSites.values()))
    pres[pres < 0] = 0
    pres = list(np.unique(pres, return_counts=True))
    
    for p, n in zip(*pres) :
        sys.stderr.write('#{2} {0} {1}\n'.format(p if p >= 1 else 'Repetitive', n, '' if p >= coreNum else '#'))

    missings = []
    coreBases = {'A':0, 'C':0, 'G':0, 'T':0}
    for n in sorted(coreSites) :
        sites = coreSites[n]
        for site, num in enumerate(sites) :
            cSite = (n, site+1)
            if num < coreNum :
                if cSite in matrix :
                    indel = matrix[cSite][1]
                    if sum(x != '-' for x in indel) >= coreNum :
                        matrix[cSite][0] = []
                    else :
                        matrix.pop(cSite, None)
                if len(missings) == 0 or missings[-1][0] != n or missings[-1][2] + 1 < cSite[1] :
                    missings.append([n, cSite[1], cSite[1]])
                else :
                    missings[-1][2] = cSite[1]
            else :
                b = refSeq[n][cSite[1]-1]
                if cSite in matrix and len(matrix[cSite][0]) :
                    matrix[cSite][0] = [ (b if s == '.' else s) for s in matrix[cSite][0]]
                    t = np.unique(matrix[cSite][0])
                    t = t[t!= '-']
                    if len(t) == 1 :
                        coreBases[t[0]] = coreBases.get(t[0], 0) + 1
                        matrix[cSite][0] = []
                else :
                    coreBases[b] = coreBases.get(b, 0) + 1

    gff_info = get_mut_gff_info(genes, matrix, refSeq, gcode) if len(genes) else {}

    outputs = {}
    if matrixOut :
        outputs['matrix'] = prefix + '.matrix.gz'
        with uopen(prefix + '.matrix.gz', 'w') as fout :
            fout.write('## Constant_bases: {A} {C} {G} {T}\n'.format(**coreBases))
            for n in refSeq :
                fout.write('## Sequence_length: {0} {1}\n'.format(n, len(refSeq[n])))
            for region in missings :
                fout.write('## Missing_region: {0} {1} {2}\n'.format(*region))
            fout.write('\t'.join(['#Seq', '#Site'] + [ mTag for mTag, mFile in alignments + lowq_aligns ] + (['#Genes:Site:Category'] if len(genes) else []) ) + '\n')
            for site in sorted(matrix) :
                bases = matrix[site]
                if len(bases[0]) :
                    fout.write('{0}\t{1}\t{2}\t{3}\n'.format(site[0], site[1], '\t'.join(bases[0]), gff_info.get((site[0], site[1]), '')))
                if len(bases[1]) :
                    fout.write('{0}\t{1}\t{2}\t{3}\n'.format(site[0], site[1], '\t'.join(bases[1]), re.sub('(Nonsyn|Syn|Nonsense):[A-Z,]+', 'Indel', gff_info.get((site[0], site[1]), '')) ))
    if alignmentOut :
        outputs['alignment'] = prefix + '.fasta.gz'
        sequences = []
        low_seq = []
        for sequence, aln, r in ((sequences, alignments, res), (low_seq, lowq_aligns, low_res)) :
            for (mTag, mFile), (presences, absences, mutations) in zip(aln, r) :
                j = alnId[mTag]
                seq = { n:['-']*len(s) for n, s in refSeq.items() } if j > 0 else { n:list(s) for n, s in refSeq.items() }
                if j :
                    for n, s, e in presences :
                        seq[n][s-1:e] = refSeq[n][s-1:e]
                    for n, s, e, c in absences :
                        seq[n][s-1:e] = '-' * (e-s+1)
                    for n, s, e in missings :
                        seq[n][s - 1:e] = '-' * (e - s + 1)
                for site in matrix :
                    bases = matrix[site]
                    if len(bases[0]) :
                        seq[site[0]][site[1]-1] = bases[0][j]
                sequence.append(seq)
        with uopen(prefix + '.fasta.gz', 'w') as fout :
            for id, n in enumerate(sorted(refSeq)) :
                if id :
                    fout.write('\n')
                for (mTag, mFile), seq in zip(alignments, sequences) :
                    fout.write('>{0}:{1}\n{2}\n'.format(mTag, n, ''.join(seq[n])))
                for (mTag, mFile), seq in zip(lowq_aligns, low_seq) :
                    fout.write('>{0}:{1}\n{2}\n'.format(mTag, n, ''.join(seq[n])))
    return outputs


def get_mut_gff_info(genes, matrix, refSeq, gcode) :
    gene = genes # sorted([ [g] + info[1:] for g, info in genes.items()], key=lambda g:g[1:4])
    gi = 0
    site_genes = {}
    for (cont, site), (mut, indel) in sorted(matrix.items()) :
        site_genes[(cont, site)] = []
        while gi < len(gene) and (cont > gene[gi][1] or ((cont == gene[gi][1]) and (site > gene[gi][3] + 120))) :
            gi += 1
        for i in np.arange(gi, len(gene)) :
            g = gene[i]
            if g[1] > cont or g[2] - 120 > site :
                break
            if site < g[2] :
                outside = site - g[2] if g[4] == '+' else g[2] - site
                if outside < 0 :
                    site_genes[(cont, site)].append(['o', outside, g[0], 'Upstream'])
                continue
            elif site > g[3] :
                outside = site - g[3] if g[4] == '+' else g[3] - site
                if outside < 0 :
                    site_genes[(cont, site)].append(['o', outside, g[0], 'Upstream'])
                continue
            in_gene = False
            for s, e in zip(g[7::2], g[8::2]) :
                if s <= site and site <= e :
                    in_gene = True
                    t = get_dNdS(refSeq, cont, s, e, g[4], site, mut, gcode) if len(mut) else 'Indel'
                    site_genes[(cont, site)].append(['g', site - s+1 if g[4] == '+' else e - site+1, g[0], t])
                    break
            if not in_gene :
                site_genes[(cont, site)].append(['i', site - g[2]+1 if g[4] == '+' else g[3] - site+1, g[0], 'Intergenic'])
    res = {}
    for cSite, gg in site_genes.items():
        gg = sorted(gg)
        if len(gg) > 0:
            gg = [g for g in gg if g[0] == gg[0][0]]
        res[cSite] = ';'.join(sorted(['{2}:{1}:{3}'.format(*g) for g in gg]))
    return res

def get_dNdS(refSeq, cont, s, e, d, site, mut, gcode) :
    if gcode == 4 :
        gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSWCWCLFLF-'))
    else :
        gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF-'))

    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    encode = {'A': 0, 'T': 3, 'G': 2, 'C': 1, 'N': -1000}
    seq = list(refSeq[cont][s-1:e])
    if d == '-' :
        seq = [comp.get(c, 'N') for c in seq[::-1]]
        mut = [comp.get(m, 'N') for m in mut]
    mut = _collections.OrderedDict([[encode.get(m, -1000), 1] for m in mut])
    i = (site - s) if d == '+' else (e - site)
    i0 = int(i/3)*3
    ci = i - i0
    codon = [encode.get(c, -1000) for c in np.array(seq)[i0:i0+3]]
    aas = _collections.OrderedDict()
    for m in mut :
        codon[ci] = m
        if len(codon) == 3 :
            x = codon[0]*16 + codon[1]*4 + codon[2]
            if x >= 0 and x < len(gtable) :
                aas[gtable[x]] = 1
    return ('Nonsense:{0}'.format(','.join(aas.keys())) if 'X' in aas else 'Nonsyn:{0}'.format(','.join(aas.keys()))) if len(aas) > 1 else 'Syn:{0}'.format(','.join(aas.keys()))


def runAlignment(prefix, reference, queries, core, aligner) :
    #alignments = list(map(alignAgainst, [[prefix +'.' + query[0].rsplit('.', 1)[0] + '.' + str(id+1), aligner, prefix + '.mmi', reference, query] for id, query in enumerate(queries)]))
    alignments = pool.map(alignAgainst, [[prefix +'.' + query[0].rsplit('.', 1)[0] + '.' + str(id+1), aligner, prefix + '.mmi', reference, query] for id, query in enumerate(queries)])

    try :
        os.unlink(reference + '.mmi')
    except :
        pass
    return alignments


def prepReference(prefix, ref_tag, reference, aligner, pilercr, trf, **args) :
    def mask_tandem(fasta_file) :
        cmd = '{0} {1} 2 4 7 80 10 60 2000 -d -h -ngs'.format(trf, fasta_file)
        trf_run = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, universal_newlines=True)
    
        region = []
        for line in iter(trf_run.stdout.readline, r'') :
            if line[0] == '@' :
                cont_name = line[1:].strip().split()[0]
            else :
                part = line.split(' ',2)[:2]
                region.append([cont_name, int(part[0])-2, int(part[1])+2])
        return region
    
    def mask_crispr(fasta_file, prefix) :
        cmd = '{0} -in {1} -out {2}.crispr'.format(pilercr, fasta_file, prefix)
        subprocess.Popen(cmd.split(), stderr=subprocess.PIPE).communicate()
        summary_trigger = 0
    
        region = []
        with open('{0}.crispr'.format(prefix)) as fin :
            for line in fin :
                if line.startswith('SUMMARY BY POSITION') :
                    summary_trigger = 1
                elif summary_trigger :
                    if line[0] == '>' :
                        cont_name = line[1:].strip().split()[0]
                    elif len(line) > 10 and line.strip()[0] in '0123456789' :
                        part = line[24:].strip().split()
                        region.append([cont_name, int(part[0]), int(part[0]) + int(part[1]) -1])
        os.unlink('{0}.crispr'.format(prefix))
        return region
    # prepare reference
    if reference :
        if not isinstance(aligner, list) :
            subprocess.Popen('{0} -k15 -w5 -d {2}.mmi {1}'.format(aligner, reference, prefix).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        else :
            subprocess.Popen('{0} -cR01 {2}.mmi {1}'.format(aligner[0], reference, prefix).split()).communicate()
        import tempfile
        with tempfile.NamedTemporaryFile(dir='.') as tf :
            seq, _ = readFastq(reference)
            tf_fas = '{0}.fasta'.format(tf.name)
            with open(tf_fas, 'wt') as fout:
                for n, s in seq.items() :
                    fout.write('>{0}\n{1}\n'.format(n, s))
            repeats = mask_tandem(tf_fas) + mask_crispr(tf_fas, tf.name)
            os.unlink(tf_fas)
        alignments = alignAgainst([prefix +'.' + ref_tag.rsplit('.', 1)[0] + '.0', aligner, prefix + '.mmi', [ref_tag, reference], [ref_tag, reference]])
        with uopen(alignments[1], 'a') as fout :
            for r in repeats :
                fout.write('{0}\trefMapper\tunsure\t{1}\t{2}\t.\t+\t.\t/inference="repetitive_regions"\n'.format(
                    r[0], r[1], r[2], 
                ))
    return alignments

divergent = False

def align(argv) :
    args = parseArgs(argv)
    #print(divergent)
    if args.qry_list == None :
        queries = args.queries
    else :
        queries = []
        with open(args.qry_list, 'rt') as fin :
            for line in fin :
                queries.append(line.strip().split()[0])
    queries = sorted([ [qt, qf] for qt, qf in [ qry.split(':', 1) if qry.find(':')>0 else [os.path.basename(qry), qry] for qry in queries ] if qt != args.reference[0] ])


    global pool
    pool = Pool(args.n_proc)
    refMask = prepReference(args.prefix, args.reference[0], args.reference[1], args.aligner, **externals)
    alignments = runAlignment(args.prefix, args.reference, queries, args.core, args.aligner)
    lowq_aligns = runAlignment(args.prefix, args.reference, args.lowq, args.core, args.aligner)
    alignments = [refMask] + alignments
    outputs = {'mappings': dict(alignments), 'low_qual_map': dict(lowq_aligns)}
    if args.matrix or args.alignment :
        outputs.update(getMatrix(args.prefix, args.reference, alignments, lowq_aligns, args.core, args.matrix, args.alignment, args.gff, args.gcode))
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
    

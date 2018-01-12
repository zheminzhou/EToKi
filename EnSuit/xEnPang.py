import sys, shlex, re, marshal, math, ete3
from operator import itemgetter
import subprocess, os, numpy as np
from collections import Counter
from time import gmtime, strftime
from multiprocessing import Pool

codons = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
          "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
          "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
          "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
          "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
          "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
          "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
          "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
          "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
          "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
          "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
          "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
          "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
          "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
          "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
          "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

codons4 = dict(codons.items())
codons4['TGA'] = 'W'

def rc(seq) :
    return ''.join([complement.get(b, 'N') for b in reversed(seq.upper())])

def transeq(seq, frame=1, transl_table=11) :
    gtable = codons if transl_table != 4 else codons4
    if str(frame).upper() == 'F' :
        frames = [1,2,3]
    elif str(frame).upper() == 'R' :
        frames = [4,5,6]
    elif str(frame).upper() == '6' :
        frames = range(1,7)
    else :
        frames = [int(frame)]

    if isinstance(seq, dict) :
        trans_seq = {}
        for n,s in seq.iteritems() :
            #print n, len(s)
            for frame in frames :
                trans_name = '{0}_{1}'.format(n, frame)
                if frame <= 3 :
                    trans_seq[trans_name] = ''.join([gtable.get(c, 'X') for c in map(''.join, zip(*[iter(s[(frame-1):])]*3))])
                else :
                    trans_seq[trans_name] = ''.join([gtable.get(c, 'X') for c in map(''.join, zip(*[iter(rc(s)[(frame-4):])]*3))])
        return trans_seq
    else :
        if len(frames) == 1 :
            frame = frames[0]
            if frame <= 3 :
                return ''.join([gtable.get(c, 'X') for c in map(''.join, zip(*[iter(seq[(frame-1):])]*3))])
            else :
                return ''.join([gtable.get(c, 'X') for c in map(''.join, zip(*[iter(rc(seq)[(frame-4):])]*3))])
        else :
            trans_seq = []
            for frame in frames :
                if frame <= 3 :
                    trans_seq.append( ''.join([gtable.get(c, 'X') for c in map(''.join, zip(*[iter(seq[(frame-1):])]*3))]) )
                else :
                    trans_seq.append( ''.join([gtable.get(c, 'X') for c in map(''.join, zip(*[iter(rc(seq)[(frame-4):])]*3))]) )
            return trans_seq

def get_similar_pairs(self_bsn) :
    def get_similar(bsn, ortho_pairs) :
        if ortho_pairs.get(bsn[0][0], {}).get(bsn[0][1], 0) == 1 :
            return
        matched_aa = {}
        len_aa = int(bsn[0][12])/3
        for part in bsn : 
            s_i, e_i, s_j, e_j = [ int(x) for x in part[6:10] ]
            frame_i, frame_j = s_i % 3, s_j % 3
            if part[14].find('-') < 0 and part[15].find('-') < 0 :
                if frame_i == frame_j :
                    matched_aa.update({site:1 for site in xrange((s_i+1)/3 *3+1, e_i, 3)})
            else :
                s_i, s_j = s_i - 1, s_j - 1
                for b_i, b_j in zip(part[14], part[15]) :
                    if b_i != '-' :
                        s_i += 1
                        if b_j != '-' :
                            s_j += 1
                            if s_i % 3 == 1 and s_j % 3 == 1 :
                                matched_aa[s_i] = 1
                    elif b_j != '-' :
                        s_j += 1
            if len(matched_aa) >= 300/3 or len(matched_aa) >= 0.5 * len_aa :
                if part[1] not in ortho_pairs :
                    ortho_pairs[part[1]] = {}
                if part[0] not in ortho_pairs :
                    ortho_pairs[part[0]] = {}
                ortho_pairs[part[1]] [part[0]] = 1
                ortho_pairs[part[0]] [part[1]] = 1
                return
            
    ortho_pairs = {}
    save = []
    with open(self_bsn) as fin :
        for id, line in enumerate(fin) :
            if id % 100000 == 0 :
                sys.stderr.write('Loading self-BLASTn results. {0}\n'.format(id))
            
            part = line.strip().split()
            if part[0] == part[1] or float(part[2]) < 50 or int(part[8]) > int(part[9]) or int(part[3]) < 40 :
                continue
            if len(save) > 0 and (save[0][0] != part[0] or save[0][1] != part[1]) :
                get_similar(save, ortho_pairs)
                save = []
            save.append(part)
    if len(save) > 0 :
        get_similar(save, ortho_pairs)
    return ortho_pairs

def load_matches(bsn_file) :
    matches = []
    with open(bsn_file) as bsn :
        for id, line in enumerate(bsn) :
            if id % 100000 == 0 :
                sys.stderr.write('Loading BLASTn results. {1} : {0}\n'.format(len(matches), id))
            part = line.strip().split()
            if float(part[2]) < params['min_identity'] or ( int(part[7]) - int(part[6]) + 1 < 60 and int(part[7]) - int(part[6]) + 1 < 0.5 * float(part[12]) ) :
                continue
            strain, cont = part[1].split('__', 1)
            if int(part[8]) < int(part[9]) :
                save = [part[0], strain, cont, float(part[2]), int(part[6]), int(part[7]), int(part[8]), int(part[9]), int(part[11]), int(part[12]), int(part[13]), len(matches), 0, part[15] if part[14].find('-') < 0 else ''.join([ y for x,y in zip(part[14], part[15]) if x != '-'])]
            else :
                save = [part[0], strain, cont, float(part[2]), int(part[6]), int(part[7]), -int(part[8]), -int(part[9]), int(part[11]), int(part[12]), int(part[13]), len(matches), 0, part[15] if part[14].find('-') < 0 else ''.join([ y for x,y in zip(part[14], part[15]) if x != '-'])]
            matches.append(save)
    return matches

def combine_matches(save) :
    def resolve_edges(edges) :
        groups = []
        for m1 in edges[0] :
            for m2 in edges[1] :
                new_len = m2[5] - m1[4] + 1
                if new_len >= params['clust_match_prop'] * m1[9] or (new_len >= 0.6666* m1[9] and new_len >= 300) :
                    if m1[2] == m2[2] :
                        s1, e1 = (m1[7], m1[8]) if m1[7] > 0 else (-m1[8], -m1[7])
                        s2, e2 = (m2[7], m2[8]) if m2[7] > 0 else (-m2[8], -m2[7])
                        overlap = min(e1, e2) - max(s1, s2) + 1
                        if overlap >= 300 or overlap >= 0.3333*min(e1-s1+1, e2-s2+1) :
                            continue
                    overlap = m1[5]-m2[4]+1
                    if overlap >= 500 or overlap <= -3000 or overlap >= 0.8*min(m1[5]-m1[4]+1, m2[5]-m2[4]+1) or -overlap >= 50+min(m1[5]-m1[4]+1, m2[5]-m2[4]+1) :
                        continue
                    if overlap > 0 :
                        score = m1[8] + m2[8] - overlap * min( float(m1[8])/(m1[5]-m1[4]+1), float(m2[8])/(m2[5]-m2[4]+1) )
                        ident = (m1[3]*(m1[5]-m1[4]+1) + m2[3]*(m2[5]-m2[4]+1) - overlap*min(m1[3], m2[3]))/((m1[5]-m1[4]+1)+(m2[5]-m2[4]+1) - overlap)
                    else :
                        score = m1[8] + m2[8]
                        ident = (m1[3]*(m1[5]-m1[4]+1) + m2[3]*(m2[5]-m2[4]+1))/((m1[5]-m1[4]+1)+(m2[5]-m2[4]+1))
                    groups.append( [new_len, score, ident ] + ([m1[11], m2[11]] if m1[3] < m2[3] else [m2[11], m1[11]]) )
                    m1[12], m2[12] = 1, 1
        return groups
    
    save.sort(key=itemgetter(0,1,2,6))

    groups = {}
    prev = ['', '']
    edges =[[], []]
    len_match = len(save)
    ref_lens = {}
    for id, m1 in enumerate(save) :
        if id % 10000 == 0 :
            sys.stderr.write('Get hits: {0} / {1}\n'.format(id, len_match))
        ref_lens[m1[0]] = m1[9]
        if m1[0] not in groups :
            groups[m1[0]] = {}
        if m1[1] not in groups[m1[0]] :
            groups[m1[0]][m1[1]] = []
        if m1[5]-m1[4]+1 >= params['clust_match_prop']*m1[9] or (m1[5] - m1[4] + 1 >= 0.6666 * m1[9] and m1[5] - m1[4] + 1 >= 300) or (m1[5] - m1[4] + 1 >= 0.5 * m1[9] and m1[5] - m1[4] + 1 >= 540) :  # minimum length requirement for a non-break gene
            groups[m1[0]][m1[1]].append([m1[5] - m1[4] + 1, m1[8], m1[3], m1[11]])
            m1[12] = 1
        if prev != m1[:2] :
            if prev[0] in groups :
                groups[prev[0]][prev[1]].extend(resolve_edges(edges))       # find genes that were fragmented due to assembly gaps
            prev = m1[:2]
            edges =[[], []]
        if m1[4] > 20 and ((m1[6] > 0 and m1[6] - 1 <= 150) or (m1[6] < 0 and m1[10] + m1[6] <= 150)) :   # any hit within the last 150 bps to either end of a scaffold is a potential fragmented gene
            edges[1].append(m1)
        if m1[5] < m1[9]-20 :
            if (m1[6] > 0 and m1[10]-m1[7] <= 150) or (m1[6] < 0 and -1-m1[7] <= 150) :
                edges[0].append(m1)
            for jd in xrange(id+1, len_match) :
                m2 = save[jd]
                if m1[:3] != m2[:3] or (m1[6] < 0 and m2[6] > 0) or m2[6] - m1[7]-1 >= 300 :    # maximum 300bps between two continuous hits in the same scaffold
                    break
                if m1[7] >= m2[7] or m1[4] > m2[4] or m1[5] > m2[5] or m2[4] - m1[5]-1 >= 300 or min([m1[5]-m1[4]+1, m2[5]-m2[4]+1]) < max(m2[6]-m1[7]-1, m2[4]-m1[5]-1)*0.6666 :
                    continue
                new_len = m2[5] - m1[4] + 1
                if new_len >= params['clust_match_prop'] * m1[9] or (new_len >= 0.6666 * m1[9] and new_len >= 300) :
                    overlap = max(m1[7]-m2[6]+1, m1[5]-m2[4]+1)   ## to improve
                    if overlap >= 0.8*(m1[5]-m1[4]+1) or overlap >= 0.8*(m2[5]-m2[4]+1) :
                        continue
                    if overlap > 0 :
                        score = m1[8] + m2[8] - overlap * min( float(m1[8])/(m1[5]-m1[4]+1), float(m2[8])/(m2[5]-m2[4]+1) )
                        ident = (m1[3]*(m1[5]-m1[4]+1) + m2[3]*(m2[5]-m2[4]+1) - overlap*min(m1[3], m2[3]))/((m1[5]-m1[4]+1)+(m2[5]-m2[4]+1) - overlap)
                    else :
                        score = m1[8] + m2[8]
                        ident = (m1[3]*(m1[5]-m1[4]+1) + m2[3]*(m2[5]-m2[4]+1))/((m1[5]-m1[4]+1)+(m2[5]-m2[4]+1))
                    groups[m1[0]][m1[1]].append( [ new_len, score, ident ] + ([m1[11], m2[11]] if m1[3] < m2[3] else [m2[11], m1[11]]) )
                    m1[12], m2[12] = 1, 1
    if prev[0] in groups :
        groups[prev[0]][prev[1]].extend(resolve_edges(edges))

    max_score = {gene:max([g[:3] for genome, group in res.iteritems() for g in group] + [[0, 1, 100.0]], key=itemgetter(1))for gene, res in groups.iteritems()}
    for gene in max_score :
        max_score[gene] = [float(max_score[gene][0])/max_score[gene][1], max_score[gene][2]/100.0]
    for gene, res in groups.iteritems() :
        for genome, group in res.iteritems() :
            for grp in group :
                grp[:3] = [1, max_score[gene][0] * grp[1] * ((float(grp[0])/ref_lens[gene])**2), min(grp[2]/max_score[gene][1], 100.0)]
            group.sort(key=lambda x:x[1], reverse=True)
            group[:] = [grp for grp in group if grp[1] >= 0.3333 * group[0][1] or grp[2] >= group[0][2]-10]
    match_inuse = {mm: 1 for group in groups.itervalues() for mat in group.itervalues() for m in mat for mm in m[3:]}
    for mat in matches :
        if mat[11] not in match_inuse :
            mat[12:] = [0]
    return groups

def get_conflicts(matches, ortho_pairs) :
    conflicts = {}
    matches.sort(key=lambda x:x[6] if x[6]>0 else -x[7])
    matches.sort(key=itemgetter(1,2))
    len_match = len(matches)
    for id, m1 in enumerate(matches) :
        if id % 10000 == 0 :
            sys.stderr.write('Get conflicts: {0} / {1}\n'.format(id, len_match))
        if m1[12] == 0 : continue
        s1, e1 = m1[6:8] if m1[6] > 0 else (-m1[7], -m1[6])
        
        for jd in xrange(id+1, len_match) :
            m2 = matches[jd]
            if m1[1:3] != m2[1:3] :
                break
            if m2[12] == 0 : continue
            s2, e2 = m2[6:8] if m2[6] > 0 else (-m2[7], -m2[6])
            if s2 > e1-10 :
                break
            overlap = min(e1, e2) - s2 + 1
            if overlap >= 200 or overlap >= 0.5 * (e1-s1+1) or overlap >= 0.5 * (e2-s2+1) :
                if m1[11] not in conflicts :
                    conflicts[m1[11]] = {}
                if m2[11] not in conflicts :
                    conflicts[m2[11]] = {}
                if m1[0] == m2[0] :
                    conflicts[m1[11]][m2[11]] = 0
                    conflicts[m2[11]][m1[11]] = 0
                else :
                    if m2[0] in ortho_pairs.get(m1[0], {}) :
                        conflicts[m1[11]][m2[11]] = 2
                    else :
                        conflicts[m1[11]][m2[11]] = 1
                        
                    if m1[0] in ortho_pairs.get(m2[0], {}) :
                        conflicts[m2[11]][m1[11]] = 2
                    else :
                        conflicts[m2[11]][m1[11]] = 1
        if m1[11] in conflicts :
            conflicts[m1[11]] = ' '.join( ['{0}={1}'.format(x, y) for x, y in conflicts[m1[11]].iteritems()] )
    return conflicts

def get_seq(matches, m) :
    seq = np.array(['-' for id in xrange(matches[m[3]][9])])
    for id in m[3:] :
        mat = matches[id]
        seq[mat[4]-1:mat[5]] = list(mat[13])
    return seq
decode = {'A':1, 'C':2, 'G':3, 'T':4, '-':0, 'N':0}
def compare_seq(seqs, mod = 1) :
    genome = seqs.keys()
    seq = np.array( [[decode.get(b, 0)for b in seqs[g]] for g in genome] )
    diff = {}
    for id, m in enumerate(genome) :
        comparable = (seq[id] != 0) * (seq[(id+1):] != 0)
        n_comparable = np.sum(comparable, 1)
        n_comparable[n_comparable < 1] = 1
        n_diff = np.sum(( seq[id] != seq[(id+1):] ) & comparable, 1).astype(float)/n_comparable
        
        if mod == 1 :
            n_diff[n_comparable < seq.shape[1]*0.5] = 2
            for jd, nd in enumerate(n_diff) :
                if nd <= 1.0 :
                    n = genome[id+jd+1]
                    key = (m, n) if m < n else (n, m)
                    diff[key] = nd
        else :
            n_diff[n_comparable < min(200, seq.shape[1]*0.2)] = 2
            for jd, nd in enumerate(n_diff) :
                n = genome[id+jd+1]
                key = (m, n) if m < n else (n, m)
                diff[key] = nd
    return diff

def global_difference(matches, groups, premature=None) :
    uclust_ref = readFasta(params['uclust'])
    global_differences = {}
    grp_order = sorted([ (grp, len([m for genome, m in mat.iteritems() if len(m) == 1])) for grp, mat in groups.iteritems() ], key=itemgetter(1), reverse=True)[:10000]
    for id, (grp, cnt) in enumerate(grp_order) :
        mat = groups[grp]
        if id % 1000 == 0 :
            sys.stderr.write('Get global similarities: {0} / {1}\n'.format(id, len(grp_order)))
            if premature is not None and id >= premature :
                break
        seqs = { genome:get_seq(matches, m[0]) for genome, m in mat.iteritems() if len(m) == 1 }
        diff = compare_seq(seqs)
        len_ref = len(uclust_ref[grp])
        for pair, difference in diff.iteritems() :
            if pair not in global_differences: 
                global_differences[pair] = [0, 0]
            global_differences[pair][0] += len_ref * difference
            global_differences[pair][1] += len_ref
    return {pair:info[0]/info[1] for pair, info in global_differences.iteritems()}

def filt_per_group(data) :
    seqs, grp, mat, tags, global_differences = data
    diff = compare_seq(seqs, mod=2)
    incompatible, distances, cleanup = {}, {}, {}
    dc = math.pow(1-params['clust_difference'], 1.0/params['distance_cut'])
    for r1, r2 in diff :
        if diff[(r1, r2)] < 1 :
            difference = math.log(min((1.0 - diff[(r1, r2)]) / (1.0 - global_differences.get((r1[0], r2[0]), 0.01)), 1.0 )) 
            if difference < 0.0 :
                difference = difference / math.log( min(1-global_differences.get((r1[0], r2[0]), 0.01), dc) )
            distances[(r1, r2)] = difference
            distances[(r2, r1)] = difference
        else :
            difference = params['distance_cut'] + 0.5
        if difference >= params['distance_cut'] :
            k1, k2 = '{0}__{1}'.format(*r1), '{0}__{1}'.format(*r2)
            if k1 not in incompatible : incompatible[k1] = {}
            if k2 not in incompatible : incompatible[k2] = {}
            incompatible[k1][k2] = incompatible[k2][k1] = 1

            if mat[r1[0]][int(r1[1])][2] - mat[r2[0]][int(r2[1])][2] >= 10 :
                cleanup[r2[0]] = min(cleanup.get(r2[0], 99999), int(r2[1]))
            elif mat[r1[0]][int(r1[1])][2] - mat[r2[0]][int(r2[1])][2] <= -10 :
                cleanup[r1[0]] = min(cleanup.get(r1[0], 99999), int(r1[1]))
    if len(incompatible) > 0 :
        if params['orthology'] == 'similarity' :
            todel = {}
            incompatible_regions = sorted([tuple(x.rsplit('__', 1)) for x in incompatible.keys()], key=lambda x:mat[x[0]][int(x[1])][1], reverse=True)
            for r1 in incompatible_regions :
                if r1 not in todel :
                    for r2 in incompatible['{0}__{1}'.format(*r1)] :
                        todel[tuple(r2.rsplit('__', 1))] = 1
            for r in sorted(todel, key=lambda x:int(x[1]), reverse=True) :
                del mat[r[0]][int(r[1])]
        else :
            groups, group_tag = [], {}
            for n, s in sorted(seqs.iteritems(), key=lambda x:mat[x[0][0]][int(x[0][1])][1], reverse=True) :
                if n[0] in cleanup and int(n[1]) > cleanup[n[0]] + 2 :
                    continue
                novel = 1
                for tag, g in groups :
                    key = (tag, n) if tag < n else (n, tag)
                    if diff[key] <= params['clust_difference'] and len(seqs[tag]) + 10 >= len(seqs[n]) :
                        g.append(n)
                        group_tag['{0}__{1}'.format(*n)] = '{0}__{1}'.format(*tag)
                        novel = 0
                        break
                if novel :
                    tags['{0}__{1}'.format(*n)] = ''.join(s.tolist())
                    groups.append([n, [n]])
                    group_tag['{0}__{1}'.format(*n)] = '{0}__{1}'.format(*n)
            for g, id in cleanup.iteritems() :
                mat[g][:] = mat[g][:(id+3)]
            
            groups = {'{0}__{1}'.format(*g): grp for g, grp in groups}
            ic2 = {}
            for r1, rr in incompatible.iteritems() :
                if r1 in group_tag :
                    t1 = group_tag[r1]
                    if t1 not in ic2 :
                        ic2[t1] = {}
                    for r2 in rr.iterkeys() :
                        if r2 in group_tag :
                            ic2[t1][group_tag[r2]] = 1
            incompatible = ic2
    
            cmd = params[params['orthology']] if len(tags) < 1000 else params['nj']
            phy_run = subprocess.Popen(shlex.split(cmd), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            gene_phy = ete3.Tree(phy_run.communicate(input=''.join(['>{0}\n{1}\n'.format(n, s) for n, s in tags.iteritems()]))[0])
            node = gene_phy.get_midpoint_outgroup()
            if node is not None :
                gene_phy.set_outgroup(node)
            s_name = set(tags.keys()) - set(['{0}__REF'.format(grp)])
    
            for ite in xrange(3000) :
                tips = set(gene_phy.get_leaf_names()) & s_name
                branches = []
                for node in gene_phy.get_descendants('preorder') :
                    descendants = set(node.get_leaf_names()) & s_name
                    if len(descendants) == 0 :
                        continue
                    unrelated = tips - descendants
                    
                    d_anchor, u_anchor, n_ic = {}, {}, 0
                    for r1 in descendants :
                        for r2 in incompatible.get(r1, {}) :
                            if r2 in unrelated :
                                n_ic += 1
                                d_anchor[r1], u_anchor[r2] = 1, 1
                    if len(d_anchor) == 0 :
                        continue
                    branches.append([n_ic, node.dist + node.get_sisters()[0].dist if node.up.is_root() else node.dist, node, descendants, unrelated, d_anchor, u_anchor])
                branches = sorted(branches, reverse=True)
                cuts = []
                for br in (branches[:40], branches[40:100], branches[100:200], branches[200:500], branches[500:1000], branches[1000:]) :
                    for n_ic, d, node, descendants, unrelated, d_anchor, u_anchor in br :
                        inter_dist, inter_cnt = 0, 0
                        for g1, g2 in ((d_anchor, unrelated), (u_anchor, descendants - set(d_anchor))) :
                            for rr1 in g1 :
                                for rr2 in g2 :
                                    for r1 in groups[rr1] :
                                        for r2 in groups[rr2] :
                                            if (r1, r2) in distances :
                                                inter_dist += distances[(r1, r2)]
                                                inter_cnt += 1
                        if inter_cnt > 0 and inter_dist/inter_cnt >= params['distance_cut'] :
                            cuts.append([inter_dist/inter_cnt, d, node])
                    if len(cuts) :
                        break
                if len(cuts) > 0 :
                    cut_node = max(cuts)[2]
                    prev_node = cut_node.up
                    cut_node.detach()
                    if grp + '__REF' in cut_node.get_leaf_names() :
                        gene_phy = cut_node
                    elif prev_node.is_root() :
                        gene_phy = gene_phy.get_children()[0]
                    else :
                        prev_node.delete(preserve_branch_length=True)
                        
                    for r1 in incompatible.keys() :
                        if r1 not in tips :
                            rr = incompatible.pop(r1, None)
                            for r2 in rr :
                                incompatible.get(r2, {}).pop(r1, None)
                    if len(incompatible) == 0 :
                        break
                    
                    sys.stderr.write('     Iteration {0}. Remains {1} tips.\n'.format(ite+1, len(gene_phy.get_leaf_names())))
                else :
                    break
            if len(gene_phy.get_leaf_names()) < len(tags) :
                tips = sorted([ [nn[0], int(nn[1])] for n in gene_phy.get_leaf_names() for nn in groups.get(n, [])])
                new_mat = { g:[] for g in mat }
                for genome, id in tips :
                    new_mat[genome].append(mat[genome][id])
                mat = new_mat
    return mat

def readFasta(fname) :
    seq = {}
    with open(fname) as fin :
        for line in fin :
            if line[0] == '>' :
                name = line[1:].strip().split()[0]
                seq[name] = []
            else :
                seq[name].append(line.strip())
    for n in seq :
        seq[n] = ''.join(seq[n])
    return seq

def get_gene(groups, first_classes, cnt=1) :
    ranking = {gene:first_classes[gene] for gene in groups if gene in first_classes} if first_classes is not None else {}
    if len(ranking) > 0 :
        min_rank = min(ranking.values())
        scores = { gene:0 for gene, score in ranking.iteritems() if score == min_rank }
    else :
        min_rank = -1
        scores = { gene:0 for gene in groups }
    scores = { gene:sum([group[0][1]*group[0][0] for genome, group in groups[gene].iteritems() if len(group) > 0 ]) for gene in scores }
    
    genes = [[gene, score, min_rank] for gene, score in sorted(sorted(scores.iteritems()), key=itemgetter(1), reverse=True)[:cnt] if score > 0]
    if len(genes) <= 0 :
        for gene in scores :
            groups.pop(gene, None)
        return []
    return genes

def filt_genes(prefix, matches, groups, global_differences, conflicts, first_classes = None) :
    uclust_ref = readFasta(params['uclust'])
    
    pangenes, used, results, run = {}, {}, {}, {}
    
    group_id = 0
    with open('{0}.PopuPr'.format(prefix), 'wb') as fout :
        while len(groups) > 0 :
            genes = get_gene(groups, first_classes, cnt=10)
            if len(genes) <= 0 :
                continue
            to_run, to_run_id, min_score, min_rank = [], [], genes[-1][1], genes[0][2]
            genes = {gene:score for gene, score, min_rank in genes}
            for gene, score in genes.iteritems() :
                if gene not in run :
                    mat = groups[gene]
                    region_score = [region[1]/m[0][1] for genome, m in mat.iteritems() for region in m]
                    if params['orthology'] == 'ref_dist' :
                        for genome, m in mat.iteritems() :
                            groups[gene][genome][:] = [ region for region in m if region[1]/m[0][1] >= 1.0 - params['clust_difference'] ]
                        run[gene] = 1
                    else :
                        if len(region_score) > len(mat) * 2 :
                            used2 = {}
                            for genome, m in mat.iteritems() :
                                remove = []
                                for id, mm in enumerate(m) :
                                    for gid in mm[3:] :
                                        if gid in used2 :
                                            remove.append(id)
                                            break
                                    if id not in remove :
                                        for gid in mm[3:] :
                                            used2[gid] = 1
                                for id in reversed(remove) :
                                    m.pop(id)
                            region_score = [region[1]/m[0][1] for genome, m in mat.iteritems() for region in m]
                        if len(region_score) > len(mat) * 3 :
                            cut = sorted(region_score, reverse=True)[len(mat)*3]
                            if cut >= 1.0 - params['clust_difference'] :
                                cut = min(sorted(region_score, reverse=True)[len(mat)*5] if len(region_score) > len(mat) * 5 else 1.0 - params['clust_difference'], 1.0 - 0.6*(params['clust_difference']))
                            for genome, m in mat.iteritems() :
                                m[:] = [ region for region in m if region[1]/m[0][1] > cut ]
                        for genome, m in mat.iteritems() :
                            m.sort(key=itemgetter(1), reverse=True)
                        
                        to_run.append([{ (genome, str(i)):get_seq(matches, region) for genome, m in mat.iteritems() for i,region in enumerate(m) }, \
                                                   gene, mat, {'{0}__REF'.format(gene):uclust_ref[gene]}, global_differences])
                        to_run_id.append(gene)
            working_groups = pool.map(filt_per_group, to_run)
            #working_groups = [filt_per_group(d) for d in to_run]
            for gene, working_group in zip(to_run_id, working_groups) :
                for g, m in working_group.iteritems() :
                    for mm in m :
                        for id in mm[3:] :
                            if len(mm) == 4 and matches[id][3] >= 85 :
                                matches[id][12] = 2
                            else :
                                matches[id][12] = 1
                        m.sort(key=lambda x:x[0]*x[1], reverse=True)
                groups[gene] = working_group
                run[gene] = 1

            while len(genes) :
                score, gene = max([[sum([group[0][1]*group[0][0] for genome, group in groups[g[0]].iteritems() if len(group) > 0 ]), g[0]] for g in genes.iteritems()])
                if score < min_score :
                    break
                working_group = groups.pop(gene)
                genes.pop(gene)

                supergroup = {}
                paralog, reg_cnt = 0, 0
                todel = {}
                for genome, region in working_group.iteritems() :
                    for reg in region :
                        reg_cnt += 1
                        skip = 0
                        for gid in reg[3:] :
                            conflict = used.get(gid, None)
                            if conflict is not None :
                                skip = 1
                                if isinstance(conflict, basestring) :
                                    supergroup[conflict] = supergroup.get(conflict, 0) + 1
                                elif conflict > 0 and reg[2] > 85 and len(reg) <= 4 :
                                    paralog = 1
                        if not skip :
                            for gid in reg[3:] :
                                used[gid] = 0
                                if gid in conflicts :
                                    for jd, cat in [c.split('=') for c in conflicts[gid].split(' ')] :
                                        jd, cat = int(jd), int(cat)
                                        if matches[jd][0] not in todel :
                                            todel[matches[jd][0]] = {}
                                        if matches[jd][1] not in todel[matches[jd][0]] :
                                            todel[matches[jd][0]][matches[jd][1]] = {}
                                        todel[matches[jd][0]][matches[jd][1]][jd] = cat
                                        if cat == 0 and jd not in used :
                                            used[jd] = 0
                        elif paralog :
                            break
                        else :
                            reg[:] = []
                if paralog :
                    continue
                if len(supergroup) > 0 :
                    pangene = max(supergroup.items(), key=itemgetter(1))
                    cnt = min(reg_cnt, len([1 for m in results[pangene[0]].values() for mm in m if matches[mm[3]][0] == pangene[0]]))
                    if (pangene[1] >= 2 and pangene[1] >= 0.1 * cnt) or (pangene[1] >= 1 and pangene[1] >= 0.5 * cnt) :
                        pangene = pangene[0]
                    else :
                        pangene = gene
                else :
                    pangene = gene
        
                if pangene not in results : 
                    results[pangene] = {}
        
                sys.stderr.write('[{3}] {5} / {6}: pan gene "{4}" : "{0}" picked from rank {1} and score {2}\n'.format(gene, min_rank, score, strftime("%Y-%m-%d %H:%M:%S", gmtime()), pangene, len(results), len(groups)+len(results)))
                
                for genome, region in working_group.iteritems() :
                    if genome not in results[pangene] :
                        results[pangene][genome] = []
                    results[pangene][genome].extend([reg for reg in region if len(reg) > 0])
                
                for genome, region in sorted(working_group.iteritems()) :
                    for grp in region :
                        if len(grp) > 0 :
                            group_id += 1
                            for id in sorted(grp[3:], key=lambda x:matches[x][4]) :
                                fout.write('{0}\t{1}\t{2}\t{3}\n'.format(pangene, min_rank, group_id, '\t'.join([str(x) for x in matches[id][:-2]])))
                
                todel2 = {}
                for gene, lists in todel.iteritems() :
                    if gene not in groups :
                        continue
                    for genome, ids in lists.iteritems() :
                        for grp in groups[gene].get(genome, []) :
                            for id in grp[3:] :
                                if id in ids :
                                    grp[0] = 0
                                    if ids[id] == 1 and grp[2] >= 100.0*(1.0 - params['clust_difference']) and len(grp) <= 4 :
                                        todel2[gene] = 1
                        groups[gene][genome].sort(key=lambda x:x[0]*x[1], reverse=True)
                        for id, cat in ids.iteritems() :
                            if cat == 1 :
                                used[id] = 1
                            elif id not in used :
                                used[id] = pangene
                for gene in todel2 :
                    groups.pop(gene, None)
                    genes.pop(gene, None)
    return results

def load_priority(priority_file) :
    first_classes = {}
    with open(priority_file) as fin :
        for line in fin :
            part = line.strip().split()
            if len(part) > 0 and part[0] not in first_classes :
                first_classes[part[0]] = int(part[1]) if len(part) > 1 else 0
    return first_classes

def get_uclust(prefix, gene, priority) :
    max_prior = max(priority.values()) + 1 if priority != {} else 1
    seq = []
    with open(gene) as fin :
        for line in fin :
            if line[0] == '>' :
                name = line[1:].strip().split()[0]
                seq.append([name, priority.get(name, max_prior), 0, []])
            else :
                seq[-1][3].append(line.strip())
    for s in seq :
        s[3] = ''.join(s[3])
        s[2] = len(s[3])
    seq.sort(key=lambda x:x[2], reverse=True)
    seq.sort(key=lambda x:x[1])
    with open('{0}.gene.sorted'.format(prefix), 'wb') as fout :
        for s in seq :
            if s[1] < 0 : continue
            if s[1] == max_prior and s[2] < 200 :
                break
            ts = transeq(s[3])
            if 'X' not in ts[:-1] :
                fout.write( '>{0}\n{3}\n'.format(*s) )
            else :
                sys.stderr.write('{0} is taken out because it has internal stop codon or uncertain bases.\n'.format(s[0]))
    subprocess.Popen('{0} -sortedby other -cluster_smallmem {1}.gene.sorted -uc {1}.gene.uc -centroids {1}.gene.reference -id {2} -query_cov {3} -target_cov {3}'.format(params['usearch'], prefix, params['clust_difference'], params['clust_match_prop']).split()).wait()
    os.unlink('{0}.gene.sorted'.format(prefix))
    return '{0}.gene.uc'.format(prefix), '{0}.gene.reference'.format(prefix)

def get_self_bsn(prefix, uclust) :
    subprocess.Popen('{0}makeblastdb -dbtype nucl -in {1}'.format(params['blast_folder'], uclust).split()).communicate()
    subprocess.Popen(shlex.split('{0}blastn -num_threads {3} -query {1} -db {1} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen qseq sseq" -max_target_seqs 2000 -out {2}.self.bsn -task blastn -reward 2 -penalty -3 -evalue 0.1'.format(params['blast_folder'], uclust, prefix, params['n_thread']))).communicate()
    return '{0}.self.bsn'.format(prefix)

def get_map_bsn(prefix, uclust, genome) :
    subprocess.Popen('{0}makeblastdb -dbtype nucl -in {1}'.format(params['blast_folder'], genome).split()).communicate()
    subprocess.Popen(shlex.split('{0}blastn -num_threads {4} -query {1} -db {2} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen qseq sseq" -max_target_seqs 2000 -out {3}.map.bsn -task blastn -reward 2 -penalty -3 -evalue 0.1'.format(params['blast_folder'], uclust, genome, prefix, params['n_thread']))).communicate()
    return '{0}.map.bsn'.format(prefix)

params = dict(
    usearch = '/home/zhemin/biosoft/usearch8.0.1623_i86linux32', 
    blast_folder = '/home/zhemin/bin/', 
    prefix = 'OPP', 
    ml = '/home/zhemin/biosoft/FastTreeMP -nt -gtr -pseudo', 
    nj = '/home/zhemin/biosoft/FastTreeMP -nt -noml -gtr -pseudo', 
    n_thread = '30', 
    orthology = 'nj', # ml, nj, ref_dist, similarity
    distance_cut = '2.0', 
    min_identity = '60', 
    clust_difference = '0.05', 
    clust_match_prop = '0.9'
)

if __name__ == '__main__' :
    params.update(dict([ p.split('=', 1) for p in sys.argv[1:] ]))
    if 'priority' in params :
        first_classes = load_priority(params['priority'])
    else :
        first_classes = {}
    if params['blast_folder'] != '' and params['blast_folder'][-1] != '/' :
        params['blast_folder'] += '/'
    params['distance_cut'] = float(params['distance_cut'])
    params['min_identity'], params['clust_match_prop'] = float(params['min_identity']), float(params['clust_match_prop'])
    params['clust_difference'] = max(0.05, float(params['clust_difference']))
    
    if 'uclust' not in params :
        params['uc'], params['uclust'] = get_uclust(params['prefix'], params['genes'], first_classes)
    if 'self_bsn' not in params :
        params['self_bsn'] = get_self_bsn(params['prefix'], params['uclust'])
    if 'map_bsn' not in params :
        params['map_bsn'] = get_map_bsn(params['prefix'], params['uclust'], params['genome'])

    pool = Pool(5)

    if 'homolog' not in params :
        matches = load_matches(params['map_bsn'])
        groups = combine_matches(matches)
        conflicts = get_conflicts(matches, get_similar_pairs(params['self_bsn']))
        marshal.dump([matches, groups, conflicts], open(params['prefix'] + '.homolog.dump', 'wb'))
    else :
        matches, groups, conflicts = marshal.load(open(params['homolog'], 'rb'))
        conflicts = { int(id): cmp for id, cmp in conflicts.iteritems() }
    matches.sort(key=lambda x:x[11])
    if 'global' not in params :
        params['global'] = params['prefix'] + '.global.dump'
        global_differences = global_difference(matches, groups)
        marshal.dump({ pair:'{0:.40f}'.format(value) for pair, value in global_differences.iteritems()}, open(params['global'], 'wb'))
    global_differences = {pair:max(float(value), 0.01) for pair, value in marshal.load(open(params['global'], 'rb')).iteritems()}
    global_differences.update({(pair[1], pair[0]):value  for pair, value in global_differences.iteritems() })
    results = filt_genes(params['prefix'], matches, groups, global_differences, conflicts, first_classes)

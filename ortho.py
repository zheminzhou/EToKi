import sys, shlex, re, marshal, math, ete3, hashlib, tempfile
from operator import itemgetter
import subprocess, os, numpy as np
from collections import Counter
from time import gmtime, strftime
from multiprocessing import Pool
from configure import externals, logger, rc, transeq, readFasta

def iter_readGFF(fname) :
    seq, cds = {}, {}
    p = subprocess.Popen(['zcat', fname], stdout = subprocess.PIPE) if fname.upper().endswith('.GZ') else subprocess.Popen(['cat', fname], stdout = subprocess.PIPE)
    for line in p.stdout :
        if line.startswith('#') : continue
        if line.startswith('>') :
            name = line[1:].strip().split()[0]
            assert name not in seq, logger('Error: duplicated sequence name {0}'.format(name))
            seq[name] = [fname, []]
        else :
            part = line.strip().split('\t')
            if len(part) > 1 :
                if part[2] == 'CDS' :
                    name = re.findall(r'ID=([^;]+)', part[8])
                    if name is None :
                        name = re.findall(r'Name=([^;]+)', part[8])
                    assert name is not None, logger('Error: CDS has no name. {0}'.format(line))
                    assert name[0] not in cds, logger('Error: duplicated CDS. {0}'.format(line))
                    cds[name[0]] = [fname, part[0], int(part[3]), int(part[4]), part[6], 0, '']
            else :
                seq[name][1].append(part[0])
    p.wait()

    for n, s in seq.iteritems() :
        s[1] = ''.join(s[1]).upper()
    for n, c in cds.iteritems() :
        try:
            c[6]= seq[c[1]][1][(c[2]-1) : c[3]]
            if c[4] == '-' :
                c[6] = rc(c[6])
            if not checkCDS(n, c[6]) :
                c[6] = ''
            else :
                c[5] = int(hashlib.sha1(c[6]).hexdigest(), 16)
        except :
            pass
    
    return seq, cds
    

def readGFF(fnames) :
    if isinstance(fnames, basestring) : fnames = [fnames]
    combo = pool.map(iter_readGFF, fnames)

    seq, cds = {}, {}
    for ss, cc in combo :
        for n, s in ss.iteritems() :
            seq[n] = s[:]
    
        for n, c in cc.iteritems() :
            if n not in cds and len(c[6]) > 0 :
                cds[n] = c[:]
    return seq, cds

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
            if len(matched_aa)*3 >= min(params['match_len2'], params['match_len']) or len(matched_aa) >= (min(params['match_prop'], params['match_prop2'])-0.1) * len_aa :
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
                logger('Loading self-BLASTn results. {0}'.format(id))

            part = line.strip().split()
            if part[0] == part[1] or float(part[2]) < params['match_identity']-10. or int(part[8]) > int(part[9]) or int(part[3]) < params['match_frag_len']-10 :
                continue
            if len(save) > 0 and (save[0][0] != part[0] or save[0][1] != part[1]) :
                get_similar(save, ortho_pairs)
                save = []
            save.append(part)
    if len(save) > 0 :
        get_similar(save, ortho_pairs)
    return ortho_pairs

def load_matches(bsn_file, genomes) :
    matches = []
    with open(bsn_file) as bsn :
        for id, line in enumerate(bsn) :
            if id % 100000 == 0 :
                logger('Loading BLASTn results. {1} : {0}'.format(len(matches), id))
            part = line.strip().split()
            if float(part[2]) < 100.*params['match_identity'] or ( int(part[7]) - int(part[6]) + 1 < params['match_frag_len'] and int(part[7]) - int(part[6]) + 1 < params['match_frag_prop'] * float(part[12]) ) :
                continue
            cont = part[1]
            strain = genomes[cont][0]
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
                if new_len >= params['clust_match_prop'] * m1[9] or (new_len >= params['match_prop2']* m1[9] and new_len >= params['match_len2']) :
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
            logger('Get hits: {0} / {1}'.format(id, len_match))
        ref_lens[m1[0]] = m1[9]
        if m1[0] not in groups :
            groups[m1[0]] = {}
        if m1[1] not in groups[m1[0]] :
            groups[m1[0]][m1[1]] = []
        match_len = m1[5]-m1[4]+1
        if match_len >= params['clust_match_prop']*m1[9] or (match_len >= params['match_prop'] * m1[9] and match_len >= params['match_len']) or \
           (match_len >= params['match_prop2'] * m1[9] and match_len >= params['match_len2']) :  # minimum length requirement for a non-break gene
            groups[m1[0]][m1[1]].append([match_len, m1[8], m1[3], m1[11]])
            m1[12] = 1
        if prev != m1[:2] :
            if prev[0] in groups :
                groups[prev[0]][prev[1]].extend(resolve_edges(edges))       # find genes that were fragmented due to assembly gaps
            prev = m1[:2]
            edges =[[], []]
        if m1[4] > 20 and ((m1[6] > 0 and m1[6] - 1 <= params['edge_rescue']) or (m1[6] < 0 and m1[10] + m1[6] <= params['edge_rescue'])) :   # any hit within the last 150 bps to either end of a scaffold is a potential fragmented gene
            edges[1].append(m1)
        if m1[5] < m1[9]-20 :
            if (m1[6] > 0 and m1[10]-m1[7] <= params['edge_rescue']) or (m1[6] < 0 and -1-m1[7] <= params['edge_rescue']) :
                edges[0].append(m1)
            for jd in xrange(id+1, len_match) :
                m2 = save[jd]
                if m1[:3] != m2[:3] or (m1[6] < 0 and m2[6] > 0) or m2[6] - m1[7]-1 >= params['synteny_gap'] :    # maximum 300bps between two continuous hits in the same scaffold
                    break
                if m1[7] >= m2[7] or m1[4] > m2[4] or m1[5] > m2[5] or m2[4] - m1[5]-1 >= params['synteny_gap'] or min([m1[5]-m1[4]+1, m2[5]-m2[4]+1])*params['synteny_diff'] < max(m2[6]-m1[7]-1, m2[4]-m1[5]-1) :
                    continue
                new_len = m2[5] - m1[4] + 1
                if new_len >= params['clust_match_prop'] * m1[9] or (new_len >= params['match_prop2'] * m1[9] and new_len >= params['match_len2']) :
                    overlap = max(m1[7]-m2[6]+1, m1[5]-m2[4]+1)   ## to improve
                    if overlap >= params['synteny_ovl_len'] or overlap >= params['synteny_ovl_prop']*min(m1[5]-m1[4]+1, m2[5]-m2[4]+1) :
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
            group[:] = [grp for grp in group if grp[1]*3. >= group[0][1] or grp[2] >= group[0][2]-10]
    match_inuse = {mm: 1 for group in groups.itervalues() for mat in group.itervalues() for m in mat for mm in m[3:]}
    for mat in save :
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
            logger('Get conflicts: {0} / {1}'.format(id, len_match))
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
            if overlap >= params['match_len'] or overlap >= params['match_prop'] * min(e1-s1+1, e2-s2+1) :
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
def compare_seq(seqs) :
    if len(seqs) < 1 :
        return {}
    genome = seqs.keys()
    seq = np.vectorize(decode.get)( [seqs[g] for g in genome], 0 )
    diff = np.zeros(shape=[seq.shape[0], seq.shape[0], 2])
    res = {}
    for id, s in enumerate(seq) :
        comparable = (s > 0) * (seq[(id+1):] > 0)
        n_comparable = np.sum(comparable, 1)
        n_comparable[n_comparable < 1] = 1        
        n_diff = np.sum(( s != seq[(id+1):] ) & comparable, 1)
        diff[id, id+1:, 0] = n_diff
        n_comparable[n_comparable < min(params['match_len2'], seq.shape[1]*params['match_prop2'])] = 0
        diff[id, id+1:, 1] = n_comparable
    for id, m in enumerate(genome) :
        for jd in xrange(id+1, len(genome)) :
            n = genome[jd]
            key = (m, n) if m < n else (n, m)
            res[key] = diff[id, jd, :]
    return res

def global_difference(matches, groups) :
    global_differences = {}
    grp_order = sorted([ (grp, len([m for genome, m in mat.iteritems() if len(m) == 1])) for grp, mat in groups.iteritems() ], key=itemgetter(1), reverse=True)[:5000]
    for id, (grp, cnt) in enumerate(grp_order) :
        if cnt == 0 : break
        mat = groups[grp]
        if id % 1000 == 0 :
            logger('Get global similarities: {0} / {1}'.format(id, len(grp_order)))
        seqs = { genome:get_seq(matches, m[0]) for genome, m in mat.iteritems() if len(m) == 1 }
        diff = compare_seq(seqs)
        for pair, (mut, aln) in diff.iteritems() :
            if pair not in global_differences:
                global_differences[pair] = [0, 0]
            global_differences[pair][0] += mut
            global_differences[pair][1] += aln
    return {pair:info[0]/info[1] for pair, info in global_differences.iteritems()}

def filt_per_group(data) :
    seqs, grp, mat, tags, global_differences = data
    diff = compare_seq(seqs)
    incompatible, distances = {}, {}
    for (r1, r2), (mut, aln) in diff.iteritems() :
        if aln > 0 :
            distances[(r1, r2)] = distances[(r2, r1)] = max(0., 1-(aln - mut)/aln/(1 - global_differences.get((r1[0], r2[0]), 0.001)) )
            ave = max(params['mutation_variation']*global_differences.get((r1[0], r2[0]), 0.01)*aln, 1.)
            difference = (mut - ave)/np.sqrt(ave)/2.807
        else :
            difference = 1.5
            distances[(r1, r2)] = distances[(r2, r1)] = 0.8
        if difference > 1. :
            k1, k2 = '{0}__{1}'.format(*r1), '{0}__{1}'.format(*r2)
            if k1 not in incompatible : incompatible[k1] = {}
            if k2 not in incompatible : incompatible[k2] = {}
            incompatible[k1][k2] = incompatible[k2][k1] = 1

    if len(incompatible) > 0 :
        groups, group_tag = [], {}
        for n, s in sorted(seqs.iteritems(), key=lambda x:mat[x[0][0]][int(x[0][1])][1], reverse=True) :
            novel = 1
            for tag, g in groups :
                key = (tag, n) if tag < n else (n, tag)
                if diff[key][0] <= 0.8*params['clust_difference']*diff[key][1] and len(seqs[tag]) + 10 >= len(seqs[n]) :
                    g.append(n)
                    group_tag['{0}__{1}'.format(*n)] = '{0}__{1}'.format(*tag)
                    novel = 0
                    break
            if novel :
                tags['{0}__{1}'.format(*n)] = ''.join(s.tolist())
                groups.append([n, [n]])
                group_tag['{0}__{1}'.format(*n)] = '{0}__{1}'.format(*n)
                
        if len(tag) < 20000 :
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
    
            tmpFile = tempfile.NamedTemporaryFile(dir='.', delete=False)
            for n, s in tags.iteritems() :
                tmpFile.write('>{0}\n{1}\n'.format(n, s))
            tmpFile.close()
            cmd = params[params['orthology']].format(tmpFile.name, **params) if len(tags) < 500 else params['nj'].format(tmpFile.name, **params)
            phy_run = subprocess.Popen(shlex.split(cmd), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            gene_phy = ete3.Tree(phy_run.communicate()[0].replace("'", ''))
            os.unlink(tmpFile.name)
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
                        if inter_cnt > 0 :
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
    
                    logger('     Iteration {0}. Remains {1} tips.'.format(ite+1, len(gene_phy.get_leaf_names())))
                else :
                    break
            if len(gene_phy.get_leaf_names()) < len(tags) :
                tips = sorted([ [nn[0], int(nn[1])] for n in gene_phy.get_leaf_names() for nn in groups.get(n, [])])
                new_mat = { g:[] for g in mat }
                for genome, id in tips :
                    new_mat[genome].append(mat[genome][id])
                mat = new_mat
    return mat

def get_gene(groups, first_classes, cnt=1) :
    ranking = {gene:first_classes[gene][0] for gene in groups if gene in first_classes}
    if len(ranking) > 0 :
        min_rank = min(ranking.values())
        scores = { gene:sum( [ group[0][1]*group[0][0] for genome, group in groups[gene].iteritems() if len(group) > 0 ] ) 
                   for gene, score in ranking.iteritems() if score == min_rank }
    else :
        min_rank = -1
        scores = {}
    genes = [[gene, score, min_rank] for gene, score in sorted(sorted(scores.iteritems()), key=itemgetter(1), reverse=True)[:cnt] if score > 0]
    if len(genes) <= 0 :
        for gene in scores :
            groups.pop(gene, None)
        return []
    else :
        for gene, score in scores :
            if score <= 0 :
                groups.pop(gene, None)
    return genes

def filt_genes(prefix, matches, groups, global_differences, conflicts, first_classes = None) :
    vclust_ref = readFasta(params['vclust'])
    
    pangenes, used, results, run = {}, {}, {}, {}
    group_id = 0
    with open('{0}.Prediction'.format(prefix), 'wb') as fout :
        while len(groups) > 0 :
            genes = get_gene(groups, first_classes, cnt=50)
            if len(genes) <= 0 :
                continue
            to_run, to_run_id, min_score, min_rank = [], [], genes[-1][1], genes[0][2]
            genes = {gene:score for gene, score, min_rank in genes}
            if params['orthology'] in ('ml', 'nj') :
                for gene, score in genes.iteritems() :
                    if gene not in run :
                        mat = groups[gene]
                        for genome, m in mat.iteritems() :
                            m.sort(key=itemgetter(1), reverse=True)
                        region_score = [region[1]/m[0][1] for genome, m in mat.iteritems() for region in m]
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
                        if len(region_score) > len(mat) * 3 and len(region_score) > 500 :
                            cut = sorted(region_score, reverse=True)[len(mat)*3]
                            if cut >= 1.0 - params['clust_difference'] :
                                cut = min(sorted(region_score, reverse=True)[len(mat)*5] if len(region_score) > len(mat) * 5 else 1.0 - params['clust_difference'], 1.0 - 0.6*(params['clust_difference']))
                            for genome, m in mat.iteritems() :
                                m[:] = [ region for region in m if region[1]/m[0][1] >= cut ]
    
                        to_run.append([{ (genome, str(i)):get_seq(matches, region) for genome, m in mat.iteritems() for i,region in enumerate(m) }, \
                                                   gene, mat, {'{0}__REF'.format(gene):vclust_ref[gene]}, global_differences])
                        to_run_id.append(gene)
                working_groups = pool.map(filt_per_group, to_run)
                #working_groups = [filt_per_group(d) for d in to_run]
                for gene, working_group in zip(to_run_id, working_groups) :
                    for g, m in working_group.iteritems() :
                        m.sort(key=lambda x:x[0]*x[1], reverse=True)
                    groups[gene] = working_group
                    run[gene] = 1
            else :
                for gene, score in genes.iteritems() :
                    if gene not in run :
                        mat = groups[gene]
                        for genome, m in mat.iteritems() :
                            m.sort(key=itemgetter(1), reverse=True)
                        cut = 1.0 - params['clust_difference']
                        for genome, m in mat.iteritems() :
                            if len(m) > 1 :
                                m[:] = [ region for region in m if region[1]/m[0][1] >= cut and region[2]/m[0][2] >= cut ]

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



            while len(genes) :
                score, gene = max([[sum([group[0][1]*group[0][0] for genome, group in groups[g[0]].iteritems() if len(group) > 0 ]), g[0]] for g in genes.iteritems()])
                if score < min_score :
                    break
                working_group = groups.pop(gene)
                genes.pop(gene)

                paralog, reg_cnt = 0, 0
                supergroup, todel = {}, {}
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
                    if pangene[1] * 2 >= cnt or (pangene[1] >= 3 and pangene[1] * 3 >= cnt) :
                        pangene = pangene[0]
                    else :
                        pangene = gene
                else :
                    pangene = gene

                if pangene not in results :
                    results[pangene] = {}

                logger('{4} / {5}: pan gene "{3}" : "{0}" picked from rank {1} and score {2}'.format(gene, min_rank, score, pangene, len(results), len(groups)+len(results)))

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
    return '{0}.Prediction'.format(prefix)

def load_priority(priority_list, genes) :
    file_priority = { fn:id for id, fnames in enumerate(priority_list.split(',')) for fn in fnames.split(':') }
    unassign_id = max(file_priority.values()) + 1
    first_classes = { n:[ file_priority.get(g[0], unassign_id), -len(g[6]), g[5], n ] for n, g in genes.iteritems() }

    return first_classes

def get_vclust(prefix, genes) :
    subprocess.Popen('{0} --dbmask soft --usersort --cluster_smallmem {4} --uc {1}.gene.uc --centroids {1}.gene.reference --id {2} --query_cov {3} --target_cov {3} --threads {5}'.format(params['vsearch'], prefix, 1.0 - params['clust_difference'], params['clust_match_prop'], genes, params['n_thread']).split()).wait()
    #os.unlink('{0}.gene'.format(prefix))
    return '{0}.gene.uc'.format(prefix), '{0}.gene.reference'.format(prefix)

def get_self_bsn(prefix, vclust) :
    subprocess.Popen('{0} -dbtype nucl -in {1}'.format(params['formatdb'], vclust).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    subprocess.Popen(shlex.split('{0} -num_threads {3} -query {1} -db {1} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen qseq sseq" -max_target_seqs 2000 -out {2}.self.bsn -task blastn -reward 2 -penalty -3 -evalue 0.1'.format(params['blast'], vclust, prefix, int(params['n_thread'])))).communicate()
    return '{0}.self.bsn'.format(prefix)

def writeGenomes(fname, seqs) :
    with open(fname, 'wb') as fout :
        for g, s in seqs.iteritems() :
            fout.write('>{0} {1}\n{2}\n'.format(g, s[0], s[-1]))
    return fname

def iter_map_bsn(data) :
    prefix, vclust, id, taxon, seq = data
    gfile, bsn = '{0}.{1}.genome'.format(prefix, id), '{0}.{1}.bsn'.format(prefix, id)
    with open(gfile, 'wb') as fout :
        fout.write(''.join(seq))
    subprocess.Popen('{0} -dbtype nucl -in {1}'.format(params['formatdb'], gfile).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    subprocess.Popen(shlex.split('{0} -num_threads {4} -query {1} -db {2} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen qseq sseq" -max_target_seqs 2000 -out {3} -task blastn -reward 2 -penalty -3 -evalue 0.01'.format(\
        params['blast'], vclust, gfile, bsn, int(params['n_thread'])/10))).communicate()
    for fn in (gfile, gfile+'.nin', gfile+'.nsq', gfile+'.nhr') :
        os.unlink(fn)
    return bsn

def get_map_bsn(prefix, vclust, genomes) :
    taxa = {}
    for g, s in genomes.iteritems() :
        if s[0] not in taxa : taxa[s[0]] = []
        taxa[s[0]].append('>{0}\n{1}\n'.format(g, s[1]))
    
    bsns = pool.map(iter_map_bsn, [(prefix, vclust, id, taxon, seq) for id, (taxon, seq) in enumerate(taxa.iteritems())])
    subprocess.Popen('cat {0} > {1}.map.bsn'.format(' '.join(bsns), prefix), shell=True).communicate()
    for bsn in bsns :
        os.unlink(bsn)
    return '{0}.map.bsn'.format(prefix)

def checkCDS(n, s) :
    if len(s) < params['min_cds'] :
        logger('{0} is too short'.format(n))
        return False
    if len(s) % 3 > 0 :
        logger('{0} is discarded due to frameshifts'.format(n))
        return False
    aa = transeq({'n':s.upper()}, transl_table='starts')['n_1']
    if aa[0] != 'M' :
        logger('{0} is discarded due to lack of start codon'.format(n))
        return False
    if aa[-1] != 'X' :
        logger('{0} is discarded due to lack of stop codon'.format(n))
        return False
    if len(aa[:-1].split('X')) - 1 > params['source_premature'] :
        logger('{0} is discarded due to too many stop codons'.format(n))
        return False
    return True
    

def addGenes(genes, gene_file) :
    for gfile in gene_file.split(',') :
        if gfile == '' : continue
        ng = readFasta(gfile)
        for name, s in ng.iteritems() :
            if name not in genes and checkCDS(name, s) :
                genes[name] = [ gfile, '', 0, 0, '+', int(hashlib.sha1(s).hexdigest(), 16), s]
    return genes

def addGenomes(genomes, genome_file) :
    for gfile in genome_file.split(',') :
        if gfile == '' : continue
        ng = readFasta(gfile)
        for name, s in ng.iteritems() :
            if name not in genomes :
                #assert name not in genomes, 'Duplicated genome "{0}" in "genomes" parameter'.format(name)
                genomes[name] = [ gfile, s.upper() ]
    return genomes

def writeGenes(fname, genes, priority) :
    uniques = {}
    with open(fname, 'wb') as fout :
        for n in sorted(priority.iteritems(), key=itemgetter(1)) :
            s = genes[n[0]][6]
            len_s, hcode = len(s), genes[n[0]][5]
            if len_s not in uniques :
                uniques = { len_s:{hcode:1} }
            elif hcode in uniques[ len_s ] :
                continue
            uniques[ len_s ][ hcode ] = 1
            fout.write( '>{0}\n{1}\n'.format(n[0], s) )
    return fname

def perform_mcl(groups, matches, method, global_differences) :
    vclust_ref = readFasta(params['vclust'])
    
    species = {}
    for s in groups.itervalues() :
        for ss in s :
            if ss not in species :
                species[ss] = len(species)
    glb_diff = {(species[g1], species[g2]): d for (g1, g2), d in global_differences.iteritems()}
    with open('{prefix}.edges'.format(**params), 'wb') as fout :    
        for gene, group in groups.iteritems() :
            if method != 'pairwise' :
                for sp, grp in group.iteritems() :
                    for id, g in enumerate(grp) :
                        fout.write('{0}__REF\t{1}\t{2}\n'.format(gene, "{0}__{1}__{2}".format(gene, species[g], id), g[1]/len(vclust_ref)))
            else :
                seqs = { (species[sp], id):get_seq(matches, g) for sp, grp in group.iteritems() for id,g in enumerate(grp) }
                differences = compare_seq(seqs)
                for (g1, g2), (mut, aln) in differences.iteritems() :
                    if aln <= 0 :
                        diff = 1
                    else :
                        gd = glb_diff.get((g1[0], g2[0]), 0.015)
                        ave = max(params['mutation_variation']*gd*aln, 1.)
                        diff = (mut - ave)/np.sqrt(ave)/2.807
                    if diff < 1 :
                        score = min((1-mut/aln)/(1-gd), 1.)
                        fout.write('{0}\t{1}\t{2}\n'.format("{0}__{1}__{2}".format(gene, g1[0], g1[1]), "{0}__{1}__{2}".format(gene, g2[0], g2[1]), score))
    subprocess.Popen('{mcl} {prefix}.edges -I 5 --abc -o {prefix}.mcl'.format(**params).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
    
    mcl = {}
    with open('{prefix}.mcl'.format(**params)) as fin :
        for line in fin :
            part = line.strip().split('\t')
            for p in part :
                mcl[p] = part[0]
    for gene, group in groups.iteritems() :
        best_mcl = {}
        for sp, grp in group.iteritems() :
            for id, g in enumerate(grp) :
                key = "{0}__{1}__{2}".format(gene, species[sp], id)
                mcl_group = mcl.get(key, key)
                best_mcl[mcl_group] = max(best_mcl.get(mcl_group, [0., 0.]), [g[1], g[2]])
        if len(best_mcl) > 0 :
            best_mcl = max(best_mcl.items(), key=itemgetter(1))[0]
            for sp, grp in group.iteritems() :
                todel = []
                for id, g in enumerate(grp) :
                    key = "{0}__{1}__{2}".format(gene, species[sp], id)
                    mcl_group = mcl.get(key, key)
                    if mcl_group != best_mcl :
                        todel.append(id)
                for id in reversed(todel) :
                    del grp[id]
    return groups

def write_output(prefix, prediction, genomes, vclust_ref, old_prediction) :
    predictions, alleles = {}, {}
    
    allele_file = open('{0}.allele.fna'.format(prefix), 'wb')
    with open(prediction) as fin :
        for line in fin :
            part = line.strip().split()
            if part[0] not in alleles :
                alleles[part[0]] = {vclust_ref[part[0]]:1}
                allele_file.write('>{0}_{1}\n{2}\n'.format(part[0], 1, vclust_ref[part[0]]))
            part[7:14] = [int(p) for p in part[7:14]]
            if part[9] > 0 :
                l, r = min(part[7]-1, part[9]-1), min(part[12]-part[8], part[13]-part[10])
            else :
                l, r = min(part[7]-1, part[13]+part[9]), min(part[12]-part[8], -part[10]-1)
            if l <= 6 and part[7] - l == 1 :
                part[7], part[9] = part[7]-l, part[9]-l
            else :
                ll = (part[7]-1) % 3
                if ll > 0 :
                    part[7], part[9] = part[7]+3-ll, part[9]+3-ll
            if r <= 6 and part[8] + r == part[12] :
                part[8], part[10] = part[8]+r, part[10]+r
            else :
                rr = (part[12] - part[8]) % 3
                if rr > 0 :
                    part[8], part[10] = part[8]-3+rr, part[10]-3+rr

            if part[9] > 0 :
                part[9:12] = part[9], part[10], '+'
            else :
                part[9:12] = -part[10], -part[9], '-'
            
            if part[4] not in predictions :
                predictions[part[4]] = []
            elif predictions[part[4]][-1][2] == part[2] :
                prev = predictions[part[4]][-1]
                if prev[5] == part[5] and part[7] - prev[8] < 300 :
                    if part[11] == '+' and part[9] - prev[10] < 300 :
                        prev[8], prev[10] = part[8], part[10]
                        continue
                    elif part[11] == '-' and prev[9] - part[10] < 300 :
                        prev[8], prev[9] = part[8], part[9]
                        continue
                predictions[part[4]][-1][1], part[1] = -1, -1
            predictions[part[4]].append(part)
    
    op = ['', 0, []]
    with open('{0}.EnOrth.gff'.format(prefix), 'wb') as fout :
        for gid, (g, predict) in enumerate(predictions.iteritems()) :
            predict.sort(key=itemgetter(5, 9, 10))
            for pid, pred in enumerate(predict) :
                cds = 'CDS'
                if pred[1] == -1 or (pred[10]-pred[9]+1) <= 0.8 * pred[12] :
                    cds, allele_id = 'fragment:{0:.2f}%'.format((pred[10]-pred[9]+1)*100/pred[12]), 'uncertain'
                    start, stop = pred[9:11]
                    s, e = pred[9:11]
                else :
                    s, e = pred[9:11]
                    if pred[11] == '+' :
                        s2, e2 = s - min(3*((s - 1)/3), 30), e + min(3*((pred[13] - e)/3), 300)
                        seq = genomes[pred[5]][1][(s2-1):e2]
                        lp, rp = s - s2, e2 - e
                    else :
                        s2, e2 = s - min(3*((s - 1)/3), 300), e + min(3*((pred[13] - e)/3), 30)
                        seq = rc(genomes[pred[5]][1][(s2-1):e2])
                        rp, lp = s - s2, e2 - e
                    if pred[7] == 1 and pred[8] == pred[12] :
                        seq2 = seq[(lp):(len(seq)-rp)]
                        if seq2 not in alleles[pred[0]] :
                            alleles[pred[0]][seq2] = len(alleles[pred[0]])+1
                        allele_id = str(alleles[pred[0]][seq2])
                        allele_file.write('>{0}_{1}\n{2}\n'.format(pred[0], allele_id, seq2))
                    else :
                        allele_id = 'uncertain'
                        
                    if len(seq) % 3 > 0 :
                        cds = 'frameshift'
                        start, stop = s, e
                    else :
                        aa_seq = transeq({'n':seq}, transl_table='starts')['n_1']
                        s0, s1 = aa_seq.find('M', lp/3, lp/3+12), aa_seq.rfind('M', 0, lp/3)
                        start = s0 if s0 >= 0 else s1
                        if start < 0 :
                            cds, start = 'nostart', lp/3
                        stop = aa_seq.find('X', start)
                        if stop >= 0 and stop < lp/3+12 :
                            s0 = aa_seq.find('M', stop, lp/3+12)
                            if s0 >= 0 :
                                start = s0
                                stop = aa_seq.find('X', start)
                        if stop < 0 :
                            cds, stop = 'nostop', len(aa_seq) - rp/3 - 1
                        elif (stop - start + 1)*3 <= 0.8 * pred[12] :
                            cds = 'premature stop:{0:.2f}%'.format((stop - start + 1)*300/pred[12])
                            stop = len(aa_seq) - rp/3 - 1
                        if pred[11] == '+' :
                            start, stop = s2 + start*3, s2 + stop*3 + 2
                        else :
                            start, stop = e2 - stop*3 - 2, e2 - start*3

                if pred[5] != op[0] :
                    op = [pred[5], 0, old_prediction.get(pred[5], [])]
                old_tag = []
                for k in xrange(op[1], len(op[2])) :
                    opd = op[2][k]
                    if opd[2] < start :
                        op[1] = k + 1
                    elif opd[1] > stop :
                        break
                    elif opd[3] != pred[11] :
                        continue
                    ovl = min(opd[2], stop) - max(opd[1], start) + 1
                    if ovl >= 300 or ovl >= 0.6 * (opd[2]-opd[1]+1) or ovl >= 0.6 * (stop - start + 1) :
                        frame = min((opd[1] - start) % 3, (opd[2] - stop) % 3)
                        if frame == 0 :
                            old_tag.append('{0}:{1}-{2}'.format(*opd))

                fout.write('{0}\t{1}\tEnOrth\t{2}\t{3}\t.\t{4}\t.\tID={5};{12}inference=ortholog group:{6},allele ID:{7},matched region:{8}-{9}{10}{11}\n'.format(
                    pred[5], 'CDS' if cds == 'CDS' else 'pseudogene', start, stop, pred[11], 
                    '{0}_{1}_{2}'.format(prefix, gid, pid), pred[0], allele_id, s, e, 
                    '' if pred[0] == pred[3] else ',structure variant group:' + part[3], 
                    '' if cds == 'CDS' else ';pseudogene=' + cds, 
                    '' if len(old_tag) == 0 else 'locus_tag={0};'.format(','.join(old_tag)), 
                ))
    allele_file.close()
    return



def add_args(a) :
    import argparse
    parser = argparse.ArgumentParser(description='''
EToKi.py ortho 
(1) Retieves genes and genomic sequences from GFF files and FASTA files.
(2) Groups genes into clusters using vsearch.
(3) Maps gene clusters back to genomes. 
(4) Filters paralogous cluster alignments.
(5) identify a set of most probable non-overlapping orthologs.
(6) Re-annotate genomes using the new set of orthologs. 
''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('GFFs', metavar='N', nargs='+')
    parser.add_argument('-g', '--genes', help='Comma delimited files for additional genes. ', default='')
    parser.add_argument('-G', '--genomes', help='Comma delimited files for additional genomes. ', default='')
    parser.add_argument('-P', '--priority', help='Comma delimited filenames that contain highly confident genes. ', default='')
    
    parser.add_argument('-p', '--prefix', help='prefix for the outputs. Default: EnOrth', default='EnOrth')
    parser.add_argument('-f', '--prefilter', dest='precluster', help='Method for MCL clustering. pairwise [default] or reference', default='pairwise')
    parser.add_argument('-o', '--orthology', help='Method to define orthologous groups. nj [default], ml or rapid (for extremely large datasets)', default='nj')

    parser.add_argument('-t', '--n_thread', help='Number of threads. Default: 30', default=30, type=int)
    parser.add_argument('--source_premature', help='Number of internal stop codon in the reference genes. Default: 0.', default=0., type=float)
    parser.add_argument('--min_cds', help='Minimum length of a reference CDS. Default: 120.', default=120., type=float)

    parser.add_argument('--clust_difference', help='max distance in vsearch clusters. Default: 0.05', default=0.05, type=float)
    parser.add_argument('--clust_match_prop', help='minimum matches in vsearch clusters. Default: 0.9', default=0.9, type=float)

    parser.add_argument('--match_identity', help='minimum identities in BLAST search. Default: 0.6', default=0.6, type=float)
    parser.add_argument('--match_prop', help='minimum match proportion for short genes in BLAST search. Default: 0.7', default=0.7, type=float)
    parser.add_argument('--match_len', help='minimum match proportion for short genes in BLAST search. Default: 300', default=300., type=float)
    parser.add_argument('--match_prop2', help='minimum match proportion for long genes in BLAST search. Default: 0.5', default=0.5, type=float)
    parser.add_argument('--match_len2', help='minimum match proportion for long genes in BLAST search. Default: 500', default=500., type=float)
    parser.add_argument('--match_frag_prop', help='Min proportion of each fragment for fragmented matches. Default: 0.4', default=0.4, type=float)
    parser.add_argument('--match_frag_len', help='Min length of each fragment for fragmented matches. Default: 60', default=60., type=float)
    
    parser.add_argument('--synteny_gap', help='Consider two fragmented matches within N bases as a synteny block. Default: 150', default=150., type=float)
    parser.add_argument('--synteny_diff', help='. Default: 1.5', default=1.5, type=float)
    parser.add_argument('--synteny_ovl_prop', help='Max proportion of overlaps between two fragments in a synteny block. Default: 0.7', default=0.7, type=float)
    parser.add_argument('--synteny_ovl_len', help='Max length of overlaps between two fragments in a synteny block. Default: 300', default=300, type=float)
    parser.add_argument('--edge_rescure', help='Consider fragments that are within N bases of contig edges as part of a synteny block. Default: 150', default=150., type=float)

    parser.add_argument('--mutation_variation', help='Relative variation level in an ortholog group. Default: 2.', default=2., type=float)

    return parser.parse_args(a)

params = dict(
    ml = '{fasttree} {0} -nt -gtr -pseudo', 
    nj = '{rapidnj} -i fa -t d {0}', 
)

pool = None
def ortho(args) :
    global params
    params.update(add_args(args).__dict__)

    global pool
    pool = Pool(10)
    genomes, genes = readGFF(params['GFFs'])
    genes, genomes = addGenes(genes, params['genes']), addGenomes(genomes, params['genomes'])
    if 'old_prediction' not in params :
        params['old_prediction'] = params['prefix']+'.old_prediction.dump'
        old_predictions = {}
        for n, g in genes.iteritems() :
            if g[1] != '' :
                if g[1] not in old_predictions :
                    old_predictions[g[1]] = []
                old_predictions[g[1]].append([n, g[2], g[3], g[4]])
        for g in old_predictions.itervalues() :
            g.sort()
        marshal.dump(old_predictions, open(params['old_prediction'], 'wb'))
        del old_predictions, n, g
    
    if 'prediction' not in params :
        first_classes = load_priority( params.get('priority', ''), genes )

        if 'vclust' not in params :
            params['genes'] = writeGenes('{0}.genes'.format(params['prefix']), genes, first_classes)
            del genes
            params['uc'], params['vclust'] = get_vclust(params['prefix'], params['genes'])
            
        if 'homolog' not in params :
            if 'self_bsn' not in params :
                params['self_bsn'] = get_self_bsn(params['prefix'], params['vclust'])
            if 'map_bsn' not in params :
                params['map_bsn'] = get_map_bsn(params['prefix'], params['vclust'], genomes)
    
            matches = load_matches(params['map_bsn'], genomes)
            groups = combine_matches(matches)
            conflicts = get_conflicts(matches, get_similar_pairs(params['self_bsn']))
            marshal.dump([matches, groups, conflicts], open(params['prefix'] + '.homolog.dump', 'wb'))
        else :
            matches, groups, conflicts = marshal.load(open(params['homolog'], 'rb'))
            conflicts = { int(id): cmp for id, cmp in conflicts.iteritems() }
        matches.sort(key=lambda x:x[11])
        global_differences = {}
        if params['precluster'] in ('pairwise', 'reference') :
            if params['precluster'] == 'pairwise' :
                if 'global' not in params :
                    params['global'] = params['prefix'] + '.global.dump'
                    global_differences = global_difference(matches, groups)
                    marshal.dump({ pair:'{0:.40f}'.format(value) for pair, value in global_differences.iteritems()}, open(params['global'], 'wb'))
                global_differences = {pair:max(float(value), 0.015) for pair, value in marshal.load(open(params['global'], 'rb')).iteritems()}
                global_differences.update({(pair[1], pair[0]):value  for pair, value in global_differences.iteritems() })
            groups = perform_mcl(groups, matches, params['precluster'], global_differences)
        
        params['prediction'] = filt_genes(params['prefix'], matches, groups, global_differences, conflicts, first_classes)
        genes = readFasta(params['vclust'])
    else :
        genes = {n:s[-1] for n,s in genes.iteritems() }
    pool.close()
    old_predictions = marshal.load( open(params['old_prediction'], 'rb')) if 'old_prediction' in params else {}
    
    write_output(params['prefix'], params['prediction'], genomes, genes, old_predictions)
    
if __name__ == '__main__' :
    ortho(sys.argv[1:])

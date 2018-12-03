import sys, shlex, re, math, ete3, hashlib, tempfile, scipy, marshal
from operator import itemgetter
import subprocess, os, numpy as np
from collections import Counter
from time import gmtime, strftime
from multiprocessing import Pool
try: 
    xrange(3)
except :
    xrange = range

def iter_readGFF(fname) :
    seq, cds = {}, {}
    names = {}
    with uopen(fname) as fin :
        sequenceMode = False
        for line in fin :
            if line.startswith('#') : 
                continue
            elif line.startswith('>') :
                sequenceMode = True
                name = line[1:].strip().split()[0]
                assert name not in seq, logger('Error: duplicated sequence name {0}'.format(name))
                seq[name] = [fname, []]
            elif sequenceMode :
                seq[name][1].extend(line.strip().split())
            else :
                part = line.strip().split('\t')
                if len(part) > 2 :
                    name = re.findall(r'locus_tag=([^;]+)', part[8])
                    if len(name) == 0 :
                        parent = re.findall(r'Parent=([^;]+)', part[8])
                        if len(parent) and parent[0] in names :
                            name = names[parent[0]]
                    if len(name) == 0 :
                        name = re.findall(r'Name=([^;]+)', part[8])
                    if len(name) == 0 :
                        name = re.findall(r'ID=([^;]+)', part[8])

                    if part[2] == 'CDS' :
                        assert len(name) > 0, logger('Error: CDS has no name. {0}'.format(line))
                        #          source_file, seqName, Start,       End,      Direction, hash, Sequences
                        cds[name[0]] = [fname, part[0], int(part[3]), int(part[4]), part[6], 0, '']
                    else :
                        ids = re.findall(r'ID=([^;]+)', part[8])
                        if len(ids) :
                            names[ids[0]] = name

    for n in seq :
        seq[n][1] = ''.join(seq[n][1]).upper()
    for n in cds :
        c = cds[n]
        try:
            c[6]= seq[c[1]][1][(c[2]-1) : c[3]]
            if c[4] == '-' :
                c[6] = rc(c[6])
            if not checkCDS(n, c[6]) :
                c[6] = ''
            else :
                c[5] = int(hashlib.sha1(c[6].encode('utf-8')).hexdigest(), 16)
        except :
            pass

    return seq, cds    

def readGFF(fnames) :
    if not isinstance(fnames, list) : fnames = [fnames]
    combo = pool.map(iter_readGFF, fnames)

    seq, cds = {}, {}
    for fname, (ss, cc) in zip(fnames, combo) :
        fprefix = fname.split('.')[0]
        for n in ss :
            seq['{0}:{1}'.format(fprefix, n)] = ss[n][:]

        for n in cc :
            c = cc[n]
            if len(c[6]) > 0 :
                c[1] = '{0}:{1}'.format(fprefix, c[1])
                cds['{0}:{1}'.format(fprefix, n)] = c[:]
    return seq, cds

def get_similar_pairs(prefix, clust) :
    def get_self_bsn(prefix, clust) :
        subprocess.Popen('{0} -dbtype nucl -in {1}'.format(params['formatdb'], clust).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        subprocess.Popen('{0} -num_threads {3} -query {1} -db {1} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen qseq sseq" -max_target_seqs 2000 -task blastn -reward 2 -penalty -3 -evalue 0.1|pigz > {2}.self.bsn.gz'.format(params['blastn'], clust, prefix, int(params['n_thread'])), shell=True ).communicate()
        return '{0}.self.bsn.gz'.format(prefix)
    self_bsn = get_self_bsn(prefix, clust)
    def get_similar(bsn, ortho_pairs) :
        key = tuple(sorted([bsn[0][0], bsn[0][1]]))
        if key in ortho_pairs :
            return
        matched_aa = {}
        len_aa = int(int(bsn[0][12])/3)
        for part in bsn :
            s_i, e_i, s_j, e_j = [ int(x) for x in part[6:10] ]
            frame_i, frame_j = s_i % 3, s_j % 3
            if part[14].find('-') < 0 and part[15].find('-') < 0 :
                if frame_i == frame_j or params['incompleteCDS'] :
                    matched_aa.update({site:1 for site in xrange(int((s_i+1)/3 *3+1), e_i, 3)})
            else :
                s_i, s_j = s_i - 1, s_j - 1
                for b_i, b_j in zip(part[14], part[15]) :
                    if b_i != '-' : s_i += 1
                    if b_j != '-' : s_j += 1
                    
                    if '-' not in {b_i, b_j} and s_i % 3 == 1 and (s_j % 3 == 1 or params['incompleteCDS']) :
                        matched_aa[s_i] = 1
                    elif b_j != '-' :
                        s_j += 1
            if len(matched_aa)*3 >= min(params['match_len2'], params['match_len']) or len(matched_aa) >= (min(params['match_prop'], params['match_prop2'])-0.1) * len_aa :
                ortho_pairs[key] = 1
                return
    
    presence = {}
    ortho_pairs = {}
    save = []
    with uopen(self_bsn) as fin :
        for id, line in enumerate(fin) :
            if id % 100000 == 0 :
                logger('Loading self-BLASTn results. {0}'.format(id))
    
            part = line.strip().split()
            if part[0] not in presence :
                presence[part[0]] = 1
            elif presence[part[0]] == 0 :
                continue
            iden, qs, qe, ss, se, ql, sl = float(part[2]), float(part[6]), float(part[7]), float(part[8]), float(part[9]), float(part[12]), float(part[13])
            if part[1] not in presence :
                if iden >= 100.*params['clust_identity'] and (qe-qs+1)*(se-ss+1) >= params['clust_match_prop']*ql*sl :
                    presence[part[1]] = 0
                    continue
            elif presence[part[1]] == 0 :
                continue
            if iden < 100.*params['match_identity']-10. or ss > se or int(part[3]) < params['match_frag_len']-10 :
                continue
            if len(save) > 0 and save[0][:2] != part[:2] :
                if len(save) >= 80 :
                    presence[part[0][1]] = 0
                elif part[0][0] != part[0][1] :
                    get_similar(save, ortho_pairs)
                save = []
            save.append(part)
        if len(save) > 0 :
            if len(save) >= 80 :
                presence[part[0][1]] = 0
            elif part[0][0] != part[0][1] :
                get_similar(save, ortho_pairs)
    os.unlink(self_bsn)
    
    toWrite = []
    with uopen(params['clust'], 'r') as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                write= True if name in presence else False
            if write :
                toWrite.append(line)
    with open(params['clust'], 'w') as fout :
        for line in toWrite :
            fout.write(line)
    return ortho_pairs

def load_matches(bsn_file, genomes) :
    matches = []
    with uopen(bsn_file) as bsn :
        for id, line in enumerate(bsn) :
            if id % 100000 == 0 :
                logger('Loading BLASTn results. {0}'.format(id))
            part = line.strip().split()
            cont = part[1]
            strain = genomes[cont][0]
            if int(part[8]) < int(part[9]) :
                save = [part[0], strain, cont, float(part[2]), int(part[6]), int(part[7]),  int(part[8]),  int(part[9]), int(part[11]), int(part[12]), int(part[13]), id, 0, part[14]]
            else :
                save = [part[0], strain, cont, float(part[2]), int(part[6]), int(part[7]), -int(part[8]), -int(part[9]), int(part[11]), int(part[12]), int(part[13]), id, 0, part[14]]
            matches.append(save)
    return matches

def combine_matches(save) :
    groups = {}
    for m in save :
        if m[0] not in groups :
            groups[m[0]] = [m]
        else :
            groups[m[0]].append(m)
    return dict(list(map(combine_matches2, groups.values())))

def combine_matches2(save) :
    gene = save[0][0]
    ref_len = save[0][9]
    min_match_len = [min(params['clust_match_prop']*ref_len, max(params['match_prop'] * ref_len, params['match_len']), max(params['match_prop2'] * ref_len, params['match_len2'])), 
                     min(params['clust_match_prop']*ref_len, max(params['match_prop2'] * ref_len, params['match_len2'])) ]
    
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

    save.sort(key=itemgetter(1,2,6))
    groups = {}
    prev, edges = save[0][1], [[], []]
    nSave = len(save)
    
    for id, m1 in enumerate(save) :
        if m1[1] not in groups :
            groups[m1[1]] = []
        match_len = m1[5]-m1[4]+1
        if match_len >= min_match_len[0] :  # minimum length requirement for a non-break gene
            groups[m1[1]].append([match_len, m1[8], m1[3], m1[11]])
            m1[12] = 1
        if prev != m1[1] :
            groups[prev].extend(resolve_edges(edges))       # find genes that were fragmented due to assembly gaps
            prev = m1[1]
            edges =[[], []]
        if m1[4] > 20 and ((m1[6] > 0 and m1[6] - 1 <= params['edge_rescue']) or (m1[6] < 0 and m1[10] + m1[6] <= params['edge_rescue'])) :   # any hit within the last 150 bps to either end of a scaffold is a potential fragmented gene
            edges[1].append(m1)
        if m1[5] < ref_len-20 :
            if (m1[6] > 0 and m1[10]-m1[7] <= params['edge_rescue']) or (m1[6] < 0 and -1-m1[7] <= params['edge_rescue']) :
                edges[0].append(m1)
            for jd in xrange(id+1, nSave) :
                m2 = save[jd]
                match_len2 = m2[5]-m2[4]+1
                if m1[:3] != m2[:3] or (m1[6] < 0 and m2[6] > 0) or m2[6] - m1[7]-1 >= params['synteny_gap'] :    # maximum 300bps between two continuous hits in the same scaffold
                    break
                if m1[7] >= m2[7] or m1[4] > m2[4] or m1[5] > m2[5] or m2[4] - m1[5]-1 >= params['synteny_gap'] or min([match_len, match_len2])*params['synteny_diff'] < max(m2[6]-m1[7]-1, m2[4]-m1[5]-1) :
                    continue
                new_len = m2[5] - m1[4] + 1
                if new_len >= min_match_len[1] :
                    overlap = max(m1[7]-m2[6]+1, m1[5]-m2[4]+1)   ## to improve
                    if overlap >= params['synteny_ovl_len'] or overlap >= params['synteny_ovl_prop']*min(match_len, match_len2) :
                        continue
                    if overlap > 0 :
                        score = m1[8] + m2[8] - overlap * min( float(m1[8])/match_len, float(m2[8])/match_len2 )
                        ident = (m1[3]*match_len + m2[3]*match_len2 - overlap*min(m1[3], m2[3]))/(match_len + match_len2 - overlap)
                    else :
                        score = m1[8] + m2[8]
                        ident = (m1[3]*match_len + m2[3]*match_len2)/(match_len + match_len2)
                    groups[m1[1]].append( [ new_len, score, ident ] + ([m1[11], m2[11]] if m1[3] < m2[3] else [m2[11], m1[11]]) )
                    m1[12], m2[12] = 1, 1
    
    groups[prev].extend(resolve_edges(edges))

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
    return (gene, groups)

def get_conflicts(matches, ortho_pairs) :
    conflicts = {}
    mm = {}
    for m in matches :
        if m[2] not in mm :
            mm[m[2]] = [m]
        else :
            mm[m[2]].append(m)
    mm = sorted(list(mm.values()), key=lambda x:len(x), reverse=True)
    for iter in xrange(0, len(mm), 100) :
        logger('Working out conflicts in hits. {0}/{1}'.format(iter, len(mm)))
        conflicts2 = pool.map(get_conflicts2, [[m, ortho_pairs] for m in mm[iter:iter+100]])
        for cfl in conflicts2 :
            conflicts.update(cfl)
    return conflicts

def get_conflicts2(data) :
    matches, ortho_pairs = data
    conflicts = {}
    matches.sort(key=lambda x:x[6] if x[6]>0 else -x[7])
    matches.sort(key=itemgetter(1))
    len_match = len(matches)
    for id, m1 in enumerate(matches) :
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
                    key = tuple(sorted([m1[0], m2[0]]))
                    if key in ortho_pairs :
                        conflicts[m1[11]][m2[11]] = 2
                        conflicts[m2[11]][m1[11]] = 2
                    else :
                        conflicts[m1[11]][m2[11]] = 1
                        conflicts[m2[11]][m1[11]] = 1
    for g1 in conflicts :
        conflicts[g1] = ' '.join(['{0}={1}'.format(x, y) for x, y in conflicts[g1].items()])
    return conflicts

def get_seq(m) :
    seq = np.zeros(m[3][9], dtype='c')
    seq.fill(b'-')
    for mat in m[3:] :
        seq[mat[4]-1:mat[5]] = mat[13]
    return seq
def compare_seq(seqs) :
    if len(seqs) < 1 :
        return {}
    genome = list(seqs.keys())
    seq = np.array([seqs[g] for g in genome])
    seq[seq == b'N'] = b'-'
    diff = np.zeros(shape=[seq.shape[0], seq.shape[0], 2])
    res = {}
    for id, s in enumerate(seq) :
        comparable = (s >= b'A') * (seq[(id+1):] >= b'A')
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

def global_difference2(mat) :
    seqs = { genome:get_seq(m[0]) for genome, m in mat.items() if len(m) == 1 }
    diff = compare_seq(seqs)
    return diff

def global_difference(groups, counts=3000) :
    global_differences = {}
    grp_order = [ go[0] for go in sorted([ (grp, len([m for genome, m in mat.items() if len(m) == 1])) for grp, mat in groups.items() ], key=itemgetter(1), reverse=True)[:counts] if go[1] > 0]
    for iter in xrange(0, len(grp_order), 100) :
        logger('finding ANIs between genomes. {0}/{1}'.format(iter, len(grp_order)))
        diffs = pool.map(global_difference2, [groups[g] for g in grp_order[iter:iter+100]])
        for diff in diffs :
            for pair, (mut, aln) in diff.items() :
                if pair not in global_differences :
                    global_differences[pair] = []
                if aln :
                    global_differences[pair].append(max(mut, .5)/aln)
    for pair, info in global_differences.items() :
        diff = np.log(info)
        mean_diff = max(np.mean(diff), -4.605)
        sigma = min(max(np.sqrt(np.mean((diff - mean_diff)**2))*3, 0.693), 1.386)
        global_differences[pair] = (np.exp(mean_diff), np.exp(sigma))
    return global_differences

def filt_per_group(data) :
    seqs, grp, mat, tags, global_differences = data
    diff = compare_seq(seqs)
    incompatible, distances = {}, {}
    for (r1, r2), (mut, aln) in diff.items() :
        if aln > 0 :
            gd = global_differences.get((r1[0], r2[0]), (0.01, 4.))
            distances[(r1, r2)] = distances[(r2, r1)] = max(0., 1-(aln - mut)/aln/(1 - gd[0]) )
            difference = mut/aln/gd[0]/gd[1]
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
        for n, s in sorted(seqs.items(), key=lambda x:mat[x[0][0]][int(x[0][1])][1], reverse=True) :
            novel = 1
            for tag, g in groups :
                key = (tag, n) if tag < n else (n, tag)
                if diff[key][0] <= 0.6*(1.0-params['clust_identity'])*diff[key][1] and len(seqs[tag]) + 10 >= len(seqs[n]) :
                    g.append(n)
                    group_tag['{0}__{1}'.format(*n)] = '{0}__{1}'.format(*tag)
                    novel = 0
                    break
            if novel :
                tags['{0}__{1}'.format(*n)] = ''.join(s.astype(str).tolist())
                groups.append([n, [n]])
                group_tag['{0}__{1}'.format(*n)] = '{0}__{1}'.format(*n)
                
        if len(tag) < 20000 :
            groups = {'{0}__{1}'.format(*g): grp for g, grp in groups}
            ic2 = {}
            for r1, rr in incompatible.items() :
                if r1 in group_tag :
                    t1 = group_tag[r1]
                    if t1 not in ic2 :
                        ic2[t1] = {}
                    for r2 in rr :
                        if r2 in group_tag and t1 != group_tag[r2] :
                            ic2[t1][group_tag[r2]] = 1
                    if not len(ic2[t1]) :
                        ic2.pop(t1, None)
            if not len(ic2) :
                return mat
            incompatible = ic2
    
            tmpFile = tempfile.NamedTemporaryFile(dir='.', delete=False)
            for n, s in tags.items() :
                tmpFile.write('>{0}\n{1}\n'.format(n, s).encode('utf-8'))
            tmpFile.close()
            cmd = params[params['orthology']].format(tmpFile.name, **params) if len(tags) < 500 else params['nj'].format(tmpFile.name, **params)
            phy_run = subprocess.Popen(shlex.split(cmd), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            gene_phy = ete3.Tree(phy_run.communicate()[0].replace("'", ''))
            os.unlink(tmpFile.name)
            node = gene_phy.get_midpoint_outgroup()
            if node is not None :
                gene_phy.set_outgroup(node)
            s_name = set(tags.keys()) - set(['{0}__REF'.format(grp.split(':', 1)[-1])])
    
            for ite in xrange(3000) :
                gene_phy.ic, gene_phy.dist = {}, 0.
                rdist = sum([c.dist for c in gene_phy.get_children()])
                for c in gene_phy.get_children() :
                    c.dist = rdist
                for node in gene_phy.iter_descendants('postorder') :
                    if node.is_leaf() :
                        node.ic = {tuple(sorted([node.name, n2])):1 for n2 in incompatible.get(node.name, {}) }
                    else :
                        node.ic = {}
                        for c in node.get_children() :
                            for x in c.ic :
                                if x in node.ic :
                                    node.ic.pop(x)
                                else :
                                    node.ic[x] = 1
                cut_node = max([[len(n.ic), n.dist, n] for n in gene_phy.iter_descendants('postorder')], key=lambda x:(x[0], x[1]))
                if cut_node[0] > 0 :
                    cut_node = cut_node[2]
                    prev_node = cut_node.up
                    cut_node.detach()
                    if grp.split(':')[-1] + '__REF' in cut_node.get_leaf_names() :
                        gene_phy = cut_node
                    elif prev_node.is_root() :
                        gene_phy = gene_phy.get_children()[0]
                    else :
                        prev_node.delete(preserve_branch_length=True)
                    
                    tips = set(gene_phy.get_leaf_names())
                    for r1 in list(incompatible.keys()) :
                        if r1 not in tips :
                            rr = incompatible.pop(r1, None)
                            for r2 in rr :
                                incompatible.get(r2, {}).pop(r1, None)
                    for r1 in list(incompatible.keys()) :
                        if len(incompatible[r1]) == 0 :
                            incompatible.pop(r1, None)
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
        scores = { gene:sum( [ group[0][1]*group[0][0] for genome, group in groups[gene].items() if len(group) > 0 ] ) 
                   for gene, score in ranking.items() if score == min_rank }
    else :
        min_rank = -1
        scores = {}
    genes = [[gene, score, min_rank] for gene, score in sorted(sorted(scores.items()), key=itemgetter(1), reverse=True)[:cnt] if score > 0]
    if len(genes) <= 0 :
        for gene in scores :
            groups.pop(gene, None)
        return []
    else :
        for gene, score in scores.items() :
            if score <= 0 :
                groups.pop(gene, None)
    return genes

def filt_genes(prefix, matches, groups, global_differences, conflicts, first_classes = None) :
    clust_ref = readFasta(params['clust'])
    
    pangenes, used, results, run = {}, {}, {}, {}
    group_id = 0
    with open('{0}.Prediction'.format(prefix), 'w') as fout :
        while len(groups) > 0 :
            genes = get_gene(groups, first_classes, cnt=50)
            if len(genes) <= 0 :
                continue
            to_run, to_run_id, min_score, min_rank = [], [], genes[-1][1], genes[0][2]
            genes = {gene:score for gene, score, min_rank in genes}
            if params['orthology'] in ('ml', 'nj') :
                for gene, score in genes.items() :
                    if gene not in run :
                        mat = groups[gene]
                        for genome, m in mat.items() :
                            m.sort(key=itemgetter(1), reverse=True)
                        region_score = [region[1]/m[0][1] for genome, m in mat.items() for region in m]
                        if len(region_score) > len(mat) * 2 :
                            used2 = {}
                            for genome, m in mat.items() :
                                remove = []
                                for id, mm in enumerate(m) :
                                    for mmm in mm[3:] :
                                        if mmm[11] in used2 :
                                            remove.append(id)
                                            break
                                    if id not in remove :
                                        for mmm in mm[3:] :
                                            used2[mmm[11]] = 1
                                for id in reversed(remove) :
                                    m.pop(id)
                            region_score = [region[1]/m[0][1] for genome, m in mat.items() for region in m]
                        if len(region_score) > len(mat) * 3 and len(region_score) > 500 :
                            cut = sorted(region_score, reverse=True)[len(mat)*3]
                            if cut >= params['clust_identity'] :
                                cut = min(sorted(region_score, reverse=True)[len(mat)*5] if len(region_score) > len(mat) * 5 else params['clust_identity'], 1.0 - 0.6*(1.0-params['clust_identity']))
                            for genome, m in mat.items() :
                                m[:] = [ region for region in m if region[1]/m[0][1] >= cut ]
    
                        to_run.append([{ (genome, str(i)):get_seq(region) for genome, m in mat.items() for i,region in enumerate(m) }, \
                                                   gene, mat, {'{0}__REF'.format(gene.split(':')[-1]):clust_ref[gene]}, global_differences])
                        to_run_id.append(gene)
                working_groups = pool.map(filt_per_group, to_run)
                #working_groups = [filt_per_group(d) for d in to_run]
                for gene, working_group in zip(to_run_id, working_groups) :
                    for g, m in working_group.items() :
                        m.sort(key=lambda x:x[0]*x[1], reverse=True)
                    groups[gene] = working_group
                    run[gene] = 1
            else :
                for gene, score in genes.items() :
                    if gene not in run :
                        mat = groups[gene]
                        for genome, m in mat.items() :
                            m.sort(key=itemgetter(1), reverse=True)
                        cut = params['clust_identity']
                        for genome, m in mat.items() :
                            if len(m) > 1 :
                                m[:] = [ region for region in m if region[1]/m[0][1] >= cut and region[2]/m[0][2] >= cut ]

                        used2 = {}
                        for genome, m in mat.items() :
                            remove = []
                            for id, mm in enumerate(m) :
                                for mmm in mm[3:] :
                                    if mmm[11] in used2 :
                                        remove.append(id)
                                        break
                                if id not in remove :
                                    for mmm in mm[3:] :
                                        used2[mmm[11]] = 1
                            for id in reversed(remove) :
                                m.pop(id)



            while len(genes) :
                score, gene = max([[sum([group[0][1]*group[0][0] for genome, group in groups[g[0]].items() if len(group) > 0 ]), g[0]] for g in genes.items()])
                if score < min_score :
                    break
                working_group = groups.pop(gene)
                genes.pop(gene)

                paralog, reg_cnt = 0, 0
                supergroup, todel = {}, {}
                for genome, region in working_group.items() :
                    for reg in region :
                        reg_cnt += 1
                        skip = 0
                        for mat in reg[3:] :
                            gid = mat[11]
                            conflict = used.get(gid, None)
                            if conflict is not None :
                                skip = 1
                                if not isinstance(conflict, int) :
                                    supergroup[conflict] = supergroup.get(conflict, 0) + 1
                                elif conflict > 0 and reg[2] > 100.0*params['clust_identity'] and len(reg) <= 4 :
                                    paralog = 1
                        if not skip :
                            for mat in reg[3:] :
                                gid = mat[11]
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
                    cnt = min(reg_cnt, len([1 for m in results[pangene[0]].values() for mm in m if mm[3][0] == pangene[0]]))
                    if pangene[1] * 2 >= cnt or (pangene[1] >= 3 and pangene[1] * 3 >= cnt) :
                        pangene = pangene[0]
                    else :
                        pangene = gene
                else :
                    pangene = gene

                if pangene not in results :
                    results[pangene] = {}

                logger('{4} / {5}: pan gene "{3}" : "{0}" picked from rank {1} and score {2}'.format(gene, min_rank, score, pangene, len(results), len(groups)+len(results)))

                for genome, region in working_group.items() :
                    if genome not in results[pangene] :
                        results[pangene][genome] = []
                    results[pangene][genome].extend([reg for reg in region if len(reg) > 0])

                for genome, region in sorted(working_group.items()) :
                    for grp in region :
                        if len(grp) > 0 :
                            group_id += 1
                            for g in sorted(grp[3:], key=lambda x:x[4]) :
                                fout.write('{0}\t{1}\t{2}\t{3}\n'.format(pangene, min_rank, group_id, '\t'.join([str(x) for x in g[:-2]])))

                todel2 = {}
                for gene, lists in todel.items() :
                    if gene not in groups :
                        continue
                    for genome, ids in lists.items() :
                        for grp in groups[gene].get(genome, []) :
                            for g in grp[3:] :
                                id = g[11]
                                if id in ids :
                                    grp[0] = 0
                                    if ids[id] == 1 and grp[2] >= 100.0*params['clust_identity'] and len(grp) <= 4 :
                                        todel2[gene] = 1
                        groups[gene][genome].sort(key=lambda x:x[0]*x[1], reverse=True)
                        for id, cat in ids.items() :
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
    first_classes = { n:[ file_priority.get(g[0], unassign_id), -len(g[6]), g[5], n ] for n, g in genes.items() }

    return first_classes

def get_clust(prefix, genes) :
    subprocess.Popen('{0} --dbmask soft --usersort --cluster_smallmem {4} --uc {1}.gene.uc --centroids {1}.gene.reference --id {2} --query_cov {3} --target_cov {3} --threads {5}'.format(params['vsearch'], prefix, params['clust_identity'], params['clust_match_prop'], genes, params['n_thread']).split()).wait()
    #os.unlink('{0}.gene'.format(prefix))
    return '{0}.gene.uc'.format(prefix), '{0}.gene.reference'.format(prefix)

def get_mmseq(prefix, genes, cnt) :
    import shutil, glob

    nRef = 999999999999999
    subprocess.Popen('{0} createdb {2} {1}.db -v 0'.format(params['mmseqs'], prefix, genes).split()).communicate()
    for ite in xrange(5) :
        if not os.path.isdir(prefix + '.tmp') :
            os.makedirs(prefix + '.tmp')
        if os.path.isfile(prefix + '.lc') :
            for fname in glob.glob(prefix+'.lc*') :
                try : 
                    os.unlink(fname)
                except :
                    pass
        subprocess.Popen('{0} linclust {1}.db {1}.lc {1}.tmp --min-seq-id {2} -c {3} --threads {4} -v 0'.format(params['mmseqs'], prefix, params['clust_identity'], params['clust_match_prop'], params['n_thread']).split(), stdout=subprocess.PIPE).communicate()
        subprocess.Popen('{0} createtsv {1}.db {1}.db {1}.lc {1}.tab'.format(params['mmseqs'], prefix).split(), stdout = subprocess.PIPE).communicate()
        grps = {}
        with open('{0}.tab'.format(prefix)) as fin :
            for line in fin :
                part = line.strip().split()
                grps[part[1]] = part[0]
        with open(genes) as fin, open('{0}.gene.reference'.format(prefix), 'w') as fout :
            write, used_grps = False, {None:1}
            
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    grp = grps.get(name, None)
                    write = False if grp in used_grps else True
                    if write :
                        used_grps[grp] = 1
                if write :
                    fout.write(line)
        shutil.rmtree(prefix + '.tmp')
        for fn in (prefix+'.db*', prefix+'.lc*') :
            for fname in glob.glob(fn) :
                try :
                    os.unlink(fname)
                except :
                    pass
        if nRef <= len(used_grps) :
            break
        nRef = len(used_grps)
        if ite < 4 :
            subprocess.Popen('{0} createdb {1}.gene.reference {1}.db -v 0'.format(params['mmseqs'], prefix).split()).communicate()
    return '{0}.gene.reference'.format(prefix)
    
def writeGenomes(fname, seqs) :
    with open(fname, 'w') as fout :
        for g, s in seqs.items() :
            fout.write('>{0} {1}\n{2}\n'.format(g, s[0], s[-1]))
    return fname

def iter_map_bsn(data) :
    prefix, clust, id, taxon, seq = data
    gfile, bsn = '{0}.{1}.genome'.format(prefix, id), '{0}.{1}.bsn.gz'.format(prefix, id)
    ortho_bsn = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ortho_bsn.py')
    with open(gfile, 'w') as fout :
        fout.write(''.join(seq))
    subprocess.Popen('{0} -dbtype nucl -in {1}'.format(params['formatdb'], gfile).split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    subprocess.Popen('{0} -num_threads {4} -query {1} -db {2} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen qseq sseq" -max_target_seqs 2000 -task blastn -reward 2 -penalty -3 -qcov_hsp_perc 10 -evalue 0.001|{5}|pigz > {3}'.format(\
        params['blastn'], clust, gfile, bsn, int(params['n_thread']/10), \
        '{0} {1} {2} {3} {4}'.format(sys.executable, ortho_bsn, params['match_identity'], params['match_frag_len'], params['match_frag_prop'])), shell=True).communicate()
    for fn in (gfile, gfile+'.nin', gfile+'.nsq', gfile+'.nhr') :
        os.unlink(fn)
    return bsn

def get_map_bsn(prefix, clust, genomes) :
    if len(genomes) == 0 :
        sys.exit(1)
    taxa = {}
    for g, s in genomes.items() :
        if s[0] not in taxa : taxa[s[0]] = []
        taxa[s[0]].append('>{0}\n{1}\n'.format(g, s[1]))
    
    bsns = pool.map(iter_map_bsn, [(prefix, clust, id, taxon, seq) for id, (taxon, seq) in enumerate(taxa.items())])
    #bsns = list(map(iter_map_bsn, [(prefix, clust, id, taxon, seq) for id, (taxon, seq) in enumerate(taxa.items())]))
    subprocess.Popen('cat {0} > {1}.map.bsn.gz'.format(' '.join(bsns), prefix), shell=True).communicate()
    for bsn in bsns :
        os.unlink(bsn)
    return '{0}.map.bsn.gz'.format(prefix)

def checkCDS(n, s) :
    if len(s) < params['min_cds'] :
        logger('{0} is too short'.format(n))
        return False
    if params['incompleteCDS'] :
        return True

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
    if len(aa[:-1].split('X')) > 1 :
        logger('{0} is discarded due to internal stop codons'.format(n))
        return False
    return True
    

def addGenes(genes, gene_file) :
    for gfile in gene_file.split(',') :
        if gfile == '' : continue
        gprefix = gfile.split('.')[0]
        ng = readFasta(gfile)
        for name in ng :
            s = ng[name]
            if checkCDS(name, s) :
                genes['{0}:{1}'.format(gprefix,name)] = [ gfile, '', 0, 0, '+', int(hashlib.sha1(s.encode('utf-8')).hexdigest(), 16), s]
    return genes

def writeGenes(fname, genes, priority) :
    uniques = {}
    cnt = 0
    with open(fname, 'w') as fout :
        for n in sorted(priority.items(), key=itemgetter(1)) :
            s = genes[n[0]][6]
            len_s, hcode = len(s), genes[n[0]][5]
            if len_s not in uniques :
                uniques = { len_s:{hcode:1} }
            elif hcode in uniques[ len_s ] :
                continue
            uniques[ len_s ][ hcode ] = 1
            fout.write( '>{0}\n{1}\n'.format(n[0], s) )
            cnt += 1
    return fname, cnt

def precluster2(data) :
    gene, group, global_differences = data
    matches = sorted([ [sp, id, (1.-g[2]/100.)] for sp, grp in group.items() for id, g in enumerate(grp) ], key=lambda x:x[2])
    nMatch = len(matches)
    inGroup = {0:1}
    newGroup = {sp:[] for sp in group}
    for i1, m1 in enumerate(matches) :
        if i1 in inGroup :
            newGroup[m1[0]].append(group[m1[0]][m1[1]])

            for i2 in xrange(i1+1, nMatch) :
                if i2 in inGroup : 
                    continue
                m2 = matches[i2]
                diff = m2[2] - m1[2]
                if diff > 0 :
                    gd = global_differences.get((m1[0], m2[0]), (0.01, 4.))
                    diff = diff/gd[0]/gd[1]
                if diff < 1 :
                    inGroup[i2] = 1
    return (gene, newGroup)

def precluster(groups, global_differences) :
    return dict(pool.map(precluster2, [(k, v, global_differences) for k,v in groups.items()]))

def write_output(prefix, prediction, genomes, clust_ref, old_prediction) :
    predictions, alleles = {}, {}
    
    allele_file = open('{0}.allele.fna'.format(prefix), 'w')
    with open(prediction) as fin :
        for line in fin :
            part = line.strip().split()
            if part[0] not in alleles :
                alleles[part[0]] = {clust_ref[part[0]]:1}
                allele_file.write('>{0}_{1}\n{2}\n'.format(part[0], 1, clust_ref[part[0]]))
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
    with open('{0}.EToKi.gff'.format(prefix), 'w') as fout :
        for gid, (g, predict) in enumerate(predictions.items()) :
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
                        s2, e2 = s - min(int(3*((s - 1)/3)), 30), e + min(3*int((pred[13] - e)/3), 300)
                        seq = genomes[pred[5]][1][(s2-1):e2]
                        lp, rp = s - s2, e2 - e
                    else :
                        s2, e2 = s - min(int(3*((s - 1)/3)), 300), e + min(3*int((pred[13] - e)/3), 30)
                        seq = rc(genomes[pred[5]][1][(s2-1):e2])
                        rp, lp = s - s2, e2 - e
                        
                    seq2 = seq[(lp):(len(seq)-rp)]
                    if seq2 not in alleles[pred[0]] :
                        if pred[3] == pred[0] and pred[7] == 1 and pred[8] == pred[12] :
                            alleles[pred[0]][seq2] = len(alleles[pred[0]])+1
                        else :
                            alleles[pred[0]][seq2] = 'LowQ{0}'.format(len(alleles[pred[0]])+1)
                        allele_id = str(alleles[pred[0]][seq2])
                        allele_file.write('>{0}_{1}\n{2}\n'.format(pred[0], allele_id, seq2))
                    else :
                        allele_id = str(alleles[pred[0]][seq2])
                        
                    if len(seq) % 3 > 0 :
                        cds = 'frameshift'
                        start, stop = s, e
                    else :
                        aa_seq = transeq({'n':seq}, transl_table='starts')['n_1']
                        s0, s1 = aa_seq.find('M', int(lp/3), int(lp/3+12)), aa_seq.rfind('M', 0, int(lp/3))
                        start = s0 if s0 >= 0 else s1
                        if start < 0 :
                            cds, start = 'nostart', int(lp/3)
                        stop = aa_seq.find('X', start)
                        if stop >= 0 and stop < lp/3+12 :
                            s0 = aa_seq.find('M', stop, int(lp/3+12))
                            if s0 >= 0 :
                                start = s0
                                stop = aa_seq.find('X', start)
                        if stop < 0 :
                            cds, stop = 'nostop', len(aa_seq) - int(rp/3) - 1
                        elif (stop - start + 1)*3 <= 0.8 * pred[12] :
                            cds = 'premature stop:{0:.2f}%'.format((stop - start + 1)*300/pred[12])
                            stop = len(aa_seq) - int(rp/3) - 1
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

                fout.write('{0}\t{1}\tEToKi-ortho\t{2}\t{3}\t.\t{4}\t.\tID={5};{12}inference=ortholog group:{6},allele ID:{7},matched region:{8}-{9}{10}{11}\n'.format(
                    pred[5], 'CDS' if cds == 'CDS' else 'pseudogene', start, stop, pred[11], 
                    '{0}_{1}_{2}'.format(prefix, gid, pid), pred[0], allele_id, s, e, 
                    '' if pred[0] == pred[3] else ',structure variant group:' + pred[3], 
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
(2) Groups genes into clusters using mmseq.
(3) Maps gene clusters back to genomes. 
(4) Filters paralogous cluster alignments.
(5) identify a set of most probable non-overlapping orthologs.
(6) Re-annotate genomes using the new set of orthologs. 
''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('GFFs', metavar='N', help='GFF files containing both annotations and sequences.', nargs='*')
    parser.add_argument('-g', '--genes', help='Comma delimited files for additional genes. ', default='')
    parser.add_argument('-P', '--priority', help='Comma delimited filenames that contain highly confident genes. ', default='')
    
    parser.add_argument('-p', '--prefix', help='prefix for the outputs. Default: EToKi', default='EToKi')
    parser.add_argument('-o', '--orthology', help='Method to define orthologous groups. nj [default], ml or rapid (for extremely large datasets)', default='nj')

    parser.add_argument('-t', '--n_thread', help='Number of threads. Default: 30', default=30, type=int)
    parser.add_argument('--min_cds', help='Minimum length of a reference CDS. Default: 150.', default=150., type=float)

    parser.add_argument('--clust_identity', help='minimum identities in mmseq clusters. Default: 0.9', default=0.9, type=float)
    parser.add_argument('--clust_match_prop', help='minimum matches in mmseq clusters. Default: 0.9', default=0.9, type=float)

    parser.add_argument('--match_identity', help='minimum identities in BLAST search. Default: 0.6', default=0.6, type=float)
    parser.add_argument('--match_prop', help='minimum match proportion for short genes in BLAST search. Default: 0.7', default=0.7, type=float)
    parser.add_argument('--match_len', help='minimum match proportion for short genes in BLAST search. Default: 300', default=300., type=float)
    parser.add_argument('--match_prop2', help='minimum match proportion for long genes in BLAST search. Default: 0.5', default=0.5, type=float)
    parser.add_argument('--match_len2', help='minimum match proportion for long genes in BLAST search. Default: 500', default=500., type=float)
    parser.add_argument('--match_frag_prop', help='Min proportion of each fragment for fragmented matches. Default: 0.4', default=0.4, type=float)
    parser.add_argument('--match_frag_len', help='Min length of each fragment for fragmented matches. Default: 90', default=90., type=float)
    
    parser.add_argument('--synteny_gap', help='Consider two fragmented matches within N bases as a synteny block. Default: 200', default=200., type=float)
    parser.add_argument('--synteny_diff', help='. Default: 1.2', default=1.2, type=float)
    parser.add_argument('--synteny_ovl_prop', help='Max proportion of overlaps between two fragments in a synteny block. Default: 0.7', default=0.7, type=float)
    parser.add_argument('--synteny_ovl_len', help='Max length of overlaps between two fragments in a synteny block. Default: 300', default=300, type=float)
    parser.add_argument('--edge_rescue', help='Consider fragments that are within N bases of contig edges as part of a synteny block. Default: 150', default=150., type=float)

    parser.add_argument('--mutation_variation', help='Relative variation level in an ortholog group. Default: 2.', default=2., type=float)
    parser.add_argument('--incompleteCDS', help='Do not do CDS checking for the reference genes. Default: False.', default=False, action='store_true')
    parser.add_argument('--metagenome', help='Set to metagenome mode. equals to "--incompleteCDS --clust_identity 0.99 --clust_match_prop 0.8 --match_identity 0.98 --prefilter reference"', default=False, action='store_true')

    parser.add_argument('--old_prediction', help='development param', default=None)
    parser.add_argument('--clust', help='development param', default=None)
    parser.add_argument('--map_bsn', help='development param', default=None)
    parser.add_argument('--self_bsn', help='development param', default=None)
    parser.add_argument('--conflicts', help='development param', default=None)
    parser.add_argument('--global', help='development param', default=None)
    parser.add_argument('--prediction', help='development param', default=None)

    params = parser.parse_args(a)
    if params.metagenome :
        params.incompleteCDS = True
        params.clust_identity = 0.99
        params.clust_match_prop = 0.8
        params.match_identity = 0.98
        params.prefilter = 'reference'
    return params

params = dict(
    ml = '{fasttree} {0} -nt -gtr -pseudo', 
    nj = '{rapidnj} -i fa -t d {0}', 
)

pool = None
def ortho(args) :
    global params
    params.update(add_args(args).__dict__)
    params.update(externals)

    global pool
    pool = Pool(params['n_thread'])
    genomes, genes = readGFF(params['GFFs'])
    genes = addGenes(genes, params['genes'])
    if params.get('old_prediction', None) is None :
        params['old_prediction'] = params['prefix']+'.old_prediction.dump'
        old_predictions = {}
        for n, g in genes.items() :
            if g[1] != '' :
                if g[1] not in old_predictions :
                    old_predictions[g[1]] = []
                old_predictions[g[1]].append([n, g[2], g[3], g[4]])
        for g in old_predictions.values() :
            g.sort()
        marshal.dump(old_predictions, open(params['old_prediction'], 'wb'))
        del old_predictions, n, g
    
    if params.get('prediction', None) is None :
        first_classes = load_priority( params.get('priority', ''), genes )

        if params.get('clust', None) is None :
            params['genes'], cnt = writeGenes('{0}.genes'.format(params['prefix']), genes, first_classes)
            del genes
            params['clust'] = get_mmseq(params['prefix'], params['genes'], cnt)
            
        if params.get('self_bsn', None) is None :
            orthoGroup = get_similar_pairs(params['prefix'], params['clust'])
            marshal.dump(orthoGroup, open(params['prefix']+'.self.bsn.dump', 'wb'))
        else :
            orthoGroup = marshal.load(open(params['self_bsn'], 'rb'))

        if params.get('map_bsn', None) is None :
            params['map_bsn'] = get_map_bsn(params['prefix'], params['clust'], genomes)
        matches = load_matches(params['map_bsn'], genomes)
        groups = combine_matches(matches)
        if params.get('conflicts', None) is None :
            conflicts = get_conflicts(matches, orthoGroup)
            marshal.dump(conflicts, open(params['prefix'] + '.conflicts.dump', 'wb'))
        else :
            conflicts = marshal.load(open(params['conflicts'], 'rb'))
            conflicts = { int(id): cmp for id, cmp in conflicts.items() }

        matches.sort(key=lambda x:x[11])
        for gene, group in groups.items() :
            for genome, grp in group.items() :
                for g in grp :
                    for id in xrange(3, len(g)) :
                        g[id] = matches[g[id]]

        if params.get('global', None) is None :
            params['global'] = params['prefix'] + '.global.dump'
            global_differences = global_difference(groups, 3000)
            marshal.dump({pair:['{:.20f}'.format(v) for v in value] for pair, value in global_differences.items()}, open(params['global'], 'wb'))
        global_differences = {pair:[float(v) for v in value] for pair, value in marshal.load(open(params['global'], 'rb')).items()}
        global_differences.update({(pair[1], pair[0]):value for pair, value in global_differences.items() })

        groups = precluster(groups, global_differences)
        params['prediction'] = filt_genes(params['prefix'], matches, groups, global_differences, conflicts, first_classes)
        genes = readFasta(params['clust'])
    else :
        genes = {n:s[-1] for n,s in genes.items() }
    pool.close()
    old_predictions = marshal.load( open(params['old_prediction'], 'rb')) if 'old_prediction' in params else {}
    
    write_output(params['prefix'], params['prediction'], genomes, genes, old_predictions)
    
if __name__ == '__main__' :
    from configure import externals, logger, rc, transeq, readFasta, uopen    
    ortho(sys.argv[1:])
else :
    from .configure import externals, logger, rc, transeq, readFasta, uopen
    

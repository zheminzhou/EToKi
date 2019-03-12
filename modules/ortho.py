import sys, shlex, re, ete3, tempfile, hashlib, shutil
from operator import itemgetter
import subprocess, os, numpy as np, pandas as pd, numba as nb
from multiprocessing import Pool
try:
    from configure import externals, logger, rc, transeq, readFasta, uopen, xrange, asc2int
    from clust import getClust
    from uberBlast import uberBlast
except :
    from .configure import externals, logger, rc, transeq, readFasta, uopen, xrange, asc2int
    from .clust import getClust
    from .uberBlast import uberBlast
    

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
            c[6] = ''

    return seq, cds    

def readGFF(fnames) :
    if not isinstance(fnames, list) : fnames = [fnames]
    combo = pool.map(iter_readGFF, fnames)

    seq, cds = {}, {}
    for fname, (ss, cc) in zip(fnames, combo) :
        fprefix = os.path.basename(fname).split('.')[0]
        for n in ss :
            seq['{0}:{1}'.format(fprefix, n)] = ss[n][:]

        for n in cc :
            c = cc[n]
            if len(c[6]) > 0 :
                c[1] = '{0}:{1}'.format(fprefix, c[1])
                cds['{0}:{1}'.format(fprefix, n)] = c[:]
    return seq, cds


    
def get_similar_pairs(prefix, clust, priorities, params) :
    def get_similar(bsn, ortho_pairs) :
        key = tuple(sorted([bsn[0][0], bsn[0][1]]))
        if key in ortho_pairs :
            return
        matched_aa = {}
        len_aa = int(int(bsn[0][12])/3)
        for part in bsn :
            s_i, e_i, s_j, e_j = [ int(x) for x in part[6:10] ]
            for s, t in re.findall(r'(\d+)([A-Z])', part[14]) :
                frame_i, frame_j = s_i % 3, s_j % 3
                s = int(s)
                if t == 'M' :
                    if frame_i == frame_j or params['incompleteCDS'] :
                        matched_aa.update({ (s_i+x): 1 for x in xrange( (3 - (frame_i - 1))%3, s )})
                    s_i += s
                    s_j += s
                    if len(matched_aa)*3 >= min(params['match_len2'], params['match_len']) or len(matched_aa) >= (min(params['match_prop'], params['match_prop2'])-0.1) * len_aa :
                        ortho_pairs[key] = 1
                        return
                elif t == 'I' :
                    s_i += s
                else :
                    s_j += s
                    
    self_bsn = uberBlast('-r {0} -q {0} --blastn --ublastSELF -f --min_id {1} --min_cov {2} -t {3} -s 2 -e 3,3'.format(clust, params['match_identity'] - 0.1, params['match_frag_len'], params['n_thread']).split())
    presence, ortho_pairs = {}, {}
    save = []
    for part in self_bsn :
        if part[0] not in presence :
            presence[part[0]] = 1
        elif presence[part[0]] == 0 :
            continue
        iden, qs, qe, ss, se, ql, sl = float(part[2]), float(part[6]), float(part[7]), float(part[8]), float(part[9]), float(part[12]), float(part[13])
        if presence.get(part[1], 1) == 0 or ss > se :
            continue
        if len(save) > 0 and np.any(save[0][:2] != part[:2]) :
            if len(save) >= 50 :
                presence[save[0][1]] = 0
            elif save[0][0] != save[0][1] :
                get_similar(save, ortho_pairs)
            save = []
        save.append(part)
    if len(save) > 0 :
        if len(save) >= 50 :
            presence[save[0][1]] = 0
        elif save[0][0] != save[0][1] :
            get_similar(save, ortho_pairs)
    
    toWrite = []
    with uopen(params['clust'], 'r') as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                write= True if presence.get(int(name), 0) > 0 else False
            if write :
                toWrite.append(line)
    with open(params['clust'], 'w') as fout :
        for line in toWrite :
            fout.write(line)
    return np.array(list(ortho_pairs.keys())).astype(int)

@nb.jit('i8[:,:,:](u1[:,:], i8[:,:,:])', nopython=True)
def compare_seq(seqs, diff) :
    for id in np.arange(seqs.shape[0]) :
        s = seqs[id]
        c = (s > 45) * (seqs[(id+1):] > 45)
        n_comparable = np.sum(c, 1)
        n_comparable[n_comparable < 1] = 1        
        n_diff = np.sum(( s != seqs[(id+1):] ) & c, 1)
        diff[id, id+1:, 0] = n_diff
        diff[id, id+1:, 1] = n_comparable
    return diff

def global_difference2(data) :
    fname, g = data
    _, idx, cnt = np.unique(g.T[1], return_counts=True, return_index=True)
    idx = idx[cnt==1]
    names, seqs = g[idx, 1], np.vstack(g[idx, 4])
    seqs[np.in1d(seqs, [65, 67, 71, 84], invert=True).reshape(seqs.shape)] = 45
    diff = compare_seq(seqs, np.zeros(shape=[seqs.shape[0], seqs.shape[0], 2], dtype=np.int64))
    res = {}
    for i, n1 in enumerate(names) :
        for j in xrange(i+1, len(names)) :
            n2 = names[j]
            if diff[i, j, 1] >= min(params['match_len2'], seqs.shape[1]*params['match_prop2']) :
                res[(n1, n2)] = diff[i, j, :]
    np.save(fname, np.array(list(res.items()), dtype=object))
    return fname

def global_difference(bsn_file, prefix, orthoGroup, counts=3000) :
    groups = np.load(bsn_file)
    genes = []
    for gene, g in groups.items() :
        _, idx, cnt = np.unique(g.T[1], return_counts=True, return_index=True)
        score = (np.sum(cnt==1)-1)*(2**41)-np.sum(g[idx[cnt==1], 2], dtype=int)
        if score > 0 :
            genes.append([score, gene])
    genes = sorted(genes, reverse=True)
    
    og = np.array(list(orthoGroup.keys()))
    grp_order, all_useds = [], set([])
    for score, gene in genes :
        tag = groups[gene][0][0]
        if tag not in all_useds :
            grp_order.append(gene)
            used = og[og.T[0] == tag, 1]
            all_useds |= set(used.tolist())
    genes = grp_order[:counts]

    global_differences = {}
    for iter in xrange(0, len(genes), 100) :
        logger('finding ANIs between genomes. {0}/{1}'.format(iter, len(genes)))
        #diffs = list(map(global_difference2, [['{0}.{1}.npy'.format(prefix, i2), groups[i]] for i2, i in enumerate(genes[iter:iter+100])])
        diffs = pool2.map(global_difference2, [['{0}.{1}.npy'.format(prefix, i2), groups[i]] for i2, i in enumerate(genes[iter:iter+100])])
        for diff_file in diffs :
            diff ={ tuple(a):b for a, b in np.load(diff_file).tolist()}
            os.unlink(diff_file)
            for pair, (mut, aln) in diff.items() :
                if pair not in global_differences :
                    global_differences[pair] = []
                if aln :
                    global_differences[pair].append(max(float(mut), .5)/aln)
    for pair, info in global_differences.items() :
        diff = np.log(info)
        mean_diff = max(np.mean(diff), -4.605)
        sigma = min(max(np.sqrt(np.mean((diff - mean_diff)**2))*3, 0.693), 1.386)
        global_differences[pair] = (np.exp(mean_diff), np.exp(sigma))
    return pd.DataFrame(list(global_differences.items())).values

def filt_per_group(data) :
    mat, ref, global_file = data
    global_differences = dict(np.load(global_file))
    nMat = mat.shape[0]
    seqs = np.vstack([np.vstack(mat.T[4]), np.array(list(ref)).view(asc2int).astype(np.uint8)[np.newaxis, :]])
    seqs[np.in1d(seqs, [65, 67, 71, 84], invert=True).reshape(seqs.shape)] = 45    
    diff = compare_seq(seqs, np.zeros(shape=[seqs.shape[0], seqs.shape[0], 2], dtype=int)).astype(float)
    incompatible, distances = {}, np.zeros(shape=[seqs.shape[0], seqs.shape[0]], dtype=float)
    for i1, m1 in enumerate(mat) :
        for i2 in xrange(i1+1, nMat) :
            m2 = mat[i2]
            mut, aln = diff[i1, i2]
            if aln > 0 :
                gd = global_differences.get(tuple(sorted([m1[1], m2[1]])), (0.01, 4))
                distances[i1, i2] = distances[i2, i1] = max(0., 1-(aln - mut)/aln/(1 - gd[0]) )
                difference = mut/aln/gd[0]/gd[1]
            else :
                distances[i1, i2] = distances[i2, i1] = 0.8
                difference = 1.5
            if difference > 1. :
                incompatible[(i1, i2)] = 1

    if len(incompatible) > 0 :
        groups = []
        for j, m in enumerate(mat) :
            novel = 1
            for g in groups :
                if diff[g[0], j, 0] <= 0.6*(1.0-params['clust_identity'])*diff[g[0], j, 1] :
                    g.append(j)
                    novel = 0
                    break
            if novel :
                groups.append([j])
        group_tag = {gg:g[0] for g in groups for gg in g}
        try :
            tags = {g[0]:mat[g[0]][4].tostring().decode('ascii') for g in groups}
        except :
            tags = {g[0]:mat[g[0]][4].tostring() for g in groups}
            
        tags.update({'REF':ref})

        ic2 = {}
        for i1, i2 in incompatible :
            t1, t2 = group_tag[i1], group_tag[i2]
            if t1 != t2 :
                t1, t2 = str(t1), str(t2)
                if t1 not in ic2 : ic2[t1] = {}
                if t2 not in ic2 : ic2[t2] = {}
                ic2[t1][t2] = ic2[t2][t1] = 1
        incompatible = ic2
        if len(incompatible) == 0 :
            return mat

        for ite in xrange(3) :
            try :
                tmpFile = tempfile.NamedTemporaryFile(dir='.', delete=False)
                for n, s in tags.items() :
                    tmpFile.write('>X{0}\n{1}\n{2}'.format(n, s, '\n'*ite).encode('utf-8'))
                tmpFile.close()
                cmd = params[params['orthology']].format(tmpFile.name, **params) if len(tags) < 500 else params['nj'].format(tmpFile.name, **params)
                phy_run = subprocess.Popen(shlex.split(cmd), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                gene_phy = ete3.Tree(phy_run.communicate()[0].replace("'", ''))
                break
            except :
                if ite == 2 :
                    return mat
            finally:
                os.unlink(tmpFile.name)
        for n in gene_phy.get_leaves() :
            if len(n.name) :
                n.name = n.name[1:] 
        
        node = gene_phy.get_midpoint_outgroup()
        if node is not None :
            gene_phy.set_outgroup(node)

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
                if 'REF' in cut_node.get_leaf_names() :
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
                logger('     Iteration {0}. Remains {1} tips with {2} conflicts.'.format(ite+1, len(gene_phy.get_leaf_names()), len(incompatible)))
                if len(incompatible) == 0 :
                    break
            else :
                break
        if len(gene_phy.get_leaf_names()) < len(tags) :
            groups = {str(g[0]):g for g in groups}
            tips = sorted([ nn for n in gene_phy.get_leaf_names() for nn in groups.get(n, [])])
            mat = mat[tips]
    return mat

def get_gene(groups, first_classes, cnt=1) :
    ranking = {gene:first_classes[gene][0] for gene in groups if gene in first_classes}
    if len(ranking) > 0 :
        min_rank = min(ranking.values())
        scores = { gene:np.sum(groups[gene][np.unique(groups[gene].T[1], return_index=True)[1]].T[2])
                   for gene, score in ranking.items() if score == min_rank }
    else :
        min_rank = -1
        scores = {}
    genes = [[gene, score, min_rank] for gene, score in sorted(sorted(scores.items()), key=itemgetter(1), reverse=True)[:cnt] if score > 0]
    if len(genes) <= 0 :
        for gene in scores :
            groups.pop(gene)
        return []
    elif min(scores.values()) <= 0 :
        for gene, score in scores.items() :
            if score <= 0 :
                groups.pop(gene)
    return genes

def filt_genes(prefix, groups, global_file, conflicts, first_classes = None, encodes = None) :
    encodes = np.array([n for i, n in sorted([[i, n] for n, i in encodes.items()])])
    outPos = np.ones(16, dtype=bool)
    outPos[[3,4,5,10,15]] = False
    
    c2 = { c:{} for c in np.unique(conflicts.T[:2]) }
    for c in conflicts :
        c2[c[0]][c[1]] = c2[c[1]][c[0]] = c[2]
    conflicts = c2
    
    clust_ref = { int(n):s for n, s in readFasta(params['clust']).items()}
    
    for gene, g in groups.items() : 
        g.T[2] *= g.T[3]
        g[:] = g[np.argsort(-g.T[2], kind='mergesort')]
    used, results, run = {}, {}, {}
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
                        _, bestPerGenome, matInGenome = np.unique(mat.T[1], return_index=True, return_inverse=True)
                        region_score = mat.T[2]/mat[bestPerGenome[matInGenome], 2]
                        if region_score.size >= bestPerGenome.size * 2 :
                            used2, kept = set([]), np.ones(mat.shape[0], dtype=bool)
                            for id, m in enumerate(mat) :
                                if m[5] in used2 :
                                    kept[id] = False
                                else :
                                    used2.update(conflicts.get(m[5], {}))
                            mat = mat[kept]
                            _, bestPerGenome, matInGenome = np.unique(mat.T[1], return_index=True, return_inverse=True)
                            region_score = mat.T[2]/mat[bestPerGenome[matInGenome], 2]
                        if region_score.size > bestPerGenome.size * 3 and len(region_score) > 500 :
                            region_score2 = sorted(region_score, reverse=True)
                            cut = region_score2[bestPerGenome.size*3-1]
                            if cut >= params['clust_identity'] :
                                cut = min(region_score2[bestPerGenome.size*5] if len(region_score) > bestPerGenome.size * 5 else params['clust_identity'], 1.0 - 0.6*(1.0-params['clust_identity']))
                            mat = mat[region_score>=cut]
    
                        to_run.append([mat, clust_ref[ mat[0][0] ], global_file])
                        to_run_id.append(gene)
                working_groups = pool2.map(filt_per_group, to_run)
                #working_groups = [filt_per_group(d) for d in to_run]
                for gene, working_group in zip(to_run_id, working_groups) :
                    groups[gene] = working_group
                    run[gene] = 1
            else :
                for gene, score in genes.items() :
                    if gene not in run :
                        mat = groups[gene]
                        _, bestPerGenome, matInGenome = np.unique(mat.T[1], return_index=True, return_inverse=True)
                        region_score = mat.T[2]/mat[bestPerGenome[matInGenome], 2]
                        mat = mat[region_score>=params['clust_identity']]
                        used2, kept = set([]), np.ones(mat.shape[0], dtype=bool)
                        for id, m in enumerate(mat) :
                            if m[5] in used2 :
                                kept[id] = False
                            else :
                                used2.update(conflicts.get(m[5], {}))
                        groups[gene] = mat[kept]

            while len(genes) :
                score, gene = max([[np.sum(groups[gene][np.unique(groups[gene].T[1], return_index=True)[1]].T[2]), gene] for gene in genes])
                if score < min_score :
                    break
                mat = groups.pop(gene, [])
                genes.pop(gene)

                paralog, paralog2 = 0, 0
                supergroup = {}
                used2 = {}
                for m in mat :
                    gid = m[5]
                    conflict = used.get(gid, None) 
                    if conflict is not None :
                        if not isinstance(conflict, int) :
                            superC = results[int(conflict)]
                            supergroup[superC] = supergroup.get(superC, 0) + 1
                        elif conflict >0 :
                            if m[6].shape[0] <= 1 and m[3] >= params['clust_identity'] :
                                paralog = 1
                                break
                            else :
                                paralog2 += 1
                        m[3] = -1
                    else :
                        for g2, gs in conflicts.get(gid, {}).items() :
                            if gs == 1 :
                                if g2 not in used :
                                    used2[g2] = str(m[0])
                            elif gs == 2 :
                                used2[g2] = 1
                            else :
                                used[g2] = 0
                if paralog or paralog2*3 >= mat.shape[0] :
                    continue
                else :
                    used.update(used2)

                pangene = mat[0][0]
                if len(supergroup) :
                    pg, pid = max(supergroup.items(), key=itemgetter(1))
                    if pid*3 >= mat.shape[0] or (pid*5 >= mat.shape[0] and pid>1) :
                        pangene = pg

                results[mat[0][0]] = pangene
                if len(results) % 100 == 0 :
                    logger('{4} / {5}: pan gene "{3}" : "{0}" picked from rank {1} and score {2}'.format(encodes[mat[0][0]], min_rank, score, encodes[pangene], len(results), len(groups)+len(results)))

                for grp in mat[mat.T[3] > 0] :
                    group_id += 1
                    grp[6].T[:2] = encodes[grp[6].T[:2].astype(int)]
                    for g in grp[6] :
                        gg = g[outPos].astype(str).tolist()
                        fout.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(encodes[pangene], min_rank, group_id, encodes[grp[1]], '\t'.join(gg)))
    return '{0}.Prediction'.format(prefix)

def load_priority(priority_list, genes, encodes) :
    file_priority = { encodes.get(fn, ''):id for id, fnames in enumerate(priority_list.split(',')) for fn in fnames.split(':') }
    unassign_id = max(file_priority.values()) + 1
    first_classes = { n:[ file_priority.get(g[0], unassign_id), -len(g[6]), g[5], n ] for n, g in genes.items() }

    return first_classes

    
def writeGenomes(fname, seqs) :
    with open(fname, 'w') as fout :
        for g, s in seqs.items() :
            fout.write('>{0} {1}\n{2}\n'.format(g, s[0], s[-1]))
    return fname

def iter_map_bsn(data) :
    prefix, clust, id, taxon, seq, params = data
    gfile, out_prefix = '{0}.{1}.genome'.format(prefix, id), '{0}.{1}'.format(prefix, id)
    with open(gfile, 'w') as fout :
        for n, s in seq :
            fout.write('>{0}\n{1}\n'.format(n, s) )

    blastab, overlap = uberBlast('-r {0} -q {1} -f -m -O --blastn --ublast --min_id {2} --min_cov {3} -t 2 -s 2 -e 0,3'.format(gfile, clust, params['match_identity']-0.1, params['match_frag_len'] ).split())
    os.unlink(gfile)
    
    groups = []
    groups2 = {}
    ids = np.zeros(np.max(blastab.T[15])+1, dtype=bool)
    for tab in blastab :
        if tab[16][1] >= params['match_identity'] and tab[16][2] >= max(params['match_prop']*tab[12], params['match_len']) and tab[16][2] >= max(params['match_prop2']*tab[12], params['match_len2']) :
            ids[tab[15]] = True
            if len(tab[16]) <= 4 :
                groups.append(tab[:2].tolist() + tab[16][:2] + [None, 0, [tab[:16]]])
            else :
                length = tab[7]-tab[6]+1
                if tab[2] >= params['match_identity'] and length >= max(params['match_prop']*tab[12], params['match_len']) and length >= max(params['match_prop2']*tab[12], params['match_len2']) :
                    groups.append(tab[:2].tolist() + [tab[11], tab[2], None, 0, [tab[:16]]])
                if tab[16][3] not in groups2 :
                    groups2[tab[16][3]] = tab[:2].tolist() + tab[16][:2] + [None, 0, [[]]*(len(tab[16])-3)]
                x = [i for i, t in enumerate(tab[16][3:]) if t == tab[15]][0]
                groups2[tab[16][3]][6][x] = tab[:16]
        else :
            tab[2] = -1
    groups.extend(list(groups2.values()))
    overlap = overlap[ids[overlap.T[0]] & ids[overlap.T[1]], :2]
    convA, convB = np.tile(-1, np.max(blastab.T[15])+1), np.tile(-1, np.max(blastab.T[15])+1)
    seq = dict(seq)
    for id, group in enumerate(groups) :
        group[4] = np.zeros(group[6][0][12], dtype=np.uint8)
        group[4].fill(45)
        group[5] = id
        group[6] = np.array(group[6])
        if group[6].shape[0] == 1 :
            convA[group[6].T[15].astype(int)] = id
        else :
            convB[group[6].T[15].astype(int)] = id
        max_sc = 0
        for tab in group[6] :
            matchedSeq = seq[tab[1]][tab[8]-1:tab[9]] if tab[8] < tab[9] else rc(seq[tab[1]][tab[9]-1:tab[8]])
            ms, i, f, sc = [], 0, 0, [0, 0, 0]
            for s, t in re.findall(r'(\d+)([A-Z])', tab[14]) :
                s = int(s)
                if t == 'M' :
                    ms.append(matchedSeq[i:i+s])
                    i += s
                    sc[f] += s
                elif t == 'D' :
                    i += s
                    f = (f-s)%3
                else :
                    ms.append('-'*s)
                    f = (f+s)%3
            group[4][tab[6]-1:tab[7]] = np.array(list(''.join(ms))).view(asc2int).astype(np.uint8)
            max_sc += max(sc[0], sc[f])
        group[2] = max_sc
    overlap = np.vstack([np.vstack([m, n]).T[(m>=0) & (n >=0)] for m in (convA[overlap.T[0]], convB[overlap.T[0]]) \
                         for n in (convA[overlap.T[1]], convB[overlap.T[1]]) ] + [np.vstack([convA, convB]).T[(convA >= 0) & (convB >=0)]])
            
    np.savez_compressed(out_prefix+'.bsn.npz', bsn=np.array(groups, dtype=object), ovl=overlap)
    return out_prefix

def get_map_bsn(prefix, clust, genomes, orthoGroup) :
    if len(genomes) == 0 :
        sys.exit(1)
    taxa = {}
    for g, s in genomes.items() :
        if s[0] not in taxa : taxa[s[0]] = []
        taxa[s[0]].append([g, s[1]])
    
    #bsns = list(map(iter_map_bsn, [(prefix, clust, id, taxon, seq, params) for id, (taxon, seq) in enumerate(taxa.items())]))
    bsns = pool.map(iter_map_bsn, [(prefix, clust, id, taxon, seq, params) for id, (taxon, seq) in enumerate(taxa.items())])
    
    blastab, overlaps = [], []
    ids = 0
    for bsnPrefix in bsns :
        tmp = np.load(bsnPrefix + '.bsn.npz')
        bsn, ovl = tmp['bsn'], tmp['ovl']

        bsn.T[5] += ids
        ovl += ids
        ids += bsn.shape[0]
        blastab.append(bsn)
        overlaps.append(ovl)

        os.unlink(bsnPrefix + '.bsn.npz')
    blastab = np.vstack(blastab)
    blastab.T[1] = np.vectorize(lambda x:genomes.get(x, ['-'])[0])(blastab.T[1])
    overlaps = np.vstack(overlaps)
    ovl_score = np.vectorize(lambda m,n:orthoGroup.get((m,n), 2))(blastab[overlaps.T[0], 0], blastab[overlaps.T[1], 0])
    overlaps = np.hstack([overlaps, ovl_score[:, np.newaxis]])
    
    blastab = blastab[np.argsort(-blastab.T[2])]
    blastab = blastab[np.argsort(blastab.T[0], kind='mergesort')]
    return blastab, overlaps

def checkCDS(n, s) :
    if len(s) < params['min_cds'] :
        #logger('{0} is too short'.format(n))
        return False
    if params['incompleteCDS'] :
        return True

    if len(s) % 3 > 0 :
        #logger('{0} is discarded due to frameshifts'.format(n))
        return False
    aa = transeq({'n':s.upper()}, frame=1, transl_table='starts')['n'][0]
    if aa[0] != 'M' :
        #logger('{0} is discarded due to lack of start codon'.format(n))
        return False
    if aa[-1] != 'X' :
        #logger('{0} is discarded due to lack of stop codon'.format(n))
        return False
    if len(aa[:-1].split('X')) > 1 :
        #logger('{0} is discarded due to internal stop codons'.format(n))
        return False
    return True
    

def addGenes(genes, gene_file) :
    for gfile in gene_file.split(',') :
        if gfile == '' : continue
        gprefix = os.path.basename(gfile).split('.')[0]
        ng = readFasta(gfile)
        for name in ng :
            s = ng[name]
            if checkCDS(name, s) :
                genes['{0}:{1}'.format(gprefix,name)] = [ gfile, '', 0, 0, '+', int(hashlib.sha1(s.encode('utf-8')).hexdigest(), 16), s]
    return genes

def writeGenes(fname, genes, priority) :
    uniques = {}
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
    return fname

def precluster2(data) :
    bsn_file, gene, global_file = data
    matches = np.load(bsn_file)[gene]
    if len(matches) <= 1 :
        return [int(gene), matches]
    global_differences = dict(np.load(global_file))
    gIden = np.hstack([matches[:, [1, 3]], np.arange(matches.shape[0])[:, np.newaxis]])
    ingroup = np.zeros(matches.shape[0], dtype=bool)
    ingroup[0] = True

    for i1, m1 in enumerate(gIden) :
        if ingroup[m1[2]] :
            m2 = gIden[i1+1:][ingroup[gIden[i1+1:, 2].astype(int)] != True]
            if m2.size :
                gs = np.vectorize(lambda g1, g2:global_differences.get(tuple(sorted([g1, g2])), (0.01, 4.) ))(m2.T[0], m1[0])
                sc = -(m2.T[1] - m1[1])/gs[0]/gs[1]
                ingroup[m2[sc < 1, 2].astype(int)] = True
            else :
                break
    return [int(gene), matches[ingroup]]

def precluster(bsn_file, global_file) :
    return dict(pool2.map(precluster2, [(bsn_file, gene, global_file) for gene in np.load(bsn_file).keys()]))
    #return dict(list(map(precluster2, [(bsn_file, gene, global_file) for gene in np.load(bsn_file).keys()])))

def write_output(prefix, prediction, genomes, clust_ref, encodes, old_prediction) :
    predictions, alleles = {}, {}
    
    allele_file = open('{0}.allele.fna'.format(prefix), 'w')
    prediction = pd.read_csv(prediction, sep='\t', header=None)
    prediction = prediction.assign(s=np.min([prediction[9], prediction[10]], 0)).sort_values(by=[3, 5, 's']).drop('s', axis=1).values
    for part in prediction :
    #with open(prediction) as fin :
        #for line in fin :
            #part = line.strip().split()
        gId = encodes[part[0]]
        if part[0] not in alleles :
            alleles[part[0]] = {clust_ref[gId]:1}
            allele_file.write('>{0}_{1}\n{2}\n'.format(part[0], 1, clust_ref[gId]))

        if part[9] < part[10] :
            l, r, d = min(part[7]-1, part[9]-1), min(part[12]-part[8], part[13]-part[10]), 1
        else :
            l, r, d = min(part[7]-1, part[13]-part[9]), min(part[12]-part[8], part[10]-1), -1
        if l <= 6 and part[7] - l == 1 :
            part[7], part[9] = part[7]-l, part[9]-l*d
        else :
            ll = (part[7]-1) % 3
            if ll > 0 :
                part[7], part[9] = part[7]+3-ll, part[9]+(3-ll)*d
        if r <= 6 and part[8] + r == part[12] :
            part[8], part[10] = part[8]+r, part[10]+r*d
        else :
            rr = (part[12] - part[8]) % 3
            if rr > 0 :
                part[8], part[10] = part[8]-3+rr, part[10]-(3+rr)*d

        if part[9] < part[10] :
            part[9:12] = part[9], part[10], '+'
        else :
            part[9:12] = part[10], part[9], '-'
        
        if part[3] not in predictions :
            predictions[part[3]] = []
        elif predictions[part[3]][-1][2] == part[2] :
            prev = predictions[part[3]][-1]
            if prev[5] == part[5] and part[9] - prev[10] < 500 :
                if part[11] == '+' and part[7] - prev[8] < 500 :
                    prev[8], prev[10] = part[8], part[10]
                    continue
                elif part[11] == '-' and prev[7] - part[8] < 500 :
                    prev[7], prev[10] = part[7], part[10]
                    continue
            predictions[part[3]][-1][1], part[1] = -1, -1
        predictions[part[3]].append(part)
    
    op = ['', 0, []]
    with open('{0}.EToKi.gff'.format(prefix), 'w') as fout :
        for gid, (g, predict) in enumerate(predictions.items()) :
            for pid, pred in enumerate(predict) :
                if pred[1] == -1 or (pred[10]-pred[9]+1) <= 0.8 * pred[12] :
                    cds, allele_id = 'fragment:{0:.2f}%'.format((pred[10]-pred[9]+1)*100/pred[12]), 'uncertain'
                    start, stop = pred[9:11]
                else :
                    s, e = pred[9:11]
                    if pred[11] == '+' :
                        s2, e2 = s - min(int(3*((s - 1)/3)), 60), e + min(3*int((pred[13] - e)/3), 600)
                        seq = genomes[encodes[pred[5]]][1][(s2-1):e2]
                        lp, rp = s - s2, e2 - e
                    else :
                        s2, e2 = s - min(int(3*((s - 1)/3)), 600), e + min(3*int((pred[13] - e)/3), 60)
                        seq = rc(genomes[encodes[pred[5]]][1][(s2-1):e2])
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
                    
                    frames = sorted(set([0, len(seq)%3]))
                    for frame, aa_seq in zip(frames, transeq({'n':seq}, transl_table='starts', frame=','.join([str(f+1) for f in frames]))['n']) :
                        cds = 'CDS'
                        s0, s1 = aa_seq.find('M', int(lp/3), int(lp/3+30)), aa_seq.rfind('M', 0, int(lp/3))
                        start = s0 if s0 >= 0 else s1
                        if start < 0 :
                            cds, start = 'nostart', int(lp/3)
                        stop = aa_seq.find('X', start)
                        if 0 <= stop < lp/3+30 :
                            s0 = aa_seq.find('M', stop, int(lp/3+30))
                            if s0 >= 0 :
                                start = s0
                                stop = aa_seq.find('X', start)
                        if stop < 0 :
                            cds = 'nostop'
                        elif (stop - start + 1)*3 <= 0.8 * pred[12] :
                            cds = 'premature stop:{0:.2f}%'.format((stop - start + 1)*300/pred[12])
                            
                        if cds == 'CDS' :
                            if pred[11] == '+' :
                                start, stop = s2 + start*3 + frame, s2 + stop*3 + 2 + frame
                            else :
                                start, stop = e2 - stop*3 - 2 - frame, e2 - start*3 - frame
                            break
                        else :
                            start, stop = s, e
                            if frame > 0 :
                                cds = 'frameshift'

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
                    '' if pred[0] == pred[4] else ',structure variant group:' + pred[4], 
                    '' if cds == 'CDS' else ';pseudogene=' + cds, 
                    '' if len(old_tag) == 0 else 'locus_tag={0};'.format(','.join(old_tag)), 
                ))
    allele_file.close()
    logger('Pan genome annotations have been saved in {0}'.format('{0}.EToKi.gff'.format(prefix)))
    logger('Gene allelic sequences have been saved in {0}'.format('{0}.allele.fna'.format(prefix)))
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

    parser.add_argument('--match_identity', help='minimum identities in BLAST search. Default: 0.6', default=0.5, type=float)
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
    parser.add_argument('--metagenome', help='Set to metagenome mode. equals to "--incompleteCDS --clust_identity 0.99 --clust_match_prop 0.8 --match_identity 0.98 --orthology rapid"', default=False, action='store_true')

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
        params.orthology = 'rapid'
    return params

def encodeNames(genomes, genes) :
    taxon = {g[0] for g in genomes.values()}
    labels = {label:labelId for labelId, label in enumerate(sorted(list(taxon) + list(genomes.keys()) + list(genes.keys())))}
    genes = { labels[gene]:[labels[info[0]], labels[info[1]]] + info[2:] for gene, info in genes.items() }
    genomes = { labels[genome]:[labels[info[0]]] + info[1:] for genome, info in genomes.items() }
    return genomes, genes, labels

params = dict(
    ml = '{fasttree} {0} -nt -gtr -pseudo', 
    nj = '{rapidnj} -i fa -t d {0}', 
)

pool, pool2 = None, None
def ortho(args) :
    global params
    params.update(add_args(args).__dict__)
    params.update(externals)

    global pool, pool2
    pool = Pool(params['n_thread'])
    pool2 = Pool(params['n_thread'])
    
    genomes, genes = readGFF(params['GFFs'])
    genes = addGenes(genes, params['genes'])
    if params.get('old_prediction', None) is None :
        params['old_prediction'] = params['prefix']+'.old_prediction.npz'
        old_predictions = {}
        for n, g in genes.items() :
            if g[1] != '' :
                if g[1] not in old_predictions :
                    old_predictions[g[1]] = []
                old_predictions[g[1]].append([n, g[2], g[3], g[4]])
        for gene, g in old_predictions.items() :
            old_predictions[gene] = np.array(sorted(g), dtype=object)
        np.savez_compressed(params['old_prediction'], **old_predictions)
        del old_predictions, n, g
    
    genomes, genes, encodes = encodeNames(genomes, genes)
    if params.get('prediction', None) is None :
        first_classes = load_priority( params.get('priority', ''), genes, encodes )

        if params.get('clust', None) is None :
            params['genes'] = writeGenes('{0}.genes'.format(params['prefix']), genes, first_classes)
            del genes
            params['clust'], params['uc'] = getClust(params['prefix'], params['genes'], dict(identity=params['clust_identity'], coverage=params['clust_match_prop'], n_thread=params['n_thread']))
        genes = { int(n):s for n, s in readFasta(params['clust']).items()}
        
        if params.get('self_bsn', None) is None :
            params['self_bsn'] = params['prefix']+'.self_bsn.npy'
            orthoGroup = get_similar_pairs(params['prefix'], params['clust'], first_classes, params)
            np.save(params['self_bsn'], orthoGroup)
        else :
            orthoGroup = np.load(params['self_bsn'])
        orthoGroup = dict([[tuple(g),1] for g in orthoGroup] + [[(g[1], g[0]),1] for g in orthoGroup] + [[(g, g), 0] for g in genes])
        
        if params.get('map_bsn', None) is None or params.get('conflicts', None) is None :
            blastab, conflicts = get_map_bsn(params['prefix'], params['clust'], genomes, orthoGroup)
            blastab = np.split(blastab, np.cumsum(np.unique(blastab.T[0], return_counts=True)[1])[:-1])
            
            params['map_bsn'], params['conflicts'] = params['prefix']+'.map_bsn.npz', params['prefix']+'.conflicts.npz'
            np.savez_compressed(params['map_bsn'], **{str(b[0,0]):b for b in blastab})
            np.savez_compressed(params['conflicts'], conflicts=conflicts)
            del blastab, conflicts
        pool.close()
        pool.join()
        
        if params.get('orthology', 'rapid') != 'rapid' :
            if params.get('global', None) is None :
                params['global'] = params['prefix']+'.global.npy'
                global_differences = global_difference(params['map_bsn'], params['prefix']+'.global', orthoGroup, 3000)
                np.save(params['global'], global_differences)
                del global_differences
            
            blastab = precluster(params['map_bsn'], params['global'])
        else :
            blastab = { int(k):v for k, v in np.load(params['map_bsn']).items()}
        params['prediction'] = filt_genes(params['prefix'], blastab, params['global'], np.load(params['conflicts'])['conflicts'], first_classes, encodes)
    else :
        genes = {n:s[-1] for n,s in genes.items() }
    pool2.close()
    pool2.join()
    old_predictions = dict(np.load(params['old_prediction'])) if 'old_prediction' in params else {}
    
    write_output(params['prefix'], params['prediction'], genomes, genes, encodes, old_predictions)
    
if __name__ == '__main__' :
    ortho(sys.argv[1:])

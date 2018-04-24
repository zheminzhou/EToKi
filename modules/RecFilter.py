import numpy as np, pandas as pd, sys, math, argparse, gzip
import phylo 
from ete3 import Tree

def parse_arg(a) :
    parser = argparse.ArgumentParser(description='Generate a matrix of only vertically inherited SNPs. ', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--prefix', '-p', help='prefix for the output', required=True)
    parser.add_argument('--snp', '-s', dest='matrix', help='SNP matrix', required=True)
    parser.add_argument('--tree', '-t', help='Labeled tree', required=True)
    parser.add_argument('--rec', '-r', help='Recombinant sketches', required=True)
    parser.add_argument('--prob', '-b', help='Minimum probability for rec sketches. Default: 0.5', default=0.5, type=float)
    parser.add_argument('--clonalframeml', help='The recombinant sketches are in ClonalFrameML format. ', default=False, action="store_true")
    parser.add_argument('--simbac', help='The recombinant sketches are in SimBac break format. ', default=False, action="store_true")
    args = parser.parse_args(a)
    return args

def read_clonalframe(fname, nodes) :
    rec = {}
    with open(fname) as fin :
        fin.readline()
        for line in fin :
            branch, s, e = line.strip().split('\t')
            if branch not in rec :
                rec[branch] = []
            rec[branch].append([int(s), int(e), '', nodes[branch]])
    return rec

def read_simbac(fname, nodes) :
    node_map = {'[{0}]'.format(','.join(sorted(n, key=lambda x:int(x)))) : b for b,n in nodes.iteritems()}
    node_map['[EXTERNAL]'] = 'External'
    rec = {}
    with open(fname) as fin :
        fin.readline()
        fin.readline()
        for line in fin :
            s, e, b1, b2 = line.strip().split('\t')
            if b1 == b2 : 
                continue
            branch = node_map[b1]
            if branch not in rec :
                rec[branch] = []
            rec[branch].append([int(s), int(e), node_map[b2], nodes[branch]])
    return {b:sorted(r) for b, r in rec.iteritems() }


def read_RecHMM(fname, nodes, prob) :
    rec = {}
    with open(fname) as fin :
        for line in fin :
            if line.startswith('\tRecomb') :
                part = line.strip().split('\t')
                if len(part) > 5 and float(part[5]) < prob :
                    continue
                if part[1] not in rec :
                    rec[part[1]] = []
                rec[part[1]].append([int(part[2]), int(part[3]), part[4], nodes[part[1]]])
    return rec

def write_filtered_matrix(fname, names, sites, snps, masks, m_weight) :
    invariants = { snp[0]:[base, snp[2]] for base, snp in snps.iteritems() if snp[1] == 0 and base[0] != '-' }
    bases = {}
    for inv in invariants.values() :
        bases[inv[0]] = bases.get(inv[0], 0) + inv[1]
    sv = {ss[0]:s.split('\t') for s, ss in snps.iteritems()}
    name_map = {name:id for id, name in enumerate(names)}
    with gzip.open(fname, 'w') as fout :
        fout.write('## Constant_bases: ' + ' '.join([str(inv[1]) for inv in sorted(bases.iteritems())]) + '\n')
        fout.write('#seq\t#site\t' + '\t'.join(names) + '\t#!W[RecFilter]\n')
        for site in sites :
            if not len(m_weight[site[1]]) : continue
            weight = np.mean(m_weight[site[1]].values()) 
            snvs = np.array(sv[site[2]])
            snv_x = []
            p = np.zeros(snvs.shape, dtype=bool)
            for m in masks.get(site[1], []) :
                pp = np.ones(snvs.shape, dtype=bool)
                pp[ [name_map[mm] for mm in m] ] = False
                p = (p | (~pp))
                snv_x.append(np.copy(snvs))
                snv_x[-1][pp] = '-'
            snv_x.append(snvs)
            snv_x[-1][p] = '-'
            
            for snv in snv_x :
                snv_type = np.unique(snv)
                if snv_type[~ np.in1d(snv_type, ['-', 'N', 'n'])].size > 1 :
                    fout.write('{2}\t{3}\t{0}\t{1:.5f}\n'.format('\t'.join(snv.tolist()), weight, *site[:2]))
    return fname


def RecFilter(argv) :
    args = parse_arg(argv)
    names, sites, snps = phylo.read_matrix(args.matrix)
    snp_list = sorted([[info[0], int(math.ceil(info[2])), line.split('\t'), info[1]] for line, info in snps.iteritems() ])
    tree = Tree(args.tree, format=1)
    
    final_tree, node_names, states = phylo.infer_ancestral(args.tree, names, snp_list, sites, infer='viterbi')
    nodes = {}
    for node in final_tree.traverse('postorder') :
        if node.is_leaf() :
            nodes[node.name] = [node.name]
        else :
            nodes[node.name] = [ d for c in node.children for d in nodes[c.name] ]

    mutations = phylo.get_mut(final_tree, node_names, np.array(states), sites)
    if args.clonalframeml :
        rec_regions = read_clonalframe(args.rec, nodes)
    elif args.simbac :
        rec_regions = read_simbac(args.rec, nodes)
    else :
        rec_regions = read_RecHMM(args.rec, nodes, args.prob)
        
    n_base = np.sum([v[2] for v in snps.values()])
    br_weight = { b:(n_base+.1)/(n_base-np.sum([ rr[1]-rr[0]+1 for rr in r ])+.1) for b, r in rec_regions.iteritems() }
    curr = ['', [], 0]
    masks = {}
    m_weight = {}
    for m in mutations :
        if curr[0] != m[0] :
            curr = [m[0], rec_regions.get(m[0], []), 0]
        if m[2] not in m_weight :
            m_weight[m[2]] = {}
        m_weight[m[2]][m[0]] = br_weight.get(m[0], 1.)
        for id in range(curr[2], len(curr[1])) :
            r = curr[1][id]
            if m[2] > r[1] :
                curr[2] += 1
            elif m[2] < r[0] :
                break
            elif m[2] not in masks :
                masks[m[2]] = [r[3]]
                m_weight[m[2]].pop(m[0])
            else :
                masks[m[2]].append(r[3])
                m_weight[m[2]].pop(m[0])
    write_filtered_matrix(args.prefix+'.filtered.gz', names, sites, snps, masks, m_weight)
    return

if __name__ == '__main__' :
    RecFilter(sys.argv[1:])

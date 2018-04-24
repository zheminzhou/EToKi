from ete3 import Tree
import sys, numpy as np, os, glob, math, gzip, re, argparse
from subprocess import Popen, PIPE
from configure import externals

raxml = externals['raxml']

def readFasta(fasta_file) :
    seqs = []
    
    fin = Popen(['zcat', fasta_file],stdout=PIPE).stdout if fasta_file[-3:].lower() == '.gz' else open(fasta_file)
    for line in fin :
        if line.startswith('>') :
            name = line[1:].strip().split()[0]
            seqs.append([name, []])
        else :
            seqs[-1][1].extend(line.strip().split())
    fin.close()
    return [ [n, ''.join(s)] for n, s in seqs ]


def parse_snps(seq, core=0.95) :
    names = np.array([ s[0] for s in seq ])
    seqs = np.array([ list(re.sub(r'[^ACGT]', r'-', s[1].upper())) for s in seq ])
    
    snps ={}
    sites = []
    const_sites = {
        'A' : '\t'.join(['A' for b in seqs.T[0]]), 
        'C' : '\t'.join(['C' for b in seqs.T[0]]), 
        'G' : '\t'.join(['G' for b in seqs.T[0]]), 
        'T' : '\t'.join(['T' for b in seqs.T[0]]), 
    }
    
    ref_site = 0
    for bases in seqs.T :
        if bases[0] != '-' :
            ref_site += 1

        types = dict(zip(*np.unique(bases, return_counts=True)))
        types.pop('-', 0)
        n_type = len(types)

        if n_type > 0 and np.sum(types.values()) >= core * bases.size :
            if n_type > 1 :
                b_key = '\t'.join(bases)
            else :
                b_key = const_sites[types.keys()[0]]

            if b_key not in snps :
                snps[b_key] = [len(snps), n_type-1, 1.]
            else :
                snps[b_key][2] += 1.
            if n_type > 1 :
                sites.append([names[0], ref_site, snps[b_key][0]])
    return names, sites, snps

def write_phylip(prefix, names, snp_list) :
    invariants = {'A':0, 'C':0, 'G':0, 'T':0, '-':0}
    valid_values = {'A':'A', 'C':'C', 'G':'G', 'T':'T', '-':'-', 'N':'-'}
    for snp in snp_list :
        if snp[3] == 0 and snp[2][0].upper() in invariants :
            invariants[ snp[2][0] ] += snp[1]
    
    snp2 = [snp for snp in snp_list if snp[3] > 0 and snp[2][0] in invariants]
    weights = [ snp[1] for snp in snp2 ]
    with open(prefix+'.phy.weight', 'w') as fout :
        fout.write(' '.join([str(x) for x in weights]))

    snp_array = np.array([ [ valid_values.get(s, '-') for s in snp[2]] for snp in snp2 ]).T
    n_tax, n_seq = snp_array.shape
    with open(prefix + '.phy', 'w') as fout :
        fout.write('\t{0} {1}\n'.format(n_tax, n_seq))
        for id, n in enumerate(names) :
            fout.write('{0} {1}\n'.format(n, ''.join(snp_array[id])))

    if sum(invariants.values()) > 0 :
        asc_file = prefix + '.asc'
        constant_file = prefix + '.constant'
        with open(asc_file, 'w') as fout :
            fout.write('[asc~{0}], ASC_DNA,p1=1-{1}\n'.format(constant_file,n_seq))
        with open(constant_file, 'w') as fout :
            fout.write(' '.join([str(int(x+0.5)) for y,x in sorted(invariants.items())[1:]]) + '\n')
    else : 
        asc_file = None
        constant_file = None

    return prefix+'.phy' , prefix + '.phy.weight', asc_file
    
def run_raxml(prefix, phy, weights, asc, model='CAT') :
    for fname in glob.glob('RAxML_*.{0}'.format(prefix)) :
        os.unlink(fname)
    if asc is None :
        cmd = '{0} -m GTR{4} -n {1} -s {2} -a {3} -T 7 -p 1234'.format(raxml, prefix, phy, weights, model)
    else :
        cmd = '{0} -m ASC_GTR{5} -n {1} -s {2} -a {3} -T 7 -p 1234 --asc-corr stamatakis -q {4}'.format(raxml, prefix, phy, weights, asc, model)
    run = Popen(cmd.split())
    run.communicate()
    if model == 'CAT' and not os.path.isfile('RAxML_bestTree.{0}'.format(prefix)) :
        return run_raxml(prefix, phy, weights, asc, 'GAMMA')
    else :
        fname = '{0}.unrooted.nwk'.format(prefix)
        os.rename('RAxML_bestTree.{0}'.format(prefix), fname)
        for fn in glob.glob('RAxML_*.{0}'.format(prefix)) :
            try:
                os.unlink(fn)
            except :
                pass
        return fname

def get_root(prefix, tree_file) :
    tree = Tree(tree_file, format=1)
    try:
        tree.set_outgroup( tree.get_midpoint_outgroup() )
    except :
        pass
    tree.write(outfile='{0}.rooted.nwk'.format(prefix), format=1)
    return '{0}.rooted.nwk'.format(prefix)

def write_matrix(fname, names, sites, snps) :
    invariants = { snp[0]:[base, snp[2]] for base, snp in snps.iteritems() if snp[1] == 0 and base[0] != '-' }
    bases = {}
    for inv in invariants.values() :
        bases[inv[0]] = bases.get(inv[0], 0) + inv[1]
    sv = {ss[0]:s for s, ss in snps.iteritems()}
    with gzip.open(fname, 'w') as fout :
        fout.write('## Constant_bases: ' + ' '.join([str(inv[1]) for inv in sorted(bases.iteritems())]) + '\n')
        fout.write('#seq\t#site\t' + '\t'.join(names) + '\n')
        for site in sites :
            fout.write('{1}\t{2}\t{0}\n'.format(sv[site[2]], *site[:2]))
    return fname

def read_matrix(fname) :
    sites, snps = [], {}
    invariant = []
    fin = Popen(['zcat', fname],stdout=PIPE).stdout if fname[-3:].lower() == '.gz' else open(fname)
    
    for line_id, line in enumerate(fin) :
        if line.startswith('##'): 
            part = line[2:].strip().split()
            invariant = zip(['A', 'C', 'G', 'T'], [float(v) for v in part[1:]])
        elif line.startswith('#') :
            part = np.array(line.strip().split('\t'))
            cols = (1 - np.char.startswith(part, '#')).astype(bool)
            w_cols = np.char.startswith(part, '#!W')
            names = part[cols]
        elif line_id == 0 :
            part = np.array(line.strip().split('\t'))
            cols = np.ones(part.shape, dtype=bool)
            cols[:2] = False
            w_cols = np.char.startswith(part, '#!W')
            names = part[cols]
        else :
            p2 = np.array(line.strip().split('\t'))
            part = np.char.upper(p2[cols])
            types = dict(zip(*np.unique(part, return_index=True)))
            types.pop('-', None)
            b_key = '\t'.join(part)
            if np.sum(w_cols) :
                w = np.multiply.reduce(p2[w_cols].astype(float))
            else :
                w = 1.
            if b_key not in snps :
                snps[b_key] = [len(snps), len(types)-1, w]
            else :
                snps[b_key][2] += w
            if len(types) > 1 :
                sites.append([ p2[0], int(p2[1]), snps[b_key][0] ])
    fin.close()
    for inv in invariant :
        b_key = '\t'.join([inv[0]] * len(names))
        if b_key not in snps :
            snps[b_key] = [len(snps), 0, float(inv[1])]
        else :
            snps[b_key][2] += float(inv[1])
    return names, sites, snps

def read_ancestor(fname, names, snp_list) :
    snp_array = np.array([snp[2] for snp in snp_list]).T
    branches = names[:]
    with open(fname) as fin :
        for line in fin :
            name, seq = line.strip().split()
            branches.append(name)
            snp_array = np.concatenate([snp_array, np.array([list(seq)])])
    return dict(zip(branches, snp_array)), [n[0] for n in snp_list]

def get_mut(final_tree, names, states, sites) :
    mutations = {}
    branches = []
    name_ids = { n:id for id, n in enumerate(names) }
    if final_tree.name == '' :
        root = {n:1 for n in names}
        for n in final_tree.traverse() : root.pop(n.name, None)
        if len(root) == 1 :
            final_tree.name = root.keys()[0]
    states = states.T
    for node in final_tree.iter_descendants('postorder') :
        for id, (m, n) in enumerate(zip(states[name_ids[node.name]], states[name_ids[node.up.name]])) :
            if m != n and m != '-' and n != '-' :
                if id not in mutations :
                    mutations[id] = []
                mutations[id].append([node.name, '{0}->{1}'.format(n,m)])
    outputs = []
    for c, p, i in sites :
        m = mutations.get(i, [])
        len_m = {}
        for mut in m :
            s = tuple(sorted([mut[1][0], mut[1][-1]]))
            mut.append(s)
            len_m[s] = len_m.get(s, 0) + 1
        for mut in m :
            outputs.append([mut[0], c, p, len_m[mut[-1]], mut[1]])
    return sorted(outputs)

def write_states(fname, names, states, sites) :
    with gzip.open(fname, 'w') as fout :
        fout.write('#Seq\t#Site\t' + '\t'.join(names) + '\n')
        for site in sites :
            fout.write('{0}\t{1}\t{2}\n'.format(site[0], site[1], '\t'.join(states[site[2]]) ))

def write_ancestral_proportion(fname, names, states, sites) :
    with gzip.open(fname, 'w') as fout :
        fout.write('#Seq\t#Site\t#Type:Proportion\n')
        for c, p, i in sites :
            tag, state = states[i]
            for n, ss in zip(names, state) :
                fout.write( '{0}\t{1}\t{2}\t{3}\n'.format(c, p, n, '\t'.join([ '{0}:{1:.5f}'.format(t, s) for t, s in zip(tag, ss)]) ))


def read_states(fname) :
    names, ss, sites = [], {}, []
    with gzip.open(fname) as fin :
        names = fin.readline().strip().split('\t')[2:]
        for line in fin :
            seq, site, snp_str = line.strip().split('\t', 2)
            if snp_str not in ss :
                ss[snp_str] = len(ss)
            sites.append([seq, int(site), ss[snp_str]])
    states = []
    for s, id in sorted(ss.iteritems(), key=lambda x:x[1]) :
        states.append(s.split('\t'))
    return names, np.array(states), sites
    

def add_args(a) :
    parser = argparse.ArgumentParser(description='Parameters for phylogeny_workflow.py', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--tasks', '-t', help='''Tasks to call. Allowed tasks are:
matrix: generate SNP matrix from alignment.
phylogeny: generate phylogeny from SNP matrix.
ancestral: generate AS (ancestral state) matrix from SNP matrix and phylogeny
ancestral_proportion: generate possibilities of AS for each site
mutation: assign SNPs into branches from AS matrix

You can run multiple tasks by sending a comma delimited task list.
There are also some pre-defined task combo:
all: matrix,phylogeny,ancestral,mutation
aln2phy: matrix,phylogeny [default]
snp2anc: phylogeny,ancestral
mat2mut: ancestral,mutation''', default='aln2phy')
    
    parser.add_argument('--prefix', '-p', help='prefix for all outputs.', required=True)
    parser.add_argument('--alignment', '-m', help='aligned sequences in either fasta format or Xmfa format. Required for "matrix" task.', default='')
    parser.add_argument('--snp', '-s', help='SNP matrix in specified format. Required for "phylogeny" and "ancestral" if alignment is not given', default='')
    parser.add_argument('--tree', '-z', help='phylogenetic tree. Required for "ancestral" task', default='')
    parser.add_argument('--ancestral', '-a', help='Inferred ancestral states in a specified format. Required for "mutation" task', default='')
    parser.add_argument('--core', '-c', help='Core genome proportion. Default: 0.95', type=float, default=0.95)
    parser.add_argument('--n_proc', '-n', help='Number of processes. Default: 5. ', type=int, default=5)
    
    args = parser.parse_args(a)
    
    args.tasks = dict(
        all = 'matrix,phylogeny,ancestral,mutation',
        aln2phy = 'matrix,phylogeny',
        snp2anc = 'phylogeny,ancestral', 
        mat2mut = 'ancestral,mutation',
    ).get(args.tasks, args.tasks).split(',')
    
    return args

def infer_ancestral(tree, names, snps, sites, infer='margin', rescale=1.0) :
    tree = Tree(tree, format=1)
    node_names = {}
    for id,branch in enumerate(tree.traverse('postorder')) :
        if branch.name == '' :
            branch.name = 'N_' + str(len(node_names))
        if not branch.up and len(branch.children) == 2 :
            branch.children[1].dist += branch.children[0].dist - 1e-8
            branch.children[0].dist = 1e-8
        node_names[str(branch.name)] = id
    states, branches = [ [ '-' for snp in snps ] for br in node_names ], []
    
    for n, s in zip(names, np.array([ [ b.upper() for b in snp[2] ] for snp in snps ]).T) :
        states[ node_names[n] ] = s
    states = np.array(states).T
    for branch in tree.traverse('postorder') :
        if branch.up :
            branches.append([ node_names[branch.up.name], node_names[branch.name], np.exp(-max(branch.dist, 1e-8)) ])
        else :
            branches.append([ None, node_names[branch.name], 1e-8 ])

    n_node = len(node_names)
    transitions = {}
    retvalue = []
    for state in states :
        tag, code = np.unique(state, return_inverse=True)
        n_state = tag.size
        missing = np.where(tag == '-')[0]
        if missing.size > 0 :
            tag, n_state = tag[tag != '-'], n_state - 1
            code[code == missing[0]] = -1
            code[code  > missing[0]] = code[code  > missing[0]] - 1
        if np.sum(np.in1d(tag, ['A', 'C', 'G', 'T'])) == n_state :
            n_state = 4
        
        if n_state not in transitions :
            transitions[n_state] = np.zeros(shape=[n_node, n_state, n_state])
            for tr, (s, t, v) in zip(transitions[n_state], branches) :
                tr.fill((1.0-v)/n_state)
                np.fill_diagonal(tr, (1.+(n_state-1.)*v)/n_state)
        
        transition = transitions[n_state]

        if infer == 'margin' :
            alpha = np.ones(shape=[n_node, n_state])/n_state
            alpha[code >= 0] = 0        
            alpha[code >= 0, code[code >= 0]] = 1
            beta = np.ones(alpha.shape)
            for (s, t, v), tr in zip(branches, transition) :
                alpha[t] = alpha[t]/np.sum(alpha[t])
                if s :
                    beta[t] = np.dot(alpha[t], tr)
                    alpha[s] *= beta[t]
            
            for (s, t, v), tr in reversed(zip(branches, transition)) :
                if s :
                    alpha[t] *= np.dot(alpha[s]/beta[t], tr)
            retvalue.append([tag, alpha])
        else :
            pt = np.log(transition)
            ids = np.arange(n_state)
            alpha = np.zeros(shape=[n_node, n_state, n_state])
            path = np.zeros(shape=[n_node, n_state], dtype=int)
            alpha[code >= 0] = -9999
            alpha[code >= 0, :, code[code >=0]] = 0
            
            for (s, t, v), tr in zip(branches, pt) :
                x = alpha[t] + tr
                path[t] = np.argmax(x, 1)
                alpha[s] += x[ids, path[t]]

            r = np.zeros(shape=[n_node], dtype=int)
            for s, t, v in reversed(branches) :
                if not s :
                    r[t] = np.argmax(alpha[t, 0])
                else :
                    r[t] = path[t][r[s]]
            retvalue.append(tag[r])
    return tree, [ k for k, v in sorted(node_names.iteritems(), key=lambda x:x[1])], retvalue

def phylo(args) :
    args = add_args(args)

    prefix = args.prefix
    
    if 'matrix' in args.tasks :
        assert os.path.isfile( args.alignment )
        seq = readFasta( args.alignment )
        names, sites, snps = parse_snps(seq, args.core)
        args.snp = write_matrix(prefix+'.matrix.gz', names, sites, snps)
        if len(names) < 4 :
            raise ValueError('Taxa too few.')
        snp_list = sorted([[info[0], int(math.ceil(info[2])), line.split('\t'), info[1]] for line, info in snps.iteritems() ])
    elif 'phylogeny' in args.tasks or 'ancestral' in args.tasks or 'ancestral_proportion' in args.tasks :
        assert os.path.isfile( args.snp )
        names, sites, snps = read_matrix(args.snp)
        if len(names) < 4 :
            raise ValueError('Taxa too few.')
        snp_list = sorted([[info[0], int(math.ceil(info[2])), line.split('\t'), info[1]] for line, info in snps.iteritems() ])
    
        
    # build tree
    if 'phylogeny' in args.tasks :
        phy, weights, asc = write_phylip(prefix+'.tre', names, snp_list)
        if phy != '' :
            args.tree = run_raxml(prefix +'.tre', phy, weights, asc)
        else :
            args.tree = prefix + '.tre'
            with open(args.tree, 'w') as fout :
                fout.write('({0}:0.0);'.format(':0.0,'.join(names)))
        args.tree = get_root(prefix, args.tree)
    elif 'ancestral' in args.tasks or 'ancestral_proportion' in args.tasks :
        tree = Tree(args.tree, format=1)
    
    # map snp
    if 'ancestral' in args.tasks :
        final_tree, node_names, states = infer_ancestral(args.tree, names, snp_list, sites, infer='viterbi')
        states = np.array(states)
        final_tree.write(format=1, outfile=prefix + '.labelled.nwk')
        write_states(prefix+'.ancestral_states.gz', node_names, states, sites)
    elif 'mutation' in args.tasks :
        final_tree = Tree(args.tree, format=1)
        node_names, states, sites = read_states(args.ancestral)
        
    if 'ancestral_proportion' in args.tasks :
        final_tree, node_names, states = infer_ancestral(args.tree, names, snp_list, sites, infer='margin')
        final_tree.write(format=1, outfile=prefix + '.labelled.nwk')
        write_ancestral_proportion(prefix+'.ancestral_proportion.gz', node_names, states, sites)

    if 'mutation' in args.tasks :
        mutations = get_mut(final_tree, node_names, states, sites)
        with gzip.open(prefix + '.mutations.gz', 'w') as fout :
            fout.write('#Node\t#Seq\t#Site\t#Homoplasy\t#Mutation\n')
            for mut in mutations :
                fout.write('\t'.join([str(m) for m in mut]) + '\n')

if __name__ == '__main__' :
    phylo(sys.argv[1:])

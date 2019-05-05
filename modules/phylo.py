from ete3 import Tree
import sys, numpy as np, os, glob, math, re, argparse
from subprocess import Popen, PIPE
from multiprocessing import Pool
import random
import pandas as pd

rint = random.randint(0, 262144)

try :
    from configure import externals, uopen, asc2int
except :
    from .configure import externals, uopen, asc2int

raxml = externals['raxml']

def fillMissingSeq(seqs, block_id) :
    check = False
    for s in seqs :
        if len(s) > 2 :
            s[2] = ''.join(s[2]).upper()
            check = True
    if not check :
        return False
    alnSize = max([ len(s[2]) for s in seqs ])
    for s in seqs :
        if len(s) == 1 :
            s.extend(['s{0}'.format(block_id), ''])
        s[2] += '-' * (alnSize - len(s[2]))
    return True
    

def xFasta2Matrix(prefix, fasta_file, core=0.95) :
    seqs = []
    snp_data = []
    nameMap = {}
    with uopen(fasta_file) as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                seqName = name.split(':', 1)[0]
                contName = name.split(':', 1)[-1]
                if seqName in nameMap :
                    seqId = nameMap[seqName]
                    seqs[seqId] = [seqName, contName, []]
                else :
                    seqId = nameMap[seqName] = len(seqs)
                    seqs.append([seqName, contName, []])
            elif line.startswith('=') :
                if fillMissingSeq(seqs, len(snp_data)) :
                    snp_data.append(parse_snps(prefix, len(snp_data), seqs, core))
                seqs = [ [n] for n,i in sorted(nameMap.items(), key=lambda x:x[1]) ]
            else :
                seqs[seqId][2].extend(line.strip().split())
    if fillMissingSeq(seqs, len(snp_data)) :
        snp_data.append(parse_snps(prefix, len(snp_data), seqs, core))
    
    const_sites = np.sum([ snp[2] for snp in snp_data ], axis=0)
    names = [n for n,i in sorted(nameMap.items(), key=lambda x:x[1])]
    with uopen(prefix+'.matrix.gz', 'w') as fout :
        fout.write('## Constant_bases: ' + ' '.join(const_sites.astype(str)) + '\n')
        for snp in snp_data :
            fout.write('## Sequence_length: {0} {1}\n'.format(*snp[:2]))
        for snp in snp_data :
            for s, e in snp[3] :
                fout.write('## Missing_region: {0} {1} {2}\n'.format(snp[0], s, e))
        fout.write('#seq\t#site\t' + '\t'.join(names) + '\n')
        for snp in snp_data :
            d = np.load(snp[4])
            sites, sv = d['sites'], d['snps']
            for s in sites :
                fout.write('{0}\t{1}\t{2}\n'.format(snp[0], s[0], '\t'.join(sv[s[1]].astype(str))))
            os.unlink(snp[4])
    return prefix+'.matrix.gz'
    
def parse_snps(prefix, id, seq, core=0.95) :
    #names = np.array([ s[0] for s in seq ])
    missing = []
    snps, sites ={}, []

    const_sites = {
        b'A' : 0., b'C' : 0.,
        b'G' : 0., b'T' : 0.,
    }
    type_id = 0

    contName = seq[0][1]
    seqs = np.array([ (re.sub(r'[^ACGT]', r'-', s[2].upper())) for s in seq ], dtype='c', order='F')

    for ref_site, bases in enumerate(seqs.T) :
        b_key = tuple(bases)
        if b_key in snps :
            s = snps[b_key]
            if s[0] == -2 :
                missing.append(ref_site+1)
            elif s[0] == -1 :
                const_sites[s[2]] += s[3]
            else :
                s[1] += 1
                sites.append([ref_site+1, s[0]])
            continue
        
        types, counts = np.unique(bases, return_counts=True)
        pType = (types != b'-')
        types, counts = types[pType], counts[pType]
        countSum = np.sum(counts)/bases.size

        if types.size > 0 and countSum >= core :
            if types.size == 1 :
                snps[b_key] = [-1, 0., types[0], 1]  #countSum]
                const_sites[types[0]] += 1  #countSum
            else :
                snps[b_key] = [type_id, 1.]
                sites.append([ref_site+1, type_id])
                type_id += 1
        else :
            snps[b_key] = [-2]
            missing.append(ref_site+1)
    if len(missing) : 
        missing2 = [[missing[0], missing[0]]]
        for m in missing[1:] :
            if missing2[-1][1] +1 == m :
                missing2[-1][1] = m
            else :
                missing2.append([m, m])
        missing = []
        del missing
    else :
        missing2 = []
    snps = np.array([ k for k,v in sorted([[k,v] for k,v in snps.items() if v[0]>=0], key=lambda v:v[1][0]) ])
    #snps = pd.DataFrame([ [k, v[1]] for k,v in sorted([[k,v] for k,v in snps.items() if v[0]>=0], key=lambda v:v[1][0]) ]).values
    const_sites = np.array([const_sites[b'A'],const_sites[b'C'],const_sites[b'G'],const_sites[b'T']])
    outputs = dict(sites = np.array(sites), snps = snps)
    np.savez_compressed('{0}.{1}.npz'.format(prefix, id), **outputs)
    return contName, seqs.shape[1], const_sites, missing2, '{0}.{1}.npz'.format(prefix, id)

def write_phylip(prefix, names, snp_list) :
    invariants = {'A':0, 'C':0, 'G':0, 'T':0, '-':0}
    for snp in snp_list :
        if snp[3] == 0 and snp[2][0].upper() in invariants :
            invariants[ snp[2][0] ] += snp[1]

    snp2 = [snp for snp in snp_list if snp[3] > 0 and snp[2][0] in invariants and max([len(s) for s in snp[2] ]) == 1]
    weights = [ snp[1] for snp in snp2 ]
    with open(prefix+'.phy.weight', 'w') as fout :
        fout.write(' '.join([str(x) for x in weights]))

    snp_array = np.array([s[2] for s in snp2]).T
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

def run_raxml(prefix, phy, weights, asc, model='CAT', n_proc=5) :
    for fname in glob.glob('RAxML_*.{0}'.format(prefix)) :
        os.unlink(fname)
    if asc is None :
        if model == 'CAT' :
            cmd = '{0} -m GTR{4} -n {1} -f D -V -s {2} -a {3} -T {5} --no-bfgs -p {6}'.format(raxml, prefix, phy, weights, model, n_proc, rint)
        else :
            cmd = '{0} -m GTR{4} -n {1} -f D -s {2} -a {3} -T {5} -p {6} --no-bfgs'.format(raxml, prefix, phy, weights, model, n_proc, rint)
    else :
        if model == 'CAT' :
            cmd = '{0} -m ASC_GTR{5} -n {1} -f D -V -s {2} -a {3} -T {6} -p {7} --asc-corr stamatakis --no-bfgs -q {4}'.format(raxml, prefix, phy, weights, asc, model, n_proc, rint)
        else :
            cmd = '{0} -m ASC_GTR{5} -n {1} -f D -s {2} -a {3} -T {6} -p {7} --asc-corr stamatakis --no-bfgs -q {4}'.format(raxml, prefix, phy, weights, asc, model, n_proc, rint)
    run = Popen(cmd.split())
    run.communicate()
    if model == 'CAT' and not os.path.isfile('RAxML_bestTree.{0}'.format(prefix)) :
        return run_raxml(prefix, phy, weights, asc, 'GAMMA', n_proc)
    else :
        cmd = '{0} -m GTRCAT -n 2.{1} -f b -z RAxML_rellBootstrap.{1} -t RAxML_bestTree.{1}'.format(raxml, prefix)
        Popen(cmd.split()).communicate()
        fname = '{0}.unrooted.nwk'.format(prefix)
        os.rename('RAxML_bipartitions.2.{0}'.format(prefix), fname)
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

def read_matrix(fname) :
    sites, snps = [], {}
    invariant = []
    seqLens, missing = [], []
    with uopen(fname) as fin :
        for line_id, line in enumerate(fin) :
            if line.startswith('##'):
                if line.startswith('## Constant_bases') :
                    part = line[2:].strip().split()
                    invariant = zip(['A', 'C', 'G', 'T'], [float(v) for v in part[1:]])
                elif line.startswith('## Sequence_length:') :
                    part = line[2:].strip().split()
                    seqLens.append([part[1], int(part[2])])
                elif line.startswith('## Missing_region:') :
                    part = line[2:].strip().split()
                    missing.append([part[1], int(part[2]), int(part[3])])
            elif line.startswith('#') :
                part = np.array(line.strip().split('\t'))
                cols = np.where((1 - np.char.startswith(part, '#')).astype(bool))[0]
                w_cols = np.where(np.char.startswith(part, '#!W'))[0]
                names = part[cols]
                break
            else :
                part = np.array(line.strip().split('\t'))
                cols = np.ones(part.shape, dtype=bool)
                cols[:2] = False
                w_cols = np.char.startswith(part, '#!W')
                names = part[cols]
                break
        mat = pd.read_csv(fin, header=None, sep='\t', usecols=cols.tolist() + w_cols.tolist() + [0,1]).values
    
    val = {'A':'A', 'C':'C', 'G':'G', 'T':'T', '-':'-', 'N':'-', '.':'.'}
    validate = np.vectorize(lambda b: b if len(b) > 1 else val.get(b, '-'))
    
    for p2 in mat :
        part = validate(np.char.upper(p2[cols].astype(str)))
        b_key = tuple(part)
        w = np.multiply.reduce(p2[w_cols].astype(float)) if w_cols.size else 1.
        
        if b_key in snps :
            snps[b_key][2] += w
        else :
            types = dict(zip(*np.unique(part, return_index=True)))
            types.pop('-', None)
            snps[b_key] = [len(snps), len(types)-1, w]
            
        if snps[b_key][1] > 0 :
            sites.append([ p2[0], int(p2[1]), snps[b_key][0] ])

    for inv in invariant :
        b_key = tuple([inv[0]] * len(names))
        if b_key not in snps :
            snps[b_key] = [len(snps), 0, float(inv[1])]
        else :
            snps[b_key][2] += float(inv[1])
    return names, sites, snps, seqLens, missing

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

def write_states(fname, names, states, sites, seqLens, missing) :
    with uopen(fname, 'w') as fout :
        for sl in seqLens :
            fout.write('## Sequence_length: {0} {1}\n'.format(*sl))
        for ms in missing :
            fout.write('## Missing_region: {0} {1} {2}\n'.format(*ms))
        fout.write('#Seq\t#Site\t' + '\t'.join(names) + '\n')
        for site in sites :
            fout.write('{0}\t{1}\t{2}\n'.format(site[0], site[1], '\t'.join(states[site[2]]) ))

def write_ancestral_proportion(fname, names, states, sites, seqLens, missing) :
    with uopen(fname, 'w') as fout :
        for sl in seqLens :
            fout.write('## Sequence_length: {0} {1}\n'.format(*sl))
        for ms in missing :
            fout.write('## Missing_region: {0} {1} {2}\n'.format(*ms))
        
        fout.write('#Seq\t#Site\t#Type:Proportion\n')
        for c, p, i in sites :
            tag, state = states[i]
            for n, ss in zip(names, state) :
                fout.write( '{0}\t{1}\t{2}\t{3}\n'.format(c, p, n, '\t'.join([ '{0}:{1:.5f}'.format(t, s) for t, s in zip(tag, ss)]) ))


def read_states(fname) :
    names, ss, sites = [], {}, []
    with uopen(fname) as fin :
        names = fin.readline().strip().split('\t')[2:]
        for line in fin :
            seq, site, snp_str = line.strip().split('\t', 2)
            if snp_str not in ss :
                ss[snp_str] = len(ss)
            sites.append([seq, int(site), ss[snp_str]])
    states = []
    for s, id in sorted(ss.items(), key=lambda x:x[1]) :
        states.append(s.split('\t'))
    return names, np.array(states), sites


def add_args(a) :
    parser = argparse.ArgumentParser(description='''
EToKi phylo runs to:
(1) Generate SNP matrix from alignment (-t matrix)
(2) Calculate ML phylogeny from SNP matrix using RAxML (-t phylogeny)
(3) Workout the nucleotide sequences of internal nodes in the tree using ML estimation (-t ancestral or -t ancestral_proportion for ratio frequencies)
(4) Place mutations onto branches of the tree (-t mutation)
''', formatter_class=argparse.RawTextHelpFormatter)
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
    parser.add_argument('--n_proc', '-n', help='Number of processes. Default: 7. ', type=int, default=7)

    args = parser.parse_args(a)

    args.tasks = dict(
        all = 'matrix,phylogeny,ancestral,mutation',
        aln2phy = 'matrix,phylogeny',
        snp2anc = 'phylogeny,ancestral',
        mat2mut = 'ancestral,mutation',
    ).get(args.tasks, args.tasks).split(',')

    return args

def infer_ancestral2(data) :
    state, branches, n_node, infer = data
    transitions = {}
    tag, code = np.unique(state, return_inverse=True)
    n_state = tag.size
    missing = np.where(tag == '-')[0]
    if missing.size > 0 :
        tag, n_state = tag[tag != '-'], n_state - 1
        code[code == missing[0]] = -1
        code[code  > missing[0]] = code[code  > missing[0]] - 1
    if len(tag) == 0 :
        tag = np.array(['-'])
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
        return [tag, alpha]
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
        return tag[r]

def infer_ancestral(tree, names, snps, sites, infer='margin', rescale=1.0) :
    tree = Tree(tree, format=1)
    node_names = {}
    for id,branch in enumerate(tree.traverse('postorder')) :
        digit = ''
        if not branch.is_leaf() :
            try :
                float(branch.name)
                digit = branch.name
                branch.name = ''
            except :
                pass
        if branch.name == '' :
            if digit == '' :
                branch.name = 'N_' + str(len(node_names))
            else :
                branch.name = 'N_' + str(len(node_names)) + '__{0}'.format(digit)
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
    #retvalue = [infer_ancestral2([state, branches, n_node, infer]) for state in states]
    try :
        retvalue = pool.map(infer_ancestral2, [[state, branches, n_node, infer] for state in states])
    except :
        pool = Pool(5)
        retvalue = pool.map(infer_ancestral2, [[state, branches, n_node, infer] for state in states])
    return tree, [ k for k, v in sorted(node_names.items(), key=lambda x:x[1])], retvalue

def phylo(args) :
    args = add_args(args)
    global pool
    pool = Pool(args.n_proc)
    
    if 'matrix' in args.tasks :
        assert os.path.isfile( args.alignment )
        args.snp = xFasta2Matrix( args.prefix, args.alignment )

    assert os.path.isfile( args.snp )
    names, sites, snps, seqLens, missing = read_matrix(args.snp)
    if len(names) < 4 :
        raise ValueError('Taxa too few.')
    snp_list = sorted([[info[0], int(math.ceil(info[2])), list(line), info[1]] for line, info in snps.items() ])


    # build tree
    if 'phylogeny' in args.tasks :
        phy, weights, asc = write_phylip(args.prefix+'.tre', names, snp_list)
        if phy != '' :
            args.tree = run_raxml(args.prefix +'.tre', phy, weights, asc, 'CAT', args.n_proc)
        else :
            args.tree = args.prefix + '.tre'
            with open(args.tree, 'w') as fout :
                fout.write('({0}:0.0);'.format(':0.0,'.join(names)))
        args.tree = get_root(args.prefix, args.tree)
    elif 'ancestral' in args.tasks or 'ancestral_proportion' in args.tasks :
        tree = Tree(args.tree, format=1)

    # map snp
    if 'ancestral' in args.tasks :
        final_tree, node_names, states = infer_ancestral(args.tree, names, snp_list, sites, infer='viterbi')
        states = np.array(states)
        final_tree.write(format=1, outfile=args.prefix + '.labelled.nwk')
        write_states(args.prefix+'.ancestral_states.gz', node_names, states, sites, seqLens, missing)
    elif 'mutation' in args.tasks :
        final_tree = Tree(args.tree, format=1)
        node_names, states, sites = read_states(args.ancestral)

    if 'ancestral_proportion' in args.tasks :
        final_tree, node_names, states = infer_ancestral(args.tree, names, snp_list, sites, infer='margin')
        final_tree.write(format=1, outfile=args.prefix + '.labelled.nwk')
        write_ancestral_proportion(args.prefix+'.ancestral_proportion.gz', node_names, states, sites, seqLens, missing)

    if 'mutation' in args.tasks :
        mutations = get_mut(final_tree, node_names, states, sites)
        with uopen(args.prefix + '.mutations.gz', 'w') as fout :
            for sl in seqLens :
                fout.write('## Sequence_length: {0} {1}\n'.format(*sl))
            for ms in missing :
                fout.write('## Missing_region: {0} {1} {2}\n'.format(*ms))
            
            fout.write('#Node\t#Seq\t#Site\t#Homoplasy\t#Mutation\n')
            for mut in mutations :
                fout.write('\t'.join([str(m) for m in mut]) + '\n')

pool = None
if __name__ == '__main__' :
    phylo(sys.argv[1:])

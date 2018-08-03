from ete3 import Tree
import numpy as np, pandas as pd, sys, os, re
try:
    from configure import uopen
except :
    from .configure import uopen


def parseTree(tre) :
    leafNames = set(tre.get_leaf_names())
    descendents, solo, duo = {}, {}, {}
    for nodeId, node in enumerate(tre.traverse('postorder')) :
        node.id = nodeId
        if node.is_leaf() :
            descendents[node.id] = [node.name]
            solo[node.name] = []
        else :
            children = node.get_children()
            descendents[node.id] = sorted([d for child in children for d in descendents[child.id]])
            if not node.name :
                node.name = ','.join(sorted(descendents[node.id]))

            for i1, c1 in enumerate(children) :
                for i2 in range(i1+1, len(children)) :
                    c2 = children[i2]
                    for d1 in descendents[c1.id] :
                        for d2 in descendents[c2.id] :
                            key = tuple(sorted([d1, d2]))
                            duo[key] = solo[d1] + solo[d2]
        node.bipartition = ':'.join([ ','.join(n) for n in sorted([sorted(descendents[node.id]), sorted(leafNames - set(descendents[node.id]))])])
        
        if node.up :
            for n in descendents[node.id] :
                solo[n].append(node.name)
        
    return tre, descendents, duo

mode = ('branch', 'pairwse')[0]
if __name__ == '__main__' :
    treeFile, recFile = sys.argv[1:]
    # read tree
    tre = Tree(treeFile, format=1)
    tre, descendents, duo = parseTree(tre)
    
    # read rec
    rec = []
    with uopen(recFile) as fin :
        for line in fin :
            part = line.strip().split('\t')
            if part[0] == 'Importation' :
                if float(part[-1]) >= 0.9 :
                    rec.append([part[1], None, int(part[3]), int(part[4])])
            elif part[0] == 'SEQUENCE' :
                n = re.findall(r'->(.+)";neg', part[-1])[0]
                rec.append([n, None, int(part[3]), int(part[4])])
            elif len(part) == 4 :
                try :
                    rec.append([part[2][1:-1], None, int(part[0]), int(part[1])])
                except :
                    pass
            elif len(part) == 3 :
                try :
                    rec.append([part[0], None, int(part[1]), int(part[2])])
                except :
                    pass
    recomb = {}
    for br, c, s, e in rec :
        if br not in recomb :
            recomb[br] = [[s, e]]
        else :
            recomb[br].append([s, e])

    tipRec = {}
    for tips, path in sorted(duo.items()) :
        rec = []
        for r in sorted([r for br in path for r in recomb.get(br, [])]) :
            if len(rec) and rec[-1][1] >= r[0] :
                if rec[-1][1] < r[1] :
                    rec[-1][1] = r[1]
            else :
                rec.append(r)
        tipRec[tips] = rec
        for r in rec :
            print('{0}\t{1}\t{2}\t{3}'.format(tips[0], tips[1], r[0], r[1]))

    ## load SNP matrix
    #with uopen(snpFile) as fin :
        #for line in fin :
            #if line.startswith('##') : 
                #continue
            #else :
                #header = line.strip().split('\t')
                #header = { h:i for i, h in enumerate(header) }
                #mat = pd.read_csv(fin, sep='\t', header=None, dtype=str).values.T
                #break
    #tipMut = {}
    #for (t1, t2) in tipRec :
        #snps = mat[1, (mat[header[t1]] != mat[header[t2]]) & (mat[header[t1]] != '-') & (mat[header[t1]] != '-')].astype(int)
        #tipMut[(t1, t2)] = snps

import argparse, tempfile, glob, os, subprocess, sys, shutil
try:
    from configure import externals, uopen, xrange, logger
except :
    from .configure import externals, uopen, xrange, logger

def clust(argv) :
    parser = argparse.ArgumentParser(description='Get clusters and exemplar sequences of the clusters. Run mmseqs linclust iteratively. ')
    parser.add_argument('-i', '--input', help='[INPUT; REQUIRED] name of the gene fasta file.', required=True)
    parser.add_argument('-p', '--prefix', help='[OUTPUT; REQUIRED] prefix of the outputs.', required=True)
    parser.add_argument('-d', '--identity', help='[PARAM; DEFAULT: 0.9] minimum intra-cluster identity.', default=0.9, type=float)
    parser.add_argument('-c', '--coverage', help='[PARAM; DEFAULT: 0.9] minimum intra-cluster coverage.', default=0.9, type=float)
    parser.add_argument('-t', '--n_thread', help='[PARAM; DEFAULT: 8]   number of threads to use.', default=8, type=int)
    args = parser.parse_args(argv)
    exemplar, clust = getClust(args.prefix, args.input, args.__dict__)
    logger('Exemplar sequences in {0}'.format(exemplar))
    logger('Clusters in {0}'.format(clust))
    return exemplar, clust
def getClust(prefix, genes, params) :
    groups = {}
    dirPath = tempfile.mkdtemp(prefix='NS_', dir='.')
    try:
        geneFile = genes
        seqDb = os.path.join(dirPath, 'seq.db')
        tmpDb = os.path.join(dirPath, 'tmp')
        lcDb = os.path.join(dirPath, 'seq.lc')
        tabFile = os.path.join(dirPath, 'clust.tab')
        refFile = os.path.join(dirPath, 'seq.ref')
        
        nRef = 999999999999999
        for ite in xrange(5) :
            if os.path.isdir(tmpDb) :
                shutil.rmtree(tmpDb)
            os.makedirs(tmpDb)
            if os.path.isfile(seqDb) :
                list(map(os.unlink, glob.glob(seqDb + '*')))
            if os.path.isfile(lcDb) :
                list(map(os.unlink, glob.glob(lcDb + '*')))

            subprocess.Popen('{0} createdb {2} {1} -v 0'.format(externals['mmseqs'], seqDb, geneFile).split()).communicate()
            subprocess.Popen('{0} linclust {1} {2} {3} --min-seq-id {4} -c {5} --threads {6} -v 0'.format( \
                externals['mmseqs'], seqDb, lcDb, tmpDb, params['identity'], params['coverage'], params['n_thread']).split(), stdout=subprocess.PIPE).communicate()
            subprocess.Popen('{0} createtsv {1} {1} {2} {3}'.format(\
                externals['mmseqs'], seqDb, lcDb, tabFile).split(), stdout = subprocess.PIPE).communicate()
            with open(tabFile) as fin :
                for line in fin :
                    part = line.strip().split()
                    groups[part[1]] = part[0]
            tmp = []
            with open(geneFile) as fin :
                toWrite, used_grps = False, {None:1}
                for line in fin :
                    if line.startswith('>') :
                        name = line[1:].strip().split()[0]
                        grp = groups.get(name, None)
                        toWrite = False if grp in used_grps else True
                        if toWrite :
                            used_grps[grp] = name
                    if toWrite :
                        tmp.append(line)
                for gene, grp in groups.items() :
                    if grp in used_grps :
                        groups[gene] = used_grps[grp]
            with open(refFile, 'w') as fout :
                for line in tmp :
                    fout.write(line)
            if nRef <= len(used_grps) :
                break
            nRef = len(used_grps)
            geneFile = refFile
        shutil.copy2(refFile, '{0}.clust.exemplar'.format(prefix))
    finally :
        shutil.rmtree(dirPath)
    with open('{0}.clust.tab'.format(prefix), 'w') as fout :
        for gene, grp in sorted(groups.items()) :
            g = gene
            while g != grp :
                g, grp = grp, groups[grp]
            groups[gene] = grp
            fout.write('{0}\t{1}\n'.format(gene, grp))
    
    return '{0}.clust.exemplar'.format(prefix), '{0}.clust.tab'.format(prefix)

if __name__ == '__main__' :
    clust(sys.argv[1:])

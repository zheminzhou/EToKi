import os, sys, shutil, pandas as pd
try :
    from configure import externals, logger, uopen, xrange, StringIO, get_md5
except :
    from .configure import externals, logger, uopen, xrange, StringIO, get_md5
import subprocess, tempfile, time

mmseqs = externals['mmseqs']
minimap2 = externals['minimap2']

def minimapFilter(sourceFna, targetFna, targetFiltFna, max_iden, min_iden, coverage, paralog, relaxEnd, orderedLoci) :
    p = subprocess.Popen('{0} -ct8 -k13 -w5 -A2 -B4 -O8,16 -E2,1 -r50 -p.2 -N500 -f2000,10000 -n1 -m19 -s40 -g200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(minimap2, sourceFna, targetFna).split(), stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines=True)
    tooClose, goodCandidates, crossLoci = {}, {}, {}
    for line in p.stdout :
        part = line.strip().split('\t')
        q, r = part[0], part[5]
        qLoc, rLoc = q.rsplit('_', 1)[0], r.rsplit('_', 1)[0]
        if q == r :
            goodCandidates[r] = 1.
            continue
        elif rLoc == qLoc :
            if part[4] == '-' :
                continue
            tl, ts, te = [int(p) for p in part[1:4]]
            rl, rs, re, ri = [int(p) for p in part[6:10]]
            if relaxEnd or (ts != rs or tl-te != rl - re) :
                continue
            iden, cov = float(ri)/(re-rs), float(re-rs)/rl
            if cov > coverage :
                if iden >= max_iden :
                    tooClose[part[0]] = 1
                elif iden >= min_iden :
                    goodCandidates[part[0]] = max(goodCandidates.get(part[0], 0), iden)
        elif orderedLoci[qLoc] > orderedLoci[rLoc] :
            rl, rs, re, ri = [int(p) for p in part[6:10]]
            iden, cov = float(ri)/(re-rs), float(re-rs)/rl
            if cov > coverage and iden >= min_iden and crossLoci.get(part[0], 0) < iden :
                crossLoci[part[0]] = iden
    with open(targetFna) as fin, open(targetFiltFna+'.fas', 'w') as fout :
        writable = False
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                writable = False if (name in tooClose or crossLoci.get(name, 0) > 1. - paralog) else True
            if writable :
                fout.write(line)
    return targetFiltFna, goodCandidates, crossLoci

def buildReference(targets, sources, max_iden=0.9,  min_iden=0.6, coverage=0.7, paralog=0.1, relaxEnd=False) :
    orderedLoci = { t['fieldname']:i for i, t in reversed(list(enumerate(sources))) }
    refsets = []
    dirPath = tempfile.mkdtemp(prefix='NS_', dir='.')
    try:
        tmpDir = os.path.join(dirPath, 'tmp')
        sourceFna = os.path.join(dirPath, 'sourceFna')
        sourceFaa = os.path.join(dirPath, 'sourceFaa')
        targetFna = os.path.join(dirPath, 'targetFna')
        targetFiltFna = os.path.join(dirPath, 'targetFiltFna')
        targetFiltFaa = os.path.join(dirPath, 'targetFiltFaa')
        targetClsFna = os.path.join(dirPath, 'targetClsFna')
        targetResFna = os.path.join(dirPath, 'targetResFna')
        alnFaa = os.path.join(dirPath, 'alnFaa')
        
        with open(sourceFna+'.fas', 'w') as fout :
            fout.write('\n'.join(['>{fieldname}_{value_id}\n{value}'.format(**s) for s in sources]))
        with open(targetFna+'.fas', 'w') as fout :
            fout.write('\n'.join(['>{fieldname}_{value_id}\n{value}'.format(**t) for t in targets]))
        targetFiltFna, goodCandidates, crossSites = minimapFilter(sourceFna+'.fas', targetFna+'.fas', targetFiltFna, max_iden, min_iden, coverage, paralog, relaxEnd, orderedLoci)
        logger('identifed {0} good exemplar alleles after nucleic search'.format(len(goodCandidates)))
        
        subprocess.Popen('{0} createdb {1}.fas {1} --dont-split-seq-by-len'.format(mmseqs, sourceFna).split(), stdout = subprocess.PIPE).communicate()
        subprocess.Popen('{0} createdb {1}.fas {1} --dont-split-seq-by-len'.format(mmseqs, targetFiltFna).split(), stdout = subprocess.PIPE).communicate()
        subprocess.Popen('{0} translatenucs {1} {2}'.format(mmseqs, sourceFna, sourceFaa).split(), stdout = subprocess.PIPE).communicate()
        subprocess.Popen('{0} translatenucs {1} {2}'.format(mmseqs, targetFiltFna, targetFiltFaa).split(), stdout = subprocess.PIPE).communicate()
        for ite in range(9) :
            if os.path.isdir(tmpDir) :
                shutil.rmtree(tmpDir)
            p=subprocess.Popen('{0} search {2} {1} {3} {4} -c {6} --min-seq-id {5} --threads {t}'.format(mmseqs, sourceFaa, targetFiltFaa, alnFaa, tmpDir, min_iden, coverage, t=9-4*int(ite/3) ).split(), stdout = subprocess.PIPE)
            p.communicate()
            if p.returncode == 0 :
                break
            time.sleep(1)
        subprocess.Popen('{0} convertalis {2} {1} {3} {3}.tab --format-mode 2 --threads 8'.format(mmseqs, sourceFaa, targetFiltFaa, alnFaa).split(), stdout = subprocess.PIPE).communicate()
        
        with open(alnFaa + '.tab') as fin :
            for line in fin :
                part = line.strip().split('\t')
                qLoc, rLoc = part[0].rsplit('_', 1)[0], part[1].rsplit('_', 1)[0]
                if qLoc == rLoc:
                    if relaxEnd or (int(part[8]) == int(part[6]) and int(part[12]) - int(part[7]) == int(part[13]) - int(part[9])) :
                        goodCandidates[part[0]] = max(goodCandidates.get(part[0], 0), float(part[2]))
                elif orderedLoci[qLoc] > orderedLoci[rLoc] and crossSites.get(part[0], 0) < float(part[2]) :
                    crossSites[part[0]] = float(part[2])
        logger('identifed a total of {0} good exemplar alleles after amino search'.format(len(goodCandidates)))
        
        nLoci = len(orderedLoci)
        for s in sources :
            key = '{0}_{1}'.format(s['fieldname'], s['value_id'])
            if crossSites.get(key, 0) > 1-paralog :
                orderedLoci.pop(s['fieldname'], None)
                #logger(key)
        if nLoci > len(orderedLoci) :
            logger('Total of {0} loci are not suitable for MLST scheme [due to paralog setting]. There are {1} left'.format(nLoci - len(orderedLoci), len(orderedLoci)))

        with open(targetFna+'.fas') as fin, open(targetClsFna+'.fas', 'w') as fout :
            writable = False
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    locus, id = name.rsplit('_', 1)
                    writable = True if locus in orderedLoci and goodCandidates.get(name, 0) - crossSites.get(name, 0) > paralog else False
                if writable :
                    fout.write(line)
        subprocess.Popen('{0} createdb {1}.fas {1} --dont-split-seq-by-len'.format(mmseqs, targetClsFna).split(), stdout = subprocess.PIPE).communicate()
        for ite in range(9) :
            if os.path.isdir(tmpDir) :
                shutil.rmtree(tmpDir)            
            p = subprocess.Popen('{0} cluster {1} {2} {3} -c {5} --min-seq-id {4} --threads {t}'.format(mmseqs, targetClsFna, targetResFna, tmpDir, max_iden, max_iden, t=9-4*int(ite/3)).split(), stdout = subprocess.PIPE)
            p.communicate()
            if p.returncode == 0 :
                break
            time.sleep(1)
        subprocess.Popen('{0} createtsv {1} {2} {3} {3}.tab'.format(mmseqs, targetClsFna, targetClsFna, targetResFna).split(), stdout = subprocess.PIPE).communicate()
        goodCandidates = {}
        with open(targetResFna + '.tab') as fin :
            for line in fin :
                goodCandidates[line.split('\t', 1)[0]] = 1
        logger('There are {0} good exemplar alleles left after final clustering'.format(len(goodCandidates)))

        with open(targetClsFna+'.fas') as fin:
            writable = False
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    writable = True if name in goodCandidates else False
                if writable :
                    refsets.append(line.strip())
    except :
        pass
    finally:
        shutil.rmtree(dirPath)
        return '\n'.join(refsets)

def readFasta(fastaText) :
    sequence = []
    for line in fastaText :
        if line.startswith('>') :
            name = line[1:].strip().split()[0]
            part = name.rsplit('_', 1)
            if len(part) > 1 :
                sequence.append({'fieldname':part[0], 'value_id':part[1], 'value':[]})
            else :
                sequence.append({'fieldname':name, 'value_id':None, 'value':[]})
        elif len(line) > 0 :
            sequence[-1]['value'].extend(line.strip().split())
    for s in sequence :
        s['value'] = ''.join(s['value'])
    return sequence


def MLSTdb(args) :
    params = getParams(args)
    database, refset, alleleFasta, refstrain, max_iden, min_iden, coverage, paralog, relaxEnd=params['database'], params['refset'], params['alleleFasta'], params['refstrain'], params['max_iden'], params['min_iden'], params['coverage'], params['paralog'], params['relaxEnd']
    if os.path.isfile(alleleFasta) :
        alleles = readFasta(uopen(alleleFasta))
    else :
        alleles = readFasta(StringIO(alleleFasta))
    alleles = [allele for allele in alleles \
                   if allele['value_id'].isdigit() and int(allele['value_id']) > 0]
    refAlleles = ''
    if refset is not None :
        if refstrain :
            if os.path.isfile(refstrain) :
                references = readFasta(uopen(refstrain))
            else :
                references = readFasta(StringIO(refstrain))
        else :
            loci, references = {}, []
            for allele in alleles :
                if allele['fieldname'] not in loci :
                    loci[allele['fieldname']] = 1
                    references.append(allele)
        
        refAlleles = buildReference(alleles, references, max_iden, min_iden, coverage, paralog, relaxEnd)
        if refset :
            with open(str(refset), 'w') as fout :
                fout.write(refAlleles + '\n')
        logger('A file of reference alleles has been generated:  {0}'.format(refset))
    if database :
        conversion = [[], []]
        with open(database, 'w') as fout :
            for allele in alleles :
                conversion[0].append(get_md5(allele['value']))
                conversion[1].append([allele['fieldname'], int(allele['value_id'])])
        conversion = pd.DataFrame(conversion[1], index=conversion[0])
        conversion.to_csv(database, header=False)
        logger('A lookup table of all alleles has been generated:  {0}'.format(database))
    return refAlleles

def getParams(args) :
    import argparse
    parser = argparse.ArgumentParser(description='MLSTdb. Create reference sets of alleles for nomenclature. ')
    parser.add_argument('-i', '--input', dest='alleleFasta', help='[REQUIRED] A single file contains all known alleles in a MLST scheme. ', required=True)
    parser.add_argument('-r', '--refset',                    help='[DEFAULT: No ref allele] Output - Reference alleles used for MLSType. ', default=None)
    parser.add_argument('-d', '--database',                  help='[DEFAULT: No allele DB] Output - A lookup table of all alleles. ', default=None)
    parser.add_argument('-s', '--refstrain',                 help='[DEFAULT: None] A single file contains alleles from the reference genome. ', default=None)
    parser.add_argument('-x', '--max_iden',                  help='[DEFAULT: 0.9 ] Maximum identities between resulting refAlleles. ', type=float, default=0.9)
    parser.add_argument('-m', '--min_iden',                  help='[DEFAULT: 0.4 ] Minimum identities between refstrain and resulting refAlleles. ', type=float, default=0.4)
    parser.add_argument('-p', '--paralog',                   help='[DEFAULT: 0.1 ] Minimum differences between difference loci. ', type=float, default=0.2)
    parser.add_argument('-c', '--coverage',                  help='[DEFAULT: 0.7 ] Proportion of aligned regions between alleles. ', type=float, default=0.7)
    parser.add_argument('-e', '--relaxEnd',                  help='[DEFAULT: False ] Allow changed ends (for pubmlst). ', action='store_true', default=False)

    return parser.parse_args(args).__dict__
    
if __name__ == '__main__' :
    sys.stdout.write(MLSTdb(sys.argv[1:])+'\n')

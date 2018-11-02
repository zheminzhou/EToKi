import os, sys
from configure import externals, logger
import subprocess, tempfile, shutil, time

mmseqs = externals['mmseqs']
minimap2 = externals['minimap2']

def minimapFilter(sourceFna, targetFna, targetFiltFna, max_iden, min_iden, coverage) :
    p = subprocess.Popen('{0} -ct8 -k13 -w5 -A2 -B4 -O8,16 -E2,1 -r50 -p.2 -N500 -f2000,10000 -n1 -m19 -s40 -g200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(minimap2, sourceFna, targetFna).split(), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    tooClose = {}
    goodCandidates = {}
    for line in p.stdout :

        part = line.strip().split('\t')
        l, s, e, i = [float(p) for p in part[6:10]]
        iden, cov = i/(e-s), (e-s)/l
        if cov > coverage :
            if iden >= max_iden :
                tooClose[part[0]] = 1
            elif iden >= min_iden :
                goodCandidates[part[0]] = 1
    with open(targetFna) as fin, open(targetFiltFna+'.fas', 'w') as fout :
        writable = False
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                writable = False if name in tooClose else True
            if writable :
                fout.write(line)
    return targetFiltFna, goodCandidates

def selfReference(targets, max_iden, coverage) :
    dirPath = tempfile.mkdtemp(prefix='NS_', dir='.')
    try:
        tmpDir = os.path.join(dirPath, 'tmp')
        targetFna = os.path.join(dirPath, 'targetFna')
        targetClsFna = os.path.join(dirPath, 'targetClsFna')
        targetClsFaa = os.path.join(dirPath, 'targetClsFaa')
        targetResFaa = os.path.join(dirPath, 'targetResFaa')
        with open(targetFna+'.fas', 'w') as fout :
            fout.write('\n'.join(['>{0}\n{1}'.format(*t) for t in targets]))
        subprocess.Popen('{0} createdb {1}.fas {1} --dont-split-seq-by-len'.format(mmseqs, targetFna).split(), stdout = subprocess.PIPE).communicate()
        for ite in range(9) :
            if os.path.isdir(tmpDir) :
                shutil.rmtree(tmpDir)
            p = subprocess.Popen('{0} cluster {1} {2} {3} -c {5} --min-seq-id {4} -e 0.01 --threads 8'.format(mmseqs, targetFna, targetClsFna, tmpDir, max_iden, coverage).split(), stdout = subprocess.PIPE)
            p.communicate()
            if p.returncode == 0 :
                break
            time.sleep(1)
        
        subprocess.Popen('{0} createtsv {1} {2} {3} {3}.tab --threads 8'.format(mmseqs, targetFna, targetFna, targetClsFna).split(), stdout = subprocess.PIPE).communicate()
        goodCandidates = {}
        with open(targetClsFna + '.tab') as fin :
            for line in fin :
                goodCandidates[line.split('\t', 1)[0]] = 1

        with open(targetFna+'.fas') as fin, open(targetClsFna+'.fas', 'w') as fout :
            writable = False
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    writable = True if name in goodCandidates else False
                if writable :
                    fout.write(line)
        subprocess.Popen('{0} createdb {1}.fas {1} --dont-split-seq-by-len'.format(mmseqs, targetClsFna).split(), stdout = subprocess.PIPE).communicate()
        subprocess.Popen('{0} translatenucs {1} {2} --threads 8'.format(mmseqs, targetClsFna, targetClsFaa).split(), stdout = subprocess.PIPE).communicate()
        for ite in range(9) :
            if os.path.isdir(tmpDir) :
                shutil.rmtree(tmpDir)
            p = subprocess.Popen('{0} cluster {1} {2} {3} -c {5} --min-seq-id {4} -e 0.01 --threads 8'.format(mmseqs, targetClsFaa, targetResFaa, tmpDir, max_iden, coverage).split(), stdout = subprocess.PIPE)
            p.communicate()
            if p.returncode == 0 :
                break
            time.sleep(1)
        
        subprocess.Popen('{0} createtsv {1} {2} {3} {3}.tab --threads 8'.format(mmseqs, targetClsFaa, targetClsFaa, targetResFaa).split(), stdout = subprocess.PIPE).communicate()
        goodCandidates = {}
        with open(targetResFaa + '.tab') as fin :
            for line in fin :
                goodCandidates[line.split('\t', 1)[0]] = 1
                
        refsets = []
        with open(targetClsFna+'.fas') as fin, open(targetResFaa+'.fas', 'w') as fout :
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
    
def buildReference(targets, sources=None, max_iden=0.95,  min_iden=0.5, coverage=0.9) :
    if sources is None :
        return selfReference(targets, max_iden, coverage)
    logger('Running {0}\n'.format(sources[0][0]))
    refsets = ['>{0}\n{1}'.format(*s) for s in sources]
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
            fout.write('\n'.join(['>{0}\n{1}'.format(*s) for s in sources]))
        with open(targetFna+'.fas', 'w') as fout :
            fout.write('\n'.join(['>{0}\n{1}'.format(*t) for t in targets]))
        targetFiltFna, goodCandidates = minimapFilter(sourceFna+'.fas', targetFna+'.fas', targetFiltFna, max_iden, min_iden, coverage)
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
                goodCandidates[line.split('\t', 1)[0]] = 1
        with open(targetFiltFna+'.fas') as fin, open(targetClsFna+'.fas', 'w') as fout :
            writable = False
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip().split()[0]
                    writable = True if name in goodCandidates else False
                if writable :
                    fout.write(line)
        subprocess.Popen('{0} createdb {1}.fas {1} --dont-split-seq-by-len'.format(mmseqs, targetClsFna).split(), stdout = subprocess.PIPE).communicate()
        for ite in range(9) :
            if os.path.isdir(tmpDir) :
                shutil.rmtree(tmpDir)            
            p = subprocess.Popen('{0} cluster {1} {2} {3} -c {5} --min-seq-id {4} --threads {t}'.format(mmseqs, targetClsFna, targetResFna, tmpDir, max_iden, coverage, t=9-4*int(ite/3)).split(), stdout = subprocess.PIPE)
            p.communicate()
            if p.returncode == 0 :
                break
            time.sleep(1)
        subprocess.Popen('{0} createtsv {1} {2} {3} {3}.tab'.format(mmseqs, targetClsFna, targetClsFna, targetResFna).split(), stdout = subprocess.PIPE).communicate()
        goodCandidates = {}
        with open(targetResFna + '.tab') as fin :
            for line in fin :
                goodCandidates[line.split('\t', 1)[0]] = 1

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


def MLSTdb(refAllele, alleleFasta, refstrain, max_iden, min_iden, coverage) :
    with open(alleleFasta) as fin :
        alleles = readFasta(fin)
    loci = {}
    for allele in alleles :
        locus = allele['fieldname']
        if locus not in loci :
            loci[locus] = [('{0}_{1}'.format(locus, allele['value_id']), allele['value'], )]
        else :
            loci[locus].append(('{0}_{1}'.format(locus, allele['value_id']), allele['value'], ))
    references = {}
    if refstrain is not None :
        with open(refstrain) as fin :
            alleles = readFasta(fin)
        for allele in alleles :
            locus = allele['fieldname']
            if locus not in references :
                references[locus] = [('{0}_{1}'.format(locus, allele['value_id']), allele['value'], )]
            else :
                references[locus].append(('{0}_{1}'.format(locus, allele['value_id']), allele['value'], ))
    with open(refAllele, 'w') as fout :
        for locus, alleles in loci.items() :
            refAlleles = buildReference(alleles, references.get(locus, None), max_iden, min_iden, coverage)
            fout.write(refAlleles + '\n')

def getParams() :
    import argparse
    parser = argparse.ArgumentParser(description='MLSTdb. Create reference sets of alleles for nomenclature. ')
    parser.add_argument('-i', '--input', dest='alleleFasta', help='[REQUIRED] A single file contains all known alleles in a MLST scheme. ', required=True)
    parser.add_argument('-o', '--output', dest='refAllele',  help='[REQUIRED] Output - A single file used for MLSType. ', required=True)
    parser.add_argument('-r', '--refstrain',                 help='[DEFAULT: None] A single file contains alleles from the reference genome. ', default=None)
    parser.add_argument('-x', '--max_iden',                  help='[DEFAULT: 0.9 ] Maximum identities between resulting refAlleles. ', type=float, default=0.9)
    parser.add_argument('-m', '--min_iden',                  help='[DEFAULT: 0.6 ] Minimum identities between refstrain and resulting refAlleles. ', type=float, default=0.6)
    parser.add_argument('-c', '--coverage',                  help='[DEFAULT: 0.9 ] Proportion of aligned regions between alleles. ', type=float, default=0.9)
    
    return parser.parse_args().__dict__
    
if __name__ == '__main__' :
    MLSTdb(**getParams())
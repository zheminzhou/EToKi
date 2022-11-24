import os, sys, subprocess, numpy as np, pandas as pd, argparse, glob, gzip, io, re
from datetime import datetime

if sys.version_info[0] < 3:
    from collections import OrderedDict
    xrange = xrange
    from cStringIO import StringIO
    asc2int = np.uint8
else :
    from _collections import OrderedDict
    from io import StringIO
    xrange = range
    asc2int = np.uint32

import hashlib, uuid
def get_md5(value, dtype=str) :
    m = hashlib.md5(str(value).encode()).hexdigest()
    if dtype == str :
        return str(uuid.UUID(m))
    else :
        return int(m, 16)


# * is designated as U; index is : (ord(r)-65)*32 + ord(q)-65
blosum62 = np.array([  4., -2.,  0., -2., -1., -2.,  0., -2., -1.,  0., -1., -1., -1., -2.,  0., -1., -1., -1.,  1.,  0., -4.,  0.,
                       -3.,  0., -2., -1.,  0.,  0.,  0.,  0.,  0.,  0., -2., 4., -3.,  4.,  1., -3., -1.,  0., -3.,  0.,  0., -4.,
                       -3.,  3.,  0., -2.,  0., -1.,  0., -1., -4., -3., -4., -1., -3.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -3.,
                       9., -3., -4., -2., -3., -3., -1.,  0., -3., -1., -1., -3.,  0., -3., -3., -3., -1., -1., -4., -1., -2., -2.,
                       -2., -3.,  0.,  0.,  0.,  0.,  0.,  0., -2.,  4., -3., 6.,  2., -3., -1., -1., -3.,  0., -1., -4., -3.,  1.,
                       0., -1.,  0., -2.,  0., -1., -4., -3., -4., -1., -3., 1.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  1., -4.,  2.,
                       5., -3., -2.,  0., -3.,  0.,  1., -3., -2.,  0.,  0., -1.,  2.,  0.,  0., -1., -4., -2., -3., -1., -2.,  4.,
                       0.,  0.,  0.,  0.,  0.,  0., -2., -3., -2., -3., -3., 6., -3., -1.,  0.,  0., -3.,  0.,  0., -3.,  0., -4.,
                       -3., -3., -2., -2., -4., -1.,  1., -1.,  3., -3.,  0., 0.,  0.,  0.,  0.,  0.,  0., -1., -3., -1., -2., -3.,
                       6., -2., -4.,  0., -2., -4., -3.,  0.,  0., -2., -2., -2.,  0., -2., -4., -3., -2., -1., -3., -2.,  0.,  0.,
                       0.,  0.,  0.,  0., -2.,  0., -3., -1.,  0., -1., -2., 8., -3.,  0., -1., -3., -2.,  1.,  0., -2.,  0.,  0.,
                       -1., -2., -4., -3., -2., -1.,  2.,  0.,  0.,  0.,  0., 0.,  0.,  0., -1., -3., -1., -3., -3.,  0., -4., -3.,
                       4.,  0., -3.,  2.,  1., -3.,  0., -3., -3., -3., -2., -1., -4.,  3., -3., -1., -1., -3.,  0.,  0.,  0.,  0.,
                       0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                       0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0., -1.,  0., -3., -1.,  1., -3., -2., -1., -3.,  0.,
                       5., -2., -1.,  0.,  0., -1.,  1.,  2.,  0., -1., -4., -2., -3., -1., -2.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,
                       -1., -4., -1., -4., -3.,  0., -4., -3.,  2.,  0., -2., 4.,  2., -3.,  0., -3., -2., -2., -2., -1., -4.,  1.,
                       -2., -1., -1., -3.,  0.,  0.,  0.,  0.,  0.,  0., -1., -3., -1., -3., -2.,  0., -3., -2.,  1.,  0., -1.,  2.,
                       5., -2.,  0., -2.,  0., -1., -1., -1., -4.,  1., -1., -1., -1., -1.,  0.,  0.,  0.,  0.,  0.,  0., -2.,  3.,
                       -3.,  1.,  0., -3.,  0.,  1., -3.,  0.,  0., -3., -2., 6.,  0., -2.,  0.,  0.,  1.,  0., -4., -3., -4., -1.,
                       -2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                       0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  0.,  0., -1., -2., -3., -1.,
                       -1., -4., -2., -2., -3.,  0., -1., -3., -2., -2.,  0., 7., -1., -2., -1., -1., -4., -2., -4., -2., -3., -1.,
                       0.,  0.,  0.,  0.,  0.,  0., -1.,  0., -3.,  0.,  2., -3., -2.,  0., -3.,  0.,  1., -2.,  0.,  0.,  0., -1.,
                       5.,  1.,  0., -1., -4., -2., -2., -1., -1.,  3.,  0., 0.,  0.,  0.,  0.,  0., -1., -1., -3., -2.,  0., -3.,
                       -2.,  0., -3.,  0.,  2., -2., -1.,  0.,  0., -2.,  1., 5., -1., -1., -4., -3., -3., -1., -2.,  0.,  0.,  0.,
                       0.,  0.,  0.,  0.,  1.,  0., -1.,  0.,  0., -2.,  0., -1., -2.,  0.,  0., -2., -1.,  1.,  0., -1.,  0., -1.,
                       4.,  1., -4., -2., -3.,  0., -2.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  0., -1., -1., -1., -1., -2., -2., -2.,
                       -1.,  0., -1., -1., -1.,  0.,  0., -1., -1., -1.,  1., 5., -4.,  0., -2.,  0., -2., -1.,  0.,  0.,  0.,  0.,
                       0.,  0., -4., -4., -4., -4., -4., -4., -4., -4., -4., 0., -4., -4., -4., -4.,  0., -4., -4., -4., -4., -4.,
                       1., -4., -4., -4., -4., -4.,  0.,  0.,  0.,  0.,  0., 0.,  0., -3., -1., -3., -2., -1., -3., -3.,  3.,  0.,
                       -2.,  1.,  1., -3.,  0., -2., -2., -3., -2.,  0., -4., 4., -3., -1., -1., -2.,  0.,  0.,  0.,  0.,  0.,  0.,
                       -3., -4., -2., -4., -3.,  1., -2., -2., -3.,  0., -3., -2., -1., -4.,  0., -4., -2., -3., -3., -2., -4., -3.,
                       11., -2.,  2., -3.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -1., -2., -1., -1., -1., -1., -1., -1.,  0., -1., -1.,
                       -1., -1.,  0., -2., -1., -1.,  0.,  0., -4., -1., -2., -1., -1., -1.,  0.,  0.,  0.,  0.,  0.,  0., -2., -3.,
                       -2., -3., -2.,  3., -3.,  2., -1.,  0., -2., -1., -1., -2.,  0., -3., -1., -2., -2., -2., -4., -1.,  2., -1.,
                       7., -2.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  1., -3., 1.,  4., -3., -2.,  0., -3.,  0.,  1., -3., -1.,  0.,
                       0., -1.,  3.,  0.,  0., -1., -4., -2., -3., -1., -2., 4.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                       0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])


class uopen(object) :
    def __init__(self, fname, label='r') :
        self.fout = None
        if label.find('r')>=0 :
            self.fstream = subprocess.Popen([externals['pigz'], '-cd', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout if fname.lower().endswith('gz') else open(fname)
        elif label.find('w') >= 0 :
            self.fout = open(fname, 'wb')
            p = subprocess.Popen([externals['pigz']], stdin=subprocess.PIPE, stdout=self.fout, universal_newlines=True)
            self.fstream = p.stdin
            
        elif label.find('a') >= 0 :
            if sys.version.startswith('3') :
                self.fout = gzip.open(fname, 'ab')
                self.fstream = io.TextIOWrapper(self.fout, encoding='utf-8')
            else :
                self.fstream = gzip.open(fname, 'ab')
    def __enter__(self) :
        return self.fstream
    def __exit__(self, type, value, traceback) :
        self.fstream.close()
        if self.fout :
            self.fout.close()
        return 
    def __iter__(self) :
        return self.fstream
    def __next__(self) :
        return self.fstream


def readFasta(fasta, headOnly=False) :
    sequence = OrderedDict()
    with uopen(fasta) as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                sequence[name] = []
            elif len(line) > 0 and not line.startswith('#') and not headOnly :
                sequence[name].extend(line.strip().split())
    for s in sequence :
        sequence[s] = (''.join(sequence[s])).upper()
    return sequence

def readFastq(fastq) :
    sequence, qual = OrderedDict(), OrderedDict()
    with uopen(fastq) as fin :
        line = fin.readline()
        if not line.startswith('@') :
            sequence = readFasta(fastq)
            return sequence, OrderedDict( [n, re.sub(r'[^!]', 'I', re.sub(r'[^ACGTacgt]', '!', s))] for n, s in sequence.items() )
    with uopen(fastq) as fin :
        for lineId, line in enumerate(fin) :
            if lineId % 4 == 0 :
                name = line[1:].strip().split()[0]
                sequence[name] = []
                qual[name] = []
            elif lineId % 4 == 1 :
                sequence[name].extend(line.strip().split())
            elif lineId % 4 == 3 :
                qual[name].extend(line.strip().split())
    for s in sequence :
        sequence[s] = (''.join(sequence[s])).upper()
        qual[s] = ''.join(qual[s])
    return sequence, qual

complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
def rc(seq, missingValue='N') :
    return ''.join([complement.get(s, missingValue) for s in reversed(seq.upper())])


def rev_transeq(s, transl_table=11) :
    if transl_table == 4 :
        gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXY*YSSSSWCWCLFLF-'))
    else :
        gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXY*YSSSSXCWCLFLF-'))
    rev_seq = {}
    na = ['A', 'C', 'G', 'T', '-']
    for aa, codon in zip(*np.unique(gtable, return_index=True)) :
        rev_seq[aa] = ''.join([na[int(codon/16)], na[int((codon%16)/4)], na[codon%4]])
    return ''.join([ rev_seq.get(x, '---') for x in s ])

baseConv = np.repeat(-100, 255)
baseConv[(np.array(['-', 'A', 'C', 'G', 'T']).view(asc2int),)] = (-100000, 0, 1, 2, 3)
def transeq(seq, frame=7, transl_table=None, markStarts=False) :
    frames = {'F': [1,2,3],
              'R': [4,5,6],
              '7': [1,2,3,4,5,6]}.get( str(frame).upper() , None)
    if frames is None :
        frames = [int(f) for f in str(frame).split(',')]
    
    if transl_table == 4 :
        gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSWCWCLFLF-'))
    else :
        gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF-'))
    if markStarts :
        gtable[(np.array([46, 62]),)] = 'M'

    seqs = seq.items() if isinstance(seq, dict) else seq
    nFrame = (max(frames) > 3)
    trans_seq = []
    for n,s in seqs :
        s = baseConv[np.array(list(s.upper())).view(asc2int)]
        if nFrame :
            rs = (3 - s)[::-1]
            rs[rs >= 100] *= -1            
        sf = s.size % 3
        tseq = []
        for f in frames :
            codons = s[f-1:] if f <= 3 else rs[f-4:]
            if codons.size % 3 :
                codons = np.concatenate([codons, [-100]*(3 - codons.size % 3)])
            codons = codons.reshape(-1, 3)
            codon2 = np.sum(codons << [4, 2, 0], 1)
            codon2[codon2 < -50000] = 64
            codon2[codon2 < 0] = 50
            tseq.append(''.join(gtable[codon2].tolist()))
        trans_seq.append([n, tseq])
    return dict(trans_seq) if isinstance(seq, dict) else trans_seq


def logger(log, pipe=sys.stderr) :
    pipe.write('{0}\t{1}\n'.format(str(datetime.now()), log))
    pipe.flush()


def getExecutable(commands) :
    def check_sys_path(cmd) :
        for path in [''] + os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, cmd)
            if os.path.exists(exe_file):
                return exe_file
        return None
    try :
        cmd = check_sys_path(commands[-1])
        if not cmd  :
            return None
        commands[-1] = cmd
        if subprocess.Popen(commands+['-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait() <= 1 or \
           subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait() <= 1 or \
           subprocess.Popen(commands+['-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait() <= 1 :
            return commands
        else :
            return None
    except  :
        return None

def download_krakenDB() :
    curdir = os.path.abspath(os.curdir)
    moveTo = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'externals')
    os.chdir(moveTo)
    if not os.path.exists('minikraken2') :
        os.makedirs('minikraken2')
    os.chdir(os.path.join(moveTo, 'minikraken2'))
    minikraken_url = 'ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz'
    logger('Downloading minikraken2 from {0}. This might take a long time.'.format(minikraken_url))    
    subprocess.Popen('curl -Lo minikraken2_v2_8GB.tgz {0}'.format(minikraken_url).split(), stderr=subprocess.PIPE).communicate()
    logger('Unpackaging minikraken2.')
    subprocess.Popen('tar -xzf minikraken2_v2_8GB.tgz'.split()).communicate()
    subprocess.Popen('mv minikraken2_v2_8GB_*/* ./', shell=True).communicate()
    os.unlink('minikraken2_v2_8GB.tgz')
    
    os.chdir(curdir)


def install_externals() :
    curdir = os.path.abspath(os.curdir)
    moveTo = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'externals')
    os.chdir(moveTo)

    if not getExecutable(['java']) :
        logger('You have not installed Java runtime. Please install it first. ')
        sys.exit(1)

    if not getExecutable([externals['treetime']]) :
        url = 'https://github.com/neherlab/treetime/archive/refs/tags/v0.9.0.tar.gz'
        logger('Downloading treetime package from {0}'.format(url))
        subprocess.Popen('curl -Lo treetime.tar.gz {0}'.format(url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging treetime package')
        subprocess.Popen('tar -xzf treetime.tar.gz'.split()).communicate()
        subprocess.Popen('cp treetime-0.9.0/bin/treetime treetime-0.9.0/treetime_cmd'.split(), stderr=subprocess.PIPE).communicate()
        subprocess.Popen('ln -fs treetime-0.9.0/treetime_cmd ./treetime'.split(), stderr=subprocess.PIPE).communicate()
        subprocess.Popen('chmod 755 ./treetime'.split(), stderr=subprocess.PIPE).communicate()
        os.unlink('treetime.tar.gz')
        logger('Done\n')

    if not getExecutable([externals['hapog']]) :
        url = 'https://github.com/institut-de-genomique/HAPO-G/archive/refs/tags/1.2.tar.gz'
        logger('Downloading Hapo-G package from {0}'.format(url))
        subprocess.Popen('curl -Lo hapog.tar.gz {0}'.format(url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging Hapo-G package')
        subprocess.Popen('tar -xzf hapog.tar.gz'.split()).communicate()
        subprocess.Popen('bash build.sh', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='HAPO-G-1.2').communicate()
        subprocess.Popen('ln -fs HAPO-G-1.2/hapog.py ./hapog.py'.split(), stderr=subprocess.PIPE).communicate()
        subprocess.Popen('chmod 755 ./hapog.py'.split(), stderr=subprocess.PIPE).communicate()
        os.unlink('hapog.tar.gz')
        logger('Done\n')


    if not getExecutable([externals['raxml_ng']]) :
        url = 'https://github.com/amkozlov/raxml-ng/releases/download/1.0.1/raxml-ng_v1.0.1_linux_x86_64.zip'
        logger('Downloading raxml-ng package from {0}'.format(url))
        subprocess.Popen('curl -Lo raxml-ng_v1.0.1_linux_x86_64.zip {0}'.format(url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging raxml-ng package')
        subprocess.Popen('unzip raxml-ng_v1.0.1_linux_x86_64.zip -d raxml-ng_v1.0.1'.split(), stderr=subprocess.PIPE).communicate()
        subprocess.Popen('ln -fs raxml-ng_v1.0.1/raxml-ng ./raxml-ng'.split(), stderr=subprocess.PIPE).communicate()
        os.unlink('raxml-ng_v1.0.1_linux_x86_64.zip')
        logger('Done\n')

    if not getExecutable([externals['diamond']]) :
        url = 'https://github.com/bbuchfink/diamond/releases/download/v2.0.11/diamond-linux64.tar.gz'
        logger('Downloading diamond package from {0}'.format(url))
        subprocess.Popen('curl -Lo diamond-linux64.tar.gz {0}'.format(url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging diamond package'.format(url))
        subprocess.Popen('tar -xzf diamond-linux64.tar.gz'.split()).communicate()
        os.unlink('diamond-linux64.tar.gz')
        try :
            os.unlink('diamond_manual.pdf')
        except :
            pass
        logger('Done\n')

    if not getExecutable(externals['flye'].split()) :
        url = 'https://github.com/fenderglass/Flye/archive/refs/tags/2.8.3.tar.gz'
        logger('Downloading Flye from {0}'.format(url))
        subprocess.Popen('curl -Lo Flye.2.8.3.tar.gz {0}'.format(url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging Flye'.format(url))
        subprocess.Popen('tar -xzf Flye.2.8.3.tar.gz'.split()).communicate()
        os.unlink('Flye.2.8.3.tar.gz')
        subprocess.Popen('make', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd='Flye-2.8.3').communicate()
        subprocess.Popen('ln -fs Flye-2.8.3/bin/flye ./flye', shell=True).communicate()
        logger('Done\n')

    if not getExecutable([externals['spades']]) :
        url = 'https://github.com/ablab/spades/releases/download/v3.15.2/SPAdes-3.15.2-Linux.tar.gz'
        logger('Downloading SPAdes-3.15.2 package from {0}'.format(url))
        subprocess.Popen('curl -Lo SPAdes-3.15.2-Linux.tar.gz {0}'.format(url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging SPAdes-3.15.2-Linux package'.format(url))
        subprocess.Popen('tar -xzf SPAdes-3.15.2-Linux.tar.gz'.split()).communicate()
        os.unlink('SPAdes-3.15.2-Linux.tar.gz')
        subprocess.Popen('ln -fs SPAdes-3.15.2-Linux/bin/spades.py ./spades.py'.split()).communicate()
        logger('Done\n')

    if not getExecutable([externals['bbduk']]) or not getExecutable([externals['bbmerge']]) or not getExecutable([externals['repair']]):
        url = 'https://netcologne.dl.sourceforge.net/project/bbmap/BBMap_38.90.tar.gz'
        logger('Downloading BBmap package from {0}'.format(url))
        subprocess.Popen('curl -Lo BBMap_38.90.tar.gz {0}'.format(url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging BBmap package'.format(url))
        subprocess.Popen('tar -xzf BBMap_38.90.tar.gz'.split()).communicate()
        os.unlink('BBMap_38.90.tar.gz')
        subprocess.Popen('ln -fs bbmap/bbduk.sh ./bbduk.sh'.split()).communicate()
        subprocess.Popen('ln -fs bbmap/bbmerge.sh ./bbmerge.sh'.split()).communicate()
        subprocess.Popen('ln -fs bbmap/repair.sh ./repair.sh'.split()).communicate()
        logger('Done\n')

    if not getExecutable([externals['megahit']]) :
        megahit_url = 'https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz'
        # megahit_url = 'https://github.com/voutcn/megahit/releases/download/v1.1.4/megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz'
        logger('Downloading megahit package from {0}'.format(megahit_url))
        subprocess.Popen('curl -Lo MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz {0}'.format(megahit_url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging megahit package'.format(megahit_url))
        subprocess.Popen('tar -xzf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz'.split()).communicate()
        os.unlink('MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz')
        subprocess.Popen('ln -fs MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit ./megahit'.split()).communicate()
        logger('Done\n')

    if not getExecutable([externals['bowtie2']]) :
        bowtie2_url = 'https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip'
        logger('Downloading bowtie2 package from {0}'.format(bowtie2_url))
        subprocess.Popen('curl -Lo bowtie2-2.3.4.3-linux-x86_64.zip {0}'.format(bowtie2_url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging bowtie2 package')
        subprocess.Popen('unzip -o bowtie2-2.3.4.3-linux-x86_64.zip'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        os.unlink('bowtie2-2.3.4.3-linux-x86_64.zip')
        subprocess.Popen('ln -fs bowtie2-2.3.4.3-linux-x86_64/bowtie2 ./bowtie2'.split()).communicate()
        subprocess.Popen('ln -fs bowtie2-2.3.4.3-linux-x86_64/bowtie2-build ./bowtie2-build'.split()).communicate()
        logger('Done\n')

    if not getExecutable([externals['mmseqs']]) :
        if int(subprocess.Popen('cat /proc/cpuinfo | grep avx2|wc', shell=True, stdout=subprocess.PIPE).communicate()[0].strip().split()[0]) :
            mmseqs_url = 'https://github.com/soedinglab/MMseqs2/releases/download/13-45111/mmseqs-linux-avx2.tar.gz'
        elif int(subprocess.Popen('cat /proc/cpuinfo | grep SSE4|wc', shell=True, stdout=subprocess.PIPE).communicate()[0].strip().split()[0]) :
            mmseqs_url = 'https://github.com/soedinglab/MMseqs2/releases/download/13-45111/mmseqs-linux-sse41.tar.gz'
        else :
            mmseqs_url = 'https://github.com/soedinglab/MMseqs2/releases/download/13-45111/mmseqs-linux-arm64.tar.gz'
        logger('Downloading mmseqs package from {0}'.format(mmseqs_url))
        subprocess.Popen('curl -Lo MMseqs2.tar.gz {0}'.format(mmseqs_url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging mmseqs package')
        subprocess.Popen('tar -xzf MMseqs2.tar.gz'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        os.unlink('MMseqs2.tar.gz')
        os.rename('mmseqs', 'mmseqs2')
        subprocess.Popen('ln -fs mmseqs2/bin/mmseqs ./mmseqs'.split()).communicate()
        logger('Done\n')

    if not getExecutable([externals['lastal']]) or not getExecutable([externals['lastdb']]) :
                
        last_url = 'https://gitlab.com/mcfrith/last/-/archive/05ffb83b103f5b6b117807b90f82353e8c27490c/last-05ffb83b103f5b6b117807b90f82353e8c27490c.tar.gz'
        logger('Downloading LAST package from {0}'.format(last_url))
        subprocess.Popen('curl -Lo last-1021.tar.gz {0}'.format(last_url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging LAST package'.format(last_url))
        subprocess.Popen('tar -xf last-1021.tar.gz'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        os.unlink('last-1021.tar.gz')
        os.rename('last-05ffb83b103f5b6b117807b90f82353e8c27490c', 'last-1021')
        logger('Installing LAST package'.format(last_url))
        subprocess.Popen('make'.split(), cwd='last-1021/src', stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        subprocess.Popen('ln -fs last-1021/src/lastal ./lastal'.split()).communicate()
        subprocess.Popen('ln -fs last-1021/src/lastdb ./lastdb'.split()).communicate()
        logger('Done\n')

    if not getExecutable(externals['kraken2'].split()) :
        kraken2_url = 'https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.tar.gz'
        logger('Downloading kraken2 package from {0}'.format(kraken2_url))
        subprocess.Popen('curl -Lo kraken2.1.2.tar.gz {0}'.format(kraken2_url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging kraken2 package')
        subprocess.Popen('tar -xzf kraken2.1.2.tar.gz'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        os.unlink('kraken2.1.2.tar.gz')
        subprocess.Popen('cd kraken2-2.1.2 && bash install_kraken2.sh ./', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()
        subprocess.Popen('ln -fs kraken2-2.1.2/kraken2 ./kraken2'.split()).communicate()
        logger('Done\n')

    if not getExecutable(externals['samtools'].split()) :
        samtools_url = 'https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2'
        logger('Downloading samtools from {0}'.format(samtools_url))
        subprocess.Popen('curl -Lo samtools-1.15.tar.bz2 {0}'.format(samtools_url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging samtools package')
        subprocess.Popen('tar -xjf samtools-1.15.tar.bz2'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        os.unlink('samtools-1.15.tar.bz2')
        subprocess.Popen('cd samtools-1.15 && ./configure --disable-bz2 --disable-lzma --without-curses && make', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()
        subprocess.Popen('ln -fs samtools-1.15/samtools ./samtools'.split()).communicate()
        logger('Done\n')

    if not getExecutable([externals['nextpolish']]):
        url = 'https://github.com/Nextomics/NextPolish/releases/download/v1.4.1/NextPolish.tgz'
        logger('Downloading nextPolish package from {0}'.format(url))
        subprocess.Popen('curl -Lo nextpolish-1.4.1.tgz {0}'.format(url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging nextPolish package')
        subprocess.Popen('tar -vxzf nextpolish-1.4.1.tgz'.split()).communicate()
        os.chdir('NextPolish')
        subprocess.Popen('make', stderr=subprocess.PIPE).communicate()
        os.chdir(moveTo)
        os.unlink('nextpolish-1.4.1.tgz')
        logger('Done\n')
        os.chdir(curdir)

    if not getExecutable([externals['blastn']]) or not getExecutable([externals['makeblastdb']]) :
        blast_url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-x64-linux.tar.gz'
        logger('Downloading ncbi-blast package from {0}'.format(blast_url))
        subprocess.Popen('curl -Lo ncbi-blast-2.8.1+-x64-linux.tar.gz {0}'.format(blast_url).split(), stderr=subprocess.PIPE).communicate()
        logger('Unpackaging ncbi-blast package'.format(blast_url))
        subprocess.Popen('tar -xzf ncbi-blast-2.8.1+-x64-linux.tar.gz'.split()).communicate()
        os.unlink('ncbi-blast-2.8.1+-x64-linux.tar.gz')
        subprocess.Popen('ln -fs ncbi-blast-2.8.1+/bin/blastn ./blastn'.split()).communicate()
        subprocess.Popen('ln -fs ncbi-blast-2.8.1+/bin/makeblastdb ./makeblastdb'.split()).communicate()
        logger('Done\n')

    if not getExecutable([externals['usearch']]) :
        logger('The 32-bit version of USEARCH is licensed at no charge for individual use. \nPlease download it at    https://www.drive5.com/usearch/download.html\nAnd copy it into the externals/usearch')
    logger('')
    os.chdir(curdir)
# -------------------------------------------------------------- #
ETOKI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def configure(args) :
    configs = load_configure()
    args = add_args(args)

    for key, value in args.__dict__.items() :
        if value is not None :
            configs[configs.T[0] == key, 1] = value
    externals = prepare_externals(conf=configs)
    if args.install :
        install_externals()
    if args.download_krakenDB :
        download_krakenDB()
    for fname, flink in sorted(externals.items()) :
        flinks = flink.split()
        if fname not in {'kraken_database', 'enbler_filter', 'pigz'} :
            if not getExecutable(flinks) :
                logger('ERROR - {0} ("{1}") is not present. '.format(fname, flinks[-1]))
                #sys.exit(0)
            else :
                logger('{0} ("{1}") is present. '.format(fname, flinks[-1]))
    if not os.path.exists(externals['kraken_database']) :
        logger('''WARNING - kraken_database is not present. 
You can still use EToKi except the parameter "--kraken" in EToKi assemble will not work.
Alternatively you can download minikraken2 database using --download_krakenDB or pass an pre-installed database into EToKi using --link_krakenDB.''')
    
    write_configure(configs)
    logger('Configuration complete.')

def prepare_externals(conf=None) :
    if conf is None :
        conf = load_configure()
    externals = {k.strip():v.split('#')[0].strip().format(ETOKI=ETOKI) for k,v in conf.tolist()}
    externals['treetime'] = sys.executable + ' ' + externals.get('treetime', '')
    externals['enbler_filter'] = sys.executable + ' {ETOKI}/modules/_EnFlt.py'.format(ETOKI=ETOKI)
    externals['pigz'] = getExecutable(['pigz'])[0] if getExecutable(['pigz']) else getExecutable(['gzip'])[0]
    return externals

def add_args(a) :
    parser = argparse.ArgumentParser(description='''Install or modify the 3rd party programs.''')
    parser.add_argument('--install', help='install 3rd party programs', default=False, action='store_true')
    parser.add_argument('--usearch', dest='usearch', help='usearch is required for ortho and MLSType. A 32-bit version of usearch can be downloaded from https://www.drive5.com/usearch/', default=None)
    parser.add_argument('--download_krakenDB', help='When specified, miniKraken2 (8GB) will be downloaded into the EToKi folder. You can also use --link_krakenDB to use a pre-installed kraken2 database.', default=False, action='store_true')
    parser.add_argument('--link_krakenDB', dest='kraken_database', help='Kraken is optional in the assemble module. You can specify your own database here', default=None)
    parser.add_argument('--path', '-p', help='Specify path to the 3rd party programs manually. Format: <program>=<path>. This parameter can be specified multiple times', default=[], action='append')
    args = parser.parse_args(a)
    for ps in args.path :
        k, v = ps.split('=')
        args.__dict__[k] = v
    return args


def load_configure() :
    EnConf_file = os.path.realpath(__file__).rsplit('.', 1)[0] + '.ini'
    try :
        return pd.read_csv(EnConf_file, sep='=', header=None).values
    except :
        return np.array([0, 2], dtype=str)
    

def write_configure(configs) :
    EnConf_file = os.path.realpath(__file__).rsplit('.', 1)[0] + '.ini'
    pd.DataFrame(configs).to_csv(EnConf_file, sep='=', index=False, header=False)

externals = prepare_externals()
if __name__ == '__main__' :
    configure(sys.argv[1:])

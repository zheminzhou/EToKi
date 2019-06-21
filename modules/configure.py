import os, sys, subprocess, numpy as np, pandas as pd, argparse, glob, gzip, io, re, requests
import multiprocessing
import multiprocessing.pool

from datetime import datetime
if sys.version_info[0] < 3:
    from cStringIO import StringIO
    xrange = xrange
    asc2int = np.uint8
else :
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
            if sys.version.startswith('3') :
                self.fout = gzip.open(fname, 'wb')
                self.fstream = io.TextIOWrapper(self.fout, encoding='utf-8')
            else :
                self.fstream = gzip.open(fname, 'wb')
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
    sequence = {}
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
    sequence, qual = {}, {}
    with uopen(fastq) as fin :
        line = fin.readline()
        if not line.startswith('@') :
            sequence = readFasta(fastq)
            return sequence, { n: re.sub(r'[^!]', 'I', re.sub(r'[^ACGTacgt]', '!', s)) for n, s in sequence.items() }
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
def rc(seq) :
    return ''.join([complement.get(s, 'N') for s in reversed(seq.upper())])


baseConv = np.empty(255, dtype=int)
baseConv.fill(-100)
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


def checkExecutable(commands) :
    try :
        return subprocess.Popen(commands+['-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait() <= 1 or subprocess.Popen(commands, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait() <= 1
    except  :
        return False

def download_krakenDB() :
    curdir = os.path.abspath(os.curdir)
    moveTo = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'externals')
    os.chdir(moveTo)
    if not os.path.exists('minikraken2') :
        os.makedirs('minikraken2')
    os.chdir(os.path.join(moveTo, 'minikraken2'))
    minikraken_url = 'ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz'
    logger('Downloading minikraken2 from {0}. This might take a long time.'.format(minikraken_url))    
    subprocess.Popen('curl -Lo minikraken2_v2_8GB.tgz {0}'.format(minikraken_url).split(), stderr=subprocess.PIPE).wait()
    logger('Unpackaging minikraken2.')
    subprocess.Popen('tar -xzf minikraken2_v2_8GB.tgz'.split()).wait()
    subprocess.Popen('mv minikraken2_v2_8GB_*/* ./', shell=True).wait()
    os.unlink('minikraken2_v2_8GB.tgz')
    
    os.chdir(curdir)


def install_externals() :
    curdir = os.path.abspath(os.curdir)
    moveTo = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'externals')
    os.chdir(moveTo)
    
    if not checkExecutable([externals['blastn']]) or not checkExecutable([externals['makeblastdb']]) :
        blast_url = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-x64-linux.tar.gz'
        logger('Downloading ncbi-blast package from {0}'.format(blast_url))
        subprocess.Popen('curl -Lo ncbi-blast-2.8.1+-x64-linux.tar.gz {0}'.format(blast_url).split(), stderr=subprocess.PIPE).wait()
        logger('Unpackaging ncbi-blast package'.format(blast_url))
        subprocess.Popen('tar -xzf ncbi-blast-2.8.1+-x64-linux.tar.gz'.split()).wait()
        os.unlink('ncbi-blast-2.8.1+-x64-linux.tar.gz')
        subprocess.Popen('ln -fs ncbi-blast-2.8.1+/bin/blastn ./blastn'.split()).wait()
        subprocess.Popen('ln -fs ncbi-blast-2.8.1+/bin/makeblastdb ./makeblastdb'.split()).wait()
        logger('Done\n')
    
    if not checkExecutable([externals['megahit']]) :
        megahit_url = 'https://github.com/voutcn/megahit/releases/download/v1.1.4/megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz'
        logger('Downloading megahit package from {0}'.format(megahit_url))
        subprocess.Popen('curl -Lo megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz {0}'.format(megahit_url).split(), stderr=subprocess.PIPE).wait()
        logger('Unpackaging megahit package'.format(megahit_url))
        subprocess.Popen('tar -xzf megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz'.split()).wait()
        os.unlink('megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz')
        subprocess.Popen('ln -fs megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin/megahit ./megahit'.split()).wait()
        logger('Done\n')

    if not checkExecutable([externals['bowtie2']]) :
        bowtie2_url = 'https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip'
        logger('Downloading bowtie2 package from {0}'.format(bowtie2_url))
        subprocess.Popen('curl -Lo bowtie2-2.3.4.3-linux-x86_64.zip {0}'.format(bowtie2_url).split(), stderr=subprocess.PIPE).wait()
        logger('Unpackaging bowtie2 package')
        subprocess.Popen('unzip -o bowtie2-2.3.4.3-linux-x86_64.zip'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
        os.unlink('bowtie2-2.3.4.3-linux-x86_64.zip')
        subprocess.Popen('ln -fs bowtie2-2.3.4.3-linux-x86_64/bowtie2 ./bowtie2'.split()).wait()
        subprocess.Popen('ln -fs bowtie2-2.3.4.3-linux-x86_64/bowtie2-build ./bowtie2-build'.split()).wait()
        logger('Done\n')

    if not checkExecutable([externals['mmseqs']]) :
        mmseqs_url = 'https://github.com/soedinglab/MMseqs2/releases/download/7-4e23d/MMseqs2-Linux-SSE4_1.tar.gz'
        logger('Downloading mmseqs package from {0}'.format(mmseqs_url))
        subprocess.Popen('curl -Lo MMseqs2-Linux-SSE4_1.tar.gz {0}'.format(mmseqs_url).split(), stderr=subprocess.PIPE).wait()
        logger('Unpackaging mmseqs package')
        subprocess.Popen('tar -xzf MMseqs2-Linux-SSE4_1.tar.gz'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
        os.unlink('MMseqs2-Linux-SSE4_1.tar.gz')
        subprocess.Popen('ln -fs mmseqs2/bin/mmseqs ./mmseqs'.split()).wait()
        logger('Done\n')

    if not checkExecutable(externals['gatk'].split()) :
        gatk_url = 'https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip'
        logger('Downloading gatk package from {0}'.format(gatk_url))
        subprocess.Popen('curl -Lo gatk-4.1.0.0.zip {0}'.format(gatk_url).split(), stderr=subprocess.PIPE).wait()
        logger('Unpackaging gatk package')
        subprocess.Popen('unzip -o gatk-4.1.0.0.zip'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
        os.unlink('gatk-4.1.0.0.zip')
        subprocess.Popen('ln -fs gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar ./gatk-package-4.1.0.0-local.jar'.split()).wait()
        logger('Done\n')

    if not checkExecutable(externals['kraken2'].split()) :
        kraken2_url = 'https://github.com/DerrickWood/kraken2/archive/v2.0.7-beta.tar.gz'
        logger('Downloading kraken2 package from {0}'.format(kraken2_url))
        subprocess.Popen('curl -Lo v2.0.7-beta.tar.gz {0}'.format(kraken2_url).split(), stderr=subprocess.PIPE).wait()
        logger('Unpackaging kraken2 package')
        subprocess.Popen('tar -xzf v2.0.7-beta.tar.gz'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
        os.unlink('v2.0.7-beta.tar.gz')
        subprocess.Popen('cd kraken2-2.0.7-beta && bash install_kraken2.sh ./', shell=True).wait()
        subprocess.Popen('ln -fs kraken2-2.0.7-beta/kraken2 ./kraken2'.split()).wait()
        logger('Done\n')

    if not checkExecutable(externals['diamond'].split()) :
        diamond_url = 'https://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz'
        logger('Downloading diamond package from {0}'.format(diamond_url))
        subprocess.Popen('curl -Lo diamond-linux64.tar.gz {0}'.format(diamond_url).split(), stderr=subprocess.PIPE).wait()
        logger('Unpackaging diamond package')
        subprocess.Popen('tar -xzf diamond-linux64.tar.gz'.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()
        os.unlink('diamond-linux64.tar.gz')
        os.unlink('diamond_manual.pdf')
        logger('Done\n')

    if not checkExecutable([externals['usearch']]) :
        logger('The 32-bit version of USEARCH is licensed at no charge for individual use. \nPlease download it at    https://www.drive5.com/usearch/download.html')
        logger('')
    os.chdir(curdir)
# -------------------------------------------------------------- #
ETOKI = os.path.dirname(os.path.dirname(__file__))
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
        if fname not in {'kraken_database', 'enbler_filter'} :
            if not checkExecutable(flinks) :
                logger('ERROR - {0} ("{1}") is not present. '.format(fname, flinks[-1]))
                sys.exit(0)
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
    externals['gatk']  = 'java -Xmx31g -jar ' + externals.get('gatk', '')
    externals['pilon'] = 'java -Xmx63g -jar ' + externals.get('pilon', '')
    externals['enbler_filter'] = sys.executable + ' {ETOKI}/modules/_EnFlt.py'.format(ETOKI=ETOKI)
    externals['pigz'] = 'pigz' if checkExecutable(['pigz']) else 'gzip'
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

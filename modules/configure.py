import os, sys, subprocess, numpy as np, pandas as pd, argparse, glob, gzip, io, re
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
            self.fstream = subprocess.Popen(['pigz', '-cd', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout if fname.lower().endswith('gz') else open(fname)
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


def readFasta(fasta) :
    sequence = {}
    with uopen(fasta) as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                sequence[name] = []
            elif len(line) > 0 and not line.startswith('#') :
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
def transeq(seq, frame=7, transl_table=None) :
    frames = {'F': [1,2,3],
              'R': [4,5,6],
              '7': [1,2,3,4,5,6]}.get( str(frame).upper() , None)
    if frames is None :
        frames = [int(f) for f in str(frame).split(',')]
    
    if transl_table == 'starts' :
        gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVMVXYXYSSSSXCWCLFMF-'))
    else :
        gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF-'))
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


# -------------------------------------------------------------- #
ETOKI = os.path.dirname(os.path.dirname(__file__))
def configure(args) :
    configs = load_configure()
    args = add_args(args)

    for key, value in args.__dict__.items() :
        if value is not None :
            configs[configs.T[0] == key, 1] = value
    externals = prepare_externals(conf=configs)
    for fname, flink in sorted(externals.items()) :
        flinks = flink.split()
        #
        try:        
            if fname in {'kraken_database', 'adapters'} :
                if not os.path.exists(flinks[-1]) :
                    raise FileNotFoundError
            else :
                subprocess.Popen(flinks, stderr=subprocess.PIPE, stdout=subprocess.PIPE).wait()
            logger('{0} ("{1}") is present. '.format(fname, flinks[-1]))
        except FileNotFoundError as e:
            logger('ERROR - {0} ("{1}") is not present. '.format(fname, flinks[-1]))
            sys.exit(0)
    write_configure(configs)
    logger('Configuration complete.')

def prepare_externals(conf=None) :
    if conf is None :
        conf = load_configure()
    externals = {k.strip():v.split('#')[0].strip().format(ETOKI=ETOKI) for k,v in conf.tolist()}
    externals['gatk']  = 'java -Xmx31g -jar ' + externals.get('gatk', '')
    externals['pilon'] = 'java -Xmx63g -jar ' + externals.get('pilon', '')
    externals['enbler_filter'] = sys.executable + ' ' + externals.get('enbler_filter', '')
    return externals

def add_args(a) :
    parser = argparse.ArgumentParser(description='''Specify links to kraken database and usearch program.''')
    parser.add_argument('--usearch', dest='usearch', help='usearch is required for ortho and MLSType. Download the 32-bit version from https://www.drive5.com/usearch/.', default=None)
    parser.add_argument('--kraken_database', dest='kraken_database', help='Kraken is optional in the assemble module. You can specify your own database or use MiniKraken2: https://ccb.jhu.edu/software/kraken2/dl/minikraken2_v2_8GB.tgz', default=None)
    return parser.parse_args(a)


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

import os, sys, subprocess, numpy as np, argparse, glob, gzip, io, re
from datetime import datetime
if sys.version_info[0] < 3:
    xrange = xrange
    asc2int = np.uint8
else :
    xrange = range
    asc2int = np.uint32
    
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
            return sequence, None
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


conv = np.empty(255, dtype=int)
conv.fill(-100)
conv[(np.array(['-', 'A', 'C', 'G', 'T']).view(asc2int),)] = (-100000, 0, 1, 2, 3)
def transeq(seq, frame=7) :
    frames = {'F': [1,2,3],
              'R': [4,5,6],
              '7': [1,2,3,4,5,6]}.get( str(frame).upper() , [frame])
    
    gtable = np.array(list('KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVXYXYSSSSXCWCLFLF-'))
    seqs = seq.items() if isinstance(seq, dict) else seq
    nFrame = (max(frames) > 3)
    trans_seq = []
    for n,s in seqs :
        s = conv[np.array(list(s.upper())).view(asc2int)]
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

search_file = lambda x: [os.path.abspath(x)] if os.path.exists(x) else [os.path.abspath(fname) for path in os.environ["PATH"].split(os.pathsep) for fname in sorted(glob.glob(os.path.join(path, x)))]
def configure(args) :
    global externals
    externals = load_configure()
    defaults = dict(
        adapters='adapters.fa', 
        kraken_database='minikraken*', 
        bbduk='bbduk2.sh', 
        bbmerge='bbmerge.sh',
        spades='spades.py', 
        megahit='megahit', 
        minimap2='minimap2', 
        bwa='bwa', 
        bowtie2='bowtie2', 
        bowtie2build='bowtie2-build', 
        samtools='samtools', 
        gatk='gatk-package-4*', 
        pilon='pilon*', 
        kraken_program='kraken', 
        kraken_report='kraken-report', 
        mmseqs='mmseqs', 
        fasttree='?fast?ree*',
        rapidnj='rapidnj',
        ublast='usearch*',
        blastn='blastn',
        formatdb='makeblastdb',
        raxml='raxml*', 
        enbler_filter=os.path.join(os.path.dirname(os.path.realpath(__file__)), '_EnFlt.py')
    )
    args = add_args(args)

    for param in defaults :
        value = defaults[param]
        if args.__dict__.get(param, None) :
            fn = search_file(args.__dict__[param])
            assert len(fn), 'The specified "{0}" is not found.'.format(param)
            logger('Found {0}'.format(fn[0]))
            externals[param] = fn[0]
        elif not os.path.exists(externals.get(param, '')) :
            fn = search_file(value)
            if not len(fn) :
                logger('{0} is not found. Be aware that some functions may not able to run.'.format(param))
            else :
                logger('Found {0}'.format(fn[0]))
                externals[param] = fn[0]

    write_configure()
    logger('Configuration complete.')

def prepare_externals() :
    externals['gatk']  = 'java -Xmx63g -jar ' + externals.get('gatk', '')
    externals['pilon'] = 'java -Xmx63g -jar ' + externals.get('pilon', '')
    externals['enbler_filter'] = sys.executable + ' ' + externals.get('enbler_filter', '')

def add_args(a) :
    parser = argparse.ArgumentParser(description='''Configure external dependencies for EToKi (Enterobase Tool Kit).
Will search executable files from system paths by default.
The path to two databases for Kraken and Illumina adapters are required for the assembler. ''', formatter_class=argparse.RawTextHelpFormatter)
    # EnBler database
    parser.add_argument('--adapters', help='adapter database to be used in bbduk.\nIt can be found as resources/adapters.fa in the BBmap package')
    parser.add_argument('--krakenDB', dest='kraken_database', help='database to be used in kraken based species prediction.\nIt can be downloaded from \n"https://ccb.jhu.edu/software/kraken/"')
    # EnBler executables
    parser.add_argument('--bbduk')
    parser.add_argument('--bbmerge')
    parser.add_argument('--spades')
    parser.add_argument('--bwa')
    parser.add_argument('--minimap2')

    parser.add_argument('--megahit')
    parser.add_argument('--bowtie2')
    parser.add_argument('--bowtie2-build', dest='bowtie2build')

    parser.add_argument('--samtools')
    parser.add_argument('--gatk4', dest='gatk')
    parser.add_argument('--pilon')
    parser.add_argument('--kraken', dest='kraken_program')
    parser.add_argument('--kraken-report', dest='kraken_report')
    # EnOrth executables
    parser.add_argument('--mmseqs')
    parser.add_argument('--fasttree')
    # EnSign executables
    parser.add_argument('--usearch', dest='ublast')
    parser.add_argument('--blastn')
    parser.add_argument('--makeblastdb', dest='formatdb')
    # EnPhyl executables
    parser.add_argument('--raxml')

    return parser.parse_args(a)


def load_configure() :
    EnConf_file = os.path.realpath(__file__).rsplit('.', 1)[0] + '.ini'
    try :
        mat = np.genfromtxt(EnConf_file, dtype=str, delimiter='=')
        return dict(mat.tolist())
    except :
        return {}
    

def write_configure() :
    EnConf_file = os.path.realpath(__file__).rsplit('.', 1)[0] + '.ini'
    with open(EnConf_file, 'w') as fout :
        for k,v in externals.items() :
            fout.write('{0}={1}\n'.format(k, v))

externals = load_configure()
if __name__ == '__main__' :
    configure(sys.argv[1:])
else :
    prepare_externals()

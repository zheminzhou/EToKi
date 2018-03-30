import os, sys, subprocess, numpy as np, argparse, glob

# packages:
#   bbmap
#   SPAdes-3.9.0
#   megahit
#   bowtie2
#   samtools
#   gatk
#   pilon
#   kraken

#   usearch
#   blast


def readFasta(filename, qual=None) :
    seq = {}
    try:
        fin = subprocess.Popen(['zcat', filename], stdout=subprocess.PIPE) if filename[-3:].lower() == '.gz' else subprocess.Popen(['cat', filename], stdout=subprocess.PIPE).stdout
    except :
        fin = open(filename)

    for line in fin:
        if line[0] == '>' :
            name = line[1:].strip().split()[0]
            seq[name] = []
        else :
            seq[name].append( line.strip() )
    fin.close()

    if qual == None :
        for n in seq:
            seq[n] = ''.join( seq[n] )
    else :
        for n in seq:
            seq[n] = [''.join(seq[n]), []]
            seq[n][1] = [chr(33+qual)] * len(seq[n][0])
            if qual > 0 :
                for id in xrange(len(seq[n][0])) :
                    if seq[n][0] in ('n', 'N') :
                        seq[n][1] = '!'
            seq[n][1] = ''.join(seq[n][1])
    return seq


complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
def rc(seq) :
    return ''.join([complement.get(s, 'N') for s in reversed(seq.upper())])

def transeq(seq, frame=1, transl_table=11) :
    gtable = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
              "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",    "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
              "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
              "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
              "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
              "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
              "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
              "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

    if transl_table == 'starts' :
        gtable.update({'GTG':'M', 'TTG':'M'})
    frames = {'F': [1,2,3],
              'R': [4,5,6],
              '7': [1,2,3,4,5,6,7]}.get( str(frame).upper(), [int(frame)] )
    trans_seq = {}
    for n,s in seq.iteritems() :
        for frame in frames :
            trans_name = '{0}_{1}'.format(n, frame)
            if frame <= 3 :
                trans_seq[trans_name] = ''.join([gtable.get(c, 'X') for c in map(''.join, zip(*[iter(s[(frame-1):])]*3))])
            else :
                trans_seq[trans_name] = ''.join([gtable.get(c, 'X') for c in map(''.join, zip(*[iter(rc(s)[(frame-4):])]*3))])
    return trans_seq


def logger(log) :
    from datetime import datetime
    sys.stderr.write('{0}\t{1}\n'.format(str(datetime.now()), log))

search_file = lambda x: [x] if os.path.exists(x) else [fname for path in os.environ["PATH"].split(os.pathsep) for fname in sorted(glob.glob(os.path.join(path, x)))]
def EnConf(args) :
    global externals
    externals = load_configure()
    defaults = dict(
        adapters='adapters.fa',
        kraken_database='minikraken*',
        bbduk='bbduk2.sh',
        spades='spades.py',
        megahit='megahit',
        bwa='bwa',
        bowtie2='bowtie2',
        bowtie2build='bowtie2-build',
        samtools='samtools',
        gatk='gatk-package-4*',
        pilon='pilon*',
        kraken_program='kraken',
        kraken_report='kraken-report',
        vsearch='vsearch',
        mcl='mcl',
        fasttree='?fast?ree*',
        ublast='usearch*',
        blastn='blastn',
        formatdb='makeblastdb',
        raxml='raxml*', 
        enbler_filter=os.path.join(os.path.dirname(os.path.realpath(__file__)), '_EnFlt.py')
    )
    args = add_args(args)

    for param, value in defaults.iteritems() :
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
    externals['gatk']  = 'java -Xmx30g -jar ' + externals['gatk']
    externals['pilon'] = 'java -Xmx30g -jar ' + externals['pilon']
    externals['enbler_filter'] = 'python ' + externals['enbler_filter']

def add_args(a) :
    parser = argparse.ArgumentParser(description='''Configure external dependencies for EToKi (Enterobase Tool Kit).
Will search executable files from system paths by default.
The path to two databases for Kraken and Illumina adapters are required for the assembler. ''', formatter_class=argparse.RawTextHelpFormatter)
    # EnBler database
    parser.add_argument('--adapters', help='adapter database to be used in bbduk.\nIt can be found as resources/adapters.fa in the BBmap package')
    parser.add_argument('--krakenDB', dest='kraken_database', help='database to be used in kraken based species prediction.\nIt can be downloaded from \n"https://ccb.jhu.edu/software/kraken/"')
    # EnBler executables
    parser.add_argument('--bbduk')
    parser.add_argument('--spades')
    parser.add_argument('--bwa')

    parser.add_argument('--megahit')
    parser.add_argument('--bowtie2')
    parser.add_argument('--bowtie2-build', dest='bowtie2build')

    parser.add_argument('--samtools')
    parser.add_argument('--gatk4', dest='gatk')
    parser.add_argument('--pilon')
    parser.add_argument('--kraken', dest='kraken_program')
    parser.add_argument('--kraken-report', dest='kraken_report')
    # EnOrth executables
    parser.add_argument('--vsearch')
    parser.add_argument('--mcl')
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
    except :
        return {}
    return dict(mat.tolist())

def write_configure() :
    EnConf_file = os.path.realpath(__file__).rsplit('.', 1)[0] + '.ini'
    with open(EnConf_file, 'wb') as fout :
        for k,v in externals.iteritems() :
            fout.write('{0}={1}\n'.format(k, v))

externals = load_configure()
if __name__ == '__main__' :
    EnConf(sys.argv[1:])
else :
    prepare_externals()
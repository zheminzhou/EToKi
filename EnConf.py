import os, sys, subprocess

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

externals = dict(
    # EnBler
    bbduk = '{HOME}/CRobot/source/bbmap/bbduk2.sh',
    adapters = '{HOME}/CRobot/source/bbmap/resources/adapters.fa',
    spades = '{HOME}/CRobot/source/SPAdes-3.9.0/spades.py',
    megahit = '{HOME}/CRobot/source/megahit/megahit',
    bowtie2 = '{HOME}/bin/bowtie2',
    bowtie2build = '{HOME}/bin/bowtie2-build',
    enbler_filter= 'python {HOME}/CRobot/pipelines/EnSuit/_EnFlt.py',
    samtools = '{HOME}/CRobot/source/samtools-1.2/samtools',
    gatk = 'java -Xmx30g -jar {HOME}/CRobot/source/gatk-4.beta.6/gatk-package-4.beta.6-local.jar',
    pilon = 'java -Xmx30g -jar {HOME}/CRobot/source/pilon-1.22.jar',
    kraken_program = '{HOME}/CRobot/source/kraken/kraken',
    kraken_report = '{HOME}/CRobot/source/kraken/kraken-report',
    kraken_database = '{HOME}/minikraken_20141208',

    # EnSign
    ublast='{HOME}/NServ/utils/usearch8.0.1623_i86linux32',
    blast='{HOME}/NServ/utils/ncbi-blast-2.2.31+/bin/blastn',
    formatdb='{HOME}/NServ/utils/ncbi-blast-2.2.31+/bin/makeblastdb',
    #
    fasttree = '{HOME}/biosoft/FastTreeMP',
    mcl = '/usr/local/bin/mcl',
)
externals = {k:v.format(HOME=os.path.expanduser('~')) for k, v in externals.iteritems()}

if __name__ == '__main__' :
    for k, p in externals.iteritems() :
        fp = p.split()[-1]
        assert os.path.exists(fp), 'ERROR - {0} : {1} is not present'.format(k, p)
        print 'Found - {1}'.format(k, p)

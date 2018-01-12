import os, sys

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

def logger(log) :
    from datetime import datetime
    sys.stderr.write('{0}\t{1}\n'.format(str(datetime.now()), log))

externals = dict(
    # EnBler
    bbduk = '{HOME}/CRobot/source/bbmap/bbduk2.sh',
    adapters = '{HOME}/CRobot/source/bbmap/resources/adapters.fa',
    spades = '{HOME}/CRobot/source/SPAdes-3.9.0/spades.py',
    megahit = '{HOME}/software/megahit/megahit',
    bowtie2 = '{HOME}/bin/bowtie2',
    bowtie2build = '{HOME}/bin/bowtie2-build',
    enbler_filter= 'python {HOME}/CRobot/pipelines/EnTools/_EnFlt.py',    
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
)
externals = {k:v.format(HOME=os.path.expanduser('~')) for k, v in externals.iteritems()}

if __name__ == '__main__' :
    for k, p in externals.iteritems() :
        fp = p.split()[-1]
        assert os.path.exists(fp), 'ERROR - {0} : {1} is not present'.format(k, p)
        print 'Found - {1}'.format(k, p)
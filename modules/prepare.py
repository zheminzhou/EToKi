import os, io, sys, re, shutil, numpy as np, signal, psutil, argparse
from glob import glob
from subprocess import Popen, PIPE, STDOUT
from time import sleep
from configure import externals, logger, readFasta
from threading import Timer

def kill_child_proc (p) :
    for child in psutil.Process(p.pid).children() :
        child.terminate()

# preprocessing
class preprocess(object) :
    def launch (self, reads) :
        reads = self.init_cleanup(reads)
        reads = self.reduce_depth(reads)
        return reads

    def init_cleanup(self, reads) :
        prefix = parameters['prefix']
        new_reads = []
        for lib_id, library in enumerate(reads) :
            library_file = ['{0}.0.{1}.1.fastq.gz'.format(prefix, lib_id)]
            Popen('cat {0} > {1}'.format(' '.join([run[0] for run in library]), library_file[0]), shell=True).wait()
            if len(library[0]) > 1 :
                library_file.append('{0}.0.{1}.2.fastq.gz'.format(prefix, lib_id))
                Popen('cat {0} > {1}'.format(' '.join([run[1] for run in library]), library_file[1]), shell=True).wait()
            if len(library_file) == 1 :
                reads = 'in=' + library_file[0]
                library_file2 = ['{0}.1.{1}.1.fastq.gz'.format(prefix, lib_id)]
                outputs = 'out=' + library_file2[0]
            else :
                reads = 'in=' + library_file[0] + ' in2=' + library_file[1]
                library_file2 = ['{0}.1.{1}.1.fastq.gz'.format(prefix, lib_id), '{0}.1.{1}.2.fastq.gz'.format(prefix, lib_id), '{0}.1.{1}.3.fastq.gz'.format(prefix, lib_id)]
                outputs = 'out=' + library_file2[0] + ' out2=' + library_file2[1] + ' outs=' + library_file2[2]

            if parameters['noTrim'] == False :
                bb_run = Popen('{bbduk} -Xmx{memory} threads=8 rref={adapters} overwrite=t qout=33 k=23 mink=13 minlength=23 tbo=t entropy=0.75 entropywindow=25 mininsert=23 maxns=2 ktrim=r trimq={read_qual} {read} {outputs}'.format( \
                                  read=reads, outputs=outputs, **parameters).split(), stdout=PIPE, stderr=PIPE)
                timer = Timer(3600, kill_child_proc, [bb_run])
                try:
                    timer.start()
                    bb_out = bb_run.communicate()
                finally: 
                    timer.cancel()
                if bb_run.returncode == 0 :
                    new_reads.append(library_file2)
                    try:
                        for fname in library_file :
                            os.unlink(fname)
                        stat = re.findall('Result:\s+(\d+) reads .+\s+(\d+) bases', bb_out[1])[0]
                        logger('Obtained {1} bases in {0} reads after BBDuk2'.format(*stat))
                    except :
                        pass
                else :
                    new_reads.append(library_file)
                    try:
                        stat = re.findall('Input:\s+(\d+) reads .+\s+(\d+) bases', bb_out[1])[0]
                        logger('BBDuk2 failed! Use original reads with {1} bases in {0} reads'.format(*stat))
                        for fname in library_file2 :
                            os.unlink(fname)
                    except :
                        pass
            else :
                new_reads.append(library_file)
        return new_reads
    def reduce_depth(self, reads) :
        encode = {'A':0, 'C':1, 'G':2, 'T':3}
        read_stats = [[] for library in reads]
        new_reads = [[] for library in reads]
        for lib_id, (library, stat, new_lib) in enumerate(zip(reads, read_stats, new_reads)) :
            for fname in library :
                p = Popen("zcat {0}|awk 'NR%4==2'|wc".format(fname), shell=True, stdout=PIPE).communicate()[0].strip().split()
                n_base = int(p[2]) - int(p[1])
                bcomp = [[0, 0, 0, 0, 0] for i in range(10)]
                p = Popen("zcat {0}|head -400000|awk 'NR%20==2'".format(fname), shell=True, stdout=PIPE, stderr=PIPE)
                for line in p.stdout :
                    for b, bc in zip(line[:10], bcomp) :
                        bc[encode.get(b, 4)] += 1
                seq_start = 0
                for c in range(9, -1, -1) :
                    bc = bcomp[c]
                    if max(bc) / 0.8 >= sum(bc) or (c < 2 and bc[4] > 0.1*sum(bc)) :
                        seq_start = c + 1
                        break
                stat.append([n_base, seq_start])
            n_base = sum([s[0] for s in stat])
            sample_freq = float(parameters['max_base'])/n_base if parameters['max_base'] > 0 else 1.
            if sample_freq >= 1 or len(stat) < 3 :
                sample_freqs = [ sample_freq for s in stat ]
            else :
                n_base2 = sum([s[0] for s in stat[:2]])
                if float(parameters['max_base']) <= n_base2 :
                    sample_freqs = [ float(parameters['max_base'])/n_base2, float(parameters['max_base'])/n_base2, 0. ]
                else :
                    sample_freqs = [ 1., 1., (float(parameters['max_base']) - n_base2)/stat[2][0] ]
            if sample_freqs[0] < 1 and sample_freqs[0] > 0 :
                logger('Read depth too high. Subsample to every {0:.2f} read'.format(1./sample_freqs[0]))

            for f_id, (lib, s, sample_freq) in enumerate(zip(library, stat, sample_freqs)) :
                    if sample_freq > 0 :
                        new_lib.append('{0}.2.{1}.{2}.fastq.gz'.format(parameters['prefix'], lib_id, f_id+1))
                        if parameters['noRename'] == False :
                            if s[1] > 0 :
                                logger('Remove potential barcode bases at the beginning {0} bps of reads in {1}'.format( s[1], lib ))
                                Popen("zcat {0}|awk '{{nr = int((NR-1)/4)}} {{id=(NR-1)%4}} int(nr*{2}) > int((nr-1)*{2}) {{if (id==1 || id == 3) {{print substr($0, {3}, 9999999)}} else {{if(id==0) {{print \"@{4}_\"nr}} else {{print \"+\"}} }} }}'|gzip > {1}".format(
                                    lib, new_lib[-1], min(sample_freq, 1.), s[1]+1, lib_id), shell=True).wait()
                            else :
                                Popen("zcat {0}|awk '{{nr = int((NR-1)/4)}} {{id=(NR-1)%4}} int(nr*{2}) > int((nr-1)*{2}) {{if (id==1 || id == 3) {{print $0}} else {{ if(id==0){{print \"@{4}_\"nr}} else {{print \"+\"}} }} }}'|gzip > {1}".format(
                                    lib, new_lib[-1], min(sample_freq, 1.), s[1]+1, lib_id), shell=True).wait()
                        else :
                            if s[1] > 0 :
                                logger('Remove potential barcode bases at the beginning {0} bps of reads in {1}'.format( s[1], lib ))
                                Popen("zcat {0}|awk '{{nr = int((NR-1)/4)}} {{id=(NR-1)%4}} int(nr*{2}) > int((nr-1)*{2}) {{if (id==1 || id == 3) {{print substr($0, {3}, 9999999)}} else {{if(id==0) {{print $0}} else {{print \"+\"}} }} }}'|gzip > {1}".format(
                                    lib, new_lib[-1], min(sample_freq, 1.), s[1]+1, lib_id), shell=True).wait()
                            else :
                                Popen("zcat {0}|awk '{{nr = int((NR-1)/4)}} {{id=(NR-1)%4}} int(nr*{2}) > int((nr-1)*{2}) {{if (id==1 || id == 3) {{print $0}} else {{ if(id==0){{print $0}} else {{print \"+\"}} }} }}'|gzip > {1}".format(
                                    lib, new_lib[-1], min(sample_freq, 1.), s[1]+1, lib_id), shell=True).wait()

                    os.unlink(lib)
        return new_reads

reads, prefix, parameters = None, None, None
def prepare(args) :
    global reads, prefix, parameters
    parameters = add_args(args).__dict__
    parameters.update(externals)
    
    reads = []
    for k, vs in zip(('pe', 'se'), (parameters['pe'], parameters['se'])) :
        for v in vs :
            if k == 'pe' :
                reads.append([])
                rnames = v.split(',')
                for r1, r2 in zip(rnames[::2], rnames[1::2]) :
                    reads[-1].append([r1, r2])
            elif k == 'se' :
                reads.append([])
                rnames = v.split(',')
                for rn in rnames :
                    reads[-1].append([rn])
    prefix = parameters['prefix']

    logger('Load in {0} read files from {1} libraries'.format(sum([ len(rr) for lib in reads for rr in lib ]), len(reads)))
    reads = preprocess().launch(reads)
    for lib_id, lib in enumerate(reads) :
        for r_id, read in enumerate(lib) :
            new_read = '{0}_L{1}_R{2}.fastq.gz'.format(parameters['prefix'], lib_id+1, r_id+1)
            lib[r_id] = new_read
            os.rename(read, new_read)
    #import json
    #print json.dumps(reads, sort_keys=True, indent=2)
    report = []
    for lib in reads :
        if len(lib) > 1 : 
            report.extend(['--pe', '{0},{1}'.format(*lib)])
        if len(lib) % 2 == 1 :
            report.extend(['--se', '{0}'.format(lib[-1])])
    print ' '.join(report)

def add_args(a) :
    parser = argparse.ArgumentParser(description='''
EToKi.py prepare 
(1) Concatenates reads of the same library together.
(2) Trims sequences based on base-qualities.
(3) Removes potential adapters and barcodes. 
(4) Limits total amount of reads to be used.
(5) Renames reads using their indexes.
''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--pe', action='append', help='sets of PE reads from the same library, delimited by commas. \ne.g.\nIf you have two sets of reads from a same library, it writes like:\n--pe a_R1.fq.gz,a_R2.fq.gz,b_R1.fq.gz,b_R2.fq.gz\nSpecify this multiple times if there are multiple libraries. ', default=[])
    parser.add_argument('--se', action='append', help='sets of SE reads from the same library, delimited by commas. \ne.g.\nIf you have two sets of reads from a same library, it writes like:\n--se a.fq.gz,b.fq.gz\nSpecify this multiple times if there are multiple libraries. ', default=[])
    parser.add_argument('-p', '--prefix', help='prefix for the outputs. Default: EnPrep', default='EnPrep')
    parser.add_argument('-q', '--read_qual', help='Minimum quality for trimming. Default: 6', type=int, default=6)
    parser.add_argument('-b', '--max_base', help='Total amount of bases to be kept. Set -1 to no restriction. Default: 600000000', type=int, default=600000000)
    parser.add_argument('-m', '--memory', help='maximum amount of memory to be used in bbduk2. Default: 30g', default='30g')
    parser.add_argument('--noTrim', help='Do not do quality trim using bbduk2', action='store_true', default=False)
    parser.add_argument('--noRename', help='Do not rename reads', action='store_true', default=False)

    return parser.parse_args(a)

if __name__ == '__main__' :
    prepare(sys.argv[1:])

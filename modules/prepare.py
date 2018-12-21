import os, io, sys, re, shutil, numpy as np, signal, psutil, argparse
from glob import glob
from subprocess import Popen, PIPE, STDOUT
from time import sleep
from threading import Timer
try:
    from .configure import externals, logger, readFasta
except :
    from configure import externals, logger, readFasta

def kill_child_proc (p) :
    for child in psutil.Process(p.pid).children() :
        child.terminate()

def monitor_proc (p) :
    timer = Timer(7200, kill_child_proc, [p])
    output = None
    try:
        timer.start()
        output = p.communicate()
    finally: 
        timer.cancel()
    return p, output
    

# preprocessing
class preprocess(object) :
    def launch (self, reads) :
        reads = self.init_cleanup(reads)
        reads = self.reduce_depth(reads)
        return reads

    def init_cleanup(self, read_libraries) :
        prefix = parameters['prefix']
        new_reads = []
        for lib_id, library in enumerate(read_libraries) :
            if len(library[0]) > 1 :
                library_file = {'PE':['{0}.0.{1}.1.fastq.gz'.format(prefix, lib_id), '{0}.0.{1}.2.fastq.gz'.format(prefix, lib_id)]}
            else :
                library_file = {'SE':['{0}.0.{1}.1.fastq.gz'.format(prefix, lib_id)]}

            if len(library_file.get('PE', [])) :
                Popen('cat {0} > {1}'.format(' '.join([run[0] for run in library]), library_file['PE'][0]), shell=True).wait()
                Popen('cat {0} > {1}'.format(' '.join([run[1] for run in library]), library_file['PE'][1]), shell=True).wait()
                if parameters['merge'] :
                    library_file2 = {'MP':['{0}.1.{1}.m.fastq.gz'.format(prefix, lib_id)], 'PE':['{0}.1.{1}.1.fastq.gz'.format(prefix, lib_id), '{0}.1.{1}.2.fastq.gz'.format(prefix, lib_id)]}
                    reads = 'in={0} in2={1}'.format(*library_file['PE'])
                    outputs = 'out={0} outu1={1} outu2={2}'.format(library_file2['MP'][0], *library_file2['PE'])
                    bb_run, bb_out = monitor_proc(
                        Popen('{bbmerge} -Xmx{memory} threads=8 loose=t mininsert=25 mininsert0=23 qtrim2=t overwrite=t qout=33 entropy=t maxns=2 trimq={read_qual} {read} {outputs}'.format( \
                            read=reads, outputs=outputs, **parameters).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True)
                    )
                    if bb_run.returncode == 0 :
                        for fname in library_file['PE'] :
                            try:
                                os.unlink(fname)
                            except :
                                pass
                        library_file.update(library_file2)
                    else :
                        for fname in library_file2['MP'] + library_file2['PE'] :
                            try:
                                os.unlink(fname)
                            except :
                                pass                

                # do PE bbduk
                if not parameters['noTrim'] :
                    library_file2 = {'PE':['{0}.1.{1}.r1.fastq.gz'.format(prefix, lib_id), '{0}.1.{1}.r2.fastq.gz'.format(prefix, lib_id)], 'SE':['{0}.1.{1}.3.fastq.gz'.format(prefix, lib_id)]}
                    reads = 'in={0} in2={1}'.format(*library_file['PE'])
                    outputs = 'out={1} out2={2} outs={0}'.format(library_file2['SE'][0], *library_file2['PE'])
                    bb_run, bb_out = monitor_proc(
                        Popen('{bbduk} -Xmx{memory} threads=8 ref={adapters} ktrim=r overwrite=t qout=33 k=23 mink=13 minlength=23 tbo=t entropy=0.75 entropywindow=25 mininsert=23 maxns=2 trimq={read_qual} {read} {outputs}'.format( \
                            read=reads, outputs=outputs, **parameters).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True)
                    )
                    if bb_run.returncode == 0 :
                        for fname in library_file['PE'] :
                            try:
                                os.unlink(fname)
                            except :
                                pass
                        library_file.update(library_file2)
                    else :
                        for fname in library_file2['PE'] + library_file2['SE'] :
                            try:
                                os.unlink(fname)
                            except :
                                pass
            else :
                Popen('cat {0} > {1}'.format(' '.join([run[0] for run in library]), library_file['SE'][0]), shell=True).wait()
                if parameters['noTrim'] == False :
                    library_file2 = {'SE':['{0}.1.{1}.s.fastq.gz'.format(prefix, lib_id)]}
                reads = 'in=' + library_file['SE'][0]
                outputs = 'out=' + library_file2['SE'][0]

                bb_run, bb_out = monitor_proc(
                    Popen('{bbduk} -Xmx{memory} threads=8 ref={adapters} ktrim=r overwrite=t qout=33 k=23 mink=13 minlength=23 tbo=t entropy=0.75 entropywindow=25 mininsert=23 maxns=2 trimq={read_qual} {read} {outputs}'.format( \
                        read=reads, outputs=outputs, **parameters).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True)
                )
                if bb_run.returncode == 0 :
                    for fname in library_file['SE'] :
                        try:
                            os.unlink(fname)
                        except :
                            pass
                    library_file.update(library_file2)
                else :
                    
                    for fname in library_file2['SE'] :
                        try:
                            os.unlink(fname)
                        except :
                            pass
            new_reads.append(library_file)
        return new_reads
    def reduce_depth(self, reads) :
        encode = {'A':0, 'C':1, 'G':2, 'T':3}
        read_stats = [{} for library in reads]
        new_reads = [{} for library in reads]
        for lib_id, (libraries, stat, new_libs) in enumerate(zip(reads, read_stats, new_reads)) :
            read_information = [0, 0]
            for lib_type, library in libraries.items() :
                stat[lib_type] = []
                for fname in library :
                    p = Popen("pigz -cd {0}|awk 'NR%4==2'|wc".format(fname), shell=True, stdout=PIPE, universal_newlines=True).communicate()[0].strip().split()
                    n_base, n_read = int(p[2]) - int(p[1]), int(p[0])
                    read_information[0] += n_base
                    read_information[1] += n_read
                    bcomp = [[0, 0, 0, 0, 0] for i in range(10)]
                    p = Popen("pigz -cd {0}|head -200000|awk 'NR%20==2'".format(fname), shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True)
                    for line in p.stdout :
                        for b, bc in zip(line[:10], bcomp) :
                            bc[encode.get(b, 4)] += 1
                    seq_start = 0
                    for c in range(9, -1, -1) :
                        bc = bcomp[c]
                        if max(bc) / 0.8 >= sum(bc) or (c < 2 and bc[4] > 0.1*sum(bc)) :
                            seq_start = c + 1
                            break
                    stat[lib_type].append([n_base, seq_start])
            logger('Obtained {1} bases in {2} reads after Trimming in Lib {0}'.format(lib_id, *read_information))
            n_base = read_information[0]
            sample_freq2 = float(parameters['max_base'])/n_base if parameters['max_base'] > 0 else 1.
            if sample_freq2 >= 1 :
                for ss in stat.values() :
                    for s in ss :
                        s.append(sample_freq2)
            else :
                max_base = float(parameters['max_base'])
                for lib_type in ('MP', 'PE', 'SE') :
                    ss = stat[lib_type]
                    n_base = sum([s[0] for s in ss])
                    sample_freq = float(max_base)/n_base
                    for s in ss :
                        s.append(sample_freq)
                    max_base = 0. if n_base >= max_base else max_base - n_base
            if 0 < sample_freq2 < 1 :
                logger('Read depth too high. Subsampling.')
            
            for lib_type, library in libraries.items() :
                if stat[lib_type][0][-1] > 0 :
                    if lib_type == 'MP' :
                        new_libs[lib_type] = ['{0}.2.{1}.m.fastq.gz'.format(parameters['prefix'], lib_id)]
                    elif lib_type == 'PE' :
                        new_libs[lib_type] = ['{0}.2.{1}.r1.fastq.gz'.format(parameters['prefix'], lib_id), 
                                              '{0}.2.{1}.r2.fastq.gz'.format(parameters['prefix'], lib_id)]
                    else :
                        new_libs[lib_type] = ['{0}.2.{1}.s.fastq.gz'.format(parameters['prefix'], lib_id)]
                    for f_id, (lib, s, nlib) in enumerate(zip(library, stat[lib_type], new_libs[lib_type])) :
                        sample_freq = s[-1]
                        if parameters['noRename'] == False :
                            if s[1] > 0 :
                                logger('Remove potential barcode bases at the beginning {0} bps of reads in {1}'.format( s[1], lib ))
                                Popen("pigz -cd {0}|awk '{{nr = int((NR-1)/4)}} {{id=(NR-1)%4}} int(nr*{2}) > int((nr-1)*{2}) {{if (id==1 || id == 3) {{print substr($0, {3}, 9999999)}} else {{if(id==0) {{print \"@{4}_\"nr}} else {{print \"+\"}} }} }}'|pigz > {1}".format(
                                    lib, nlib, min(sample_freq, 1.), s[1]+1, lib_id), shell=True).wait()
                            else :
                                Popen("pigz -cd {0}|awk '{{nr = int((NR-1)/4)}} {{id=(NR-1)%4}} int(nr*{2}) > int((nr-1)*{2}) {{if (id==1 || id == 3) {{print $0}} else {{ if(id==0){{print \"@{4}_\"nr}} else {{print \"+\"}} }} }}'|pigz > {1}".format(
                                    lib, nlib, min(sample_freq, 1.), s[1]+1, lib_id), shell=True).wait()
                        else :
                            if s[1] > 0 :
                                logger('Remove potential barcode bases at the beginning {0} bps of reads in {1}'.format( s[1], lib ))
                                Popen("pigz -cd {0}|awk '{{nr = int((NR-1)/4)}} {{id=(NR-1)%4}} int(nr*{2}) > int((nr-1)*{2}) {{if (id==1 || id == 3) {{print substr($0, {3}, 9999999)}} else {{if(id==0) {{print $0}} else {{print \"+\"}} }} }}'|pigz > {1}".format(
                                    lib, nlib, min(sample_freq, 1.), s[1]+1, lib_id), shell=True).wait()
                            else :
                                Popen("pigz -cd {0}|awk '{{nr = int((NR-1)/4)}} {{id=(NR-1)%4}} int(nr*{2}) > int((nr-1)*{2}) {{if (id==1 || id == 3) {{print $0}} else {{ if(id==0){{print $0}} else {{print \"+\"}} }} }}'|pigz > {1}".format(
                                    lib, nlib, min(sample_freq, 1.), s[1]+1, lib_id), shell=True).wait()
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
    for lib_id, libary in enumerate(reads) :
        for lib_type, lib in libary.items() :
            for r_id, read in enumerate(lib) :
                if lib_type == 'SE' :
                    new_read = '{0}_L{1}_SE.fastq.gz'.format(parameters['prefix'], lib_id+1)
                elif lib_type == 'MP' :
                    new_read = '{0}_L{1}_MP.fastq.gz'.format(parameters['prefix'], lib_id+1)
                else :
                    new_read = '{0}_L{1}_R{2}.fastq.gz'.format(parameters['prefix'], lib_id+1, r_id+1)
                lib[r_id] = new_read
                os.rename(read, new_read)
    report = []
    for libary in reads :
        for lib_type, lib in libary.items() :
            if lib_type == 'PE' : 
                report.extend(['--pe', '{0},{1}'.format(*lib)])
            else :
                report.extend(['--se', '{0}'.format(lib[-1])])
    print(' '.join(report))

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
    parser.add_argument('-p', '--prefix', help='prefix for the outputs. Default: EToKiPrepare', default='EToKiPrepare')
    parser.add_argument('-q', '--read_qual', help='Minimum quality for trimming. Default: 6', type=int, default=6)
    parser.add_argument('-b', '--max_base', help='Total amount of bases (in BPs) to be kept. Set -1 to no restriction. Use ~100X coverage when you know the size of genome. Default: -1', type=int, default=-1)
    parser.add_argument('-m', '--memory', help='maximum amount of memory to be used in bbduk2. Default: 30g', default='30g')
    parser.add_argument('--noTrim', help='Do not do quality trim using bbduk2', action='store_true', default=False)
    parser.add_argument('--merge', help='Try to merge PE reads by their overlaps', action='store_true', default=False)
    parser.add_argument('--noRename', help='Do not rename reads', action='store_true', default=False)

    args = parser.parse_args(a)
    return args

if __name__ == '__main__' :
    prepare(sys.argv[1:])

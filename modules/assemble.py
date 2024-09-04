import os, sys, re, shutil, numpy as np, pandas as pd
from collections import OrderedDict
from time import sleep
from subprocess import Popen, PIPE
try:
    from .configure import externals, logger, readFasta, xrange, uopen, ETOKI
except :
    from configure import externals, logger, readFasta, xrange, uopen, ETOKI

# mainprocess
class mainprocess(object) :
    def launch (self, reads, longReads=[]) :
        result = parameters.pop('reference', None)
        result = os.path.abspath(result) if result else None

        self.cwd = os.getcwd()
        self.folder = prefix
        if not os.path.exists(self.folder) :
            os.makedirs(self.folder)
        assert os.path.isdir(self.folder), '"{0}" exists and is not a folder'.format(self.folder)
        os.chdir(self.folder)

        self.snps = None
        if not result and len(reads) > 0 :
            if parameters['assembler'] == 'spades' :
                result = self.do_spades(reads, parameters['nohammer'])
            else :
                result = self.do_megahit(reads)
        if sum([len(lr) for lr in longReads]) > 0 :
            longReads = self.rename_longReads(longReads)
            result = self.do_flye(longReads, result, parameters['metagenome'])

        if parameters['outgroup'] or parameters['excluded'] :
            parameters['excluded'] = self.identify_outgroups(result, reads, parameters['ingroup'], parameters['outgroup'], parameters['excluded'], parameters)
        for ite in np.arange(int(parameters['numPolish'])) :
            result, n_changes = self.do_polish(result, reads, parameters['onlySNP'])
            if not n_changes :
                break

        if not parameters['noQuality'] :
            result = self.get_quality(result, reads)
        
        result = os.path.abspath(result)
        os.chdir(self.cwd)
        logger('Final result is written into {0}'.format('{0}.result.fastq'.format(prefix)))
        shutil.move(result, '{0}.result.fastq'.format(prefix))
        return '{0}.result.fastq'.format(prefix)

    def rename_longReads(self, longReads) :
        for i0, libs in enumerate(longReads) :
            for i1, lib in enumerate(libs) :
                for i2, longRead in enumerate(lib) :
                    revised = 'LR_{0}_{1}_{2}.fastq.gz'.format(i0, i1, i2)
                    with uopen(longRead, 'rt') as fin, uopen(revised, 'wt') as fout :
                        for lid, line in enumerate(fin) :
                            if lid % 4 == 0 :
                                fout.write('@{0}\n'.format(int(lid/4)))
                            else :
                                fout.write(line)
                    lib[i2] = revised
        return longReads

    def __get_read_len(self, reads) :
        read_len = 0.
        for lib_id, lib in enumerate(reads) :
            rl = [0, 0]
            for rname in lib :
                p = Popen("{pigz} -cd {0}|head -30000000|awk 'NR%4 == 2'|wc".format(rname, **parameters), shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()[0].split()
                rl[0] += int(p[0])
                rl[1] += int(p[2]) - int(p[0])
            read_len = max(rl[1]/float(rl[0]), read_len) if float(rl[0]) > 0 else read_len
        logger('Estimated read length: {0}'.format(read_len))
        return read_len
    
    def __run_minimap(self, prefix, reference, reads, clean=True) :
        if not os.path.isfile(reference+'.mmi') or (os.path.getmtime(reference+'.mmi') < os.path.getmtime(reference)) :
            Popen('{minimap2} -k13 -w5 -d {0}.mmi {0}'.format(reference, **parameters).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()
        else :
            sleep(1)
    
        outputs = []
        for lib_id, lib in enumerate(reads) :
            se = [ lib[-1], '{0}.mapping.{1}.se1.bam'.format(prefix, lib_id), '{0}.mapping.{1}.se.bam'.format(prefix, lib_id) ] if len(lib) != 2 else [None, None, None]
            pe = [ '{0} {1}'.format(*lib[:2]), '{0}.mapping.{1}.pe1.bam'.format(prefix, lib_id), '{0}.mapping.{1}.pe.bam'.format(prefix, lib_id) ] if len(lib) >= 2 else [None, None, None]

            for r, o1, o in (se, pe) :
                if r is not None :
                    logger('Run minimap2 with: {0}'.format(r))
                    if clean :
                        cmd = '{minimap2} -t{n_cpu} -ax sr --sr --frag=yes -A2 -B4 -O8,16 -E2,1 -r50 -p.6 -N 1 -f2000,10000 -Y -n1 -m19 -s40 -g200 -2K10m --heap-sort=yes --secondary=yes {reference}.mmi {r} |{enbler_filter} {max_diff} {excluded} | {samtools} fixmate -m -@{n_cpu} - - | {samtools} sort -m 1G -@ {n_cpu} -O bam -l 0 -T {prefix} - | {samtools} markdup -r -@ {n_cpu} -O BAM - {o}'.format(
                                r = r, o = o, reference = reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        for line in st_run[1].split('\n') :
                            logger(line.rstrip())
                        outputs.append(o) #self.__markDuplicates(o1, o))
                        Popen('{samtools} index {o}'.format(o=o, **parameters).split(), stdout=PIPE).communicate()
                    else :
                        cmd = '''{minimap2} -t{n_cpu} -ax sr --sr --frag=yes -A2 -B4 -O8,16 -E2,1 -r50 -p.6 -N 1 -f2000,10000 -Y -n1 -m19 -s40 -g200 -2K10m --heap-sort=yes --secondary=yes {reference}.mmi {r} | {enbler_filter} -1 | {samtools} view -bo {o} -'''.format(
                                r = r, o = o, reference = reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        outputs.append(o)
        return outputs

    def __run_bowtie(self, prefix, reference, reads, clean=True) :
        if not os.path.isfile(reference + '.4.bt2') or (os.path.getmtime(reference + '.4.bt2') < os.path.getmtime(reference)) :
            Popen('{bowtie2build} {reference} {reference}'.format(reference=reference, bowtie2build=parameters['bowtie2build']).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
        else :
            sleep(1)

        outputs = []
        for lib_id, lib in enumerate(reads) :
            se = [ '-U ' + lib[-1], '{0}.mapping.{1}.se1.bam'.format(prefix, lib_id), '{0}.mapping.{1}.se.bam'.format(prefix, lib_id) ] if len(lib) != 2 else [None, None, None]
            pe = [ '-1 {0} -2 {1}'.format(*lib[:2]), '{0}.mapping.{1}.pe1.bam'.format(prefix, lib_id), '{0}.mapping.{1}.pe.bam'.format(prefix, lib_id) ] if len(lib) >= 2 else [None, None, None]
            for r, o1, o in (se, pe) :
                if r is not None :
                    logger('Run Bowtie2 with: {0}'.format(r))
                    if clean :
                        cmd= '{bowtie2} -p {n_cpu} --no-unal --mp 4,4 --np 4 --sensitive-local -q -I 25 -X 800 -x {reference} {r} | {enbler_filter} {max_diff} {excluded} | {samtools} fixmate -m -@{n_cpu} - - | {samtools} sort -m 4G -@{n_cpu} -O bam -l 0 -T {prefix} - | {samtools} markdup -r -@{n_cpu} -O BAM - {o}'.format(
                            r=r, o=o, reference=reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        for line in st_run[1].split('\n') :
                            logger(line.rstrip())
                        outputs.append(o)#self.__markDuplicates(o1, o))
                        Popen('{samtools} index {o}'.format(o=o, **parameters).split(), stdout=PIPE).communicate()
                    else :
                        cmd = '''{bowtie2} -p {n_cpu} --no-unal --mp 4,4 --np 4 --sensitive-local -q -I 25 -X 800 -x {reference} {r} | awk '$2 %8 < 4'| {samtools} view -bo {o} -'''.format(
                                r = r, o = o, reference = reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        outputs.append(o)

        return outputs

    def __run_bwa(self, prefix, reference, reads, clean=True) :
        if not os.path.isfile(reference + '.bwt') or (os.path.getmtime(reference + '.bwt') < os.path.getmtime(reference)) :
            Popen('{bwa} index {reference}'.format(reference=reference, bwa=parameters['bwa']).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
        else :
            sleep(1)

        outputs = []
        for lib_id, lib in enumerate(reads) :
            se = [ lib[-1], '{0}.mapping.{1}.se1.bam'.format(prefix, lib_id), '{0}.mapping.{1}.se.bam'.format(prefix, lib_id) ] if len(lib) != 2 else [None, None, None]
            pe = [ '{0} {1}'.format(*lib[:2]), '{0}.mapping.{1}.pe1.bam'.format(prefix, lib_id), '{0}.mapping.{1}.pe.bam'.format(prefix, lib_id) ] if len(lib) >= 2 else [None, None, None]
            for r, o1, o in (se, pe) :
                if r is not None :
                    logger('Run bwa with: {0}'.format(r))
                    if clean :
                        cmd= '{bwa} mem -A 2 -B 4 -T 40 -t {n_cpu} -m 40 {reference} {r} | {enbler_filter} {max_diff} {excluded} | | {samtools} fixmate -m -@{n_cpu} - - | {samtools} sort -m 4G -@{n_cpu} -O bam -l 0 -T {prefix} - | {samtools} markdup -r -@{n_cpu} -O BAM - {o}'.format(
                            r=r, o=o, reference=reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        for line in st_run[1].split('\n') :
                            logger(line.rstrip())
                        outputs.append(o) #self.__markDuplicates(o1, o))
                        Popen('{samtools} index {o}'.format(o=o, **parameters).split(), stdout=PIPE).communicate()
                    else :
                        cmd = '''{bwa} mem -A 2 -B 4 -T 40 -t {n_cpu} -m 40 {reference} {r} | awk '$2 %8 < 4'| {samtools} view -bo {o} -'''.format(
                                r = r, o = o, reference = reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        outputs.append(o)

        return outputs


    def do_megahit(self, reads) :
        outdir = 'megahit'
        output_file = 'megahit.fasta'
        if os.path.isdir(outdir) :
            shutil.rmtree(outdir)
        read_input = [[], [], []]
        for lib in reads :
            if len(lib) % 2 == 1 :
                read_input[2].append(lib[-1])
            else :
                read_input[0].append(lib[0])
                read_input[1].append(lib[1])

        if len(read_input[0]) > 0 :
            if len(read_input[2]) > 0 :
                read_input = '-1 {0} -2 {1} -r {2}'.format(','.join(read_input[0]), ','.join(read_input[1]), ','.join(read_input[2]))
            else :
                read_input = '-1 {0} -2 {1}'.format(','.join(read_input[0]), ','.join(read_input[1]), ','.join(read_input[2]))
        else :
            read_input = '-r {2}'.format(','.join(read_input[0]), ','.join(read_input[1]), ','.join(read_input[2]))
        cmd = '{megahit} {read_input} --k-min 21 --k-max 201 --k-step 18 -t {n_cpu} -m 0.9 -o {outdir}'.format(
              read_input=read_input, outdir=outdir, **parameters)
        run = Popen( cmd.split(), stdout=PIPE, bufsize=0, universal_newlines=True )
        run.communicate()
        if run.returncode != 0 :
            sys.exit(7351685)
        shutil.copyfile( '{outdir}/final.contigs.fa'.format(outdir=outdir), output_file )
        logger('MEGAHIT contigs in {0}'.format(output_file))
        return output_file

    def _flye_dedup(self, input_name, output_name, assembly_info='') :
        # special: remove large duplicates generated by flye
        circular = {}
        if os.path.isfile(assembly_info) :
            with open(assembly_info) as fin :
                fin.readline()
                for line in fin :
                    part = line.strip().split()
                    if part[3] in ('+', 'Y') :
                        circular[part[0]] = 1
        cmd = '{makeblastdb} -dbtype nucl -in {0}'.format(input_name, **parameters)
        Popen(cmd, shell=True).communicate()
        cmd = '{blastn} -num_threads {n_cpu} -db {0} -query {0} -outfmt "6 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue score qlen slen"'.format(input_name, **parameters)
        p = Popen(cmd, shell=True, universal_newlines=True, stdout=PIPE)
        matches = pd.read_csv(p.stdout, sep='\t', header=None).values
        matches = matches[(matches.T[2] >= 98.5) & (matches.T[3] >= 1000)]
        matches = matches[(matches.T[0] != matches.T[1]) | (matches.T[6] != matches.T[8]) | (matches.T[7] != matches.T[9])]
        matches.T[12] = np.min([matches.T[12] - matches.T[7], matches.T[6] - 1], 0)
        matches = matches[matches.T[12] < 500]
        matches = matches[np.argsort(-matches.T[3])]
        matches = matches[np.argsort(matches.T[12])]
        
        toRemove = []
        for mat in matches :
            if mat[0] in circular :
                continue
            r0 = [mat[0], mat[6], mat[7]]
            r1 = [mat[1]] + sorted([mat[8], mat[9]])
            if r0[0] == r1[0] :
                s, e = max(r1[1], r0[1]), min(r1[2], r0[2])
                if (e - s + 1) > 0 :
                    continue
                
            for r2 in toRemove :
                if r1[0] == r2[0] :
                    s, e = max(r1[1], r2[1]), min(r1[2], r2[2])
                    if (e-s+1) >= 0.8 * (r1[2] -r1[1]+1) :
                        r0 = []
                        break
                if r0[0] == r2[0] :
                    s, e = max(r0[1], r2[1]), min(r0[2], r2[2])
                    if (e-s+1) > 0 :
                        r2[1:] = [min(r0[1], r2[1]), max(r0[2], r2[2])]
                        r0 = []
                        break
            if r0 :
                toRemove.append(r0)
        toRemove.sort(reverse=True)

        seq = OrderedDict()
        with open(input_name) as fin :
            for line in fin :
                if line.startswith('>') :
                    name = line[1:].strip()
                    seq[name] = []
                else :
                    seq[name].extend(line.strip().split())
        for n, s in seq.items() :
            seq[n] = list(''.join(s))
        for i, (n, s, e) in enumerate(toRemove) :
            if e < len(seq[n]) :
                seq[n+'.{0}'.format(i)] = seq[n][e:]
            seq[n][(s-1):] = []
 
        for n in circular :
            s = seq.pop(n, None)
            if s :
                seq[n+'_circular'] = s

        with open(output_name, 'w') as fout :
            for n, s in sorted(seq.items()) :
                if len(s) > 100 :
                    fout.write('>{0}\n{1}\n'.format(n, ''.join(s)))
        return output_name

    def do_flye(self, reads, contigs, isMetagenome) :
        isMetagenomes = ['--meta'] if isMetagenome else ['', '--meta']
        if contigs :
            isMetagenomes.append('{0} --meta'.format(contigs))
        outdir = 'flye'
        output_file = 'flye.fasta'
        
        if os.path.isdir(outdir) :
            shutil.rmtree(outdir)
        
        finished = False
        for isMetagenome in isMetagenomes :
            if len(reads[0]) and len(reads[1]) :
                read_inputs = [ '--pacbio-raw {0}'.format(' '.join([r for read in reads for r2 in read for r in r2])), \
                              '--pacbio-raw {0}'.format(' '.join([r for read in reads[0] for r in read])) ]
                cmds = ['{flye} -t {n_cpu} -g 5m --iterations 0 --plasmids {read_input} {isMetagenome} -o {outdir}'.format( \
                        read_input=read_inputs[0], outdir=outdir, isMetagenome=isMetagenome, **parameters), \
                        '{flye} -t {n_cpu} -g 5m --resume-from polishing --plasmids {read_input} {isMetagenome} -o {outdir}'.format( \
                        read_input=read_inputs[1], outdir=outdir, isMetagenome=isMetagenome, **parameters) ]
                finished = True
                for cmd in cmds :
                    flye_run = Popen( cmd.split(), stdout=PIPE, bufsize=0, universal_newlines=True )
                    flye_run.communicate()
                    if flye_run.returncode != 0:
                        finished = False
                        break
            else :
                if len(reads[0]) and not len(reads[1]) :
                    read_input = '--pacbio-raw {0}'.format(' '.join([r for read in reads[0] for r in read]))
                elif len(reads[1]) and not len(reads[0]) :
                    read_input = '--nano-raw {0}'.format(' '.join([r for read in reads[1] for r in read]))
    
                cmd = '{flye} -t {n_cpu} -g 5m --plasmids {read_input} {isMetagenome} -o {outdir}'.format(
                      read_input=read_input, outdir=outdir, isMetagenome=isMetagenome, **parameters)
    
                flye_run = Popen( cmd.split(), stdout=PIPE, bufsize=0, universal_newlines=True )
                flye_run.communicate()
                if flye_run.returncode == 0 :
                    finished = True
            if finished :
                break
        if not finished :
            sys.exit(20155)
                
        flye_file = self._flye_dedup('{outdir}/assembly.fasta'.format(outdir=outdir), '{outdir}/assembly.dedup.fasta'.format(outdir=outdir), '{outdir}/assembly_info.txt'.format(outdir=outdir))
        
        if contigs :
            # polish by two means. 
            #(1) try combine two assemblies together
            outdir2 = 'flye.hybrid'
            if not os.path.isdir(outdir2) :
                os.makedirs(outdir2)

            asm1 = flye_file
            asm2 = '{outdir2}/assembly.dedup.fasta'.format(outdir2=outdir2)
            n1 = 0
            with open(asm1, 'r') as fin, open(asm2, 'w') as fout :
                for line in fin :
                    if line.startswith('>') :
                        fout.write('>x_{0}'.format(line[1:]))
                        n1 += 1
                    else :
                        fout.write(line)

            cmd = '{flye} -t {n_cpu} -g 5m --plasmids --subassemblies {asm} -o {outdir}'.format(
                  asm=' '.join([contigs, asm1, asm2]), outdir=outdir2, **parameters)

            flye_run = Popen( cmd.split(), stdout=PIPE, bufsize=0, universal_newlines=True)
            flye_run.communicate()
            if flye_run.returncode != 0 :
                sys.exit(20123)
            asm2 = self._flye_dedup('{outdir2}/assembly.fasta'.format(outdir2=outdir2), '{outdir2}/assembly.dedup.fasta'.format(outdir2=outdir2), '{outdir2}/assembly_info.txt'.format(outdir2=outdir2))
            
            with open(asm2, 'r') as fin :
                for line in fin :
                    if line.startswith('>') :
                        n1 -= 1
            if n1 >= 0 :
                asm1 = asm2

            #(1) polish flye assembly with contig sequences 
            outdir3 = 'flye.polish'
            if not os.path.isdir(outdir3) :
                os.makedirs(outdir3)

            asm2 = '{outdir3}/contigs.fasta'.format(outdir3=outdir3)
            asm3 = '{outdir3}/flye_contig.fasta'.format(outdir3=outdir3)
            with open(asm1, 'r') as fin, open(asm3, 'w') as fout :
                for line in fin :
                    if line.startswith('>') :
                        fout.write('>x_{0}'.format(line[1:]))
                    else :
                        fout.write(line)
            
            with open(asm2, 'w') as fout :
                for ite in (0,1,2,) :
                    with open(contigs, 'r') as fin :
                        for line in fin :
                            if line.startswith('>') :
                                fout.write('>x{1}_{0}'.format(line[1:], ite))
                            else :
                                fout.write(line)
            
            cmd = '{flye} -t {n_cpu} -g 5m --plasmids --subassemblies {asm} --polish-target {asm1} -o {outdir3}'.format(
                  asm=' '.join([contigs, asm1, asm2, asm3]), asm1=asm1, outdir3=outdir3, **parameters)
            flye_run = Popen( cmd.split(), stdout=PIPE, bufsize=0, universal_newlines=True)
            flye_run.communicate()
            if flye_run.returncode != 0 :
                sys.exit(20123)
            flye_file = '{outdir3}/polished_1.fasta'.format(outdir3=outdir3)
        shutil.copyfile( flye_file, output_file )
        logger('Flye assembly in {0}'.format(output_file))
        return output_file

    def do_spades(self, reads, nohammer) :
        outdir = 'spades'
        output_file = 'spades.fasta'
        if os.path.isdir(outdir) :
            shutil.rmtree(outdir)
        #read_len = min(self.__get_read_len(reads), 140)
        #kmer = ','.join([str(max(min(int(read_len*float(x)/200)*2+1, 127),17)) for x in parameters['kmers'].split(',')])
        read_input = []
        for lib_id, lib in enumerate(reads) :
            if len(lib) == 1 :
                read_input.append('--pe-s {0} {1}'.format(lib_id+1, lib[0]))
            elif len(lib) == 2 :
                read_input.append('--pe-1 {0} {1} --pe-2 {0} {2}'.format(lib_id+1, lib[0], lib[1]))
            elif len(lib) == 3 :
                read_input.append('--pe-1 {0} {1} --pe-2 {0} {2} --pe-s {0} {3}'.format(lib_id+1, lib[0], lib[1], lib[2]))

        cmd = '{python} {spades} -t {n_cpu} {read_input} -o {outdir}'.format(
              python=sys.executable, read_input=' '.join(read_input), outdir=outdir, **parameters)
        if nohammer :
            cmd += ' --only-assembler'
        print(cmd)
        spades_run = Popen( cmd.split(' '), stdout=PIPE, bufsize=0, universal_newlines=True)
        spades_run.communicate()
        if spades_run.returncode != 0 :
            sys.exit(20123)
        try :
            shutil.copyfile( '{outdir}/scaffolds.fasta'.format(outdir=outdir), output_file )
        except :
            shutil.copyfile( '{outdir}/contigs.fasta'.format(outdir=outdir), output_file )
        logger('SPAdes scaffolds in {0}'.format(output_file))
        return output_file

    def do_polish_with_SNPs(self, reference, snp_file) :
        sequence = readFasta(reference)
        snps = { n:[] for n in sequence }
        if snp_file != '' :
            with open(snp_file) as fin :
                for line in fin :
                    part = line.strip().split()
                    snps[part[0]].append([int(part[1]), part[-1]])
            self.snps = snps

        for n, s in sequence.items() :
            sequence[n] = list(s)

        for cont, sites in sorted(snps.items()) :
            for site,base in reversed(sites) :
                if base.startswith('+') :
                    sequence[cont][site-1:site-1] = base[1:]
                elif base.startswith('-') :
                    sequence[cont][site-1:(site+len(base)-2)] = []
                else :
                    sequence[cont][site-1] = base

        with open('etoki.fasta', 'w') as fout :
            for n, s in sorted(sequence.items()) :
                s = ''.join(s)
                fout.write('>{0}\n{1}\n'.format(n, '\n'.join([ s[site:(site+100)] for site in xrange(0, len(s), 100)])))
        return 'etoki.fasta', 0
    
    def identify_outgroups(self, reference, reads, ingroup='', outgroup='', excluded='', parameters={}) :
        excludedReads = {}
        ingroups = [ os.path.abspath(os.path.join(self.cwd, fn)) for fn in ingroup.split(',') if fn != '' ]
        outgroups = [ os.path.abspath(os.path.join(self.cwd, fn)) for fn in outgroup.split(',') if fn != '' ]
        ingroups, outgroups = list( set(ingroups + [reference]) ), list( set(outgroups) - set([reference]) )

        assert np.all([os.path.isfile(fn) for fn in ingroups]), 'Some filenames in ingroup are not valid.'
        assert np.all([os.path.isfile(fn) for fn in outgroups]), 'Some filenames in outgroup are not valid.'
        
        inRef, outRef = 'etoki.inRef', 'etoki.outRef'
        
        if len(outgroups) :
            Popen('cat {0} > {1}'.format(' '.join(outgroups), outRef), shell=True).communicate()
            if parameters['mapper'] == 'minimap2' :
                bams = self.__run_minimap('etoki', outRef, reads, clean=False )
                os.unlink(outRef+'.mmi')
            elif parameters['mapper'] != 'bwa' :
                bams = self.__run_bowtie('etoki', outRef, reads, clean=False )
            else :
                bams = self.__run_bwa('etoki', outRef, reads, clean=False )
            os.unlink(outRef) 
            for bam in bams :
                p = Popen('{samtools} view {0}'.format(bam, **parameters).split(), stdout=PIPE, universal_newlines=True)
                for line in p.stdout :
                    try :
                        rname, score = line.split('\t', 1)[0] , int(re.findall('AS:i:(\d+)', line)[0])
                        if rname not in excludedReads or score > excludedReads[rname] :
                            excludedReads[rname] = score
                    except :
                        continue
        
            if len(ingroups) :
                Popen('cat {0} > {1}'.format(' '.join(ingroups), inRef), shell=True).communicate()
                if parameters['mapper'] == 'minimap2' :
                    bams = self.__run_minimap('etoki', inRef, reads, clean=False)
                    os.unlink(inRef+'.mmi')
                elif parameters['mapper'] != 'bwa' :
                    bams = self.__run_bowtie('etoki', inRef, reads, clean=False )
                else :
                    bams = self.__run_bwa('etoki', inRef, reads, clean=False )
                os.unlink(inRef)
                for bam in bams :
                    p = Popen('{samtools} view {0}'.format(bam, **parameters).split(), stdout=PIPE, universal_newlines=True)
                    for line in p.stdout :
                        rname = line.split('\t', 1)[0]
                        if rname in excludedReads:
                            try :
                                score = int(re.findall('AS:i:(\d+)', line)[0])
                                if score >= excludedReads[rname] :
                                    excludedReads.pop(rname, None)
                            except :
                                continue
        
        if os.path.isfile(excluded) :
            excludedReads.update({ r:9999999999 for r in pd.read_csv(excluded, sep='\t').values.T[0]})
        excludedFile = 'etoki.excludedReads'
        np.savetxt(excludedFile, list(excludedReads.items()), delimiter='\t', fmt='%s', header=[])
        return excludedFile


    def do_polish(self, reference, reads, onlySNP=False):
        if parameters.get('SNP', None) is not None:
            return self.do_polish_with_SNPs(reference, parameters['SNP'])
        else:
            if os.path.isdir('etoki.hapog') :
                shutil.rmtree('etoki.hapog')
            seq = readFasta(reference)
            names = []
            with open('etoki.hapog.ref', 'wt') as fout :
                for id, (n, s) in enumerate(seq.items()) :
                    fout.write('>{0}\n{1}\n'.format(id, s))
                    seq[n] = list(s)
                    names.append(n)
            reference = 'etoki.hapog.ref'
            if parameters['mapper'] == 'minimap2':
                bams = self.__run_minimap('etoki', reference, reads)
            elif parameters['mapper'] != 'bwa':
                bams = self.__run_bowtie('etoki', reference, reads)
            else:
                bams = self.__run_bwa('etoki', reference, reads)

            merged_bam = 'etoki.mapping.merged.bam'
            Popen('{samtools} merge -f {merged_bam} {bams}'.format(
                merged_bam=merged_bam, bams=' '.join(bams), **parameters
            ).split(), stdout=PIPE, universal_newlines=True).communicate()
            Popen('{samtools} index {bam}'.format(bam=merged_bam, **parameters).split()).communicate()

            Popen('{samtools} faidx {0}'.format(reference, **parameters).split()).communicate()

            my_env = os.environ.copy()
            my_env["PATH"] = ETOKI+'/externals' + ':' + my_env["PATH"]
            out, err = Popen('{hapog} -g {0} -u -t {n_cpu} -b {bam} -o etoki.hapog'.format(reference,
                bam=merged_bam, **parameters).split(),
                  universal_newlines=True, stdout=PIPE, stderr=PIPE, env=my_env).communicate()
            if len(err):
                logger(err)
            changes = 0
            try :
                n = Popen('grep read etoki.hapog/hapog_results/hapog.changes'.split(), stdout=PIPE,
                          universal_newlines=True).communicate()
                # Hapog is designed to produce a new assembly with both major and minor alleles
                # we just want an updated assembly which is consistent with the major allele, so update
                # assembly using the changes file which records where there were specific nucleotides
                # found more in the reads than the reference, implying that the reference is the
                # minor allele.  Hapog outputs two ratios, just use the value from the first.
                diffs = [ [p for p in nn.split('\t')] for nn in (n[0].split('\n')) if len(nn) ]
                diffs = [ [names[int(p[0])], int(p[1]), p[2][4:], p[3][5:].upper().replace('-', '')] for p in diffs if float(p[7][7:]) > 0.5 ]
                for n, i, o, r in diffs[::-1] :
                    if not onlySNP or len(r) == 1 :
                        seq[n][i] = r
            except Exception as e:
                logger('Problem while parsing Hapog output {0}'.format(str(e)))
            with open('etoki.fasta', 'wt') as fout :
                for n, s in seq.items():
                    fout.write('>{0}\n{1}\n'.format(n, ''.join(s)))
            return 'etoki.fasta', changes



    def get_ave_depth(self, sites, accurate_depth=False, isMetagenome=False) :
        if True :
            glob_depth = np.zeros(65536, dtype=int)
            for n, s in sites.items() :
                ss = np.bincount(s)[:65536]
                glob_depth[:ss.size] += ss
            total_site = np.sum(glob_depth)
            acc_site, iv = 0, [-1, -1]
            for d, c in enumerate(glob_depth) :
                acc_site += c
                if acc_site >= 0.25 * total_site and iv[0] < 0 :
                    iv[0] = d
                if acc_site >= 0.75 * total_site :
                    iv[1] = d
                    break
            div = max(3, 3*(iv[1] - iv[0]))
            iv = [iv[0]-div, iv[1]+div] if not isMetagenome else [0, 65536]
            sites2 = {}
            for n, s in sites.items() :
                x = s[ (iv[0] <= s) & (s <= iv[1]) ]
                if len(x) :
                    sites2[n] = [s.size, np.mean(s[ (iv[0] <= s) & (s <= iv[1]) ]), 0.]
                else :
                    sites2[n] = [s.size, 0., 0.]

        depth = np.array(list(sites2.values()))
        depth = depth[np.argsort(-depth.T[0])]
        size = np.sum(depth.T[0])
        acc = [0, 0]
        for d in depth :
            acc[0], acc[1] = acc[0] + d[0], acc[1] + d[0]*d[1]
            if acc[0] *2 >= size :
                break
        ave_depth = acc[1]/acc[0]

        return sites2, ave_depth

    def __calc_qual(self, n_ref, n_alt):
        if n_ref < n_alt :
            return 0
        p1, p2 = 0.1, 0.5
        qual = int(10*n_ref*(np.log10(1-p1) - np.log10(1-p2)) + 10*n_alt*(np.log10(p1) - np.log10(p2))+0.5)
        return max(qual, 0)

    def get_quality(self, reference, reads ) :
        if parameters['mapper'] == 'minimap2' :
            bams = self.__run_minimap('etoki', reference, reads, )
        elif parameters['mapper'] != 'bwa' :
            bams = self.__run_bowtie('etoki', reference, reads, )
        else :
            bams = self.__run_bwa('etoki', reference, reads, )

        merged_bam = 'etoki.mapping.merged.bam'
        Popen('{samtools} merge -f {merged_bam} {bams}'.format(
            merged_bam=merged_bam, bams=' '.join(bams), **parameters
        ).split(), stdout=PIPE, universal_newlines=True).communicate()
        Popen('{samtools} index {bam}'.format(bam=merged_bam, **parameters).split()).communicate()

        Popen('{samtools} faidx {0}'.format(reference, **parameters).split()).communicate()

        sequence = readFasta(reference)
        for n, s in sequence.items() :
            q = ['!'] * len(s)
            sequence[n] = [s, q]

        sites = { n:np.array([0 for _ in s[1] ]) for n, s in sequence.items() }
        depth = Popen('{samtools} depth -aa -q 0 -Q 0 {bam}'.format(bam=merged_bam, **parameters).split(), stdout=PIPE, universal_newlines=True)
        for line in depth.stdout :
            part = line.strip().split()
            if len(part) > 2 and float(part[2]) > 0 :
                sites[part[0]][int(part[1]) - 1] = float(part[2])

        sites, ave_depth = self.get_ave_depth(sites, parameters['accurate_depth'], parameters['metagenome'])
        cont_depth = [float(d)*ave_depth for d in parameters['cont_depth'].split(',')]
        logger('Average read depth: {0}'.format(ave_depth))
        logger('Contigs with less than {0} depth will be removed from the assembly'.format(cont_depth[0]))
        for n, s in sorted(sites.items()) :
            s[2] = s[1]/ave_depth

        cmd = '{samtools} mpileup -ABf {0} {1}'.format(reference, merged_bam, **parameters)
        q = Popen( cmd.split(), stdout=PIPE, universal_newlines=True )
        for line in q.stdout:
            p = line.strip().split()
            depth = int(p[3])
            n_ref = p[4].count('.') + p[4].count(',')
            n_alt = depth - n_ref
            if depth > cont_depth[1] :
                n_ref -= depth - cont_depth[1]
            q = self.__calc_qual(n_ref, n_alt)
            sequence[p[0]][1][int(p[1])-1] = chr( min(q, 40) + 33 )


        with open('etoki.result.fastq', 'w') as fout :
            p = prefix.rsplit('/', 1)[-1]
            for n, (s, q) in sorted(sequence.items()) :
                if sites[n][1] >= cont_depth[0] :
                    fout.write( '@{0} {3} {4} {5}\n{1}\n+\n{2}\n'.format( p+'_'+n, s, ''.join(q), *sites[n] ) )
        return 'etoki.result.fastq'
    
# postprocess
class postprocess(object) :
    def launch (self, assembly) :
        try:
            seq, fasfile = self.__readAssembly(assembly)
        except :
            return {}
        evaluation = dict(assembly=assembly, fasta=fasfile)
        evaluation.update(self.do_evaluation(seq))
        if parameters['runKraken'] and os.path.exists(parameters['kraken_database']) :
            evaluation.update({'kraken': self.do_kraken(assembly, seq)})
        return evaluation
    def __readAssembly(self, assembly) :
        seq = {}
        with uopen(assembly) as fin :
            header = fin.read(1)
        with uopen(assembly) as fin :
            if header == '@' :
                for id, line in enumerate(fin) :
                    if id % 4 == 0 :
                        part = line[1:].strip().split()
                        name = part[0]
                        seq[name]= [0, float(part[2]) if len(part) > 2 else 0., None, None]
                    elif id % 4 == 1 :
                        seq[name][2] = line.strip()
                        seq[name][0] = len(seq[name][2])
                    elif id % 4 == 3 :
                        seq[name][3] = np.array(list(line.strip()))
                fasfile = assembly.rsplit('.', 1)[0] + '.fasta'
                logger('Write fasta sequences into {0}'.format(fasfile))
                with open(fasfile, 'w') as fout :
                    for n, s in sorted(seq.items()) :
                        fout.write('>{0}\n{1}\n'.format(n, '\n'.join([ s[2][site:(site+100)] for site in xrange(0, len(s[2]), 100)])))
            else :
                fasfile = assembly
                for id, line in enumerate(fin) :
                    if line.startswith('>') :
                        name = line[1:].strip().split()[0]
                        seq[name] = [0, 0., []]
                    else :
                        seq[name][2].extend( line.strip().split() )
                for n, s in sorted(seq.items()) :
                    s[2] = ''.join(s[2])
                    s[0] = len(s[2])
        return seq, fasfile

    def do_kraken(self, assembly, seq) :
        with open(assembly+'.filter', 'w') as fout :
            for n, s in sorted(seq.items()) :
                if s[0] > 1000 :
                    fout.write('>{0}\n{1}\n'.format(n, s[2]))
        cmd = '{kraken2} -db {kraken_database} --threads {n_cpu} --output - --report {assembly}.kraken {assembly}.filter'.format(
            assembly=assembly, **parameters
        )
            
        Popen(cmd, stderr=PIPE, stdout=PIPE, shell=True, universal_newlines=True).communicate()
        species = {}
        with open('{0}.kraken'.format(assembly)) as fin :
            for line in fin :
                part = line.strip().split('\t')
                if len(part) < 5 :
                    continue
                if part[3] == 'S' :
                    s_name = part[5].strip(' ')
                    if s_name.startswith( 'Shigella') or s_name.startswith('Escherichia coli') :
                        s_name = 'Escherichia coli / Shigella'
                    species[s_name] = species.get(s_name, 0) + float(part[0])
        species = sorted(species.items(), key=lambda x:x[1], reverse=True)
        species_sum = sum([x[1] for x in species])
        species = [[x[0], int(10000.0*x[1]/species_sum + 0.5)/100.0] for x in species]
        try :
            os.unlink('{0}.kraken'.format(assembly))
        except :
            pass
        return species

    def do_evaluation(self, fastq) :
        seq = sorted([s for s in fastq.values() if s[0] >= 300], key=lambda x:-x[0])
        n_seq = len(seq)
        n_base = sum([s[0] for s in seq])
        n50, acc = 0, [0, 0]
        l50, ave_depth, n_low = 0, 0, 0
        if n_seq > 0:
            for l50, s in enumerate(seq) :
                acc[0], acc[1] = acc[0] + s[0], acc[1] + s[0]*s[1]
                if acc[0] * 2 >= n_base :
                    n50 = s[0]
                    break
            l50 += 1
            ave_depth = acc[1]/acc[0]
            for s in seq :
                if len(s) > 3 :
                    n_low += np.sum( (np.vectorize(ord)(s[3]) < 43) | (np.array(list(s[2])) == 'N') )
                else :
                    n_low += np.sum( (np.array(list(s[2])) == 'N') )
        return dict(n_contig = n_seq,
                    n_base = n_base,
                    ave_depth = ave_depth,
                    n_lowQual = float(n_low),
                    N50 = n50,
                    L50 = l50)

reads, prefix, parameters = None, None, None
def assemble(args) :
    global reads, prefix, parameters
    parameters = add_args(args).__dict__
    parameters.update(externals)
    prefix = parameters['prefix']
    
    reads = OrderedDict([['PE', []],['SE', []],['PacBio', []],['ONT', []]])
    for (k, d), vs in zip(reads.items(), (parameters['pe'], parameters['se'], parameters['pacbio'], parameters['ont'])) :
        for v in vs :
            rnames = [os.path.abspath(rn) for rn in v.split(',')]
            if len(rnames) > 0 :
                if k == 'PE' :
                    assert len(rnames) == 2, 'Allows 2 reads per PE library. You specified {0}'.format(len(rnames))
                elif k in ('SE', 'PacBio', 'ONT') :
                    assert len(rnames) == 1, 'Allows one file per {1} library. You specified {0}'.format(len(rnames), k)
                d.append(rnames)

    logger('Load in {0} read files from {1} libraries'.format(sum(len(l) for lib in reads.values() for l in lib), sum(len(lib) for lib in reads.values())))
    if not parameters['onlyEval'] :
        assembly = mainprocess().launch(reads['PE'] + reads['SE'], [reads['PacBio'], reads['ONT']])
    else :
        assembly = parameters['reference']
    
    report = postprocess().launch(assembly)
    import json
    print(json.dumps(report, sort_keys=True, indent=2))

def add_args(a) :
    import argparse
    parser = argparse.ArgumentParser(description='''
EToKi.py assemble
(1.1) Assembles short reads into assemblies, or 
(1.2) Maps them onto a reference. 
And
(2) Polishes consensus using polish, 
(3) Removes low level contaminations. 
(4) Estimates the base quality of the consensus. 
(5) Predicts taxonomy using Kraken.
''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--pe', action='append', help='comma delimited two files of PE reads. ', default=[])
    parser.add_argument('--se', action='append', help='one file of SE read. \n', default=[])
    parser.add_argument('--pacbio', action='append', help='one file of pacbio read. \n', default=[])
    parser.add_argument('--ont', action='append', help='one file of nanopore read. \n', default=[])
    parser.add_argument('-p', '--prefix', help='prefix for the outputs. Default: EToKi_assemble', default='EToKi_assemble')
    parser.add_argument('-a', '--assembler', help='Assembler used for de novo assembly. \nDisabled if you specify a reference. \nDefault: spades for single colony isolates, megahit for metagenome. \n Long reads will always be assembled with Flye', default='')
    parser.add_argument('-r', '--reference', help='Reference for read mapping. Specify this for reference mapping module. ', default=None)
    parser.add_argument('--nohammer', help='Flag to disable read correction in spades. ', default=False, action='store_true')
    parser.add_argument('-k', '--kmers', help='relative lengths of kmers used in SPAdes. Default: 30,50,70,90', default='30,50,70,90')
    parser.add_argument('-m', '--mapper', help='aligner used for read mapping.\noptions are: miminap (default), bwa or bowtie2', default='minimap2')
    parser.add_argument('-d', '--max_diff', help='Maximum proportion of variations allowed for a aligned reads. \nDefault: 0.1 for single isolates, 0.05 for metagenome', type=float, default=-1)
    parser.add_argument('-i', '--ingroup', help='Additional references presenting intra-population genetic diversities. ', default='')
    parser.add_argument('-o', '--outgroup', help='Additional references presenting genetic diversities outside of the studied population. \nReads that are more similar to outgroups will be excluded from analysis. ', default='')
    
    parser.add_argument('-S', '--SNP', help='Exclusive set of SNPs. This will overwrite the polish process. \nRequired format:\n<cont_name> <site> <base_type>\n...', default=None)
    parser.add_argument('-c', '--cont_depth', help='Allowed range of read depth variations relative to average value.\nDefault: 0.2,100\nContigs with read depths outside of this range will be removed from the final assembly.', default='')
    
    parser.add_argument('--excluded', help='A name of the file that contains reads to be excluded from the analysis.', default='')
    parser.add_argument('--metagenome', help='Reads are from metagenomic samples', action='store_true', default=False)
    parser.add_argument('--accurate_depth', help='A better estimation of read depths in each contig', action='store_true', default=True)

    parser.add_argument('--n_cpu', help='Number of CPUs. Default: 8', default=8, type=int)
    parser.add_argument('--numPolish', help='Number of Pilon polish iterations. Default: 3 (single) or 1 (metagenome)', default=-1, type=int)
    parser.add_argument('--onlySNP', help='Only modify substitutions during the PILON polish.', action='store_true', default=False)
    parser.add_argument('--noQuality', help='Do not estimate base qualities.', action='store_true', default=False)
    parser.add_argument('--onlyEval', help='Do not run assembly/mapping. Only evaluate assembly status.', action='store_true', default=False)
    parser.add_argument('--kraken', dest='runKraken', help='Run kmer based species predicton on contigs.', action='store_true', default=False)

    args = parser.parse_args(a)
    if args.cont_depth == '' :
        args.cont_depth = '0.2,100' if not args.metagenome else '0.0001,10000'
    if args.assembler == '' :
        args.assembler = 'spades' if not args.metagenome else 'megahit'
    if args.max_diff < 0 :
        args.max_diff = 0.1 if not args.metagenome else 0.05
    if args.numPolish < 0 :
        args.numPolish = 3 if not args.metagenome else 1

    return args

if __name__ == '__main__' :
    assemble(sys.argv[1:])

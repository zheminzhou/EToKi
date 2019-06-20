import os, io, sys, re, shutil, numpy as np, pandas as pd
from glob import glob
from time import sleep
from subprocess import Popen, PIPE, STDOUT
try:
    from .configure import externals, logger, readFasta, xrange, uopen
except :
    from configure import externals, logger, readFasta, xrange, uopen

# mainprocess
class mainprocess(object) :
    def launch (self, reads) :
        self.snps = None
        result = parameters.pop('reference', None)
        if not result :
            if parameters['assembler'] == 'spades' :
                result = self.do_spades(reads)
            else :
                result = self.do_megahit(reads)
        if parameters['outgroup'] or parameters['excluded'] :
            parameters['excluded'] = self.identify_outgroups(result, parameters['ingroup'], parameters['outgroup'], parameters['excluded'])
        if not parameters['noPolish'] :
            result = self.do_polish(result, reads, parameters['reassemble'], parameters['onlySNP'])
        if not parameters['noQuality'] :
            result = self.get_quality(result, reads)
        return result

    def __get_read_len(self, reads) :
        read_len = 0.
        for lib_id, lib in enumerate(reads) :
            rl = [0, 0]
            for rname in lib :
                p = Popen("pigz -cd {0}|head -100000|awk 'NR%4 == 2'|wc".format(rname), shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()[0].split()
                rl[0] += int(p[0])
                rl[1] += int(p[2]) - int(p[0])
            read_len = max(rl[1]/float(rl[0]), read_len) if float(rl[0]) > 0 else read_len
        logger('Estimated read length: {0}'.format(read_len))
        return read_len
    
    def __markDuplicates(self, ori_bam, tgt_bam) :
        try:
            x = Popen('{gatk} MarkDuplicates -O {output} -M {prefix}.mapping.dup --REMOVE_DUPLICATES true -I {input}'.format( \
                output=tgt_bam, input=ori_bam, \
                **parameters).split(), stdout=PIPE, stderr=PIPE, universal_newlines=True )
            x.communicate()
            assert x.returncode == 0
            os.unlink(ori_bam)
            output = tgt_bam
        except :
            output = ori_bam
        Popen('{samtools} index {output}'.format(output=output, **parameters).split(), stdout=PIPE, universal_newlines=True ).communicate()
        return output
    
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
                        cmd = '{minimap2} -t8 -ax sr --sr --frag=yes -A2 -B4 -O8,16 -E2,1 -r50 -p.6 -N 8 -f2000,10000 -Y -n1 -m19 -s40 -g200 -2K10m --heap-sort=yes --secondary=yes {reference}.mmi {r} |{enbler_filter} {max_diff} {excluded} | {samtools} sort -@ 8 -O bam -l 0 -T {prefix} - > {o}'.format(
                                r = r, o = o1, reference = reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        for line in st_run[1].split('\n') :
                            logger(line.rstrip())
                        outputs.append(self.__markDuplicates(o1, o))
                    else :
                        cmd = '''{minimap2} -t8 -ax sr --sr --frag=yes -A2 -B4 -O8,16 -E2,1 -r50 -p.6 -N 1 -f2000,10000 -Y -n1 -m19 -s40 -g200 -2K10m --heap-sort=yes --secondary=yes {reference}.mmi {r} | awk '$2 %8 < 4'| {samtools} view -bo {o} -'''.format(
                                r = r, o = o1, reference = reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        outputs.append(o1)
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
                        cmd= '{bowtie2} -p 8 --no-unal --mp 4,4 --np 4 --sensitive-local -q -I 25 -X 800 -x {reference} {r} | {enbler_filter} {max_diff} {excluded} | {samtools} sort -@ 8 -O bam -l 0 -T {prefix} - > {o}'.format(
                            r=r, o=o1, reference=reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        for line in st_run[1].split('\n') :
                            logger(line.rstrip())
                        outputs.append(self.__markDuplicates(o1, o))
                    else :
                        cmd = '''{bowtie2} -p 8 --no-unal --mp 4,4 --np 4 --sensitive-local -q -I 25 -X 800 -x {reference} {r} | awk '$2 %8 < 4'| {samtools} view -bo {o} -'''.format(
                                r = r, o = o1, reference = reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        outputs.append(o1)

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
                        cmd= '{bwa} mem -A 2 -B 4 -T 40 -t 8 -m 40 {reference} {r} | {enbler_filter} {max_diff} {excluded} | {samtools} sort -@ 8 -O bam -l 0 -T {prefix} - > {o}'.format(
                            r=r, o=o1, reference=reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        for line in st_run[1].split('\n') :
                            logger(line.rstrip())
                        outputs.append(self.__markDuplicates(o1, o))
                    else :
                        cmd = '''{bwa} mem -A 2 -B 4 -T 40 -t 8 -m 40 {reference} {r} | awk '$2 %8 < 4'| {samtools} view -bo {o} -'''.format(
                                r = r, o = o1, reference = reference, **parameters)
                        st_run = Popen( cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True ).communicate()
                        outputs.append(o1)
        return outputs


    def do_megahit(self, reads) :
        outdir = prefix + '.megahit'
        output_file = prefix + '.megahit.fasta'
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
        cmd = '{megahit} {read_input} --k-max 201 --k-step 10 -t 8 -m 0.9 -o {outdir}'.format(
              megahit=parameters['megahit'], read_input=read_input, outdir=outdir)
        run = Popen( cmd.split(), stdout=PIPE, bufsize=0, universal_newlines=True )
        run.communicate()
        if run.returncode != 0 :
            sys.exit(7351685)
        shutil.copyfile( '{outdir}/final.contigs.fa'.format(outdir=outdir), output_file )
        logger('MEGAHIT contigs in {0}'.format(output_file))
        return output_file

    def do_spades(self, reads) :
        outdir = prefix + '.spades'
        output_file = prefix + '.spades.fasta'
        if os.path.isdir(outdir) :
            shutil.rmtree(outdir)
        read_len = min(self.__get_read_len(reads), 140)
        kmer = ','.join([str(max(min(int(read_len*float(x)/200)*2+1, 127),17)) for x in parameters['kmers'].split(',')])
        read_input = []
        for lib_id, lib in enumerate(reads) :
            if len(lib) == 1 :
                read_input.append('--s{0} {1}'.format(lib_id+1, lib[0]))
            elif len(lib) == 2 :
                read_input.append('--pe{0}-1 {1} --pe{0}-2 {2}'.format(lib_id+1, lib[0], lib[1]))
            elif len(lib) == 3 :
                read_input.append('--pe{0}-1 {1} --pe{0}-2 {2} --pe{0}-s {3}'.format(lib_id+1, lib[0], lib[1], lib[2]))
        cmd = '{spades} -t 8 --only-assembler {read_input} -k {kmer} -o {outdir}'.format(
              spades=parameters['spades'], read_input=' '.join(read_input), kmer=kmer, outdir=outdir)
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

        for cont, sites in snps.items() :
            for site,base in reversed(sites) :
                if base.startswith('+') :
                    sequence[cont][site-1:site-1] = base[1:]
                elif base.startswith('-') :
                    sequence[cont][site-1:(site+len(base)-2)] = []
                else :
                    sequence[cont][site-1] = base

        with open('{0}.fasta'.format(prefix), 'w') as fout :
            for n, s in sorted(sequence.items()) :
                s = ''.join(s)
                fout.write('>{0}\n{1}\n'.format(n, '\n'.join([ s[site:(site+100)] for site in xrange(0, len(s), 100)])))
        return '{0}.fasta'.format(prefix)
    
    def identify_outgroups(self, reference, ingroup='', outgroup='', excluded='') :
        excludedReads = {}
        ingroups, outgroups = list(set(ingroup.split(',') + [reference]) - set([''])), list(set(outgroup.split(',')) - set(['']))
        assert np.all([os.path.isfile(fn) for fn in ingroups]), 'Some filenames in ingroup is not valid.'
        assert np.all([os.path.isfile(fn) for fn in outgroups]), 'Some filenames in outgroup is not valid.'
        
        inRef, outRef = '{0}.inRef'.format(parameters['prefix']), '{0}.outRef'.format(parameters['prefix'])
        
        if len(outgroups) :
            Popen('cat {0} > {1}'.format(' '.join(outgroups), outRef), shell=True).communicate()
            if parameters['mapper'] == 'minimap2' :
                bams = self.__run_minimap(prefix, outRef, reads, clean=False )
            elif parameters['mapper'] != 'bwa' :
                bams = self.__run_bowtie(prefix, outRef, reads, clean=False )
            else :
                bams = self.__run_bwa(prefix, outRef, reads, clean=False )
            
            for bam in bams :
                p = Popen('samtools view {0}'.format(bam).split(), stdout=PIPE, universal_newlines=True)
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
                    bams = self.__run_minimap(prefix, inRef, reads, clean=False )
                elif parameters['mapper'] != 'bwa' :
                    bams = self.__run_bowtie(prefix, inRef, reads, clean=False )
                else :
                    bams = self.__run_bwa(prefix, inRef, reads, clean=False )
                
                for bam in bams :
                    p = Popen('samtools view {0}'.format(bam).split(), stdout=PIPE, universal_newlines=True)
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
        excludedFile = '{0}.excludedReads'.format(parameters['prefix'])
        np.savetxt(excludedFile, list(excludedReads.items()), delimiter='\t', fmt='%s', header=[])
        return excludedFile
    
    def do_polish(self, reference, reads, reassemble=False, onlySNP=False) :
        if parameters.get('SNP', None) is not None :
            return self.do_polish_with_SNPs(reference, parameters['SNP'])
        else :
            if parameters['mapper'] == 'minimap2' :
                bams = self.__run_minimap(prefix, reference, reads )
            elif parameters['mapper'] != 'bwa' :
                bams = self.__run_bowtie(prefix, reference, reads )
            else :
                bams = self.__run_bwa(prefix, reference, reads )
            sites = {}
            for bam in bams :
                if bam is not None :
                    depth = Popen('{samtools} depth -q 0 -Q 0 {bam}'.format(bam=bam, **parameters).split(), stdout=PIPE, universal_newlines=True)
                    for line in depth.stdout :
                        part = line.strip().split()
                        if len(part) > 2 and float(part[2]) > 0 :
                            sites[part[0]] = 1
            sequence = readFasta(reference)
            sequence = {n:s for n,s in sequence.items() if n in sites}

            with open('{0}.mapping.reference.fasta'.format(prefix), 'w') as fout :
                for n, s in sorted(sequence.items()) :
                    fout.write('>{0}\n{1}\n'.format(n, '\n'.join([ s[site:(site+100)] for site in xrange(0, len(s), 100)])))

            bam_opt = ' '.join(['--bam {0}'.format(b) for b in bams if b is not None])
            if reassemble :
                pilon_cmd = '{pilon} --fix all,breaks --vcf --output {prefix}.mapping --genome {prefix}.mapping.reference.fasta {bam_opt}'.format(bam_opt=bam_opt, **parameters)
                Popen( pilon_cmd.split(), stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()
            else :
                pilon_cmd = '{pilon} --fix all --vcf --output {prefix}.mapping --genome {prefix}.mapping.reference.fasta {bam_opt}'.format(bam_opt=bam_opt, **parameters)
                Popen( pilon_cmd.split(), stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()
            
            if not os.path.isfile('{0}.mapping.vcf'.format(prefix)) :
                pilon_cmd = '{pilon} --fix snps,indels,gaps --vcf --output {prefix}.mapping --genome {prefix}.mapping.reference.fasta {bam_opt}'.format(bam_opt=bam_opt, **parameters)
                Popen( pilon_cmd.split(), stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()                    
            
            snps = []
            with open('{0}.mapping.vcf'.format(prefix)) as fin, open('{0}.mapping.changes'.format(prefix), 'w') as fout :
                for line in fin :
                    if line.startswith('#') : continue
                    part = line.strip().split('\t')
                    if part[-1] != '0/0':
                        try :
                            if (part[6] == 'PASS' or float(part[7][-4:]) >= 0.75) and re.match(r'^[ACGTN]+$', part[4]):
                                if (not onlySNP) or (len(part[3]) == 1 and len(part[4]) == 1 ) :
                                    snps.append( [ part[0], int(part[1])-1, part[3], part[4] ] )
                                    fout.write(line)
                        except :
                            pass

            os.unlink('{0}.mapping.vcf'.format(prefix))
            for n in sequence.keys() :
                sequence[n] = list(sequence[n])
            for n, site, ori, alt in reversed(snps) :
                s = sequence[n]
                end = site + len(ori)
                s[site:end] = alt
            logger('Observed and corrected {0} changes using PILON'.format(len(snps)))
            with open('{0}.fasta'.format(prefix), 'w') as fout :
                for n, s in sorted(sequence.items()) :
                    s = ''.join(s)
                    fout.write('>{0}\n{1}\n'.format(n, '\n'.join([ s[site:(site+100)] for site in xrange(0, len(s), 100)])))
            return '{0}.fasta'.format(prefix)

    def get_quality(self, reference, reads ) :
        if parameters['mapper'] == 'minimap2' :
            bams = self.__run_minimap(prefix, reference, reads, )
        elif parameters['mapper'] != 'bwa' :
            bams = self.__run_bowtie(prefix, reference, reads, )
        else :
            bams = self.__run_bwa(prefix, reference, reads, )
        
        sequence = readFasta(reference)
        for n, s in sequence.items() :
            q = ['!'] * len(s)
            sequence[n] = [s, q]

        sites = { n:np.array([0 for ss in s[1] ]) for n, s in sequence.items() }
        for bam in bams :
            if bam is not None :
                depth = Popen('{samtools} depth -q 0 -Q 0 {bam}'.format(bam=bam, **parameters).split(), stdout=PIPE, universal_newlines=True)
                for line in depth.stdout :
                    part = line.strip().split()
                    if len(part) > 2 and float(part[2]) > 0 :
                        sites[part[0]][int(part[1]) - 1] += float(part[2])
        sites = {n:[s.size, np.mean(s), 0.] for n, s in sites.items()}
        depth = np.array(list(sites.values()))
        depth = depth[np.argsort(-depth.T[0])]
        size = np.sum(depth.T[0])
        acc = [0, 0]
        for d in depth :
            acc[0], acc[1] = acc[0] + d[0], acc[1] + d[0]*d[1]
            if acc[0] *2 >= size :
                break
        ave_depth = acc[1]/acc[0]
        exp_mut_depth = max(ave_depth * 0.2, 2.)
        for n, s in sites.items() :
            s[2] = s[1]/ave_depth
        logger('Average read depth: {0}'.format(ave_depth))
        sequence = {n:s for n, s in sequence.items() if sites[n][1]>0.}
        with open('{0}.mapping.reference.fasta'.format(prefix), 'w') as fout :
            for n, s in sorted(sequence.items()) :
                fout.write('>{0}\n{1}\n'.format(n, '\n'.join([ s[0][site:(site+100)] for site in xrange(0, len(s[0]), 100)])))
        bam_opt = ' '.join(['--bam {0}'.format(b) for b in bams if b is not None])
        pilon_cmd = '{pilon} --fix all,breaks --vcf --output {prefix}.mapping --genome {prefix}.mapping.reference.fasta {bam_opt}'.format(bam_opt=bam_opt, **parameters)
        Popen( pilon_cmd.split(), stdout=PIPE, universal_newlines=True ).communicate()
        if not os.path.isfile('{0}.mapping.vcf'.format(prefix)) :
            pilon_cmd = '{pilon} --fix snps,indels,gaps,breaks --vcf --output {prefix}.mapping --genome {prefix}.mapping.reference.fasta {bam_opt}'.format(bam_opt=bam_opt, **parameters)
            Popen( pilon_cmd.split(), stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate()                    

        cont_depth = [float(d) for d in parameters['cont_depth'].split(',')]
        logger('Contigs with less than {0} depth will be removed from the assembly'.format(cont_depth[0]*ave_depth))
        logger('Contigs with more than {0} depth will be treated as duplicates'.format(cont_depth[1]*ave_depth))

        with open('{0}.mapping.vcf'.format(prefix)) as fin, open('{0}.mapping.difference'.format(prefix), 'w') as fout :
            for line in fin :
                if line.startswith('#') : continue
                part = line.strip().split('\t')
                if sites[part[0]][2] < cont_depth[0] or sites[part[0]][2] >= cont_depth[1] :
                    continue
                try:
                    if part[-1] == '0/0' and len(part[3]) == 1 and len(part[4]) == 1 :
                        pp = part[7].split(';')
                        dp = float(pp[0][3:])
                        af = 100 - sorted([float(af) for af in pp[6][3:].split(',')])[-1]
                        if af <= 15 and dp >= 3 and dp * af/100. <= exp_mut_depth :
                            if part[6] == 'PASS' or (part[6] == 'LowCov' and parameters['metagenome']) :
                                site = int(part[1])-1
                                qual = chr(int(pp[4][3:])+33)
                                sequence[part[0]][1][site] = qual
                        else :
                            fout.write(line)
                    else :
                        fout.write(line)
                except :
                    fout.write(line)
        if self.snps is not None :
            for n, snvs in self.snps.items() :
                for site, snv in snvs :
                    if snv.find('N') >= 0 : continue
                    if snv.startswith('+') :
                        s, e = site-4, site+3+len(snv)
                    else :
                        s, e = site-4, site+4
                    for k in xrange(s, e) :
                        sequence[n][1][k] = max(chr(40+33), sequence[n][1][k])

        with open('{0}.result.fastq'.format(prefix), 'w') as fout :
            p = prefix.rsplit('/', 1)[-1]
            for n, (s, q) in sequence.items() :
                if sites[n][2] >= cont_depth[0] :
                    fout.write( '@{0} {3} {4} {5}\n{1}\n+\n{2}\n'.format( p+'_'+n, s, ''.join(q), *sites[n] ) )
        os.unlink( '{0}.mapping.vcf'.format(prefix) )
        logger('Final result is written into {0}'.format('{0}.result.fastq'.format(prefix)))
        return '{0}.result.fastq'.format(prefix)
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
            evaluation.update({'kraken': self.do_kraken(fasfile)})
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
                for n, s in seq.items() :
                    s[2] = ''.join(s[2])
                    s[0] = len(s[2])
        return seq, fasfile

    def do_kraken(self, assembly) :
        if assembly.lower().endswith('.gz') :
            cmd = '{kraken2} -db {kraken_database} --threads 8 --output - --gzip-compressed --report {assembly}.kraken {assembly}'.format(
                assembly=assembly, **parameters
            )
        else :
            cmd = '{kraken2} -db {kraken_database} --threads 8 --output - --report {assembly}.kraken {assembly}'.format(
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
        for l50, s in enumerate(seq) :
            acc[0], acc[1] = acc[0] + s[0], acc[1] + s[0]*s[1]
            if acc[0] * 2 >= n_base :
                n50 = s[0]
                break
        l50 += 1
        ave_depth = acc[1]/acc[0]
        n_low = 0
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
    
    reads = []
    for k, vs in zip(('pe', 'se'), (parameters['pe'], parameters['se'])) :
        for v in vs :
            if k == 'pe' :
                rnames = v.split(',')
                if len(rnames) > 0 :
                    assert len(rnames) == 2, 'Allows 2 reads per PE library. You specified {0}'.format(len(rnames))
                    reads.append(rnames)
            elif k == 'se' :
                rnames = v.split(',')
                if len(rnames) > 0 :
                    assert len(rnames) == 1, 'Allows one file per SE library. You specified {0}'.format(len(rnames))
                    reads.append(rnames)

    logger('Load in {0} read files from {1} libraries'.format(sum([ len(lib) for lib in reads ]), len(reads)))
    if not parameters['onlyEval'] :
        assembly = mainprocess().launch(reads)
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
    parser.add_argument('-p', '--prefix', help='prefix for the outputs. Default: EToKi_assemble', default='EToKi_assemble')
    parser.add_argument('-a', '--assembler', help='Assembler used for de novo assembly. \nDisabled if you specify a reference. \nDefault: spades for single colony isolates, megahit for metagenome', default='')
    parser.add_argument('-r', '--reference', help='Reference for read mapping. Specify this for reference mapping module. ', default=None)
    parser.add_argument('-k', '--kmers', help='relative lengths of kmers used in SPAdes. Default: 30,50,70,90', default='30,50,70,90')
    parser.add_argument('-m', '--mapper', help='aligner used for read mapping.\noptions are: miminap (default), bwa and bowtie2', default='minimap2')
    parser.add_argument('-d', '--max_diff', help='Maximum proportion of variations allowed for a aligned reads. \nDefault: 0.1 for single isolates, 0.05 for metagenome', type=float, default=-1)
    parser.add_argument('-i', '--ingroup', help='Additional references presenting intra-population genetic diversities. ', default='')
    parser.add_argument('-o', '--outgroup', help='Additional references presenting genetic diversities outside of the studied population. \nReads that are more similar to outgroups will be excluded from analysis. ', default='')
    
    parser.add_argument('-S', '--SNP', help='Exclusive set of SNPs. This will overwrite the polish process. \nRequired format:\n<cont_name> <site> <base_type>\n...', default=None)
    parser.add_argument('-c', '--cont_depth', help='Allowed range of read depth variations relative to average value.\nDefault: 0.2,2.5\nContigs with read depths outside of this range will be removed from the final assembly.', default='')
    
    parser.add_argument('--excluded', help='A name of the file that contains reads to be excluded from the analysis.', default='')
    parser.add_argument('--metagenome', help='Reads are from metagenomic samples', action='store_true', default=False)

    parser.add_argument('--noPolish', help='Do not do PILON polish.', action='store_true', default=False)
    parser.add_argument('--reassemble', help='Do local re-assembly in PILON', action='store_true', default=False)
    parser.add_argument('--onlySNP', help='Only modify substitutions during the PILON polish.', action='store_true', default=False)
    parser.add_argument('--noQuality', help='Do not estimate base qualities.', action='store_true', default=False)
    parser.add_argument('--onlyEval', help='Do not run assembly/mapping. Only evaluate assembly status.', action='store_true', default=False)
    parser.add_argument('--kraken', dest='runKraken', help='Run kmer based species predicton on contigs.', action='store_true', default=False)

    args = parser.parse_args(a)
    if args.cont_depth == '' :
        args.cont_depth = '0.2,2.5' if not args.metagenome else '0.0001,10000'
    if args.assembler == '' :
        args.assembler = 'spades' if not args.metagenome else 'megahit'
    if args.max_diff < 0 :
        args.max_diff = 0.1 if not args.metagenome else 0.05
    if args.metagenome :
        args.reassemble = False
        
    return args

if __name__ == '__main__' :
    assemble(sys.argv[1:])

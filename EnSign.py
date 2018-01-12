import os, sys, numpy as np
from subprocess import Popen, PIPE
from operator import itemgetter
from EnStal import externals, logger

class dualBlast(object) :
    def fastaLength(self, filename) :
        seq = {}
        with open(filename) as fin :
            for line in fin:
                if line[0] == '>' :
                    name = line[1:].strip().split()[0]
                    seq[name] = 0
                else :
                    seq[name] += len(line.strip())
        return seq
    def load_ublast(self, filename) :
        blasttab = []
        min_iden = parameters['min_iden']*100
        with open(filename) as fin :
            for line in fin :
                part = line.strip().split()
                if float(part[2]) > min_iden - 10. and int(part[7]) - int(part[6]) + 1 >= parameters['min_frag_len'] :
                    blasttab.append(part[:2] + [float(x) for x in part[2:14]])
        return blasttab
    def transform_blastp(self, blasttab, ql, tl) :
        for part in blasttab :
            part[0], f0 = part[0].rsplit('_', 1)
            part[1], f1 = part[1].rsplit('_', 1)
            frames = [int(f0), int(f1)]
            part[12], part[13] = ql[part[0]], tl[part[1]]
    
            if (part[7]-1)*3 + frames[0] %3 + 2 > part[12] or \
               (part[9]-1)*3 + frames[1] %3 + 2 > part[13] :
                part[7], part[9], part[3] = part[7]-1, part[9]-1, part[3]-1
                
            if frames[0] <= 3 :
                part[6:8] = [(int(part[6])-1)*3 + frames[0], (int(part[7])-1)*3 + frames[0] + 2]
            else :
                part[6:8] = [part[12] - ( (int(part[6])-1)*3 + frames[0] -4 ), part[12] - ( (int(part[7])-1)*3 + frames[0] -2 )]
            if frames[1] <= 3 :
                part[8:10] = [(int(part[8])-1)*3 + frames[1], (int(part[9])-1)*3 + frames[1] + 2]
            else :
                part[8:10] = [part[13] - ( (int(part[8])-1)*3 + frames[1] -4 ), part[13] - ( (int(part[9])-1)*3 + frames[1] -2 )]
            part[3], part[4] = 3 * int(part[3]), 3 * int(part[4])
        return blasttab
    
    def run_ublast(self, fna_target, faa_target, fna_query, faa_query=None) :
        tgt_len, qry_len = self.fastaLength(fna_target), self.fastaLength(fna_query)
        
        format_cmd = '{formatdb} -dbtype nucl -in {na_db}'.format(na_db=fna_target, na_input=fna_query, **parameters)
        Popen(format_cmd.split(), stderr=PIPE, stdout=PIPE).communicate()
        blast_cmd = '{blast} -db {na_db} -query {na_input} -out {na_db}.b6a -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen slen" -num_threads 6 -task blastn -evalue 1e-3 -dbsize 5000000 -reward 2 -penalty -2 -gapopen 6 -gapextend 2'.format(
            na_db=fna_target, na_input=fna_query,  **parameters
        )
        Popen(blast_cmd, stderr=PIPE, stdout=PIPE, shell=True).communicate()
        na_blast = self.load_ublast(fna_target + '.b6a')
        logger('Obtain {0} hits with BLASTn'.format(len(na_blast)))
        if faa_query :
            ublast_cmd = '{ublast} -threads 6 -db {aa_db} -ublast {aa_input} -evalue 1e-3 -accel 0.9 -maxhits 6 -userout {aa_db}.b6a -ka_dbsize 5000000 -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+raw+ql+tl'.format(
                aa_db=faa_target, aa_input=faa_query,  **parameters
            )
            Popen(ublast_cmd.split(), stderr=PIPE, stdout=PIPE).communicate()
            na_blast.extend( self.transform_blastp(self.load_ublast(faa_target + '.b6a'), qry_len, tgt_len ) )
            logger('Obtain {0} hits with BLASTn & uBLASTp'.format(len(na_blast)))            
        return na_blast
    
class blastParser(object) :
    def linear_merge(self, blasttab, min_iden, min_frag_prop, min_frag_len, max_dist, diag_diff, max_diff, **params) :
        for part in blasttab :
            if part[8] > part[9] :
                part[8], part[9] = -part[8], -part[9]

        blasttab.sort(key=itemgetter(0,1,6,8))
        
        syntenies = []
        for id, p1 in enumerate(blasttab) :
            for jd in xrange(id+1, len(blasttab))  :
                p2 = blasttab[jd]
                
                if p1[0] != p2[0] or p1[1] != p2[1] or p2[6] - p1[7] > max_dist : 
                    break
                elif (p1[8] < 0 and p2[8] > 0) or p2[8] < p1[8] + 15 or p2[9] < p1[9] + 15 or p2[7]<p1[7]+15 or p2[6] < p1[6]+15 or p2[8] - p1[9] > max_dist : 
                    continue
                m = max(p2[7], p1[7]) - min(p2[6], p1[6]) + 1
                n = max(p2[9], p1[9]) - min(p2[8], p1[8]) + 1
                if m < min_frag_len or m < min_frag_prop*p1[12] or max(m, n) - min(m, n) > max_diff or max(m,n) > diag_diff * min(m, n) : 
                    continue
                p1_len = min(p1[7], p2[6]-1) - p1[6] + 1
                p2_len = p2[7] - max(p2[6], p1[7]+1) + 1 if p2[7] > p1[7] else 0
                o_len = 0 if p2[6] > p1[7] else p1[7] - p2[6] + 1
                iden = (p1[2]*p1_len + p2[2]*p2_len +max(p1[2], p2[2])*o_len)/(p1_len+p2_len+o_len)
                if iden < min_iden : 
                    continue
                
                p1s, p2s = p1[11]/(p1[7]-p1[6]+1), p2[11]/(p2[7]-p2[6]+1)
                dist = max(p2[6]-p1[7]-1, p2[8]-p1[9]-1, 0)
                score = p1s*p1_len +p2s*p2_len+max(p1s, p2s)*o_len - dist
                if score > 0 :
                    syntenies.append([id, jd, iden, score])
        syn_score = {}
        for id, syn in enumerate(syntenies) :
            if syn[0] not in syn_score or syntenies[syn_score[syn[0]]][3] < syn[3] :
                syn_score[syn[0]] = id
            if syn[1] not in syn_score or syntenies[syn_score[syn[1]]][3] < syn[3] :
                syn_score[syn[1]] = id
        for id , syn in enumerate(syntenies) :
            if syn_score[ syn[0] ] == id and syn_score[ syn[1] ] == id :
                p1, p2 = blasttab[syn[0]][:], blasttab[syn[1]]
                part[7], part[9], part[2], part[11] = p2[7], p2[9], syn[2], syn[3]
                blasttab.append(part)
        for part in blasttab :
            if part[8] < 0 :
                part[8], part[9] = -part[8], -part[9]
        return blasttab
    def parse_ublast(self, hits, parameters) :
        for part in hits:
            allele = part[0].rsplit('_', 1)
            r_start, r_end = 1, part[12]
            direct = 1 if part[8] < part[9] else -1
            tailing = min(part[6]-r_start, part[7] - r_end)
    
            ## raw estimation of q_start and q_end
            if part[6] < r_start + 7 :
                while part[6] > 1 and part[8] > 1 and part[8] < part[13] :
                    part[6], part[8] = part[6]-1, part[8]-direct
            if part[7] > r_end - 7 :
                while part[7] < r_end and part[9] > 1 and part[9] < part[13] :
                    part[7], part[9] = part[7]+1, part[9]+direct
    
            q_start, q_end = part[8], part[9]
            
            ## ignore the hit if it is too small
            if direct*(q_end - q_start) + 1 < parameters['min_frag_prop'] * (r_end - r_start + 1) :
                part[0] = ''
            else :
                part.extend([r_start, r_end, q_start, q_end, float(part[11]), allele[0], allele[1], tailing])
                #              14       15      16      17        18            19          20        21
        hits = sorted([hit for hit in hits if hit[0] != ''], key = itemgetter(21, 18), reverse=True)
        alleles = {'__non_specific__':[]}
    
        for part in hits:
            ## check whether the allele has overlapped with other better hits
            if part[18] <= 0 : continue
            #direct = 1 if part[8] < part[9] else -1
            allele_full = '{0}_{1}'.format(*part[19:])
            r_start, r_end, q_start, q_end = part[14:18]
            overlap, to_move = -1, []

            if part[19] not in alleles :
                alleles[part[19]] = []
            else :
                s, e = sorted([q_start, q_end])
                for id, region in enumerate( alleles[part[19]] ):
                    if region[0] == '' : continue
                    if region[2] == part[1] and (min(e, region[4]) - max(s, region[3]) + 1) >= parameters['merging_prop']*(e-s+1) :
                        if region[1] < part[2] and (min(e, region[4]) - max(s, region[3]) + 1) >= parameters['merging_prop']*(region[4]-region[3]+1) :
                            if region[1] < parameters['min_iden']*100 or (part[21] >= 0 and part[2] - region[1] >= parameters['merging_error']*100):
                                to_move.append(id)
                                continue
                            else :
                                region[1] = part[2]
                        overlap = id
    
                    elif region[1] - part[2] >= parameters['merging_error']*100 :
                        overlap = 999999999
                        break
                    elif part[2] - region[1] >= parameters['merging_error'] * 100 :
                        to_move.append(id)
                if overlap < 999999999 :
                    for id in reversed(to_move) :
                        alleles['__non_specific__'].append( alleles[part[19]].pop(id)[:-1] + [''] )
            
            

            if overlap < 0 :
                # insert a new region if there is no overlap
                if q_start < q_end :
                    alleles[part[19]].append([part[19], float(part[2]), part[1], q_start, q_end, '+', r_start - part[6], part[7]-r_end, '', allele_full])
                else :
                    alleles[part[19]].append([part[19], float(part[2]), part[1], q_end, q_start, '-', r_start - part[6], part[7]-r_end, '', allele_full])
            elif overlap == 999999999 :
                if q_start < q_end :
                    alleles['__non_specific__'].append([part[19], float(part[2]), part[1], q_start, q_end, '+', r_start - part[6], part[7]-r_end, '', ''])#part[0]])
                else :
                    alleles['__non_specific__'].append([part[19], float(part[2]), part[1], q_end, q_start, '-', r_start - part[6], part[7]-r_end, '', ''])#part[0]])
        return alleles
        
    def inter_loci_overlap(self, alleles, parameters) :
        regions = [reg for region in alleles.values() for reg in region]
        # sort with contig name and start points
        regions.sort(key=itemgetter(2,3))
        
        for id, regi in enumerate(regions) :
            if regi[0] == '' : continue
            todel, deleted = [], 0
            for jd in xrange(id+1, len(regions)) :
                regj = regions[jd]
                if regj[0] == '' or regi[0] == regj[0]: continue
                if regi[2] != regj[2] or regj[3] > regi[4] :
                    break
                overlap = min(regi[4], regj[4]) - regj[3] + 1
                if regi[-1] != '' and float(overlap) >= parameters['merging_prop'] * (regi[4]-regi[3]+1) :
                    delta = regi[1] - regj[1]
                            
                    if delta > parameters['merging_error']*100 :
                        todel.append(jd)
                    elif delta < - parameters['merging_error']*100 or (delta < 0 and regj[-1] == '') :
                        deleted = 1
                        break
                    else :
                        regi[8] += '{{M:{0}:{1}}}'.format(*regj)
                if regj[-1] != '' and float(overlap) >= parameters['merging_prop'] * (regj[4]-regj[3]+1) :
                    delta = regj[1] - regi[1]
                    if delta > parameters['merging_error']*100 :
                        deleted = 1
                        break
                    elif delta < - parameters['merging_error']*100 or (delta < 0 and regi[-1] == '') :
                        todel.append(jd)
                    else :
                        regj[8] += '{{M:{0}:{1}}}'.format(*regi)
            if deleted == 0 :
                for jd in todel:
                    regions[jd][0] = ''
            else :
                regi[0] = ''
                
        return [{'locus':reg[0], 'identity':reg[1], 'coordinates':[reg[2], int(reg[3]), int(reg[4]), reg[5]], 'flanking':reg[6:8], 'status':reg[8], 'accepted':(0 if reg[8] == '' else 128)} for reg in regions if reg[0] != '' and reg[-1] != '']
    
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}    
    def get_seq(self, fastq, header, start, end, direction) :
        fasta = fastq[header][0]
        if direction == '+' :
            seq = fasta[(start-1):end].upper()
        else :
            seq = ''.join([self.complement.get(x.upper(), 'N') for x in reversed(fasta[(start-1):end])])
        return seq
    def get_qual(self, fastq, header, start, end, direction='+') :
        if fastq == None:
            return 20
        s = max(1, min(start, end))
        e = min(len(fastq[header][0]), max(start, end))
        if 'N' in fastq[header][0][(start-1):end].upper() :
            return 0
        else :
            try:
                return min(fastq[header][1][(s-1):e])
            except :
                return 0
    
    def form_alleles(self, regions, fastq, genome_id, accepted, argument) :
        alleles = {}
        regions.sort(key=lambda x:x['identity'], reverse=True)
        regions.sort(key=lambda x:min(x['flanking'] + [0]), reverse=True)
        for region in regions:
            region['seq'] = self.get_seq(fastq, *region['coordinates'])
            if min(region['flanking']) >= 0 and 'N' not in region['seq'] :  ## add proportional check
                region['accepted'] = region['accepted'] | 1
            else :
                region['status'] = '{Fragmented}'
                region['accepted'] = region['accepted'] | 64
                region['allele_id'] = -1

            if region['locus'] in alleles :
                if region['accepted'] & 64 > 0 :
                    if alleles[ region['locus'] ]['accepted'] & 64 > 0 :
                        if 'secondary' not in alleles[ region['locus'] ] :
                            alleles[ region['locus'] ]['secondary'] = []
                        alleles[ region['locus'] ]['secondary'].append( dict(coordinates =region['coordinates'], seq=region['seq'], identity=region['identity']) )
                elif alleles[ region['locus'] ] ['accepted'] & 32 == 0 :
                    alleles[ region['locus'] ] ['status'] += '{Duplicated}'
                    alleles[ region['locus'] ] ['accepted'] = (alleles[ region['locus'] ] ['accepted'] | 32) & (~1)
                    alleles[ region['locus'] ] ['allele_id'] = -1
                    if 'secondary' not in alleles[ region['locus'] ] :
                        alleles[ region['locus'] ]['secondary'] = []
                    alleles[ region['locus'] ]['secondary'].append( dict(coordinates =region['coordinates'], seq=region['seq'], identity=region['identity']) )
            else :
                if accepted == 0 or self.get_qual(fastq, *region['coordinates']) < 10:
                    region['accepted'] = (region['accepted'] & (~ 1)) | 2
                region['reference'] = {'source':'enSign', 'barcode': genome_id}
                alleles[region['locus']] = region
        for locus, allele in alleles.items() :
            if allele['accepted'] & 65 == 64 :
                allele_len = allele['coordinates'][2] - allele['coordinates'][1] + 1
                for ale in allele.get('secondary', {}) :
                    allele_len += ale['coordinates'][2] - ale['coordinates'][1] + 1
                if allele_len < argument['min_frag_prop'] :
                    alleles.pop(locus)
            if 'identity' in allele and allele['identity'] < argument['min_iden']* 100 :
                allele['allele_id'] = -1
                allele['accepted'] = (allele['accepted'] & (~1)) | 256
                allele['status'] += '{Low identities}'
                if allele['accepted'] & 224 > 0 :
                    alleles.pop(locus, None)
        return alleles

class seqOperation(object) :
    def readSequence(self, fname) :
        seq = {}
        with open(fname) as fin :
            header = fin.read(1)
            fin.seek(0, 0)
            if header == '@' :
                for id, line in enumerate(fin) :
                    if id % 4 == 0 :
                        name = line[1:].strip().split()[0]
                        seq[name]= [None, None]
                    elif id % 4 == 1 :
                        seq[name][0] = line.strip()
                    elif id % 4 == 3 :
                        seq[name][1] = [ord(q)-33 for q in line.strip()]
            else :
                for line in fin :
                    if line.startswith('>') :
                        name = line[1:].strip().split()[0]
                        seq[name] = [[], None]
                    else :
                        seq[name][0].extend( line.strip().split() )
                for n, s in seq.iteritems() :
                    s[0] = ''.join(s[0])
                    s[1] = [0 for ss in s[0]]
        return seq
    
    def transeq(self, seq, frame=1, transl_table=11) :
        gtable = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
                  "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",    "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
                  "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                  "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                  "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                  "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
                  "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                  "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
        
        complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
        def rc(seq) :
            return ''.join([complement.get(s, 'N') for s in reversed(seq.upper())])
        
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
    
    def write_refsets(self, reference) :
        ref_aa = '{0}.refset.aa'.format(parameters['unique_key'])
        refseq = self.readSequence(reference)
        refamino = self.transeq({n:s[0] for n,s in refseq.iteritems() }, 1)
        with open(ref_aa, 'w') as fout :
            for n, s in refamino.iteritems() :
                if s[:-1].find('X') == -1 :
                    fout.write('>{0}\n{1}\n'.format(n, s))
        return ref_aa
    
    def write_query(self, query) :
        fna, faa = '{0}.query.na'.format(parameters['unique_key']), '{0}.query.aa'.format(parameters['unique_key'])
        qryseq = self.readSequence(query)
        qryamino = self.transeq({n:s[0] for n,s in qryseq.iteritems()}, frame=7)
        with open(fna, 'w') as fout:
            for n, s in qryseq.iteritems() :
                fout.write('>{0}\n{1}\n'.format(n, s[0]))
        with open(faa, 'w') as fout :
            for n, s in qryamino.iteritems() :
                fout.write('>{0}\n{1}\n'.format(n, s))
        return qryseq, fna, faa

def nomenclature(query, reference, ref_aa='', **params) :
    # write query
    logger('EnSign starts')
    sequence, qry_fna, qry_faa = seqOperation().write_query(query)
    logger('Read in {0} bases as query'.format(sum([len(s[0]) for s in sequence.itervalues()])))
    # write refset
    if not os.path.isfile(str(ref_aa)) :
        ref_aa = seqOperation().write_refsets(reference)
        logger('Prepare translated references')
    # do comparison
    blasttab = dualBlast().run_ublast(fna_target=qry_fna, faa_target=qry_faa, fna_query=reference, faa_query=ref_aa)
    # filter
    blasttab_parser = blastParser()
    blasttab = blasttab_parser.linear_merge(blasttab, **parameters)
    logger('Merge closely located hits. {0} hits'.format(len(blasttab)))
    loci = blasttab_parser.parse_ublast(blasttab, parameters)
    logger('Identify homologous groups. {0} groups'.format(len([1 for lc in loci if lc != '__non_specific__'])))
    regions = blasttab_parser.inter_loci_overlap(loci, parameters)
    logger('Resolve potential paralogs. {0} regions'.format(len(regions)))
    
    # submission
    alleles = blasttab_parser.form_alleles(regions, sequence, parameters['unique_key'], parameters['high_quality'], parameters)
    logger('Generate allelic sequences. {0} remains'.format(len(alleles)))
    #results = blasttab_parser.typing(alleles, parameters, dbname, scheme, submission=submission)
    return alleles

parameters = dict(
    unique_key = 'entype', 
    ref_na='',
    ref_aa = '', 
    query = '', 
    high_quality = True,
    
    min_iden = 0.65, 
    min_frag_prop = .6, 
    min_frag_len = 50, 
    
    merging_prop = 0.5,
    merging_error = 0.05, 
    
    max_dist=300, 
    diag_diff=1.2, 
    max_diff=200, 
)
parameters.update(externals)

def enSign() :
    parameters.update(dict( arg.split('=', 1) for arg in sys.argv[1:] ))
    return nomenclature(**parameters)


if __name__ == '__main__' :
    import json
    print json.dumps(enSign(), indent=2, sort_keys=True)
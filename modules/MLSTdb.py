import os, sys, shutil, pandas as pd
try :
    from configure import externals, logger, uopen, xrange, StringIO, get_md5, readFasta
    from uberBlast import uberBlast
    from clust import clust
except :
    from .configure import externals, logger, uopen, xrange, StringIO, get_md5, readFasta
    from .uberBlast import uberBlast
    from .clust import clust
import subprocess, tempfile, time

mmseqs = externals['mmseqs']
minimap2 = externals['minimap2']

def minimapFilter(sourceFna, targetFna, targetFiltFna, max_iden, min_iden, coverage, paralog, relaxEnd, orderedLoci) :
    p = subprocess.Popen('{0} -ct8 -k13 -w5 -A2 -B4 -O8,16 -E2,1 -r50 -p.2 -N500 -f2000,10000 -n1 -m19 -s40 -g200 -2K10m --heap-sort=yes --secondary=yes {1} {2}'.format(minimap2, sourceFna, targetFna).split(), stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines=True)
    tooClose, goodCandidates, crossLoci = {}, {}, {}
    for line in p.stdout :
        part = line.strip().split('\t')
        q, r = part[0], part[5]
        qLoc, rLoc = q.rsplit('_', 1)[0], r.rsplit('_', 1)[0]
        if q == r :
            goodCandidates[r] = 1.
            continue
        elif rLoc == qLoc :
            if part[4] == '-' :
                continue
            tl, ts, te = [int(p) for p in part[1:4]]
            rl, rs, re, ri = [int(p) for p in part[6:10]]
            if relaxEnd or (ts != rs or tl-te != rl - re) :
                continue
            iden, cov = float(ri)/(re-rs), float(re-rs)/rl
            if cov > coverage :
                if iden >= max_iden :
                    tooClose[part[0]] = 1
                elif iden >= min_iden :
                    goodCandidates[part[0]] = max(goodCandidates.get(part[0], 0), iden)
        elif orderedLoci[qLoc] > orderedLoci[rLoc] :
            rl, rs, re, ri = [int(p) for p in part[6:10]]
            iden, cov = float(ri)/(re-rs), float(re-rs)/rl
            if cov > coverage and iden >= min_iden and crossLoci.get(part[0], 0) < iden :
                crossLoci[part[0]] = iden
    with open(targetFna) as fin, open(targetFiltFna+'.fas', 'w') as fout :
        writable = False
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip().split()[0]
                writable = False if (name in tooClose or crossLoci.get(name, 0) > 1. - paralog) else True
            if writable :
                fout.write(line)
    return targetFiltFna, goodCandidates, crossLoci

def buildReference(alleles, references, max_iden=0.9,  min_iden=0.6, coverage=0.7, paralog=0.1, relaxEnd=False) :
    orderedLoci = { t['fieldname']:i for i, t in reversed(list(enumerate(references))) }
    with tempfile.TemporaryDirectory(prefix='NS_', dir='.') as dirPath :
        sourceFna = os.path.join(dirPath, 'sourceFna')
        clsFna = os.path.join(dirPath, 'clsFna')
        targetFna = os.path.join(dirPath, 'targetFna')
        with open(sourceFna, 'w') as fout :
            fout.write('\n'.join(['>{fieldname}_{value_id}\n{value}'.format(**s) for s in alleles]))
        with open(targetFna, 'w') as fout :
            fout.write('\n'.join(['>{fieldname}_{value_id}\n{value}'.format(**t) for t in references]))

        if 1-paralog > max_iden :
            exampler, cluster = clust('-i {0} -p {1} -d {2} -c 1 -t 8'.format( \
                sourceFna, clsFna, 1-paralog).split())
            tooClose, goodCandidates, crossSites = {}, {}, {}
            with open(cluster) as fin:
                for line in fin:
                    part = line.strip().split()
                    locus = [p.rsplit('_', 1)[0] for p in part]
                    if locus[0] != locus[1]:
                        crossSites[part[0]] = crossSites.get(part[0], []) + [locus[1]]
                        crossSites[part[1]] = crossSites.get(part[1], []) + [locus[0]]
            exampler, cluster = clust('-i {0} -p {1} -d {2} -c 1 -t 8'.format( \
                sourceFna, clsFna, max_iden).split())
        else :
            # get cluster
            exampler, cluster = clust('-i {0} -p {1} -d {2} -c 1 -t 8'.format(\
                sourceFna, clsFna, max_iden).split())
            tooClose, goodCandidates, crossSites = {}, {}, {}
            with open(cluster) as fin :
                for line in fin :
                    part = line.strip().split()
                    locus = [ p.rsplit('_', 1)[0] for p in part ]
                    if locus[0] != locus[1] :
                        crossSites[part[0]] = crossSites.get(part[0], [])+[locus[1]]
                        crossSites[part[1]] = crossSites.get(part[1], [])+[locus[0]]
                
        # compare with references
        blastab = uberBlast('-r {0} -q {1} -f --blastn --diamondSELF --min_id {2} --min_ratio {3} -t 8 -p -s 1 -e 0,3'.format(\
            targetFna, exampler, min_iden, coverage ).split())
        #blastab = blastab[blastab.T[0] != blastab.T[1]]

    for tab in blastab :
        locus = [ p.rsplit('_', 1)[0] for p in tab[:2] ]
        c = (tab[7]-tab[6]+1)/tab[12]
        e = max(abs(tab[8] - tab[6]), abs(tab[12]-tab[7] - (tab[13]-tab[9])))
        if c >= coverage and tab[2] >= min_iden :
            if locus[0] != locus[1] and tab[2] >= 1-paralog :
                crossSites[tab[0]] = crossSites.get(tab[0], [])+[locus[1]]
                crossSites[tab[1]] = crossSites.get(tab[0], [])+[locus[0]]
            elif e <= 0 or relaxEnd :
                if tab[2] >= max_iden and tab[0] != tab[1] :
                    tooClose[tab[0]] = 1
                else :
                    goodCandidates[tab[0]] = tab[2]
    paralogous_loci = {}
    for ref in references :
        key = '{0}_{1}'.format(ref['fieldname'], ref['value_id'])
        if key in crossSites:
            conflicts = crossSites[key]
            for cfl in conflicts :
                if orderedLoci[ref['fieldname']] > orderedLoci[cfl] and cfl not in paralogous_loci :
                    paralogous_loci[ref['fieldname']]=1
                    break
    refsets = []
    for allele in alleles :
        if allele['fieldname'] in paralogous_loci :
            allele['fieldname'] = ''
        else :
            key = '{0}_{1}'.format(allele['fieldname'], allele['value_id'])
            if key in crossSites :
                allele['fieldname'] = ''
            elif key in goodCandidates and key not in tooClose :
                refsets.append('>{fieldname}_{value_id}\n{value}'.format(**allele))
    alleles = [ '>{fieldname}_{value_id}\n{value}'.format(**allele) for allele in alleles if allele['fieldname'] != '' ]
    logger('removed {0} paralogous sites.'.format(len(paralogous_loci)))
    logger('obtained {0} alleles and {1} references alleles'.format(len(alleles), len(refsets)))
    return '\n'.join(alleles), '\n'.join(refsets)


# Read a fasta into an array whose elements are a dict with keys:
#   fieldname: seqname (STRING)
#   value_id:  allele (STRING)
#   value:     sequence (STRING)
def readFasta(fastaText) :
    sequence = []
    for line in fastaText :
        # Get the defline
        if line.startswith('>') :
            # The name of this sequence starts in position 1 (after the >)
            # but then we remove anything after whitespace to just keep the unique seqid
            name = line[1:].strip().split()[0]
            # part evaluates to the whole seqid and then the part after the last underscore
            part = name.rsplit('_', 1)
            if len(part) > 1 :
                # If there is the part after the last underscore, save that as
                # the allele number into 'value_id'
                # and the seqid minus the allele into 'fieldname'
                sequence.append({'fieldname':part[0], 'value_id':part[1], 'value':[]})
            else :
                # If there is no allele, just save the seqid into fieldname
                sequence.append({'fieldname':name, 'value_id':None, 'value':[]})
        elif len(line) > 0 :
            # add the sequence onto the latest element of the sequence array
            # into 'value'
            sequence[-1]['value'].extend(line.strip().split())
    for s in sequence :
        # Turn 'value' into a string
        s['value'] = ''.join(s['value'])
    return sequence


def MLSTdb(args) :
    params = getParams(args)
    database, refset, alleleFasta, refstrain, max_iden, min_iden, coverage, paralog, relaxEnd = \
        params['database'], params['refset'], params['alleleFasta'], params['refstrain'], params['max_iden'], params['min_iden'], params['coverage'], params['paralog'], params['relaxEnd']
    # Read the fasta from either a path or a string and I guess 'alleleFasta' could be either
    if os.path.isfile(alleleFasta) :
        alleles = readFasta(uopen(alleleFasta))
    else :
        alleles = readFasta(StringIO(alleleFasta))
    # filter the alleles array for those with an allele (value_id.isdigit())
    # and whose allele is greater than 0 and whose locus name does not have '/'
    alleles = [allele for allele in alleles \
                   if allele['value_id'].isdigit() and int(allele['value_id']) > 0 and allele['fieldname'].find('/') < 0]
    refAlleles = ''
    if refset is not None :
        if refstrain :
            if os.path.isfile(refstrain) :
                references = readFasta(uopen(refstrain))
            else :
                references = readFasta(StringIO(refstrain))
        else :
            loci, references = {}, []
            for allele in alleles :
                if allele['fieldname'] not in loci :
                    loci[allele['fieldname']] = 1
                    references.append(allele)
        
        allele_text, refAlleles = buildReference(alleles, references, max_iden, min_iden, coverage, paralog, relaxEnd)
        if refset :
            with open(str(refset), 'w') as fout :
                fout.write(refAlleles + '\n')
        logger('A file of reference alleles has been generated:  {0}'.format(refset))

    # Create a file for lookup table of all alleles
    if database :
        # Conversion is an array of arrays:
        #   0 - md5sums of allele sequences
        #   1 - array: [locus, allele]
        conversion = [[], []]
        for allele in alleles :
            # If this fasta entry has a locus name...
            if allele['fieldname'] :
                conversion[0].append(get_md5(allele['value']))
                conversion[1].append([allele['fieldname'], int(allele['value_id'])])

        conversion = pd.DataFrame(conversion[1], index=conversion[0])
        conversion.to_csv(database, header=False)
        logger('A lookup table of all alleles has been generated:  {0}'.format(database))
    return allele_text, refAlleles

def getParams(args) :
    import argparse
    parser = argparse.ArgumentParser(description='MLSTdb. Create reference sets of alleles for nomenclature. ')
    parser.add_argument('-i', '--input', dest='alleleFasta', help='[REQUIRED] A single file contains all known alleles in a MLST scheme. ', required=True)
    parser.add_argument('-r', '--refset',                    help='[DEFAULT: No ref allele] Output - Reference alleles used for MLSType. ', default=None)
    parser.add_argument('-d', '--database',                  help='[DEFAULT: No allele DB] Output - A lookup table of all alleles. ', default=None)
    parser.add_argument('-s', '--refstrain',                 help='[DEFAULT: None] A single file contains alleles from the reference genome. ', default=None)
    parser.add_argument('-x', '--max_iden',                  help='[DEFAULT: 0.9 ] Maximum identities between resulting refAlleles. ', type=float, default=0.9)
    parser.add_argument('-m', '--min_iden',                  help='[DEFAULT: 0.6 ] Minimum identities between refstrain and resulting refAlleles. ', type=float, default=0.6)
    parser.add_argument('-p', '--paralog',                   help='[DEFAULT: 0.2 ] Minimum differences between difference loci. ', type=float, default=0.2)
    parser.add_argument('-c', '--coverage',                  help='[DEFAULT: 0.7 ] Proportion of aligned regions between alleles. ', type=float, default=0.7)
    parser.add_argument('-e', '--relaxEnd',                  help='[DEFAULT: False ] Allow changed ends (for pubmlst). ', action='store_true', default=False)

    return parser.parse_args(args).__dict__
    
if __name__ == '__main__' :
    MLSTdb(sys.argv[1:])

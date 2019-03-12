import sys, subprocess, os, shlex, json
try:
    from .configure import externals, logger, readFasta, xrange
except :
    from configure import externals, logger, readFasta, xrange

def parse_bsn(save) :
    region = []
    value = 0.
    for bsn in save :
        s, e, rate = float(bsn[6]), float(bsn[7]), (float(bsn[2])/100.)**4
        for reg in region :
            if e >= reg[0] and s<= reg[1]:
                if s <= reg[0] and e >= reg[1] :
                    s = e+1
                if s >= reg[0] :
                    s = reg[1] + 1
                if e <= reg[1] :
                    e = reg[0] - 1
                if s > e :
                    break
        if e >= s :
            region.append([s, e])
            value += (e-s+1)*rate
    return value / float(save[0][-1])

def run_prediction(prefix, assembly, db) :
    subprocess.Popen('{0} -in {1} -dbtype nucl -out {2}'.format(externals['makeblastdb'], assembly, prefix).split(), stdout=subprocess.PIPE).communicate()
    antigen_genes = {'H':{}, 'O':{}}
    with open(db) as fin :
        for line in fin :
            if line[0] == '>' :
                antigen, gene = line[1:].split('__')[:2]
                category = antigen[0] if antigen[0] in antigen_genes else 'O'
                gene = gene.upper()
                if antigen not in antigen_genes[category] :
                    antigen_genes[category][antigen] = {gene:-1.}
                else :
                    antigen_genes[category][antigen][gene] = -1.
    
    subprocess.Popen(shlex.split('{0} -task blastn -db {2} -query {1} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score qlen" -reward 1 -penalty -3 -out {2}.bsn'.format(externals['blastn'], db, prefix))).communicate()
    
    save = []
    with open('{0}.bsn'.format(prefix)) as fin :
        for line in fin :
            part = line.strip().split()
            if len(save) > 0 and save[0][0] != part[0] :
                antigen, gene = save[0][0].split('__')[:2]
                category = antigen[0] if antigen[0] in antigen_genes else 'O'
                gene = gene.upper()
                score = parse_bsn(save)
                if antigen_genes[category][antigen].get(gene, 0) < score :
                    antigen_genes[category][antigen][gene] = score
                save = []
            if int(part[3]) >= 50 and float(part[2]) >= 80 :
                save.append(part)
    for fname in {'{0}.nhr', '{0}.nin', '{0}.nsq', '{0}.bsn'} :
        try :
            os.unlink(fname.format(prefix))
        except :
            pass
    if len(save) > 0 :
        antigen, gene = save[0][0].split('__')[:2]
        category = antigen[0] if antigen[0] in antigen_genes else 'O'
        gene = gene.upper()
        score = parse_bsn(save)
        if antigen_genes[category][antigen].get(gene, 0) < score :
            antigen_genes[category][antigen][gene] = score
    prediction = {'query':assembly, 'prefix':prefix}
    for category, antigens in antigen_genes.items() :
        scores = sorted([[antigen, sum(genes.values())/len(genes) if antigen not in {'H3', 'H17', 'H35', 'H36', 'H44', 'H47', 'H53', 'H54', 'H55'} else sum(genes.values())*1.1/len(genes) ] for antigen, genes in antigens.items()], key=lambda x:x[1], reverse=True)
        antigen, score = scores[0]

        if score * 2. < 1. :
            prediction[category] = '-'
        elif score < .8 :
            prediction[category] = 'uncertain'
        elif  (len(scores) > 1 and scores[1][1] >= 0.99 * scores[0][1]) :
            prediction[category] = ' or '.join(sorted([t for t, s in scores if s >= 0.99*scores[0][1]]))
        else :
            prediction[category] = antigen
    return prediction

def EBEis(args):
    import argparse
    parser = argparse.ArgumentParser(description='''EnteroBase Escherichia in silico serotyping''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-q', '--query', help='file name for the queried assembly in multi-FASTA format.', required=True)
    parser.add_argument('-t', '--taxon', help='Taxon database to compare with. Only support Escherichia', default='Escherichia')
    parser.add_argument('-p', '--prefix', help='prefix for intermediate files.', default='EBEis')
    args = parser.parse_args(args)
    
    db = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'EBEis.{0}.fas'.format(args.taxon))
    prediction = run_prediction(args.prefix, args.query, db)
    sys.stdout.write(json.dumps(prediction, sort_keys=True) + '\n')

if __name__ == '__main__' :
    EBEis(sys.argv[1:])

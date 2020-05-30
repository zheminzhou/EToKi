import numpy as np, pandas as pd, click, re
try:
    from configure import uopen, xrange
except:
    from .configure import uopen, xrange

@click.command()
@click.argument('mgs', nargs=-1)
def main(mgs) :
    for mg in mgs :
        res = []
        with uopen(mg) as fin :
            for line in fin :
                logp = re.findall('logp:\t([-eE\d\.]+)', line)
                if len(logp) :
                    logp=float(logp[0])
                    res.append([logp])
                else :
                    genotype = re.findall('Genotype (\d+):\tMean proportion:\t([eE\d\.]+)\tCI95%:\t(\[ [eE\d\.]+ - [eE\d\.]+ \])', line)
                    if len(genotype) :
                        res[-1].append([genotype[0][1], genotype[0][2], '', '', ''])
                    elif len(res) and len(res[-1])>1 and res[-1][-1][-1] == '' :
                        part = line.strip().split('\t')
                        res[-1][-1][2:] = [part[0], part[1], part[3]+' '+part[5]]
        try :
            res = max(res)
            res[1:] = sorted(res[1:], key=lambda x:-float(x[0]))
            for i in xrange(1, len(res)) :
                r = res[i]
                print('{0}\t{1}\t{2}'.format(mg, i, '\t'.join(r)))
        except :
            pass

if __name__ == '__main__' :
    main()
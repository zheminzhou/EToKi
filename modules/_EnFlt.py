#!/usr/bin/env python

import sys, re, pandas as pd

cutoff = float(sys.argv[1])
try :
    outScore = {} if len(sys.argv) <= 2 else dict(pd.read_csv(sys.argv[2], sep='\t', header=None).values.tolist())
except :
    outScore = {}

pairs = {}

sq_name = {}

for line in sys.stdin :
    if line.startswith('@') :
        if line.startswith('@SQ') :
            if line.strip() not in sq_name :
                sq_name[line.strip()] = 1
            else :
                continue
        sys.stdout.write(line)
    else :
        if cutoff < 0 :
            sys.stdout.write(line)
            sys.stdout.flush()
        part = line.strip().split('\t', 10)
        flag = int(part[1])
        if flag & 4 > 0 :
            if (flag & 9) == 1 :
                sys.stdout.write(line)
                sys.stdout.flush()
            continue
        try:
            score = int(re.findall('AS:i:(\d+)', line)[0])
            if part[0] in outScore :
                if score < outScore[part[0]] :
                    if flag < 256 :
                        part[5] = '{0}S'.format(len(part[9]))
                        sys.stdout.write('\t'.join(part) + '\n')
                    continue
                else :
                    del outScore[part[0]]
        except :
            if flag < 256:
                part[5] = '{0}S'.format(len(part[9]))
                sys.stdout.write('\t'.join(part) + '\n')
            continue
        if 'H' in part[5] :
            gap = [int(g) for g in re.findall(r'(\d+)H', part[5])]
            aln = len(part[9])
        else :
            gap = [int(g) for g in re.findall(r'(\d+)S', part[5])]
            aln = len(part[9]) - sum(gap)

        if len(gap) > 1 and min(gap) >= 8 :
            if flag < 256:
                part[5] = '{0}S'.format(len(part[9]))
                sys.stdout.write('\t'.join(part) + '\n')
            continue
        if sum(gap) >= 1.5 * aln :
            if flag < 256:
                part[5] = '{0}S'.format(len(part[9]))
                sys.stdout.write('\t'.join(part) + '\n')
            continue
        dist = (aln*2. - score)/6.
        if dist >= aln*cutoff :
            if flag < 256:
                part[5] = '{0}S'.format(len(part[9]))
                sys.stdout.write('\t'.join(part) + '\n')
            continue
        sys.stdout.write(line)
        sys.stdout.flush()

#!/usr/bin/env python

import sys, re, pandas as pd

cutoff = float(sys.argv[1])
try :
    outScore = {} if len(sys.argv) <= 2 else dict(pd.read_csv(sys.argv[2], sep='\t', header=None).values.tolist())
except :
    outScore = {}
        
for line in sys.stdin :
    if line.startswith('@') :
        sys.stdout.write(line)
    else :
        part = line.strip().split('\t', 10)
        flag = int(part[1])
        if flag & 4 > 0 :
            continue
        try:
            score = int(re.findall('AS:i:(\d+)', line)[0])
            if part[0] in outScore :
                if score < outScore[part[0]] :
                    continue
                else :
                    del outScore[part[0]]
        except :
            continue
        if flag < 1024 :
            gap = [int(g) for g in re.findall(r'(\d+)S', part[5])]
            aln = len(part[9]) - sum(gap)
        else :
            gap = [int(g) for g in re.findall(r'(\d+)H', part[5])]
            aln = len(part[9])
        if len(gap) > 1 and min(gap) >= 8 :
            continue
        if sum(gap) >= 1.5 * aln :
            continue
        dist = (aln*2. - score)/6.
        if dist >= aln*cutoff :
            continue
        sys.stdout.write(line)
        sys.stdout.flush()

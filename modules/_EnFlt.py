#!/usr/bin/env python

import sys, re

cutoff = float(sys.argv[1])

for line in sys.stdin :
    if line.startswith('@') :
        sys.stdout.write(line)
    else :
        part = line.strip().split('\t')
        flag = int(part[1])
        if flag & 4 > 0 :
            continue
        if flag < 1024 :
            gap = [int(g) for g in re.findall(r'(\d+)S', part[5])]
            aln = len(part[10]) - sum(gap)
        else :
            gap = [int(g) for g in re.findall(r'(\d+)H', part[5])]
            aln = len(part[10])
        if len(gap) > 1 and min(gap) >= 10 :
            continue
        if sum(gap) >= 4 * aln :
            continue
        try:
            score = int(re.findall('AS:i:(\d+)', line)[0])
            dist = (aln*2. - score)/8.
        except :
            continue
        if dist >= aln*cutoff :
            continue
        sys.stdout.write(line)
        sys.stdout.flush()

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
        dist = (aln*2. - float(part[11][5:]))/8.
        if dist >= aln*cutoff :
            continue
        sys.stdout.write(line)
        sys.stdout.flush()

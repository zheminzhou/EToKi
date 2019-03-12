import os, sys
try :
    xrange(3)
except :
    xrange = range

if __name__ == '__main__' :
    match_identity, match_frag_len, match_frag_prop = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    for line in sys.stdin :
        part = line.strip().split()
        if float(part[2]) < 100.*match_identity or ( int(part[7]) - int(part[6]) + 1 < match_frag_len and int(part[7]) - int(part[6]) + 1 < match_frag_len * float(part[12]) ) :
            continue
        
        part[14] = part[15] if part[14].find('-') < 0 else \
            ''.join([ y if x.find('-')<0 else ''.join([yy for xx, yy in zip(x, y) if xx != '-' ]) for i in xrange(0, len(part[14]), 100) for x, y in [[part[14][i:i+100], part[15][i:i+100]]] ])
        
        sys.stdout.write('\t'.join(part[:15])+'\n')
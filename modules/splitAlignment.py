import sys, numpy as np, os, glob, re, argparse

try :
    from configure import externals, uopen, asc2int, logger, readFasta
except :
    from .configure import externals, uopen, asc2int, logger, readFasta


    
def add_args(a) :
    parser = argparse.ArgumentParser(description='split an aligment into even fragments', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--alignment', '-a', help='[REQUIRED; INPUT] the alignment in fasta format', required=True)
    parser.add_argument('--size', '-s', help='[DEFAULT: 1000; PARAMETER] the size of fragments', default=1000, type=int)
    parser.add_argument('--dir', '-d', help='[REQUIRED; OUTPUT] folder for all outputs.', required=True)
    args = parser.parse_args(a)
    return args
    
def splitAlignment(args) :
    args = add_args(args)
    seqs = readFasta(args.alignment)
    if not os.path.isdir(args.dir) :
        os.makedirs(args.dir)
    prefix = os.path.basename(args.dir)
    seq_size = len(list(seqs.values())[0])
    for i in np.arange(0, seq_size, args.size) :
        subseqs = []
        for n, s in seqs.items() :
            ss = s[i:(i+args.size)]
            if len(ss.replace('-', '')) > args.size/2 :
                subseqs.append([n, ss])
        if len(subseqs) >= 4 :
            with open(os.path.join(args.dir, '{0}.{1}.fas'.format(prefix, i+1)), 'wt') as fout :
                for n, ss in sorted(subseqs) :
                    fout.write('>{0}\n{1}\n'.format(n, ss))

if __name__ == '__main__' :
    splitAlignment(sys.argv[1:])

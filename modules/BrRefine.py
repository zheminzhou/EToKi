import numpy as np, pandas as pd, sys, math, argparse, gzip
import phylo 
from ete3 import Tree

def parse_arg(a) :
    parser = argparse.ArgumentParser(description='Generate a matrix of only vertically inherited SNPs. ', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--prefix', '-p', help='prefix for the output', required=True)
    parser.add_argument('--tree', '-t', help='Labeled tree', required=True)
    parser.add_argument('--rec', '-r', help='Recombinant sketches', required=True)
    args = parser.parse_args(a)
    return args

def read_RecHMM(fname) :
    br = {}
    with open(fname) as fin :
        for line in fin :
            part = line.strip().split('\t')
            if part[0] == 'Branch' :
                mut = float(part[2][2:])
                br[part[1]] = mut 
    return br

def BrRefine(argv) :
    args = parse_arg(argv)
    tree = Tree(args.tree, format=1)
    mut_branch = read_RecHMM(args.rec)
    
    for node in tree.traverse('postorder') :
        node.dist = mut_branch.get(node.name, 1e-9)
    if len(node.children) == 2 :
        rdist = sum([n.dist for n in node.children])/2
        for n in node.children :
            n.dist = rdist
    tree.write(format=1, outfile='{0}.mutation.labelled.nwk'.format(args.prefix))
    return

if __name__ == '__main__' :
    BrRefine(sys.argv[1:])

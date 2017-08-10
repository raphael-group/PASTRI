from fileio import *
from methods import read_vaf_file
from scipy.stats import binom

import numpy as np
def parse_arguments():
    '''
    Parses command line arguments.

    Returns:
        data_file:      <str> Filename containing variant allele frequencies
        proposal_file:  <str> Filename containing the initial proposal distribution
        tree_file:      <str> Perfect phylogeny matrices file
        tree_poo:       <int> The position of the tree of interest, such that 1 
                              is the highest likelihood tree in the result file. 
        output_prefix:  <str> Path and prefix for output data

    '''
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file", type=str)
    parser.add_argument("result_file", type = str)
    parser.add_argument("sample_file", type = str)
    parser.add_argument("-i", "--tree_pos", type=int, default=1) 
    parser.add_argument("-o", "--output_prefix", type=str, default="PASTRI")
    
    args = parser.parse_args()
    return args.data_file, args.result_file, args.sample_file, args.tree_pos, args.output_prefix


def get_tree(filename, pos):
    with open(filename) as f:
        for i in range(pos):
            try:
                T, T_header = read_in_matrix(f, int, True)
            except ValueError:
                raise ValueError("Tree index out of bounds")
    ti = int(T_header.split(":")[1])
    
    ll = float(T_header.split(":")[2])
    if ll == float('-inf'):
        print "Selected tree has a likelihood of 0. No samples found"
        exit()

    return T, ti

def get_F(filename, tree_index):
    with open(filename) as f:
        best_ll = float('-inf')
        best_F = None
        while True:
            try:
                F, header = read_in_matrix(f, float, True)
                ll, trees = header.split(':')[-1].split(',')
                ll = float(ll)
                if ll >= best_ll:
                    trees = map(int, trees.strip().split(' '))
                    if trees[tree_index] == 1:
                        best_ll = ll
                        best_F = F
            except ValueError:
                break
        return best_F

def get_maxll_cluster(a_vec, d_vec, F, num_clusts, num_samples):
    print a_vec, d_vec
    for a, d in zip(a_vec, d_vec):
        print a*2./d,
    print
    # Calculate the likelihood of it coming from any of the clusters
    cluster_likelihoods = []
    maxll = float("-inf")
    max_clust = None
    for i in range(num_clusts):
        clust_ll = 0
        for j in range(num_samples):
            freq = min(1,F.item((i,j))+0.00001)
            likelihood = binom.logpmf(a_vec[j],d_vec[j],freq/2.)
            clust_ll += likelihood
            print freq

        print i, clust_ll

        if clust_ll >= maxll:
            maxll = clust_ll
            max_clust = i
    
    return max_clust

if __name__ == '__main__':

    input_file, tree_likelihoods, sample_file, tree_pos, prefix = parse_arguments()
    A,D,num_snvs,num_samples = read_vaf_file(input_file)
    tree, tree_index = get_tree(tree_likelihoods, tree_pos)
    F = get_F(sample_file, tree_index)

    num_clusts, num_samples = F.shape
    assgmts=[]
    for i in range(num_snvs):
        a_vec = A[i,:]
        d_vec = D[i,:]

        assgmt = get_maxll_cluster(a_vec, d_vec, F, num_clusts, num_samples)
        assgmts.append(assgmt)

    with open(prefix+"."+str(tree_pos)+".C", 'w') as out:
        nodes = {}
        for i, v in enumerate(assgmts):
            if v not in nodes:
                nodes[v] = []
            nodes[v].append(i)
        for v in nodes:
            out.write(str(v)+"\t"+"\t".join(map(str,nodes[v]))+"\n")

    with open(prefix+"."+str(tree_pos)+".F", 'w') as out:
        write_matrix(F, "F", out)


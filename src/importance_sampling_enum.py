from generate_data import write_matrix
from io import read_in_trees
from gabow_myers import *
from methods import *

import numpy as np

import sys
import random
import heapq

def main():
    if len(sys.argv) <= 3: 
        print "Usage: python importance_sampling.py VAF_file pp_matrix_file num_iters estimate prefix"
        quit()
    vaf_file = sys.argv[1] 
    pp_matrices = sys.argv[2]
    num_iters = int(sys.argv[3])
    estimate_file = sys.argv[4]
    prefix = sys.argv[5]

    print vaf_file, pp_matrices, num_iters, estimate_file, prefix

    A,D, num_snvs, num_samples = read_vaf_file(vaf_file)
    alpha, beta, num_clusters = read_estimate_file(estimate_file)
    trees, num_chars = read_in_trees(pp_matrices)

    ### 
    #   Setup -- generate tree spectra for enumeration
    #   to avoid doing this repeatedly
    ### 
    tree_spectrums = {}
    for i, tree in enumerate(trees):
        child_list = c_list_from_matrix(tree)
        tree_spectrums[child_spectrum(child_list,0)]=i

    #print tree_spectrums, len(tree_spectrums)

    tree_likelihoods = [[] for i in range(len(trees))]
    max_ll_means = None
    max_ll = float('-inf')
    for i in range(num_iters):
        f_sample, sll = generate_sample(alpha, beta, num_clusters, num_samples)
        
        tree_counts = enumerate_trees(tree_spectrums, f_sample)
        if sum(tree_counts) > 0:
            ll = get_binomial_log_likelihood(A,D,f_sample,num_snvs,num_clusters,num_samples)
            llw = ll/sll # Weighted likelihood
            for j,v in enumerate(tree_counts):
                for k in range(v):
                    tree_likelihoods[j].append(ll)
            if ll > max_ll:
                max_ll = ll
                max_ll_means = f_sample

    best_tree = None
    best_tree_ll = float('-inf')

    for i in range(len(trees)):
        lls = tree_likelihoods[i]
        if len(lls) > 0: tree_ll = logsumexp(tree_likelihoods[i])
        else: tree_ll = float('-inf')
        if tree_ll > best_tree_ll:
            best_tree = trees[i]

    with open(prefix+".result", 'w') as out:
        write_matrix(best_tree, "Max LL tree", out)
        write_matrix(max_ll_means, "Max LL means", out)
    

if __name__=="__main__":
    main()







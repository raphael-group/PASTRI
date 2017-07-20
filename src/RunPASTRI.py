#from generate_data import read_in_trees, generate_random_usage_matrix, write_matrix

import numpy as np
import random
import heapq
import time

from methods import *
from fileio import *


###
#   Convert between representations of beta distribution
###
def get_f_matrix(alpha, beta):
    num_clusters, num_samples = alpha.shape
    newF = np.zeros((num_clusters, num_samples))
    sigma = np.zeros((num_clusters, num_samples))
    for i in range(num_clusters):
        for j in range(num_samples):
            b = float(beta.item((i,j)))
            a = float(alpha.item((i,j)))
            newF[i,j] = a / (a+b)
            sigma[i,j] = a+b

    return newF, sigma

def get_alpha_beta_matrix(f, sigma):
    num_clusters, num_samples = f.shape
    newA = np.zeros((num_clusters, num_samples))
    newB = np.zeros((num_clusters, num_samples))
    for i in range(num_clusters):
        for j in range(num_samples):
            alpha = max(sigma.item((i,j)) * f.item((i,j)), 0.0000001)
            beta = max(sigma.item((i,j))-alpha, 0.0000001)
            newA[i,j] = alpha
            newB[i,j] = beta
    return newA, newB

def update_proposal(cur_size, cur_ll, cur_distribution, sigma, f_sample, largest_tree, ll, base_distribution, base_size, base_ll, epsilon):
    #print ll, cur_ll
    if largest_tree > cur_size: 
        print "Bigger tree", largest_tree, cur_size
        alpha, beta = get_alpha_beta_matrix(f_sample, sigma)
        #print "Accepting new, larger size", largest_tree, cur_size
        return alpha, beta, largest_tree, ll 
    if ll != None and ll > cur_ll :
        alpha, beta = get_alpha_beta_matrix(f_sample, sigma)
        print "Accepting new, better likelihood"
        return alpha, beta, largest_tree, ll
    val = random.random()

    if val < epsilon:
        print "Restarting"
        alpha, beta = base_distribution
        #print "Restarting"
        return alpha, beta, base_size, base_ll
    else:
        alpha, beta = cur_distribution
        #print "Keeping"
        return alpha, beta, cur_size, cur_ll

def algorithm3(A, D, trees, alpha, beta, sigma, num_clusters, num_samples, num_iters, num_snvs, track_samples=False):
    ### 
    #   Setup -- generate tree spectra for enumeration
    #   to avoid doing this repeatedly
    ### 
    biggest_tree = 0
    tree_spectrums = {}
    for i, tree in enumerate(trees):
        child_list = c_list_from_matrix(tree)
        tree_spectrums[child_spectrum(child_list,0)]=i

    tree_likelihoods = [[] for i in range(len(trees))]
    max_ll_means = None
    max_ll = float('-inf')
    samples = []

    # Current Proposal consists of alpha, beta, a likelihood value, and a tree size
    f_start, sigma = get_f_matrix(alpha, beta)
    tree_counts, largest_tree = enumerate_trees(tree_spectrums, f_start)
    print tree_counts
    ll = get_binomial_log_likelihood(A,D,f_start,num_snvs,num_clusters,num_samples)
    base_size = largest_tree
    base_distribution = (alpha,beta)
    cur_size = largest_tree
    if cur_size == num_clusters:
        cur_ll = ll
    else:
        cur_ll = float("-inf")

    base_ll = cur_ll
    biggest_tree = cur_size

    print "Starting dist", f_start
    print "Starting lll:", cur_ll
    print "Starting size:", cur_size
    f_sample = f_start
    sll = cur_ll

    epsilon = 0.01
    for i in range(num_iters):
        if i%50 == 0: print i
        # Generate sample
        
        # Enumerate all trees that fit that sample
        tree_counts, largest_tree = enumerate_trees(tree_spectrums, f_sample)
        if sum(tree_counts) > 0:
            ll = get_binomial_log_likelihood(A,D,f_sample,num_snvs,num_clusters,num_samples)
            print "Found tree", ll, sll
            llw = ll
          
            for j,v in enumerate(tree_counts):
                for k in range(v):
                    tree_likelihoods[j].append(llw)
            if ll > max_ll:
                max_ll = ll
                max_ll_means = f_sample
            samples.append((ll, f_sample, ((1 if v > 0 else 0) for v in tree_counts)))
        else:
            ll = None

        alpha, beta, cur_size, cur_ll = update_proposal(cur_size, cur_ll, (alpha, beta), sigma, f_sample, largest_tree, ll, base_distribution, base_size, base_ll, epsilon)
        biggest_tree = max(biggest_tree, cur_size)

        f_sample, sll = generate_sample(alpha, beta, num_clusters, num_samples)

    best_tree = None
    best_tree_ll = float('-inf')

    likelihoods = [] 

    for i in range(len(trees)):
        lls = tree_likelihoods[i]
        if len(lls) > 0: tree_ll = logsumexp(tree_likelihoods[i])
        else: tree_ll = float('-inf')
        heapq.heappush(likelihoods, (-tree_ll, i, trees[i]))

    return likelihoods, samples

def read_in_files(vaf_file, estimate_file, pp_matrices):
    '''
    Read in input files.

    Returns:
        A, D:           <Numpy Matrix> Matrices containing the number of variant and total reads for the input data
        alpha, beta:    <Numpy Matrix> Matrices containing parameters for beta, for the initial proposal distribution
        trees:          <list of Numpy Matrix> Binary perfect phylogeny matrices
        num_snvs:       <int> The number of observed mutations
        num_samples:    <int> The number of samples
        num_clusters:   <int> The number of clusters as given in the proposal distribution

    '''
    A,D, num_snvs, num_samples = read_vaf_file(vaf_file)
    alpha, beta, num_clusters = read_estimate_file(estimate_file)
    trees, num_chars = read_in_trees(pp_matrices)
    assert(num_clusters == num_chars)
    return A, D, alpha, beta, trees, num_snvs, num_samples, num_clusters 

def parse_arguments():
    '''
    Parses command line arguments.

    Returns:
        data_file:      <str> Filename containing variant allele frequencies
        proposal_file:  <str> Filename containing the initial proposal distribution
        tree_file:      <str> Perfect phylogeny matrices file
        num_iters:      <int> The total number of iterations
        output_prefix:  <str> Path and prefix for output data

    '''
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file", type=str)
    parser.add_argument("proposal_file", type = str)
    parser.add_argument("tree_file", type=str)
    parser.add_argument("-n", "--num_iters", type=int, default=1000)
    parser.add_argument("-o", "--output_prefix", type=str, default=None)
    
    args = parser.parse_args()
    return args.data_file, args.proposal_file, args.tree_file, args.num_iters, args.output_prefix

def main():

    ###
    #   Initialize
    ###
 
    vaf_file, estimate_file, pp_matrices, num_iters, prefix = \
            parse_arguments()
    A, D, alpha, beta, trees, num_snvs, num_samples, num_clusters = \
            read_in_files(vaf_file, estimate_file, pp_matrices)
    

    # TODO: Replace
    sigma = np.zeros((num_clusters, num_samples))
    for i in range(num_clusters):
        for j in range(num_samples):
            sigma[i,j] = 10


    ###
    #   Run Algorithm
    ###
    likelihoods, samples = algorithm3(A, D, trees, alpha, beta, sigma, num_clusters, num_samples, num_iters, num_snvs)


    ###
    #   Write out Results
    ###

    with open(prefix+".trees", 'w') as out:
        i = 0
        while likelihoods:
            i += 1
            ll, index, tree = heapq.heappop(likelihoods)
            write_matrix(tree, str(i)+":"+str(index)+":"+str(-ll), out)

    with open(prefix+".fsamples", 'w') as out:
        for sample in samples:
            ll, f, trees = sample
            nrows, ncols = f.shape
            for i in range(nrows):
                for j in range(ncols):
                    f[i,j] = f.item((i,j)) * 2

            write_matrix(f, "LL:"+str(ll)+", " + " ".join(map(str, trees)), out)

if __name__=="__main__":
    main()







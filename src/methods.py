#from generate_data import read_in_trees, generate_random_usage_matrix, write_matrix
from gabow_myers import *
from fileio import *

import numpy as np

# A is the matrix of total read counts for all SNVs. A is num_snvs x num_samples
# D is the matrix of variant read counts for all SNVs. A is num_snvs x num_samples
# F is the matrix of cluster frequencies. F is num_chars x num_samples
# B is the perfect phylogeny matrix. B is num_chars x num_clones
# U is the usage matrix. U is num_clones x num_samples

from scipy.stats import binom
import math

def get_snv_log_likelihood(a_vec, d_vec, F, num_clusts, num_samples):
    # Calculate the likelihood of it coming from any of the clusters
    cluster_likelihoods = []
    for i in range(num_clusts):
        clust_ll = 0
        for j in range(num_samples):
            freq = min(1,F.item((i,j))+0.00001)
            likelihood = binom.logpmf(a_vec[j],d_vec[j],freq)
            clust_ll += likelihood

        if not(np.isnan(clust_ll) or clust_ll == float("-inf")):
            cluster_likelihoods.append(clust_ll)

    return logsumexp(cluster_likelihoods)

def get_binomial_log_likelihood(A,D,F,num_snvs,num_chars,num_samples):
    total_log_likelihood = 0
    for i in range(num_snvs):
        snv_likelihood = get_snv_log_likelihood(A[i,:], D[i,:], F, num_chars, num_samples)
        total_log_likelihood += snv_likelihood 
    return total_log_likelihood

def read_vaf_file(filename):
    with open(filename) as f:
        A_header = f.readline()
        num_snvs, num_samples = map(int, f.readline().strip()[1:-1].split(","))
        A = np.zeros((num_snvs, num_samples))
        
        for i in range(num_snvs):
            line = f.readline()
            line = map(int, map(float, line.strip().split()))
            for j,v in enumerate(line):
                A[i,j] = v

        f.readline() #blankline
        D_header = f.readline()
        num_snvs, num_samples = map(int, f.readline().strip()[1:-1].split(","))
        D = np.zeros((num_snvs, num_samples))
        
        for i in range(num_snvs):
            line = map(int, map(float, f.readline().strip().split()))
            for j,v in enumerate(line):
                D[i,j] = v

    return A,D,num_snvs, num_samples
from scipy.misc import logsumexp
import scipy.stats

def generate_sample(alpha, beta, num_clusters, num_samples):
    # generating samples from product of betas
    newF = np.zeros((num_clusters, num_samples))
    row_permutation=np.random.permutation(num_clusters)
    ll = 1
    for i in range(num_clusters):
        row = row_permutation[i]
        #row = i
        for j in range(num_samples):
            try:
                a = alpha.item((row,j))
                b = max(beta.item((row,j)), 0.00001)
                v = np.random.beta(a, b)
                ll += math.log(scipy.stats.beta.pdf(v,a,b))
            except:
                print "ALPHA"
                print alpha
                print "BETA"
                print beta
                raise
            newF[row,j] = v
    return newF, ll


def gen_child_list(tree):
    # Given a perfect phylogeny matrix, generate a dict that maps 
    # row indices to the row indices of their children
    # Row B is the child of Row A iff B and A differ in one position and in that position, A has a 0 and B has a 1

    tree = tree.transpose()
    
    child_list = {}
    for i, row1 in enumerate(tree.tolist()):
        clist = []
        for j, row2 in enumerate(tree.tolist()):
            num_diff = sum([1 if a != b else 0 for a,b in zip(row1, row2)])
            dir_diff = sum([ a - b for a,b in zip(row1, row2)])
            if num_diff == 1 and dir_diff < 0:
                clist.append(j)
        child_list[i] = clist
    return child_list
    

def check_sum_condition(child_list, F, num_chars, num_samples):
    for i in range(num_chars):
        children = child_list[i]
        for sample in range(num_samples):
            parent = F.item((i, sample))*2
            child_sum = sum([F.item((c, sample))*2 for c in children])
            if parent < child_sum:
                return False

    return True

                




from generate_data import read_in_trees, generate_random_usage_matrix, write_matrix
from gabow_myers import *

import numpy as np



def matrix_difference(matrixA, matrixB):

    # Need to do a maximum matching
    # First calculate the costs between any pair of rows


    rows, cols = matrixA.shape
    costs = {}
    
    for i in range(rows): #matrix A
        for j in range(rows): #matrix B
            cost = 0
            for k in range(cols):
                a = matrixA.item((i,k))
                b = matrixB.item((j,k)) 
                cost += abs(a-b)
            costs[i,j] = cost

    matching = []

    vals1 = range(rows)
    vals2 = range(rows)
    total_cost = 0
    while len(vals1) > 0:
        min_cost = float('inf')
        for v in vals1:
            for w in vals2:
                if costs[v,w] <= min_cost:
                    min_cost = costs[v,w]
                    min_val = (v,w)
        total_cost += min_cost
        vals1.remove(min_val[0])
        vals2.remove(min_val[1])

    return total_cost


def read_in_true(filename):
    with open(filename) as f:
        tree = read_in_matrix(f)
        U = read_in_matrix(f)
        F = read_in_matrix(f)
    return tree, U, F

def read_in_matrix(f):
        A_header = f.readline()
        num_snvs, num_samples = map(int, f.readline().strip()[1:-1].split(","))
        A = np.zeros((num_snvs, num_samples))
        
        for i in range(num_snvs):
            line = f.readline()
            line = map(float, line.strip().split())
            for j,v in enumerate(line):
                A[i,j] = v
        f.readline()
        return A

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

def read_estimate_file(filename):
    with open(filename) as f:
        alpha_header = f.readline()
        num_clusters, num_samples = map(int, f.readline().strip()[1:-1].split(","))
        alpha = np.zeros((num_clusters, num_samples))
        
        for i in range(num_clusters):
            line = f.readline()
            line = map(float, line.strip().split())
            for j,v in enumerate(line):
                alpha[i,j] = max(0.00001,v)

        f.readline() #blankline
        beta_header = f.readline()
        num_clusters, num_samples = map(int, f.readline().strip()[1:-1].split(","))
        beta = np.zeros((num_clusters, num_samples))

        
        for i in range(num_clusters):
            line = map(float, f.readline().strip().split())
            print line
            for j,v in enumerate(line):
                beta[i,j] = max(0.00001, v)
        return alpha, beta, num_clusters


from scipy.misc import logsumexp
def get_tree_likelihood(tree, A, D, num_snvs, num_chars, num_samples, precision):
    sample_lls = []
    for i in range(precision):
        F = random_F(tree, num_chars, num_samples, i%num_samples)
        sample_lls.append(get_binomial_log_likelihood(A,D,F, num_snvs, num_chars, num_samples))
    return logsumexp(sample_lls)
import scipy.stats

def generate_sample(alpha, beta, num_clusters, num_samples):
    # generating samples from product of betas
    newF = np.zeros((num_clusters, num_samples))

    #print alpha
    #print beta
    #print alpha.shape
    #print beta.shape
    #print (num_clusters, num_samples)
    row_permutation=np.random.permutation(num_clusters)
    ll = 1
    for i in range(num_clusters):
        row = row_permutation[i]
        #row = i
        for j in range(num_samples):
            try:
                a = alpha.item((row,j))
                b = max(beta.item((row,j)) - a, 0.00001)
                f=10.
                v = np.random.beta(a/f, b/f)
                ll += math.log(scipy.stats.beta.pdf(v,a/f,b/f))
            except:
                #print alpha.item((i,j))
                #print beta.item((i,j))
                print "ALPHA"
                print alpha
                print "BETA"
                print beta
                raise
            #v = alpha.item((i,j))/(beta.item((i,j)) + alpha.item((i,j)))
            #print v
            newF[row,j] = v
    #print newF
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

                




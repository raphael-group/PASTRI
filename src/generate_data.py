import numpy as np
import random

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

            #print row1, row2, num_diff, dir_diff
            # if a is 0 and b is 1, then row2 is child. Thus, dir_diff is negative
            if num_diff == 1 and dir_diff < 0:
                clist.append(j)
        child_list[i] = clist
    #print child_list
    return child_list
    

def check_sum_condition(child_list, F, num_chars, num_samples):
    #print F
    #print num_chars
    #print num_samples
    for i in range(num_chars):
        children = child_list[i]
        #print children
        for sample in range(num_samples):
            parent = F.item((i, sample))*2
            child_sum = sum([F.item((c, sample))*2 for c in children])
            #print parent, child_sum
            if parent < child_sum:
                return False

    #print "----------------"
    return True


def read_in_trees(filename):
    all_trees = []
    with open(filename) as f:
        num_chars, num_trees = map(int, f.readline().strip().split(" "))
        trees = f.readlines()
        for j in range(num_trees):
            treestart = j*((num_chars)+1)
            treeend = treestart + num_chars
            tree = np.matrix([map(int, line.strip().split(" ")) for line in trees[treestart:treeend]])
            all_trees.append(tree)
    return all_trees, num_chars

def generate_random_usage_matrix(num_chars, num_samples):
    alpha = [1]*num_chars
    U = np.random.dirichlet(alpha,num_samples).transpose()
    rows, cols = U.shape
    for r in range(rows):
        for c in range(cols):
            U[r,c] /= 2.
 
    return U

def generate_vafs(num_snvs, num_chars, num_samples, depth, F):
    cluster_assignment_dist = np.random.dirichlet([1]*num_chars, 1)
    clust_assignments = np.random.choice(num_chars, num_snvs, p = cluster_assignment_dist[0])
    D = np.ones((num_snvs, num_samples))
    A = np.ones((num_snvs, num_samples))
    for i in range(num_snvs):
        cluster = clust_assignments[i]
        for j in range(num_samples): 
            # The min is here because otherwise I get floating point errors 
            # when the frequency is actually 1
            clust_freq = min(F.item((cluster,j)), 1) 
            D[(i,j)] = depth
            A[(i,j)] = np.random.binomial(depth, clust_freq, 1)
             
    return A, D, clust_assignments

def write_out_result(tree, U, F, A, D, clust_F, clust_assignments, out_prefix):
    with open(out_prefix+".true", 'w') as out:
        #Write out B, U, and F
        write_matrix(tree, "B", out)
        write_matrix(U, "U", out)
        write_matrix(F, "F", out)
        write_matrix(clust_F, "Cluster Frequencies", out)
        write_vector(clust_assignments, "Cluster Assignments", out)

    with open(out_prefix+".input", 'w') as out:
        write_matrix(A, "A", out)
        write_matrix(D, "D", out)

def write_out_estimate(F, depth, out_prefix):
    num_clusters, num_samples = F.shape
    alpha = F * depth
    beta = depth*np.ones((num_clusters, num_samples)) - alpha
    with open(out_prefix+".estimate", 'w') as out:
        write_matrix(alpha, "Alpha", out)
        write_matrix(beta, "Beta", out)

def write_matrix(matrix,name, stream):
    stream.write("> " + name + "\n")
    
    if matrix == None:
        stream.write(str((0,0))+"\n")
    else:
        stream.write(str(matrix.shape)+"\n")
        for line in matrix:
            stream.write(str(line).replace("\n", " ")[2:-2]+"\n")

    stream.write("\n")

def write_vector(vector, name, stream):
    stream.write("> " + name + "\n")
    stream.write(str(vector.shape)+"\n")
    stream.write(str(vector)[1:-1]+"\n")
    
def main():
    import sys
    
    if len(sys.argv) < 6:
        print "Usage: python generate_data.py pp_matrix_file num_snvs num_samples depth out_prefix (seed)"
        quit()
    
    pp_matrix_file = sys.argv[1]
    num_snvs = int(sys.argv[2])
    num_samples = int(sys.argv[3])
    depth = int(sys.argv[4])
    out_prefix = sys.argv[5]

    # 0) Set random seeds
    if len(sys.argv) > 6:
        seed = int(sys.argv[6])
        random.seed(seed)
        np.random.seed(seed)
    
    # 1) Read in perfect phylogeny matrices
    trees, num_chars = read_in_trees(pp_matrix_file)
    
    # 2) Select a perfect phylogeny matrix at random
    tree = random.choice(trees)
    
    # 3) Generate a random usage matrix
    U = generate_random_usage_matrix(num_chars, num_samples)
    
    # 4) Calculate the resultant F matrix
    F =  tree*U
    
    # 5) Generate a random cluster assignment probabilty distribution (using a dirichlet)
    # For num_snvs snvs, randomly assign to a cluster, then for num_samples samples, draw from binomial to get frequencies
    A,D,clust_assignments = generate_vafs(num_snvs, num_chars, num_samples, depth, F)

    counts = [sum([(1 if c == i else 0) for c in clust_assignments]) for i in range(num_chars)] 

    clust_F = np.zeros((num_chars, num_samples))

    for j in range(num_samples):
        for i in range(num_chars):

            var = 0
            tot = 0
            for k, c in enumerate(clust_assignments):
                if c == i:
                    var += A.item((k,j))
                    tot += D.item((k,j))
            try:
                clust_F[i,j] = var*1.0/(tot)
            except:
                clust_F[i,j] = 0

    clist = gen_child_list(tree)
    check = check_sum_condition(clist, clust_F, num_chars, num_samples)
    #

    #if not check:
    write_out_result(tree, U, F, A, D, clust_F, clust_assignments, out_prefix)
    write_out_estimate(F, depth, out_prefix)
    #else:
    #    print "Meets Sum Condition"
    #    exit(1)


if __name__ == '__main__':
    main()

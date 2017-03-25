import numpy as np
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
                beta[i,j] = max(0.00001, v - alpha.item((i,j)))

    
    print "---------------------------"
    print "ALPHA:"
    print alpha
    print "BETA:"
    print beta
    print "---------------------------"
    return alpha, beta, num_clusters

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
    stream.write(u"> " + name + "\n")
    
    if matrix == None:
        stream.write(str((0,0))+u"\n")
    else:
        stream.write(str(matrix.shape)+u"\n")
        for line in matrix:
            stream.write(str(line).replace("\n", " ")[2:-2]+u"\n")

    stream.write(u"\n")

def write_vector(vector, name, stream):
    stream.write(u"> " + name + "\n")
    stream.write(str(vector.shape)+u"\n")
    stream.write(str(vector)[1:-1]+u"\n")

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


    



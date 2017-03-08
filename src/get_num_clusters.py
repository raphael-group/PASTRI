import numpy as np
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


import sys
filename = sys.argv[1]
with open(filename) as f:
    num_clusters, num_samples = read_in_matrix(f).shape

print num_clusters
    

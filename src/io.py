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



import sys
import copy 
import numpy as np
from io import *

def space(tree_nodes):
    return "\t"*tree_nodes

def enumerate_trees(tree_spectrums, freqs):
    # Create ancestry graph with F
    # Enumerate all spanning trees
    # Report (t,l) pairs

    graph = {}
    for i, row1 in enumerate(freqs):
        graph[i] = []
        for j, row2 in enumerate(freqs):
            if i == j: continue
            all_gt = True
            for v,w in zip(row1,row2):
                if v < w:
                    all_gt = False
                    break
            if all_gt: graph[i].append(j)
    num_verts, num_samples = freqs.shape
    all_trees = []
    for i in range(num_verts):
        frontier = []
        root = i
        for edge in graph[root]:
            frontier.append((root, edge))
        tree_nodes = [0]*num_verts
        tree_nodes[i] = 1
        root_i_trees = grow([], tree_nodes, num_verts, frontier, graph, freqs, [], root)
        all_trees+=[(tree,i) for tree in root_i_trees]
    #print all_trees

    counts = [0]*len(tree_spectrums)

    # Now I have all the trees, I need to return unlabeled trees and permutations       
    #print "----------"
    #print all_trees
    # For any tree, find the topology that it corresponds to
    for (tree, i) in all_trees:
        #convert_to_pp(tree, 1, 4)
        child_list = c_list_from_edge_list(tree, num_verts)
        spectrum = child_spectrum(child_list, i)
        counts[tree_spectrums[spectrum]] += 1

    return counts
 
def child_spectrum(child_list, root):
    val = spectrum_recurse(child_list, root)
    spectrum = sorted(spectrum_recurse(child_list, root), key=len)
    return tuple(spectrum)

def spectrum_recurse(child_list, i):
    if len(child_list[i]) == 0: return ()
    else: 
        spectrum = tuple(sorted([spectrum_recurse(child_list,j) for j in child_list[i]]))
        return spectrum

def c_list_from_matrix(tree):
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
    
def c_list_from_edge_list(tree, verts):
    child_list = {}
    for i in range(verts):
        child_list[i] = []
    for edge in tree:
        s,t = edge
        child_list[s].append(t)
    return child_list 


def convert_to_pp(tree, root, verts):
    # first convert to childlist
    child_list = {}
    for i in range(verts):
        child_list[i] = []
    for edge in tree:
        s,t = edge
        child_list[s].append(t)
    #print child_list

    matrix = np.zeros((verts,verts))
    matrix[root,root] = 1
    for c in child_list[root]: recurse_convert_to_pp(matrix, child_list, c, root, verts)
    #print matrix

def recurse_convert_to_pp(matrix, child_list, i, pi, verts):
    # Add i and all the children of i to matrix
    for j in range(verts):
        matrix[i,j] = matrix[pi,j]
    matrix[i,i] = 1
    for c in child_list[i]: recurse_convert_to_pp(matrix, child_list, c, i, verts)



def canoncial_tree(matrix):
    pass     
    

def grow(tree, tree_nodes, k, F, graph, freqs, trees, root):
    z = sum(tree_nodes)
    if sum(tree_nodes) == k:
           
        #print space(z), "FOUND TREE"

        #print space(z), tree
        trees.append(tree[:])
        
    else:
        FF = []
        while len(F) > 0:
            assert(test_frontier(tree_nodes, graph, F))
            assert(test_tree(tree,tree_nodes, root))
            
            # Add a new node to the tree
            edge = F.pop()
            start, end = edge
            tree.append(edge)
            tree_nodes[end] = 1

            removed_edges = []
            
            F = construct_frontier(tree_nodes, tree, graph, freqs)
            trees = grow(tree, tree_nodes, k, F, copy.deepcopy(graph), freqs, trees, root)
            graph[start].remove(end)
    
            # Remove the node you just added from the tree

            tree_nodes[end] = 0
            popped_edge = tree.pop()
            assert(popped_edge == edge)
            F = construct_frontier(tree_nodes, tree, graph, freqs)
    return trees
    
def construct_frontier(tree_nodes, tree, graph, freqs):
    # Need to check sum condition across samples
    num_nodes, num_samples = freqs.shape
    F = []
    for s in graph.keys():
        Ts = graph[s]
        slack = [0]*num_samples
        for j in range(num_samples):
            
            par_freq = freqs.item((s,j))
            child_set = [t for (s2,t) in tree if s == s2]
            child_freq = sum([ freqs.item((t,j)) for t in child_set])
            slack[j] = par_freq - child_freq
            assert(round(slack[j],10) >= 0)

        if tree_nodes[s] == 1:
             for t in Ts:
                 # If the tree 
                 if tree_nodes[t] == 0:
                    fits = True
                    for j in range(num_samples):
                        t_freq = freqs.item((t,j))
                        if round(t_freq,10) > round(slack[j],10):
                            fits = False 
                            break
                    if fits: F.append((s,t))
    return F


#####################################
###         Testing Code          ###
#####################################
def test_frontier(tree_nodes, graph, F):
    # Check that all source nodes in frontier are in tree
    # Check that all target nodes in frontier are not in tree
    for edge in F:
        s,t = edge 
        if tree_nodes[s] == 0: return False
        if tree_nodes[t] == 1: return False
    # Check that all edges in the graph that meet these properties 
    # are in the frontier
    #for s in graph.keys():
    #    Ts = graph[s]
    #    if tree_nodes[s] == 0:
    #        for t in Ts:
    #            if (s,t) in F: return False
    #    elif tree_nodes[s] == 1:
    #         for t in Ts:
    #             if tree_nodes[t] == 0 and (s,t) not in F: return False
    #             if tree_nodes[t] == 1 and (s,t) in F: return False
    return True

def test_tree(tree, tree_nodes, root):
    # Checks that tree and tree_nodes are consistent
    #print tree, tree_nodes

    # Construct new tree_nodes
    new_tree_nodes = [0] * len(tree_nodes)
    for edge in tree:
        s,t = edge
        new_tree_nodes[s] = 1
        new_tree_nodes[t] = 1
    new_tree_nodes[root] = 1


    for v,w in zip(tree_nodes, new_tree_nodes):
        if v != w: return False
    return True


# 0 is a dummy node with edges to every other node

#filename = sys.argv[1]
#trees = read_in_trees(filename) 

#tree_spectrums = {}
# Calculate child spectrum of trees
#for i, tree in enumerate(trees[0]):
#        child_list = c_list_from_matrix(tree)
#        tree_spectrums[child_spectrum(child_list,0)]=i

#print tree_spectrums


#print enumerate_trees(tree_spectrums, freqs)





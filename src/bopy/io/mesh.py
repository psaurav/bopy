import os
import sys
import numpy as np

def read_mesh_nodes(filename):
    """
    Read the mesh nodes

    Args:
        filename
    Returns:
        np.array - size (n, 3)
    """
    f = open(filename, 'r')

    # header
    l = f.readline().rstrip().split()    # header line
    if l[0] != 'MeshData' or l[1] != '5.0':
        print >>sys.stderr, l[0], l[1]
        raise Exception('{0} does not have recognizable header'.
                         format(filename))
    l = f.readline()            # empty line

    # Nodelist 
    l = f.readline().rstrip().split()    # NodeList
    if l[0] == 'NodeList':
        nn = int(l[1])
    else:
        raise Exception('{0} could not find NodeList'.
                         format(filename))

    # read the nodes
    nodes = []
    for i in range(nn):
        l = f.readline().rstrip()
        node = [float(x) for x in l.split('[')[1].split(']')[0].split()]
        nodes.append(node)
    f.close()
    return np.array(nodes)

def read_mesh_bbox(filename=None, nodes=None):
    """
    Get the bounding box of the mesh.  
    Essentially the min and max of node positions.

    Args:
        filename - mesh file name [default: None]
        nodes - np.array (n, 3) [default: None]
    Returns:
        np.array - (2, 3) 
    """
    if filename != None:
        nodes = read_mesh_nodes(filename)
    return np.array([[np.min(nodes[:,0]), np.min(nodes[:,1]), np.min(nodes[:,2])],
                     [np.max(nodes[:,0]), np.max(nodes[:,1]), np.max(nodes[:,2])]])

def read_mesh_nodes_bbox(filename):
    """
    Get the nodes as well as the bbox

    Args:
        filename - Mesh file name

    Returns:
        nodes - the node positions (type np.array)
        bbox - the bounding box (type np.array)
    """
    nodes = read_mesh_nodes(filename)
    bbox = read_mesh_bbox(nodes=nodes)
    return nodes, bbox

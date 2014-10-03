import os
import sys
import numpy as np

def node_in_breast(nodes, bs, method='cubic'):
    """
    Returns an bool array 
    
    Args:
        nodes - mesh nodes np.array (nn, 3)
        bs - breast and chest shape 
        method - interpolation method [default: cubic]
    """
    from scipy.interpolate import griddata
    #print bs[:,0:2]
    #print bs[0:,2]
    zn = griddata(bs[:,0:2], bs[:,2], nodes[:,0:2], method=method)
    #print zn
    return nodes[:,2] < zn

def add_sphere_to_nim(meshfile, nimfile, val, c0, r, outfilename):
    """
    Add a sphere to an existing NIM file

    Args:
        meshfile - meshfile name
        nimfile - an existing nim file,
        val - new value
        c0 - center of sphere (3D tupple, x0, y0, z0)
        r - radius of sphere
        outfilename - name of output NIM file
    """
    import bopy.io as bpio
    nodes = bpio.read_mesh_nodes(meshfile) 
    nim = bpio.read_nim(nimfile)
    c0 = np.array(c0)
    r2 = r**2
    insphere = np.sum((nodes-c0)**2) <= r2
    nim[insphere] = val
    bpio.write_nim(outfilename, nim, meshname=meshfile)

def create_2region_nims(meshfile, bsfile, bval, rval, outfilename):
    """
    Create two region NIM files 

    Args:
        meshfile - meshfile name
        bsfile - breast shape filename name [.npy]
        vbal - the value in the breast/chest region
        rval - the value in outside the breast region
        outfilename - output a NIM file
    """
    bs = np.load(bsfile)
    import bopy.io as bpio
    nodes = bpio.read_mesh_nodes(meshfile) 
    nib = node_in_breast(nodes, bs)
    vect = np.ones(nib.shape)*rval
    vect[nib] = bval
    bpio.write_nim(outfilename, vect, meshname=meshfile)

def create_chestwall_nim(meshfile, zchest, cval, rval, outfilename):
    """
    Create two region (chestwall) NIM files 

    Args:
        meshfile - meshfile name
        zchest - chest position
        cval - the value in the chest region
        rval - the value in outside the chest region
        outfilename - output a NIM file
    """
    import bopy.io as bpio
    nodes = bpio.read_mesh_nodes(meshfile) 
    nic = nodes[:,2] 
    vect = np.ones(nic.shape)*rval
    vect[nic] = cval
    bpio.write_nim(outfilename, vect, meshname=meshfile)
    
    
def nim2grid(meshfile, nimfile, mm=2, method='cubic'):
    """
    Interpolate the mesh values on a regular grid, for visualization

    Args:
        meshfile - mesh file name
        nimfile - NIM file name
        mm - the grid spacing size in mm [default: 2mm]
        method - interpolation method [default: cubic] 

    Returns:
        np.array (3D)
    """
    import bopy.io as bpio
    nodes, bbox = bpio.read_mesh_nodes_bbox(meshfile)
    nx = int(bbox[1,0] - bbox[0,0])/mm
    ny = int(bbox[1,1] - bbox[0,1])/mm
    nz = int(bbox[1,2] - bbox[0,2])/mm
    xgrid = np.linspace(bbox[0,0], bbox[1,0], nx+1)
    ygrid = np.linspace(bbox[0,1], bbox[1,1], ny+1)
    zgrid = np.linspace(bbox[0,2], bbox[1,2], nz+1)
    ix, iy, iz = np.indices([len(xgrid), len(ygrid), len(zgrid)])
    xgrid = np.ravel(xgrid[ix])
    ygrid = np.ravel(ygrid[iy])
    zgrid = np.ravel(zgrid[iz])
    xyzgrid = np.vstack((xgrid, ygrid, zgrid)).T
    nim = bpio.read_nim(nimfile)
    from scipy.interpolate import griddata
    return griddata(nodes, nim, xyzgrid, method=method)

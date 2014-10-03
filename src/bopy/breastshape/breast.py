#!/usr/bin/python

import os
import sys
import numpy as np
import numpy.ma as ma
from mayavi import mlab
import nxutils.nxutils as nx

def poly33val(x,y,p):
    b00=p[0]
    b10=p[1]
    b01=p[2]
    b20=p[3]
    b11=p[4]
    b02=p[5]
    b30=p[6]
    b21=p[7]
    b12=p[8]
    b03=p[9]
    z=b00+b10*x+b01*y+b20*x**2+b11*x*y+b02*y**2+b30*x**3+b21*x**2*y+b12*x*y**2+b03*y**3
    return z

def poly33(xy,b00,b01,b10,b20,b11,b02,b30,b21,b12,b03):
    x=xy[0]
    y=xy[1]
    p=[b00,b01,b10,b20,b11,b02,b30,b21,b12,b03]
    return poly33val(x,y,p)

def poly4val(x,p):
    b0=p[0];b1=p[1];b2=p[2];b3=p[3];b4=p[4]
    y=b0+b1*x+b2*x**2+b3*x**3+b4*x**4
    return y

def poly4(x,b0,b1,b2,b3,b4):
    p=[b0,b1,b2,b3,b4]
    return poly4val(x,p)

def poly6val(x, p):
    b0 = p[0]
    b1 = p[1]
    b2 = p[2]
    b3 = p[3]
    b4 = p[4]
    b5 = p[5]
    b6 = p[6]
    y = b0 + b1*x + b2*x**2 + b3*x**3 + b4*x**4 + b5*x**5 + b6*x**6
    return y

def poly6(x, b0, b1, b2, b3, b4, b5, b6):
    p = [b0, b1, b2, b3, b4, b5, b6]
    return poly6val(x, p)

class BreastResult:

    def __init__(self, filename):
        surface_result = np.load(filename)
        shape = surface_result['shape']
        self.x_all, self.y_all, self.z_all=shape[:,0],shape[:,1],shape[:,2]
        self.right_surface = surface_result['right_surface'][0]
        self.chestwall = self.right_surface['z_min']
        self.left_surface = surface_result['left_surface'][0]
        self.bottom_surface = surface_result['bottom_surface']
        self.p_r = self.right_surface['p']
        self.p_l = self.left_surface['p']
        self.dense = np.float(surface_result['Param_0'][0]['dense'])
        self.y_min = np.float(self.right_surface['y_min'])
        self.y_max = np.float(self.right_surface['y_max'])
        self.y_grid = np.arange(self.y_min, self.y_max+self.dense, self.dense)

        self.x_all_max = np.max(self.x_all)
        self.x_all_min = np.min(self.x_all)
        self.z_all_max = np.max(self.z_all)

    def isInsideBreast(self, x0, y0, z0, zchestwall=-np.inf, z_offset=0.0):
        x0 = np.float(x0)
        y0 = np.float(y0)
        z0 = np.float(z0)
        z_offset = np.float(z_offset)
        z_pp = z0 - z_offset
        x_pp = x0          

        if z_pp < zchestwall:
            return True
        
        if x_pp > self.x_all_max or \
           x_pp < self.x_all_min or \
           z_pp < self.chestwall or \
           z_pp > self.z_all_max:  
            return False
    
        delta_y = np.abs(y0 - self.y_grid)   
        index = np.where(delta_y == delta_y.min())[0][0]
        y_pp = self.y_grid[index]
        bottom_curve = self.bottom_surface[index]
        x_min_b = bottom_curve['x_min']         
        x_max_b = bottom_curve['x_max']        
        p_b = bottom_curve['p']               
        z_max_r = poly6val(x_max_b,p_b)      
        z_max_l = poly6val(x_min_b,p_b)     

        x_pp = np.float(x_pp) 
        y_pp = np.float(y_pp) 
        z_pp = np.float(z_pp)              
        flag = False                      

        # Method 1 point in polygon
        delta = 2*self.dense
        xpoints = np.array([])
        zpoints = np.array([])
        z_tmp = np.arange(self.chestwall, z_max_l, delta)
        y_tmp = y_pp*np.ones_like(z_tmp)
        x_tmp = poly33val(y_tmp, z_tmp, self.p_l)
        x_min_all = x_tmp[0]
        xpoints = np.append(xpoints, x_tmp)
        zpoints = np.append(zpoints, z_tmp)
        x_tmp = np.arange(x_min_b, x_max_b, delta)
        z_tmp = poly6val(x_tmp, p_b)
        xpoints = np.append(xpoints, x_tmp)
        zpoints = np.append(zpoints, z_tmp)
        z_tmp = np.arange(z_max_r, self.chestwall+delta, -delta)
        y_tmp = y_pp*np.ones_like(z_tmp)
        x_tmp = poly33val(y_tmp, z_tmp, self.p_r)
        x_max_all = x_tmp[-1]
        xpoints = np.append(xpoints, x_tmp)
        zpoints = np.append(zpoints, z_tmp)
        x_tmp = np.arange(x_max_all-self.dense, x_min_all, -delta)
        z_tmp = self.chestwall*np.ones_like(x_tmp)
        xpoints = np.append(xpoints, x_tmp)
        zpoints = np.append(zpoints, z_tmp)
        xzpoints = np.vstack((xpoints, zpoints))
        point = np.array([[x_pp, z_pp]])
        flag = nx.points_inside_poly(point, np.transpose(xzpoints))
        return flag

def isInsideBreast(SurfaceResult, x0, y0, z0, z_offset=0.0):
    x0=np.float(x0); y0=np.float(y0); z0=np.float(z0)       #
    z_offset=np.float(z_offset)                             #
    z_pp=z0-z_offset                                        #
    x_pp=x0                                                 #
    #SurfaceFittingResultFilePath=r'E:/Data/Thesis Data/20130410_G3033_Patient/20130410_cam0_proj_calib_00_reconstruction.npz'
    #SurfaceResult=np.load(SurfaceFittingResultFilePath)
    shape=SurfaceResult['shape']
    x_all,y_all,z_all=shape[:,0],shape[:,1],shape[:,2]
    right_surface=SurfaceResult['right_surface'][0]
    chestwall=right_surface['z_min']

    if x_pp>x_all.max() or x_pp<x_all.min() or z_pp<chestwall or z_pp>z_all.max():  #
        return False                                                                #

    left_surface=SurfaceResult['left_surface'][0]
    bottom_surface=SurfaceResult['bottom_surface']
    p_r=right_surface['p']
    p_l=left_surface['p']
    dense=np.float(SurfaceResult['Param_0'][0]['dense'])
    y_min=np.float(right_surface['y_min'])
    y_max=np.float(right_surface['y_max'])
    y_grid=np.arange(y_min,y_max+dense,dense)
    delta_y=np.abs(y0-y_grid)                                                       #
    index=np.where(delta_y==delta_y.min())[0][0]                                    #
    y_pp=y_grid[index]                                                              #
    bottom_curve=bottom_surface[index]                                              #
    x_min_b=bottom_curve['x_min']                                                   #
    x_max_b=bottom_curve['x_max']                                                   #
    p_b=bottom_curve['p']                                                           #
    z_max_r=poly6val(x_max_b,p_b)                                                   #
    z_max_l=poly6val(x_min_b,p_b)                                                   #

    x_pp=np.float(x_pp); y_pp=np.float(y_pp); z_pp=np.float(z_pp)                   #
    flag=False                                                                      #

    # Method 1 point in polygon
    delta=2*dense
    xpoints=np.array([])
    zpoints=np.array([])
    z_tmp=np.arange(chestwall,z_max_l,delta)
    y_tmp=y_pp*np.ones_like(z_tmp)
    x_tmp=poly33val(y_tmp,z_tmp,p_l)
    x_min_all=x_tmp[0]
    xpoints=np.append(xpoints,x_tmp)
    zpoints=np.append(zpoints,z_tmp)
    x_tmp=np.arange(x_min_b,x_max_b,delta)
    z_tmp=poly6val(x_tmp,p_b)
    xpoints=np.append(xpoints,x_tmp)
    zpoints=np.append(zpoints,z_tmp)
    z_tmp=np.arange(z_max_r,chestwall+delta,-delta)
    y_tmp=y_pp*np.ones_like(z_tmp)
    x_tmp=poly33val(y_tmp,z_tmp,p_r)
    x_max_all=x_tmp[-1]
    xpoints=np.append(xpoints,x_tmp)
    zpoints=np.append(zpoints,z_tmp)
    x_tmp=np.arange(x_max_all-dense,x_min_all,-delta)
    z_tmp=chestwall*np.ones_like(x_tmp)
    xpoints=np.append(xpoints,x_tmp)
    zpoints=np.append(zpoints,z_tmp)
    xzpoints=np.vstack((xpoints,zpoints))
    point=np.array([[x_pp,z_pp]])
    flag=nx.points_inside_poly(point, np.transpose(xzpoints))

    return flag


def read_breast_shape(filename):
    return np.load(filename)['shape']

def breast_and_chestwall(filename, xx, xm, yx, ym, zm, dx, dy):
    sh = read_breast_shape(filename)
    print >>sys.stderr, 'sh.shape =', sh.shape
    print >>sys.stderr, '0min =', np.min(sh[:,0]), '  0max =', np.max(sh[:,0])
    print >>sys.stderr, '1min =', np.min(sh[:,1]), '  1max =', np.max(sh[:,1])
    print >>sys.stderr, '2min =', np.min(sh[:,2]), '  2max =', np.max(sh[:,2])
    x = sh[:,0]
    y = sh[:,1]
    z = sh[:,2]
    
    x_set = np.unique(x)
    y_set = np.unique(y)
    z_set = np.unique(z)

    print >>sys.stderr, 'len(x_set) =', len(x_set)
    print >>sys.stderr, 'len(y_set) =', len(y_set)
    print >>sys.stderr, 'len(z_set) =', len(z_set)
    print >>sys.stderr, 'len(x) =', len(x)

    # create a nx2 array of points (x, y)
    d_xy = np.vstack((x, y)).T
    # create a grid of points for interpolation 
    xn = np.linspace(xm, xx, dx+1)
    yn = np.linspace(ym, yx, dy+1)
    n_xy = np.meshgrid(xn, yn)

    #zn = griddata(d_xy, z, n_xy, method='cubic', fill_value=-50.0)
    from scipy.interpolate import griddata
    zn = griddata(d_xy, z, n_xy, method='cubic')
    xn = n_xy[0]
    yn = n_xy[1]
    zshape = zn.shape

    # correct for missing last
    #zn[dy-1,:] = zn[dy-2,:]
    zn[dy,:] = zn[dy-1,:]

    znew = []
    from scipy.optimize import curve_fit
    for i in range(len(yn)):  
        xt = ma.compressed(ma.masked_array(xn[i,:], np.isnan(zn[i,:])))
        zt = ma.compressed(ma.masked_array(zn[i,:], np.isnan(zn[i,:])))
        p, pcov = curve_fit(poly6, xt, zt)
        zt = poly6val(xn[i,:], p)
        znew = np.hstack((znew, zt))
    znew[znew < zm] = zm 
    znew = np.array(znew)

    # unravel for display 
    xn = np.ravel(n_xy[0])
    yn = np.ravel(n_xy[1])
    zn = np.ravel(zn)
    
    # display
    mlab.points3d(x, y, z, scale_factor=1.0,color=(0,1,0))
    mlab.points3d(xn, yn, zn, scale_factor=1.0,color=(1,0,0))
    mlab.points3d(xn, yn, znew, scale_factor=1.0,color=(0,0,1))
    mlab.axes()
    mlab.show()

    return (xn, yn, zn)

def driver_breast_chestwall(*agrv):
    from bopy.utils import ParamList
    pl = ParamList()
    for a in args:
        pl.read_from_file(a)

    xx = float(pl.get_val('MESH_XMAX'))
    yx = float(pl.get_val('MESH_YMAX'))
    xm = float(pl.get_val('MESH_XMIN'))
    ym = float(pl.get_val('MESH_YMIN'))
    zm - float(pl.get_val('MESH_CHESTWALL'))
    dx = int(xx-xm)
    dy = int(yx-ym)
    bc_shape = breast_and_chestwall(filename, xx, xm, yx, ym, zm, dx, dy)
    d = np.vstack((bc_shape[0], bc_shape[1], bc_shape[2])).T
    np.save('breast_and_chestwall.npy', d)
    

def main():
    try:
        filename = sys.argv[1]
    except:
        sys.exit('Usage: {0} <filename>'.format(os.path.basename(sys.argv[0])))

    xx = 130.0
    xm = -130.0
    yx = 60.0
    ym = 0.0
    zm = -50.0
    dx = int(xx-xm)
    dy = int(yx-ym)
    bc_shape = breast_and_chestwall(filename, xx, xm, yx, ym, zm, dx, dy)
    d = np.vstack((bc_shape[0], bc_shape[1], bc_shape[2])).T
    np.savez('breast_and_chestwall.npy', d)

if __name__ == '__main__':
    main()

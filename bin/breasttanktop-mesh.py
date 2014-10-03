#!/usr/bin/python

import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

class struct():
    def __init__(self, *args, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

def draw_box(axes, bbox, color, alpha=1):

    xmin = bbox[0]
    xmax = bbox[1]
    ymin = bbox[2]
    ymax = bbox[3]
    zmin = bbox[4]
    zmax = bbox[5]

    xs = [xmin, xmin, xmin, xmin, xmin, xmax, xmax, xmax, xmax, xmax]
    ys = [ymin, ymax, ymax, ymin, ymin, ymin, ymax, ymax, ymin, ymin]
    zs = [zmin, zmin, zmax, zmax, zmin, zmin, zmin, zmax, zmax, zmin]

    x2 = [xmin, xmax]
    y2 = [ymax, ymax]
    z2 = [zmax, zmax]

    x3 = [xmin, xmax]
    y3 = [ymin, ymin]
    z3 = [zmax, zmax]

    x4 = [xmin, xmax]
    y4 = [ymax, ymax]
    z4 = [zmin, zmin]

    axes.plot(xs, ys, zs, color=color, alpha=alpha)
    axes.plot(x2, y2, z2, color=color, alpha=alpha)
    axes.plot(x3, y3, z3, color=color, alpha=alpha)
    axes.plot(x4, y4, z4, color=color, alpha=alpha)


def draw_breast_tank(axes, bt):

    # breast box bounding box
    xmin = -float(bt.width)/2.0
    zmin = -float(bt.origin[1])  
    xmax = float(bt.width)/2.0
    zmax = float(bt.height)-float(bt.origin[1])
    ymin = 0
    ymax = bt.thickness
    bt_bbox = [xmin, xmax, ymin, ymax, zmin, zmax]
    draw_box(axes, bt_bbox, 'gray')

    #draw top extension
    ztopmax = zmin
    ztopmin = zmin - bt.topheight
    bt_top_bbox = [xmin, xmax, ymin, ymax, ztopmin, ztopmax]
    draw_box(axes, bt_top_bbox, 'gray')

    # breasktbox window bounding box
    xmin_win = -float(bt.window_width)/2.0
    xmax_win =  float(bt.window_width)/2.0
    zmin_win = -float(bt.origin[1])
    zmax_win = float(bt.window_height)-float(bt.origin[1])
    det_bbox = [xmin_win, xmax_win, ymax, ymax, zmin_win, zmax_win]
    draw_box(axes, det_bbox, 'blue')

    # source positions
    #hz = (bt.nrowsrc-1)/2
    #hx = (bt.ncolsrc-1)/2
    #sx = []
    #sy = []
    #sz = []
    #for z in range(-hz, hz):
    #    for x in range(-hx, hx):
    #        xx = x*float(bt.dsrc)
    #        zz = z*float(bt.dsrc)
    #        sx.append(xx)
    #        sz.append(zz)
    #        sy.append(ymin)
    #axes.plot(sx, sy, sz, '.', color='green') 

    # aspect ratio workaround
    bxmin = xmin 
    bxmax = xmax
    bymin = -(float(bt.width)/2.0 + float(bt.thickness)/2.0)
    bymax =  (float(bt.width)/2.0 - float(bt.thickness)/2.0)
    bzmax = float(bt.height)-float(bt.origin[1]) 
    bzmin = -float(bt.width)+bzmax
    b_bbox = [bxmin, bxmax, bymin, bymax, bzmin, bzmax]
    draw_box(axes, b_bbox, 'white', 0.0)

def draw_3dpoints(axes, pts, color):
    axes.plot(pts[:,0], pts[:,1], pts[:,2], '.', color=color)


print >>sys.stderr, 'Usage: {0} [.qm] [.msh]'.format(os.path.basename(sys.argv[0]))
    
files = []
if len(sys.argv) > 1:
    files = sys.argv[1:]
print files

fig = plt.figure()
ax = Axes3D(fig, azim=-50, elev=-158)

# breast box dimensions etc
bt = struct()
bt.topheight = 101.6
bt.height = 285
bt.width = 385
bt.thickness = 60
bt.window_width = 260.35
bt.window_height = 180.975
bt.origin = [130, 48]
bt.dsrc = 8
bt.nrowsrc = 11
bt.ncolsrc = 19
draw_breast_tank(ax, bt)

# mesh 
#print len(sys.argv)
#if (len(sys.argv) == 7):
#    mxmin = float(sys.argv[1])
#    mxmax = float(sys.argv[2])
#    mymin = float(sys.argv[3])
#    mymax = float(sys.argv[4])
#    mzmin = float(sys.argv[5])
#    mzmax = float(sys.argv[6])
#    mbbox = [mxmin, mxmax, mymin, mymax, mzmin, mzmax]
#    draw_box(ax, mbbox, 'red')
if len(files) > 0:
    import bopy.toast.utils as btu
    for f in files:
        if '.qm' in f:
            qm = btu.QM(f) 
            detpos = qm.get_detpos()
            draw_3dpoints(ax, detpos, 'red')
            srcpos = qm.get_srcpos()
            draw_3dpoints(ax, srcpos, 'blue')
        if '.msh' in f:
            b = btu.read_mesh_bbox(f)
            bbox = [b[0,0], b[1,0], b[0,1], b[1,1], b[0,2], b[1,2]]
            draw_box(ax, bbox, 'green')
            
         

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.axis('scaled')
plt.show()

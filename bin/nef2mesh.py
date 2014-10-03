#!/usr/bin/python 

import os
import sys

try:
    basefilename = sys.argv[1]
except:
    print >>sys.stderr, 'Usage: {0} <basename of node file>'.format(os.path.basename(sys.argv[0]))
    sys.exit(1)

print 'Basename:', basefilename
node_filename = basefilename + ".node"
ele_filename = basefilename + ".ele"
msh_filename = basefilename + ".msh"

n = open(node_filename, 'r')
m = open(msh_filename, 'w')
print >>m, "MeshData 5.0"
print >>m, ''

count = -1
nnodes = 0

for l in n:
    a = [v for v in l.strip().split()]
    if count > 0:
        if (int(a[4]) == 1):
            print >>m,"B[%s %s %s]R0" % (a[1], a[2], a[3])
        else:
            print >>m, "N[%s %s %s]R0" % (a[1], a[2], a[3])
        count = count - 1 
    elif count == -1:
        nnodes = int(a[0])
        print >>m, "NodeList %d 1" % (nnodes)
        count = nnodes
    else:
        print >>m, ""
        break
n.close()
e = open(ele_filename, 'r')
count = -1
nele = 0
for l in e:
    a = [v for v in l.strip().split()]
    if count > 0:
        print >>m, "c", a[1], a[2], a[3], a[4]
        count = count - 1
    elif count == -1:
        nele = int(a[0]) 
        print >>m, "ElementList", nele
        count = nele
    else:
        print >>m, '' 
        break
e.close()

print >>m, '[ParameterList]'
print >>m, 'Size {0}'.format(nnodes)
print >>m, 'Param1 MUA'
print >>m, 'Param2 KAPPA'
print >>m, 'Param3 N'
print >>m, 'Data'
for i in range(nnodes):
    print >>m, "0.003 0.4151100042 1.5"



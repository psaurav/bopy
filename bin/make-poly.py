#!/usr/bin/python

import os
import sys

#xmin = -130
#xmax = 130
#ymin = 0
#ymax = 60
#zmin = -70
#zmax = 50

try:
    xmin = float(sys.argv[1])
    xmax = float(sys.argv[2])
    ymin = float(sys.argv[3])
    ymax = float(sys.argv[4])
    zmin = float(sys.argv[5])
    zmax = float(sys.argv[6])
except:
    print >>sys.stderr, 'Usage: {0} <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> [output.poly]'.format(os.path.basename(sys.argv[0]))
    print >>sys.stderr, '(use mm and only integers)'
    sys.exit(1)

try:
    outputfile = sys.argv[7]
except:
    outputfile = "cuboid-{0}-{1}-{2}-{3}-{4}-{5}.poly".format(xmin, xmax, ymin, ymax, zmin, zmax)
print >>sys.stderr, 'Output:', outputfile

string = """# Part 1 - node list 
# node count, 3 dim, no attribute, no boundary marker
8  3  0  1
# Node index, node coordinates
1  {0} {1} {2} #0.0 0.0 0.0
2  {3} {1} {2} #1.0 0.0 0.0
3  {3} {4} {2} #1.0 1.0 0.0
4  {0} {4} {2} #0.0 1.0 0.0
5  {0} {1} {5} #0.0 0.0 1.0
6  {3} {1} {5} #1.0 0.0 1.0
7  {3} {4} {5} #1.0 1.0 1.0
8  {0} {4} {5} #0.0 1.0 1.0

# Part 2 - facet list
# facet count, no boundary marker
6  1
# facets
1  0  1      # 1 polygon, no hole, no boundary marker
4  1 2 3 4   # front
1  0  1
4  5 6 7 8   # back
1  0  1
4  1 2 6 5   # bottom
1  0  1
4  2 3 7 6   # right
1  0  1
4  3 4 8 7   # top
1  0  1
4  4 1 5 8   # left

# Part 3 - hole list
0            # no hole

# Part 4 - region list
0            # no region""".format(xmin, ymin, zmin, xmax, ymax, zmax)

f = open(outputfile, 'w')
print >>f, string
f.close()


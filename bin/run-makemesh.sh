#!/bin/bash

# all in mm
# get the thickness right (fix ymin=0, ymax=tank thickness).  play with others
#
default_config() {
cat > mesh.cfg<<DELIM
MESH_XMIN=-130
MESH_XMAX=130
MESH_YMIN=0
MESH_YMAX=${1}
MESH_ZMIN=-78
MESH_ZMAX=122
MESH_AVAL=10
DELIM
}

me=`basename $0`
cfgfile=mesh.cfg
# the file BREAST_TANK_THICKNESS_MM is taken from
pcfgfile=configure_calibration.cfg

if [ ! -e $cfgfile ]; then
    echo "${me}: config file ${cfgfile} not found"
    btt=`cat configure_calibration.cfg | egrep "^BREAST_TANK_THICKNESS_MM=" | awk -F= '{print $2}'`
    default_config $btt
    echo "${me}: config file ${cfgfile} written"
    echo "     * check file and make necessary changes"
    echo "     * Especially, check value of MESH_YMAX"
    exit
fi

declare -A keys
oldIFS="$IFS"
IFS="="
while read name value
do
    keys[$name]="$value"
done < $cfgfile
IFS="$oldIFS"

xmin=${keys['MESH_XMIN']}
xmax=${keys['MESH_XMAX']}
ymin=${keys['MESH_YMIN']}
ymax=${keys['MESH_YMAX']}
zmin=${keys['MESH_ZMIN']}
zmax=${keys['MESH_ZMAX']}
a=${keys['MESH_AVAL']}

make_poly=make-poly.py
tetgen=tetgen
nef2mesh=nef2mesh.py

polyfilename=cuboid-eh2-${a}.poly
bname=`basename ${polyfilename} .poly`

# make initial polygon
echo "***************************************************************"
echo "Making polygon with the following command"
echo "$make_poly $xmin $xmax $ymin $ymax $zmin $zmax ${polyfilename}"
$make_poly $xmin $xmax $ymin $ymax $zmin $zmax ${polyfilename} 
echo "Done making"
echo "***************************************************************"

echo ""
echo "***************************************************************"
echo "Making mesh with the following command"
echo "${tetgen} -p -a${a} ${polyfilename}"
# create mesh using tetgen
${tetgen} -p -a${a} ${polyfilename}
# .ele, .node, .face files created.
echo "Done making: .ele, .node, .face files created."
echo "***************************************************************"

# convert .node and .ele files to .msh file that Toast understands
echo ""
echo "***************************************************************"
echo "Making .msh file with the following command"
echo "$nef2mesh ${bname}.1"
$nef2mesh ${bname}.1
echo "Done making: .msh file created."
echo "***************************************************************"
#
# creates .msh file

# move files to a simpler naming convention
echo ""
echo "***************************************************************"
echo "Moving files"
mv ${bname}.1.node cuboid-eh2e-${a}.node
mv ${bname}.1.ele cuboid-eh2e-${a}.ele
mv ${bname}.1.face cuboid-eh2e-${a}.face
mv ${bname}.1.msh cuboid-eh2e-${a}.msh
echo "Done moving: Look at cuboid-eh2e-${a}.msh, the mesh file"
echo "You may check the file BBOX with the following command"
echo "  breasttanktop-mesh.py <qm-filename> <mesh-filename>"
echo "***************************************************************"

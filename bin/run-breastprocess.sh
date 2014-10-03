#!/bin/bash

# The previous configuration file
pcfgfile=configure_calibration.cfg
cfgfile=breasttoast.cfg

echo $pcfgfile

function breast_config() {
cat >${cfgfile}<<DELIM
BREAST_SHAPEFILE=20140307_cam0_proj_calib_00_reconstruction.npz
BREAST_MESH_FILE=cuboid-exthead2-10.msh
BREAST_BULK_CONC_HBO2_MICROMOLAR=29.0143
BREAST_BULK_CONC_HB_MICROMOLAR=17.6183
BREAST_BULK_PERCENTAGE_WATER=31
BREAST_BULK_PERCENTAGE_LIPID=57
BREAST_BULK_SCATTERING_PREFACTOR_A=10.1938
BREAST_BULK_SCATTERING_POWER_B=0.401569
BREAST_SPECT_WATER=segelstein81
BREAST_SPECT_LIPID=vanveen
BREAST_SPECT_HBO2HB=prahl
BREAST_CHESTWALLPOS_Z_MM=-50
BREAST_REFINDEX=1.44
HOMOG_REFINDEX=1.33
BREAST_GRID_DELTA_MM=2
DELIM
}

function get_val() {
echo `cat ${1} | egrep "^${2}=" | awk -F= '{print $2}'`
}

if [ ! -e "${cfgfile}" ]; then
    echo "${0}: Configuration file ${cfgfile} not found."
    echo "    ... Writing example file and exiting."
    breast_config
    echo "************************************************************************"
    echo "Configuration file written."
    echo "  * Edit the BREAST_SHAPEFILE"
    echo "  * Edit BREAST_MESH_FILE"
    echo "  * Edit the correct values for breast bulk properties (LBS?)"
    echo "  * then run this program again"
    echo "************************************************************************"
    exit
fi

if [ ! -e "${pcfgfile}" ]; then
    echo "${0}: Previous configuration file, ${pcfgfile} not found. Exiting."
    exit
fi


bshp=`get_val ${cfgfile} BREAST_SHAPEFILE`
mshf=`get_val ${cfgfile} BREAST_MESH_FILE`
zpos=`get_val ${cfgfile} BREAST_CHESTWALLPOS_Z_MM`

echo "************************************************************************"
echo "${0}: producing in1out0.nim file"
msh2inoutbreast.py $mshf $bshp $zpos
echo "    ... done producing in1out0.nim file"
echo "************************************************************************"

#wls=`get_val ${pcfgfile} WAVELENGTHS` 
#echo ${wls}

echo "************************************************************************"
echo "${0}: producing breast values for TOASt configuration"
breastvalues.py ${pcfgfile} ${cfgfile} 
echo "Done: look at the file"
echo "-----"
thisdir=`pwd`
echo "And this is your directory, to help set the NIM, and MESH directories:"
echo "${thisdir}"
echo "************************************************************************"


#!/bin/bash

toast_data_cfg() {
cat > toast.cfg << DELIM
PP_DATA_DIR = ..
TRG_DIR_NAME = ${1}-gsmooth-${3}
REF_DIR_NAME = ${2}-gsmooth-${3}
SINGLE_SPECT_MASK_DIR = ../masks
DATA_ROW_OFFSET = 0
DATA_COL_OFFSET = 0
DATA_ROW_STRIDE = 10
DATA_COL_STRIDE = 10
MUSP_CORR_SOURCE_POSITION = ${4}
LOESS_SPAN = 0.02
DATA_FILE_STRING = {0}_{1}_wl{2}_s{3}.npy
SRC_ROW_OFFSET = 0
SRC_COL_OFFSET = 0
SRC_ROW_STRIDE = 2
SRC_COL_STRIDE = 2
SRC_NROW = 11
SRC_NCOL = 19
DELIM
}

# Initialize our own variables:
target=""
reference=""
cfgfile="conf_calib_pomxrw_g3fit_params.cfg"

target=`cat $cfgfile |egrep "^TRG_NAME=" | awk -F= '{print $2}'`
reference=`cat $cfgfile |egrep "^REF_NAME=" | awk -F= '{print $2}'`

sigma=`cat $cfgfile |egrep "^SIGMA=" |awk -F= '{print $2}'`
echo "sigma: $sigma"
#
# Get the mean musp of all musp's.
mean_musp=`fgrep CALCULATED_MUSP_MM $cfgfile | awk -F= '{ total += $2; count++ } END { print total/count }'`

#
# Use mean musp to calculate position of isotropic diffuse source.
toast_data_cfg $target $reference $sigma $mean_musp

echo " "
echo "*****************************************************************************************"
echo "Created configuration file toast.cfg.  Check the file and edit according to your needs."
echo "Next run-maketoastdata.sh" 
echo "*****************************************************************************************"


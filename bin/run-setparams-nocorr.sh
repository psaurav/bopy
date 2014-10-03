#!/bin/bash

default_params_cfg() {
cat >> $1 <<DELIM
SIGMA=8.0
MASK_TOPROWS=14
MASK_BOTTOMROWS=5
MASK_LEFTCOLS=5
MASK_RIGHTCOLS=7
MASK_PICKOFF_ROW=100
MASK_PICKOFF_COL=59
SIGNAL_MAX_CUTOFF=0.03
DELIM
}

function get_val() {
echo `cat ${1} | egrep "^${2}=" | awk -F= '{print $2}'`
}


target=""
reference=""
cfgfile="conf_calib_pomxrw_g3fit.cfg"

echo "cfgfile='$cfgfile'"

if [ ! -e "${cfgfile}" ]
then
    echo "Error: Configuration file ${cfgfile} does not exist"
    exit 1
else
    cfgfile=`readlink -f "${cfgfile}"`
fi
target=`cat $cfgfile |egrep "^TRG_NAME=" | awk -F= '{print $2}'`
reference=`cat $cfgfile |egrep "^REF_NAME=" | awk -F= '{print $2}'`

ocfgfile="conf_calib_pomxrw_g3fit_params.cfg"
cat $cfgfile > $ocfgfile
default_params_cfg $ocfgfile
echo "***************************************************************************"
echo "Wrote config file: $ocfgfile"
echo "  * Make adjustments as needed"
echo "***************************************************************************"

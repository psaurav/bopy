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
CALIB_PIXELS=116 189
DKINDS=A phi DC
CALIB_SRCS=0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 209
PICKOFF_PIXELS=14 28 3 16
WLS_CALIB_AMPLITUDE=785 830
WLS_PICKOFF_AMPLITUDE=660 690 808
WLS_PICKOFF_PHASE=660 690 785 808 830
DELIM
}

target=""
reference=""
cfgfile="conf_calib_pomxrw_g3fit_pomncl_poreg.cfg"

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

ocfgfile="conf_calib_pomxrw_g3fit_pomncl_poreg_params.cfg"
cat $cfgfile > $ocfgfile
default_params_cfg $ocfgfile
echo "***************************************************************************"
echo "Wrote config file: $ocfgfile"
echo "  * Make adjustments as needed"
echo "***************************************************************************"

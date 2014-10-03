#!/bin/bash

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
target=""
reference=""
cfgfile="conf_calib_pomxrw_g3fit.cfg"

#echo "Name of previous cfgfile='$cfgfile'"

if [ ! -e "${cfgfile}" ]
then
    echo "Error: Configuration file ${cfgfile} does not exist"
    exit 1
else
    cfgfile=`readlink -f "${cfgfile}"`
fi

target=`cat $cfgfile |egrep "^TRG_NAME=" | awk -F= '{print $2}'`
reference=`cat $cfgfile |egrep "^REF_NAME=" | awk -F= '{print $2}'`

WLSS=`cat $cfgfile |egrep "^WAVELENGTHS" | awk -F= '{print $2}'`
WLS=(`echo $WLSS | tr " " "\n"`)

STUDY_NAME=`cat $cfgfile |egrep "^STUDY_NAME" | awk -F= '{print $2}'`
DATA_DIRECTORY=`cat $cfgfile |egrep "^DATA_DIRECTORY" | awk -F= '{print $2}'`
DATA_DIRECTORY=`readlink -f ${DATA_DIRECTORY}`
ddir=${DATA_DIRECTORY}/${STUDY_NAME}

#
# Fit the plain sources
#
for meas in ${target} ${reference}
do
    if [ ! -d "${meas}" ]; then
        mkdir ${meas} 
    fi
    cd ${meas} 
    for w in "${WLS[@]}"
    do
        echo "$w"
        (gen3-tsfit.py $ddir ${meas} $cfgfile $w src 2>&1 |tee fit-src-${w}.log &)
    done
    cd ..
done
exit

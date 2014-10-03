#!/bin/bash


# Initialize our own variables:
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

WLSS=`cat $cfgfile |egrep "^WAVELENGTHS" | awk -F= '{print $2}'`
WLS=(`echo $WLSS | tr " " "\n"`)

STUDY_NAME=`cat $cfgfile |egrep "^STUDY_NAME" | awk -F= '{print $2}'`
DATA_DIRECTORY=`cat $cfgfile |egrep "^DATA_DIRECTORY" | awk -F= '{print $2}'`
DATA_DIRECTORY=`readlink -f ${DATA_DIRECTORY}`
ddir=${DATA_DIRECTORY}/${STUDY_NAME}

#
# Now get the From column
#
ocfgfile=conf_calib_pomxrw_g3fit_pomncl.cfg
gen3-pomincol.py $cfgfile $ocfgfile

cfgfile="$ocfgfile"
cfgfile=`readlink -f "${cfgfile}"`
#
# Now the calib fits
#
for meas in ${target} ${reference}
do
    if [ ! -d "calibsrc-${meas}" ]; then
        mkdir calibsrc-${meas} 
    fi
    cd calibsrc-${meas} 
    for w in "${WLS[@]}"
    do
        echo "$w"
        (gen3-tsfit.py $ddir ${meas} $cfgfile $w calib 2>&1 |tee fit-src-${w}.log &)
    done
    cd ..
done
exit

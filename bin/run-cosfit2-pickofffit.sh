#!/bin/bash

# Initialize our own variables:
target=""
reference=""
cfgfile="conf_calib_pomxrw_g3fit_pomncl.cfg"

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
# Now get the pickoff region
#
ocfgfile=conf_calib_pomxrw_g3fit_pomncl_poreg.cfg
imgfile=`ls -1 ${reference}/*_phi.npy | fgrep "_wl${WLS[0]}_" | fgrep "_s199_"`
gen3-poregion.py $cfgfile $ocfgfile ${imgfile}
cfgfile="$ocfgfile"
cfgfile=`readlink -f "${cfgfile}"`

#
# Now the pickoff fits
#
for meas in ${target} ${reference}
do
    if [ ! -d "pickoff-${meas}" ]; then
        mkdir pickoff-${meas} 
    fi
    cd pickoff-${meas} 
    for w in "${WLS[@]}"
    do
        echo "$w"
        (gen3-tsfit.py $ddir ${meas} $cfgfile $w pickoff 2>&1 |tee fit-src-${w}.log &)
    done
    cd ..
done
wait


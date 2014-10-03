#!/bin/bash

default_gen3fit_cfg() {
cat > $1 <<DELIM
GEN3FIT_SAMPLING_INTERVAL_SEC=0.1
GEN3FIT_BEAT_FREQ_HZ=1.0
GEN3FIT_N_DATAPOINTS=17
GEN3FIT_ABS_ERROR=0.0
GEN3FIT_REL_ERROR=1.0e-5
GEN3FIT_MAX_ITERATIONS=200
GEN3FIT_OUTPUTS=A phi DC
GEN3FIT_FILENAME_ENCODING=GEN3
GEN3FIT_DROP_FIRST_FRAMES=0
GEN3FIT_CALL=aphidc_flat_andro_udrot2
DELIM
}

usage() { echo "Usage: $0 -d ddir" 1>&2; exit 1; }

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
ddir=""
cfgfile=""

while getopts ":d:" opt; do
    #echo "$opt $OPTARG"
    case "$opt" in
    d)  ddir=$OPTARG
        ;;
    *) 
        usage
        ;;
    esac
done

shift $((OPTIND-1))

if [ -z "${ddir}" ]; then
    usage
fi

echo "ddir='$ddir'"

if [ ! -d "${ddir}" ]
then
    echo "Error: Directory ${ddir} does not exist"
    exit 1
else
    ddir=`readlink -f "${ddir}"`
fi

root=`basename ${ddir}`
datadir=`dirname ${ddir}`
cp ${ddir}/*.cfg .
cp ${ddir}/*.txt .

echo "************************************************************"
echo "Now performing init"
gen3-init.py
echo "************************************************************"
echo "************************************************************"
echo "Now performing calibration"
echo "Do following in sequence"
echo "  * Middle button - identify middle source"
echo "  * Right button - identify top-left source"
echo "  * Right button - identify top-right source"
echo "  * 'w' key - save to file"
gen3-calib.py ${root}-homog.cfg
echo "Done init"
echo "************************************************************"
echo ""
echo "************************************************************"
echo "Now perform max row pickoff"
echo "  * Do this irrespective of whether there exists a pickoff signal"
echo "  * Middle button - click to so that the hirozontal line is between"
echo "                    the pickoff and lowest row of sources" 
echo "  * 'w' key - save to file"
gen3-pomaxrow.py configure_calibration.cfg configure_calibration-pomaxrow.cfg
echo "Done max row pickoff"
echo "************************************************************"
echo ""
echo "************************************************************"
echo "Writing final configuration files"
default_gen3fit_cfg conf_calib_pomxrw_g3fit.cfg
cat configure_calibration-pomaxrow.cfg >>conf_calib_pomxrw_g3fit.cfg
echo "Done writing final config file: conf_calib_pomxrw_g3fit.cfg"
echo "  * Check the configuration parameters in the file"
echo "************************************************************"
exit

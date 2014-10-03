#!/bin/bash

# Initialize our own variables:
target=""
reference=""
cfgfile="conf_calib_pomxrw_g3fit_params.cfg"

function get_val() {
echo `cat ${1} | egrep "^${2}=" | awk -F= '{print $2}'`
}

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
root=`get_val $cfgfile ROOT`
fname=${reference}/${root}_${reference}_wl660_s99_A.npy
echo ${fname}
if [ ! -e "${fname}" ]; then
    testf=`ls -1 ${reference}/*_wl660_s99_A.npy | head -n 1`
    echo "*********"
    echo ""
    echo "  ERROR: File $fname not found"
    echo "         listed file is $testf"
    echo "         possibly ROOT parameter is wrong.  check"
    echo "         if so, correct in file ${cfgfile} and repeat"
    echo ""
    echo "*********"
    exit
fi

#
# Gaussian smoothing
#
sigma=`cat $cfgfile | egrep "^SIGMA=" | awk -F= '{print $2}'`
echo "Gausss smoothing target"
gen3-togauss.py $target $sigma
echo "Gausss smoothing reference"
gen3-togauss.py $reference $sigma
echo "Done Gausss smoothing"

#
# masks
#
echo ""
echo "************************************************************"
echo "Creating masks..."
fname=${reference}/${root}_${reference}_wl660_s99_A.npy
nrow=`get_val $cfgfile RF_DATA_NROW`
ncol=`get_val $cfgfile RF_DATA_NCOL`
msk_tr=`get_val $cfgfile MASK_TOPROWS`
msk_br=`get_val $cfgfile MASK_BOTTOMROWS`
msk_lc=`get_val $cfgfile MASK_LEFTCOLS`
msk_rc=`get_val $cfgfile MASK_RIGHTCOLS`
rmin=${msk_tr}
cmin=${msk_lc}
rmax=`echo "${nrow} - ${msk_br}" |bc`
cmax=`echo "${ncol} - ${msk_rc}" |bc`
rpo=`get_val $cfgfile MASK_PICKOFF_ROW`
cpo=`get_val $cfgfile MASK_PICKOFF_COL`
wls=`cat $cfgfile | egrep "^WAVELENGTHS=" | awk -F= '{print $2}'`
ddir=.
refdir=${reference}-gsmooth-${sigma}
trgdir=${target}-gsmooth-${sigma}
ref=${reference}
trg=${target}
cutoff=`get_val $cfgfile SIGNAL_MAX_CUTOFF`
echo "rmin = ${rmin}"
echo "cmin = ${cmin}"
echo "rmax = ${rmax}"
echo "cmax = ${cmax}"
echo "rpo = ${rpo}"
echo "cpo = ${cpo}"
echo "cutoff = ${cutoff}"
echo "Running..."
echo "create_genmask2_masks.py $fname $rmin $rmax $cmin $cmax $rpo $cpo \"${wls}\" $ddir $refdir $trgdir $root $ref $trg $cutoff"
create_genmask2_masks.py $fname $rmin $rmax $cmin $cmax $rpo $cpo "${wls}" $ddir $refdir $trgdir $root $ref $trg $cutoff
echo "Done creating masks..."
echo "************************************************************"

echo "Doing phaseunwrap..."
tdir="${target}-gsmooth-${sigma}"
troot="${root}_${target}"
rdir="../${reference}-gsmooth-${sigma}"
rroot="${root}_${reference}"
cd $tdir
gen3-phaseunwrap.py $tdir $troot $rdir $rroot
cd ..

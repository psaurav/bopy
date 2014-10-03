#!/bin/bash

cfgfile1=conf_calib_pomxrw_g3fit_params.cfg
cfgfile2=toast.cfg

echo $cfgfile1
echo $cfgfile2

mkdir toastdata
cd toastdata
gen3-toastdata.py $cfgfile1 $cfgfile2
cd ..

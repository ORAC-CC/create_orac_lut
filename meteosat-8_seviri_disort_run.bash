#!/bin/csh
source initpath
# monochromatic calculation

# liquid-water section
setenv CREATE_ORAC_LUT_DRIVER 'input_files/driver/meteosat-8_seviri_cloud.driver'
echo "making ..." $CREATE_ORAC_LUT_DRIVER
idl -e "create_orac_lut_wrapper,srf_quad=1,mmfile='liquid-water_old.mm',lutfile='liquid-water-cloud.lut',version=7,reuse_scat=1"

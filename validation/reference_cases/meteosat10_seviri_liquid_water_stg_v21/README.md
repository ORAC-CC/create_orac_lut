# Meteosat-10 SEVIRI liquid-water STG V2

This is the initial Python reproduction target. The numerical reference remains
in the authorised read-only archive; it is not copied into this repository.

The product is a current V2 NetCDF4 LUT with all 11 SEVIRI channels, solar,
mixed, and thermal operators, and the simple `liquid-water_stg.mm` particle
model. It is a clean first test of the reader, dimensions, metadata, spectral
channel handling, and output operator layout before implementing particle optics
or DISORT.

The archived output driver and current IDL writer are the configuration evidence.
The exact historical cloud shell script is not retained in the current working
tree; `meteosat-10_seviri_run` currently contains aerosol commands only. This
case must therefore preserve the archived driver, source revision marker, and
the documented uncertainty rather than inventing a command line.

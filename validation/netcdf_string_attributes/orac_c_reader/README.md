# ORAC C-reader harness

`main.c` wraps ORAC's unmodified `common/nc_get_string_att.c` (the routine
behind `ncdf_get_string_att` in `src/read_sad_lut.F90`) so that the five LUT
axis `spacing` attributes can be read with exactly the code path that failed
in production.  `build_and_run.sh ORAC_CHECKOUT FILE...` compiles it against
the project Conda environment's NetCDF library and prints one line per axis.
The ORAC source is used in place from a read-only checkout and is not copied
here.  Results of the 2026-10-09 campaign are in `results_*.tsv`.

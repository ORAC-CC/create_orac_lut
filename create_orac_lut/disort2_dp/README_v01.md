# Double-precision DISORT for the ORAC LUT generators

This directory holds the double-precision DISORT used by default by both
LUT generators since 2026-10-10:

| File | Role |
|---|---|
| `DISORTtoIDL_DP.c` | IDL interface: procedure `DISORT_DP` (the same 44 arguments as `DISORT`) and function `DISORT_DP_INFO()` |
| `dp_support.f` | `DPREALBYTES` (kernel REAL size; the DLM refuses to load unless 8) and a fatal `BDREF` stub (the BRDF path is never used by ORAC and its C helper is single precision) |
| `DISORTtoIDL_DP.dlm`, `exports.map` | DLM descriptor; linker map exporting only `IDL_Load`, so the kernel cannot collide with the historical `../disort2/src/DISORTtoIDL.so` loaded in the same IDL session |
| `build_dlm.sh` | reproducible build; writes `DISORTtoIDL_DP.so`, `libdisort_dp.so` and `BUILD_INFO.txt` |
| `DISORTtoIDL_DP.so` | the IDL DLM (build product) |
| `libdisort_dp.so` | the same kernel objects as a plain library exporting `disort_`, loaded by the Python generator (build product) |

The kernel is the unmodified production DISORT source of `../disort2/src`
(`DISORTfunctions.f`, `PRTFIN.f`, `disort/ErrPack.f`, `disort/LINPAK.f`,
`disort/RDI1MACH.f`) compiled with ifort 2021.6.0 (module
`intel-compilers/2022`), `-O2 -fp-model precise -fpic -real-size 64
-double-size 64 -integer-size 32`.  Its object files are bitwise identical to
those validated through IDL in
`validation/v21_v22_terra_modis/REPORT_idl_double_precision.md`.

## Selection

Both generators read the environment variable `ORAC_LUT_DISORT_PRECISION`:

| Value | IDL (`create_orac_lut/call_disort.pro`) | Python (`oraclut.radiative_transfer.legacy_disort`) |
|---|---|---|
| unset or `double` | `DISORT_DP` (this directory) | `libdisort_dp.so` (this directory) |
| `single` | the historical `DISORT` DLM, the original statement unchanged | the historical single-precision kernel |

Anything else is an error.  IDL prints `call_disort: DISORT precision = ...`
once per session; the Python provenance sidecar records `disort_precision`.
`GETMOM` and `PLKAVG` always use the historical single-precision library in
both generators.  The stream count stays 60 (`setup_disort`).

## Rebuild

```shell
bash create_orac_lut/disort2_dp/build_dlm.sh
```

Check the result against the release manifests
(`python scripts/check_release_manifest.py releases/idl_lut_v21_double.release`
and the launcher's `--check 22` with `ORAC_LUT_RELEASE_MANIFEST=lut_v22_double.release`).

Do not SourceVault-version `DISORTtoIDL_DP.dlm` in place: a second `.dlm`
declaring the same module in the IDL DLM path would conflict at load time.

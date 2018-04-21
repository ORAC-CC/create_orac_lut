Northrop Grumman October 2011 RSR Release.

This data is provided by Northrop Grumman under the NPOESS contract (F04701-02-0502, NPOESS A&O).

Note: Flight model 1 ("F1") is the VIIRS sensor model flying on the Suomi NPP platform.

NG Fused RSRs for F1 VisNIR M bands:
Fused RSR is derived by combing spacecraft-level T-SIRCUS test data and sensor-level Tvac data. 
Fused RSR adopts most of all high quality T-SIRCUS data. Tvac RSR are inserted where no high quality T-SIRCUS exists. 
T-SIRCUS RSR is based on GT_F1_SC_RSR_release1.0 data files. 
Tvac Effective RSR is based on polarization-corrected effective RSR files.
Bias in Tvac RSR are removed based on comparison of Tvac and T-SIRCUS RSR in regions where 
high quality data exists for both.

NOTES ON the fused RSR:

1)  Fused RSR adopts all high quality T-SIRCUS data identified in GT data files. 
    However, some low quality are retained if standard deviation of detector response normalized to
    detector averaged response is less than average uncertainty of individual detector. 
    Additional SNR threshold flag of 2 is used to discard low quality if it pass the dispersion test.
2)  Detector dispersion test is used to determine whether to use per detector RSR valises or an
    averaged RSR value for to both T-SIRCUS and Tvac data. If the test is true, detector average value is used 
    for all individual detectors. If not, individual detector values are retained. 
3)  Electronic crosstalk flag in T-SIRCUS data is from Government released data files. Similar electronic crosstalk flag
    is derived for Tvac data. It is recommended that electronic crosstalk driven response be removed from RSR 
    if an electronic crosstalk correction applied to VIIRS observations before being calibrated to SDR.
4)  Most of angular resolved scatter (ARS) and fluorescence driven response in Tvac RSR exhibit positive bias in 
    comparison with T-SIRCUS data.  The bias was removed based on comparison of Tvac and T-SIRCUS RSR in regions 
    where high quality data exists for both, however, excluding electronic crosstalk and color-glass transmission leaks.
5)  Fused RSR has minimal polarization bias caused by test equipment artifacts, since T-SIRCUS source is unpolarized, 
    whereas polarization bias in Tvac RSR has been corrected through a correction algorithm.
6)  SIRCUS RSR uncertainty is derived by RSS of errors in dn, top and bottom radiance monitor measurements, 
    and bias in radiance monitor calibration.  The bias in radiance monitor calibration is found to be low, 
    ~1.2%, based on data presented by Tom Schwarting of NICST.

NG Effective RSRs for F1 VisNIR I bands:
Tvac Effective RSR is based on polarization-corrected effective RSR files.
I band RSRs were not updated using T-SIRCUS data due to the small differnces between the Tvac and T-SIRCUS data for these
bands and the fact that the OCC algorithm does not use these bands.
Band I2 shows a small in band wavelength shift relative to M7.  This feature was comfirmed by T-SIRCUS testing.
Tvac Effective RSR show greater detector to detector variation than T-SIRCUS data; however, since only detector average
RSRs are used for these bands, this will not impact any of the VIIRS data products.

NG RSRs for F1 SMWIR and LWIR bands:
These are true RSR. M9 RSR has been updated for water-vapor absorption
Notes on the SMWIR and LWIR RSRs:

1)  SNR thresholds are used to identify useful lower and upper wavelength limits. Outside of these limits are mostly 
    low quality data that are largely noise driven.
2)  NG recommends the use of the RSR data within the limits to define the spectral characteristics of the VIIRS F1 
    sensor for the SMWIR and LWIR bands. NG does not recommend the use of low quality data outside of the limits.
3)  Recomended wavelength limits and SNR threshold are given below:

   Band  lower_limit,micron  upper_limit, micron  SNR threshold
   M8           1.15           2.30               0.4
   M9           1.28           1.77               0.4
   M10          1.35           2.30               0.4
   M11          1.29           2.46               0.4
   M12          3.20           4.25               0.4
   M13          3.60           4.56               0.4
   M14          8.20           9.00               0.8
   M15          7.50          12.50               0.5
   M16A         8.0           13.00               0.8
   M16B         8.0           13.00               0.8
   I3           1.40           2.30               0.4
   I4           3.40           4.10               0.4
   I5          10.00          13.00               0.4

4) No good measurement of ZnSe window spectral transmission used during F1 test program.  The window will be 
   measured during F2 test program possibly resulting in changes to the RSRs that were measured through this window.

DATA FORMAT:
Column 1  - wavelength in nm
Column 2  - detector average RSR  (aka "band averaged")
Column 3 to 18(M) or 34(I)  - per detector RSR for detector 1 to detector 16 (32) in sensor order

CONTACT INFO:
Mau-Song Chou: mau-song.chou@ngc.com
20-Oct-11

Copyright 2011, Northrop Grumman Systems Corporation. All Rights Reserved.


****************************************** Government Team Notes **************************************************
NOTE 1: Because the band I5 RSR upper limit wavelength cutoff is at a response greater than 0.01 in the I5 band 
averaged RSR, the Govt Team has decided to extend the cutoff from 13.00um out to 13.036um so as to include the 0.01
response wavelength. 

NOTE 2: The Oct 2011 band I3 RSR contained anomalous low response in a portion of the out-of-band region.  The
Govt Team replaced this response using high quality response from a Nov 2011 I3 RSR update provided by NG.

NOTE 3: The band averaged RSR of all bands were normalized to a peak response of 1.0 by the Govt Team.
 
NOTE 4: The Govt Team applied the recommended NG wavelength limits to create a user-ready NG RSR product containing 
band averaged and supporting detector level RSR.  The filenames containing this modified NG RSR product are:

3243240 Jan 13 21:30 NPP_VIIRS_NG_RSR_I1_Effective_Oct2011f.dat
3243240 Jan 13 21:31 NPP_VIIRS_NG_RSR_I2_Effective_Oct2011f.dat
 343560 Feb  1 18:43 NPP_VIIRS_NG_RSR_I3_filtered_Oct2011f_Mod.dat
 308440 Jan 13 21:22 NPP_VIIRS_NG_RSR_I4_filtered_Oct2011f.dat
1336280 Jan 13 21:22 NPP_VIIRS_NG_RSR_I5_filtered_Oct2011f.dat
1443272 Jan 13 21:31 NPP_VIIRS_NG_RSR_M1_Fused_Oct2011f.dat
1549992 Jan 13 21:32 NPP_VIIRS_NG_RSR_M2_Fused_Oct2011f.dat
1508232 Jan 13 21:33 NPP_VIIRS_NG_RSR_M3_Fused_Oct2011f.dat
1538392 Jan 13 21:33 NPP_VIIRS_NG_RSR_M4_Fused_Oct2011f.dat
1547672 Jan 13 21:33 NPP_VIIRS_NG_RSR_M5_Fused_Oct2011f.dat
1473432 Jan 13 21:34 NPP_VIIRS_NG_RSR_M6_Fused_Oct2011f.dat
1556952 Jan 13 21:34 NPP_VIIRS_NG_RSR_M7_Fused_Oct2011f.dat
 267032 Jan 13 21:23 NPP_VIIRS_NG_RSR_M8_filtered_Oct2011f.dat
 113912 Jan 13 21:23 NPP_VIIRS_NG_RSR_M9_filtered_Oct2011f.dat
 220632 Jan 13 21:24 NPP_VIIRS_NG_RSR_M10_filtered_Oct2011f.dat
 271672 Jan 13 21:24 NPP_VIIRS_NG_RSR_M11_filtered_Oct2011f.dat
 243832 Jan 13 21:25 NPP_VIIRS_NG_RSR_M12_filtered_Oct2011f.dat
 222952 Jan 13 21:25 NPP_VIIRS_NG_RSR_M13_filtered_Oct2011f.dat
 185832 Jan 13 21:26 NPP_VIIRS_NG_RSR_M14_filtered_Oct2011f.dat
1160232 Jan 13 21:26 NPP_VIIRS_NG_RSR_M15_filtered_Oct2011f.dat
1160232 Jan 13 21:26 NPP_VIIRS_NG_RSR_M16A_filtered_Oct2011f.dat
1160232 Jan 13 21:27 NPP_VIIRS_NG_RSR_M16B_filtered_Oct2011f.dat

NOTE 5: From the data files listed in NOTE 3, the Govt Team has created a user-ready NG RSR product 
containing band averaged RSR only (no detector level RSR).  The filenames are:

 184275 Jan 25 22:17 NPP_VIIRS_NG_RSR_I1_Effective_Oct2011f_BA.dat
 184275 Jan 25 22:18 NPP_VIIRS_NG_RSR_I2_Effective_Oct2011f_BA.dat
  22525 Jan 25 22:18 NPP_VIIRS_NG_RSR_I3_filtered_Oct2011f_BA.dat
  19632 Feb  1 18:49 NPP_VIIRS_NG_RSR_I3_filtered_Oct2011f_BA_Mod.dat
  17525 Jan 25 22:18 NPP_VIIRS_NG_RSR_I4_filtered_Oct2011f_BA.dat
  75925 Jan 25 22:18 NPP_VIIRS_NG_RSR_I5_filtered_Oct2011f_BA.dat
  23775 Jan 25 22:18 NPP_VIIRS_NG_RSR_M10_filtered_Oct2011f_BA.dat
  29275 Jan 25 22:18 NPP_VIIRS_NG_RSR_M11_filtered_Oct2011f_BA.dat
  26275 Jan 25 22:19 NPP_VIIRS_NG_RSR_M12_filtered_Oct2011f_BA.dat
  24025 Jan 25 22:19 NPP_VIIRS_NG_RSR_M13_filtered_Oct2011f_BA.dat
  20025 Jan 25 22:19 NPP_VIIRS_NG_RSR_M14_filtered_Oct2011f_BA.dat
 125025 Jan 25 22:19 NPP_VIIRS_NG_RSR_M15_filtered_Oct2011f_BA.dat
 125025 Jan 25 22:19 NPP_VIIRS_NG_RSR_M16A_filtered_Oct2011f_BA.dat
 125025 Jan 25 22:19 NPP_VIIRS_NG_RSR_M16B_filtered_Oct2011f_BA.dat
 155525 Jan 25 22:20 NPP_VIIRS_NG_RSR_M1_Fused_Oct2011f_BA.dat
 167025 Jan 25 22:20 NPP_VIIRS_NG_RSR_M2_Fused_Oct2011f_BA.dat
 162525 Jan 25 22:20 NPP_VIIRS_NG_RSR_M3_Fused_Oct2011f_BA.dat
 165775 Jan 25 22:20 NPP_VIIRS_NG_RSR_M4_Fused_Oct2011f_BA.dat
 166775 Jan 25 22:20 NPP_VIIRS_NG_RSR_M5_Fused_Oct2011f_BA.dat
 158775 Jan 25 22:20 NPP_VIIRS_NG_RSR_M6_Fused_Oct2011f_BA.dat
 167775 Jan 25 22:20 NPP_VIIRS_NG_RSR_M7_Fused_Oct2011f_BA.dat
  28775 Jan 25 22:21 NPP_VIIRS_NG_RSR_M8_filtered_Oct2011f_BA.dat
  12275 Jan 25 22:21 NPP_VIIRS_NG_RSR_M9_filtered_Oct2011f_BA.dat

NOTE 6: The NG October 2011 band averaged RSR, with modifications of Notes 1-5, have been deemed
acceptable by the Govt Team to represent the NPP VIIRS At-Launch RSR. Identical band averaged 
RSR are contained in column 2 of the data files listed in Note 4 and Note 5.

Chris Moeller, Univ. Wisconsin
chrism@ssec.wisc.edu
Feb 1, 2012

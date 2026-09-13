#FC = ifort
FC = gfortran
clean:
		rm -f OBJECTS/*.o mod3p5.exe
tidy:
		rm -f OBJECTS/*.o
mod3p5.exe:	OBJECTS/Modtrn.o OBJECTS/ab.o     OBJECTS/abcd.o   \
		OBJECTS/abcdta.o OBJECTS/aerext.o OBJECTS/aernsm.o \
		OBJECTS/aerprf.o OBJECTS/aertmp.o OBJECTS/aitk.o   \
		OBJECTS/albtrn.o OBJECTS/altrin.o OBJECTS/asymtx.o \
		OBJECTS/atmcon.o OBJECTS/bbfn.o   OBJECTS/betabs.o \
		OBJECTS/bfh2o.o  OBJECTS/bmdata.o OBJECTS/bmerfu.o \
		OBJECTS/bmflux.o OBJECTS/bmload.o OBJECTS/bmod.o   \
		OBJECTS/bmtran.o OBJECTS/bmwids.o OBJECTS/bmxld.o  \
		OBJECTS/bo2c.o   OBJECTS/bo3hh0.o OBJECTS/bo3hh1.o \
		OBJECTS/bo3hh2.o OBJECTS/bs.o     OBJECTS/c4d.o    \
		OBJECTS/c4dta.o  OBJECTS/c6dta.o  OBJECTS/check.o  \
		OBJECTS/chekin.o OBJECTS/cirr18.o OBJECTS/cirrus.o \
		OBJECTS/ckd.o    OBJECTS/cldprf.o OBJECTS/cmpint.o \
		OBJECTS/compar.o OBJECTS/convrt.o OBJECTS/cph2o.o  \
		OBJECTS/cpo3.o   OBJECTS/cptrcg.o OBJECTS/cpumix.o \
		OBJECTS/crdriv.o OBJECTS/crmerg.o OBJECTS/crprof.o \
		OBJECTS/crspec.o OBJECTS/cruprf.o OBJECTS/cruspc.o \
		OBJECTS/cxdta.o  OBJECTS/d1mach.o OBJECTS/dasum.o  \
		OBJECTS/daxpy.o  OBJECTS/dbltx.o  OBJECTS/dchdc.o  \
		OBJECTS/ddot.o   OBJECTS/debye.o  OBJECTS/defalt.o \
		OBJECTS/del.o    OBJECTS/desatt.o OBJECTS/dgbco.o  \
		OBJECTS/dgbfa.o  OBJECTS/dgbsl.o  OBJECTS/dgeco.o  \
		OBJECTS/dgefa.o  OBJECTS/dgesl.o  OBJECTS/disbbf.o \
		OBJECTS/disort.o OBJECTS/dop.o    OBJECTS/dpandx.o \
		OBJECTS/dpexnt.o OBJECTS/dpfill.o OBJECTS/dpfish.o \
		OBJECTS/dpfnmn.o OBJECTS/dplayr.o OBJECTS/dprarf.o \
		OBJECTS/dprfpa.o OBJECTS/dpscht.o OBJECTS/dref.o   \
		OBJECTS/driver.o OBJECTS/dscal.o  OBJECTS/dstdta.o \
		OBJECTS/dswap.o  OBJECTS/errmsg.o \
		OBJECTS/exabin.o OBJECTS/expint.o OBJECTS/extdta.o \
		OBJECTS/fdbeta.o OBJECTS/flayz.o  OBJECTS/fluxes.o \
		OBJECTS/flxadd.o OBJECTS/fndpth.o OBJECTS/frn296.o \
		OBJECTS/frq5dt.o OBJECTS/ftrang.o \
		OBJECTS/gamfog.o OBJECTS/geo.o    \
		OBJECTS/geoinp.o OBJECTS/gtstrt.o OBJECTS/gmrain.o \
		OBJECTS/h2src.o  OBJECTS/hengns.o OBJECTS/hertda.o \
		OBJECTS/hno3.o   OBJECTS/idamax.o OBJECTS/indx.o   \
		OBJECTS/interp.o OBJECTS/irfxn.o  OBJECTS/jou.o    \
		OBJECTS/lay5dt.o OBJECTS/laycld.o OBJECTS/layvsa.o \
		OBJECTS/lepoly.o OBJECTS/locate.o OBJECTS/loop.o   \
		OBJECTS/mapms.o  OBJECTS/mardta.o OBJECTS/marine.o \
		OBJECTS/mdta.o   OBJECTS/mlatmb.o OBJECTS/msrad.o  \
		OBJECTS/mssolr.o OBJECTS/msthml.o OBJECTS/newh2.o  \
		OBJECTS/newsrc.o OBJECTS/no2xs.o  OBJECTS/nunit.o  \
		OBJECTS/o2cont.o OBJECTS/o2int.o  OBJECTS/o3chap.o \
		OBJECTS/o3hht0.o OBJECTS/o3hht1.o OBJECTS/o3hht2.o \
		OBJECTS/o3int.o  OBJECTS/o3uv.o   OBJECTS/o3uvfd.o \
		OBJECTS/pf.o     OBJECTS/phasef.o OBJECTS/phsdta.o \
		OBJECTS/plkavg.o OBJECTS/praltr.o OBJECTS/pravin.o \
		OBJECTS/prfdta.o OBJECTS/prtinp.o OBJECTS/prtint.o \
		OBJECTS/psi.o    OBJECTS/psieca.o OBJECTS/pslct.o  \
		OBJECTS/qgausn.o OBJECTS/rain.o   OBJECTS/random.o \
		OBJECTS/ranfun.o OBJECTS/ranset.o OBJECTS/ratio.o  \
		OBJECTS/rdexa.o  OBJECTS/rdnsm.o  OBJECTS/rdsun.o  \
		OBJECTS/reduce.o OBJECTS/rhoeps.o OBJECTS/rnscat.o \
		OBJECTS/rtbis.o  OBJECTS/schrun.o OBJECTS/sctang.o \
		OBJECTS/setdis.o OBJECTS/setmtx.o OBJECTS/sf260.o  \
		OBJECTS/sf296.o  OBJECTS/shade.o  OBJECTS/shumg.o  \
		OBJECTS/sint.o   OBJECTS/slf260.o OBJECTS/slf296.o \
		OBJECTS/slftst.o OBJECTS/smgeo.o  OBJECTS/smprep.o \
		OBJECTS/so2xs.o  OBJECTS/soleig.o OBJECTS/solve0.o \
		OBJECTS/solve1.o OBJECTS/solzen.o OBJECTS/source.o \
		OBJECTS/spaltr.o OBJECTS/ssgeo.o  OBJECTS/ssrad.o  \
		OBJECTS/stdmdl.o OBJECTS/subsol.o OBJECTS/sunbd.o  \
		OBJECTS/surfac.o OBJECTS/svsola.o OBJECTS/tab.o    \
		OBJECTS/tanht.o  OBJECTS/terpev.o OBJECTS/terpso.o \
		OBJECTS/title.o  OBJECTS/tnrain.o OBJECTS/trans.o  \
		OBJECTS/trlay.o  OBJECTS/tstbad.o OBJECTS/upbeam.o \
		OBJECTS/upisot.o OBJECTS/usrint.o OBJECTS/vsa.o    \
		OBJECTS/vsansm.o OBJECTS/watvap.o OBJECTS/wrtbad.o \
		OBJECTS/wrtdim.o OBJECTS/wvbnrg.o OBJECTS/xintrp.o \
		OBJECTS/xmlatm.o OBJECTS/xprofl.o OBJECTS/xsectb.o \
		OBJECTS/xterp.o  OBJECTS/zeroal.o OBJECTS/zeroit.o
		${FC} -O -o mod3p5.exe \
		OBJECTS/Modtrn.o OBJECTS/ab.o     OBJECTS/abcd.o   \
		OBJECTS/abcdta.o OBJECTS/aerext.o OBJECTS/aernsm.o \
		OBJECTS/aerprf.o OBJECTS/aertmp.o OBJECTS/aitk.o   \
		OBJECTS/albtrn.o OBJECTS/altrin.o OBJECTS/asymtx.o \
		OBJECTS/atmcon.o OBJECTS/bbfn.o   OBJECTS/betabs.o \
		OBJECTS/bfh2o.o  OBJECTS/bmdata.o OBJECTS/bmerfu.o \
		OBJECTS/bmflux.o OBJECTS/bmload.o OBJECTS/bmod.o   \
		OBJECTS/bmtran.o OBJECTS/bmwids.o OBJECTS/bmxld.o  \
		OBJECTS/bo2c.o   OBJECTS/bo3hh0.o OBJECTS/bo3hh1.o \
		OBJECTS/bo3hh2.o OBJECTS/bs.o     OBJECTS/c4d.o    \
		OBJECTS/c4dta.o  OBJECTS/c6dta.o  OBJECTS/check.o  \
		OBJECTS/chekin.o OBJECTS/cirr18.o OBJECTS/cirrus.o \
		OBJECTS/ckd.o    OBJECTS/cldprf.o OBJECTS/cmpint.o \
		OBJECTS/compar.o OBJECTS/convrt.o OBJECTS/cph2o.o  \
		OBJECTS/cpo3.o   OBJECTS/cptrcg.o OBJECTS/cpumix.o \
		OBJECTS/crdriv.o OBJECTS/crmerg.o OBJECTS/crprof.o \
		OBJECTS/crspec.o OBJECTS/cruprf.o OBJECTS/cruspc.o \
		OBJECTS/cxdta.o  OBJECTS/d1mach.o OBJECTS/dasum.o  \
		OBJECTS/daxpy.o  OBJECTS/dbltx.o  OBJECTS/dchdc.o  \
		OBJECTS/ddot.o   OBJECTS/debye.o  OBJECTS/defalt.o \
		OBJECTS/del.o    OBJECTS/desatt.o OBJECTS/dgbco.o  \
		OBJECTS/dgbfa.o  OBJECTS/dgbsl.o  OBJECTS/dgeco.o  \
		OBJECTS/dgefa.o  OBJECTS/dgesl.o  OBJECTS/disbbf.o \
		OBJECTS/disort.o OBJECTS/dop.o    OBJECTS/dpandx.o \
		OBJECTS/dpexnt.o OBJECTS/dpfill.o OBJECTS/dpfish.o \
		OBJECTS/dpfnmn.o OBJECTS/dplayr.o OBJECTS/dprarf.o \
		OBJECTS/dprfpa.o OBJECTS/dpscht.o OBJECTS/dref.o   \
		OBJECTS/driver.o OBJECTS/dscal.o  OBJECTS/dstdta.o \
		OBJECTS/dswap.o  OBJECTS/errmsg.o \
		OBJECTS/exabin.o OBJECTS/expint.o OBJECTS/extdta.o \
		OBJECTS/fdbeta.o OBJECTS/flayz.o  OBJECTS/fluxes.o \
		OBJECTS/flxadd.o OBJECTS/fndpth.o OBJECTS/frn296.o \
		OBJECTS/frq5dt.o OBJECTS/ftrang.o \
		OBJECTS/gamfog.o OBJECTS/geo.o    \
		OBJECTS/geoinp.o OBJECTS/gtstrt.o OBJECTS/gmrain.o \
		OBJECTS/h2src.o  OBJECTS/hengns.o OBJECTS/hertda.o \
		OBJECTS/hno3.o   OBJECTS/idamax.o OBJECTS/indx.o   \
		OBJECTS/interp.o OBJECTS/irfxn.o  OBJECTS/jou.o    \
		OBJECTS/lay5dt.o OBJECTS/laycld.o OBJECTS/layvsa.o \
		OBJECTS/lepoly.o OBJECTS/locate.o OBJECTS/loop.o   \
		OBJECTS/mapms.o  OBJECTS/mardta.o OBJECTS/marine.o \
		OBJECTS/mdta.o   OBJECTS/mlatmb.o OBJECTS/msrad.o  \
		OBJECTS/mssolr.o OBJECTS/msthml.o OBJECTS/newh2.o  \
		OBJECTS/newsrc.o OBJECTS/no2xs.o  OBJECTS/nunit.o  \
		OBJECTS/o2cont.o OBJECTS/o2int.o  OBJECTS/o3chap.o \
		OBJECTS/o3hht0.o OBJECTS/o3hht1.o OBJECTS/o3hht2.o \
		OBJECTS/o3int.o  OBJECTS/o3uv.o   OBJECTS/o3uvfd.o \
		OBJECTS/pf.o     OBJECTS/phasef.o OBJECTS/phsdta.o \
		OBJECTS/plkavg.o OBJECTS/praltr.o OBJECTS/pravin.o \
		OBJECTS/prfdta.o OBJECTS/prtinp.o OBJECTS/prtint.o \
		OBJECTS/psi.o    OBJECTS/psieca.o OBJECTS/pslct.o  \
		OBJECTS/qgausn.o OBJECTS/rain.o   OBJECTS/random.o \
		OBJECTS/ranfun.o OBJECTS/ranset.o OBJECTS/ratio.o  \
		OBJECTS/rdexa.o  OBJECTS/rdnsm.o  OBJECTS/rdsun.o  \
		OBJECTS/reduce.o OBJECTS/rhoeps.o OBJECTS/rnscat.o \
		OBJECTS/rtbis.o  OBJECTS/schrun.o OBJECTS/sctang.o \
		OBJECTS/setdis.o OBJECTS/setmtx.o OBJECTS/sf260.o  \
		OBJECTS/sf296.o  OBJECTS/shade.o  OBJECTS/shumg.o  \
		OBJECTS/sint.o   OBJECTS/slf260.o OBJECTS/slf296.o \
		OBJECTS/slftst.o OBJECTS/smgeo.o  OBJECTS/smprep.o \
		OBJECTS/so2xs.o  OBJECTS/soleig.o OBJECTS/solve0.o \
		OBJECTS/solve1.o OBJECTS/solzen.o OBJECTS/source.o \
		OBJECTS/spaltr.o OBJECTS/ssgeo.o  OBJECTS/ssrad.o  \
		OBJECTS/stdmdl.o OBJECTS/subsol.o OBJECTS/sunbd.o  \
		OBJECTS/surfac.o OBJECTS/svsola.o OBJECTS/tab.o    \
		OBJECTS/tanht.o  OBJECTS/terpev.o OBJECTS/terpso.o \
		OBJECTS/title.o  OBJECTS/tnrain.o OBJECTS/trans.o  \
		OBJECTS/trlay.o  OBJECTS/tstbad.o OBJECTS/upbeam.o \
		OBJECTS/upisot.o OBJECTS/usrint.o OBJECTS/vsa.o    \
		OBJECTS/vsansm.o OBJECTS/watvap.o OBJECTS/wrtbad.o \
		OBJECTS/wrtdim.o OBJECTS/wvbnrg.o OBJECTS/xintrp.o \
		OBJECTS/xmlatm.o OBJECTS/xprofl.o OBJECTS/xsectb.o \
		OBJECTS/xterp.o  OBJECTS/zeroal.o OBJECTS/zeroit.o

OBJECTS/Modtrn.o:	SOURCE/Modtrn.f
			${FC} -O -c       SOURCE/Modtrn.f -o OBJECTS/Modtrn.o

OBJECTS/ab.o:		SOURCE/ab.f
			${FC} -O -c       SOURCE/ab.f     -o OBJECTS/ab.o

OBJECTS/abcd.o:		SOURCE/abcd.f
			${FC} -O -c       SOURCE/abcd.f   -o OBJECTS/abcd.o

OBJECTS/abcdta.o:	SOURCE/abcdta.f
			${FC} -O -c       SOURCE/abcdta.f -o OBJECTS/abcdta.o

OBJECTS/aerext.o:	SOURCE/aerext.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/aerext.f -o OBJECTS/aerext.o

OBJECTS/aernsm.o:	SOURCE/aernsm.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/aernsm.f -o OBJECTS/aernsm.o

OBJECTS/aerprf.o:	SOURCE/aerprf.f
			${FC} -O -c       SOURCE/aerprf.f -o OBJECTS/aerprf.o

OBJECTS/aertmp.o:	SOURCE/aertmp.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/aertmp.f -o OBJECTS/aertmp.o

OBJECTS/aitk.o:		SOURCE/aitk.f
			${FC} -O -c       SOURCE/aitk.f   -o OBJECTS/aitk.o

OBJECTS/albtrn.o:	SOURCE/albtrn.f
			${FC} -O -c       SOURCE/albtrn.f -o OBJECTS/albtrn.o

OBJECTS/altrin.o:	SOURCE/altrin.f
			${FC} -O -c       SOURCE/altrin.f -o OBJECTS/altrin.o

OBJECTS/asymtx.o:	SOURCE/asymtx.f
			${FC} -O -c       SOURCE/asymtx.f -o OBJECTS/asymtx.o

OBJECTS/atmcon.o:	SOURCE/atmcon.f
			${FC} -O -c       SOURCE/atmcon.f -o OBJECTS/atmcon.o

OBJECTS/bbfn.o:		SOURCE/bbfn.f
			${FC} -O -c       SOURCE/bbfn.f   -o OBJECTS/bbfn.o

OBJECTS/betabs.o:	SOURCE/betabs.f
			${FC} -O -c       SOURCE/betabs.f -o OBJECTS/betabs.o

OBJECTS/bfh2o.o:	SOURCE/bfh2o.f
			${FC} -O -c       SOURCE/bfh2o.f  -o OBJECTS/bfh2o.o

OBJECTS/bmdata.o:	SOURCE/bmdata.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/bmdata.f -o OBJECTS/bmdata.o

OBJECTS/bmerfu.o:	SOURCE/bmerfu.f
			${FC} -O -c       SOURCE/bmerfu.f -o OBJECTS/bmerfu.o

OBJECTS/bmflux.o:	SOURCE/bmflux.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/bmflux.f -o OBJECTS/bmflux.o

OBJECTS/bmxld.o:	SOURCE/bmxld.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/bmxld.f  -o OBJECTS/bmxld.o

OBJECTS/bmload.o:	SOURCE/bmload.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/bmload.f -o OBJECTS/bmload.o

OBJECTS/bmod.o:		SOURCE/bmod.f   SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/bmod.f   -o OBJECTS/bmod.o

OBJECTS/bmtran.o:	SOURCE/bmtran.f
			${FC} -O -c       SOURCE/bmtran.f -o OBJECTS/bmtran.o

OBJECTS/bmwids.o:	SOURCE/bmwids.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/bmwids.f -o OBJECTS/bmwids.o

OBJECTS/bo2c.o:		SOURCE/bo2c.f
			${FC} -O -c       SOURCE/bo2c.f   -o OBJECTS/bo2c.o

OBJECTS/bo3hh0.o:	SOURCE/bo3hh0.f
			${FC} -O -c       SOURCE/bo3hh0.f -o OBJECTS/bo3hh0.o

OBJECTS/bo3hh1.o:	SOURCE/bo3hh1.f
			${FC} -O -c       SOURCE/bo3hh1.f -o OBJECTS/bo3hh1.o

OBJECTS/bo3hh2.o:	SOURCE/bo3hh2.f
			${FC} -O -c       SOURCE/bo3hh2.f -o OBJECTS/bo3hh2.o

OBJECTS/bs.o:		SOURCE/bs.f
			${FC} -O -c       SOURCE/bs.f     -o OBJECTS/bs.o

OBJECTS/c4d.o:		SOURCE/c4d.f
			${FC} -O -c       SOURCE/c4d.f    -o OBJECTS/c4d.o

OBJECTS/c4dta.o:	SOURCE/c4dta.f
			${FC} -O -c       SOURCE/c4dta.f  -o OBJECTS/c4dta.o

OBJECTS/c6dta.o:	SOURCE/c6dta.f
			${FC} -O -c       SOURCE/c6dta.f  -o OBJECTS/c6dta.o

OBJECTS/check.o:	SOURCE/check.f
			${FC} -O -c       SOURCE/check.f  -o OBJECTS/check.o

OBJECTS/chekin.o:	SOURCE/chekin.f
			${FC} -O -c       SOURCE/chekin.f -o OBJECTS/chekin.o

OBJECTS/cirr18.o:	SOURCE/cirr18.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/cirr18.f -o OBJECTS/cirr18.o

OBJECTS/cirrus.o:	SOURCE/cirrus.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/cirrus.f -o OBJECTS/cirrus.o

OBJECTS/ckd.o:		SOURCE/ckd.f
			${FC} -O -c       SOURCE/ckd.f    -o OBJECTS/ckd.o

OBJECTS/cldprf.o:	SOURCE/cldprf.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/cldprf.f -o OBJECTS/cldprf.o

OBJECTS/cmpint.o:	SOURCE/cmpint.f
			${FC} -O -c       SOURCE/cmpint.f -o OBJECTS/cmpint.o

OBJECTS/compar.o:	SOURCE/compar.f
			${FC} -O -c       SOURCE/compar.f -o OBJECTS/compar.o

OBJECTS/convrt.o:	SOURCE/convrt.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/convrt.f -o OBJECTS/convrt.o

OBJECTS/cph2o.o:	SOURCE/cph2o.f
			${FC} -O -c       SOURCE/cph2o.f  -o OBJECTS/cph2o.o

OBJECTS/cpo3.o:		SOURCE/cpo3.f
			${FC} -O -c       SOURCE/cpo3.f   -o OBJECTS/cpo3.o

OBJECTS/cptrcg.o:	SOURCE/cptrcg.f
			${FC} -O -c       SOURCE/cptrcg.f -o OBJECTS/cptrcg.o

OBJECTS/cpumix.o:	SOURCE/cpumix.f
			${FC} -O -c       SOURCE/cpumix.f -o OBJECTS/cpumix.o

OBJECTS/crdriv.o:	SOURCE/crdriv.f
			${FC} -O -c       SOURCE/crdriv.f -o OBJECTS/crdriv.o

OBJECTS/crmerg.o:	SOURCE/crmerg.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/crmerg.f -o OBJECTS/crmerg.o

OBJECTS/crprof.o:	SOURCE/crprof.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/crprof.f -o OBJECTS/crprof.o

OBJECTS/crspec.o:	SOURCE/crspec.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/crspec.f -o OBJECTS/crspec.o

OBJECTS/cruprf.o:	SOURCE/cruprf.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/cruprf.f -o OBJECTS/cruprf.o

OBJECTS/cruspc.o:	SOURCE/cruspc.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/cruspc.f -o OBJECTS/cruspc.o

OBJECTS/cxdta.o:	SOURCE/cxdta.f
			${FC} -O -c       SOURCE/cxdta.f  -o OBJECTS/cxdta.o

OBJECTS/d1mach.o:	SOURCE/d1mach.f
			${FC} -O -c       SOURCE/d1mach.f -o OBJECTS/d1mach.o

OBJECTS/dasum.o:	SOURCE/dasum.f
			${FC} -O -c       SOURCE/dasum.f  -o OBJECTS/dasum.o

OBJECTS/daxpy.o:	SOURCE/daxpy.f
			${FC} -O -c       SOURCE/daxpy.f  -o OBJECTS/daxpy.o

OBJECTS/dbltx.o:	SOURCE/dbltx.f
			${FC} -O -c       SOURCE/dbltx.f  -o OBJECTS/dbltx.o

OBJECTS/dchdc.o:	SOURCE/dchdc.f
			${FC} -O -c       SOURCE/dchdc.f  -o OBJECTS/dchdc.o

OBJECTS/ddot.o:		SOURCE/ddot.f
			${FC} -O -c       SOURCE/ddot.f   -o OBJECTS/ddot.o

OBJECTS/debye.o:	SOURCE/debye.f
			${FC} -O -c       SOURCE/debye.f  -o OBJECTS/debye.o

OBJECTS/defalt.o:	SOURCE/defalt.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/defalt.f -o OBJECTS/defalt.o

OBJECTS/del.o:		SOURCE/del.f
			${FC} -O -c       SOURCE/del.f    -o OBJECTS/del.o

OBJECTS/desatt.o:	SOURCE/desatt.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/desatt.f -o OBJECTS/desatt.o

OBJECTS/dgbco.o:	SOURCE/dgbco.f
			${FC} -O -c       SOURCE/dgbco.f  -o OBJECTS/dgbco.o

OBJECTS/dgbfa.o:	SOURCE/dgbfa.f
			${FC} -O -c       SOURCE/dgbfa.f  -o OBJECTS/dgbfa.o

OBJECTS/dgbsl.o:	SOURCE/dgbsl.f
			${FC} -O -c       SOURCE/dgbsl.f  -o OBJECTS/dgbsl.o

OBJECTS/dgeco.o:	SOURCE/dgeco.f
			${FC} -O -c       SOURCE/dgeco.f  -o OBJECTS/dgeco.o

OBJECTS/dgefa.o:	SOURCE/dgefa.f
			${FC} -O -c       SOURCE/dgefa.f  -o OBJECTS/dgefa.o

OBJECTS/dgesl.o:	SOURCE/dgesl.f
			${FC} -O -c       SOURCE/dgesl.f  -o OBJECTS/dgesl.o

OBJECTS/disbbf.o:	SOURCE/disbbf.f
			${FC} -O -c       SOURCE/disbbf.f -o OBJECTS/disbbf.o

OBJECTS/disort.o:	SOURCE/disort.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/disort.f -o OBJECTS/disort.o

OBJECTS/dop.o:		SOURCE/dop.f
			${FC} -O -c       SOURCE/dop.f    -o OBJECTS/dop.o

OBJECTS/dpandx.o:	SOURCE/dpandx.f
			${FC} -O -c       SOURCE/dpandx.f -o OBJECTS/dpandx.o

OBJECTS/dpexnt.o:	SOURCE/dpexnt.f
			${FC} -O -c       SOURCE/dpexnt.f -o OBJECTS/dpexnt.o

OBJECTS/dpfill.o:	SOURCE/dpfill.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/dpfill.f -o OBJECTS/dpfill.o

OBJECTS/dpfish.o:	SOURCE/dpfish.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/dpfish.f -o OBJECTS/dpfish.o

OBJECTS/dpfnmn.o:	SOURCE/dpfnmn.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/dpfnmn.f -o OBJECTS/dpfnmn.o

OBJECTS/dplayr.o:	SOURCE/dplayr.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/dplayr.f -o OBJECTS/dplayr.o

OBJECTS/dprarf.o:	SOURCE/dprarf.f
			${FC} -O -c       SOURCE/dprarf.f -o OBJECTS/dprarf.o

OBJECTS/dprfpa.o:	SOURCE/dprfpa.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/dprfpa.f -o OBJECTS/dprfpa.o

OBJECTS/dpscht.o:	SOURCE/dpscht.f
			${FC} -O -c       SOURCE/dpscht.f -o OBJECTS/dpscht.o

OBJECTS/dref.o:		SOURCE/dref.f
			${FC} -O -c       SOURCE/dref.f   -o OBJECTS/dref.o

OBJECTS/driver.o:	SOURCE/driver.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/driver.f -o OBJECTS/driver.o

OBJECTS/dscal.o:	SOURCE/dscal.f
			${FC} -O -c       SOURCE/dscal.f  -o OBJECTS/dscal.o

OBJECTS/dstdta.o:	SOURCE/dstdta.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/dstdta.f -o OBJECTS/dstdta.o

OBJECTS/dswap.o:	SOURCE/dswap.f
			${FC} -O -c       SOURCE/dswap.f  -o OBJECTS/dswap.o

OBJECTS/errmsg.o:	SOURCE/errmsg.f
			${FC} -O -c       SOURCE/errmsg.f -o OBJECTS/errmsg.o

OBJECTS/exabin.o:	SOURCE/exabin.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/exabin.f -o OBJECTS/exabin.o

OBJECTS/expint.o:	SOURCE/expint.f
			${FC} -O -c       SOURCE/expint.f -o OBJECTS/expint.o

OBJECTS/extdta.o:	SOURCE/extdta.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/extdta.f -o OBJECTS/extdta.o

OBJECTS/fdbeta.o:	SOURCE/fdbeta.f
			${FC} -O -c       SOURCE/fdbeta.f -o OBJECTS/fdbeta.o

OBJECTS/flayz.o:	SOURCE/flayz.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/flayz.f  -o OBJECTS/flayz.o

OBJECTS/fluxes.o:	SOURCE/fluxes.f
			${FC} -O -c       SOURCE/fluxes.f -o OBJECTS/fluxes.o

OBJECTS/flxadd.o:	SOURCE/flxadd.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/flxadd.f -o OBJECTS/flxadd.o

OBJECTS/fndpth.o:	SOURCE/fndpth.f
			${FC} -O -c       SOURCE/fndpth.f -o OBJECTS/fndpth.o

OBJECTS/frn296.o:	SOURCE/frn296.f
			${FC} -O -c       SOURCE/frn296.f -o OBJECTS/frn296.o

OBJECTS/frq5dt.o:	SOURCE/frq5dt.f
			${FC} -O -c       SOURCE/frq5dt.f -o OBJECTS/frq5dt.o

OBJECTS/ftrang.o:	SOURCE/ftrang.f
			${FC} -O -c       SOURCE/ftrang.f -o OBJECTS/ftrang.o

OBJECTS/gamfog.o:	SOURCE/gamfog.f
			${FC} -O -c       SOURCE/gamfog.f -o OBJECTS/gamfog.o

OBJECTS/geo.o:		SOURCE/geo.f    SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/geo.f    -o OBJECTS/geo.o

OBJECTS/geoinp.o:	SOURCE/geoinp.f
			${FC} -O -c       SOURCE/geoinp.f -o OBJECTS/geoinp.o

OBJECTS/gtstrt.o:	SOURCE/gtstrt.f
			${FC} -O -c       SOURCE/gtstrt.f -o OBJECTS/gtstrt.o

OBJECTS/gmrain.o:	SOURCE/gmrain.f
			${FC} -O -c       SOURCE/gmrain.f -o OBJECTS/gmrain.o

OBJECTS/h2src.o:	SOURCE/h2src.f
			${FC} -O -c       SOURCE/h2src.f  -o OBJECTS/h2src.o

OBJECTS/hengns.o:	SOURCE/hengns.f
			${FC} -O -c       SOURCE/hengns.f -o OBJECTS/hengns.o

OBJECTS/hertda.o:	SOURCE/hertda.f
			${FC} -O -c       SOURCE/hertda.f -o OBJECTS/hertda.o

OBJECTS/hno3.o:		SOURCE/hno3.f
			${FC} -O -c       SOURCE/hno3.f   -o OBJECTS/hno3.o

OBJECTS/idamax.o:	SOURCE/idamax.f
			${FC} -O -c       SOURCE/idamax.f -o OBJECTS/idamax.o

OBJECTS/indx.o:		SOURCE/indx.f
			${FC} -O -c       SOURCE/indx.f   -o OBJECTS/indx.o

OBJECTS/interp.o:	SOURCE/interp.f
			${FC} -O -c       SOURCE/interp.f -o OBJECTS/interp.o

OBJECTS/irfxn.o:	SOURCE/irfxn.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/irfxn.f  -o OBJECTS/irfxn.o

OBJECTS/jou.o:		SOURCE/jou.f
			${FC} -O -c       SOURCE/jou.f    -o OBJECTS/jou.o

OBJECTS/lay5dt.o:	SOURCE/lay5dt.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/lay5dt.f -o OBJECTS/lay5dt.o

OBJECTS/laycld.o:	SOURCE/laycld.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/laycld.f -o OBJECTS/laycld.o

OBJECTS/layvsa.o:	SOURCE/layvsa.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/layvsa.f -o OBJECTS/layvsa.o

OBJECTS/lepoly.o:	SOURCE/lepoly.f
			${FC} -O -c       SOURCE/lepoly.f -o OBJECTS/lepoly.o

OBJECTS/locate.o:	SOURCE/locate.f
			${FC} -O -c       SOURCE/locate.f -o OBJECTS/locate.o

OBJECTS/loop.o:		SOURCE/loop.f   SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/loop.f   -o OBJECTS/loop.o

OBJECTS/mapms.o:	SOURCE/mapms.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/mapms.f  -o OBJECTS/mapms.o

OBJECTS/mardta.o:	SOURCE/mardta.f
			${FC} -O -c       SOURCE/mardta.f -o OBJECTS/mardta.o

OBJECTS/marine.o:	SOURCE/marine.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/marine.f -o OBJECTS/marine.o

OBJECTS/mdta.o:		SOURCE/mdta.f   SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/mdta.f   -o OBJECTS/mdta.o

OBJECTS/mlatmb.o:	SOURCE/mlatmb.f
			${FC} -O -c       SOURCE/mlatmb.f -o OBJECTS/mlatmb.o

OBJECTS/msrad.o:	SOURCE/msrad.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/msrad.f  -o OBJECTS/msrad.o

OBJECTS/mssolr.o:	SOURCE/mssolr.f
			${FC} -O -c       SOURCE/mssolr.f -o OBJECTS/mssolr.o

OBJECTS/msthml.o:	SOURCE/msthml.f
			${FC} -O -c       SOURCE/msthml.f -o OBJECTS/msthml.o

OBJECTS/newh2.o:	SOURCE/newh2.f
			${FC} -O -c       SOURCE/newh2.f  -o OBJECTS/newh2.o

OBJECTS/newsrc.o:	SOURCE/newsrc.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/newsrc.f -o OBJECTS/newsrc.o

OBJECTS/no2xs.o:	SOURCE/no2xs.f
			${FC} -O -c       SOURCE/no2xs.f  -o OBJECTS/no2xs.o

OBJECTS/nunit.o:	SOURCE/nunit.f 
			${FC} -O -c       SOURCE/nunit.f -o OBJECTS/nunit.o

OBJECTS/o2cont.o:	SOURCE/o2cont.f
			${FC} -O -c       SOURCE/o2cont.f -o OBJECTS/o2cont.o

OBJECTS/o2int.o:	SOURCE/o2int.f
			${FC} -O -c       SOURCE/o2int.f  -o OBJECTS/o2int.o

OBJECTS/o3chap.o:	SOURCE/o3chap.f
			${FC} -O -c       SOURCE/o3chap.f -o OBJECTS/o3chap.o

OBJECTS/o3hht0.o:	SOURCE/o3hht0.f
			${FC} -O -c       SOURCE/o3hht0.f -o OBJECTS/o3hht0.o

OBJECTS/o3hht1.o:	SOURCE/o3hht1.f
			${FC} -O -c       SOURCE/o3hht1.f -o OBJECTS/o3hht1.o

OBJECTS/o3hht2.o:	SOURCE/o3hht2.f
			${FC} -O -c       SOURCE/o3hht2.f -o OBJECTS/o3hht2.o

OBJECTS/o3int.o:	SOURCE/o3int.f
			${FC} -O -c       SOURCE/o3int.f  -o OBJECTS/o3int.o

OBJECTS/o3uv.o:		SOURCE/o3uv.f
			${FC} -O -c       SOURCE/o3uv.f   -o OBJECTS/o3uv.o

OBJECTS/o3uvfd.o:	SOURCE/o3uvfd.f
			${FC} -O -c       SOURCE/o3uvfd.f -o OBJECTS/o3uvfd.o

OBJECTS/pf.o:		SOURCE/pf.f
			${FC} -O -c       SOURCE/pf.f     -o OBJECTS/pf.o

OBJECTS/phasef.o:	SOURCE/phasef.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/phasef.f -o OBJECTS/phasef.o

OBJECTS/phsdta.o:	SOURCE/phsdta.f
			${FC} -O -c       SOURCE/phsdta.f -o OBJECTS/phsdta.o

OBJECTS/plkavg.o:	SOURCE/plkavg.f
			${FC} -O -c       SOURCE/plkavg.f -o OBJECTS/plkavg.o

OBJECTS/praltr.o:	SOURCE/praltr.f
			${FC} -O -c       SOURCE/praltr.f -o OBJECTS/praltr.o

OBJECTS/pravin.o:	SOURCE/pravin.f
			${FC} -O -c       SOURCE/pravin.f -o OBJECTS/pravin.o

OBJECTS/prfdta.o:	SOURCE/prfdta.f
			${FC} -O -c       SOURCE/prfdta.f -o OBJECTS/prfdta.o

OBJECTS/prtinp.o:	SOURCE/prtinp.f
			${FC} -O -c       SOURCE/prtinp.f -o OBJECTS/prtinp.o

OBJECTS/prtint.o:	SOURCE/prtint.f
			${FC} -O -c       SOURCE/prtint.f -o OBJECTS/prtint.o

OBJECTS/psi.o:		SOURCE/psi.f
			${FC} -O -c       SOURCE/psi.f    -o OBJECTS/psi.o

OBJECTS/psieca.o:	SOURCE/psieca.f
			${FC} -O -c       SOURCE/psieca.f -o OBJECTS/psieca.o

OBJECTS/pslct.o:	SOURCE/pslct.f
			${FC} -O -c       SOURCE/pslct.f  -o OBJECTS/pslct.o

OBJECTS/qgausn.o:	SOURCE/qgausn.f
			${FC} -O -c       SOURCE/qgausn.f -o OBJECTS/qgausn.o

OBJECTS/rain.o:		SOURCE/rain.f   SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/rain.f   -o OBJECTS/rain.o

OBJECTS/random.o:	SOURCE/random.f
			${FC} -O -c       SOURCE/random.f -o OBJECTS/random.o

OBJECTS/ranfun.o:	SOURCE/ranfun.f
			${FC} -O -c       SOURCE/ranfun.f -o OBJECTS/ranfun.o

OBJECTS/ranset.o:	SOURCE/ranset.f
			${FC} -O -c       SOURCE/ranset.f -o OBJECTS/ranset.o

OBJECTS/ratio.o:	SOURCE/ratio.f
			${FC} -O -c       SOURCE/ratio.f  -o OBJECTS/ratio.o

OBJECTS/rdexa.o:	SOURCE/rdexa.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/rdexa.f  -o OBJECTS/rdexa.o

OBJECTS/rdnsm.o:	SOURCE/rdnsm.f
			${FC} -O -c       SOURCE/rdnsm.f  -o OBJECTS/rdnsm.o

OBJECTS/rdsun.o:	SOURCE/rdsun.f 
			${FC} -O -c       SOURCE/rdsun.f -o OBJECTS/rdsun.o

OBJECTS/reduce.o:	SOURCE/reduce.f
			${FC} -O -c       SOURCE/reduce.f -o OBJECTS/reduce.o

OBJECTS/rhoeps.o:	SOURCE/rhoeps.f
			${FC} -O -c       SOURCE/rhoeps.f -o OBJECTS/rhoeps.o

OBJECTS/rnscat.o:	SOURCE/rnscat.f
			${FC} -O -c       SOURCE/rnscat.f -o OBJECTS/rnscat.o

OBJECTS/rtbis.o:	SOURCE/rtbis.f
			${FC} -O -c       SOURCE/rtbis.f  -o OBJECTS/rtbis.o

OBJECTS/schrun.o:	SOURCE/schrun.f
			${FC} -O -c       SOURCE/schrun.f -o OBJECTS/schrun.o

OBJECTS/sctang.o:	SOURCE/sctang.f
			${FC} -O -c       SOURCE/sctang.f -o OBJECTS/sctang.o

OBJECTS/setdis.o:	SOURCE/setdis.f
			${FC} -O -c       SOURCE/setdis.f -o OBJECTS/setdis.o

OBJECTS/setmtx.o:	SOURCE/setmtx.f
			${FC} -O -c       SOURCE/setmtx.f -o OBJECTS/setmtx.o

OBJECTS/sf260.o:	SOURCE/sf260.f
			${FC} -O -c       SOURCE/sf260.f  -o OBJECTS/sf260.o

OBJECTS/sf296.o:	SOURCE/sf296.f
			${FC} -O -c       SOURCE/sf296.f  -o OBJECTS/sf296.o

OBJECTS/shade.o:	SOURCE/shade.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/shade.f  -o OBJECTS/shade.o

OBJECTS/shumg.o:	SOURCE/shumg.f
			${FC} -O -c       SOURCE/shumg.f  -o OBJECTS/shumg.o

OBJECTS/sint.o:		SOURCE/sint.f
			${FC} -O -c       SOURCE/sint.f   -o OBJECTS/sint.o

OBJECTS/slf260.o:	SOURCE/slf260.f
			${FC} -O -c       SOURCE/slf260.f -o OBJECTS/slf260.o

OBJECTS/slf296.o:	SOURCE/slf296.f
			${FC} -O -c       SOURCE/slf296.f -o OBJECTS/slf296.o

OBJECTS/slftst.o:	SOURCE/slftst.f
			${FC} -O -c       SOURCE/slftst.f -o OBJECTS/slftst.o

OBJECTS/smgeo.o:	SOURCE/smgeo.f
			${FC} -O -c       SOURCE/smgeo.f  -o OBJECTS/smgeo.o

OBJECTS/smprep.o:	SOURCE/smprep.f
			${FC} -O -c       SOURCE/smprep.f -o OBJECTS/smprep.o

OBJECTS/so2xs.o:	SOURCE/so2xs.f
			${FC} -O -c       SOURCE/so2xs.f  -o OBJECTS/so2xs.o

OBJECTS/sola1.o:	SOURCE/sola1.f
			${FC} -O -c       SOURCE/sola1.f  -o OBJECTS/sola1.o

OBJECTS/soleig.o:	SOURCE/soleig.f
			${FC} -O -c       SOURCE/soleig.f -o OBJECTS/soleig.o

OBJECTS/solve0.o:	SOURCE/solve0.f
			${FC} -O -c       SOURCE/solve0.f -o OBJECTS/solve0.o

OBJECTS/solve1.o:	SOURCE/solve1.f
			${FC} -O -c       SOURCE/solve1.f -o OBJECTS/solve1.o

OBJECTS/solzen.o:	SOURCE/solzen.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/solzen.f -o OBJECTS/solzen.o

OBJECTS/source.o:	SOURCE/source.f
			${FC} -O -c       SOURCE/source.f -o OBJECTS/source.o

OBJECTS/spaltr.o:	SOURCE/spaltr.f
			${FC} -O -c       SOURCE/spaltr.f -o OBJECTS/spaltr.o

OBJECTS/ssgeo.o:	SOURCE/ssgeo.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/ssgeo.f  -o OBJECTS/ssgeo.o

OBJECTS/ssrad.o:	SOURCE/ssrad.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/ssrad.f  -o OBJECTS/ssrad.o

OBJECTS/stdmdl.o:	SOURCE/stdmdl.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/stdmdl.f -o OBJECTS/stdmdl.o

OBJECTS/subsol.o:	SOURCE/subsol.f
			${FC} -O -c       SOURCE/subsol.f -o OBJECTS/subsol.o

OBJECTS/sunbd.o:	SOURCE/sunbd.f
			${FC} -O -c       SOURCE/sunbd.f    -o OBJECTS/sunbd.o

OBJECTS/surfac.o:	SOURCE/surfac.f
			${FC} -O -c       SOURCE/surfac.f -o OBJECTS/surfac.o

OBJECTS/svsola.o:	SOURCE/svsola.f
			${FC} -O -c       SOURCE/svsola.f -o OBJECTS/svsola.o

OBJECTS/tab.o:		SOURCE/tab.f
			${FC} -O -c       SOURCE/tab.f    -o OBJECTS/tab.o

OBJECTS/tanht.o:	SOURCE/tanht.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/tanht.f  -o OBJECTS/tanht.o

OBJECTS/terpev.o:	SOURCE/terpev.f
			${FC} -O -c       SOURCE/terpev.f -o OBJECTS/terpev.o

OBJECTS/terpso.o:	SOURCE/terpso.f
			${FC} -O -c       SOURCE/terpso.f -o OBJECTS/terpso.o

OBJECTS/title.o:	SOURCE/title.f
			${FC} -O -c       SOURCE/title.f  -o OBJECTS/title.o

OBJECTS/tnrain.o:	SOURCE/tnrain.f
			${FC} -O -c       SOURCE/tnrain.f -o OBJECTS/tnrain.o

OBJECTS/trans.o:	SOURCE/trans.f  SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/trans.f  -o OBJECTS/trans.o

OBJECTS/trlay.o:	SOURCE/trlay.f
			${FC} -O -c       SOURCE/trlay.f  -o OBJECTS/trlay.o

OBJECTS/tstbad.o:	SOURCE/tstbad.f
			${FC} -O -c       SOURCE/tstbad.f -o OBJECTS/tstbad.o

OBJECTS/upbeam.o:	SOURCE/upbeam.f
			${FC} -O -c       SOURCE/upbeam.f -o OBJECTS/upbeam.o

OBJECTS/upisot.o:	SOURCE/upisot.f
			${FC} -O -c       SOURCE/upisot.f -o OBJECTS/upisot.o

OBJECTS/usrint.o:	SOURCE/usrint.f
			${FC} -O -c       SOURCE/usrint.f -o OBJECTS/usrint.o

OBJECTS/vsa.o:		SOURCE/vsa.f
			${FC} -O -c       SOURCE/vsa.f    -o OBJECTS/vsa.o

OBJECTS/vsansm.o:	SOURCE/vsansm.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/vsansm.f -o OBJECTS/vsansm.o

OBJECTS/watvap.o:	SOURCE/watvap.f
			${FC} -O -c       SOURCE/watvap.f -o OBJECTS/watvap.o

OBJECTS/wrtbad.o:	SOURCE/wrtbad.f
			${FC} -O -c       SOURCE/wrtbad.f -o OBJECTS/wrtbad.o

OBJECTS/wrtdim.o:	SOURCE/wrtdim.f
			${FC} -O -c       SOURCE/wrtdim.f -o OBJECTS/wrtdim.o

OBJECTS/wvbnrg.o:	SOURCE/wvbnrg.f
			${FC} -O -c       SOURCE/wvbnrg.f -o OBJECTS/wvbnrg.o

OBJECTS/xintrp.o:	SOURCE/xintrp.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/xintrp.f -o OBJECTS/xintrp.o

OBJECTS/xmlatm.o:	SOURCE/xmlatm.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/xmlatm.f -o OBJECTS/xmlatm.o

OBJECTS/xprofl.o:	SOURCE/xprofl.f SOURCE/PARAM.LST
			${FC} -O -c       SOURCE/xprofl.f -o OBJECTS/xprofl.o

OBJECTS/xsectb.o:	SOURCE/xsectb.f
			${FC} -O -c       SOURCE/xsectb.f -o OBJECTS/xsectb.o

OBJECTS/xterp.o:	SOURCE/xterp.f
			${FC} -O -c       SOURCE/xterp.f  -o OBJECTS/xterp.o

OBJECTS/zeroal.o:	SOURCE/zeroal.f
			${FC} -O -c       SOURCE/zeroal.f -o OBJECTS/zeroal.o

OBJECTS/zeroit.o:	SOURCE/zeroit.f
			${FC} -O -c       SOURCE/zeroit.f -o OBJECTS/zeroit.o

#
#
#  #############################################################
#  ##                                                         ##
#  ##  linkapbs.make  --  link the TINKER programs with APBS  ##
#  ##        (Intel Fortran Compiler for Linux Version)       ##
#  ##                                                         ##
#  #############################################################
#
#
#  The libraries in APBSLIBS are generated during the standard
#  build procedure for the APBS software package, and should be
#  moved to a directory compatible with the APBSLIB definition
#
#
setenv APBSLIBS "../apbs/linux/lib/libapbsmainroutines.a ../apbs/linux/lib/libapbs.a ../apbs/linux/lib/libmaloc.a ../apbs/linux/lib/libapbsblas.a"
ifort -fast -no-ipo -static-intel -o alchemy.x alchemy.o libtinker.a $APBSLIBS ; strip alchemy.x
ifort -fast -no-ipo -static-intel -o analyze.x analyze.o libtinker.a $APBSLIBS ; strip analyze.x
ifort -fast -no-ipo -static-intel -o anneal.x anneal.o libtinker.a $APBSLIBS ; strip anneal.x
ifort -fast -no-ipo -static-intel -o archive.x archive.o libtinker.a $APBSLIBS ; strip archive.x
ifort -fast -no-ipo -static-intel -o correlate.x correlate.o libtinker.a $APBSLIBS ; strip correlate.x
ifort -fast -no-ipo -static-intel -o crystal.x crystal.o libtinker.a $APBSLIBS ; strip crystal.x
ifort -fast -no-ipo -static-intel -o diffuse.x diffuse.o libtinker.a $APBSLIBS ; strip diffuse.x
ifort -fast -no-ipo -static-intel -o distgeom.x distgeom.o libtinker.a $APBSLIBS ; strip distgeom.x
ifort -fast -no-ipo -static-intel -o document.x document.o libtinker.a $APBSLIBS ; strip document.x
ifort -fast -no-ipo -static-intel -o dynamic.x dynamic.o libtinker.a $APBSLIBS ; strip dynamic.x
ifort -fast -no-ipo -static-intel -o gda.x gda.o libtinker.a $APBSLIBS ; strip gda.x
ifort -fast -no-ipo -static-intel -o intedit.x intedit.o libtinker.a $APBSLIBS ; strip intedit.x
ifort -fast -no-ipo -static-intel -o intxyz.x intxyz.o libtinker.a $APBSLIBS ; strip intxyz.x
ifort -fast -no-ipo -static-intel -o minimize.x minimize.o libtinker.a $APBSLIBS ; strip minimize.x
ifort -fast -no-ipo -static-intel -o minirot.x minirot.o libtinker.a $APBSLIBS ; strip minirot.x
ifort -fast -no-ipo -static-intel -o minrigid.x minrigid.o libtinker.a $APBSLIBS ; strip minrigid.x
ifort -fast -no-ipo -static-intel -o monte.x monte.o libtinker.a $APBSLIBS ; strip monte.x
ifort -fast -no-ipo -static-intel -o newton.x newton.o libtinker.a $APBSLIBS ; strip newton.x
ifort -fast -no-ipo -static-intel -o newtrot.x newtrot.o libtinker.a $APBSLIBS ; strip newtrot.x
ifort -fast -no-ipo -static-intel -o nucleic.x nucleic.o libtinker.a $APBSLIBS ; strip nucleic.x
ifort -fast -no-ipo -static-intel -o optimize.x optimize.o libtinker.a $APBSLIBS ; strip optimize.x
ifort -fast -no-ipo -static-intel -o optirot.x optirot.o libtinker.a $APBSLIBS ; strip optirot.x
ifort -fast -no-ipo -static-intel -o optrigid.x optrigid.o libtinker.a $APBSLIBS ; strip optrigid.x
ifort -fast -no-ipo -static-intel -o path.x path.o libtinker.a $APBSLIBS ; strip path.x
ifort -fast -no-ipo -static-intel -o pdbxyz.x pdbxyz.o libtinker.a $APBSLIBS ; strip pdbxyz.x
ifort -fast -no-ipo -static-intel -o polarize.x polarize.o libtinker.a $APBSLIBS ; strip polarize.x
ifort -fast -no-ipo -static-intel -o poledit.x poledit.o libtinker.a $APBSLIBS ; strip poledit.x
ifort -fast -no-ipo -static-intel -o potential.x potential.o libtinker.a $APBSLIBS ; strip potential.x
ifort -fast -no-ipo -static-intel -o prmedit.x prmedit.o libtinker.a $APBSLIBS ; strip prmedit.x
ifort -fast -no-ipo -static-intel -o protein.x protein.o libtinker.a $APBSLIBS ; strip protein.x
ifort -fast -no-ipo -static-intel -o pss.x pss.o libtinker.a $APBSLIBS ; strip pss.x
ifort -fast -no-ipo -static-intel -o pssrigid.x pssrigid.o libtinker.a $APBSLIBS ; strip pssrigid.x
ifort -fast -no-ipo -static-intel -o pssrot.x pssrot.o libtinker.a $APBSLIBS ; strip pssrot.x
ifort -fast -no-ipo -static-intel -o radial.x radial.o libtinker.a $APBSLIBS ; strip radial.x
ifort -fast -no-ipo -static-intel -o saddle.x saddle.o libtinker.a $APBSLIBS ; strip saddle.x
ifort -fast -no-ipo -static-intel -o scan.x scan.o libtinker.a $APBSLIBS ; strip scan.x
ifort -fast -no-ipo -static-intel -o sniffer.x sniffer.o libtinker.a $APBSLIBS ; strip sniffer.x
ifort -fast -no-ipo -static-intel -o spacefill.x spacefill.o libtinker.a $APBSLIBS ; strip spacefill.x
ifort -fast -no-ipo -static-intel -o spectrum.x spectrum.o libtinker.a $APBSLIBS ; strip spectrum.x
ifort -fast -no-ipo -static-intel -o superpose.x superpose.o libtinker.a $APBSLIBS ; strip superpose.x
ifort -fast -no-ipo -static-intel -o sybylxyz.x sybylxyz.o libtinker.a $APBSLIBS ; strip sybylxyz.x
ifort -fast -no-ipo -static-intel -o testgrad.x testgrad.o libtinker.a $APBSLIBS ; strip testgrad.x
ifort -fast -no-ipo -static-intel -o testhess.x testhess.o libtinker.a $APBSLIBS ; strip testhess.x
ifort -fast -no-ipo -static-intel -o testpair.x testpair.o libtinker.a $APBSLIBS ; strip testpair.x
ifort -fast -no-ipo -static-intel -o testrot.x testrot.o libtinker.a $APBSLIBS ; strip testrot.x
ifort -fast -no-ipo -static-intel -o timer.x timer.o libtinker.a $APBSLIBS ; strip timer.x
ifort -fast -no-ipo -static-intel -o timerot.x timerot.o libtinker.a $APBSLIBS ; strip timerot.x
ifort -fast -no-ipo -static-intel -o torsfit.x torsfit.o libtinker.a $APBSLIBS ; strip torsfit.x
ifort -fast -no-ipo -static-intel -o valence.x valence.o libtinker.a $APBSLIBS ; strip valence.x
ifort -fast -no-ipo -static-intel -o vibbig.x vibbig.o libtinker.a $APBSLIBS ; strip vibbig.x
ifort -fast -no-ipo -static-intel -o vibrate.x vibrate.o libtinker.a $APBSLIBS ; strip vibrate.x
ifort -fast -no-ipo -static-intel -o vibrot.x vibrot.o libtinker.a $APBSLIBS ; strip vibrot.x
ifort -fast -no-ipo -static-intel -o xtalfit.x xtalfit.o libtinker.a $APBSLIBS ; strip xtalfit.x
ifort -fast -no-ipo -static-intel -o xtalmin.x xtalmin.o libtinker.a $APBSLIBS ; strip xtalmin.x
ifort -fast -no-ipo -static-intel -o xyzedit.x xyzedit.o libtinker.a $APBSLIBS ; strip xyzedit.x
ifort -fast -no-ipo -static-intel -o xyzint.x xyzint.o libtinker.a $APBSLIBS ; strip xyzint.x
ifort -fast -no-ipo -static-intel -o xyzpdb.x xyzpdb.o libtinker.a $APBSLIBS ; strip xyzpdb.x
ifort -fast -no-ipo -static-intel -o xyzsybyl.x xyzsybyl.o libtinker.a $APBSLIBS ; strip xyzsybyl.x

#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##         (Intel Fortran Compiler for Linux Version)        ##
#  ##                                                           ##
#  ###############################################################
#
#
ifort -fast -no-ipo -static-intel -o alchemy.x alchemy.o libtinker.a ; strip alchemy.x
ifort -fast -no-ipo -static-intel -o analyze.x analyze.o libtinker.a ; strip analyze.x
ifort -fast -no-ipo -static-intel -o anneal.x anneal.o libtinker.a ; strip anneal.x
ifort -fast -no-ipo -static-intel -o archive.x archive.o libtinker.a ; strip archive.x
ifort -fast -no-ipo -static-intel -o correlate.x correlate.o libtinker.a ; strip correlate.x
ifort -fast -no-ipo -static-intel -o crystal.x crystal.o libtinker.a ; strip crystal.x
ifort -fast -no-ipo -static-intel -o diffuse.x diffuse.o libtinker.a ; strip diffuse.x
ifort -fast -no-ipo -static-intel -o distgeom.x distgeom.o libtinker.a ; strip distgeom.x
ifort -fast -no-ipo -static-intel -o document.x document.o libtinker.a ; strip document.x
ifort -fast -no-ipo -static-intel -o dynamic.x dynamic.o libtinker.a ; strip dynamic.x
ifort -fast -no-ipo -static-intel -o gda.x gda.o libtinker.a ; strip gda.x
ifort -fast -no-ipo -static-intel -o intedit.x intedit.o libtinker.a ; strip intedit.x
ifort -fast -no-ipo -static-intel -o intxyz.x intxyz.o libtinker.a ; strip intxyz.x
ifort -fast -no-ipo -static-intel -o minimize.x minimize.o libtinker.a ; strip minimize.x
ifort -fast -no-ipo -static-intel -o minirot.x minirot.o libtinker.a ; strip minirot.x
ifort -fast -no-ipo -static-intel -o minrigid.x minrigid.o libtinker.a ; strip minrigid.x
ifort -fast -no-ipo -static-intel -o monte.x monte.o libtinker.a ; strip monte.x
ifort -fast -no-ipo -static-intel -o newton.x newton.o libtinker.a ; strip newton.x
ifort -fast -no-ipo -static-intel -o newtrot.x newtrot.o libtinker.a ; strip newtrot.x
ifort -fast -no-ipo -static-intel -o nucleic.x nucleic.o libtinker.a ; strip nucleic.x
ifort -fast -no-ipo -static-intel -o optimize.x optimize.o libtinker.a ; strip optimize.x
ifort -fast -no-ipo -static-intel -o optirot.x optirot.o libtinker.a ; strip optirot.x
ifort -fast -no-ipo -static-intel -o optrigid.x optrigid.o libtinker.a ; strip optrigid.x
ifort -fast -no-ipo -static-intel -o path.x path.o libtinker.a ; strip path.x
ifort -fast -no-ipo -static-intel -o pdbxyz.x pdbxyz.o libtinker.a ; strip pdbxyz.x
ifort -fast -no-ipo -static-intel -o polarize.x polarize.o libtinker.a ; strip polarize.x
ifort -fast -no-ipo -static-intel -o poledit.x poledit.o libtinker.a ; strip poledit.x
ifort -fast -no-ipo -static-intel -o potential.x potential.o libtinker.a ; strip potential.x
ifort -fast -no-ipo -static-intel -o prmedit.x prmedit.o libtinker.a ; strip prmedit.x
ifort -fast -no-ipo -static-intel -o protein.x protein.o libtinker.a ; strip protein.x
ifort -fast -no-ipo -static-intel -o pss.x pss.o libtinker.a ; strip pss.x
ifort -fast -no-ipo -static-intel -o pssrigid.x pssrigid.o libtinker.a ; strip pssrigid.x
ifort -fast -no-ipo -static-intel -o pssrot.x pssrot.o libtinker.a ; strip pssrot.x
ifort -fast -no-ipo -static-intel -o radial.x radial.o libtinker.a ; strip radial.x
ifort -fast -no-ipo -static-intel -o saddle.x saddle.o libtinker.a ; strip saddle.x
ifort -fast -no-ipo -static-intel -o scan.x scan.o libtinker.a ; strip scan.x
ifort -fast -no-ipo -static-intel -o sniffer.x sniffer.o libtinker.a ; strip sniffer.x
ifort -fast -no-ipo -static-intel -o spacefill.x spacefill.o libtinker.a ; strip spacefill.x
ifort -fast -no-ipo -static-intel -o spectrum.x spectrum.o libtinker.a ; strip spectrum.x
ifort -fast -no-ipo -static-intel -o superpose.x superpose.o libtinker.a ; strip superpose.x
ifort -fast -no-ipo -static-intel -o sybylxyz.x sybylxyz.o libtinker.a ; strip sybylxyz.x
ifort -fast -no-ipo -static-intel -o testgrad.x testgrad.o libtinker.a ; strip testgrad.x
ifort -fast -no-ipo -static-intel -o testhess.x testhess.o libtinker.a ; strip testhess.x
ifort -fast -no-ipo -static-intel -o testpair.x testpair.o libtinker.a ; strip testpair.x
ifort -fast -no-ipo -static-intel -o testrot.x testrot.o libtinker.a ; strip testrot.x
ifort -fast -no-ipo -static-intel -o timer.x timer.o libtinker.a ; strip timer.x
ifort -fast -no-ipo -static-intel -o timerot.x timerot.o libtinker.a ; strip timerot.x
ifort -fast -no-ipo -static-intel -o valence.x valence.o libtinker.a ; strip valence.x
ifort -fast -no-ipo -static-intel -o vibbig.x vibbig.o libtinker.a ; strip vibbig.x
ifort -fast -no-ipo -static-intel -o vibrate.x vibrate.o libtinker.a ; strip vibrate.x
ifort -fast -no-ipo -static-intel -o vibrot.x vibrot.o libtinker.a ; strip vibrot.x
ifort -fast -no-ipo -static-intel -o xtalfit.x xtalfit.o libtinker.a ; strip xtalfit.x
ifort -fast -no-ipo -static-intel -o xtalmin.x xtalmin.o libtinker.a ; strip xtalmin.x
ifort -fast -no-ipo -static-intel -o xyzedit.x xyzedit.o libtinker.a ; strip xyzedit.x
ifort -fast -no-ipo -static-intel -o xyzint.x xyzint.o libtinker.a ; strip xyzint.x
ifort -fast -no-ipo -static-intel -o xyzpdb.x xyzpdb.o libtinker.a ; strip xyzpdb.x
ifort -fast -no-ipo -static-intel -o xyzsybyl.x xyzsybyl.o libtinker.a ; strip xyzsybyl.x

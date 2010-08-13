#
#
#  ####################################################################
#  ##                                                                ##
#  ##  linkserial.make  --  link all the TINKER programs for serial  ##
#  ##          (Intel Fortran Compiler for Mac OSX Version)          ##
#  ##                                                                ##
#  ####################################################################
#
#
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o alchemy.x alchemy.o libtinker.a ; strip alchemy.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o analyze.x analyze.o libtinker.a ; strip analyze.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o anneal.x anneal.o libtinker.a ; strip anneal.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o archive.x archive.o libtinker.a ; strip archive.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o correlate.x correlate.o libtinker.a ; strip correlate.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o crystal.x crystal.o libtinker.a ; strip crystal.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o diffuse.x diffuse.o libtinker.a ; strip diffuse.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o distgeom.x distgeom.o libtinker.a ; strip distgeom.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o document.x document.o libtinker.a ; strip document.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o dynamic.x dynamic.o libtinker.a ; strip dynamic.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o gda.x gda.o libtinker.a ; strip gda.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o intedit.x intedit.o libtinker.a ; strip intedit.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o intxyz.x intxyz.o libtinker.a ; strip intxyz.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o minimize.x minimize.o libtinker.a ; strip minimize.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o minirot.x minirot.o libtinker.a ; strip minirot.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o minrigid.x minrigid.o libtinker.a ; strip minrigid.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o monte.x monte.o libtinker.a ; strip monte.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o newton.x newton.o libtinker.a ; strip newton.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o newtrot.x newtrot.o libtinker.a ; strip newtrot.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o nucleic.x nucleic.o libtinker.a ; strip nucleic.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o optimize.x optimize.o libtinker.a ; strip optimize.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o optirot.x optirot.o libtinker.a ; strip optirot.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o optrigid.x optrigid.o libtinker.a ; strip optrigid.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o path.x path.o libtinker.a ; strip path.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o pdbxyz.x pdbxyz.o libtinker.a ; strip pdbxyz.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o polarize.x polarize.o libtinker.a ; strip polarize.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o poledit.x poledit.o libtinker.a ; strip poledit.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o potential.x potential.o libtinker.a ; strip potential.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o prmedit.x prmedit.o libtinker.a ; strip prmedit.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o protein.x protein.o libtinker.a ; strip protein.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o pss.x pss.o libtinker.a ; strip pss.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o pssrigid.x pssrigid.o libtinker.a ; strip pssrigid.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o pssrot.x pssrot.o libtinker.a ; strip pssrot.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o radial.x radial.o libtinker.a ; strip radial.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o saddle.x saddle.o libtinker.a ; strip saddle.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o scan.x scan.o libtinker.a ; strip scan.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o sniffer.x sniffer.o libtinker.a ; strip sniffer.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o spacefill.x spacefill.o libtinker.a ; strip spacefill.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o spectrum.x spectrum.o libtinker.a ; strip spectrum.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o superpose.x superpose.o libtinker.a ; strip superpose.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o sybylxyz.x sybylxyz.o libtinker.a ; strip sybylxyz.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o testgrad.x testgrad.o libtinker.a ; strip testgrad.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o testhess.x testhess.o libtinker.a ; strip testhess.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o testpair.x testpair.o libtinker.a ; strip testpair.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o testrot.x testrot.o libtinker.a ; strip testrot.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o timer.x timer.o libtinker.a ; strip timer.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o timerot.x timerot.o libtinker.a ; strip timerot.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o torsfit.x torsfit.o libtinker.a ; strip torsfit.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o valence.x valence.o libtinker.a ; strip valence.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o vibbig.x vibbig.o libtinker.a ; strip vibbig.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o vibrate.x vibrate.o libtinker.a ; strip vibrate.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o vibrot.x vibrot.o libtinker.a ; strip vibrot.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o xtalfit.x xtalfit.o libtinker.a ; strip xtalfit.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o xtalmin.x xtalmin.o libtinker.a ; strip xtalmin.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o xyzedit.x xyzedit.o libtinker.a ; strip xyzedit.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o xyzint.x xyzint.o libtinker.a ; strip xyzint.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o xyzpdb.x xyzpdb.o libtinker.a ; strip xyzpdb.x
ifort -O3 -no-ipo -no-prec-div -static-intel -mmacosx-version-min=10.4 -o xyzsybyl.x xyzsybyl.o libtinker.a ; strip xyzsybyl.x

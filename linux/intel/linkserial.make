#
#
#  ################################################################
#  ##                                                            ##
#  ##  linkserial.make  --  link the Tinker programs for serial  ##
#  ##              (Intel Fortran for Linux Version)             ##
#  ##                                                            ##
#  ################################################################
#
#
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o alchemy.x alchemy.o libtinker.a ; strip alchemy.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o analyze.x analyze.o libtinker.a ; strip analyze.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o anneal.x anneal.o libtinker.a ; strip anneal.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o arcedit.x arcedit.o libtinker.a ; strip arcedit.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o bar.x bar.o libtinker.a ; strip bar.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o correlate.x correlate.o libtinker.a ; strip correlate.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o critical.x critical.o libtinker.a ; strip critical.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o crystal.x crystal.o libtinker.a ; strip crystal.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o diffuse.x diffuse.o libtinker.a ; strip diffuse.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o distgeom.x distgeom.o libtinker.a ; strip distgeom.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o document.x document.o libtinker.a ; strip document.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o dynamic.x dynamic.o libtinker.a ; strip dynamic.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o freefix.x freefix.o libtinker.a ; strip freefix.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o gda.x gda.o libtinker.a ; strip gda.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o intedit.x intedit.o libtinker.a ; strip intedit.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o intxyz.x intxyz.o libtinker.a ; strip intxyz.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o minimize.x minimize.o libtinker.a ; strip minimize.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o minirot.x minirot.o libtinker.a ; strip minirot.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o minrigid.x minrigid.o libtinker.a ; strip minrigid.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o mol2xyz.x mol2xyz.o libtinker.a ; strip mol2xyz.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o molxyz.x molxyz.o libtinker.a ; strip molxyz.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o monte.x monte.o libtinker.a ; strip monte.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o newton.x newton.o libtinker.a ; strip newton.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o newtrot.x newtrot.o libtinker.a ; strip newtrot.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o nucleic.x nucleic.o libtinker.a ; strip nucleic.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o optimize.x optimize.o libtinker.a ; strip optimize.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o optirot.x optirot.o libtinker.a ; strip optirot.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o optrigid.x optrigid.o libtinker.a ; strip optrigid.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o path.x path.o libtinker.a ; strip path.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o pdbxyz.x pdbxyz.o libtinker.a ; strip pdbxyz.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o polarize.x polarize.o libtinker.a ; strip polarize.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o poledit.x poledit.o libtinker.a ; strip poledit.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o potential.x potential.o libtinker.a ; strip potential.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o prmedit.x prmedit.o libtinker.a ; strip prmedit.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o protein.x protein.o libtinker.a ; strip protein.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o pss.x pss.o libtinker.a ; strip pss.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o pssrigid.x pssrigid.o libtinker.a ; strip pssrigid.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o pssrot.x pssrot.o libtinker.a ; strip pssrot.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o radial.x radial.o libtinker.a ; strip radial.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o saddle.x saddle.o libtinker.a ; strip saddle.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o scan.x scan.o libtinker.a ; strip scan.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o sniffer.x sniffer.o libtinker.a ; strip sniffer.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o spacefill.x spacefill.o libtinker.a ; strip spacefill.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o spectrum.x spectrum.o libtinker.a ; strip spectrum.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o superpose.x superpose.o libtinker.a ; strip superpose.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o testgrad.x testgrad.o libtinker.a ; strip testgrad.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o testhess.x testhess.o libtinker.a ; strip testhess.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o testpair.x testpair.o libtinker.a ; strip testpair.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o testpol.x testpol.o libtinker.a ; strip testpol.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o testrot.x testrot.o libtinker.a ; strip testrot.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o testvir.x testvir.o libtinker.a ; strip testvir.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o timer.x timer.o libtinker.a ; strip timer.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o timerot.x timerot.o libtinker.a ; strip timerot.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o torsfit.x torsfit.o libtinker.a ; strip torsfit.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o valence.x valence.o libtinker.a ; strip valence.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o vibbig.x vibbig.o libtinker.a ; strip vibbig.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o vibrate.x vibrate.o libtinker.a ; strip vibrate.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o vibrot.x vibrot.o libtinker.a ; strip vibrot.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o xtalfit.x xtalfit.o libtinker.a ; strip xtalfit.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o xtalmin.x xtalmin.o libtinker.a ; strip xtalmin.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o xyzedit.x xyzedit.o libtinker.a ; strip xyzedit.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o xyzint.x xyzint.o libtinker.a ; strip xyzint.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o xyzmol2.x xyzmol2.o libtinker.a ; strip xyzmol2.x
ifort -O3 -no-ipo -no-prec-div -recursive -static-libgcc -static-intel -o xyzpdb.x xyzpdb.o libtinker.a ; strip xyzpdb.x

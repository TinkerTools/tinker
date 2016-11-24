#
#
#  ##############################################################
#  ##                                                          ##
#  ##  linkd.make  --  link the TINKER programs for debugging  ##
#  ##             (Intel Fortran for MacOS Version)            ##
#  ##                                                          ##
#  ##############################################################
#
#
ifort -g -Wl,-no_pie -o alchemy.x alchemy.o libtinker.a
ifort -g -Wl,-no_pie -o analyze.x analyze.o libtinker.a
ifort -g -Wl,-no_pie -o anneal.x anneal.o libtinker.a
ifort -g -Wl,-no_pie -o archive.x archive.o libtinker.a
ifort -g -Wl,-no_pie -o bar.x bar.o libtinker.a
ifort -g -Wl,-no_pie -o correlate.x correlate.o libtinker.a
ifort -g -Wl,-no_pie -o crystal.x crystal.o libtinker.a
ifort -g -Wl,-no_pie -o diffuse.x diffuse.o libtinker.a
ifort -g -Wl,-no_pie -o distgeom.x distgeom.o libtinker.a
ifort -g -Wl,-no_pie -o document.x document.o libtinker.a
ifort -g -Wl,-no_pie -o dynamic.x dynamic.o libtinker.a
ifort -g -Wl,-no_pie -o gda.x gda.o libtinker.a
ifort -g -Wl,-no_pie -o intedit.x intedit.o libtinker.a
ifort -g -Wl,-no_pie -o intxyz.x intxyz.o libtinker.a
ifort -g -Wl,-no_pie -o minimize.x minimize.o libtinker.a
ifort -g -Wl,-no_pie -o minirot.x minirot.o libtinker.a
ifort -g -Wl,-no_pie -o minrigid.x minrigid.o libtinker.a
ifort -g -Wl,-no_pie -o molxyz.x molxyz.o libtinker.a
ifort -g -Wl,-no_pie -o monte.x monte.o libtinker.a
ifort -g -Wl,-no_pie -o newton.x newton.o libtinker.a
ifort -g -Wl,-no_pie -o newtrot.x newtrot.o libtinker.a
ifort -g -Wl,-no_pie -o nucleic.x nucleic.o libtinker.a
ifort -g -Wl,-no_pie -o optimize.x optimize.o libtinker.a
ifort -g -Wl,-no_pie -o optirot.x optirot.o libtinker.a
ifort -g -Wl,-no_pie -o optrigid.x optrigid.o libtinker.a
ifort -g -Wl,-no_pie -o path.x path.o libtinker.a
ifort -g -Wl,-no_pie -o pdbxyz.x pdbxyz.o libtinker.a
ifort -g -Wl,-no_pie -o polarize.x polarize.o libtinker.a
ifort -g -Wl,-no_pie -o poledit.x poledit.o libtinker.a
ifort -g -Wl,-no_pie -o potential.x potential.o libtinker.a
ifort -g -Wl,-no_pie -o prmedit.x prmedit.o libtinker.a
ifort -g -Wl,-no_pie -o protein.x protein.o libtinker.a
ifort -g -Wl,-no_pie -o pss.x pss.o libtinker.a
ifort -g -Wl,-no_pie -o pssrigid.x pssrigid.o libtinker.a
ifort -g -Wl,-no_pie -o pssrot.x pssrot.o libtinker.a
ifort -g -Wl,-no_pie -o radial.x radial.o libtinker.a
ifort -g -Wl,-no_pie -o saddle.x saddle.o libtinker.a
ifort -g -Wl,-no_pie -o scan.x scan.o libtinker.a
ifort -g -Wl,-no_pie -o sniffer.x sniffer.o libtinker.a
ifort -g -Wl,-no_pie -o spacefill.x spacefill.o libtinker.a
ifort -g -Wl,-no_pie -o spectrum.x spectrum.o libtinker.a
ifort -g -Wl,-no_pie -o superpose.x superpose.o libtinker.a
ifort -g -Wl,-no_pie -o sybylxyz.x sybylxyz.o libtinker.a
ifort -g -Wl,-no_pie -o testgrad.x testgrad.o libtinker.a
ifort -g -Wl,-no_pie -o testhess.x testhess.o libtinker.a
ifort -g -Wl,-no_pie -o testpair.x testpair.o libtinker.a
ifort -g -Wl,-no_pie -o testpol.x testpol.o libtinker.a
ifort -g -Wl,-no_pie -o testrot.x testrot.o libtinker.a
ifort -g -Wl,-no_pie -o timer.x timer.o libtinker.a
ifort -g -Wl,-no_pie -o timerot.x timerot.o libtinker.a
ifort -g -Wl,-no_pie -o torsfit.x torsfit.o libtinker.a
ifort -g -Wl,-no_pie -o valence.x valence.o libtinker.a
ifort -g -Wl,-no_pie -o vibbig.x vibbig.o libtinker.a
ifort -g -Wl,-no_pie -o vibrate.x vibrate.o libtinker.a
ifort -g -Wl,-no_pie -o vibrot.x vibrot.o libtinker.a
ifort -g -Wl,-no_pie -o xtalfit.x xtalfit.o libtinker.a
ifort -g -Wl,-no_pie -o xtalmin.x xtalmin.o libtinker.a
ifort -g -Wl,-no_pie -o xyzedit.x xyzedit.o libtinker.a
ifort -g -Wl,-no_pie -o xyzint.x xyzint.o libtinker.a
ifort -g -Wl,-no_pie -o xyzpdb.x xyzpdb.o libtinker.a
ifort -g -Wl,-no_pie -o xyzsybyl.x xyzsybyl.o libtinker.a

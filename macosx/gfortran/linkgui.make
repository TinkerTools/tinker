#
#
#  ##################################################################
#  ##                                                              ##
#  ##  linkgui.make  --  link each of the TINKER programs for FFE  ##
#  ##                (gfortran for Mac OSX Version)                ##
#  ##                                                              ##
#  ##################################################################
#
#
gfortran -o alchemy.x alchemy.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o analyze.x analyze.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o anneal.x anneal.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o archive.x archive.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o correlate.x correlate.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o crystal.x crystal.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o diffuse.x diffuse.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o distgeom.x distgeom.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o document.x document.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o dynamic.x dynamic.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o gda.x gda.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o intedit.x intedit.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o intxyz.x intxyz.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o minimize.x minimize.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o minirot.x minirot.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o minrigid.x minrigid.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o molxyz.x molxyz.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o monte.x monte.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o newton.x newton.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o newtrot.x newtrot.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o nucleic.x nucleic.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o optimize.x optimize.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o optirot.x optirot.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o optrigid.x optrigid.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o path.x path.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o pdbxyz.x pdbxyz.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o polarize.x polarize.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o potential.x potential.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o prmedit.x prmedit.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o protein.x protein.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o pss.x pss.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o pssrigid.x pssrigid.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o pssrot.x pssrot.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o radial.x radial.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o saddle.x saddle.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o scan.x scan.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o sniffer.x sniffer.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o spacefill.x spacefill.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o spectrum.x spectrum.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o superpose.x superpose.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o sybylxyz.x sybylxyz.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o testgrad.x testgrad.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o testhess.x testhess.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o testpair.x testpair.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o testpol.x testpol.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o testrot.x testrot.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o timer.x timer.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o timerot.x timerot.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o torsfit.x torsfit.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o valence.x valence.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o vibbig.x vibbig.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o vibrate.x vibrate.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o vibrot.x vibrot.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o xtalfit.x xtalfit.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o xtalmin.x xtalmin.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o xyzedit.x xyzedit.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o xyzint.x xyzint.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o xyzpdb.x xyzpdb.o libtinker.a -framework JavaVM -framework CoreFoundation
gfortran -o xyzsybyl.x xyzsybyl.o libtinker.a -framework JavaVM -framework CoreFoundation

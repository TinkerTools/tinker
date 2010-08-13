#
#
#  ##################################################################
#  ##                                                              ##
#  ##  linkgui.make  --  link each of the TINKER programs for FFE  ##
#  ##                    (Mac OS X/g95 Version)                    ##
#  ##                                                              ##
#  ##################################################################
#
#
g95 -o alchemy.x alchemy.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o analyze.x analyze.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o anneal.x anneal.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o archive.x archive.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o correlate.x correlate.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o crystal.x crystal.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o diffuse.x diffuse.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o distgeom.x distgeom.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o document.x document.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o dynamic.x dynamic.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o gda.x gda.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o intedit.x intedit.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o intxyz.x intxyz.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o minimize.x minimize.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o minirot.x minirot.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o minrigid.x minrigid.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o monte.x monte.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o newton.x newton.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o newtrot.x newtrot.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o nucleic.x nucleic.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o optimize.x optimize.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o optirot.x optirot.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o optrigid.x optrigid.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o path.x path.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o pdbxyz.x pdbxyz.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o polarize.x polarize.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o potential.x potential.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o prmedit.x prmedit.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o protein.x protein.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o pss.x pss.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o pssrigid.x pssrigid.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o pssrot.x pssrot.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o radial.x radial.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o saddle.x saddle.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o scan.x scan.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o sniffer.x sniffer.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o spacefill.x spacefill.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o spectrum.x spectrum.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o superpose.x superpose.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o sybylxyz.x sybylxyz.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o testgrad.x testgrad.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o testhess.x testhess.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o testpair.x testpair.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o testrot.x testrot.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o timer.x timer.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o timerot.x timerot.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o torsfit.x torsfit.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o valence.x valence.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o vibbig.x vibbig.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o vibrate.x vibrate.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o vibrot.x vibrot.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o xtalfit.x xtalfit.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o xtalmin.x xtalmin.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o xyzedit.x xyzedit.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o xyzint.x xyzint.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o xyzpdb.x xyzpdb.o libtinker.a -framework JavaVM -framework CoreFoundation
g95 -o xyzsybyl.x xyzsybyl.o libtinker.a -framework JavaVM -framework CoreFoundation

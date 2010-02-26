#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##                  (IBM AIX Unix Version)                   ##
#  ##                                                           ##
#  ###############################################################
#
#
xlf -o alchemy.x alchemy.o -L. -ltinker
xlf -o analyze.x analyze.o -L. -ltinker
xlf -o anneal.x anneal.o -L. -ltinker
xlf -o archive.x archive.o -L. -ltinker
xlf -o correlate.x correlate.o -L. -ltinker
xlf -o crystal.x crystal.o -L. -ltinker
xlf -o diffuse.x diffuse.o -L. -ltinker
xlf -o distgeom.x distgeom.o -L. -ltinker
xlf -o document.x document.o -L. -ltinker
xlf -o dynamic.x dynamic.o -L. -ltinker
xlf -o gda.x gda.o -L. -ltinker
xlf -o intedit.x intedit.o -L. -ltinker
xlf -o intxyz.x intxyz.o -L. -ltinker
xlf -o minimize.x minimize.o -L. -ltinker
xlf -o minirot.x minirot.o -L. -ltinker
xlf -o minrigid.x minrigid.o -L. -ltinker
xlf -o monte.x monte.o -L. -ltinker
xlf -o newton.x newton.o -L. -ltinker
xlf -o newtrot.x newtrot.o -L. -ltinker
xlf -o nucleic.x nucleic.o -L. -ltinker
xlf -o optimize.x optimize.o -L. -ltinker
xlf -o optirot.x optirot.o -L. -ltinker
xlf -o optrigid.x optrigid.o -L. -ltinker
xlf -o path.x path.o -L. -ltinker
xlf -o pdbxyz.x pdbxyz.o -L. -ltinker
xlf -o polarize.x polarize.o -L. -ltinker
xlf -o poledit.x poledit.o -L. -ltinker
xlf -o potential.x potential.o -L. -ltinker
xlf -o prmedit.x prmedit.o -L. -ltinker
xlf -o protein.x protein.o -L. -ltinker
xlf -o pss.x pss.o -L. -ltinker
xlf -o pssrigid.x pssrigid.o -L. -ltinker
xlf -o pssrot.x pssrot.o -L. -ltinker
xlf -o radial.x radial.o -L. -ltinker
xlf -o saddle.x saddle.o -L. -ltinker
xlf -o scan.x scan.o -L. -ltinker
xlf -o sniffer.x sniffer.o -L. -ltinker
xlf -o spacefill.x spacefill.o -L. -ltinker
xlf -o spectrum.x spectrum.o -L. -ltinker
xlf -o superpose.x superpose.o -L. -ltinker
xlf -o sybylxyz.x sybylxyz.o -L. -ltinker
xlf -o testgrad.x testgrad.o -L. -ltinker
xlf -o testhess.x testhess.o -L. -ltinker
xlf -o testpair.x testpair.o -L. -ltinker
xlf -o testrot.x testrot.o -L. -ltinker
xlf -o timer.x timer.o -L. -ltinker
xlf -o timerot.x timerot.o -L. -ltinker
xlf -o valence.x valence.o -L. -ltinker
xlf -o vibbig.x vibbig.o -L. -ltinker
xlf -o vibrate.x vibrate.o -L. -ltinker
xlf -o vibrot.x vibrot.o -L. -ltinker
xlf -o xtalfit.x xtalfit.o -L. -ltinker
xlf -o xtalmin.x xtalmin.o -L. -ltinker
xlf -o xyzedit.x xyzedit.o -L. -ltinker
xlf -o xyzint.x xyzint.o -L. -ltinker
xlf -o xyzpdb.x xyzpdb.o -L. -ltinker
xlf -o xyzsybyl.x xyzsybyl.o -L. -ltinker

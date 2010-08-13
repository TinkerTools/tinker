#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##              (Hewlett-Packard HP-UX Version)              ##
#  ##                                                           ##
#  ###############################################################
#
#
fort77 -O2 +E1 +U77 -o alchemy.x alchemy.o -L. -ltinker
fort77 -O2 +E1 +U77 -o analyze.x analyze.o -L. -ltinker
fort77 -O2 +E1 +U77 -o anneal.x anneal.o -L. -ltinker
fort77 -O2 +E1 +U77 -o archive.x archive.o -L. -ltinker
fort77 -O2 +E1 +U77 -o correlate.x correlate.o -L. -ltinker
fort77 -O2 +E1 +U77 -o crystal.x crystal.o -L. -ltinker
fort77 -O2 +E1 +U77 -o diffuse.x distgeom.o -L. -ltinker
fort77 -O2 +E1 +U77 -o distgeom.x distgeom.o -L. -ltinker
fort77 -O2 +E1 +U77 -o document.x document.o -L. -ltinker
fort77 -O2 +E1 +U77 -o dynamic.x dynamic.o -L. -ltinker
fort77 -O2 +E1 +U77 -o gda.x gda.o -L. -ltinker
fort77 -O2 +E1 +U77 -o intedit.x intedit.o -L. -ltinker
fort77 -O2 +E1 +U77 -o intxyz.x intxyz.o -L. -ltinker
fort77 -O2 +E1 +U77 -o minimize.x minimize.o -L. -ltinker
fort77 -O2 +E1 +U77 -o minirot.x minirot.o -L. -ltinker
fort77 -O2 +E1 +U77 -o minrigid.x minrigid.o -L. -ltinker
fort77 -O2 +E1 +U77 -o newton.x newton.o -L. -ltinker
fort77 -O2 +E1 +U77 -o newtrot.x newtrot.o -L. -ltinker
fort77 -O2 +E1 +U77 -o nucleic.x nucleic.o -L. -ltinker
fort77 -O2 +E1 +U77 -o optimize.x optimize.o -L. -ltinker
fort77 -O2 +E1 +U77 -o optirot.x optirot.o -L. -ltinker
fort77 -O2 +E1 +U77 -o optrigid.x optrigid.o -L. -ltinker
fort77 -O2 +E1 +U77 -o path.x path.o -L. -ltinker
fort77 -O2 +E1 +U77 -o pdbxyz.x pdbxyz.o -L. -ltinker
fort77 -O2 +E1 +U77 -o polarize.x polarize.o -L. -ltinker
fort77 -O2 +E1 +U77 -o poledit.x poledit.o -L. -ltinker
fort77 -O2 +E1 +U77 -o potential.x potential.o -L. -ltinker
fort77 -O2 +E1 +U77 -o prmedit.x prmedit.o -L. -ltinker
fort77 -O2 +E1 +U77 -o protein.x protein.o -L. -ltinker
fort77 -O2 +E1 +U77 -o pss.x pss.o -L. -ltinker
fort77 -O2 +E1 +U77 -o pssrigid.x pssrigid.o -L. -ltinker
fort77 -O2 +E1 +U77 -o pssrot.x pssrot.o -L. -ltinker
fort77 -O2 +E1 +U77 -o radial.x radial.o -L. -ltinker
fort77 -O2 +E1 +U77 -o saddle.x saddle.o -L. -ltinker
fort77 -O2 +E1 +U77 -o scan.x scan.o -L. -ltinker
fort77 -O2 +E1 +U77 -o sniffer.x sniffer.o -L. -ltinker
fort77 -O2 +E1 +U77 -o spacefill.x spacefill.o -L. -ltinker
fort77 -O2 +E1 +U77 -o spectrum.x spectrum.o -L. -ltinker
fort77 -O2 +E1 +U77 -o superpose.x superpose.o -L. -ltinker
fort77 -O2 +E1 +U77 -o sybylxyz.x sybylxyz.o -L. -ltinker
fort77 -O2 +E1 +U77 -o testgrad.x testgrad.o -L. -ltinker
fort77 -O2 +E1 +U77 -o testhess.x testhess.o -L. -ltinker
fort77 -O2 +E1 +U77 -o testpair.x testpair.o -L. -ltinker
fort77 -O2 +E1 +U77 -o testrot.x testrot.o -L. -ltinker
fort77 -O2 +E1 +U77 -o timer.x timer.o -L. -ltinker
fort77 -O2 +E1 +U77 -o timerot.x timerot.o -L. -ltinker
fort77 -O2 +E1 +U77 -o torsfit.x torsfit.o -L. -ltinker
fort77 -O2 +E1 +U77 -o valence.x valence.o -L. -ltinker
fort77 -O2 +E1 +U77 -o vibbig.x vibbig.o -L. -ltinker
fort77 -O2 +E1 +U77 -o vibrate.x vibrate.o -L. -ltinker
fort77 -O2 +E1 +U77 -o vibrot.x vibrot.o -L. -ltinker
fort77 -O2 +E1 +U77 -o xtalfit.x xtalfit.o -L. -ltinker
fort77 -O2 +E1 +U77 -o xtalmin.x xtalmin.o -L. -ltinker
fort77 -O2 +E1 +U77 -o xyzedit.x xyzedit.o -L. -ltinker
fort77 -O2 +E1 +U77 -o xyzint.x xyzint.o -L. -ltinker
fort77 -O2 +E1 +U77 -o xyzpdb.x xyzpdb.o -L. -ltinker
fort77 -O2 +E1 +U77 -o xyzsybyl.x xyzsybyl.o -L. -ltinker

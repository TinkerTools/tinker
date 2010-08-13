#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##                   (SGI/GNU g77Version)                   ##
#  ##                                                           ##
#  ###############################################################
#
#
g77 -o alchemy.x alchemy.o libtinker.a
g77 -o analyze.x analyze.o libtinker.a
g77 -o anneal.x anneal.o libtinker.a
g77 -o archive.x archive.o libtinker.a
g77 -o correlate.x correlate.o libtinker.a
g77 -o crystal.x crystal.o libtinker.a
g77 -o diffuse.x diffuse.o libtinker.a
g77 -o distgeom.x distgeom.o libtinker.a
g77 -o document.x document.o libtinker.a
g77 -o dynamic.x dynamic.o libtinker.a
g77 -o gda.x gda.o libtinker.a
g77 -o intedit.x intedit.o libtinker.a
g77 -o intxyz.x intxyz.o libtinker.a
g77 -o minimize.x minimize.o libtinker.a
g77 -o minirot.x minirot.o libtinker.a
g77 -o minrigid.x minrigid.o libtinker.a
g77 -o monte.x monte.o libtinker.a
g77 -o newton.x newton.o libtinker.a
g77 -o newtrot.x newtrot.o libtinker.a
g77 -o nucleic.x nucleic.o libtinker.a
g77 -o optimize.x optimize.o libtinker.a
g77 -o optirot.x optirot.o libtinker.a
g77 -o optrigid.x optrigid.o libtinker.a
g77 -o path.x path.o libtinker.a
g77 -o pdbxyz.x pdbxyz.o libtinker.a
g77 -o polarize.x polarize.o libtinker.a
g77 -o poledit.x poledit.o libtinker.a
g77 -o potential.x potential.o libtinker.a
g77 -o prmedit.x prmedit.o libtinker.a
g77 -o protein.x protein.o libtinker.a
g77 -o pss.x pss.o libtinker.a
#g77-o pssrigid.x pssrigid.o libtinker.a
g77 -o pssrot.x pssrot.o libtinker.a
g77 -o radial.x radial.o libtinker.a
g77 -o saddle.x saddle.o libtinker.a
g77 -o scan.x scan.o libtinker.a
g77 -o sniffer.x sniffer.o libtinker.a
g77 -o spacefill.x spacefill.o libtinker.a
g77 -o spectrum.x spectrum.o libtinker.a
g77 -o superpose.x superpose.o libtinker.a
g77 -o sybylxyz.x sybylxyz.o libtinker.a
g77 -o testgrad.x testgrad.o libtinker.a
g77 -o testhess.x testhess.o libtinker.a
g77 -o testpair.x testpair.o libtinker.a
g77 -o testrot.x testrot.o libtinker.a
g77 -o timer.x timer.o libtinker.a
g77 -o timerot.x timerot.o libtinker.a
g77 -o torsfit.x torsfit.o libtinker.a
g77 -o valence.x valence.o libtinker.a
g77 -o vibbig.x vibbig.o libtinker.a
g77 -o vibrate.x vibrate.o libtinker.a
g77 -o vibrot.x vibrot.o libtinker.a
g77 -o xtalfit.x xtalfit.o libtinker.a
g77 -o xtalmin.x xtalmin.o libtinker.a
g77 -o xyzedit.x xyzedit.o libtinker.a
g77 -o xyzint.x xyzint.o libtinker.a
g77 -o xyzpdb.x xyzpdb.o libtinker.a
g77 -o xyzsybyl.x xyzsybyl.o libtinker.a

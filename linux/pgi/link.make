#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##              (PGI Fortran for Linux Version)              ##
#  ##                                                           ##
#  ###############################################################
#
#
pgf77 -s -o alchemy.x alchemy.o libtinker.a
pgf77 -s -o analyze.x analyze.o libtinker.a
pgf77 -s -o anneal.x anneal.o libtinker.a
pgf77 -s -o archive.x archive.o libtinker.a
pgf77 -s -o correlate.x correlate.o libtinker.a
pgf77 -s -o crystal.x crystal.o libtinker.a
pgf77 -s -o diffuse.x diffuse.o libtinker.a
pgf77 -s -o distgeom.x distgeom.o libtinker.a
pgf77 -s -o document.x document.o libtinker.a
pgf77 -s -o dynamic.x dynamic.o libtinker.a
pgf77 -s -o gda.x gda.o libtinker.a
pgf77 -s -o intedit.x intedit.o libtinker.a
pgf77 -s -o intxyz.x intxyz.o libtinker.a
pgf77 -s -o minimize.x minimize.o libtinker.a
pgf77 -s -o minirot.x minirot.o libtinker.a
pgf77 -s -o minrigid.x minrigid.o libtinker.a
pgf77 -s -o molxyz.x molxyz.o libtinker.a
pgf77 -s -o monte.x monte.o libtinker.a
pgf77 -s -o newton.x newton.o libtinker.a
pgf77 -s -o newtrot.x newtrot.o libtinker.a
pgf77 -s -o nucleic.x nucleic.o libtinker.a
pgf77 -s -o optimize.x optimize.o libtinker.a
pgf77 -s -o optirot.x optirot.o libtinker.a
pgf77 -s -o optrigid.x optrigid.o libtinker.a
pgf77 -s -o path.x path.o libtinker.a
pgf77 -s -o pdbxyz.x pdbxyz.o libtinker.a
pgf77 -s -o polarize.x polarize.o libtinker.a
pgf77 -s -o poledit.x poledit.o libtinker.a
pgf77 -s -o potential.x potential.o libtinker.a
pgf77 -s -o prmedit.x prmedit.o libtinker.a
pgf77 -s -o protein.x protein.o libtinker.a
pgf77 -s -o pss.x pss.o libtinker.a
pgf77 -s -o pssrigid.x pssrigid.o libtinker.a
pgf77 -s -o pssrot.x pssrot.o libtinker.a
pgf77 -s -o radial.x radial.o libtinker.a
pgf77 -s -o saddle.x saddle.o libtinker.a
pgf77 -s -o scan.x scan.o libtinker.a
pgf77 -s -o sniffer.x sniffer.o libtinker.a
pgf77 -s -o spacefill.x spacefill.o libtinker.a
pgf77 -s -o spectrum.x spectrum.o libtinker.a
pgf77 -s -o superpose.x superpose.o libtinker.a
pgf77 -s -o sybylxyz.x sybylxyz.o libtinker.a
pgf77 -s -o testgrad.x testgrad.o libtinker.a
pgf77 -s -o testhess.x testhess.o libtinker.a
pgf77 -s -o testpair.x testpair.o libtinker.a
pgf77 -s -o testrot.x testrot.o libtinker.a
pgf77 -s -o timer.x timer.o libtinker.a
pgf77 -s -o timerot.x timerot.o libtinker.a
pgf77 -s -o torsfit.x torsfit.o libtinker.a
pgf77 -s -o valence.x valence.o libtinker.a
pgf77 -s -o vibbig.x vibbig.o libtinker.a
pgf77 -s -o vibrate.x vibrate.o libtinker.a
pgf77 -s -o vibrot.x vibrot.o libtinker.a
pgf77 -s -o xtalfit.x xtalfit.o libtinker.a
pgf77 -s -o xtalmin.x xtalmin.o libtinker.a
pgf77 -s -o xyzedit.x xyzedit.o libtinker.a
pgf77 -s -o xyzint.x xyzint.o libtinker.a
pgf77 -s -o xyzpdb.x xyzpdb.o libtinker.a
pgf77 -s -o xyzsybyl.x xyzsybyl.o libtinker.a

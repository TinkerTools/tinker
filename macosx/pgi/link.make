#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##              (PGI pgf95 for Mac OSX Version)              ##
#  ##                                                           ##
#  ###############################################################
#
#
pgf95 -fast -mp -o alchemy.x alchemy.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o analyze.x analyze.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o anneal.x anneal.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o archive.x archive.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o correlate.x correlate.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o crystal.x crystal.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o diffuse.x diffuse.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o distgeom.x distgeom.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o document.x document.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o dynamic.x dynamic.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o gda.x gda.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o intedit.x intedit.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o intxyz.x intxyz.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o minimize.x minimize.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o minirot.x minirot.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o minrigid.x minrigid.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o molxyz.x molxyz.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o monte.x monte.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o newton.x newton.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o newtrot.x newtrot.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o nucleic.x nucleic.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o optimize.x optimize.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o optirot.x optirot.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o optrigid.x optrigid.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o path.x path.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o polarize.x polarize.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o poledit.x poledit.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o potential.x potential.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o prmedit.x prmedit.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o protein.x protein.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o pss.x pss.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o pssrigid.x pssrigid.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o pssrot.x pssrot.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o radial.x radial.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o saddle.x saddle.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o scan.x scan.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o sniffer.x sniffer.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o spacefill.x spacefill.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o spectrum.x spectrum.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o superpose.x superpose.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o sybylxyz.x sybylxyz.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o testgrad.x testgrad.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o testhess.x testhess.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o testpair.x testpair.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o testrot.x testrot.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o timer.x timer.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o timerot.x timerot.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o torsfit.x torsfit.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o valence.x valence.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o vibbig.x vibbig.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o vibrate.x vibrate.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o vibrot.x vibrot.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o xtalfit.x xtalfit.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o xtalmin.x xtalmin.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o xyzedit.x xyzedit.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o xyzint.x xyzint.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_omp.a libfftw3.a
pgf95 -fast -mp -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3_omp.a libfftw3.a

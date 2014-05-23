#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##             (GNU gfortran for Mac OSX Version)            ##
#  ##                                                           ##
#  ###############################################################
#
#
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o alchemy.x alchemy.o libtinker.a libfftw3_threads.a libfftw3.a ; strip alchemy.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o analyze.x analyze.o libtinker.a libfftw3_threads.a libfftw3.a ; strip analyze.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o anneal.x anneal.o libtinker.a libfftw3_threads.a libfftw3.a ; strip anneal.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o archive.x archive.o libtinker.a libfftw3_threads.a libfftw3.a ; strip archive.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o bar.x bar.o libtinker.a libfftw3_threads.a libfftw3.a ; strip bar.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o correlate.x correlate.o libtinker.a libfftw3_threads.a libfftw3.a ; strip correlate.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o crystal.x crystal.o libtinker.a libfftw3_threads.a libfftw3.a ; strip crystal.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o diffuse.x diffuse.o libtinker.a libfftw3_threads.a libfftw3.a ; strip diffuse.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o distgeom.x distgeom.o libtinker.a libfftw3_threads.a libfftw3.a ; strip distgeom.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o document.x document.o libtinker.a libfftw3_threads.a libfftw3.a ; strip document.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o dynamic.x dynamic.o libtinker.a libfftw3_threads.a libfftw3.a ; strip dynamic.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o gda.x gda.o libtinker.a libfftw3_threads.a libfftw3.a ; strip gda.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o intedit.x intedit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip intedit.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o intxyz.x intxyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip intxyz.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o minimize.x minimize.o libtinker.a libfftw3_threads.a libfftw3.a ; strip minimize.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o minirot.x minirot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip minirot.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o minrigid.x minrigid.o libtinker.a libfftw3_threads.a libfftw3.a ; strip minrigid.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o molxyz.x molxyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip molxyz.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o monte.x monte.o libtinker.a libfftw3_threads.a libfftw3.a ; strip monte.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o newton.x newton.o libtinker.a libfftw3_threads.a libfftw3.a ; strip newton.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o newtrot.x newtrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip newtrot.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o nucleic.x nucleic.o libtinker.a libfftw3_threads.a libfftw3.a ; strip nucleic.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o optimize.x optimize.o libtinker.a libfftw3_threads.a libfftw3.a ; strip optimize.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o optirot.x optirot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip optirot.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o optrigid.x optrigid.o libtinker.a libfftw3_threads.a libfftw3.a ; strip optrigid.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o path.x path.o libtinker.a libfftw3_threads.a libfftw3.a ; strip path.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o pdbxyz.x pdbxyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pdbxyz.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o polarize.x polarize.o libtinker.a libfftw3_threads.a libfftw3.a ; strip polarize.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o poledit.x poledit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip poledit.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o potential.x potential.o libtinker.a libfftw3_threads.a libfftw3.a ; strip potential.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o prmedit.x prmedit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip prmedit.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o protein.x protein.o libtinker.a libfftw3_threads.a libfftw3.a ; strip protein.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o pss.x pss.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pss.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o pssrigid.x pssrigid.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pssrigid.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o pssrot.x pssrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip pssrot.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o radial.x radial.o libtinker.a libfftw3_threads.a libfftw3.a ; strip radial.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o saddle.x saddle.o libtinker.a libfftw3_threads.a libfftw3.a ; strip saddle.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o scan.x scan.o libtinker.a libfftw3_threads.a libfftw3.a ; strip scan.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o sniffer.x sniffer.o libtinker.a libfftw3_threads.a libfftw3.a ; strip sniffer.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o spacefill.x spacefill.o libtinker.a libfftw3_threads.a libfftw3.a ; strip spacefill.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o spectrum.x spectrum.o libtinker.a libfftw3_threads.a libfftw3.a ; strip spectrum.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o superpose.x superpose.o libtinker.a libfftw3_threads.a libfftw3.a ; strip superpose.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o sybylxyz.x sybylxyz.o libtinker.a libfftw3_threads.a libfftw3.a ; strip sybylxyz.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o testgrad.x testgrad.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testgrad.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o testhess.x testhess.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testhess.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o testpair.x testpair.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testpair.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o testpol.x testpol.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testpol.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o testrot.x testrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip testrot.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o timer.x timer.o libtinker.a libfftw3_threads.a libfftw3.a ; strip timer.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o timerot.x timerot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip timerot.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o torsfit.x torsfit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip torsfit.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o valence.x valence.o libtinker.a libfftw3_threads.a libfftw3.a ; strip valence.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o vibbig.x vibbig.o libtinker.a libfftw3_threads.a libfftw3.a ; strip vibbig.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o vibrate.x vibrate.o libtinker.a libfftw3_threads.a libfftw3.a ; strip vibrate.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o vibrot.x vibrot.o libtinker.a libfftw3_threads.a libfftw3.a ; strip vibrot.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o xtalfit.x xtalfit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xtalfit.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o xtalmin.x xtalmin.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xtalmin.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o xyzedit.x xyzedit.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzedit.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o xyzint.x xyzint.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzint.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o xyzpdb.x xyzpdb.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzpdb.x
gfortran -O3 -ffast-math -fopenmp -static-libgcc -o xyzsybyl.x xyzsybyl.o libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzsybyl.x

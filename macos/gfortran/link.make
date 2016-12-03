#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##              (GNU gfortran for MacOS Version)             ##
#  ##                                                           ##
#  ###############################################################
#
#
#  if the -static-libgcc flag is not enforced, then copy static
#  versions of the main gcc libraries (libgcc.a, libgfortran.a,
#  libgomp.a and libquadmath.a) into the present build directory;
#  with Linux, all four libraries must usually be present locally,
#  while for MacOS, libquadmath.a is the only one usually needed
#
#
gfortran -Ofast -fopenmp -static-libgcc -o alchemy.x alchemy.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip alchemy.x
gfortran -Ofast -fopenmp -static-libgcc -o analyze.x analyze.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip analyze.x
gfortran -Ofast -fopenmp -static-libgcc -o anneal.x anneal.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip anneal.x
gfortran -Ofast -fopenmp -static-libgcc -o archive.x archive.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip archive.x
gfortran -Ofast -fopenmp -static-libgcc -o bar.x bar.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip bar.x
gfortran -Ofast -fopenmp -static-libgcc -o correlate.x correlate.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip correlate.x
gfortran -Ofast -fopenmp -static-libgcc -o crystal.x crystal.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip crystal.x
gfortran -Ofast -fopenmp -static-libgcc -o diffuse.x diffuse.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip diffuse.x
gfortran -Ofast -fopenmp -static-libgcc -o distgeom.x distgeom.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip distgeom.x
gfortran -Ofast -fopenmp -static-libgcc -o document.x document.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip document.x
gfortran -Ofast -fopenmp -static-libgcc -o dynamic.x dynamic.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip dynamic.x
gfortran -Ofast -fopenmp -static-libgcc -o gda.x gda.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip gda.x
gfortran -Ofast -fopenmp -static-libgcc -o intedit.x intedit.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip intedit.x
gfortran -Ofast -fopenmp -static-libgcc -o intxyz.x intxyz.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip intxyz.x
gfortran -Ofast -fopenmp -static-libgcc -o minimize.x minimize.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip minimize.x
gfortran -Ofast -fopenmp -static-libgcc -o minirot.x minirot.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip minirot.x
gfortran -Ofast -fopenmp -static-libgcc -o minrigid.x minrigid.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip minrigid.x
gfortran -Ofast -fopenmp -static-libgcc -o molxyz.x molxyz.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip molxyz.x
gfortran -Ofast -fopenmp -static-libgcc -o monte.x monte.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip monte.x
gfortran -Ofast -fopenmp -static-libgcc -o newton.x newton.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip newton.x
gfortran -Ofast -fopenmp -static-libgcc -o newtrot.x newtrot.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip newtrot.x
gfortran -Ofast -fopenmp -static-libgcc -o nucleic.x nucleic.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip nucleic.x
gfortran -Ofast -fopenmp -static-libgcc -o optimize.x optimize.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip optimize.x
gfortran -Ofast -fopenmp -static-libgcc -o optirot.x optirot.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip optirot.x
gfortran -Ofast -fopenmp -static-libgcc -o optrigid.x optrigid.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip optrigid.x
gfortran -Ofast -fopenmp -static-libgcc -o path.x path.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip path.x
gfortran -Ofast -fopenmp -static-libgcc -o pdbxyz.x pdbxyz.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip pdbxyz.x
gfortran -Ofast -fopenmp -static-libgcc -o polarize.x polarize.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip polarize.x
gfortran -Ofast -fopenmp -static-libgcc -o poledit.x poledit.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip poledit.x
gfortran -Ofast -fopenmp -static-libgcc -o potential.x potential.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip potential.x
gfortran -Ofast -fopenmp -static-libgcc -o prmedit.x prmedit.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip prmedit.x
gfortran -Ofast -fopenmp -static-libgcc -o protein.x protein.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip protein.x
gfortran -Ofast -fopenmp -static-libgcc -o pss.x pss.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip pss.x
gfortran -Ofast -fopenmp -static-libgcc -o pssrigid.x pssrigid.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip pssrigid.x
gfortran -Ofast -fopenmp -static-libgcc -o pssrot.x pssrot.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip pssrot.x
gfortran -Ofast -fopenmp -static-libgcc -o radial.x radial.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip radial.x
gfortran -Ofast -fopenmp -static-libgcc -o saddle.x saddle.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip saddle.x
gfortran -Ofast -fopenmp -static-libgcc -o scan.x scan.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip scan.x
gfortran -Ofast -fopenmp -static-libgcc -o sniffer.x sniffer.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip sniffer.x
gfortran -Ofast -fopenmp -static-libgcc -o spacefill.x spacefill.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip spacefill.x
gfortran -Ofast -fopenmp -static-libgcc -o spectrum.x spectrum.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip spectrum.x
gfortran -Ofast -fopenmp -static-libgcc -o superpose.x superpose.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip superpose.x
gfortran -Ofast -fopenmp -static-libgcc -o sybylxyz.x sybylxyz.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip sybylxyz.x
gfortran -Ofast -fopenmp -static-libgcc -o testgrad.x testgrad.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip testgrad.x
gfortran -Ofast -fopenmp -static-libgcc -o testhess.x testhess.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip testhess.x
gfortran -Ofast -fopenmp -static-libgcc -o testpair.x testpair.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip testpair.x
gfortran -Ofast -fopenmp -static-libgcc -o testpol.x testpol.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip testpol.x
gfortran -Ofast -fopenmp -static-libgcc -o testrot.x testrot.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip testrot.x
gfortran -Ofast -fopenmp -static-libgcc -o timer.x timer.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip timer.x
gfortran -Ofast -fopenmp -static-libgcc -o timerot.x timerot.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip timerot.x
gfortran -Ofast -fopenmp -static-libgcc -o torsfit.x torsfit.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip torsfit.x
gfortran -Ofast -fopenmp -static-libgcc -o valence.x valence.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip valence.x
gfortran -Ofast -fopenmp -static-libgcc -o vibbig.x vibbig.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip vibbig.x
gfortran -Ofast -fopenmp -static-libgcc -o vibrate.x vibrate.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip vibrate.x
gfortran -Ofast -fopenmp -static-libgcc -o vibrot.x vibrot.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip vibrot.x
gfortran -Ofast -fopenmp -static-libgcc -o xtalfit.x xtalfit.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip xtalfit.x
gfortran -Ofast -fopenmp -static-libgcc -o xtalmin.x xtalmin.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip xtalmin.x
gfortran -Ofast -fopenmp -static-libgcc -o xyzedit.x xyzedit.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzedit.x
gfortran -Ofast -fopenmp -static-libgcc -o xyzint.x xyzint.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzint.x
gfortran -Ofast -fopenmp -static-libgcc -o xyzpdb.x xyzpdb.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzpdb.x
gfortran -Ofast -fopenmp -static-libgcc -o xyzsybyl.x xyzsybyl.o -L. libtinker.a libfftw3_threads.a libfftw3.a ; strip xyzsybyl.x

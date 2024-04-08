#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the Tinker package programs  ##
#  ##              (GNU gfortran for macOS Version)             ##
#  ##                                                           ##
#  ###############################################################
#
#
#  if the -static-libgcc flag is not enforced, then copy static
#  versions of the main gcc libraries (libgcc.a, libgfortran.a,
#  libgomp.a and libquadmath.a) into the present build directory;
#  with Linux, all four libraries must be present locally, while
#  for macOS, libquadmath.a is the only one needed
#
#
gfortran -Ofast -fopenmp -static-libgcc -o alchemy.x alchemy.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip alchemy.x
gfortran -Ofast -fopenmp -static-libgcc -o analyze.x analyze.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip analyze.x
gfortran -Ofast -fopenmp -static-libgcc -o anneal.x anneal.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip anneal.x
gfortran -Ofast -fopenmp -static-libgcc -o arcedit.x arcedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip arcedit.x
gfortran -Ofast -fopenmp -static-libgcc -o bar.x bar.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip bar.x
gfortran -Ofast -fopenmp -static-libgcc -o correlate.x correlate.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip correlate.x
gfortran -Ofast -fopenmp -static-libgcc -o critical.x critical.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip critical.x
gfortran -Ofast -fopenmp -static-libgcc -o crystal.x crystal.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip crystal.x
gfortran -Ofast -fopenmp -static-libgcc -o diffuse.x diffuse.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip diffuse.x
gfortran -Ofast -fopenmp -static-libgcc -o distgeom.x distgeom.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip distgeom.x
gfortran -Ofast -fopenmp -static-libgcc -o document.x document.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip document.x
gfortran -Ofast -fopenmp -static-libgcc -o dynamic.x dynamic.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip dynamic.x
gfortran -Ofast -fopenmp -static-libgcc -o freefix.x freefix.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip freefix.x
gfortran -Ofast -fopenmp -static-libgcc -o gda.x gda.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip gda.x
gfortran -Ofast -fopenmp -static-libgcc -o intedit.x intedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip intedit.x
gfortran -Ofast -fopenmp -static-libgcc -o intxyz.x intxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip intxyz.x
gfortran -Ofast -fopenmp -static-libgcc -o minimize.x minimize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minimize.x
gfortran -Ofast -fopenmp -static-libgcc -o minirot.x minirot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minirot.x
gfortran -Ofast -fopenmp -static-libgcc -o minrigid.x minrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minrigid.x
gfortran -Ofast -fopenmp -static-libgcc -o mol2xyz.x mol2xyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip mol2xyz.x
gfortran -Ofast -fopenmp -static-libgcc -o molxyz.x molxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip molxyz.x
gfortran -Ofast -fopenmp -static-libgcc -o monte.x monte.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip monte.x
gfortran -Ofast -fopenmp -static-libgcc -o newton.x newton.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip newton.x
gfortran -Ofast -fopenmp -static-libgcc -o newtrot.x newtrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip newtrot.x
gfortran -Ofast -fopenmp -static-libgcc -o nucleic.x nucleic.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip nucleic.x
gfortran -Ofast -fopenmp -static-libgcc -o optimize.x optimize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optimize.x
gfortran -Ofast -fopenmp -static-libgcc -o optirot.x optirot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optirot.x
gfortran -Ofast -fopenmp -static-libgcc -o optrigid.x optrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optrigid.x
gfortran -Ofast -fopenmp -static-libgcc -o path.x path.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip path.x
gfortran -Ofast -fopenmp -static-libgcc -o pdbxyz.x pdbxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pdbxyz.x
gfortran -Ofast -fopenmp -static-libgcc -o polarize.x polarize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip polarize.x
gfortran -Ofast -fopenmp -static-libgcc -o poledit.x poledit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip poledit.x
gfortran -Ofast -fopenmp -static-libgcc -o potential.x potential.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip potential.x
gfortran -Ofast -fopenmp -static-libgcc -o prmedit.x prmedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip prmedit.x
gfortran -Ofast -fopenmp -static-libgcc -o protein.x protein.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip protein.x
gfortran -Ofast -fopenmp -static-libgcc -o pss.x pss.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pss.x
gfortran -Ofast -fopenmp -static-libgcc -o pssrigid.x pssrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pssrigid.x
gfortran -Ofast -fopenmp -static-libgcc -o pssrot.x pssrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pssrot.x
gfortran -Ofast -fopenmp -static-libgcc -o radial.x radial.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip radial.x
gfortran -Ofast -fopenmp -static-libgcc -o saddle.x saddle.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip saddle.x
gfortran -Ofast -fopenmp -static-libgcc -o scan.x scan.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip scan.x
gfortran -Ofast -fopenmp -static-libgcc -o sniffer.x sniffer.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip sniffer.x
gfortran -Ofast -fopenmp -static-libgcc -o spacefill.x spacefill.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip spacefill.x
gfortran -Ofast -fopenmp -static-libgcc -o spectrum.x spectrum.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip spectrum.x
gfortran -Ofast -fopenmp -static-libgcc -o superpose.x superpose.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip superpose.x
gfortran -Ofast -fopenmp -static-libgcc -o testgrad.x testgrad.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testgrad.x
gfortran -Ofast -fopenmp -static-libgcc -o testhess.x testhess.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testhess.x
gfortran -Ofast -fopenmp -static-libgcc -o testpair.x testpair.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testpair.x
gfortran -Ofast -fopenmp -static-libgcc -o testpol.x testpol.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testpol.x
gfortran -Ofast -fopenmp -static-libgcc -o testrot.x testrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testrot.x
gfortran -Ofast -fopenmp -static-libgcc -o testsurf.x testsurf.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testsurf.x
gfortran -Ofast -fopenmp -static-libgcc -o testvir.x testvir.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testvir.x
gfortran -Ofast -fopenmp -static-libgcc -o timer.x timer.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip timer.x
gfortran -Ofast -fopenmp -static-libgcc -o timerot.x timerot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip timerot.x
gfortran -Ofast -fopenmp -static-libgcc -o torsfit.x torsfit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip torsfit.x
gfortran -Ofast -fopenmp -static-libgcc -o valence.x valence.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip valence.x
gfortran -Ofast -fopenmp -static-libgcc -o vibbig.x vibbig.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibbig.x
gfortran -Ofast -fopenmp -static-libgcc -o vibrate.x vibrate.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibrate.x
gfortran -Ofast -fopenmp -static-libgcc -o vibrot.x vibrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibrot.x
gfortran -Ofast -fopenmp -static-libgcc -o xtalfit.x xtalfit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xtalfit.x
gfortran -Ofast -fopenmp -static-libgcc -o xtalmin.x xtalmin.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xtalmin.x
gfortran -Ofast -fopenmp -static-libgcc -o xyzedit.x xyzedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzedit.x
gfortran -Ofast -fopenmp -static-libgcc -o xyzint.x xyzint.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzint.x
gfortran -Ofast -fopenmp -static-libgcc -o xyzmol2.x xyzmol2.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzmol2.x
gfortran -Ofast -fopenmp -static-libgcc -o xyzpdb.x xyzpdb.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzpdb.x

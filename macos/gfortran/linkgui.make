#
#
#  ##################################################################
#  ##                                                              ##
#  ##  linkgui.make  --  link each of the Tinker programs for FFE  ##
#  ##               (GNU gfortran for MacOS Version)               ##
#  ##                                                              ##
#  ##################################################################
#
#
gfortran -Ofast -fopenmp -static-libgcc -o alchemy.x alchemy.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o analyze.x analyze.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o anneal.x anneal.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o archive.x archive.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o bar.x bar.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o correlate.x correlate.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o crystal.x crystal.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o diffuse.x diffuse.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o distgeom.x distgeom.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o document.x document.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o dynamic.x dynamic.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o gda.x gda.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o intedit.x intedit.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o intxyz.x intxyz.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o minimize.x minimize.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o minirot.x minirot.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o minrigid.x minrigid.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o mol2xyz.x mol2xyz.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o molxyz.x molxyz.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o monte.x monte.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o newton.x newton.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o newtrot.x newtrot.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o nucleic.x nucleic.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o optimize.x optimize.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o optirot.x optirot.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o optrigid.x optrigid.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o path.x path.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o pdbxyz.x pdbxyz.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o polarize.x polarize.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o potential.x potential.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o prmedit.x prmedit.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o protein.x protein.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o pss.x pss.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o pssrigid.x pssrigid.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o pssrot.x pssrot.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o radial.x radial.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o saddle.x saddle.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o scan.x scan.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o sniffer.x sniffer.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o spacefill.x spacefill.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o spectrum.x spectrum.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o superpose.x superpose.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o testgrad.x testgrad.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o testhess.x testhess.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o testpair.x testpair.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o testpol.x testpol.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o testrot.x testrot.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o timer.x timer.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o timerot.x timerot.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o torsfit.x torsfit.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o valence.x valence.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o vibbig.x vibbig.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o vibrate.x vibrate.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o vibrot.x vibrot.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o xtalfit.x xtalfit.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o xtalmin.x xtalmin.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o xyzedit.x xyzedit.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o xyzint.x xyzint.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o xyzmol2.x xyzmol2.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation
gfortran -Ofast -fopenmp -static-libgcc -o xyzpdb.x xyzpdb.o -L. libtinker.a libfftw3_threads.a libfftw3.a -framework JavaVM -framework CoreFoundation

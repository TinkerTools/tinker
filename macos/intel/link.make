#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the Tinker package programs  ##
#  ##             (Intel Fortran for macOS Version)             ##
#  ##                                                           ##
#  ###############################################################
#
#
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o alchemy.x alchemy.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip alchemy.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o analyze.x analyze.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip analyze.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o anneal.x anneal.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip anneal.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o arcedit.x arcedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip arcedit.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o bar.x bar.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip bar.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o correlate.x correlate.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip correlate.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o critical.x critical.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip critical.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o crystal.x crystal.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip crystal.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o diffuse.x diffuse.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip diffuse.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o distgeom.x distgeom.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip distgeom.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o document.x document.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip document.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o dynamic.x dynamic.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip dynamic.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o freefix.x freefix.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip freefix.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o gda.x gda.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip gda.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o intedit.x intedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip intedit.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o intxyz.x intxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip intxyz.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o minimize.x minimize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minimize.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o minirot.x minirot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minirot.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o minrigid.x minrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minrigid.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o mol2xyz.x mol2xyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip mol2xyz.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o molxyz.x molxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip molxyz.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o monte.x monte.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip monte.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o newton.x newton.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip newton.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o newtrot.x newtrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip newtrot.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o nucleic.x nucleic.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip nucleic.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o optimize.x optimize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optimize.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o optirot.x optirot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optirot.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o optrigid.x optrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optrigid.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o path.x path.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip path.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o pdbxyz.x pdbxyz.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pdbxyz.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o polarize.x polarize.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip polarize.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o poledit.x poledit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip poledit.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o potential.x potential.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip potential.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o prmedit.x prmedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip prmedit.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o protein.x protein.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip protein.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o pss.x pss.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pss.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o pssrigid.x pssrigid.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pssrigid.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o pssrot.x pssrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pssrot.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o radial.x radial.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip radial.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o saddle.x saddle.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip saddle.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o scan.x scan.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip scan.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o sniffer.x sniffer.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip sniffer.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o spacefill.x spacefill.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip spacefill.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o spectrum.x spectrum.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip spectrum.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o superpose.x superpose.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip superpose.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o testgrad.x testgrad.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testgrad.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o testhess.x testhess.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testhess.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o testpair.x testpair.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testpair.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o testpol.x testpol.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testpol.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o testrot.x testrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testrot.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o testsurf.x testsurf.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testsurf.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o testvir.x testvir.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testvir.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o timer.x timer.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip timer.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o timerot.x timerot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip timerot.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o torsfit.x torsfit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip torsfit.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o valence.x valence.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip valence.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o vibbig.x vibbig.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibbig.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o vibrate.x vibrate.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibrate.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o vibrot.x vibrot.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibrot.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o xtalfit.x xtalfit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xtalfit.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o xtalmin.x xtalmin.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xtalmin.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o xyzedit.x xyzedit.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzedit.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o xyzint.x xyzint.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzint.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o xyzmol2.x xyzmol2.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzmol2.x
ifort -O3 -no-ipo -no-prec-div -inline -mdynamic-no-pic -qopenmp -static-intel -Wl,-stack_size,0x10000000 -o xyzpdb.x xyzpdb.o -L. -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzpdb.x

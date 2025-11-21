#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the Tinker package programs  ##
#  ##               (LLVM flang for macOS Version)              ##
#  ##                                                           ##
#  ###############################################################
#
#
flang -O3 -fopenmp -o alchemy.x alchemy.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip alchemy.x
flang -O3 -fopenmp -o analyze.x analyze.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip analyze.x
flang -O3 -fopenmp -o anneal.x anneal.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip anneal.x
flang -O3 -fopenmp -o arcedit.x arcedit.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip arcedit.x
flang -O3 -fopenmp -o bar.x bar.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip bar.x
flang -O3 -fopenmp -o correlate.x correlate.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip correlate.x
flang -O3 -fopenmp -o critical.x critical.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip critical.x
flang -O3 -fopenmp -o crystal.x crystal.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip crystal.x
flang -O3 -fopenmp -o diffuse.x diffuse.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip diffuse.x
flang -O3 -fopenmp -o distgeom.x distgeom.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip distgeom.x
flang -O3 -fopenmp -o document.x document.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip document.x
flang -O3 -fopenmp -o dynamic.x dynamic.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip dynamic.x
flang -O3 -fopenmp -o freefix.x freefix.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip freefix.x
flang -O3 -fopenmp -o gda.x gda.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip gda.x
flang -O3 -fopenmp -o intedit.x intedit.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip intedit.x
flang -O3 -fopenmp -o intxyz.x intxyz.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip intxyz.x
flang -O3 -fopenmp -o mdavg.x mdavg.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip mdavg.x
flang -O3 -fopenmp -o minimize.x minimize.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minimize.x
flang -O3 -fopenmp -o minirot.x minirot.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minirot.x
flang -O3 -fopenmp -o minrigid.x minrigid.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip minrigid.x
flang -O3 -fopenmp -o mol2xyz.x mol2xyz.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip mol2xyz.x
flang -O3 -fopenmp -o molxyz.x molxyz.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip molxyz.x
flang -O3 -fopenmp -o monte.x monte.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip monte.x
flang -O3 -fopenmp -o newton.x newton.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip newton.x
flang -O3 -fopenmp -o newtrot.x newtrot.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip newtrot.x
flang -O3 -fopenmp -o nucleic.x nucleic.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip nucleic.x
flang -O3 -fopenmp -o optimize.x optimize.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optimize.x
flang -O3 -fopenmp -o optirot.x optirot.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optirot.x
flang -O3 -fopenmp -o optrigid.x optrigid.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip optrigid.x
flang -O3 -fopenmp -o path.x path.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip path.x
flang -O3 -fopenmp -o pdbxyz.x pdbxyz.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pdbxyz.x
flang -O3 -fopenmp -o polarize.x polarize.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip polarize.x
flang -O3 -fopenmp -o poledit.x poledit.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip poledit.x
flang -O3 -fopenmp -o potential.x potential.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip potential.x
flang -O3 -fopenmp -o prmedit.x prmedit.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip prmedit.x
flang -O3 -fopenmp -o protein.x protein.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip protein.x
flang -O3 -fopenmp -o pss.x pss.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pss.x
flang -O3 -fopenmp -o pssrigid.x pssrigid.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pssrigid.x
flang -O3 -fopenmp -o pssrot.x pssrot.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip pssrot.x
flang -O3 -fopenmp -o radial.x radial.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip radial.x
flang -O3 -fopenmp -o saddle.x saddle.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip saddle.x
flang -O3 -fopenmp -o scan.x scan.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip scan.x
flang -O3 -fopenmp -o sniffer.x sniffer.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip sniffer.x
flang -O3 -fopenmp -o spacefill.x spacefill.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip spacefill.x
flang -O3 -fopenmp -o spectrum.x spectrum.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip spectrum.x
flang -O3 -fopenmp -o superpose.x superpose.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip superpose.x
flang -O3 -fopenmp -o testgrad.x testgrad.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testgrad.x
flang -O3 -fopenmp -o testhess.x testhess.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testhess.x
flang -O3 -fopenmp -o testpair.x testpair.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testpair.x
flang -O3 -fopenmp -o testpol.x testpol.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testpol.x
flang -O3 -fopenmp -o testrot.x testrot.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testrot.x
flang -O3 -fopenmp -o testsurf.x testsurf.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testsurf.x
flang -O3 -fopenmp -o testvir.x testvir.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip testvir.x
flang -O3 -fopenmp -o timer.x timer.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip timer.x
flang -O3 -fopenmp -o timerot.x timerot.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip timerot.x
flang -O3 -fopenmp -o torsfit.x torsfit.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip torsfit.x
flang -O3 -fopenmp -o valence.x valence.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip valence.x
flang -O3 -fopenmp -o vibbig.x vibbig.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibbig.x
flang -O3 -fopenmp -o vibrate.x vibrate.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibrate.x
flang -O3 -fopenmp -o vibrot.x vibrot.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip vibrot.x
flang -O3 -fopenmp -o xtalfit.x xtalfit.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xtalfit.x
flang -O3 -fopenmp -o xtalmin.x xtalmin.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xtalmin.x
flang -O3 -fopenmp -o xyzedit.x xyzedit.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzedit.x
flang -O3 -fopenmp -o xyzint.x xyzint.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzint.x
flang -O3 -fopenmp -o xyzmol2.x xyzmol2.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzmol2.x
flang -O3 -fopenmp -o xyzpdb.x xyzpdb.o -L. -L/opt/homebrew/opt/llvm/lib -L../lib/macos -L../fftw/lib libtinker.a -lfftw3_threads -lfftw3 ; strip xyzpdb.x

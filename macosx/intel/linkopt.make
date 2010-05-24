#
#
#  ###################################################################
#  ##                                                               ##
#  ##  linkopt.make  --  optimized link of all the TINKER programs  ##
#  ##          (Intel Fortran Compiler for Mac OSX Version)         ##
#  ##                                                               ##
#  ###################################################################
#
#
ifort -fast -static-intel -vec-report0 -o alchemy.x alchemy.f *.o ; strip alchemy.x
ifort -fast -static-intel -vec-report0 -o analyze.x analyze.f *.o ; strip analyze.x
ifort -fast -static-intel -vec-report0 -o anneal.x anneal.f *.o ; strip anneal.x
ifort -fast -static-intel -vec-report0 -o archive.x archive.f *.o ; strip archive.x
ifort -fast -static-intel -vec-report0 -o correlate.x correlate.f *.o ; strip correlate.x
ifort -fast -static-intel -vec-report0 -o crystal.x crystal.f *.o ; strip crystal.x
ifort -fast -static-intel -vec-report0 -o diffuse.x diffuse.f *.o ; strip diffuse.x
ifort -fast -static-intel -vec-report0 -o distgeom.x distgeom.f *.o ; strip distgeom.x
ifort -fast -static-intel -vec-report0 -o document.x document.f *.o ; strip document.x
ifort -fast -static-intel -vec-report0 -o dynamic.x dynamic.f *.o ; strip dynamic.x
ifort -fast -static-intel -vec-report0 -o gda.x gda.f *.o ; strip gda.x
ifort -fast -static-intel -vec-report0 -o intedit.x intedit.f *.o ; strip intedit.x
ifort -fast -static-intel -vec-report0 -o intxyz.x intxyz.f *.o ; strip intxyz.x
ifort -fast -static-intel -vec-report0 -o minimize.x minimize.f *.o ; strip minimize.x
ifort -fast -static-intel -vec-report0 -o minirot.x minirot.f *.o ; strip minirot.x
ifort -fast -static-intel -vec-report0 -o minrigid.x minrigid.f *.o ; strip minrigid.x
ifort -fast -static-intel -vec-report0 -o monte.x monte.f *.o ; strip monte.x
ifort -fast -static-intel -vec-report0 -o newton.x newton.f *.o ; strip newton.x
ifort -fast -static-intel -vec-report0 -o newtrot.x newtrot.f *.o ; strip newtrot.x
ifort -fast -static-intel -vec-report0 -o nucleic.x nucleic.f *.o ; strip nucleic.x
ifort -fast -static-intel -vec-report0 -o optimize.x optimize.f *.o ; strip optimize.x
ifort -fast -static-intel -vec-report0 -o optirot.x optirot.f *.o ; strip optirot.x
ifort -fast -static-intel -vec-report0 -o optrigid.x optrigid.f *.o ; strip optrigid.x
ifort -fast -static-intel -vec-report0 -o path.x path.f *.o ; strip path.x
ifort -fast -static-intel -vec-report0 -o pdbxyz.x pdbxyz.f *.o ; strip pdbxyz.x
ifort -fast -static-intel -vec-report0 -o polarize.x polarize.f *.o ; strip polarize.x
ifort -fast -static-intel -vec-report0 -o poledit.x poledit.f *.o ; strip poledit.x
ifort -fast -static-intel -vec-report0 -o potential.x potential.f *.o ; strip potential.x
ifort -fast -static-intel -vec-report0 -o prmedit.x prmedit.f *.o ; strip prmedit.x
ifort -fast -static-intel -vec-report0 -o protein.x protein.f *.o ; strip protein.x
ifort -fast -static-intel -vec-report0 -o pss.x pss.f *.o ; strip pss.x
ifort -fast -static-intel -vec-report0 -o pssrigid.x pssrigid.f *.o ; strip pssrigid.x
ifort -fast -static-intel -vec-report0 -o pssrot.x pssrot.f *.o ; strip pssrot.x
ifort -fast -static-intel -vec-report0 -o radial.x radial.f *.o ; strip radial.x
ifort -fast -static-intel -vec-report0 -o saddle.x saddle.f *.o ; strip saddle.x
ifort -fast -static-intel -vec-report0 -o scan.x scan.f *.o ; strip scan.x
ifort -fast -static-intel -vec-report0 -o sniffer.x sniffer.f *.o ; strip sniffer.x
ifort -fast -static-intel -vec-report0 -o spacefill.x spacefill.f *.o ; strip spacefill.x
ifort -fast -static-intel -vec-report0 -o spectrum.x spectrum.f *.o ; strip spectrum.x
ifort -fast -static-intel -vec-report0 -o superpose.x superpose.f *.o ; strip superpose.x
ifort -fast -static-intel -vec-report0 -o sybylxyz.x sybylxyz.f *.o ; strip sybylxyz.x
ifort -fast -static-intel -vec-report0 -o testgrad.x testgrad.f *.o ; strip testgrad.x
ifort -fast -static-intel -vec-report0 -o testhess.x testhess.f *.o ; strip testhess.x
ifort -fast -static-intel -vec-report0 -o testpair.x testpair.f *.o ; strip testpair.x
ifort -fast -static-intel -vec-report0 -o testrot.x testrot.f *.o ; strip testrot.x
ifort -fast -static-intel -vec-report0 -o timer.x timer.f *.o ; strip timer.x
ifort -fast -static-intel -vec-report0 -o timerot.x timerot.f *.o ; strip timerot.x
ifort -fast -static-intel -vec-report0 -o valence.x valence.f *.o ; strip valence.x
ifort -fast -static-intel -vec-report0 -o vibbig.x vibbig.f *.o ; strip vibbig.x
ifort -fast -static-intel -vec-report0 -o vibrate.x vibrate.f *.o ; strip vibrate.x
ifort -fast -static-intel -vec-report0 -o vibrot.x vibrot.f *.o ; strip vibrot.x
ifort -fast -static-intel -vec-report0 -o xtalfit.x xtalfit.f *.o ; strip xtalfit.x
ifort -fast -static-intel -vec-report0 -o xtalmin.x xtalmin.f *.o ; strip xtalmin.x
ifort -fast -static-intel -vec-report0 -o xyzedit.x xyzedit.f *.o ; strip xyzedit.x
ifort -fast -static-intel -vec-report0 -o xyzint.x xyzint.f *.o ; strip xyzint.x
ifort -fast -static-intel -vec-report0 -o xyzpdb.x xyzpdb.f *.o ; strip xyzpdb.x
ifort -fast -static-intel -vec-report0 -o xyzsybyl.x xyzsybyl.f *.o ; strip xyzsybyl.x

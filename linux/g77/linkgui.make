#
#
#  ##################################################################
#  ##                                                              ##
#  ##  linkgui.make  --  link each of the TINKER programs for FFE  ##
#  ##                   (Linux/GNU g77 Version)                    ##
#  ##                                                              ##
#  ##################################################################
#
#
#  Copy the Java libjvm.so library from ../jre/lib/i386/client
#  in the standard Java distribution to use the script below,
#  or add it to the library search path
#
#
g77 -pthread -s -o alchemy.x alchemy.o libtinker.a libjvm.so 
g77 -pthread -s -o analyze.x analyze.o libtinker.a libjvm.so 
g77 -pthread -s -o anneal.x anneal.o libtinker.a libjvm.so 
g77 -pthread -s -o archive.x archive.o libtinker.a libjvm.so 
g77 -pthread -s -o correlate.x correlate.o libtinker.a libjvm.so 
g77 -pthread -s -o crystal.x crystal.o libtinker.a libjvm.so 
g77 -pthread -s -o diffuse.x diffuse.o libtinker.a libjvm.so 
g77 -pthread -s -o distgeom.x distgeom.o libtinker.a libjvm.so 
g77 -pthread -s -o document.x document.o libtinker.a libjvm.so 
g77 -pthread -s -o dynamic.x dynamic.o libtinker.a libjvm.so 
g77 -pthread -s -o gda.x gda.o libtinker.a libjvm.so 
g77 -pthread -s -o intedit.x intedit.o libtinker.a libjvm.so 
g77 -pthread -s -o intxyz.x intxyz.o libtinker.a libjvm.so 
g77 -pthread -s -o minimize.x minimize.o libtinker.a libjvm.so 
g77 -pthread -s -o minirot.x minirot.o libtinker.a libjvm.so 
g77 -pthread -s -o minrigid.x minrigid.o libtinker.a libjvm.so 
g77 -pthread -s -o monte.x monte.o libtinker.a libjvm.so 
g77 -pthread -s -o newton.x newton.o libtinker.a libjvm.so 
g77 -pthread -s -o newtrot.x newtrot.o libtinker.a libjvm.so 
g77 -pthread -s -o nucleic.x nucleic.o libtinker.a libjvm.so 
g77 -pthread -s -o optimize.x optimize.o libtinker.a libjvm.so 
g77 -pthread -s -o optirot.x optirot.o libtinker.a libjvm.so 
g77 -pthread -s -o optrigid.x optrigid.o libtinker.a libjvm.so 
g77 -pthread -s -o path.x path.o libtinker.a libjvm.so 
g77 -pthread -s -o pdbxyz.x pdbxyz.o libtinker.a libjvm.so 
g77 -pthread -s -o polarize.x polarize.o libtinker.a libjvm.so 
g77 -pthread -s -o poledit.x poledit.o libtinker.a libjvm.so 
g77 -pthread -s -o potential.x potential.o libtinker.a libjvm.so 
g77 -pthread -s -o prmedit.x prmedit.o libtinker.a libjvm.so 
g77 -pthread -s -o protein.x protein.o libtinker.a libjvm.so 
g77 -pthread -s -o pss.x pss.o libtinker.a libjvm.so 
g77 -pthread -s -o pssrigid.x pssrigid.o libtinker.a libjvm.so 
g77 -pthread -s -o pssrot.x pssrot.o libtinker.a libjvm.so 
g77 -pthread -s -o radial.x radial.o libtinker.a libjvm.so 
g77 -pthread -s -o saddle.x saddle.o libtinker.a libjvm.so 
g77 -pthread -s -o scan.x scan.o libtinker.a libjvm.so 
g77 -pthread -s -o sniffer.x sniffer.o libtinker.a libjvm.so 
g77 -pthread -s -o spacefill.x spacefill.o libtinker.a libjvm.so 
g77 -pthread -s -o spectrum.x spectrum.o libtinker.a libjvm.so 
g77 -pthread -s -o superpose.x superpose.o libtinker.a libjvm.so 
g77 -pthread -s -o sybylxyz.x sybylxyz.o libtinker.a libjvm.so 
g77 -pthread -s -o testgrad.x testgrad.o libtinker.a libjvm.so 
g77 -pthread -s -o testhess.x testhess.o libtinker.a libjvm.so 
g77 -pthread -s -o testpair.x testpair.o libtinker.a libjvm.so 
g77 -pthread -s -o testrot.x testrot.o libtinker.a libjvm.so 
g77 -pthread -s -o timer.x timer.o libtinker.a libjvm.so 
g77 -pthread -s -o timerot.x timerot.o libtinker.a libjvm.so 
g77 -pthread -s -o torsfit.x torsfit.o libtinker.a libjvm.so 
g77 -pthread -s -o valence.x valence.o libtinker.a libjvm.so 
g77 -pthread -s -o vibbig.x vibbig.o libtinker.a libjvm.so 
g77 -pthread -s -o vibrate.x vibrate.o libtinker.a libjvm.so 
g77 -pthread -s -o vibrot.x vibrot.o libtinker.a libjvm.so 
g77 -pthread -s -o xtalfit.x xtalfit.o libtinker.a libjvm.so 
g77 -pthread -s -o xtalmin.x xtalmin.o libtinker.a libjvm.so 
g77 -pthread -s -o xyzedit.x xyzedit.o libtinker.a libjvm.so 
g77 -pthread -s -o xyzint.x xyzint.o libtinker.a libjvm.so 
g77 -pthread -s -o xyzpdb.x xyzpdb.o libtinker.a libjvm.so 
g77 -pthread -s -o xyzsybyl.x xyzsybyl.o libtinker.a libjvm.so 

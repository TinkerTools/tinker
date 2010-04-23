#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##                  (Linux/GNU gcc Version)                  ##
#  ##                                                           ##
#  ###############################################################
#
#
gcc -o alchemy.x alchemy.o libtinker.a libf2c.a -lm
gcc -o analyze.x analyze.o libtinker.a libf2c.a -lm
gcc -o anneal.x anneal.o libtinker.a libf2c.a -lm
gcc -o archive.x archive.o libtinker.a libf2c.a -lm
gcc -o correlate.x correlate.o libtinker.a libf2c.a -lm
gcc -o crystal.x crystal.o libtinker.a libf2c.a -lm
gcc -o diffuse.x diffuse.o libtinker.a libf2c.a -lm
gcc -o distgeom.x distgeom.o libtinker.a libf2c.a -lm
gcc -o document.x document.o libtinker.a libf2c.a -lm
gcc -o dynamic.x dynamic.o libtinker.a libf2c.a -lm
gcc -o gda.x gda.o libtinker.a libf2c.a -lm
gcc -o intedit.x intedit.o libtinker.a libf2c.a -lm
gcc -o intxyz.x intxyz.o libtinker.a libf2c.a -lm
gcc -o minimize.x minimize.o libtinker.a libf2c.a -lm
gcc -o minirot.x minirot.o libtinker.a libf2c.a -lm
gcc -o minrigid.x minrigid.o libtinker.a libf2c.a -lm
gcc -o monte.x monte.o libtinker.a libf2c.a -lm
gcc -o newton.x newton.o libtinker.a libf2c.a -lm
gcc -o newtrot.x newtrot.o libtinker.a libf2c.a -lm
gcc -o nucleic.x nucleic.o libtinker.a libf2c.a -lm
gcc -o optimize.x optimize.o libtinker.a libf2c.a -lm
gcc -o optirot.x optirot.o libtinker.a libf2c.a -lm
gcc -o optrigid.x optrigid.o libtinker.a libf2c.a -lm
gcc -o path.x path.o libtinker.a libf2c.a -lm
gcc -o pdbxyz.x pdbxyz.o libtinker.a libf2c.a -lm
gcc -o polarize.x polarize.o libtinker.a libf2c.a -lm
gcc -o poledit.x poledit.o libtinker.a libf2c.a -lm
gcc -o potential.x potential.o libtinker.a libf2c.a -lm
gcc -o prmedit.x prmedit.o libtinker.a libf2c.a -lm
gcc -o protein.x protein.o libtinker.a libf2c.a -lm
gcc -o pss.x pss.o libtinker.a libf2c.a -lm
gcc -o pssrigid.x pssrigid.o libtinker.a libf2c.a -lm
gcc -o pssrot.x pssrot.o libtinker.a libf2c.a -lm
gcc -o radial.x radial.o libtinker.a libf2c.a -lm
gcc -o saddle.x saddle.o libtinker.a libf2c.a -lm
gcc -o scan.x scan.o libtinker.a libf2c.a -lm
gcc -o sniffer.x sniffer.o libtinker.a libf2c.a -lm
gcc -o spacefill.x spacefill.o libtinker.a libf2c.a -lm
gcc -o spectrum.x spectrum.o libtinker.a libf2c.a -lm
gcc -o superpose.x superpose.o libtinker.a libf2c.a -lm
gcc -o sybylxyz.x sybylxyz.o libtinker.a libf2c.a -lm
gcc -o testgrad.x testgrad.o libtinker.a libf2c.a -lm
gcc -o testhess.x testhess.o libtinker.a libf2c.a -lm
gcc -o testpair.x testpair.o libtinker.a libf2c.a -lm
gcc -o testrot.x testrot.o libtinker.a libf2c.a -lm
gcc -o timer.x timer.o libtinker.a libf2c.a -lm
gcc -o timerot.x timerot.o libtinker.a libf2c.a -lm
gcc -o valence.x valence.o libtinker.a libf2c.a -lm
gcc -o vibbig.x vibbig.o libtinker.a libf2c.a -lm
gcc -o vibrate.x vibrate.o libtinker.a libf2c.a -lm
gcc -o vibrot.x vibrot.o libtinker.a libf2c.a -lm
gcc -o xtalfit.x xtalfit.o libtinker.a libf2c.a -lm
gcc -o xtalmin.x xtalmin.o libtinker.a libf2c.a -lm
gcc -o xyzedit.x xyzedit.o libtinker.a libf2c.a -lm
gcc -o xyzint.x xyzint.o libtinker.a libf2c.a -lm
gcc -o xyzpdb.x xyzpdb.o libtinker.a libf2c.a -lm
gcc -o xyzsybyl.x xyzsybyl.o libtinker.a libf2c.a -lm

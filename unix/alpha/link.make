#
#
#  ###############################################################
#  ##                                                           ##
#  ##  link.make  --  link each of the TINKER package programs  ##
#  ##                   (HP Tru64 Unix Version)                 ##
#  ##                                                           ##
#  ###############################################################
#
#
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o alchemy.x alchemy.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o analyze.x analyze.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o anneal.x anneal.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o archive.x archive.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o correlate.x correlate.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o crystal.x crystal.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o diffuse.x diffuse.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o distgeom.x distgeom.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o document.x document.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o dynamic.x dynamic.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o gda.x gda.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o intedit.x intedit.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o intxyz.x intxyz.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o minimize.x minimize.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o minirot.x minirot.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o minrigid.x minrigid.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o monte.x monte.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o newton.x newton.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o newtrot.x newtrot.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o nucleic.x nucleic.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o optimize.x optimize.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o optirot.x optirot.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o optrigid.x optrigid.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o path.x path.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o pdbxyz.x pdbxyz.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o polarize.x polarize.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o poledit.x poledit.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o potential.x potential.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o prmedit.x prmedit.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o protein.x protein.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o pss.x pss.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o pssrigid.x pssrigid.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o pssrot.x pssrot.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o radial.x radial.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o saddle.x saddle.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o scan.x scan.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o sniffer.x sniffer.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o spacefill.x spacefill.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o spectrum.x spectrum.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o superpose.x superpose.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o sybylxyz.x sybylxyz.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o testgrad.x testgrad.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o testhess.x testhess.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o testpair.x testpair.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o testrot.x testrot.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o timer.x timer.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o timerot.x timerot.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o valence.x valence.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o vibbig.x vibbig.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o vibrate.x vibrate.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o vibrot.x vibrot.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o xtalfit.x xtalfit.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o xtalmin.x xtalmin.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o xyzedit.x xyzedit.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o xyzint.x xyzint.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o xyzpdb.x xyzpdb.o libtinker.a
f77 -fast -s -non_shared -om -WL,-om_no_inst_sched -o xyzsybyl.x xyzsybyl.o libtinker.a

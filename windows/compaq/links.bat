@echo off
rem
rem
rem  ###############################################################
rem  ##                                                           ##
rem  ##  links.bat  --  link each of the TINKER package programs  ##
rem  ##        (Compaq Visual Fortran for Windows Version)        ##
rem  ##                                                           ##
rem  ###############################################################
rem
rem
@echo on
DF /optimize /fast alchemy.obj tinker.lib
DF /optimize /fast analyze.obj tinker.lib
DF /optimize /fast anneal.obj tinker.lib
DF /optimize /fast archive.obj tinker.lib
DF /optimize /fast correlate.obj tinker.lib
DF /optimize /fast crystal.obj tinker.lib
DF /optimize /fast diffuse.obj tinker.lib
DF /optimize /fast distgeom.obj tinker.lib
DF /optimize /fast document.obj tinker.lib
DF /optimize /fast dynamic.obj tinker.lib
DF /optimize /fast gda.obj tinker.lib
DF /optimize /fast intedit.obj tinker.lib
DF /optimize /fast intxyz.obj tinker.lib
DF /optimize /fast minimize.obj tinker.lib
DF /optimize /fast minirot.obj tinker.lib
DF /optimize /fast minrigid.obj tinker.lib
DF /optimize /fast monte.obj tinker.lib
DF /optimize /fast newton.obj tinker.lib
DF /optimize /fast newtrot.obj tinker.lib
DF /optimize /fast nucleic.obj tinker.lib
DF /optimize /fast optimize.obj tinker.lib
DF /optimize /fast optirot.obj tinker.lib
DF /optimize /fast optrigid.obj tinker.lib
DF /optimize /fast path.obj tinker.lib
DF /optimize /fast pdbxyz.obj tinker.lib
DF /optimize /fast polarize.obj tinker.lib
DF /optimize /fast poledit.obj tinker.lib
DF /optimize /fast potential.obj tinker.lib
DF /optimize /fast prmedit.obj tinker.lib
DF /optimize /fast protein.obj tinker.lib
DF /optimize /fast pss.obj tinker.lib
DF /optimize /fast pssrigid.obj tinker.lib
DF /optimize /fast pssrot.obj tinker.lib
DF /optimize /fast radial.obj tinker.lib
DF /optimize /fast saddle.obj tinker.lib
DF /optimize /fast scan.obj tinker.lib
DF /optimize /fast sniffer.obj tinker.lib
DF /optimize /fast spacefill.obj tinker.lib
DF /optimize /fast spectrum.obj tinker.lib
DF /optimize /fast superpose.obj tinker.lib
DF /optimize /fast sybylxyz.obj tinker.lib
DF /optimize /fast testgrad.obj tinker.lib
DF /optimize /fast testhess.obj tinker.lib
DF /optimize /fast testpair.obj tinker.lib
DF /optimize /fast testrot.obj tinker.lib
DF /optimize /fast timer.obj tinker.lib
DF /optimize /fast timerot.obj tinker.lib
DF /optimize /fast valence.obj tinker.lib
DF /optimize /fast vibbig.obj tinker.lib
DF /optimize /fast vibrate.obj tinker.lib
DF /optimize /fast vibrot.obj tinker.lib
DF /optimize /fast xtalfit.obj tinker.lib
DF /optimize /fast xtalmin.obj tinker.lib
DF /optimize /fast xyzedit.obj tinker.lib
DF /optimize /fast xyzint.obj tinker.lib
DF /optimize /fast xyzpdb.obj tinker.lib
DF /optimize /fast xyzsybyl.obj tinker.lib

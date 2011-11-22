@echo off
rem
rem
rem  ###############################################################
rem  ##                                                           ##
rem  ##  links.bat  --  link each of the TINKER package programs  ##
rem  ##        (Intel Fortran Compiler for Windows Version)       ##
rem  ##                                                           ##
rem  ###############################################################
rem
rem
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib alchemy.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib analyze.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib anneal.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib archive.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib correlate.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib crystal.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib diffuse.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib distgeom.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib document.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib dynamic.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib gda.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib intedit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib intxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib minimize.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib minirot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib minrigid.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib monte.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib newton.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib newtrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib nucleic.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib optimize.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib optirot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib optrigid.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib path.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib pdbxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib polarize.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib poledit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib potential.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib prmedit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib protein.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib pss.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib pssrigid.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib pssrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib radial.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib saddle.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib scan.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib sniffer.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib spacefill.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib spectrum.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib superpose.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib sybylxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib testgrad.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib testhess.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib testpair.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib testrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib timer.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib timerot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib torsfit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib valence.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib vibbig.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib vibrate.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib vibrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib xtalfit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib xtalmin.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib xyzedit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib xyzint.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib xyzpdb.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /4Yportlib xyzsybyl.obj tinker.lib

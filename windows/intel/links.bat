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
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static alchemy.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static analyze.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static anneal.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static archive.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static correlate.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static crystal.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static diffuse.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static distgeom.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static document.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static dynamic.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static gda.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static intedit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static intxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static minimize.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static minirot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static minrigid.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static molxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static monte.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static newton.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static newtrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static nucleic.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static optimize.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static optirot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static optrigid.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static path.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static pdbxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static polarize.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static poledit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static potential.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static prmedit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static protein.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static pss.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static pssrigid.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static pssrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static radial.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static saddle.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static scan.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static sniffer.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static spacefill.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static spectrum.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static superpose.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static sybylxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static testgrad.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static testhess.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static testpair.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static testrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static timer.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static timerot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static torsfit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static valence.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static vibbig.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static vibrate.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static vibrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static xtalfit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static xtalmin.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static xyzedit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static xyzint.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static xyzpdb.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /libs:static /Qopenmp-link:static xyzsybyl.obj tinker.lib

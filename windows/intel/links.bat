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
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static alchemy.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static analyze.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static anneal.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static archive.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static bar.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static correlate.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static crystal.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static diffuse.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static distgeom.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static document.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static dynamic.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static gda.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static intedit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static intxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static minimize.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static minirot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static minrigid.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static molxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static monte.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static newton.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static newtrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static nucleic.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static optimize.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static optirot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static optrigid.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static path.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static pdbxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static polarize.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static poledit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static potential.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static prmedit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static protein.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static pss.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static pssrigid.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static pssrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static radial.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static saddle.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static scan.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static sniffer.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static spacefill.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static spectrum.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static superpose.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static sybylxyz.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static testgrad.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static testhess.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static testpair.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static testpol.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static testrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static timer.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static timerot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static torsfit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static valence.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static vibbig.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static vibrate.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static vibrot.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static xtalfit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static xtalmin.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static xyzedit.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static xyzint.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static xyzpdb.obj tinker.lib
ifort /O3 /Qprec-div- /Qopenmp /recursive /libs:static xyzsybyl.obj tinker.lib

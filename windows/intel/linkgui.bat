@echo off
rem
rem
rem  #################################################################
rem  ##                                                             ##
rem  ##  linkgui.bat  --  link each of the TINKER programs for FFE  ##
rem  ##         (Intel Fortran Compiler for Windows Version)        ##
rem  ##                                                             ##
rem  #################################################################
rem
rem
rem  Create a variable with the static libraries to link against;
rem  Change "javalib" to reflect location of the Java library
rem
set javalib="C:\Program Files\Java\jdk1.6.0_22\lib\jvm.lib"
rem
rem  Link using Intel compiler suite and compatibility library
rem
ifort /libs:static alchemy.obj tinker.lib %javalib%
ifort /libs:static analyze.obj tinker.lib %javalib%
ifort /libs:static anneal.obj tinker.lib %javalib%
ifort /libs:static archive.obj tinker.lib %javalib%
ifort /libs:static bar.obj tinker.lib %javalib%
ifort /libs:static correlate.obj tinker.lib %javalib%
ifort /libs:static crystal.obj tinker.lib %javalib%
ifort /libs:static diffuse.obj tinker.lib %javalib%
ifort /libs:static distgeom.obj tinker.lib %javalib%
ifort /libs:static document.obj tinker.lib %javalib%
ifort /libs:static dynamic.obj tinker.lib %javalib%
ifort /libs:static gda.obj tinker.lib %javalib%
ifort /libs:static intedit.obj tinker.lib %javalib%
ifort /libs:static intxyz.obj tinker.lib %javalib%
ifort /libs:static minimize.obj tinker.lib %javalib%
ifort /libs:static minirot.obj tinker.lib %javalib%
ifort /libs:static minrigid.obj tinker.lib %javalib%
ifort /libs:static molxyz.obj tinker.lib %javalib%
ifort /libs:static monte.obj tinker.lib %javalib%
ifort /libs:static newton.obj tinker.lib %javalib%
ifort /libs:static newtrot.obj tinker.lib %javalib%
ifort /libs:static nucleic.obj tinker.lib %javalib%
ifort /libs:static optimize.obj tinker.lib %javalib%
ifort /libs:static optirot.obj tinker.lib %javalib%
ifort /libs:static optrigid.obj tinker.lib %javalib%
ifort /libs:static path.obj tinker.lib %javalib%
ifort /libs:static pdbxyz.obj tinker.lib %javalib%
ifort /libs:static polarize.obj tinker.lib %javalib%
ifort /libs:static prmedit.obj tinker.lib %javalib%
ifort /libs:static protein.obj tinker.lib %javalib%
ifort /libs:static pss.obj tinker.lib %javalib%
ifort /libs:static pssrigid.obj tinker.lib %javalib%
ifort /libs:static pssrot.obj tinker.lib %javalib%
ifort /libs:static radial.obj tinker.lib %javalib%
ifort /libs:static saddle.obj tinker.lib %javalib%
ifort /libs:static scan.obj tinker.lib %javalib%
ifort /libs:static sniffer.obj tinker.lib %javalib%
ifort /libs:static spacefill.obj tinker.lib %javalib%
ifort /libs:static spectrum.obj tinker.lib %javalib%
ifort /libs:static superpose.obj tinker.lib %javalib%
ifort /libs:static sybylxyz.obj tinker.lib %javalib%
ifort /libs:static testgrad.obj tinker.lib %javalib%
ifort /libs:static testhess.obj tinker.lib %javalib%
ifort /libs:static testpair.obj tinker.lib %javalib%
ifort /libs:static testpol.obj tinker.lib %javalib%
ifort /libs:static testrot.obj tinker.lib %javalib%
ifort /libs:static timer.obj tinker.lib %javalib%
ifort /libs:static timerot.obj tinker.lib %javalib%
ifort /libs:static torsfit.obj tinker.lib %javalib%
ifort /libs:static valence.obj tinker.lib %javalib%
ifort /libs:static vibbig.obj tinker.lib %javalib%
ifort /libs:static vibrate.obj tinker.lib %javalib%
ifort /libs:static vibrot.obj tinker.lib %javalib%
ifort /libs:static xtalfit.obj tinker.lib %javalib%
ifort /libs:static xtalmin.obj tinker.lib %javalib%
ifort /libs:static xyzedit.obj tinker.lib %javalib%
ifort /libs:static xyzint.obj tinker.lib %javalib%
ifort /libs:static xyzpdb.obj tinker.lib %javalib%
ifort /libs:static xyzsybyl.obj tinker.lib %javalib%

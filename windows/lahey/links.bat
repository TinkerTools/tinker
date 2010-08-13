@echo off
rem
rem
rem  ###############################################################
rem  ##                                                           ##
rem  ##  links.bat  --  link each of the TINKER package programs  ##
rem  ##          (Lahey/Fujitsu LF95 for Windows Version)         ##
rem  ##                                                           ##
rem  ###############################################################
rem
rem
@echo on
LF95 -o1 -tpp -nstchk -ntrace -nomap alchemy.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap analyze.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap anneal.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap archive.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap correlate.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap crystal.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap diffuse.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap distgeom.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap document.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap dynamic.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap gda.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap intedit.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap intxyz.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap minimize.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap minirot.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap minrigid.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap monte.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap newton.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap newtrot.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap nucleic.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap optimize.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap optirot.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap optrigid.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap path.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap pdbxyz.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap polarize.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap poledit.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap potential.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap prmedit.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap protein.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap pss.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap pssrigid.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap pssrot.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap radial.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap saddle.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap scan.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap sniffer.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap spacefill.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap spectrum.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap superpose.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap sybylxyz.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap testgrad.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap testhess.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap testpair.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap testrot.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap timer.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap timerot.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap torsfit.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap valence.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap vibbig.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap vibrate.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap vibrot.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap xtalfit.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap xtalmin.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap xyzedit.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap xyzint.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap xyzpdb.obj -lib tinker
LF95 -o1 -tpp -nstchk -ntrace -nomap xyzsybyl.obj -lib tinker

##
###################################################################
##                                                               ##
##  Makefile for Building the TINKER Molecular Modeling Package  ##
##                                                               ##
###################################################################
##
##  Invocation Options:
##
##   1. make all              Build all the TINKER executables
##   2. make install          Move the executables to BINDIR
##   3. make create_links     Create soft links in LINKDIR
##   4. make remove_links     Remove soft links from LINKDIR
##   6. make listing          Concatenate source to tinker.txt
##   5. make clean            Delete objects and executables
##
##  Original version of this file is due to Peter Happersberger
##  and Jochen Buehler of the University of Konstanz, January 1998.
##  Modifications by Reece Hart & Jay Ponder, Washington University
##
##  Building TINKER from source requires libraries from the FFTW
##  Fast Fourier Transform package and (optionally) APBS Poisson-
##  Boltzmann solver. These libraries are assumed to already be
##  present, and are not built by this Makefile
##
###################################################################

###################################################################
##  Master Directory Locations; Change as Needed for Local Site  ##
###################################################################

##  TINKERDIR    TINKER Distribution Directory
##  BINDIR       Hard Copies of TINKER Executables
##  LINKDIR      Linked Copies of TINKER Executables

TINKERDIR = $(HOME)/tinker
BINDIR = $(TINKERDIR)/bin
LINKDIR = /usr/local/bin

####################################################################
##  Known Machine Types; Uncomment One of the Following Sections  ##
##  May Need Editing to Match Your Desired OS & Compiler Version  ##
####################################################################

##  Machine:  Generic Linux or Macintosh
##  CPU Type: Intel x86 Compatible
##  Compiler: GNU gfortran

#F77 = gfortran
#LIBS =
#F77FLAGS = -c
#OPTFLAGS = -O3 -ffast-math
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-libgcc

##  Machine:  Generic Linux or Macintosh
##  CPU Type: Intel x86 Compatible
##  Compiler: GNU gfortran
##  Parallel: OpenMP

F77 = gfortran
LIBS = -L$(TINKERDIR)/fftw/lib -lfftw3_threads -lfftw3
F77FLAGS = -c
OPTFLAGS = -O3 -ffast-math -fopenmp
LIBFLAGS = -crusv
RANLIB = ranlib -c
LINKFLAGS = $(OPTFLAGS) -static-libgcc

##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible (also AMD)
##  Oper Sys: Fedora Core
##  Compiler: Intel Fortran for Linux 14.0

#F77 = ifort
#LIBS =
#F77FLAGS = -c -xHost
#OPTFLAGS = -O3 -no-ipo -no-prec-div -recursive
#LIBFLAGS = -crusv
#RANLIB = echo
#LINKFLAGS = -static-libgcc -static-intel

##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible (also AMD)
##  Oper Sys: Fedora Core
##  Compiler: Intel Fortran for Linux 14.0
##  Parallel: OpenMP

#F77 = ifort
#LIBS = -L$(TINKERDIR)/fftw/lib -lfftw3_threads -lfftw3
#F77FLAGS = -c -xHost
#OPTFLAGS = -O3 -no-ipo -no-prec-div -recursive -openmp
#LIBFLAGS = -crusv
#RANLIB = echo
#LINKFLAGS = $(OPTFLAGS) -static-libgcc -static-intel

##  Machine:  Macintosh
##  CPU Type: Intel Xeon
##  Oper Sys: OS X 10.9 (Mavericks)
##  Compiler: Intel Fortran for Mac 14.0

#F77 = ifort
#LIBS =
#F77FLAGS = -c -axSSSE3
#OPTFLAGS = -O3 -no-ipo -no-prec-div -mdynamic-no-pic
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-intel -mmacosx-version-min=10.6 \
             -Wl,-stack_size,0x10000000

##  Machine:  Macintosh
##  CPU Type: Intel Xeon
##  Oper Sys: OS X 10.9 (Mavericks)
##  Compiler: Intel Fortran for Mac 14.0
##  Parallel: OpenMP

#F77 = ifort
#LIBS = -L$(TINKERDIR)/fftw/lib -lfftw3_threads -lfftw3
#F77FLAGS = -c -axSSSE3
#OPTFLAGS = -O3 -no-ipo -no-prec-div -mdynamic-no-pic -openmp
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-intel -mmacosx-version-min=10.6 \
#            -Wl,-stack_size,0x10000000

#################################################################
##  Should not be Necessary to Change Things Below this Point  ##
#################################################################

OBJS = action.o \
       active.o \
       alchemy.o \
       align.o \
       analysis.o \
       analyz.o \
       analyze.o \
       angang.o \
       angbnd.o \
       angles.o \
       angpot.o \
       angtor.o \
       anneal.o \
       archive.o \
       argue.o \
       ascii.o \
       atmlst.o \
       atomid.o \
       atoms.o \
       attach.o \
       bar.o \
       basefile.o \
       bath.o \
       beeman.o \
       bicubic.o \
       bitor.o \
       bitors.o \
       bndpot.o \
       bndstr.o \
       bonds.o \
       born.o \
       bound.o \
       bounds.o \
       boxes.o \
       bussi.o \
       calendar.o \
       cell.o \
       center.o \
       charge.o \
       chgpot.o \
       chkpole.o \
       chkring.o \
       chkxyz.o \
       cholesky.o \
       chrono.o \
       chunks.o \
       clock.o \
       cluster.o \
       column.o \
       command.o \
       connect.o \
       connolly.o \
       control.o \
       correlate.o \
       couple.o \
       crystal.o \
       cspline.o \
       cutoffs.o \
       deflate.o \
       delete.o \
       deriv.o \
       diagq.o \
       diffeq.o \
       diffuse.o \
       dipole.o \
       disgeo.o \
       distgeom.o \
       dma.o \
       document.o \
       domega.o \
       dynamic.o \
       eangang.o \
       eangang1.o \
       eangang2.o \
       eangang3.o \
       eangle.o \
       eangle1.o \
       eangle2.o \
       eangle3.o \
       eangtor.o \
       eangtor1.o \
       eangtor2.o \
       eangtor3.o \
       ebond.o \
       ebond1.o \
       ebond2.o \
       ebond3.o \
       ebuck.o \
       ebuck1.o \
       ebuck2.o \
       ebuck3.o \
       echarge.o \
       echarge1.o \
       echarge2.o \
       echarge3.o \
       echgdpl.o \
       echgdpl1.o \
       echgdpl2.o \
       echgdpl3.o \
       edipole.o \
       edipole1.o \
       edipole2.o \
       edipole3.o \
       egauss.o \
       egauss1.o \
       egauss2.o \
       egauss3.o \
       egeom.o \
       egeom1.o \
       egeom2.o \
       egeom3.o \
       ehal.o \
       ehal1.o \
       ehal2.o \
       ehal3.o \
       eimprop.o \
       eimprop1.o \
       eimprop2.o \
       eimprop3.o \
       eimptor.o \
       eimptor1.o \
       eimptor2.o \
       eimptor3.o \
       elj.o \
       elj1.o \
       elj2.o \
       elj3.o \
       embed.o \
       emetal.o \
       emetal1.o \
       emetal2.o \
       emetal3.o \
       emm3hb.o \
       emm3hb1.o \
       emm3hb2.o \
       emm3hb3.o \
       empole.o \
       empole1.o \
       empole2.o \
       empole3.o \
       energi.o \
       energy.o \
       eopbend.o \
       eopbend1.o \
       eopbend2.o \
       eopbend3.o \
       eopdist.o \
       eopdist1.o \
       eopdist2.o \
       eopdist3.o \
       epitors.o \
       epitors1.o \
       epitors2.o \
       epitors3.o \
       erf.o \
       erxnfld.o \
       erxnfld1.o \
       erxnfld2.o \
       erxnfld3.o \
       esolv.o \
       esolv1.o \
       esolv2.o \
       esolv3.o \
       estrbnd.o \
       estrbnd1.o \
       estrbnd2.o \
       estrbnd3.o \
       estrtor.o \
       estrtor1.o \
       estrtor2.o \
       estrtor3.o \
       etors.o \
       etors1.o \
       etors2.o \
       etors3.o \
       etortor.o \
       etortor1.o \
       etortor2.o \
       etortor3.o \
       eurey.o \
       eurey1.o \
       eurey2.o \
       eurey3.o \
       evcorr.o \
       ewald.o \
       extra.o \
       extra1.o \
       extra2.o \
       extra3.o \
       faces.o \
       fatal.o \
       fft.o \
       fft3d.o \
       fftpack.o \
       field.o \
       fields.o \
       files.o \
       final.o \
       flatten.o \
       fracs.o \
       freeunit.o \
       freeze.o \
       gda.o \
       geometry.o \
       getint.o \
       getkey.o \
       getmol.o \
       getmol2.o \
       getnumb.o \
       getpdb.o \
       getprm.o \
       getref.o \
       getstring.o \
       gettext.o \
       getword.o \
       getxyz.o \
       ghmcstep.o \
       gkstuf.o \
       gradient.o \
       gradrgd.o \
       gradrot.o \
       group.o \
       groups.o \
       grpline.o \
       gyrate.o \
       hescut.o \
       hessian.o \
       hessn.o \
       hessrgd.o \
       hessrot.o \
       hpmf.o \
       hybrid.o \
       image.o \
       impose.o \
       improp.o \
       imptor.o \
       induce.o \
       inertia.o \
       inform.o \
       initatom.o \
       initial.o \
       initprm.o \
       initres.o \
       initrot.o \
       insert.o \
       intedit.o \
       inter.o \
       intxyz.o \
       invbeta.o \
       invert.o \
       iounit.o \
       jacobi.o \
       kanang.o \
       kangang.o \
       kangle.o \
       kangs.o \
       kangtor.o \
       kantor.o \
       katom.o \
       katoms.o \
       kbond.o \
       kbonds.o \
       kcharge.o \
       kchrge.o \
       kdipol.o \
       kdipole.o \
       kewald.o \
       kextra.o \
       keys.o \
       kgeom.o \
       khbond.o \
       kimprop.o \
       kimptor.o \
       kinetic.o \
       kiprop.o \
       kitors.o \
       kmetal.o \
       kmpole.o \
       kmulti.o \
       kopbend.o \
       kopbnd.o \
       kopdist.o \
       kopdst.o \
       korbit.o \
       korbs.o \
       kpitor.o \
       kpitors.o \
       kpolar.o \
       kpolr.o \
       ksolv.o \
       kstbnd.o \
       kstrbnd.o \
       kstrtor.o \
       ksttor.o \
       ktors.o \
       ktorsn.o \
       ktortor.o \
       ktrtor.o \
       kurey.o \
       kurybr.o \
       kvdw.o \
       kvdwpr.o \
       kvdws.o \
       lattice.o \
       lbfgs.o \
       light.o \
       lights.o \
       limits.o \
       linmin.o \
       makeint.o \
       makeref.o \
       makexyz.o \
       math.o \
       maxwell.o \
       mdinit.o \
       mdrest.o \
       mdsave.o \
       mdstat.o \
       mdstuf.o \
       mechanic.o \
       merck.o \
       merge.o \
       minima.o \
       minimize.o \
       minirot.o \
       minrigid.o \
       molcul.o \
       moldyn.o \
       molecule.o \
       molxyz.o \
       moment.o \
       moments.o \
       monte.o \
       mplpot.o \
       mpole.o \
       mutant.o \
       mutate.o \
       nblist.o \
       neigh.o \
       newton.o \
       newtrot.o \
       nextarg.o \
       nexttext.o \
       nonpol.o \
       nose.o \
       nspline.o \
       nucleic.o \
       nucleo.o \
       number.o \
       numeral.o \
       numgrad.o \
       ocvm.o \
       omega.o \
       opbend.o \
       opdist.o \
       openend.o \
       openmp.o \
       optimize.o \
       optirot.o \
       optrigid.o \
       optsave.o \
       orbital.o \
       orbits.o \
       orient.o \
       orthog.o \
       output.o \
       overlap.o \
       params.o \
       path.o \
       paths.o \
       pbstuf.o \
       pdb.o \
       pdbxyz.o \
       phipsi.o \
       picalc.o \
       piorbs.o \
       pistuf.o \
       pitors.o \
       pme.o \
       pmestuff.o \
       pmpb.o \
       polar.o \
       polarize.o \
       poledit.o \
       polgrp.o \
       polpot.o \
       polymer.o \
       potent.o \
       potential.o \
       potfit.o \
       precis.o \
       precise.o \
       pressure.o \
       prmedit.o \
       prmkey.o \
       promo.o \
       protein.o \
       prtdyn.o \
       prterr.o \
       prtint.o \
       prtmol2.o \
       prtpdb.o \
       prtprm.o \
       prtseq.o \
       prtxyz.o \
       pss.o \
       pssrigid.o \
       pssrot.o \
       ptable.o \
       qmstuf.o \
       qrfact.o \
       quatfit.o \
       radial.o \
       random.o \
       rattle.o \
       readdyn.o \
       readgau.o \
       readint.o \
       readmol.o \
       readmol2.o \
       readpdb.o \
       readprm.o \
       readseq.o \
       readxyz.o \
       refer.o \
       replica.o \
       resdue.o \
       respa.o \
       restrn.o \
       rgddyn.o \
       rgdstep.o \
       rigid.o \
       ring.o \
       rings.o \
       rmsfit.o \
       rotbnd.o \
       rotlist.o \
       rotpole.o \
       rxnfld.o \
       rxnpot.o \
       saddle.o \
       scales.o \
       scan.o \
       sdstep.o \
       search.o \
       sequen.o \
       server.o \
       shakeup.o \
       shunt.o \
       sigmoid.o \
       sizes.o \
       sktstuff.o \
       sniffer.o \
       socket.o \
       solute.o \
       sort.o \
       spacefill.o \
       spectrum.o \
       square.o \
       stodyn.o \
       strbnd.o \
       strtor.o \
       suffix.o \
       superpose.o \
       surface.o \
       surfatom.o \
       switch.o \
       sybylxyz.o \
       syntrn.o \
       tarray.o \
       temper.o \
       testgrad.o \
       testhess.o \
       testpair.o \
       testpol.o \
       testrot.o \
       timer.o \
       timerot.o \
       titles.o \
       tncg.o \
       torphase.o \
       torpot.o \
       torque.o \
       tors.o \
       torsfit.o \
       torsions.o \
       tortor.o \
       tree.o \
       trimtext.o \
       unitcell.o \
       units.o \
       uprior.o \
       urey.o \
       urypot.o \
       usage.o \
       usolve.o \
       valence.o \
       valfit.o \
       vdw.o \
       vdwpot.o \
       verlet.o \
       version.o \
       vibbig.o \
       vibrate.o \
       vibrot.o \
       vibs.o \
       virial.o \
       volume.o \
       warp.o \
       xtalfit.o \
       xtalmin.o \
       xtals.o \
       xyzatm.o \
       xyzedit.o \
       xyzint.o \
       xyzpdb.o \
       xyzsybyl.o \
       zatom.o \
       zclose.o \
       zcoord.o

EXEFILES = alchemy.x \
           analyze.x \
           anneal.x \
           archive.x \
           bar.x \
           correlate.x \
           crystal.x \
           diffuse.x \
           distgeom.x \
           document.x \
           dynamic.x \
           gda.x \
           intedit.x \
           intxyz.x \
           minimize.x \
           minirot.x \
           minrigid.x \
           molxyz.x \
           monte.x \
           newton.x \
           newtrot.x \
           nucleic.x \
           optimize.x \
           optirot.x \
           optrigid.x \
           path.x \
           pdbxyz.x \
           polarize.x \
           poledit.x \
           potential.x \
           prmedit.x \
           protein.x \
           pss.x \
           pssrigid.x \
           pssrot.x \
           radial.x \
           saddle.x \
           scan.x \
           sniffer.x \
           spacefill.x \
           spectrum.x \
           superpose.x \
           sybylxyz.x \
           testgrad.x \
           testhess.x \
           testpair.x \
           testpol.x \
           testrot.x \
           timer.x \
           timerot.x \
           torsfit.x \
           valence.x \
           vibbig.x \
           vibrate.x \
           vibrot.x \
           xtalfit.x \
           xtalmin.x \
           xyzedit.x \
           xyzint.x \
           xyzpdb.x \
           xyzsybyl.x

%.o: %.f
	${F77} ${F77FLAGS} ${OPTFLAGS} $< -o $@ 

%.o: %.c
	${CC} ${CFLAGS} ${INCLUDEDIR} ${OPTFLAGS} $<

%.x: %.o libtinker.a
	${F77} ${LINKFLAGS} -o $@ $^ ${LIBS}; strip $@

all:	${EXEFILES}

install:
	rename

clean:
	rm -f *.o *.mod *.a *.x

listing:
	cat *.i *.f *.c > tinker.txt

rename:
	mv  alchemy.x    $(BINDIR)/alchemy
	mv  analyze.x    $(BINDIR)/analyze
	mv  anneal.x     $(BINDIR)/anneal
	mv  archive.x    $(BINDIR)/archive
	mv  bar.x        $(BINDIR)/bar
	mv  correlate.x  $(BINDIR)/correlate
	mv  crystal.x    $(BINDIR)/crystal
	mv  diffuse.x    $(BINDIR)/diffuse
	mv  distgeom.x   $(BINDIR)/distgeom
	mv  document.x   $(BINDIR)/document
	mv  dynamic.x    $(BINDIR)/dynamic
	mv  gda.x        $(BINDIR)/gda
	mv  intedit.x    $(BINDIR)/intedit
	mv  intxyz.x     $(BINDIR)/intxyz
	mv  minimize.x   $(BINDIR)/minimize
	mv  minirot.x    $(BINDIR)/minirot
	mv  minrigid.x   $(BINDIR)/minrigid
	mv  molxyz.x     $(BINDIR)/molxyz
	mv  monte.x      $(BINDIR)/monte
	mv  newton.x     $(BINDIR)/newton
	mv  newtrot.x    $(BINDIR)/newtrot
	mv  nucleic.x    $(BINDIR)/nucleic
	mv  optimize.x   $(BINDIR)/optimize
	mv  optirot.x    $(BINDIR)/optirot
	mv  optrigid.x   $(BINDIR)/optrigid
	mv  path.x       $(BINDIR)/path
	mv  pdbxyz.x     $(BINDIR)/pdbxyz
	mv  polarize.x   $(BINDIR)/polarize
	mv  poledit.x    $(BINDIR)/poledit
	mv  potential.x  $(BINDIR)/potential
	mv  prmedit.x    $(BINDIR)/prmedit
	mv  protein.x    $(BINDIR)/protein
	mv  pss.x        $(BINDIR)/pss
	mv  pssrigid.x   $(BINDIR)/pssrigid
	mv  pssrot.x     $(BINDIR)/pssrot
	mv  radial.x     $(BINDIR)/radial
	mv  saddle.x     $(BINDIR)/saddle
	mv  scan.x       $(BINDIR)/scan
	mv  sniffer.x    $(BINDIR)/sniffer
	mv  spacefill.x  $(BINDIR)/spacefill
	mv  spectrum.x   $(BINDIR)/spectrum
	mv  superpose.x  $(BINDIR)/superpose
	mv  sybylxyz.x   $(BINDIR)/sybylxyz
	mv  testgrad.x   $(BINDIR)/testgrad
	mv  testhess.x   $(BINDIR)/testhess
	mv  testpair.x   $(BINDIR)/testpair
	mv  testpol.x    $(BINDIR)/testpol
	mv  testrot.x    $(BINDIR)/testrot
	mv  timer.x      $(BINDIR)/timer
	mv  timerot.x    $(BINDIR)/timerot
	mv  torsfit.x    $(BINDIR)/torsfit
	mv  valence.x    $(BINDIR)/valence
	mv  vibbig.x     $(BINDIR)/vibbig
	mv  vibrate.x    $(BINDIR)/vibrate
	mv  vibrot.x     $(BINDIR)/vibrot
	mv  xtalfit.x    $(BINDIR)/xtalfit
	mv  xtalmin.x    $(BINDIR)/xtalmin
	mv  xyzedit.x    $(BINDIR)/xyzedit
	mv  xyzint.x     $(BINDIR)/xyzint
	mv  xyzpdb.x     $(BINDIR)/xyzpdb
	mv  xyzsybyl.x   $(BINDIR)/xyzsybyl

remove_links:
	rm -f $(LINKDIR)/alchemy
	rm -f $(LINKDIR)/analyze
	rm -f $(LINKDIR)/anneal
	rm -f $(LINKDIR)/archive
	rm -f $(LINKDIR)/bar
	rm -f $(LINKDIR)/correlate
	rm -f $(LINKDIR)/crystal
	rm -f $(LINKDIR)/diffuse
	rm -f $(LINKDIR)/distgeom
	rm -f $(LINKDIR)/document
	rm -f $(LINKDIR)/dynamic
	rm -f $(LINKDIR)/gda
	rm -f $(LINKDIR)/intedit
	rm -f $(LINKDIR)/intxyz
	rm -f $(LINKDIR)/minimize
	rm -f $(LINKDIR)/minirot
	rm -f $(LINKDIR)/minrigid
	rm -f $(LINKDIR)/molxyz
	rm -f $(LINKDIR)/monte
	rm -f $(LINKDIR)/newton
	rm -f $(LINKDIR)/newtrot
	rm -f $(LINKDIR)/nucleic
	rm -f $(LINKDIR)/optimize
	rm -f $(LINKDIR)/optirot
	rm -f $(LINKDIR)/optrigid
	rm -f $(LINKDIR)/path
	rm -f $(LINKDIR)/pdbxyz
	rm -f $(LINKDIR)/polarize
	rm -f $(LINKDIR)/poledit
	rm -f $(LINKDIR)/potential
	rm -f $(LINKDIR)/prmedit
	rm -f $(LINKDIR)/protein
	rm -f $(LINKDIR)/pss
	rm -f $(LINKDIR)/pssrigid
	rm -f $(LINKDIR)/pssrot
	rm -f $(LINKDIR)/radial
	rm -f $(LINKDIR)/saddle
	rm -f $(LINKDIR)/scan
	rm -f $(LINKDIR)/sniffer
	rm -f $(LINKDIR)/spacefill
	rm -f $(LINKDIR)/spectrum
	rm -f $(LINKDIR)/superpose
	rm -f $(LINKDIR)/sybylxyz
	rm -f $(LINKDIR)/testgrad
	rm -f $(LINKDIR)/testhess
	rm -f $(LINKDIR)/testpair
	rm -f $(LINKDIR)/testpol
	rm -f $(LINKDIR)/testrot
	rm -f $(LINKDIR)/timer
	rm -f $(LINKDIR)/timerot
	rm -f $(LINKDIR)/torsfit
	rm -f $(LINKDIR)/valence
	rm -f $(LINKDIR)/vibbig
	rm -f $(LINKDIR)/vibrate
	rm -f $(LINKDIR)/vibrot
	rm -f $(LINKDIR)/xtalfit
	rm -f $(LINKDIR)/xtalmin
	rm -f $(LINKDIR)/xyzedit
	rm -f $(LINKDIR)/xyzint
	rm -f $(LINKDIR)/xyzpdb
	rm -f $(LINKDIR)/xyzsybyl

create_links:
	ln -s $(BINDIR)/alchemy    $(LINKDIR)/alchemy
	ln -s $(BINDIR)/analyze    $(LINKDIR)/analyze
	ln -s $(BINDIR)/anneal     $(LINKDIR)/anneal
	ln -s $(BINDIR)/archive    $(LINKDIR)/archive
	ln -s $(BINDIR)/bar        $(LINKDIR)/bar
	ln -s $(BINDIR)/correlate  $(LINKDIR)/correlate
	ln -s $(BINDIR)/crystal    $(LINKDIR)/crystal
	ln -s $(BINDIR)/diffuse    $(LINKDIR)/diffuse
	ln -s $(BINDIR)/distgeom   $(LINKDIR)/distgeom
	ln -s $(BINDIR)/document   $(LINKDIR)/document
	ln -s $(BINDIR)/dynamic    $(LINKDIR)/dynamic
	ln -s $(BINDIR)/gda        $(LINKDIR)/gda
	ln -s $(BINDIR)/intedit    $(LINKDIR)/intedit
	ln -s $(BINDIR)/intxyz     $(LINKDIR)/intxyz
	ln -s $(BINDIR)/minimize   $(LINKDIR)/minimize
	ln -s $(BINDIR)/minirot    $(LINKDIR)/minirot
	ln -s $(BINDIR)/minrigid   $(LINKDIR)/minrigid
	ln -s $(BINDIR)/molxyz     $(LINKDIR)/molxyz
	ln -s $(BINDIR)/monte      $(LINKDIR)/monte
	ln -s $(BINDIR)/newton     $(LINKDIR)/newton
	ln -s $(BINDIR)/newtrot    $(LINKDIR)/newtrot
	ln -s $(BINDIR)/nucleic    $(LINKDIR)/nucleic
	ln -s $(BINDIR)/optimize   $(LINKDIR)/optimize
	ln -s $(BINDIR)/optirot    $(LINKDIR)/optirot
	ln -s $(BINDIR)/optrigid   $(LINKDIR)/optrigid
	ln -s $(BINDIR)/path       $(LINKDIR)/path
	ln -s $(BINDIR)/pdbxyz     $(LINKDIR)/pdbxyz
	ln -s $(BINDIR)/polarize   $(LINKDIR)/polarize
	ln -s $(BINDIR)/poledit    $(LINKDIR)/poledit
	ln -s $(BINDIR)/potential  $(LINKDIR)/potential
	ln -s $(BINDIR)/prmedit    $(LINKDIR)/prmedit
	ln -s $(BINDIR)/protein    $(LINKDIR)/protein
	ln -s $(BINDIR)/pss        $(LINKDIR)/pss
	ln -s $(BINDIR)/pssrigid   $(LINKDIR)/pssrigid
	ln -s $(BINDIR)/pssrot     $(LINKDIR)/pssrot
	ln -s $(BINDIR)/radial     $(LINKDIR)/radial
	ln -s $(BINDIR)/saddle     $(LINKDIR)/saddle
	ln -s $(BINDIR)/scan       $(LINKDIR)/scan
	ln -s $(BINDIR)/sniffer    $(LINKDIR)/sniffer
	ln -s $(BINDIR)/spacefill  $(LINKDIR)/spacefill
	ln -s $(BINDIR)/spectrum   $(LINKDIR)/spectrum
	ln -s $(BINDIR)/superpose  $(LINKDIR)/superpose
	ln -s $(BINDIR)/sybylxyz   $(LINKDIR)/sybylxyz
	ln -s $(BINDIR)/testgrad   $(LINKDIR)/testgrad
	ln -s $(BINDIR)/testhess   $(LINKDIR)/testhess
	ln -s $(BINDIR)/testpair   $(LINKDIR)/testpair
	ln -s $(BINDIR)/testpol    $(LINKDIR)/testpol
	ln -s $(BINDIR)/testrot    $(LINKDIR)/testrot
	ln -s $(BINDIR)/timer      $(LINKDIR)/timer
	ln -s $(BINDIR)/timerot    $(LINKDIR)/timerot
	ln -s $(BINDIR)/torsfit    $(LINKDIR)/torsfit
	ln -s $(BINDIR)/valence    $(LINKDIR)/valence
	ln -s $(BINDIR)/vibbig     $(LINKDIR)/vibbig
	ln -s $(BINDIR)/vibrate    $(LINKDIR)/vibrate
	ln -s $(BINDIR)/vibrot     $(LINKDIR)/vibrot
	ln -s $(BINDIR)/xtalfit    $(LINKDIR)/xtalfit
	ln -s $(BINDIR)/xtalmin    $(LINKDIR)/xtalmin
	ln -s $(BINDIR)/xyzedit    $(LINKDIR)/xyzedit
	ln -s $(BINDIR)/xyzint     $(LINKDIR)/xyzint
	ln -s $(BINDIR)/xyzpdb     $(LINKDIR)/xyzpdb
	ln -s $(BINDIR)/xyzsybyl   $(LINKDIR)/xyzsybyl

libtinker.a: ${OBJS} 
	ar ${LIBFLAGS} libtinker.a \
        action.o \
        active.o \
        align.o \
        analysis.o \
        analyz.o \
        angang.o \
        angbnd.o \
        angles.o \
        angpot.o \
        angtor.o \
        argue.o \
        ascii.o \
        atmlst.o \
        atomid.o \
        atoms.o \
        attach.o \
        basefile.o \
        bath.o \
        beeman.o \
        bicubic.o \
        bitor.o \
        bitors.o \
        bndpot.o \
        bndstr.o \
        bonds.o \
        born.o \
        bound.o \
        bounds.o \
        boxes.o \
        bussi.o \
        calendar.o \
        cell.o \
        center.o \
        charge.o \
        chgpot.o \
        chkpole.o \
        chkring.o \
        chkxyz.o \
        cholesky.o \
        chrono.o \
        chunks.o \
        clock.o \
        cluster.o \
        column.o \
        command.o \
        connect.o \
        connolly.o \
        control.o \
        couple.o \
        cspline.o \
        cutoffs.o \
        deflate.o \
        delete.o \
        deriv.o \
        diagq.o \
        diffeq.o \
        dipole.o \
        disgeo.o \
        dma.o \
        domega.o \
        eangang.o \
        eangang1.o \
        eangang2.o \
        eangang3.o \
        eangle.o \
        eangle1.o \
        eangle2.o \
        eangle3.o \
        eangtor.o \
        eangtor1.o \
        eangtor2.o \
        eangtor3.o \
        ebond.o \
        ebond1.o \
        ebond2.o \
        ebond3.o \
        ebuck.o \
        ebuck1.o \
        ebuck2.o \
        ebuck3.o \
        echarge.o \
        echarge1.o \
        echarge2.o \
        echarge3.o \
        echgdpl.o \
        echgdpl1.o \
        echgdpl2.o \
        echgdpl3.o \
        edipole.o \
        edipole1.o \
        edipole2.o \
        edipole3.o \
        egauss.o \
        egauss1.o \
        egauss2.o \
        egauss3.o \
        egeom.o \
        egeom1.o \
        egeom2.o \
        egeom3.o \
        ehal.o \
        ehal1.o \
        ehal2.o \
        ehal3.o \
        eimprop.o \
        eimprop1.o \
        eimprop2.o \
        eimprop3.o \
        eimptor.o \
        eimptor1.o \
        eimptor2.o \
        eimptor3.o \
        elj.o \
        elj1.o \
        elj2.o \
        elj3.o \
        embed.o \
        emetal.o \
        emetal1.o \
        emetal2.o \
        emetal3.o \
        emm3hb.o \
        emm3hb1.o \
        emm3hb2.o \
        emm3hb3.o \
        empole.o \
        empole1.o \
        empole2.o \
        empole3.o \
        energi.o \
        energy.o \
        eopbend.o \
        eopbend1.o \
        eopbend2.o \
        eopbend3.o \
        eopdist.o \
        eopdist1.o \
        eopdist2.o \
        eopdist3.o \
        epitors.o \
        epitors1.o \
        epitors2.o \
        epitors3.o \
        erf.o \
        erxnfld.o \
        erxnfld1.o \
        erxnfld2.o \
        erxnfld3.o \
        esolv.o \
        esolv1.o \
        esolv2.o \
        esolv3.o \
        estrbnd.o \
        estrbnd1.o \
        estrbnd2.o \
        estrbnd3.o \
        estrtor.o \
        estrtor1.o \
        estrtor2.o \
        estrtor3.o \
        etors.o \
        etors1.o \
        etors2.o \
        etors3.o \
        etortor.o \
        etortor1.o \
        etortor2.o \
        etortor3.o \
        eurey.o \
        eurey1.o \
        eurey2.o \
        eurey3.o \
        evcorr.o \
        ewald.o \
        extra.o \
        extra1.o \
        extra2.o \
        extra3.o \
        faces.o \
        fatal.o \
        fft.o \
        fft3d.o \
        fftpack.o \
        field.o \
        fields.o \
        files.o \
        final.o \
        flatten.o \
        fracs.o \
        freeunit.o \
        freeze.o \
        geometry.o \
        getint.o \
        getkey.o \
        getmol.o \
        getmol2.o \
        getnumb.o \
        getpdb.o \
        getprm.o \
        getref.o \
        getstring.o \
        gettext.o \
        getword.o \
        getxyz.o \
        ghmcstep.o \
        gkstuf.o \
        gradient.o \
        gradrgd.o \
        gradrot.o \
        group.o \
        groups.o \
        grpline.o \
        gyrate.o \
        hescut.o \
        hessian.o \
        hessn.o \
        hessrgd.o \
        hessrot.o \
        hpmf.o \
        hybrid.o \
        image.o \
        impose.o \
        improp.o \
        imptor.o \
        induce.o \
        inertia.o \
        inform.o \
        initatom.o \
        initial.o \
        initprm.o \
        initres.o \
        initrot.o \
        insert.o \
        inter.o \
        invbeta.o \
        invert.o \
        iounit.o \
        jacobi.o \
        kanang.o \
        kangang.o \
        kangle.o \
        kangs.o \
        kangtor.o \
        kantor.o \
        katom.o \
        katoms.o \
        kbond.o \
        kbonds.o \
        kcharge.o \
        kchrge.o \
        kdipol.o \
        kdipole.o \
        kewald.o \
        kextra.o \
        keys.o \
        kgeom.o \
        khbond.o \
        kimprop.o \
        kimptor.o \
        kinetic.o \
        kiprop.o \
        kitors.o \
        kmetal.o \
        kmpole.o \
        kmulti.o \
        kopbend.o \
        kopbnd.o \
        kopdist.o \
        kopdst.o \
        korbit.o \
        korbs.o \
        kpitor.o \
        kpitors.o \
        kpolar.o \
        kpolr.o \
        ksolv.o \
        kstbnd.o \
        kstrbnd.o \
        kstrtor.o \
        ksttor.o \
        ktors.o \
        ktorsn.o \
        ktortor.o \
        ktrtor.o \
        kurey.o \
        kurybr.o \
        kvdw.o \
        kvdwpr.o \
        kvdws.o \
        lattice.o \
        lbfgs.o \
        light.o \
        lights.o \
        limits.o \
        linmin.o \
        makeint.o \
        makeref.o \
        makexyz.o \
        math.o \
        maxwell.o \
        mdinit.o \
        mdrest.o \
        mdsave.o \
        mdstat.o \
        mdstuf.o \
        mechanic.o \
        merck.o \
        merge.o \
        minima.o \
        molcul.o \
        moldyn.o \
        molecule.o \
        moment.o \
        moments.o \
        mplpot.o \
        mpole.o \
        mutant.o \
        mutate.o \
        nblist.o \
        neigh.o \
        nextarg.o \
        nexttext.o \
        nonpol.o \
        nose.o \
        nspline.o \
        nucleo.o \
        number.o \
        numeral.o \
        numgrad.o \
        ocvm.o \
        omega.o \
        opbend.o \
        opdist.o \
        openend.o \
        openmp.o \
        optsave.o \
        orbital.o \
        orbits.o \
        orient.o \
        orthog.o \
        output.o \
        overlap.o \
        params.o \
        paths.o \
        pbstuf.o \
        pdb.o \
        phipsi.o \
        picalc.o \
        piorbs.o \
        pistuf.o \
        pitors.o \
        pme.o \
        pmestuff.o \
        pmpb.o \
        polar.o \
        polgrp.o \
        polpot.o \
        polymer.o \
        potent.o \
        potfit.o \
        precis.o \
        precise.o \
        pressure.o \
        prmkey.o \
        promo.o \
        prtdyn.o \
        prterr.o \
        prtint.o \
        prtmol2.o \
        prtpdb.o \
        prtprm.o \
        prtseq.o \
        prtxyz.o \
        ptable.o \
        qmstuf.o \
        qrfact.o \
        quatfit.o \
        random.o \
        rattle.o \
        readdyn.o \
        readgau.o \
        readint.o \
        readmol.o \
        readmol2.o \
        readpdb.o \
        readprm.o \
        readseq.o \
        readxyz.o \
        refer.o \
        replica.o \
        resdue.o \
        respa.o \
        restrn.o \
        rgddyn.o \
        rgdstep.o \
        rigid.o \
        ring.o \
        rings.o \
        rmsfit.o \
        rotbnd.o \
        rotlist.o \
        rotpole.o \
        rxnfld.o \
        rxnpot.o \
        scales.o \
        sdstep.o \
        search.o \
        sequen.o \
        server.o \
        shakeup.o \
        shunt.o \
        sigmoid.o \
        sizes.o \
        sktstuff.o \
        socket.o \
        solute.o \
        sort.o \
        square.o \
        stodyn.o \
        strbnd.o \
        strtor.o \
        suffix.o \
        surface.o \
        surfatom.o \
        switch.o \
        syntrn.o \
        tarray.o \
        temper.o \
        titles.o \
        tncg.o \
        torphase.o \
        torpot.o \
        torque.o \
        tors.o \
        torsions.o \
        tortor.o \
        tree.o \
        trimtext.o \
        unitcell.o \
        units.o \
        uprior.o \
        urey.o \
        urypot.o \
        usage.o \
        usolve.o \
        valfit.o \
        vdw.o \
        vdwpot.o \
        verlet.o \
        version.o \
        vibs.o \
        virial.o \
        volume.o \
        warp.o \
        xtals.o \
        xyzatm.o \
        zatom.o \
        zclose.o \
        zcoord.o
	${RANLIB} libtinker.a

###############################################################
##  Next Section has Explicit Dependencies on Include Files  ##
###############################################################

action.o:
active.o: atoms.o inform.o iounit.o keys.o sizes.o usage.o
alchemy.o: analyz.o atoms.o energi.o files.o inform.o iounit.o katoms.o mutant.o potent.o sizes.o units.o usage.o
align.o:
analysis.o: analyz.o atoms.o bound.o energi.o group.o inter.o iounit.o limits.o potent.o sizes.o vdwpot.o
analyz.o:
analyze.o: action.o analyz.o angang.o angbnd.o angpot.o angtor.o atomid.o atoms.o bitor.o bndstr.o bound.o boxes.o charge.o chgpot.o couple.o dipole.o energi.o ewald.o fields.o files.o improp.o imptor.o inform.o inter.o iounit.o korbs.o ktrtor.o kvdws.o limits.o math.o molcul.o moment.o mpole.o opbend.o opdist.o piorbs.o pistuf.o pitors.o pme.o polar.o polgrp.o potent.o sizes.o solute.o strbnd.o strtor.o tors.o tortor.o units.o urey.o vdw.o vdwpot.o virial.o
angang.o:
angbnd.o:
angles.o: angbnd.o atmlst.o atoms.o couple.o iounit.o sizes.o
angpot.o:
angtor.o:
anneal.o: atomid.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o mdstuf.o potent.o sizes.o solute.o usage.o warp.o
archive.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o sizes.o titles.o usage.o
argue.o:
ascii.o:
atmlst.o:
atomid.o: sizes.o
atoms.o: sizes.o
attach.o: atoms.o couple.o iounit.o sizes.o
bar.o: atoms.o energi.o files.o inform.o iounit.o keys.o sizes.o units.o
basefile.o: ascii.o files.o
bath.o:
beeman.o: atomid.o atoms.o freeze.o mdstuf.o moldyn.o sizes.o units.o usage.o
bicubic.o:
bitor.o:
bitors.o: angbnd.o atoms.o bitor.o couple.o iounit.o sizes.o
bndpot.o:
bndstr.o:
bonds.o: atmlst.o atoms.o bndstr.o couple.o iounit.o sizes.o
born.o: atomid.o atoms.o bath.o chgpot.o couple.o deriv.o inform.o iounit.o math.o pbstuf.o sizes.o solute.o virial.o
bound.o:
bounds.o: atomid.o atoms.o boxes.o molcul.o sizes.o
boxes.o:
bussi.o: atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o sizes.o units.o usage.o
calendar.o:
cell.o:
center.o: align.o
charge.o:
chgpot.o:
chkpole.o: atoms.o mpole.o sizes.o
chkring.o: couple.o
chkxyz.o: atoms.o iounit.o sizes.o
cholesky.o:
chrono.o:
chunks.o:
clock.o: chrono.o
cluster.o: atomid.o atoms.o bound.o group.o inform.o iounit.o keys.o limits.o molcul.o sizes.o
column.o:
command.o: argue.o
connect.o: atoms.o couple.o sizes.o zclose.o zcoord.o
connolly.o: atoms.o faces.o inform.o iounit.o math.o sizes.o
control.o: argue.o inform.o keys.o output.o
correlate.o: ascii.o atomid.o atoms.o files.o inform.o iounit.o sizes.o
couple.o: sizes.o
crystal.o: atomid.o atoms.o bound.o boxes.o couple.o files.o iounit.o math.o molcul.o sizes.o
cspline.o: iounit.o
cutoffs.o: atoms.o bound.o hescut.o keys.o limits.o neigh.o polpot.o sizes.o tarray.o
deflate.o: iounit.o
delete.o: atomid.o atoms.o couple.o inform.o iounit.o sizes.o
deriv.o:
diagq.o:
diffeq.o: atoms.o iounit.o math.o sizes.o warp.o
diffuse.o: atomid.o atoms.o bound.o inform.o iounit.o molcul.o sizes.o usage.o
dipole.o:
disgeo.o:
distgeom.o: angbnd.o atomid.o atoms.o bndstr.o couple.o disgeo.o files.o inform.o iounit.o kvdws.o math.o refer.o restrn.o sizes.o tors.o
dma.o:
document.o: iounit.o
domega.o:
dynamic.o: atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o sizes.o solute.o stodyn.o usage.o
eangang.o: angang.o angbnd.o angpot.o atoms.o bound.o energi.o group.o math.o sizes.o usage.o
eangang1.o: angang.o angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o usage.o virial.o
eangang2.o: angang.o angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o sizes.o
eangang3.o: action.o analyz.o angang.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o usage.o
eangle.o: angbnd.o angpot.o atoms.o bound.o energi.o group.o math.o sizes.o usage.o
eangle1.o: angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o usage.o virial.o
eangle2.o: angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o sizes.o
eangle3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o usage.o
eangtor.o: angbnd.o angtor.o atoms.o bound.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o
eangtor1.o: angbnd.o angtor.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o virial.o
eangtor2.o: angbnd.o angtor.o atoms.o bound.o group.o hessn.o math.o sizes.o torpot.o tors.o
eangtor3.o: action.o analyz.o angbnd.o angtor.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o torpot.o tors.o usage.o
ebond.o: atoms.o bndpot.o bndstr.o bound.o energi.o group.o sizes.o usage.o
ebond1.o: atoms.o bndpot.o bndstr.o bound.o deriv.o energi.o group.o sizes.o usage.o virial.o
ebond2.o: atmlst.o atoms.o bndpot.o bndstr.o bound.o couple.o group.o hessn.o sizes.o
ebond3.o: action.o analyz.o atomid.o atoms.o bndpot.o bndstr.o bound.o energi.o group.o inform.o iounit.o sizes.o usage.o
ebuck.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o iounit.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
ebuck1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
ebuck2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o iounit.o math.o shunt.o sizes.o vdw.o vdwpot.o warp.o
ebuck3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
echarge.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o energi.o ewald.o group.o iounit.o light.o limits.o math.o neigh.o pme.o shunt.o sizes.o usage.o warp.o
echarge1.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o light.o limits.o math.o molcul.o neigh.o pme.o shunt.o sizes.o usage.o virial.o warp.o
echarge2.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o ewald.o group.o hessn.o limits.o math.o shunt.o sizes.o warp.o
echarge3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o energi.o ewald.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o warp.o
echgdpl.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o energi.o group.o shunt.o sizes.o units.o usage.o
echgdpl1.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o deriv.o dipole.o energi.o group.o inter.o molcul.o shunt.o sizes.o units.o usage.o virial.o
echgdpl2.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o group.o hessn.o shunt.o sizes.o units.o
echgdpl3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o energi.o group.o inform.o inter.o iounit.o molcul.o shunt.o sizes.o units.o usage.o
edipole.o: atoms.o bound.o cell.o chgpot.o dipole.o energi.o group.o shunt.o sizes.o units.o usage.o
edipole1.o: atoms.o bound.o cell.o chgpot.o deriv.o dipole.o energi.o group.o inter.o molcul.o shunt.o sizes.o units.o usage.o virial.o
edipole2.o: atoms.o bound.o cell.o chgpot.o dipole.o group.o hessn.o shunt.o sizes.o units.o
edipole3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o chgpot.o dipole.o energi.o group.o inform.o inter.o iounit.o molcul.o shunt.o sizes.o units.o usage.o
egauss.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
egauss1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
egauss2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o warp.o
egauss3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
egeom.o: atomid.o atoms.o bound.o energi.o group.o math.o molcul.o restrn.o sizes.o usage.o
egeom1.o: atomid.o atoms.o bound.o deriv.o energi.o group.o inter.o math.o molcul.o restrn.o sizes.o usage.o virial.o
egeom2.o: atomid.o atoms.o bound.o deriv.o group.o hessn.o math.o molcul.o restrn.o sizes.o
egeom3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o inter.o iounit.o math.o molcul.o restrn.o sizes.o usage.o
ehal.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
ehal1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o
ehal2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o
ehal3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
eimprop.o: atoms.o bound.o energi.o group.o improp.o math.o sizes.o torpot.o usage.o
eimprop1.o: atoms.o bound.o deriv.o energi.o group.o improp.o math.o sizes.o torpot.o usage.o virial.o
eimprop2.o: atoms.o bound.o group.o hessn.o improp.o math.o sizes.o torpot.o
eimprop3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o improp.o inform.o iounit.o math.o sizes.o torpot.o usage.o
eimptor.o: atoms.o bound.o energi.o group.o imptor.o sizes.o torpot.o usage.o
eimptor1.o: atoms.o bound.o deriv.o energi.o group.o imptor.o sizes.o torpot.o usage.o virial.o
eimptor2.o: atoms.o bound.o group.o hessn.o imptor.o sizes.o torpot.o
eimptor3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o imptor.o inform.o iounit.o math.o sizes.o torpot.o usage.o
elj.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
elj1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
elj2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o math.o shunt.o sizes.o vdw.o vdwpot.o warp.o
elj3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
embed.o: angbnd.o atoms.o bndstr.o couple.o disgeo.o files.o inform.o iounit.o keys.o light.o math.o minima.o output.o refer.o restrn.o sizes.o tors.o units.o
emetal.o: atomid.o atoms.o couple.o energi.o kchrge.o sizes.o
emetal1.o: atomid.o atoms.o couple.o deriv.o energi.o kchrge.o sizes.o
emetal2.o:
emetal3.o: action.o analyz.o atomid.o atoms.o energi.o kchrge.o sizes.o
emm3hb.o: atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o energi.o group.o light.o limits.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
emm3hb1.o: atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o
emm3hb2.o: atmlst.o atomid.o atoms.o bndstr.o bound.o cell.o chgpot.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o
emm3hb3.o: action.o analyz.o atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
empole.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o energi.o ewald.o group.o limits.o math.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o
empole1.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
empole2.o: atoms.o bound.o boxes.o chgpot.o couple.o deriv.o group.o hessn.o limits.o molcul.o mplpot.o mpole.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o
empole3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o chgpot.o couple.o energi.o ewald.o group.o inform.o inter.o iounit.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o
energi.o:
energy.o: bound.o energi.o iounit.o limits.o potent.o rigid.o sizes.o vdwpot.o
eopbend.o: angbnd.o angpot.o atoms.o bound.o energi.o fields.o group.o math.o opbend.o sizes.o usage.o
eopbend1.o: angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o opbend.o sizes.o usage.o virial.o
eopbend2.o: angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o opbend.o sizes.o
eopbend3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o opbend.o sizes.o usage.o
eopdist.o: angpot.o atoms.o bound.o energi.o group.o opdist.o sizes.o usage.o
eopdist1.o: angpot.o atoms.o bound.o deriv.o energi.o group.o opdist.o sizes.o usage.o virial.o
eopdist2.o: angpot.o atoms.o bound.o group.o hessn.o opdist.o sizes.o usage.o
eopdist3.o: action.o analyz.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o opdist.o sizes.o usage.o
epitors.o: atoms.o bound.o energi.o group.o pitors.o sizes.o torpot.o usage.o
epitors1.o: atoms.o bound.o deriv.o energi.o group.o pitors.o sizes.o torpot.o usage.o virial.o
epitors2.o: angbnd.o atoms.o bound.o group.o hessn.o pitors.o sizes.o torpot.o usage.o
epitors3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o pitors.o sizes.o torpot.o usage.o
erf.o: iounit.o math.o
erxnfld.o: atoms.o chgpot.o energi.o mpole.o rxnfld.o rxnpot.o shunt.o sizes.o usage.o
erxnfld1.o: atoms.o deriv.o energi.o sizes.o
erxnfld2.o:
erxnfld3.o: action.o analyz.o atomid.o atoms.o chgpot.o energi.o inform.o iounit.o mpole.o shunt.o sizes.o usage.o
esolv.o: atomid.o atoms.o bound.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o kvdws.o limits.o math.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o warp.o
esolv1.o: atomid.o atoms.o bound.o boxes.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o inter.o iounit.o kvdws.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o virial.o warp.o
esolv2.o: atoms.o charge.o chgpot.o hessn.o math.o potent.o shunt.o sizes.o solute.o warp.o
esolv3.o: action.o analyz.o atomid.o atoms.o bound.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o inform.o inter.o iounit.o kvdws.o limits.o math.o molcul.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o warp.o
estrbnd.o: angbnd.o angpot.o atoms.o bndstr.o bound.o energi.o group.o math.o sizes.o strbnd.o usage.o
estrbnd1.o: angbnd.o angpot.o atoms.o bndstr.o bound.o deriv.o energi.o group.o math.o sizes.o strbnd.o usage.o virial.o
estrbnd2.o: angbnd.o angpot.o atoms.o bndstr.o bound.o group.o hessn.o math.o sizes.o strbnd.o
estrbnd3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bndstr.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o strbnd.o usage.o
estrtor.o: atoms.o bndstr.o bound.o energi.o group.o sizes.o strtor.o torpot.o tors.o usage.o
estrtor1.o: atoms.o bndstr.o bound.o deriv.o energi.o group.o sizes.o strtor.o torpot.o tors.o usage.o virial.o
estrtor2.o: atoms.o bndstr.o bound.o group.o hessn.o sizes.o strtor.o torpot.o tors.o
estrtor3.o: action.o analyz.o atomid.o atoms.o bndstr.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o strtor.o torpot.o tors.o usage.o
etors.o: atoms.o bound.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o warp.o
etors1.o: atoms.o bound.o deriv.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o virial.o warp.o
etors2.o: atoms.o bound.o group.o hessn.o math.o sizes.o torpot.o tors.o warp.o
etors3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o torpot.o tors.o usage.o warp.o
etortor.o: atomid.o atoms.o bitor.o bound.o couple.o energi.o group.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o
etortor1.o: atoms.o bitor.o bound.o deriv.o energi.o group.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o virial.o
etortor2.o: atoms.o bitor.o bound.o group.o hessn.o ktrtor.o math.o sizes.o torpot.o tortor.o units.o
etortor3.o: action.o analyz.o atoms.o bitor.o bound.o energi.o group.o inform.o iounit.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o
eurey.o: atoms.o bound.o energi.o group.o sizes.o urey.o urypot.o usage.o
eurey1.o: atoms.o bound.o deriv.o energi.o group.o sizes.o urey.o urypot.o usage.o virial.o
eurey2.o: atoms.o bound.o couple.o group.o hessn.o sizes.o urey.o urypot.o
eurey3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o sizes.o urey.o urypot.o usage.o
evcorr.o: atoms.o bound.o boxes.o limits.o math.o mutant.o shunt.o sizes.o vdw.o vdwpot.o
ewald.o:
extra.o: energi.o
extra1.o: atoms.o deriv.o energi.o sizes.o
extra2.o: atoms.o hessn.o sizes.o
extra3.o: action.o analyz.o atoms.o energi.o sizes.o
faces.o:
fatal.o: iounit.o
fft.o:
fft3d.o: fft.o openmp.o pme.o
fftpack.o: math.o
field.o: keys.o potent.o
fields.o:
files.o:
final.o: align.o analyz.o angang.o angbnd.o angtor.o atmlst.o bitor.o bndstr.o cell.o charge.o chunks.o couple.o deriv.o dipole.o disgeo.o domega.o faces.o fracs.o freeze.o group.o hessn.o hpmf.o improp.o imptor.o inform.o iounit.o katoms.o kchrge.o kpolr.o kvdws.o light.o merck.o molcul.o moldyn.o mpole.o mutant.o neigh.o nonpol.o omega.o opbend.o opdist.o orbits.o paths.o pbstuf.o pdb.o piorbs.o pistuf.o pitors.o pme.o polar.o polgrp.o qmstuf.o refer.o restrn.o rgddyn.o rigid.o ring.o rotbnd.o sizes.o socket.o solute.o stodyn.o strbnd.o strtor.o syntrn.o tarray.o tors.o tortor.o uprior.o urey.o usage.o usolve.o vdw.o vibs.o warp.o
flatten.o: atoms.o fields.o inform.o iounit.o keys.o sizes.o warp.o
fracs.o:
freeunit.o: iounit.o
freeze.o:
gda.o: atoms.o files.o iounit.o minima.o potent.o sizes.o vdwpot.o warp.o
geometry.o: atoms.o math.o sizes.o
getint.o: atoms.o inform.o iounit.o output.o sizes.o
getkey.o: argue.o files.o iounit.o keys.o openmp.o
getmol.o: files.o iounit.o
getmol2.o: files.o iounit.o
getnumb.o: ascii.o
getpdb.o: iounit.o
getprm.o: files.o iounit.o keys.o params.o
getref.o: atomid.o atoms.o boxes.o couple.o files.o refer.o sizes.o titles.o
getstring.o: ascii.o
gettext.o: ascii.o
getword.o: ascii.o
getxyz.o: inform.o iounit.o output.o
ghmcstep.o: atomid.o atoms.o bath.o freeze.o iounit.o mdstuf.o moldyn.o sizes.o stodyn.o units.o usage.o virial.o
gkstuf.o: sizes.o
gradient.o: atoms.o bound.o couple.o deriv.o energi.o inter.o iounit.o limits.o potent.o rigid.o sizes.o vdwpot.o virial.o
gradrgd.o: atoms.o group.o rigid.o sizes.o
gradrot.o: atoms.o deriv.o domega.o omega.o potent.o rotbnd.o sizes.o
group.o:
groups.o: group.o
grpline.o: atomid.o atoms.o group.o rgddyn.o sizes.o
gyrate.o: atoms.o sizes.o usage.o
hescut.o:
hessian.o: atoms.o bound.o couple.o hescut.o hessn.o inform.o iounit.o limits.o mpole.o potent.o rigid.o sizes.o usage.o vdw.o vdwpot.o
hessn.o:
hessrgd.o: atoms.o group.o rigid.o sizes.o
hessrot.o: math.o omega.o sizes.o zcoord.o
hpmf.o:
hybrid.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o charge.o couple.o dipole.o imptor.o inform.o iounit.o kangs.o katoms.o kbonds.o kchrge.o kdipol.o kitors.o kstbnd.o ksttor.o ktorsn.o kvdws.o math.o mutant.o sizes.o strbnd.o strtor.o tors.o vdw.o vdwpot.o
image.o: boxes.o cell.o sizes.o
impose.o: align.o inform.o iounit.o sizes.o
improp.o:
imptor.o:
induce.o: atoms.o bound.o boxes.o cell.o couple.o ewald.o gkstuf.o group.o inform.o iounit.o limits.o math.o mpole.o neigh.o openmp.o pbstuf.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o tarray.o units.o uprior.o usolve.o
inertia.o: atomid.o atoms.o iounit.o math.o sizes.o
inform.o:
initatom.o: ptable.o
initial.o: align.o atoms.o bath.o bound.o boxes.o cell.o files.o group.o inform.o iounit.o keys.o linmin.o minima.o molcul.o mutant.o neigh.o openmp.o output.o params.o pdb.o precis.o rigid.o scales.o sequen.o sizes.o socket.o warp.o zclose.o
initprm.o: angpot.o bndpot.o chgpot.o fields.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o math.o merck.o mplpot.o polpot.o rxnpot.o solute.o torpot.o units.o urypot.o vdwpot.o
initres.o: resdue.o
initrot.o: atoms.o couple.o group.o inform.o iounit.o math.o omega.o potent.o restrn.o rotbnd.o sizes.o usage.o zcoord.o
insert.o: atomid.o atoms.o couple.o inform.o iounit.o sizes.o
intedit.o: atomid.o atoms.o files.o iounit.o katoms.o sizes.o zcoord.o
inter.o:
intxyz.o: files.o iounit.o titles.o
invbeta.o:
invert.o: iounit.o
iounit.o:
jacobi.o: iounit.o
kanang.o: sizes.o
kangang.o: angang.o angbnd.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o kanang.o keys.o potent.o sizes.o tors.o
kangle.o: angbnd.o angpot.o atomid.o atoms.o bndstr.o couple.o fields.o inform.o iounit.o kangs.o keys.o merck.o potent.o ring.o sizes.o usage.o
kangs.o:
kangtor.o: angtor.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o kantor.o keys.o potent.o sizes.o tors.o
kantor.o:
katom.o: atomid.o atoms.o couple.o inform.o iounit.o katoms.o keys.o sizes.o
katoms.o: sizes.o
kbond.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o couple.o fields.o inform.o iounit.o kbonds.o keys.o merck.o potent.o sizes.o tors.o usage.o
kbonds.o:
kcharge.o: atomid.o atoms.o charge.o chgpot.o couple.o fields.o inform.o iounit.o kchrge.o keys.o merck.o potent.o sizes.o
kchrge.o: sizes.o
kdipol.o:
kdipole.o: atmlst.o atoms.o bndstr.o couple.o dipole.o inform.o iounit.o kdipol.o keys.o potent.o sizes.o
kewald.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o math.o openmp.o pme.o sizes.o
kextra.o:
keys.o:
kgeom.o: atomid.o atoms.o bound.o couple.o group.o iounit.o keys.o molcul.o potent.o restrn.o sizes.o
khbond.o:
kimprop.o: atomid.o atoms.o couple.o improp.o inform.o iounit.o keys.o kiprop.o potent.o sizes.o tors.o
kimptor.o: atomid.o atoms.o couple.o imptor.o inform.o iounit.o keys.o kitors.o math.o potent.o sizes.o tors.o
kinetic.o: atomid.o atoms.o bath.o group.o mdstuf.o moldyn.o rgddyn.o sizes.o units.o usage.o
kiprop.o:
kitors.o:
kmetal.o:
kmpole.o: atoms.o couple.o inform.o iounit.o keys.o kmulti.o mpole.o polar.o polgrp.o potent.o sizes.o units.o
kmulti.o:
kopbend.o: angbnd.o angpot.o atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o kopbnd.o merck.o opbend.o potent.o sizes.o usage.o
kopbnd.o:
kopdist.o: angbnd.o angpot.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o keys.o kopdst.o opdist.o potent.o sizes.o
kopdst.o:
korbit.o: atomid.o atoms.o bndstr.o inform.o iounit.o keys.o korbs.o orbits.o piorbs.o pistuf.o sizes.o tors.o units.o
korbs.o: sizes.o
kpitor.o:
kpitors.o: atomid.o atoms.o bndstr.o couple.o inform.o iounit.o keys.o kpitor.o pitors.o potent.o sizes.o tors.o
kpolar.o: atoms.o couple.o inform.o iounit.o keys.o kpolr.o mpole.o polar.o polgrp.o polpot.o potent.o sizes.o usolve.o
kpolr.o: sizes.o
ksolv.o: angbnd.o atmlst.o atomid.o atoms.o bath.o bndstr.o chgpot.o couple.o gkstuf.o hpmf.o inform.o iounit.o keys.o kvdws.o math.o nonpol.o pbstuf.o potent.o ptable.o sizes.o solute.o
kstbnd.o:
kstrbnd.o: angbnd.o angpot.o atmlst.o atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o kstbnd.o merck.o potent.o ring.o sizes.o strbnd.o
kstrtor.o: atmlst.o atomid.o atoms.o couple.o inform.o iounit.o keys.o ksttor.o potent.o sizes.o strtor.o tors.o
ksttor.o:
ktors.o: atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o ktorsn.o math.o merck.o potent.o ring.o sizes.o tors.o usage.o
ktorsn.o:
ktortor.o: atomid.o atoms.o bitor.o inform.o iounit.o keys.o ktrtor.o potent.o sizes.o tortor.o
ktrtor.o:
kurey.o: angbnd.o atomid.o atoms.o inform.o iounit.o keys.o kurybr.o potent.o sizes.o urey.o
kurybr.o:
kvdw.o: atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o khbond.o kvdwpr.o kvdws.o math.o merck.o potent.o sizes.o vdw.o vdwpot.o
kvdwpr.o:
kvdws.o: sizes.o
lattice.o: boxes.o cell.o inform.o iounit.o math.o
lbfgs.o: inform.o iounit.o keys.o linmin.o math.o minima.o output.o scales.o
light.o:
lights.o: bound.o boxes.o cell.o iounit.o light.o sizes.o
limits.o:
linmin.o:
makeint.o: atoms.o couple.o inform.o iounit.o math.o sizes.o zclose.o zcoord.o
makeref.o: atomid.o atoms.o boxes.o couple.o files.o refer.o sizes.o titles.o
makexyz.o: atoms.o sizes.o zcoord.o
math.o:
maxwell.o: units.o
mdinit.o: atomid.o atoms.o bath.o bound.o couple.o files.o freeze.o group.o inform.o iounit.o keys.o mdstuf.o molcul.o moldyn.o mpole.o rgddyn.o rigid.o sizes.o stodyn.o units.o uprior.o usage.o
mdrest.o: atomid.o atoms.o bound.o group.o inform.o iounit.o mdstuf.o moldyn.o rgddyn.o sizes.o units.o
mdsave.o: atomid.o atoms.o bound.o boxes.o energi.o files.o group.o inform.o iounit.o mdstuf.o moldyn.o mpole.o output.o polar.o potent.o rgddyn.o sizes.o socket.o titles.o units.o
mdstat.o: atoms.o bath.o bound.o boxes.o inform.o inter.o iounit.o limits.o mdstuf.o molcul.o sizes.o units.o usage.o warp.o
mdstuf.o:
mechanic.o: inform.o iounit.o limits.o potent.o vdwpot.o
merck.o: sizes.o
merge.o: atomid.o atoms.o couple.o iounit.o refer.o sizes.o
minima.o:
minimize.o: atoms.o files.o inform.o iounit.o keys.o scales.o sizes.o usage.o
minirot.o: files.o inform.o iounit.o keys.o math.o omega.o scales.o sizes.o zcoord.o
minrigid.o: files.o group.o inform.o iounit.o keys.o math.o output.o rigid.o sizes.o
molcul.o:
moldyn.o:
molecule.o: atomid.o atoms.o couple.o molcul.o sizes.o
molxyz.o: files.o iounit.o titles.o
moment.o:
moments.o: atomid.o atoms.o bound.o charge.o dipole.o limits.o moment.o mpole.o polar.o potent.o rigid.o sizes.o solute.o units.o usage.o
monte.o: atoms.o files.o inform.o iounit.o omega.o output.o sizes.o units.o usage.o zcoord.o
mplpot.o:
mpole.o:
mutant.o:
mutate.o: atomid.o atoms.o charge.o inform.o iounit.o katoms.o keys.o mpole.o mutant.o polar.o potent.o sizes.o
nblist.o: atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o
neigh.o:
newton.o: atoms.o files.o inform.o iounit.o keys.o sizes.o usage.o
newtrot.o: files.o hescut.o inform.o iounit.o keys.o math.o omega.o sizes.o zcoord.o
nextarg.o: argue.o
nexttext.o:
nonpol.o:
nose.o: atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o sizes.o units.o usage.o virial.o
nspline.o:
nucleic.o: atoms.o couple.o files.o group.o inform.o iounit.o katoms.o math.o molcul.o nucleo.o output.o potent.o resdue.o restrn.o rigid.o sequen.o sizes.o titles.o usage.o
nucleo.o: sizes.o
number.o: inform.o iounit.o
numeral.o:
numgrad.o: atoms.o sizes.o
ocvm.o: inform.o iounit.o keys.o linmin.o math.o minima.o output.o potent.o scales.o
omega.o:
opbend.o:
opdist.o:
openend.o:
openmp.o:
optimize.o: atoms.o files.o inform.o iounit.o keys.o scales.o sizes.o usage.o
optirot.o: files.o inform.o iounit.o keys.o math.o omega.o scales.o sizes.o zcoord.o
optrigid.o: files.o group.o inform.o iounit.o keys.o math.o output.o rigid.o
optsave.o: atoms.o files.o iounit.o math.o omega.o output.o scales.o sizes.o socket.o usage.o zcoord.o
orbital.o: atomid.o atoms.o bndstr.o couple.o iounit.o keys.o piorbs.o potent.o sizes.o tors.o
orbits.o:
orient.o: atomid.o atoms.o group.o math.o rigid.o sizes.o
orthog.o:
output.o:
overlap.o: units.o
params.o:
path.o: align.o atomid.o atoms.o files.o inform.o iounit.o linmin.o minima.o output.o paths.o sizes.o
paths.o:
pbstuf.o:
pdb.o:
pdbxyz.o: atomid.o atoms.o couple.o fields.o files.o inform.o iounit.o katoms.o pdb.o resdue.o sequen.o sizes.o titles.o
phipsi.o: sizes.o
picalc.o: atomid.o atoms.o bndstr.o couple.o inform.o iounit.o orbits.o piorbs.o pistuf.o sizes.o tors.o units.o
piorbs.o:
pistuf.o:
pitors.o:
pme.o:
pmestuff.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
pmpb.o: iounit.o
polar.o:
polarize.o: atoms.o inform.o iounit.o molcul.o mpole.o polar.o polgrp.o polpot.o potent.o sizes.o units.o
poledit.o: atomid.o atoms.o couple.o dma.o files.o iounit.o keys.o kpolr.o mpole.o polar.o polgrp.o polpot.o potent.o sizes.o units.o
polgrp.o:
polpot.o:
polymer.o: atoms.o bndstr.o bound.o boxes.o iounit.o keys.o sizes.o
potent.o:
potential.o: atomid.o atoms.o bndstr.o charge.o chgpot.o couple.o dipole.o files.o inform.o iounit.o katoms.o kchrge.o kdipol.o keys.o kmulti.o kpolr.o math.o minima.o moment.o mpole.o neigh.o output.o polar.o potent.o potfit.o ptable.o refer.o sizes.o titles.o units.o
potfit.o: sizes.o
precis.o:
precise.o:
pressure.o: atomid.o atoms.o bath.o bound.o boxes.o group.o iounit.o math.o mdstuf.o molcul.o sizes.o units.o usage.o virial.o
prmedit.o: angpot.o bndpot.o iounit.o math.o params.o sizes.o urypot.o vdwpot.o
prmkey.o: angpot.o bndpot.o chgpot.o fields.o mplpot.o polpot.o potent.o rxnpot.o torpot.o urypot.o vdwpot.o
promo.o: iounit.o
protein.o: atomid.o atoms.o couple.o files.o group.o inform.o iounit.o katoms.o math.o molcul.o output.o phipsi.o potent.o resdue.o restrn.o rigid.o sequen.o sizes.o titles.o usage.o
prtdyn.o: atoms.o boxes.o files.o group.o mdstuf.o moldyn.o rgddyn.o sizes.o titles.o
prterr.o: files.o output.o
prtint.o: atomid.o atoms.o files.o inform.o sizes.o titles.o zclose.o zcoord.o
prtmol2.o: atomid.o atoms.o bndstr.o couple.o files.o iounit.o sizes.o titles.o
prtpdb.o: files.o pdb.o sequen.o sizes.o titles.o
prtprm.o: angpot.o bndpot.o chgpot.o fields.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o mplpot.o polpot.o sizes.o urypot.o vdwpot.o
prtseq.o: files.o sequen.o sizes.o
prtxyz.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o sizes.o titles.o
pss.o: atoms.o files.o hescut.o inform.o iounit.o math.o omega.o refer.o sizes.o tree.o warp.o zcoord.o
pssrigid.o: atoms.o files.o group.o inform.o iounit.o math.o minima.o molcul.o refer.o rigid.o sizes.o warp.o
pssrot.o: atoms.o files.o inform.o iounit.o math.o minima.o omega.o refer.o sizes.o warp.o zcoord.o
ptable.o:
qmstuf.o:
qrfact.o:
quatfit.o: align.o sizes.o
radial.o: argue.o atomid.o atoms.o bound.o boxes.o files.o inform.o iounit.o limits.o math.o molcul.o potent.o sizes.o
random.o: inform.o iounit.o keys.o math.o
rattle.o: atomid.o atoms.o freeze.o group.o inform.o iounit.o moldyn.o sizes.o units.o usage.o virial.o
readdyn.o: atoms.o boxes.o files.o group.o iounit.o mdstuf.o moldyn.o rgddyn.o sizes.o
readgau.o: ascii.o iounit.o qmstuf.o sizes.o units.o
readint.o: atomid.o atoms.o files.o inform.o iounit.o sizes.o titles.o zclose.o zcoord.o
readmol.o: atomid.o atoms.o couple.o files.o iounit.o ptable.o sizes.o titles.o
readmol2.o: atomid.o atoms.o couple.o files.o iounit.o sizes.o titles.o
readpdb.o: files.o inform.o iounit.o pdb.o resdue.o sequen.o titles.o
readprm.o: fields.o iounit.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o merck.o params.o sizes.o
readseq.o: files.o iounit.o resdue.o sequen.o
readxyz.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o sizes.o titles.o
refer.o: sizes.o
replica.o: bound.o boxes.o cell.o inform.o iounit.o
resdue.o:
respa.o: atomid.o atoms.o freeze.o limits.o moldyn.o potent.o sizes.o units.o usage.o virial.o
restrn.o:
rgddyn.o:
rgdstep.o: atomid.o atoms.o bound.o group.o iounit.o rgddyn.o sizes.o units.o virial.o
rigid.o:
ring.o:
rings.o: angbnd.o atoms.o bitor.o bndstr.o couple.o inform.o iounit.o ring.o tors.o
rmsfit.o: align.o sizes.o
rotbnd.o:
rotlist.o: atoms.o couple.o iounit.o molcul.o rotbnd.o sizes.o zclose.o
rotpole.o: atoms.o mpole.o sizes.o
rxnfld.o:
rxnpot.o:
saddle.o: atoms.o inform.o iounit.o keys.o linmin.o minima.o sizes.o syntrn.o titles.o zcoord.o
scales.o: sizes.o
scan.o: atoms.o files.o inform.o iounit.o math.o minima.o omega.o output.o sizes.o zcoord.o
sdstep.o: atomid.o atoms.o bath.o couple.o freeze.o kvdws.o math.o mdstuf.o moldyn.o sizes.o stodyn.o units.o usage.o virial.o
search.o: linmin.o math.o sizes.o
sequen.o: sizes.o
server.o:
shakeup.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o bound.o couple.o freeze.o keys.o math.o molcul.o ring.o sizes.o usage.o
shunt.o:
sigmoid.o:
sizes.o:
sktstuff.o: atomid.o atoms.o charge.o couple.o deriv.o fields.o files.o inform.o iounit.o keys.o moldyn.o mpole.o polar.o potent.o sizes.o socket.o
sniffer.o: atoms.o files.o inform.o iounit.o linmin.o math.o minima.o output.o scales.o sizes.o usage.o
socket.o:
solute.o:
sort.o:
spacefill.o: atomid.o atoms.o files.o inform.o iounit.o kvdws.o math.o sizes.o usage.o
spectrum.o: files.o iounit.o math.o units.o
square.o: inform.o iounit.o keys.o minima.o sizes.o
stodyn.o:
strbnd.o:
strtor.o:
suffix.o: ascii.o
superpose.o: align.o atomid.o atoms.o bound.o files.o inform.o iounit.o sizes.o titles.o
surface.o: atoms.o inform.o iounit.o math.o sizes.o usage.o
surfatom.o: atoms.o iounit.o math.o sizes.o
switch.o: limits.o nonpol.o shunt.o
sybylxyz.o: files.o iounit.o titles.o
syntrn.o:
tarray.o:
temper.o: atomid.o atoms.o bath.o group.o mdstuf.o molcul.o moldyn.o rgddyn.o sizes.o units.o usage.o
testgrad.o: atoms.o deriv.o energi.o files.o inform.o inter.o iounit.o sizes.o solute.o usage.o
testhess.o: atoms.o files.o hescut.o inform.o iounit.o sizes.o usage.o
testpair.o: atoms.o deriv.o energi.o inform.o iounit.o light.o limits.o neigh.o polpot.o potent.o sizes.o tarray.o vdwpot.o
testpol.o: atoms.o bound.o inform.o iounit.o limits.o polar.o polpot.o potent.o rigid.o sizes.o units.o
testrot.o: domega.o energi.o inform.o iounit.o math.o omega.o sizes.o zcoord.o
timer.o: atoms.o hescut.o inform.o iounit.o limits.o sizes.o
timerot.o: iounit.o limits.o omega.o sizes.o
titles.o:
tncg.o: atoms.o hescut.o inform.o iounit.o keys.o linmin.o math.o minima.o output.o piorbs.o potent.o sizes.o
torphase.o:
torpot.o:
torque.o: atoms.o deriv.o mpole.o sizes.o
tors.o:
torsfit.o: atomid.o atoms.o files.o inform.o iounit.o keys.o ktorsn.o math.o output.o potent.o qmstuf.o restrn.o scales.o sizes.o tors.o usage.o
torsions.o: atoms.o bndstr.o couple.o iounit.o sizes.o tors.o
tortor.o:
tree.o:
trimtext.o:
unitcell.o: bound.o boxes.o iounit.o keys.o
units.o:
uprior.o:
urey.o:
urypot.o:
usage.o:
usolve.o:
valence.o: angbnd.o angpot.o atomid.o atoms.o bndpot.o bndstr.o couple.o files.o hescut.o inform.o iounit.o kangs.o kbonds.o keys.o kopbnd.o kstbnd.o ktorsn.o kurybr.o kvdws.o linmin.o math.o minima.o opbend.o output.o potent.o qmstuf.o scales.o sizes.o strbnd.o torpot.o tors.o units.o urey.o urypot.o usage.o valfit.o vdwpot.o
valfit.o:
vdw.o:
vdwpot.o:
verlet.o: atomid.o atoms.o freeze.o moldyn.o sizes.o units.o usage.o
version.o: iounit.o output.o
vibbig.o: atomid.o atoms.o bound.o couple.o files.o hescut.o hessn.o inform.o iounit.o keys.o limits.o mpole.o potent.o rigid.o sizes.o units.o usage.o vdw.o vdwpot.o vibs.o
vibrate.o: atomid.o atoms.o files.o hescut.o iounit.o math.o sizes.o units.o usage.o
vibrot.o: iounit.o omega.o sizes.o
vibs.o:
virial.o:
volume.o: atoms.o iounit.o math.o sizes.o
warp.o:
xtalfit.o: atomid.o atoms.o bound.o boxes.o charge.o dipole.o energi.o files.o fracs.o inform.o iounit.o kvdws.o limits.o math.o molcul.o mpole.o polar.o potent.o sizes.o vdw.o vdwpot.o xtals.o
xtalmin.o: atoms.o boxes.o files.o inform.o iounit.o keys.o math.o scales.o sizes.o
xtals.o:
xyzatm.o: atoms.o inform.o iounit.o math.o sizes.o
xyzedit.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o limits.o math.o molcul.o ptable.o refer.o sizes.o titles.o units.o usage.o
xyzint.o: files.o iounit.o titles.o
xyzpdb.o: atomid.o atoms.o couple.o fields.o files.o inform.o molcul.o pdb.o resdue.o sequen.o sizes.o
xyzsybyl.o: files.o iounit.o sizes.o titles.o
zatom.o: angbnd.o atomid.o atoms.o bndstr.o fields.o iounit.o kangs.o katoms.o kbonds.o sizes.o zclose.o zcoord.o
zclose.o: sizes.o
zcoord.o: sizes.o
Last login: Thu Feb 11 21:22:49 on console
sleigh 1>ls
Applications/			ice2/
Applications (Parallels)/	info/
Caches/				intel/
Calibre Library/		jac/
Desktop/			jogamp/
Documents/			josh/
Downloads/			lectures/
Google Drive/			lee-ping/
Library/			letterhead/
MolecularTools/			manual.pdf
Movies/				mdvry-1.0/
Music/				meetings.txt
Pictures/			mike-solv/
Public/				network.txt
Sites/				nih-report/
amber-ff03/			nobel-2013/
amber10.ffparms/		nobel-legacy/
amber9/				ocean/
aminos/				openmm-master/
amoeba/				openmm-new/
apbs-1.3/			openmm-test/
apbs-1.4/			orca/
apbs-build/			papers/
apbs-test/			pdb2pqr/
autoconf/			penetration/
bachscore/			plugins/
backgrounds/			pmemd/
basketball/			ponder/
books/				posters/
boss/				process.run*
brooke/				proteins/
build_openmm/			proteins-old/
cb7/				psi4.0b5/
charmm34/			qchem/
chemspider/			reprints/
crambin/			salt/
cucu/				sampl5/
dhfr/				signature/
diala/				talks/
diis-eastman.txt		tests/
diis-sherrill.pdf		tinker/
dynamic_new/			tinker-6.2/
dynamic_omm/			tinker-git/
examples/			tinker-iso/
family/				tinker-lfmm/
ffe/				tinker-osrw/
ffe-win/			tinker-osrw-old/
ffe-win.tar.gz			tinker-osrw-older/
ffe.tar.gz			tinker-svn/
ffx/				tinker-svn.tar.gz
ffx-13apr/			tinker.tar.gz
finance/			tirex.tar.gz
fitpen/				torsion-fit/
forms/				usolve/
grc-2014/			water/
hbond/				water2/
hdac/				wca/
iamoeba/
sleigh 2>dir
total 8620280
drwxr-xr-x@  168 ponder  staff        5712 Feb 11 21:22 ./
drwxr-xr-x     6 root    admin         204 Feb  4 00:49 ../
-r--------     1 ponder  staff           7 Nov  5  2014 .CFUserTextEncoding
-rw-r--r--@    1 ponder  staff       32772 Feb  9 18:24 .DS_Store
drwxr-xr-x     3 ponder  staff         102 Sep 17  2014 .DownloadManager/
drwxr-xr-x     2 root    staff          68 Jun 25  2015 .InstallAnywhere/
drwxr-xr-x     3 ponder  staff         102 Oct  7 16:26 .RSA/
-rw-r--r--     1 ponder  staff           0 Sep 17  2014 .Rapp.history
drwx------     2 ponder  staff          68 Feb 11 21:22 .Trash/
-rw-------     1 ponder  staff          51 Apr  9  2015 .Xauthority
drwxr-xr-x     2 ponder  staff          68 Sep 17  2014 .Xcode/
drwxr-x--x     4 ponder  staff         136 Sep 17  2014 .adobe/
-rw-r--r--     1 ponder  staff         140 Sep 17  2014 .appletviewer
-rw-------     1 ponder  staff        2191 Sep 17  2014 .bash_history
drwxr-xr-x     3 ponder  staff         102 Apr  9  2015 .cache/
drwxr-xr-x     4 ponder  staff         136 Jan  2  2015 .config/
drwx------     3 ponder  staff         102 Sep 17  2014 .cups/
-rw-------     1 ponder  staff          58 Sep 17  2014 .cvspass
-rw-r--r--     1 ponder  staff          32 Sep 25 17:20 .deskmetrics
drwx------    24 ponder  staff         816 Feb 11 21:22 .dropbox/
drwxr-xr-x     4 ponder  staff         136 Oct 19 22:51 .dvdcss/
-rwxr-xr-x     1 ponder  staff         158 Jun 25  2015 .flexlmrc*
drwxr-xr-x     3 ponder  staff         102 Sep 17  2014 .fontconfig/
-rw-------     1 ponder  staff          43 Sep 17  2014 .gdb_history
-rw-------     1 ponder  staff       30823 Aug 29 20:07 .gnuplot_history
drwxr-xr-x     3 ponder  staff         102 Sep 17  2014 .groovy/
drwxr-xr-x     3 ponder  staff         102 Sep 17  2014 .idlerc/
-rw-rw-r--     1 ponder  staff         636 Aug 12  2015 .install4j
drwxr-xr-x     3 ponder  staff         102 Sep 17  2014 .install4j5/
drwxr-xr-x     4 ponder  staff         136 Aug 12  2015 .install4j6/
-rw-r--r--     1 ponder  staff        1443 Nov 13 10:45 .iqmol.log
drwxr-xr-x    18 ponder  staff         612 Sep 17  2014 .jbidwatcher/
drwxr-xr-x     3 ponder  staff         102 Sep 17  2014 .jmol/
-rw-------     1 ponder  staff          35 Aug 14 11:17 .lesshst
-rw-------     1 ponder  staff        1822 Sep 17  2014 .lldb-history
drwxr-xr-x     3 ponder  staff         102 Jan  2  2015 .local/
drwxr-xr-x     2 ponder  staff          68 Feb  8 08:36 .maplesoft/
drwxr-xr-x     4 ponder  staff         136 Jan 13 23:23 .matlab/
drwxr-xr-x     4 ponder  staff         136 Apr  9  2015 .mgltools/
drwxr-xr-x     3 ponder  staff         102 Sep 17  2014 .mobac/
drwxr-xr-x     5 ponder  staff         170 Nov 11 13:49 .oracle_jre_usage/
drwxr-xr-x     4 ponder  staff         136 Sep 17  2014 .parallels/
-rw-r--r--     1 ponder  staff          47 Sep 17  2014 .parallels_settings
-rw-------     1 ponder  staff        1024 Apr  9  2015 .rnd
drwx------     9 ponder  staff         306 Sep 17  2014 .sage/
-rw-------     1 ponder  staff         133 Sep 17  2014 .sh_history
drwx------     7 ponder  staff         238 Dec  9  2014 .ssh/
drwxr-xr-x     6 ponder  staff         204 Sep 17  2014 .subversion/
-rw-r--r--     1 ponder  staff         252 Sep 17  2014 .swfinfo
drwxr-xr-x     3 ponder  staff         102 Sep 22  2014 .sybyl/
-rw-r--r--@    1 ponder  staff        5010 Jan 10 16:33 .tcshrc
-rw-------     1 ponder  staff        8150 Oct 26 18:18 .viminfo
-rw-r--r--     1 root    staff         996 Apr  9  2015 .xinitrc
drwxr-xr-x    13 ponder  staff         442 Aug 27 21:40 Applications/
drwxr-xr-x@   12 ponder  staff         408 Feb 10 22:20 Applications (Parallels)/
drwxr-xr-x     2 ponder  staff          68 Jan  3 00:19 Caches/
drwxr-xr-x     5 ponder  staff         170 Jan 12 22:39 Calibre Library/
drwx------+   36 ponder  staff        1224 Feb 11 22:01 Desktop/
drwx------@   17 ponder  staff         578 Jan 26 02:21 Documents/
drwx------+    4 ponder  staff         136 Jul 17  2015 Downloads/
drwxr-xr-x@    3 ponder  staff         102 Jan  1  2013 Google Drive/
drwx------@   73 ponder  staff        2482 Feb  4 01:02 Library/
drwxr-xr-x     9 ponder  staff         306 Aug 25 19:42 MolecularTools/
drwx------+   11 ponder  staff         374 Jul 12  2015 Movies/
drwx------+   95 ponder  staff        3230 Jan 23 23:38 Music/
drwx------+    7 ponder  staff         238 Aug 21 14:53 Pictures/
drwxr-xr-x+    5 ponder  staff         170 Dec 22  2013 Public/
drwxr-xr-x+    3 ponder  staff         102 Nov 25  2010 Sites/
drwxr-xr-x    12 ponder  staff         408 Sep 17  2014 amber-ff03/
drwx------     7 ponder  staff         238 Sep 17  2014 amber10.ffparms/
drwxr-xr-x    15 ponder  staff         510 Sep 17  2014 amber9/
drwxr-xr-x    11 ponder  staff         374 Sep 17  2014 aminos/
drwxr-xr-x    51 ponder  staff        1734 Jan 16 18:11 amoeba/
drwxr-xr-x    35 ponder  staff        1190 Sep 17  2014 apbs-1.3/
drwxr-xr-x@   16 ponder  staff         544 Sep 17  2014 apbs-1.4/
drwxr-xr-x     9 ponder  staff         306 Sep 17  2014 apbs-build/
drwxr-xr-x    11 ponder  staff         374 Sep 17  2014 apbs-test/
drwxr-xr-x     6 ponder  staff         204 May 13  2014 autoconf/
drwxr-xr-x     8 ponder  staff         272 Feb 13  2014 bachscore/
drwxr-xr-x     9 ponder  staff         306 Feb 12  2014 backgrounds/
drwxr-xr-x    26 ponder  staff         884 Nov 21  2011 basketball/
drwxr-xr-x   111 ponder  staff        3774 Feb  4 12:44 books/
drwxr-xr-x    32 ponder  staff        1088 Aug 21  2009 boss/
drwxr-xr-x    14 ponder  staff         476 May 22  2015 brooke/
drwxr-xr-x   201 ponder  staff        6834 Apr 13  2015 build_openmm/
drwxr-xr-x     5 ponder  staff         170 Apr 28  2015 cb7/
drwxr-xr-x    20 ponder  staff         680 Mar 22  2014 charmm34/
drwxr-xr-x    29 ponder  staff         986 Jan 30  2014 chemspider/
drwxr-xr-x    48 ponder  staff        1632 Jul  2  2015 crambin/
drwxr-xr-x    64 ponder  staff        2176 Sep  5  2013 cucu/
drwxr-xr-x     8 ponder  staff         272 Dec  7  2014 dhfr/
drwxr-xr-x     5 ponder  staff         170 Oct 24  2011 diala/
-rw-r--r--     1 ponder  staff        6952 Nov 27  2014 diis-eastman.txt
-rw-r--r--@    1 ponder  staff       49357 Nov 27  2014 diis-sherrill.pdf
drwxr-xr-x  1159 ponder  staff       39406 Mar 19  2015 dynamic_new/
drwxr-xr-x  1190 ponder  staff       40460 May  1  2015 dynamic_omm/
drwxr-xr-x     6 ponder  staff         204 Mar 12  2015 examples/
drwxr-xr-x    98 ponder  staff        3332 Feb  4 11:47 family/
drwxr-xr-x     5 ponder  staff         170 Aug 14 13:24 ffe/
drwxr-xr-x     5 ponder  staff         170 Aug 17 17:45 ffe-win/
-rw-rw-rw-     1 ponder  staff   784177476 Aug 17 17:50 ffe-win.tar.gz
-rw-r--r--     1 ponder  staff  3350419282 Aug 15 00:16 ffe.tar.gz
drwxr-xr-x    16 ponder  staff         544 Aug 10  2010 ffx/
drwxr-xr-x    14 ponder  staff         476 Apr 19  2013 ffx-13apr/
drwxr-xr-x    10 ponder  staff         340 Jan 29  2014 finance/
drwxr-xr-x    10 ponder  staff         340 Sep  2  2014 fitpen/
drwxr-xr-x     7 ponder  staff         238 Feb  1 19:14 forms/
drwxr-xr-x    23 ponder  staff         782 Nov 13 10:27 grc-2014/
drwxr-xr-x    25 ponder  staff         850 Oct  2  2009 hbond/
drwxr-xr-x     6 ponder  staff         204 Aug  8  2013 hdac/
drwxr-xr-x    16 ponder  staff         544 Feb  8  2012 iamoeba/
drwxr-xr-x     5 ponder  staff         170 Mar  3  2013 ice2/
drwxr-xr-x   104 ponder  staff        3536 Feb  9 12:55 info/
drwxr-xr-x    10 ponder  staff         340 May  9  2015 intel/
drwxr-xr-x    52 ponder  staff        1768 Jan 27  2014 jac/
drwxr-xr-x     6 ponder  staff         204 Jun  3  2015 jogamp/
drwxr-xr-x    24 ponder  staff         816 Jun  7  2015 josh/
drwxr-xr-x     6 ponder  staff         204 Aug 24 23:56 lectures/
drwxr-xr-x    61 ponder  staff        2074 Jul 14  2014 lee-ping/
drwxr-xr-x    13 ponder  staff         442 Jul  2  2015 letterhead/
-rw-r--r--@    1 ponder  staff      434651 Feb  9 13:03 manual.pdf
drwxr-xr-x+    9 ponder  staff         306 Oct 24  2011 mdvry-1.0/
-rw-r--r--@    1 ponder  staff         750 Jan  8 12:06 meetings.txt
drwxr-xr-x     5 ponder  staff         170 Apr 28  2009 mike-solv/
-rw-r--r--     1 ponder  staff        7593 Feb  7 17:02 network.txt
drwxr-xr-x    15 ponder  staff         510 Feb 27  2014 nih-report/
drwxr-xr-x     4 ponder  staff         136 Oct 17  2013 nobel-2013/
drwxr-xr-x     7 ponder  staff         238 Oct 20 10:42 nobel-legacy/
drwxr-xr-x    22 ponder  staff         748 Jan 10  2015 ocean/
drwxr-xr-x@   20 ponder  staff         680 Mar 19  2015 openmm-master/
drwxr-xr-x@   20 ponder  staff         680 Mar 15  2015 openmm-new/
drwxr-xr-x    10 ponder  staff         340 Feb 11  2015 openmm-test/
drwxr-xr-x    10 ponder  staff         340 Apr 23  2008 orca/
drwxr-xr-x    65 ponder  staff        2210 Oct 22 12:59 papers/
drwxr-xr-x    32 ponder  staff        1088 Apr 30  2012 pdb2pqr/
drwxr-xr-x    16 ponder  staff         544 Oct 27  2014 penetration/
drwxr-xr-x     3 ponder  staff         102 Apr 21  2015 plugins/
drwxr-xr-x    14 ponder  staff         476 Mar  3  2011 pmemd/
drwxr-xr-x     2 ponder  staff          68 Jan  3 00:19 ponder/
drwxr-xr-x    22 ponder  staff         748 Oct  6 21:26 posters/
-rwxr-xr-x@    1 ponder  staff        1060 Oct 27 17:17 process.run*
drwxr-xr-x    18 ponder  staff         612 Aug 11  2010 proteins/
drwxr-xr-x    15 ponder  staff         510 Mar  8  2010 proteins-old/
drwxr-xr-x     8 ponder  staff         272 May 19  2014 psi4.0b5/
drwxr-xr-x    15 ponder  staff         510 Oct 22 13:48 qchem/
drwxr-xr-x   330 ponder  staff       11220 Feb  5 23:54 reprints/
drwxr-xr-x     8 ponder  staff         272 Oct 23  2014 salt/
drwxr-xr-x    13 ponder  staff         442 Feb  4 21:08 sampl5/
drwxr-xr-x    12 ponder  staff         408 Aug 22 20:44 signature/
drwxr-xr-x    58 ponder  staff        1972 Jan 12 14:47 talks/
drwxr-xr-x     6 ponder  staff         204 Mar 12  2015 tests/
drwxr-xr-x    27 ponder  staff         918 Jan 29 23:15 tinker/
drwxr-xr-x    18 ponder  staff         612 Jun  5  2015 tinker-6.2/
drwxr-xr-x    22 ponder  staff         748 Dec 26 00:00 tinker-git/
drwxr-xr-x    22 ponder  staff         748 Feb  8 13:41 tinker-iso/
drwxr-xr-x     8 ponder  staff         272 Apr 30  2015 tinker-lfmm/
drwxr-xr-x@    5 ponder  staff         170 Mar  4  2014 tinker-osrw/
drwxr-xr-x     5 ponder  staff         170 Mar  1  2014 tinker-osrw-old/
drwxr-xr-x   526 ponder  staff       17884 Feb 11  2014 tinker-osrw-older/
drwxr-xr-x    19 ponder  staff         646 Jun 23  2015 tinker-svn/
-rw-r--r--     1 ponder  staff    68560883 Jun 23  2015 tinker-svn.tar.gz
-rw-r--r--     1 ponder  staff   209240453 Feb  6 10:17 tinker.tar.gz
-rw-r--r--     1 ponder  staff      500996 Sep 29  2010 tirex.tar.gz
drwxr-xr-x    12 ponder  staff         408 Aug 24  2011 torsion-fit/
drwxr-xr-x    11 ponder  staff         374 Oct  9 00:02 usolve/
drwxr-xr-x    93 ponder  staff        3162 Aug 12  2015 water/
drwxr-xr-x    83 ponder  staff        2822 Dec  1  2011 water2/
drwxr-xr-x    15 ponder  staff         510 May 24  2008 wca/
sleigh 3>cd build_openmm
sleigh 4>ls
CMakeCache.txt
CMakeFiles/
CTestTestfile.cmake
DartConfiguration.tcl
DoxyfileC++*
DoxyfilePython*
HelloArgon*
HelloArgonInC*
HelloEthane*
HelloSodiumChloride*
HelloSodiumChlorideInC*
HelloWaterBox*
Makefile
TestCpuCustomGBForce*
TestCpuCustomManyParticleForce*
TestCpuCustomNonbondedForce*
TestCpuEwald*
TestCpuGBSAOBCForce*
TestCpuLangevinIntegrator*
TestCpuNeighborList*
TestCpuNonbondedForce*
TestCpuPeriodicTorsionForce*
TestCpuRBTorsionForce*
TestCpuSettle*
TestCudaAmoebaAngleForce*
TestCudaAmoebaBondForce*
TestCudaAmoebaGeneralizedKirkwoodForce*
TestCudaAmoebaInPlaneAngleForce*
TestCudaAmoebaMultipoleForce*
TestCudaAmoebaOutOfPlaneBendForce*
TestCudaAmoebaPiTorsionForce*
TestCudaAmoebaStretchBendForce*
TestCudaAmoebaTorsionTorsionForce*
TestCudaAmoebaVdwForce*
TestCudaAndersenThermostat*
TestCudaBrownianIntegrator*
TestCudaCMAPTorsionForce*
TestCudaCMMotionRemover*
TestCudaCheckpoints*
TestCudaCustomAngleForce*
TestCudaCustomBondForce*
TestCudaCustomCompoundBondForce*
TestCudaCustomExternalForce*
TestCudaCustomGBForce*
TestCudaCustomHbondForce*
TestCudaCustomIntegrator*
TestCudaCustomManyParticleForce*
TestCudaCustomNonbondedForce*
TestCudaCustomTorsionForce*
TestCudaDrudeForce*
TestCudaDrudeLangevinIntegrator*
TestCudaDrudeSCFIntegrator*
TestCudaEwald*
TestCudaGBSAOBCForce*
TestCudaHarmonicAngleForce*
TestCudaHarmonicBondForce*
TestCudaLangevinIntegrator*
TestCudaLocalEnergyMinimizer*
TestCudaMonteCarloAnisotropicBarostat*
TestCudaMonteCarloBarostat*
TestCudaMultipleForces*
TestCudaNonbondedForce*
TestCudaPeriodicTorsionForce*
TestCudaRBTorsionForce*
TestCudaRandom*
TestCudaRpmd*
TestCudaSettle*
TestCudaSort*
TestCudaVariableLangevinIntegrator*
TestCudaVariableVerletIntegrator*
TestCudaVerletIntegrator*
TestCudaVirtualSites*
TestCudaWcaDispersionForce*
TestFindExclusions*
TestFindMolecules*
TestParser*
TestReferenceAmoebaAngleForce*
TestReferenceAmoebaBondForce*
TestReferenceAmoebaGeneralizedKirkwoodForce*
TestReferenceAmoebaInPlaneAngleForce*
TestReferenceAmoebaMultipoleForce*
TestReferenceAmoebaOutOfPlaneBendForce*
TestReferenceAmoebaPiTorsionForce*
TestReferenceAmoebaStretchBendForce*
TestReferenceAmoebaTorsionTorsionForce*
TestReferenceAmoebaVdwForce*
TestReferenceAndersenThermostat*
TestReferenceBrownianIntegrator*
TestReferenceCMAPTorsionForce*
TestReferenceCMMotionRemover*
TestReferenceCheckpoints*
TestReferenceCustomAngleForce*
TestReferenceCustomBondForce*
TestReferenceCustomCompoundBondForce*
TestReferenceCustomExternalForce*
TestReferenceCustomGBForce*
TestReferenceCustomHbondForce*
TestReferenceCustomIntegrator*
TestReferenceCustomManyParticleForce*
TestReferenceCustomNonbondedForce*
TestReferenceCustomTorsionForce*
TestReferenceDrudeForce*
TestReferenceDrudeLangevinIntegrator*
TestReferenceDrudeSCFIntegrator*
TestReferenceEwald*
TestReferenceGBSAOBCForce*
TestReferenceGBVIForce*
TestReferenceHarmonicAngleForce*
TestReferenceHarmonicBondForce*
TestReferenceKineticEnergy*
TestReferenceLangevinIntegrator*
TestReferenceLocalEnergyMinimizer*
TestReferenceMonteCarloAnisotropicBarostat*
TestReferenceMonteCarloBarostat*
TestReferenceMonteCarloMembraneBarostat*
TestReferenceNeighborList*
TestReferenceNonbondedForce*
TestReferencePeriodicTorsionForce*
TestReferenceRBTorsionForce*
TestReferenceRandom*
TestReferenceRpmd*
TestReferenceSettle*
TestReferenceVariableLangevinIntegrator*
TestReferenceVariableVerletIntegrator*
TestReferenceVerletIntegrator*
TestReferenceVirtualSites*
TestReferenceWcaDispersionForce*
TestSerializationNode*
TestSerializeAmoebaAngleForce*
TestSerializeAmoebaBondForce*
TestSerializeAmoebaGeneralizedKirkwoodForce*
TestSerializeAmoebaInPlaneAngleForce*
TestSerializeAmoebaMultipoleForce*
TestSerializeAmoebaOutOfPlaneBendForce*
TestSerializeAmoebaPiTorsionForce*
TestSerializeAmoebaStretchBendForce*
TestSerializeAmoebaTorsionTorsionForce*
TestSerializeAmoebaVdwForce*
TestSerializeAmoebaWcaDispersionForce*
TestSerializeAndersenThermostat*
TestSerializeCMAPTorsion*
TestSerializeCMMotionRemover*
TestSerializeCustomAngleForce*
TestSerializeCustomBondForce*
TestSerializeCustomCompoundBondForce*
TestSerializeCustomExternalForce*
TestSerializeCustomGBForce*
TestSerializeCustomHbondForce*
TestSerializeCustomManyParticleForce*
TestSerializeCustomNonbondedForce*
TestSerializeCustomTorsionForce*
TestSerializeDrudeForce*
TestSerializeDrudeLangevinIntegrator*
TestSerializeGBSAOBCForce*
TestSerializeGBVIForce*
TestSerializeHarmonicAngleForce*
TestSerializeHarmonicBondForce*
TestSerializeIntegrator*
TestSerializeMonteCarloAnisotropicBarostat*
TestSerializeMonteCarloBarostat*
TestSerializeMonteCarloMembraneBarostat*
TestSerializeNonbondedForce*
TestSerializePeriodicTorsionForce*
TestSerializeRBTorsionForce*
TestSerializeState*
TestSerializeSystem*
TestSerializeTabulatedFunctions*
TestSplineFitter*
TestSystem*
TestVectorize*
Testing/
api-c++/
api-python/
build.log
cmake_install.cmake
docs-source/
examples/
install-python.log
install.log
install_manifest.txt
libOpenMM.dylib*
libOpenMMAmoeba.dylib*
libOpenMMAmoebaCUDA.dylib*
libOpenMMAmoebaReference.dylib*
libOpenMMCPU.dylib*
libOpenMMCUDA.dylib*
libOpenMMDrude.dylib*
libOpenMMDrudeCUDA.dylib*
libOpenMMDrudeReference.dylib*
libOpenMMRPMD.dylib*
libOpenMMRPMDCUDA.dylib*
libOpenMMRPMDReference.dylib*
platforms/
plugins/
python/
serialization/
sphinx-docs/
tests/
wrappers/
sleigh 5>cd ..
sleigh 6>ls
Applications/			ice2/
Applications (Parallels)/	info/
Caches/				intel/
Calibre Library/		jac/
Desktop/			jogamp/
Documents/			josh/
Downloads/			lectures/
Google Drive/			lee-ping/
Library/			letterhead/
MolecularTools/			manual.pdf
Movies/				mdvry-1.0/
Music/				meetings.txt
Pictures/			mike-solv/
Public/				network.txt
Sites/				nih-report/
amber-ff03/			nobel-2013/
amber10.ffparms/		nobel-legacy/
amber9/				ocean/
aminos/				openmm-master/
amoeba/				openmm-new/
apbs-1.3/			openmm-test/
apbs-1.4/			orca/
apbs-build/			papers/
apbs-test/			pdb2pqr/
autoconf/			penetration/
bachscore/			plugins/
backgrounds/			pmemd/
basketball/			ponder/
books/				posters/
boss/				process.run*
brooke/				proteins/
build_openmm/			proteins-old/
cb7/				psi4.0b5/
charmm34/			qchem/
chemspider/			reprints/
crambin/			salt/
cucu/				sampl5/
dhfr/				signature/
diala/				talks/
diis-eastman.txt		tests/
diis-sherrill.pdf		tinker/
dynamic_new/			tinker-6.2/
dynamic_omm/			tinker-git/
examples/			tinker-iso/
family/				tinker-lfmm/
ffe/				tinker-osrw/
ffe-win/			tinker-osrw-old/
ffe-win.tar.gz			tinker-osrw-older/
ffe.tar.gz			tinker-svn/
ffx/				tinker-svn.tar.gz
ffx-13apr/			tinker.tar.gz
finance/			tirex.tar.gz
fitpen/				torsion-fit/
forms/				usolve/
grc-2014/			water/
hbond/				water2/
hdac/				wca/
iamoeba/
sleigh 7>ssh dasher
Last login: Thu Feb 11 14:27:05 2016 from 99-7-216-113.lightspeed.stlsmo.sbcglobal.net
dasher 1>ls
Applications/			marie/
Desktop/			max/
Documents/			methanol-water/
Downloads/			mike-molecules/
Dropbox/			mike-tests/
Library/			mmcc/
MolecularTools/			motto.txt
Movies/				nanotube/
Music/				network.txt
Pictures/			nma-test/
Public/				nobel-2013/
Sites/				nucleic-acid/
Templates/			ocean/
Videos/				openmm-master/
aib/				openmm-test/
alamet/				openmpi-1.6.3/
amber-darden/			opep/
aminos/				other-3.pdf
amoeba/				papers/
amoeba-bio.out			penetration/
amoeba-bio.out2			perl/
amoeba-bio.out3			photos/
amoeba.out3			posters/
amoeba2.out			potential/
amoebabio.out			pressure/
amoebabio.out2			proteins/
amoebabio09-melys.prm		proteins-old/
archive/			psi4.0b5/
bachscore/			qmtinker/
backgrounds/			qsar/
basketball/			report.pdf
biophysj-106-a403-14.txt	reprints/
biophysj-96-a5-09.pdf		salt/
books/				sampl-info/
brooke/				sampl5/
build_openmm/			sas-pro/
clip/				scans/
clipart/			science-302-1364-03.pdf
clu/				sda7/
cluster/			sept/
crambin/			signature/
dasher.cer			sleuth/
dimer/				small-molecules/
dynamic_omm/			talks/
dynamic_omm-linux/		taylor-nucleosome/
embo/				tinker/
enkeph/				tinker-arc/
evfold/				tinker-bin/
exam-02/			tinker-extra/
ffe/				tinker-gpu.tar.gz
ffe.tar.gz*			tinker-iso/
ffx/				tinker-iso.tar.gz
forcefield/			tinker-lfmm/
forms/				tinker-mike/
gamess/				tinker-new/
gaussian/			tinker-new.tar.gz
graphite/			tinker-osrw/
grc-2014/			tinker-osrw-pren/
hbond/				tinker-quotes/
hints/				tinker-refs/
info/				tinker-src/
integrator/			tinker-svn.tar.gz
intel/				tinker.tar.gz
ionic-liquids/			tirex/
iwater/				to-be-filed/
jac/				torsion-fit/
jcc-37-494-16.pdf		trpcage/
jcc-8-1016-87.pdf		typer/
josh/				usage/
letterhead/			waterbox/
loos/				wca/
maple/				woodward/
dasher 2>cd dynamic_omm
dasher 3>ls
0README			etors2.f		nose.f
Makefile		etors2.o		nose.o
action.f		etors3.f		nspline.f
action.mod		etors3.o		nspline.o
action.o		etortor.f		nucleic.f
active.f		etortor.o		nucleo.f
active.o		etortor1.f		nucleo.mod
alchemy.f		etortor1.o		nucleo.o
align.f			etortor2.f		number.f
align.mod		etortor2.o		number.o
align.o			etortor3.f		numeral.f
amoebabio09.prm		etortor3.o		numeral.o
amoebapro13.prm		eurey.f			numgrad.f
analysis.f		eurey.o			numgrad.o
analysis.o		eurey1.f		ocvm.f
analyz.f		eurey1.o		ocvm.o
analyz.mod		eurey2.f		omega.f
analyz.o		eurey2.o		omega.mod
analyze.f		eurey3.f		omega.o
angang.f		eurey3.o		ommstuff.cpp
angang.mod		evcorr.f		ommstuff.o
angang.o		evcorr.o		opbend.f
angbnd.f		ewald.f			opbend.mod
angbnd.mod		ewald.mod		opbend.o
angbnd.o		ewald.o			opdist.f
angles.f		extra.f			opdist.mod
angles.o		extra.o			opdist.o
angpot.f		extra1.f		openend.f
angpot.mod		extra1.o		openend.o
angpot.o		extra2.f		openmm.f
angtor.f		extra2.o		openmm.mod
angtor.mod		extra3.f		openmm.o
angtor.o		extra3.o		openmp.f
anneal.f		faces.f			openmp.mod
archive.f		faces.mod		openmp.o
argue.f			faces.o			optimize.f
argue.mod		fatal.f			optirot.f
argue.o			fatal.o			optrigid.f
ascii.f			fft.f			optsave.f
ascii.mod		fft.mod			optsave.o
ascii.o			fft.o			orbital.f
atmlst.f		fft3d.f			orbital.o
atmlst.mod		fft3d.o			orbits.f
atmlst.o		fftpack.f		orbits.mod
atomid.f		fftpack.o		orbits.o
atomid.mod		field.f			orient.f
atomid.o		field.o			orient.o
atoms.f			fields.f		orthog.f
atoms.mod		fields.mod		orthog.o
atoms.o			fields.o		output.f
attach.f		files.f			output.mod
attach.o		files.mod		output.o
bar.f			files.o			overlap.f
basefile.f		final.f			overlap.o
basefile.o		final.o			params.f
bath.f			flatten.f		params.mod
bath.mod		flatten.o		params.o
bath.o			fracs.f			path.f
beeman.f		fracs.mod		paths.f
beeman.o		fracs.o			paths.mod
bicubic.f		freeunit.f		paths.o
bicubic.o		freeunit.o		pbstuf.f
bitor.f			freeze.f		pbstuf.mod
bitor.mod		freeze.mod		pbstuf.o
bitor.o			freeze.o		pdb.f
bitors.f		gda.f			pdb.mod
bitors.o		geometry.f		pdb.o
bndpot.f		geometry.o		pdbxyz.f
bndpot.mod		getint.f		phipsi.f
bndpot.o		getint.o		phipsi.mod
bndstr.f		getkey.f		phipsi.o
bndstr.mod		getkey.o		picalc.f
bndstr.o		getmol.f		picalc.o
bonds.f			getmol.o		piorbs.f
bonds.o			getmol2.f		piorbs.mod
born.f			getmol2.o		piorbs.o
born.o			getnumb.f		pistuf.f
bound.f			getnumb.o		pistuf.mod
bound.mod		getpdb.f		pistuf.o
bound.o			getpdb.o		pitors.f
bounds.f		getprm.f		pitors.mod
bounds.o		getprm.o		pitors.o
boxes.f			getref.f		pme.f
boxes.mod		getref.o		pme.mod
boxes.o			getstring.f		pme.o
bussi.f			getstring.o		pmestuff.f
bussi.o			gettext.f		pmestuff.o
calendar.f		gettext.o		pmpb.f
calendar.o		getword.f		pmpb.o
cell.f			getword.o		polar.f
cell.mod		getxyz.f		polar.mod
cell.o			getxyz.o		polar.o
center.f		ghmcstep.f		polarize.f
center.o		ghmcstep.o		poledit.f
charge.f		gkstuf.f		polgrp.f
charge.mod		gkstuf.mod		polgrp.mod
charge.o		gkstuf.o		polgrp.o
chgpot.f		gradient.f		polpot.f
chgpot.mod		gradient.o		polpot.mod
chgpot.o		gradrgd.f		polpot.o
chkpole.f		gradrgd.o		polymer.f
chkpole.o		gradrot.f		polymer.o
chkring.f		gradrot.o		potent.f
chkring.o		group.f			potent.mod
chkxyz.f		group.mod		potent.o
chkxyz.o		group.o			potential.f
cholesky.f		groups.f		potfit.f
cholesky.o		groups.o		potfit.mod
chrono.f		grpline.f		potfit.o
chrono.mod		grpline.o		precis.f
chrono.o		gyrate.f		precis.mod
chunks.f		gyrate.o		precis.o
chunks.mod		hescut.f		precise.f
chunks.o		hescut.mod		precise.o
clock.f			hescut.o		pressure.f
clock.o			hessian.f		pressure.o
cluster.f		hessian.o		prmedit.f
cluster.o		hessn.f			prmkey.f
column.f		hessn.mod		prmkey.o
column.o		hessn.o			promo.f
command.f		hessrgd.f		promo.o
command.o		hessrgd.o		protein.f
connect.f		hessrot.f		prtdyn.f
connect.o		hessrot.o		prtdyn.o
connolly.f		hpmf.f			prterr.f
connolly.o		hpmf.mod		prterr.o
control.f		hpmf.o			prtint.f
control.o		hybrid.f		prtint.o
correlate.f		hybrid.o		prtmol2.f
couple.f		image.f			prtmol2.o
couple.mod		image.o			prtpdb.f
couple.o		impose.f		prtpdb.o
crambin2.key		impose.o		prtprm.f
crambin2.xyz		improp.f		prtprm.o
crambin3.arc-750m	improp.mod		prtseq.f
crambin3.dyn-orig	improp.o		prtseq.o
crambin3.key		imptor.f		prtxyz.f
crambin3.log-750m	imptor.mod		prtxyz.o
crambin3.log-970	imptor.o		pss.f
crambin3.run*		induce.f		pssrigid.f
crambin3.xyz		induce.o		pssrot.f
crystal.f		inertia.f		ptable.f
cspline.f		inertia.o		ptable.mod
cspline.o		inform.f		ptable.o
cutoffs.f		inform.mod		qmstuf.f
cutoffs.o		inform.o		qmstuf.mod
deflate.f		initatom.f		qmstuf.o
deflate.o		initatom.o		qrfact.f
delete.f		initial.f		qrfact.o
delete.o		initial.o		quatfit.f
deriv.f			initprm.f		quatfit.o
deriv.mod		initprm.o		radial.f
deriv.o			initres.f		random.f
dhfr2.arc		initres.o		random.o
dhfr2.dyn		initrot.f		rattle.f
dhfr2.dyn-orig		initrot.o		rattle.o
dhfr2.key		insert.f		readdyn.f
dhfr2.log-750m		insert.o		readdyn.o
dhfr2.log-970		intedit.f		readgau.f
dhfr2.run*		inter.f			readgau.o
dhfr2.seq		inter.mod		readint.f
dhfr2.xyz		inter.o			readint.o
diagq.f			intxyz.f		readmol.f
diagq.o			invbeta.f		readmol.o
diffeq.f		invbeta.o		readmol2.f
diffeq.o		invert.f		readmol2.o
diffuse.f		invert.o		readpdb.f
dipole.f		iounit.f		readpdb.o
dipole.mod		iounit.mod		readprm.f
dipole.o		iounit.o		readprm.o
disgeo.f		jacobi.f		readseq.f
disgeo.mod		jacobi.o		readseq.o
disgeo.o		kanang.f		readxyz.f
distgeom.f		kanang.mod		readxyz.o
dma.f			kanang.o		refer.f
dma.mod			kangang.f		refer.mod
dma.o			kangang.o		refer.o
document.f		kangle.f		replica.f
domega.f		kangle.o		replica.o
domega.mod		kangs.f			resdue.f
domega.o		kangs.mod		resdue.mod
dynamic.f		kangs.o			resdue.o
dynamic_omm.f		kangtor.f		respa.f
dynamic_omm.f-6.3.3	kangtor.o		respa.o
dynamic_omm.o		kantor.f		restrn.f
dynamic_omm.x*		kantor.mod		restrn.mod
eangang.f		kantor.o		restrn.o
eangang.o		katom.f			rgddyn.f
eangang1.f		katom.o			rgddyn.mod
eangang1.o		katoms.f		rgddyn.o
eangang2.f		katoms.mod		rgdstep.f
eangang2.o		katoms.o		rgdstep.o
eangang3.f		kbond.f			rigid.f
eangang3.o		kbond.o			rigid.mod
eangle.f		kbonds.f		rigid.o
eangle.o		kbonds.mod		ring.f
eangle1.f		kbonds.o		ring.mod
eangle1.o		kcharge.f		ring.o
eangle2.f		kcharge.o		rings.f
eangle2.o		kchrge.f		rings.o
eangle3.f		kchrge.mod		rmsfit.f
eangle3.o		kchrge.o		rmsfit.o
eangtor.f		kdipol.f		rotbnd.f
eangtor.o		kdipol.mod		rotbnd.mod
eangtor1.f		kdipol.o		rotbnd.o
eangtor1.o		kdipole.f		rotlist.f
eangtor2.f		kdipole.o		rotlist.o
eangtor2.o		kewald.f		rotpole.f
eangtor3.f		kewald.o		rotpole.o
eangtor3.o		kextra.f		rxnfld.f
ebond.f			kextra.o		rxnfld.mod
ebond.o			keys.f			rxnfld.o
ebond1.f		keys.mod		rxnpot.f
ebond1.o		keys.o			rxnpot.mod
ebond2.f		kgeom.f			rxnpot.o
ebond2.o		kgeom.o			saddle.f
ebond3.f		khbond.f		scales.f
ebond3.o		khbond.mod		scales.mod
ebuck.f			khbond.o		scales.o
ebuck.o			kimprop.f		scan.f
ebuck1.f		kimprop.o		sdstep.f
ebuck1.o		kimptor.f		sdstep.o
ebuck2.f		kimptor.o		search.f
ebuck2.o		kinetic.f		search.o
ebuck3.f		kinetic.o		sequen.f
ebuck3.o		kiprop.f		sequen.mod
echarge.f		kiprop.mod		sequen.o
echarge.o		kiprop.o		server.f
echarge1.f		kitors.f		server.o
echarge1.o		kitors.mod		shakeup.f
echarge2.f		kitors.o		shakeup.o
echarge2.o		kmetal.f		shunt.f
echarge3.f		kmetal.o		shunt.mod
echarge3.o		kmpole.f		shunt.o
echgdpl.f		kmpole.o		sigmoid.f
echgdpl.o		kmulti.f		sigmoid.o
echgdpl1.f		kmulti.mod		sizes.f
echgdpl1.o		kmulti.o		sizes.mod
echgdpl2.f		kopbend.f		sizes.o
echgdpl2.o		kopbend.o		sktstuff.f
echgdpl3.f		kopbnd.f		sktstuff.o
echgdpl3.o		kopbnd.mod		sniffer.f
edipole.f		kopbnd.o		socket.f
edipole.o		kopdist.f		socket.mod
edipole1.f		kopdist.o		socket.o
edipole1.o		kopdst.f		solute.f
edipole2.f		kopdst.mod		solute.mod
edipole2.o		kopdst.o		solute.o
edipole3.f		korbit.f		sort.f
edipole3.o		korbit.o		sort.o
egauss.f		korbs.f			spacefill.f
egauss.o		korbs.mod		spectrum.f
egauss1.f		korbs.o			square.f
egauss1.o		kpitor.f		square.o
egauss2.f		kpitor.mod		stodyn.f
egauss2.o		kpitor.o		stodyn.mod
egauss3.f		kpitors.f		stodyn.o
egauss3.o		kpitors.o		strbnd.f
egeom.f			kpolar.f		strbnd.mod
egeom.o			kpolar.o		strbnd.o
egeom1.f		kpolr.f			strtor.f
egeom1.o		kpolr.mod		strtor.mod
egeom2.f		kpolr.o			strtor.o
egeom2.o		ksolv.f			suffix.f
egeom3.f		ksolv.o			suffix.o
egeom3.o		kstbnd.f		superpose.f
ehal.f			kstbnd.mod		surface.f
ehal.o			kstbnd.o		surface.o
ehal1.f			kstrbnd.f		surfatom.f
ehal1.o			kstrbnd.o		surfatom.o
ehal2.f			kstrtor.f		switch.f
ehal2.o			kstrtor.o		switch.o
ehal3.f			ksttor.f		sybylxyz.f
ehal3.o			ksttor.mod		syntrn.f
eimprop.f		ksttor.o		syntrn.mod
eimprop.o		ktors.f			syntrn.o
eimprop1.f		ktors.o			tarray.f
eimprop1.o		ktorsn.f		tarray.mod
eimprop2.f		ktorsn.mod		tarray.o
eimprop2.o		ktorsn.o		temper.f
eimprop3.f		ktortor.f		temper.o
eimprop3.o		ktortor.o		testgrad.f
eimptor.f		ktrtor.f		testhess.f
eimptor.o		ktrtor.mod		testpair.f
eimptor1.f		ktrtor.o		testpol.f
eimptor1.o		kurey.f			testrot.f
eimptor2.f		kurey.o			timer.f
eimptor2.o		kurybr.f		timerot.f
eimptor3.f		kurybr.mod		titles.f
eimptor3.o		kurybr.o		titles.mod
elj.f			kvdw.f			titles.o
elj.o			kvdw.o			tncg.f
elj1.f			kvdwpr.f		tncg.o
elj1.o			kvdwpr.mod		torphase.f
elj2.f			kvdwpr.o		torphase.o
elj2.o			kvdws.f			torpot.f
elj3.f			kvdws.mod		torpot.mod
elj3.o			kvdws.o			torpot.o
embed.f			lattice.f		torque.f
embed.o			lattice.o		torque.o
emetal.f		lbfgs.f			tors.f
emetal.o		lbfgs.o			tors.mod
emetal1.f		libfftw3.a		tors.o
emetal1.o		libfftw3_threads.a	torsfit.f
emetal2.f		libtinker.a		torsions.f
emetal2.o		light.f			torsions.o
emetal3.f		light.mod		tortor.f
emetal3.o		light.o			tortor.mod
emm3hb.f		lights.f		tortor.o
emm3hb.o		lights.o		tree.f
emm3hb1.f		limits.f		tree.mod
emm3hb1.o		limits.mod		tree.o
emm3hb2.f		limits.o		trimtext.f
emm3hb2.o		linmin.f		trimtext.o
emm3hb3.f		linmin.mod		unitcell.f
emm3hb3.o		linmin.o		unitcell.o
empole.f		makeint.f		units.f
empole.o		makeint.o		units.mod
empole1.f		makeref.f		units.o
empole1.o		makeref.o		uprior.f
empole2.f		makexyz.f		uprior.mod
empole2.o		makexyz.o		uprior.o
empole3.f		math.f			urey.f
empole3.o		math.mod		urey.mod
energi.f		math.o			urey.o
energi.mod		maxwell.f		urypot.f
energi.o		maxwell.o		urypot.mod
energy.f		mdinit.f		urypot.o
energy.o		mdinit.o		usage.f
eopbend.f		mdrest.f		usage.mod
eopbend.o		mdrest.o		usage.o
eopbend1.f		mdsave.f		usolve.f
eopbend1.o		mdsave.o		usolve.mod
eopbend2.f		mdstat.f		usolve.o
eopbend2.o		mdstat.o		valence.f
eopbend3.f		mdstuf.f		valfit.f
eopbend3.o		mdstuf.mod		valfit.mod
eopdist.f		mdstuf.o		valfit.o
eopdist.o		mechanic.f		vdw.f
eopdist1.f		mechanic.o		vdw.mod
eopdist1.o		merck.f			vdw.o
eopdist2.f		merck.mod		vdwpot.f
eopdist2.o		merck.o			vdwpot.mod
eopdist3.f		merge.f			vdwpot.o
eopdist3.o		merge.o			verlet.f
epitors.f		minima.f		verlet.o
epitors.o		minima.mod		version.f
epitors1.f		minima.o		version.o
epitors1.o		minimize.f		vibbig.f
epitors2.f		minirot.f		vibrate.f
epitors2.o		minrigid.f		vibrot.f
epitors3.f		molcul.f		vibs.f
epitors3.o		molcul.mod		vibs.mod
erf.f			molcul.o		vibs.o
erf.o			moldyn.f		virial.f
erxnfld.f		moldyn.mod		virial.mod
erxnfld.o		moldyn.o		virial.o
erxnfld1.f		molecule.f		volume.f
erxnfld1.o		molecule.o		volume.o
erxnfld2.f		molxyz.f		warp.f
erxnfld2.o		moment.f		warp.mod
erxnfld3.f		moment.mod		warp.o
erxnfld3.o		moment.o		water03.prm
esolv.f			moments.f		watermin.key
esolv.o			moments.o		watermin.xyz
esolv1.f		monte.f			watersmall.dyn-orig
esolv1.o		mplpot.f		watersmall.key
esolv2.f		mplpot.mod		watersmall.log-750m
esolv2.o		mplpot.o		watersmall.log-970
esolv3.f		mpole.f			watersmall.run*
esolv3.o		mpole.mod		watersmall.xyz
estrbnd.f		mpole.o			xtalfit.f
estrbnd.o		mutant.f		xtalmin.f
estrbnd1.f		mutant.mod		xtals.f
estrbnd1.o		mutant.o		xtals.mod
estrbnd2.f		mutate.f		xtals.o
estrbnd2.o		mutate.o		xyzatm.f
estrbnd3.f		nblist.f		xyzatm.o
estrbnd3.o		nblist.o		xyzedit.f
estrtor.f		neigh.f			xyzint.f
estrtor.o		neigh.mod		xyzpdb.f
estrtor1.f		neigh.o			xyzsybyl.f
estrtor1.o		newton.f		zatom.f
estrtor2.f		newtrot.f		zatom.o
estrtor2.o		nextarg.f		zclose.f
estrtor3.f		nextarg.o		zclose.mod
estrtor3.o		nexttext.f		zclose.o
etors.f			nexttext.o		zcoord.f
etors.o			nonpol.f		zcoord.mod
etors1.f		nonpol.mod		zcoord.o
etors1.o		nonpol.o
dasher 4>ls *.make
ls: No match.
dasher 5>ls *.make
ls: No match.
dasher 6>dir Makefile
-rw-r--r--  1 ponder  staff  53643 Mar 16  2015 Makefile
dasher 7>cat Makefile
##
######################################################################
##                                                                  ##
##  Makefile for Building an OpenMM-Enabled TINKER DYNAMIC Program  ##
##                                                                  ##
######################################################################
##
##  Invocation Options:
##
##   1. make all              Build TINKER dynamic executable
##   2. make install          Move the executable to BINDIR
##   3. make create_links     Create soft link in LINKDIR
##   4. make remove_links     Remove soft link from LINKDIR
##   6. make listing          Concatenate source to tinker.txt
##   5. make clean            Delete objects and executable
##
##  Modified from the standard TINKER Makefile in March 2015 by
##  Jay Ponder by build an OpenMM-accelerated dynamic executable.
##  Requires TINKER 7.1.2 or later, and assumes an AMOEBA-capable
##  version of OpenMM is already built and present
##
##  Original version of this file is due to Peter Happersberger
##  and Jochen Buehler of the University of Konstanz, January 1998.
##  Modifications by Reece Hart & Jay Ponder, Washington University
##
##  Building TINKER from source requires libraries from the FFTW
##  Fast Fourier Transform package and (optionally) APBS Poisson-
##  Boltzmann solver. These libraries are assumed to already be
##  present, and are not built by this Makefile
##
###################################################################

###################################################################
##  Master Directory Locations; Change as Needed for Local Site  ##
###################################################################

##  TINKERDIR    TINKER Distribution Directory
##  BINDIR       Hard Copies of TINKER Executables
##  LINKDIR      Linked Copies of TINKER Executables

TINKERDIR = $(HOME)/tinker
BINDIR = $(TINKERDIR)/bin
LINKDIR = /usr/local/bin
OPENMMDIR = /usr/local/openmm

OPENMM_LIB_DIR = -L$(OPENMMDIR)/lib
OPENMM_INCLUDE_DIR = $(OPENMMDIR)/include
OPENMM_LIBS = -lOpenMM -lOpenMMAmoeba

####################################################################
##  Known Machine Types; Uncomment One of the Following Sections  ##
##  May Need Editing to Match Your Desired OS & Compiler Version  ##
####################################################################

##  Machine:  Generic Linux or Macintosh
##  CPU Type: Intel x86 Compatible
##  Compiler: GNU gfortran

#CC = gcc
#CXX = g++
#F77 = gfortran
#LIBS = -lstdc++
#CFLAGS = -c
#F77FLAGS = -c
#OPTFLAGS = -O3 -ffast-math
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-libgcc

##  Machine:  Generic Linux or Macintosh
##  CPU Type: Intel x86 Compatible
##  Compiler: GNU gfortran
##  Parallel: OpenMP

#CC = gcc
#CXX = g++
#F77 = gfortran
#LIBS = -lstdc++ -L. -lfftw3_threads -lfftw3
#CFLAGS = -c
#F77FLAGS = -c
#OPTFLAGS = -O3 -ffast-math -fopenmp
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-libgcc

##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible (also AMD)
##  Oper Sys: Fedora Core
##  Compiler: Intel Fortran for Linux 14.0

#CC = icc
#CXX = icc
#F77 = ifort
#LIBS = -lstdc++
#CFLAGS = -c -xHost
#F77FLAGS = -c -xHost
#OPTFLAGS = -O3 -no-ipo -no-prec-div -recursive
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-libgcc -static-intel

##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible (also AMD)
##  Oper Sys: Fedora Core
##  Compiler: Intel Fortran for Linux 14.0
##  Parallel: OpenMP

#CC = icc
#CXX = icc
#F77 = ifort
#LIBS = -L. -lfftw3_threads -lfftw3
#CFLAGS = -c -xHost
#F77FLAGS = -c -xHost
#OPTFLAGS = -O3 -no-ipo -no-prec-div -recursive -openmp
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-libgcc -static-intel

##  Machine:  Macintosh
##  CPU Type: Intel Xeon
##  Oper Sys: OS X 10.10 (Yosemite)
##  Compiler: Intel Fortran for Mac 14.0

#CC = icc
#CXX = icc
#F77 = ifort
#LIBS = -lstdc++
#CFLAGS = -c -axSSSE3
#F77FLAGS = -c -axSSSE3
#OPTFLAGS = -O3 -no-ipo -no-prec-div -mdynamic-no-pic
#LIBFLAGS = -crusv
#RANLIB = ranlib -c
#LINKFLAGS = $(OPTFLAGS) -static-intel -mmacosx-version-min=10.6 \
#            -Wl,-stack_size,0x10000000

##  Machine:  Macintosh
##  CPU Type: Intel Xeon
##  Oper Sys: OS X 10.10 (Yosemite)
##  Compiler: Intel Fortran for Mac 14.0
##  Parallel: OpenMP

CC = icc
CXX = icc
F77 = ifort
LIBS = -lstdc++ -L. -lfftw3_threads -lfftw3
CFLAGS = -c -axSSSE3
F77FLAGS = -c -axSSSE3
OPTFLAGS = -O3 -no-ipo -no-prec-div -mdynamic-no-pic -openmp
LIBFLAGS = -crusv
RANLIB = ranlib -c
LINKFLAGS = $(OPTFLAGS) -static-intel -mmacosx-version-min=10.6 \
            -Wl,-stack_size,0x10000000

#################################################################
##  Should not be Necessary to Change Things Below this Point  ##
#################################################################

OBJS = action.o \
       active.o \
       align.o \
       analysis.o \
       analyz.o \
       angang.o \
       angbnd.o \
       angles.o \
       angpot.o \
       angtor.o \
       argue.o \
       ascii.o \
       atmlst.o \
       atomid.o \
       atoms.o \
       attach.o \
       basefile.o \
       bath.o \
       beeman.o \
       bicubic.o \
       bitor.o \
       bitors.o \
       bndpot.o \
       bndstr.o \
       bonds.o \
       born.o \
       bound.o \
       bounds.o \
       boxes.o \
       bussi.o \
       calendar.o \
       cell.o \
       center.o \
       charge.o \
       chgpot.o \
       chkpole.o \
       chkring.o \
       chkxyz.o \
       cholesky.o \
       chrono.o \
       chunks.o \
       clock.o \
       cluster.o \
       column.o \
       command.o \
       connect.o \
       connolly.o \
       control.o \
       couple.o \
       cspline.o \
       cutoffs.o \
       deflate.o \
       delete.o \
       deriv.o \
       diagq.o \
       diffeq.o \
       dipole.o \
       disgeo.o \
       dma.o \
       domega.o \
       dynamic_omm.o \
       eangang.o \
       eangang1.o \
       eangang2.o \
       eangang3.o \
       eangle.o \
       eangle1.o \
       eangle2.o \
       eangle3.o \
       eangtor.o \
       eangtor1.o \
       eangtor2.o \
       eangtor3.o \
       ebond.o \
       ebond1.o \
       ebond2.o \
       ebond3.o \
       ebuck.o \
       ebuck1.o \
       ebuck2.o \
       ebuck3.o \
       echarge.o \
       echarge1.o \
       echarge2.o \
       echarge3.o \
       echgdpl.o \
       echgdpl1.o \
       echgdpl2.o \
       echgdpl3.o \
       edipole.o \
       edipole1.o \
       edipole2.o \
       edipole3.o \
       egauss.o \
       egauss1.o \
       egauss2.o \
       egauss3.o \
       egeom.o \
       egeom1.o \
       egeom2.o \
       egeom3.o \
       ehal.o \
       ehal1.o \
       ehal2.o \
       ehal3.o \
       eimprop.o \
       eimprop1.o \
       eimprop2.o \
       eimprop3.o \
       eimptor.o \
       eimptor1.o \
       eimptor2.o \
       eimptor3.o \
       elj.o \
       elj1.o \
       elj2.o \
       elj3.o \
       embed.o \
       emetal.o \
       emetal1.o \
       emetal2.o \
       emetal3.o \
       emm3hb.o \
       emm3hb1.o \
       emm3hb2.o \
       emm3hb3.o \
       empole.o \
       empole1.o \
       empole2.o \
       empole3.o \
       energi.o \
       energy.o \
       eopbend.o \
       eopbend1.o \
       eopbend2.o \
       eopbend3.o \
       eopdist.o \
       eopdist1.o \
       eopdist2.o \
       eopdist3.o \
       epitors.o \
       epitors1.o \
       epitors2.o \
       epitors3.o \
       epolar.o \
       epolar1.o \
       epolar2.o \
       epolar3.o \
       erf.o \
       erxnfld.o \
       erxnfld1.o \
       erxnfld2.o \
       erxnfld3.o \
       esolv.o \
       esolv1.o \
       esolv2.o \
       esolv3.o \
       estrbnd.o \
       estrbnd1.o \
       estrbnd2.o \
       estrbnd3.o \
       estrtor.o \
       estrtor1.o \
       estrtor2.o \
       estrtor3.o \
       etors.o \
       etors1.o \
       etors2.o \
       etors3.o \
       etortor.o \
       etortor1.o \
       etortor2.o \
       etortor3.o \
       eurey.o \
       eurey1.o \
       eurey2.o \
       eurey3.o \
       evcorr.o \
       ewald.o \
       extra.o \
       extra1.o \
       extra2.o \
       extra3.o \
       faces.o \
       fatal.o \
       fft.o \
       fft3d.o \
       fftpack.o \
       field.o \
       fields.o \
       files.o \
       final.o \
       flatten.o \
       fracs.o \
       freeunit.o \
       freeze.o \
       geometry.o \
       getint.o \
       getkey.o \
       getmol.o \
       getmol2.o \
       getnumb.o \
       getpdb.o \
       getprm.o \
       getref.o \
       getstring.o \
       gettext.o \
       getword.o \
       getxyz.o \
       ghmcstep.o \
       gkstuf.o \
       gradient.o \
       gradrgd.o \
       gradrot.o \
       group.o \
       groups.o \
       grpline.o \
       gyrate.o \
       hescut.o \
       hessian.o \
       hessn.o \
       hessrgd.o \
       hessrot.o \
       hpmf.o \
       hybrid.o \
       image.o \
       impose.o \
       improp.o \
       imptor.o \
       induce.o \
       inertia.o \
       inform.o \
       initatom.o \
       initial.o \
       initprm.o \
       initres.o \
       initrot.o \
       insert.o \
       inter.o \
       invbeta.o \
       invert.o \
       iounit.o \
       jacobi.o \
       kanang.o \
       kangang.o \
       kangle.o \
       kangs.o \
       kangtor.o \
       kantor.o \
       katom.o \
       katoms.o \
       kbond.o \
       kbonds.o \
       kcharge.o \
       kchrge.o \
       kdipol.o \
       kdipole.o \
       kewald.o \
       kextra.o \
       keys.o \
       kgeom.o \
       khbond.o \
       kimprop.o \
       kimptor.o \
       kinetic.o \
       kiprop.o \
       kitors.o \
       kmetal.o \
       kmpole.o \
       kmulti.o \
       kopbend.o \
       kopbnd.o \
       kopdist.o \
       kopdst.o \
       korbit.o \
       korbs.o \
       kpitor.o \
       kpitors.o \
       kpolar.o \
       kpolr.o \
       ksolv.o \
       kstbnd.o \
       kstrbnd.o \
       kstrtor.o \
       ksttor.o \
       ktors.o \
       ktorsn.o \
       ktortor.o \
       ktrtor.o \
       kurey.o \
       kurybr.o \
       kvdw.o \
       kvdwpr.o \
       kvdws.o \
       lattice.o \
       lbfgs.o \
       light.o \
       lights.o \
       limits.o \
       linmin.o \
       makeint.o \
       makeref.o \
       makexyz.o \
       math.o \
       maxwell.o \
       mdinit.o \
       mdrest.o \
       mdsave.o \
       mdstat.o \
       mdstuf.o \
       mechanic.o \
       merck.o \
       merge.o \
       minima.o \
       molcul.o \
       moldyn.o \
       molecule.o \
       moment.o \
       moments.o \
       mplpot.o \
       mpole.o \
       mrecip.o \
       mutant.o \
       mutate.o \
       nblist.o \
       neigh.o \
       nextarg.o \
       nexttext.o \
       nonpol.o \
       nose.o \
       nspline.o \
       nucleo.o \
       number.o \
       numeral.o \
       numgrad.o \
       ocvm.o \
       omega.o \
       opbend.o \
       opdist.o \
       openend.o \
       openmm.o \
       openmp.o \
       optsave.o \
       orbital.o \
       orbits.o \
       orient.o \
       orthog.o \
       output.o \
       overlap.o \
       params.o \
       paths.o \
       pbstuf.o \
       pdb.o \
       phipsi.o \
       picalc.o \
       piorbs.o \
       pistuf.o \
       pitors.o \
       pme.o \
       pmestuff.o \
       pmpb.o \
       polar.o \
       polgrp.o \
       polpot.o \
       polymer.o \
       potent.o \
       potfit.o \
       precis.o \
       precise.o \
       pressure.o \
       prmkey.o \
       promo.o \
       prtdyn.o \
       prterr.o \
       prtint.o \
       prtmol2.o \
       prtpdb.o \
       prtprm.o \
       prtseq.o \
       prtxyz.o \
       ptable.o \
       qmstuf.o \
       qrfact.o \
       quatfit.o \
       random.o \
       rattle.o \
       readdyn.o \
       readgau.o \
       readint.o \
       readmol.o \
       readmol2.o \
       readpdb.o \
       readprm.o \
       readseq.o \
       readxyz.o \
       refer.o \
       replica.o \
       resdue.o \
       respa.o \
       restrn.o \
       rgddyn.o \
       rgdstep.o \
       rigid.o \
       ring.o \
       rings.o \
       rmsfit.o \
       rotbnd.o \
       rotlist.o \
       rotpole.o \
       rxnfld.o \
       rxnpot.o \
       scales.o \
       sdstep.o \
       search.o \
       sequen.o \
       server.o \
       shakeup.o \
       shunt.o \
       sigmoid.o \
       sizes.o \
       sktstuff.o \
       socket.o \
       solute.o \
       sort.o \
       square.o \
       stodyn.o \
       strbnd.o \
       strtor.o \
       suffix.o \
       surface.o \
       surfatom.o \
       switch.o \
       syntrn.o \
       tarray.o \
       temper.o \
       titles.o \
       tncg.o \
       torphase.o \
       torpot.o \
       torque.o \
       tors.o \
       torsions.o \
       tortor.o \
       tree.o \
       trimtext.o \
       unitcell.o \
       units.o \
       uprior.o \
       urey.o \
       urypot.o \
       usage.o \
       usolve.o \
       valfit.o \
       vdw.o \
       vdwpot.o \
       verlet.o \
       version.o \
       vibs.o \
       virial.o \
       volume.o \
       warp.o \
       xtals.o \
       xyzatm.o \
       zatom.o \
       zclose.o \
       zcoord.o

EXEFILES = dynamic_omm.x

%.o: %.f
	${F77} ${F77FLAGS} ${OPTFLAGS} $< -o $@ 

%.o: %.c
	${CC} ${CFLAGS} ${INCLUDEDIR} ${OPTFLAGS} $<

%.x: %.o libtinker.a
	${F77} ${LINKFLAGS} -o $@ $^ ${LIBS}; strip $@

all:	${EXEFILES}

clean:
	rm -f *.o *.mod *.a *.x

listing:
	cat *.f *.cpp > tinker.txt

ommstuff.o: ommstuff.cpp
	${CXX} ${CFLAGS} ${OPTFLAGS} -I${OPENMM_INCLUDE_DIR} ommstuff.cpp

dynamic_omm.x: dynamic_omm.o ommstuff.o libtinker.a
	${F77} ${LINKFLAGS} ${OPENMM_LIB_DIR} -o $@ $^ ${LIBS} ${OPENMM_LIBS}

rename:
	mv  dynamic_omm.x   $(BINDIR)/dynamic_omm

remove_links:
	rm -f $(LINKDIR)/dynamic_omm

create_links:
	ln -s $(BINDIR)/dynamic_omm   $(LINKDIR)/dynamic_omm

libtinker.a: ${OBJS} 
	ar ${LIBFLAGS} libtinker.a \
        action.o \
        active.o \
        align.o \
        analysis.o \
        analyz.o \
        angang.o \
        angbnd.o \
        angles.o \
        angpot.o \
        angtor.o \
        argue.o \
        ascii.o \
        atmlst.o \
        atomid.o \
        atoms.o \
        attach.o \
        basefile.o \
        bath.o \
        beeman.o \
        bicubic.o \
        bitor.o \
        bitors.o \
        bndpot.o \
        bndstr.o \
        bonds.o \
        born.o \
        bound.o \
        bounds.o \
        boxes.o \
        bussi.o \
        calendar.o \
        cell.o \
        center.o \
        charge.o \
        chgpot.o \
        chkpole.o \
        chkring.o \
        chkxyz.o \
        cholesky.o \
        chrono.o \
        chunks.o \
        clock.o \
        cluster.o \
        column.o \
        command.o \
        connect.o \
        connolly.o \
        control.o \
        couple.o \
        cspline.o \
        cutoffs.o \
        deflate.o \
        delete.o \
        deriv.o \
        diagq.o \
        diffeq.o \
        dipole.o \
        disgeo.o \
        dma.o \
        domega.o \
        eangang.o \
        eangang1.o \
        eangang2.o \
        eangang3.o \
        eangle.o \
        eangle1.o \
        eangle2.o \
        eangle3.o \
        eangtor.o \
        eangtor1.o \
        eangtor2.o \
        eangtor3.o \
        ebond.o \
        ebond1.o \
        ebond2.o \
        ebond3.o \
        ebuck.o \
        ebuck1.o \
        ebuck2.o \
        ebuck3.o \
        echarge.o \
        echarge1.o \
        echarge2.o \
        echarge3.o \
        echgdpl.o \
        echgdpl1.o \
        echgdpl2.o \
        echgdpl3.o \
        edipole.o \
        edipole1.o \
        edipole2.o \
        edipole3.o \
        egauss.o \
        egauss1.o \
        egauss2.o \
        egauss3.o \
        egeom.o \
        egeom1.o \
        egeom2.o \
        egeom3.o \
        ehal.o \
        ehal1.o \
        ehal2.o \
        ehal3.o \
        eimprop.o \
        eimprop1.o \
        eimprop2.o \
        eimprop3.o \
        eimptor.o \
        eimptor1.o \
        eimptor2.o \
        eimptor3.o \
        elj.o \
        elj1.o \
        elj2.o \
        elj3.o \
        embed.o \
        emetal.o \
        emetal1.o \
        emetal2.o \
        emetal3.o \
        emm3hb.o \
        emm3hb1.o \
        emm3hb2.o \
        emm3hb3.o \
        empole.o \
        empole1.o \
        empole2.o \
        empole3.o \
        energi.o \
        energy.o \
        eopbend.o \
        eopbend1.o \
        eopbend2.o \
        eopbend3.o \
        eopdist.o \
        eopdist1.o \
        eopdist2.o \
        eopdist3.o \
        epitors.o \
        epitors1.o \
        epitors2.o \
        epitors3.o \
        epolar.o \
        epolar1.o \
        epolar2.o \
        epolar3.o \
        erf.o \
        erxnfld.o \
        erxnfld1.o \
        erxnfld2.o \
        erxnfld3.o \
        esolv.o \
        esolv1.o \
        esolv2.o \
        esolv3.o \
        estrbnd.o \
        estrbnd1.o \
        estrbnd2.o \
        estrbnd3.o \
        estrtor.o \
        estrtor1.o \
        estrtor2.o \
        estrtor3.o \
        etors.o \
        etors1.o \
        etors2.o \
        etors3.o \
        etortor.o \
        etortor1.o \
        etortor2.o \
        etortor3.o \
        eurey.o \
        eurey1.o \
        eurey2.o \
        eurey3.o \
        evcorr.o \
        ewald.o \
        extra.o \
        extra1.o \
        extra2.o \
        extra3.o \
        faces.o \
        fatal.o \
        fft.o \
        fft3d.o \
        fftpack.o \
        field.o \
        fields.o \
        files.o \
        final.o \
        flatten.o \
        fracs.o \
        freeunit.o \
        freeze.o \
        geometry.o \
        getint.o \
        getkey.o \
        getmol.o \
        getmol2.o \
        getnumb.o \
        getpdb.o \
        getprm.o \
        getref.o \
        getstring.o \
        gettext.o \
        getword.o \
        getxyz.o \
        ghmcstep.o \
        gkstuf.o \
        gradient.o \
        gradrgd.o \
        gradrot.o \
        group.o \
        groups.o \
        grpline.o \
        gyrate.o \
        hescut.o \
        hessian.o \
        hessn.o \
        hessrgd.o \
        hessrot.o \
        hpmf.o \
        hybrid.o \
        image.o \
        impose.o \
        improp.o \
        imptor.o \
        induce.o \
        inertia.o \
        inform.o \
        initatom.o \
        initial.o \
        initprm.o \
        initres.o \
        initrot.o \
        insert.o \
        inter.o \
        invbeta.o \
        invert.o \
        iounit.o \
        jacobi.o \
        kanang.o \
        kangang.o \
        kangle.o \
        kangs.o \
        kangtor.o \
        kantor.o \
        katom.o \
        katoms.o \
        kbond.o \
        kbonds.o \
        kcharge.o \
        kchrge.o \
        kdipol.o \
        kdipole.o \
        kewald.o \
        kextra.o \
        keys.o \
        kgeom.o \
        khbond.o \
        kimprop.o \
        kimptor.o \
        kinetic.o \
        kiprop.o \
        kitors.o \
        kmetal.o \
        kmpole.o \
        kmulti.o \
        kopbend.o \
        kopbnd.o \
        kopdist.o \
        kopdst.o \
        korbit.o \
        korbs.o \
        kpitor.o \
        kpitors.o \
        kpolar.o \
        kpolr.o \
        ksolv.o \
        kstbnd.o \
        kstrbnd.o \
        kstrtor.o \
        ksttor.o \
        ktors.o \
        ktorsn.o \
        ktortor.o \
        ktrtor.o \
        kurey.o \
        kurybr.o \
        kvdw.o \
        kvdwpr.o \
        kvdws.o \
        lattice.o \
        lbfgs.o \
        light.o \
        lights.o \
        limits.o \
        linmin.o \
        makeint.o \
        makeref.o \
        makexyz.o \
        math.o \
        maxwell.o \
        mdinit.o \
        mdrest.o \
        mdsave.o \
        mdstat.o \
        mdstuf.o \
        mechanic.o \
        merck.o \
        merge.o \
        minima.o \
        molcul.o \
        moldyn.o \
        molecule.o \
        moment.o \
        moments.o \
        mplpot.o \
        mpole.o \
        mrecip.o \
        mutant.o \
        mutate.o \
        nblist.o \
        neigh.o \
        nextarg.o \
        nexttext.o \
        nonpol.o \
        nose.o \
        nspline.o \
        nucleo.o \
        number.o \
        numeral.o \
        numgrad.o \
        ocvm.o \
        omega.o \
        opbend.o \
        opdist.o \
        openend.o \
        openmm.o \
        openmp.o \
        optsave.o \
        orbital.o \
        orbits.o \
        orient.o \
        orthog.o \
        output.o \
        overlap.o \
        params.o \
        paths.o \
        pbstuf.o \
        pdb.o \
        phipsi.o \
        picalc.o \
        piorbs.o \
        pistuf.o \
        pitors.o \
        pme.o \
        pmestuff.o \
        pmpb.o \
        polar.o \
        polgrp.o \
        polpot.o \
        polymer.o \
        potent.o \
        potfit.o \
        precis.o \
        precise.o \
        pressure.o \
        prmkey.o \
        promo.o \
        prtdyn.o \
        prterr.o \
        prtint.o \
        prtmol2.o \
        prtpdb.o \
        prtprm.o \
        prtseq.o \
        prtxyz.o \
        ptable.o \
        qmstuf.o \
        qrfact.o \
        quatfit.o \
        random.o \
        rattle.o \
        readdyn.o \
        readgau.o \
        readint.o \
        readmol.o \
        readmol2.o \
        readpdb.o \
        readprm.o \
        readseq.o \
        readxyz.o \
        refer.o \
        replica.o \
        resdue.o \
        respa.o \
        restrn.o \
        rgddyn.o \
        rgdstep.o \
        rigid.o \
        ring.o \
        rings.o \
        rmsfit.o \
        rotbnd.o \
        rotlist.o \
        rotpole.o \
        rxnfld.o \
        rxnpot.o \
        scales.o \
        sdstep.o \
        search.o \
        sequen.o \
        server.o \
        shakeup.o \
        shunt.o \
        sigmoid.o \
        sizes.o \
        sktstuff.o \
        socket.o \
        solute.o \
        sort.o \
        square.o \
        stodyn.o \
        strbnd.o \
        strtor.o \
        suffix.o \
        surface.o \
        surfatom.o \
        switch.o \
        syntrn.o \
        tarray.o \
        temper.o \
        titles.o \
        tncg.o \
        torphase.o \
        torpot.o \
        torque.o \
        tors.o \
        torsions.o \
        tortor.o \
        tree.o \
        trimtext.o \
        unitcell.o \
        units.o \
        uprior.o \
        urey.o \
        urypot.o \
        usage.o \
        usolve.o \
        valfit.o \
        vdw.o \
        vdwpot.o \
        verlet.o \
        version.o \
        vibs.o \
        virial.o \
        volume.o \
        warp.o \
        xtals.o \
        xyzatm.o \
        zatom.o \
        zclose.o \
        zcoord.o
	${RANLIB} libtinker.a

###############################################################
##  Next Section has Explicit Dependencies on Include Files  ##
###############################################################

action.o:
active.o: atoms.o inform.o iounit.o keys.o sizes.o usage.o
alchemy.o: analyz.o atoms.o energi.o files.o inform.o iounit.o katoms.o mutant.o potent.o sizes.o units.o usage.o
align.o:
analysis.o: analyz.o atoms.o bound.o energi.o group.o inter.o iounit.o limits.o potent.o sizes.o vdwpot.o
analyz.o:
analyze.o: action.o analyz.o angang.o angbnd.o angpot.o angtor.o atomid.o atoms.o bitor.o bndstr.o bound.o boxes.o charge.o chgpot.o couple.o dipole.o energi.o ewald.o fields.o files.o improp.o imptor.o inform.o inter.o iounit.o korbs.o ktrtor.o kvdws.o limits.o math.o molcul.o moment.o mpole.o opbend.o opdist.o piorbs.o pistuf.o pitors.o pme.o polar.o polgrp.o potent.o sizes.o solute.o strbnd.o strtor.o tors.o tortor.o units.o urey.o vdw.o vdwpot.o virial.o
angang.o:
angbnd.o:
angles.o: angbnd.o atmlst.o atoms.o couple.o iounit.o sizes.o
angpot.o:
angtor.o:
anneal.o: atomid.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o mdstuf.o potent.o sizes.o solute.o usage.o warp.o
archive.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o sizes.o titles.o usage.o
argue.o:
ascii.o:
atmlst.o:
atomid.o: sizes.o
atoms.o: sizes.o
attach.o: atoms.o couple.o iounit.o sizes.o
bar.o: atoms.o boxes.o energi.o files.o inform.o iounit.o keys.o sizes.o units.o
basefile.o: ascii.o files.o
bath.o:
beeman.o: atomid.o atoms.o freeze.o mdstuf.o moldyn.o sizes.o units.o usage.o
bicubic.o:
bitor.o:
bitors.o: angbnd.o atoms.o bitor.o couple.o iounit.o sizes.o
bndpot.o:
bndstr.o:
bonds.o: atmlst.o atoms.o bndstr.o couple.o iounit.o sizes.o
born.o: atomid.o atoms.o bath.o chgpot.o couple.o deriv.o inform.o iounit.o math.o pbstuf.o sizes.o solute.o virial.o
bound.o:
bounds.o: atomid.o atoms.o boxes.o molcul.o sizes.o
boxes.o:
bussi.o: atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o sizes.o units.o usage.o
calendar.o:
cell.o:
center.o: align.o
charge.o:
chgpot.o:
chkpole.o: atoms.o mpole.o sizes.o
chkring.o: couple.o
chkxyz.o: atoms.o iounit.o sizes.o
cholesky.o:
chrono.o:
chunks.o:
clock.o: chrono.o
cluster.o: atomid.o atoms.o bound.o group.o inform.o iounit.o keys.o limits.o molcul.o sizes.o
column.o:
command.o: argue.o
connect.o: atoms.o couple.o sizes.o zclose.o zcoord.o
connolly.o: atoms.o faces.o inform.o iounit.o math.o sizes.o
control.o: argue.o inform.o keys.o output.o
correlate.o: ascii.o atomid.o atoms.o files.o inform.o iounit.o sizes.o
couple.o: sizes.o
crystal.o: atomid.o atoms.o bound.o boxes.o couple.o files.o iounit.o math.o molcul.o sizes.o
cspline.o: iounit.o
cutoffs.o: atoms.o bound.o hescut.o keys.o limits.o neigh.o polpot.o sizes.o tarray.o
deflate.o: iounit.o
delete.o: atomid.o atoms.o couple.o inform.o iounit.o sizes.o
deriv.o:
diagq.o:
diffeq.o: atoms.o iounit.o math.o sizes.o warp.o
diffuse.o: atomid.o atoms.o bound.o inform.o iounit.o molcul.o sizes.o usage.o
dipole.o:
disgeo.o:
distgeom.o: angbnd.o atomid.o atoms.o bndstr.o couple.o disgeo.o files.o inform.o iounit.o kvdws.o math.o refer.o restrn.o sizes.o tors.o
dma.o:
document.o: iounit.o
domega.o:
dynamic.o: atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o sizes.o solute.o stodyn.o usage.o
eangang.o: angang.o angbnd.o angpot.o atoms.o bound.o energi.o group.o math.o sizes.o usage.o
eangang1.o: angang.o angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o usage.o virial.o
eangang2.o: angang.o angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o sizes.o
eangang3.o: action.o analyz.o angang.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o usage.o
eangle.o: angbnd.o angpot.o atoms.o bound.o energi.o group.o math.o sizes.o usage.o
eangle1.o: angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o usage.o virial.o
eangle2.o: angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o sizes.o
eangle3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o usage.o
eangtor.o: angbnd.o angtor.o atoms.o bound.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o
eangtor1.o: angbnd.o angtor.o atoms.o bound.o deriv.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o virial.o
eangtor2.o: angbnd.o angtor.o atoms.o bound.o group.o hessn.o math.o sizes.o torpot.o tors.o
eangtor3.o: action.o analyz.o angbnd.o angtor.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o torpot.o tors.o usage.o
ebond.o: atoms.o bndpot.o bndstr.o bound.o energi.o group.o sizes.o usage.o
ebond1.o: atoms.o bndpot.o bndstr.o bound.o deriv.o energi.o group.o sizes.o usage.o virial.o
ebond2.o: atmlst.o atoms.o bndpot.o bndstr.o bound.o couple.o group.o hessn.o sizes.o
ebond3.o: action.o analyz.o atomid.o atoms.o bndpot.o bndstr.o bound.o energi.o group.o inform.o iounit.o sizes.o usage.o
ebuck.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o iounit.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
ebuck1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
ebuck2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o iounit.o math.o shunt.o sizes.o vdw.o vdwpot.o warp.o
ebuck3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
echarge.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o energi.o ewald.o group.o iounit.o light.o limits.o math.o neigh.o pme.o shunt.o sizes.o usage.o warp.o
echarge1.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o light.o limits.o math.o molcul.o neigh.o pme.o shunt.o sizes.o usage.o virial.o warp.o
echarge2.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o ewald.o group.o hessn.o limits.o math.o shunt.o sizes.o warp.o
echarge3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o energi.o ewald.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o warp.o
echgdpl.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o energi.o group.o shunt.o sizes.o units.o usage.o
echgdpl1.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o deriv.o dipole.o energi.o group.o inter.o molcul.o shunt.o sizes.o units.o usage.o virial.o
echgdpl2.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o group.o hessn.o shunt.o sizes.o units.o
echgdpl3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o energi.o group.o inform.o inter.o iounit.o molcul.o shunt.o sizes.o units.o usage.o
edipole.o: atoms.o bound.o cell.o chgpot.o dipole.o energi.o group.o shunt.o sizes.o units.o usage.o
edipole1.o: atoms.o bound.o cell.o chgpot.o deriv.o dipole.o energi.o group.o inter.o molcul.o shunt.o sizes.o units.o usage.o virial.o
edipole2.o: atoms.o bound.o cell.o chgpot.o dipole.o group.o hessn.o shunt.o sizes.o units.o
edipole3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o chgpot.o dipole.o energi.o group.o inform.o inter.o iounit.o molcul.o shunt.o sizes.o units.o usage.o
egauss.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
egauss1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
egauss2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o warp.o
egauss3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
egeom.o: atomid.o atoms.o bound.o energi.o group.o math.o molcul.o restrn.o sizes.o usage.o
egeom1.o: atomid.o atoms.o bound.o deriv.o energi.o group.o inter.o math.o molcul.o restrn.o sizes.o usage.o virial.o
egeom2.o: atomid.o atoms.o bound.o deriv.o group.o hessn.o math.o molcul.o restrn.o sizes.o
egeom3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o inter.o iounit.o math.o molcul.o restrn.o sizes.o usage.o
ehal.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
ehal1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o
ehal2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o
ehal3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o mutant.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
eimprop.o: atoms.o bound.o energi.o group.o improp.o math.o sizes.o torpot.o usage.o
eimprop1.o: atoms.o bound.o deriv.o energi.o group.o improp.o math.o sizes.o torpot.o usage.o virial.o
eimprop2.o: atoms.o bound.o group.o hessn.o improp.o math.o sizes.o torpot.o
eimprop3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o improp.o inform.o iounit.o math.o sizes.o torpot.o usage.o
eimptor.o: atoms.o bound.o energi.o group.o imptor.o sizes.o torpot.o usage.o
eimptor1.o: atoms.o bound.o deriv.o energi.o group.o imptor.o sizes.o torpot.o usage.o virial.o
eimptor2.o: atoms.o bound.o group.o hessn.o imptor.o sizes.o torpot.o
eimptor3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o imptor.o inform.o iounit.o math.o sizes.o torpot.o usage.o
elj.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o math.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
elj1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o inter.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o warp.o
elj2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o math.o shunt.o sizes.o vdw.o vdwpot.o warp.o
elj3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o warp.o
embed.o: angbnd.o atoms.o bndstr.o couple.o disgeo.o files.o inform.o iounit.o keys.o light.o math.o minima.o output.o refer.o restrn.o sizes.o tors.o units.o
emetal.o: atomid.o atoms.o couple.o energi.o kchrge.o sizes.o
emetal1.o: atomid.o atoms.o couple.o deriv.o energi.o kchrge.o sizes.o
emetal2.o:
emetal3.o: action.o analyz.o atomid.o atoms.o energi.o kchrge.o sizes.o
emm3hb.o: atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o energi.o group.o light.o limits.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
emm3hb1.o: atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o group.o inter.o iounit.o light.o limits.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o virial.o
emm3hb2.o: atmlst.o atomid.o atoms.o bndstr.o bound.o cell.o chgpot.o couple.o group.o hessn.o shunt.o sizes.o vdw.o vdwpot.o
emm3hb3.o: action.o analyz.o atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o neigh.o shunt.o sizes.o usage.o vdw.o vdwpot.o
empole.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o energi.o ewald.o group.o limits.o math.o mplpot.o mpole.o mrecip.o neigh.o pme.o potent.o shunt.o sizes.o usage.o
empole1.o: atoms.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o ewald.o group.o inter.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o virial.o
empole2.o: atoms.o bound.o boxes.o chgpot.o couple.o deriv.o group.o hessn.o limits.o molcul.o mplpot.o mpole.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o usage.o
empole3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o chgpot.o couple.o energi.o ewald.o group.o inform.o inter.o iounit.o limits.o math.o molcul.o mplpot.o mpole.o mrecip.o neigh.o pme.o potent.o shunt.o sizes.o usage.o
energi.o:
energy.o: bound.o energi.o iounit.o limits.o potent.o rigid.o sizes.o vdwpot.o
eopbend.o: angbnd.o angpot.o atoms.o bound.o energi.o fields.o group.o math.o opbend.o sizes.o usage.o
eopbend1.o: angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o opbend.o sizes.o usage.o virial.o
eopbend2.o: angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o opbend.o sizes.o
eopbend3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o opbend.o sizes.o usage.o
eopdist.o: angpot.o atoms.o bound.o energi.o group.o opdist.o sizes.o usage.o
eopdist1.o: angpot.o atoms.o bound.o deriv.o energi.o group.o opdist.o sizes.o usage.o virial.o
eopdist2.o: angpot.o atoms.o bound.o group.o hessn.o opdist.o sizes.o usage.o
eopdist3.o: action.o analyz.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o opdist.o sizes.o usage.o
epitors.o: atoms.o bound.o energi.o group.o pitors.o sizes.o torpot.o usage.o
epitors1.o: atoms.o bound.o deriv.o energi.o group.o pitors.o sizes.o torpot.o usage.o virial.o
epitors2.o: angbnd.o atoms.o bound.o deriv.o group.o hessn.o pitors.o sizes.o torpot.o usage.o
epitors3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o pitors.o sizes.o torpot.o usage.o
epolar.o: atoms.o boxes.o chgpot.o energi.o ewald.o limits.o math.o mpole.o polar.o potent.o sizes.o units.o
epolar1.o:
epolar2.o:
epolar3.o: action.o analyz.o atomid.o atoms.o boxes.o chgpot.o energi.o ewald.o inform.o iounit.o limits.o math.o mpole.o polar.o potent.o sizes.o units.o
erf.o: iounit.o math.o
erxnfld.o: atoms.o chgpot.o energi.o mpole.o rxnfld.o rxnpot.o shunt.o sizes.o usage.o
erxnfld1.o: atoms.o deriv.o energi.o sizes.o
erxnfld2.o:
erxnfld3.o: action.o analyz.o atomid.o atoms.o chgpot.o energi.o inform.o iounit.o mpole.o shunt.o sizes.o usage.o
esolv.o: atomid.o atoms.o bound.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o kvdws.o limits.o math.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o warp.o
esolv1.o: atomid.o atoms.o bound.o boxes.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o inter.o iounit.o kvdws.o limits.o math.o molcul.o mplpot.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o virial.o warp.o
esolv2.o: atoms.o charge.o chgpot.o deriv.o hessn.o math.o mpole.o potent.o shunt.o sizes.o solute.o warp.o
esolv3.o: action.o analyz.o atomid.o atoms.o bound.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o inform.o inter.o iounit.o kvdws.o limits.o math.o molcul.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o usage.o vdw.o warp.o
estrbnd.o: angbnd.o angpot.o atoms.o bndstr.o bound.o energi.o group.o math.o sizes.o strbnd.o usage.o
estrbnd1.o: angbnd.o angpot.o atoms.o bndstr.o bound.o deriv.o energi.o group.o math.o sizes.o strbnd.o usage.o virial.o
estrbnd2.o: angbnd.o angpot.o atoms.o bndstr.o bound.o group.o hessn.o math.o sizes.o strbnd.o
estrbnd3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bndstr.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o strbnd.o usage.o
estrtor.o: atoms.o bndstr.o bound.o energi.o group.o sizes.o strtor.o torpot.o tors.o usage.o
estrtor1.o: atoms.o bndstr.o bound.o deriv.o energi.o group.o sizes.o strtor.o torpot.o tors.o usage.o virial.o
estrtor2.o: atoms.o bndstr.o bound.o group.o hessn.o sizes.o strtor.o torpot.o tors.o
estrtor3.o: action.o analyz.o atomid.o atoms.o bndstr.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o strtor.o torpot.o tors.o usage.o
etors.o: atoms.o bound.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o warp.o
etors1.o: atoms.o bound.o deriv.o energi.o group.o math.o sizes.o torpot.o tors.o usage.o virial.o warp.o
etors2.o: atoms.o bound.o group.o hessn.o math.o sizes.o torpot.o tors.o warp.o
etors3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o sizes.o torpot.o tors.o usage.o warp.o
etortor.o: atomid.o atoms.o bitor.o bound.o couple.o energi.o group.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o
etortor1.o: atoms.o bitor.o bound.o deriv.o energi.o group.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o virial.o
etortor2.o: atoms.o bitor.o bound.o group.o hessn.o ktrtor.o math.o sizes.o torpot.o tortor.o units.o
etortor3.o: action.o analyz.o atoms.o bitor.o bound.o energi.o group.o inform.o iounit.o ktrtor.o math.o sizes.o torpot.o tortor.o usage.o
eurey.o: atoms.o bound.o energi.o group.o sizes.o urey.o urypot.o usage.o
eurey1.o: atoms.o bound.o deriv.o energi.o group.o sizes.o urey.o urypot.o usage.o virial.o
eurey2.o: atoms.o bound.o couple.o group.o hessn.o sizes.o urey.o urypot.o
eurey3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o sizes.o urey.o urypot.o usage.o
evcorr.o: atoms.o bound.o boxes.o limits.o math.o mutant.o shunt.o sizes.o vdw.o vdwpot.o
ewald.o:
extra.o: energi.o
extra1.o: atoms.o deriv.o energi.o sizes.o
extra2.o: atoms.o hessn.o sizes.o
extra3.o: action.o analyz.o atoms.o energi.o sizes.o
faces.o:
fatal.o: iounit.o
fft.o:
fft3d.o: fft.o openmp.o pme.o
fftpack.o: math.o
field.o: keys.o potent.o
fields.o:
files.o:
final.o: align.o analyz.o angang.o angbnd.o angtor.o atmlst.o bitor.o bndstr.o cell.o charge.o chunks.o couple.o deriv.o dipole.o disgeo.o domega.o faces.o fracs.o freeze.o group.o hessn.o hpmf.o improp.o imptor.o inform.o iounit.o katoms.o kchrge.o kpolr.o kvdws.o light.o merck.o molcul.o moldyn.o mpole.o mrecip.o mutant.o neigh.o nonpol.o omega.o opbend.o opdist.o orbits.o paths.o pbstuf.o pdb.o piorbs.o pistuf.o pitors.o pme.o polar.o polgrp.o potfit.o qmstuf.o refer.o restrn.o rgddyn.o rigid.o ring.o rotbnd.o sizes.o socket.o solute.o stodyn.o strbnd.o strtor.o syntrn.o tarray.o tors.o tortor.o uprior.o urey.o usage.o usolve.o vdw.o vibs.o warp.o
flatten.o: atoms.o fields.o inform.o iounit.o keys.o sizes.o warp.o
fracs.o:
freeunit.o: iounit.o
freeze.o:
gda.o: atoms.o files.o iounit.o minima.o potent.o sizes.o vdwpot.o warp.o
geometry.o: atoms.o math.o sizes.o
getint.o: atoms.o inform.o iounit.o output.o sizes.o
getkey.o: argue.o files.o iounit.o keys.o openmp.o
getmol.o: files.o iounit.o
getmol2.o: files.o iounit.o
getnumb.o: ascii.o
getpdb.o: iounit.o
getprm.o: files.o iounit.o keys.o params.o
getref.o: atomid.o atoms.o boxes.o couple.o files.o refer.o sizes.o titles.o
getstring.o: ascii.o
gettext.o: ascii.o
getword.o: ascii.o
getxyz.o: inform.o iounit.o output.o
ghmcstep.o: atomid.o atoms.o bath.o freeze.o iounit.o mdstuf.o moldyn.o sizes.o stodyn.o units.o usage.o virial.o
gkstuf.o: sizes.o
gradient.o: atoms.o bound.o couple.o deriv.o energi.o inter.o iounit.o limits.o potent.o rigid.o sizes.o vdwpot.o virial.o
gradrgd.o: atoms.o group.o rigid.o sizes.o
gradrot.o: atoms.o deriv.o domega.o omega.o potent.o rotbnd.o sizes.o
group.o:
groups.o: group.o
grpline.o: atomid.o atoms.o group.o rgddyn.o sizes.o
gyrate.o: atoms.o sizes.o usage.o
hescut.o:
hessian.o: atoms.o bound.o couple.o hescut.o hessn.o inform.o iounit.o limits.o mpole.o potent.o rigid.o sizes.o usage.o vdw.o vdwpot.o
hessn.o:
hessrgd.o: atoms.o group.o rigid.o sizes.o
hessrot.o: math.o omega.o sizes.o zcoord.o
hpmf.o:
hybrid.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o charge.o couple.o dipole.o imptor.o inform.o iounit.o kangs.o katoms.o kbonds.o kchrge.o kdipol.o kitors.o kstbnd.o ksttor.o ktorsn.o kvdws.o math.o mutant.o sizes.o strbnd.o strtor.o tors.o vdw.o vdwpot.o
image.o: boxes.o cell.o sizes.o
impose.o: align.o inform.o iounit.o sizes.o
improp.o:
imptor.o:
induce.o: atoms.o bound.o boxes.o cell.o couple.o ewald.o gkstuf.o group.o inform.o iounit.o limits.o math.o mpole.o neigh.o openmp.o pbstuf.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o sizes.o solute.o tarray.o units.o uprior.o usolve.o
inertia.o: atomid.o atoms.o iounit.o math.o sizes.o
inform.o:
initatom.o: ptable.o
initial.o: align.o atoms.o bath.o bound.o boxes.o cell.o files.o group.o inform.o iounit.o keys.o linmin.o minima.o molcul.o mutant.o neigh.o openmp.o output.o params.o pdb.o precis.o rigid.o scales.o sequen.o sizes.o socket.o warp.o zclose.o
initprm.o: angpot.o bndpot.o chgpot.o fields.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o math.o merck.o mplpot.o polpot.o rxnpot.o sizes.o solute.o torpot.o units.o urypot.o vdwpot.o
initres.o: resdue.o
initrot.o: atoms.o couple.o group.o inform.o iounit.o math.o omega.o potent.o restrn.o rotbnd.o sizes.o usage.o zcoord.o
insert.o: atomid.o atoms.o couple.o inform.o iounit.o sizes.o
intedit.o: atomid.o atoms.o files.o iounit.o katoms.o sizes.o zcoord.o
inter.o:
intxyz.o: files.o iounit.o titles.o
invbeta.o:
invert.o: iounit.o
iounit.o:
jacobi.o: iounit.o
kanang.o: sizes.o
kangang.o: angang.o angbnd.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o kanang.o keys.o potent.o sizes.o tors.o
kangle.o: angbnd.o angpot.o atomid.o atoms.o bndstr.o couple.o fields.o inform.o iounit.o kangs.o keys.o merck.o potent.o ring.o sizes.o usage.o
kangs.o:
kangtor.o: angtor.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o kantor.o keys.o potent.o sizes.o tors.o
kantor.o:
katom.o: atomid.o atoms.o couple.o inform.o iounit.o katoms.o keys.o sizes.o
katoms.o:
kbond.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o couple.o fields.o inform.o iounit.o kbonds.o keys.o merck.o potent.o sizes.o tors.o usage.o
kbonds.o:
kcharge.o: atomid.o atoms.o charge.o chgpot.o couple.o fields.o inform.o iounit.o kchrge.o keys.o merck.o potent.o sizes.o
kchrge.o:
kdipol.o:
kdipole.o: atmlst.o atoms.o bndstr.o couple.o dipole.o inform.o iounit.o kdipol.o keys.o potent.o sizes.o
kewald.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o math.o openmp.o pme.o sizes.o
kextra.o:
keys.o:
kgeom.o: atomid.o atoms.o bound.o couple.o group.o iounit.o keys.o molcul.o potent.o restrn.o sizes.o
khbond.o:
kimprop.o: atomid.o atoms.o couple.o improp.o inform.o iounit.o keys.o kiprop.o potent.o sizes.o tors.o
kimptor.o: atomid.o atoms.o couple.o imptor.o inform.o iounit.o keys.o kitors.o math.o potent.o sizes.o tors.o
kinetic.o: atomid.o atoms.o bath.o group.o mdstuf.o moldyn.o rgddyn.o sizes.o units.o usage.o
kiprop.o:
kitors.o:
kmetal.o:
kmpole.o: atoms.o couple.o inform.o iounit.o keys.o kmulti.o mpole.o polar.o polgrp.o potent.o sizes.o units.o
kmulti.o:
kopbend.o: angbnd.o angpot.o atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o kopbnd.o merck.o opbend.o potent.o sizes.o usage.o
kopbnd.o:
kopdist.o: angbnd.o angpot.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o keys.o kopdst.o opdist.o potent.o sizes.o
kopdst.o:
korbit.o: atomid.o atoms.o bndstr.o inform.o iounit.o keys.o korbs.o orbits.o piorbs.o pistuf.o sizes.o tors.o units.o
korbs.o: sizes.o
kpitor.o:
kpitors.o: atomid.o atoms.o bndstr.o couple.o inform.o iounit.o keys.o kpitor.o pitors.o potent.o sizes.o tors.o
kpolar.o: atoms.o couple.o inform.o iounit.o keys.o kpolr.o mpole.o neigh.o polar.o polgrp.o polpot.o potent.o sizes.o usolve.o
kpolr.o:
ksolv.o: angbnd.o atmlst.o atomid.o atoms.o bath.o bndstr.o chgpot.o couple.o gkstuf.o hpmf.o inform.o iounit.o keys.o kvdws.o math.o nonpol.o pbstuf.o polar.o potent.o ptable.o sizes.o solute.o
kstbnd.o:
kstrbnd.o: angbnd.o angpot.o atmlst.o atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o kstbnd.o merck.o potent.o ring.o sizes.o strbnd.o
kstrtor.o: atmlst.o atomid.o atoms.o couple.o inform.o iounit.o keys.o ksttor.o potent.o sizes.o strtor.o tors.o
ksttor.o:
ktors.o: atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o ktorsn.o math.o merck.o potent.o ring.o sizes.o tors.o usage.o
ktorsn.o:
ktortor.o: atomid.o atoms.o bitor.o inform.o iounit.o keys.o ktrtor.o potent.o sizes.o tortor.o
ktrtor.o:
kurey.o: angbnd.o atomid.o atoms.o inform.o iounit.o keys.o kurybr.o potent.o sizes.o urey.o
kurybr.o:
kvdw.o: atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o khbond.o kvdwpr.o kvdws.o math.o merck.o potent.o sizes.o vdw.o vdwpot.o
kvdwpr.o:
kvdws.o:
lattice.o: boxes.o cell.o inform.o iounit.o math.o
lbfgs.o: inform.o iounit.o keys.o linmin.o math.o minima.o output.o scales.o
light.o:
lights.o: bound.o boxes.o cell.o iounit.o light.o sizes.o
limits.o:
linmin.o:
makeint.o: atoms.o couple.o inform.o iounit.o math.o sizes.o zclose.o zcoord.o
makeref.o: atomid.o atoms.o boxes.o couple.o files.o refer.o sizes.o titles.o
makexyz.o: atoms.o sizes.o zcoord.o
math.o:
maxwell.o: units.o
mdinit.o: atomid.o atoms.o bath.o bound.o couple.o files.o freeze.o group.o inform.o iounit.o keys.o mdstuf.o molcul.o moldyn.o mpole.o rgddyn.o rigid.o sizes.o stodyn.o units.o uprior.o usage.o
mdrest.o: atomid.o atoms.o bound.o group.o inform.o iounit.o mdstuf.o moldyn.o rgddyn.o sizes.o units.o
mdsave.o: atomid.o atoms.o bound.o boxes.o files.o group.o inform.o iounit.o mdstuf.o moldyn.o mpole.o output.o polar.o potent.o rgddyn.o sizes.o socket.o titles.o units.o
mdstat.o: atoms.o bath.o bound.o boxes.o inform.o inter.o iounit.o limits.o mdstuf.o molcul.o sizes.o units.o usage.o warp.o
mdstuf.o:
mechanic.o: inform.o iounit.o limits.o potent.o vdwpot.o
merck.o: sizes.o
merge.o: atomid.o atoms.o couple.o iounit.o refer.o sizes.o
minima.o:
minimize.o: atoms.o files.o inform.o iounit.o keys.o scales.o sizes.o usage.o
minirot.o: files.o inform.o iounit.o keys.o math.o omega.o scales.o sizes.o zcoord.o
minrigid.o: files.o group.o inform.o iounit.o keys.o math.o output.o rigid.o sizes.o
molcul.o:
moldyn.o:
molecule.o: atomid.o atoms.o couple.o molcul.o sizes.o
molxyz.o: files.o iounit.o titles.o
moment.o:
moments.o: atomid.o atoms.o bound.o charge.o dipole.o limits.o moment.o mpole.o polar.o potent.o rigid.o sizes.o solute.o units.o usage.o
monte.o: atoms.o files.o inform.o iounit.o omega.o output.o sizes.o units.o usage.o zcoord.o
mplpot.o:
mpole.o:
mrecip.o:
mutant.o:
mutate.o: atomid.o atoms.o charge.o inform.o iounit.o katoms.o keys.o mpole.o mutant.o polar.o potent.o sizes.o
nblist.o: atoms.o bound.o boxes.o cell.o charge.o iounit.o light.o limits.o mpole.o neigh.o potent.o sizes.o vdw.o
neigh.o:
newton.o: atoms.o files.o inform.o iounit.o keys.o sizes.o usage.o
newtrot.o: files.o hescut.o inform.o iounit.o keys.o math.o omega.o sizes.o zcoord.o
nextarg.o: argue.o
nexttext.o:
nonpol.o:
nose.o: atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o sizes.o units.o usage.o virial.o
nspline.o:
nucleic.o: atoms.o couple.o files.o group.o inform.o iounit.o katoms.o math.o molcul.o nucleo.o output.o potent.o resdue.o restrn.o rigid.o sequen.o sizes.o titles.o usage.o
nucleo.o: sizes.o
number.o: inform.o iounit.o
numeral.o:
numgrad.o: atoms.o sizes.o
ocvm.o: inform.o iounit.o keys.o linmin.o math.o minima.o output.o potent.o scales.o
omega.o:
opbend.o:
opdist.o:
openend.o:
openmp.o:
optimize.o: atoms.o files.o inform.o iounit.o keys.o scales.o sizes.o usage.o
optirot.o: files.o inform.o iounit.o keys.o math.o omega.o scales.o sizes.o zcoord.o
optrigid.o: files.o group.o inform.o iounit.o keys.o math.o output.o rigid.o
optsave.o: atoms.o files.o iounit.o math.o omega.o output.o scales.o sizes.o socket.o usage.o zcoord.o
orbital.o: atomid.o atoms.o bndstr.o couple.o iounit.o keys.o piorbs.o potent.o sizes.o tors.o
orbits.o:
orient.o: atomid.o atoms.o group.o math.o rigid.o sizes.o
orthog.o:
output.o:
overlap.o: units.o
params.o:
path.o: align.o atomid.o atoms.o files.o inform.o iounit.o linmin.o minima.o output.o paths.o sizes.o
paths.o:
pbstuf.o:
pdb.o:
pdbxyz.o: atomid.o atoms.o couple.o fields.o files.o inform.o iounit.o katoms.o pdb.o resdue.o sequen.o sizes.o titles.o
phipsi.o: sizes.o
picalc.o: atomid.o atoms.o bndstr.o couple.o inform.o iounit.o orbits.o piorbs.o pistuf.o sizes.o tors.o units.o
piorbs.o:
pistuf.o:
pitors.o:
pme.o:
pmestuff.o: atoms.o boxes.o charge.o chunks.o mpole.o pme.o potent.o sizes.o
pmpb.o: iounit.o
polar.o:
polarize.o: atoms.o inform.o iounit.o molcul.o mpole.o polar.o polgrp.o polpot.o potent.o sizes.o units.o
poledit.o: atomid.o atoms.o couple.o dma.o files.o iounit.o keys.o kpolr.o mpole.o polar.o polgrp.o polpot.o potent.o sizes.o units.o
polgrp.o:
polpot.o:
polymer.o: atoms.o bndstr.o bound.o boxes.o iounit.o keys.o sizes.o
potent.o:
potential.o: atomid.o atoms.o charge.o chgpot.o dipole.o files.o inform.o iounit.o katoms.o keys.o math.o minima.o moment.o mpole.o neigh.o output.o polar.o potent.o potfit.o ptable.o refer.o sizes.o titles.o units.o
potfit.o: sizes.o
precis.o:
precise.o:
pressure.o: atomid.o atoms.o bath.o bound.o boxes.o group.o iounit.o math.o mdstuf.o molcul.o sizes.o units.o usage.o virial.o
prmedit.o: angpot.o bndpot.o iounit.o math.o params.o sizes.o urypot.o vdwpot.o
prmkey.o: angpot.o bndpot.o chgpot.o fields.o mplpot.o polpot.o potent.o rxnpot.o torpot.o urypot.o vdwpot.o
promo.o: iounit.o
protein.o: atomid.o atoms.o couple.o files.o group.o inform.o iounit.o katoms.o math.o molcul.o output.o phipsi.o potent.o resdue.o restrn.o rigid.o sequen.o sizes.o titles.o usage.o
prtdyn.o: atoms.o boxes.o files.o group.o mdstuf.o moldyn.o rgddyn.o sizes.o titles.o
prterr.o: files.o output.o
prtint.o: atomid.o atoms.o files.o inform.o sizes.o titles.o zclose.o zcoord.o
prtmol2.o: atomid.o atoms.o bndstr.o couple.o files.o iounit.o sizes.o titles.o
prtpdb.o: files.o pdb.o sequen.o sizes.o titles.o
prtprm.o: angpot.o bndpot.o chgpot.o fields.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o mplpot.o polpot.o sizes.o urypot.o vdwpot.o
prtseq.o: files.o sequen.o sizes.o
prtxyz.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o sizes.o titles.o
pss.o: atoms.o files.o hescut.o inform.o iounit.o math.o omega.o refer.o sizes.o tree.o warp.o zcoord.o
pssrigid.o: atoms.o files.o group.o inform.o iounit.o math.o minima.o molcul.o refer.o rigid.o sizes.o warp.o
pssrot.o: atoms.o files.o inform.o iounit.o math.o minima.o omega.o refer.o sizes.o warp.o zcoord.o
ptable.o:
qmstuf.o:
qrfact.o:
quatfit.o: align.o sizes.o
radial.o: argue.o atomid.o atoms.o bound.o boxes.o files.o inform.o iounit.o limits.o math.o molcul.o potent.o sizes.o
random.o: inform.o iounit.o keys.o math.o
rattle.o: atomid.o atoms.o freeze.o group.o inform.o iounit.o moldyn.o sizes.o units.o usage.o virial.o
readdyn.o: atoms.o boxes.o files.o group.o iounit.o mdstuf.o moldyn.o rgddyn.o sizes.o
readgau.o: ascii.o iounit.o qmstuf.o sizes.o units.o
readint.o: atomid.o atoms.o files.o inform.o iounit.o sizes.o titles.o zclose.o zcoord.o
readmol.o: atomid.o atoms.o couple.o files.o iounit.o ptable.o sizes.o titles.o
readmol2.o: atomid.o atoms.o couple.o files.o iounit.o ptable.o sizes.o titles.o
readpdb.o: files.o inform.o iounit.o pdb.o resdue.o sequen.o titles.o
readprm.o: fields.o iounit.o kanang.o kangs.o kantor.o katoms.o kbonds.o kchrge.o kdipol.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o merck.o params.o sizes.o
readseq.o: files.o iounit.o resdue.o sequen.o
readxyz.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o sizes.o titles.o
refer.o: sizes.o
replica.o: bound.o boxes.o cell.o inform.o iounit.o
resdue.o:
respa.o: atomid.o atoms.o freeze.o limits.o moldyn.o potent.o sizes.o units.o usage.o virial.o
restrn.o:
rgddyn.o:
rgdstep.o: atomid.o atoms.o bound.o group.o iounit.o rgddyn.o sizes.o units.o virial.o
rigid.o:
ring.o:
rings.o: angbnd.o atoms.o bitor.o bndstr.o couple.o inform.o iounit.o ring.o tors.o
rmsfit.o: align.o sizes.o
rotbnd.o:
rotlist.o: atoms.o couple.o iounit.o molcul.o rotbnd.o sizes.o zclose.o
rotpole.o: atoms.o mpole.o sizes.o
rxnfld.o:
rxnpot.o:
saddle.o: atoms.o inform.o iounit.o keys.o linmin.o minima.o sizes.o syntrn.o titles.o zcoord.o
scales.o: sizes.o
scan.o: atoms.o files.o inform.o iounit.o math.o minima.o omega.o output.o sizes.o zcoord.o
sdstep.o: atomid.o atoms.o bath.o couple.o freeze.o kvdws.o math.o mdstuf.o moldyn.o sizes.o stodyn.o units.o usage.o virial.o
search.o: linmin.o math.o sizes.o
sequen.o: sizes.o
server.o:
shakeup.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o bound.o couple.o freeze.o keys.o math.o molcul.o ring.o sizes.o usage.o
shunt.o:
sigmoid.o:
sizes.o:
sktstuff.o: atomid.o atoms.o charge.o couple.o deriv.o fields.o files.o inform.o iounit.o keys.o moldyn.o mpole.o polar.o potent.o sizes.o socket.o
sniffer.o: atoms.o files.o inform.o iounit.o linmin.o math.o minima.o output.o scales.o sizes.o usage.o
socket.o:
solute.o:
sort.o:
spacefill.o: atomid.o atoms.o files.o inform.o iounit.o kvdws.o math.o sizes.o usage.o
spectrum.o: files.o iounit.o math.o units.o
square.o: inform.o iounit.o keys.o minima.o sizes.o
stodyn.o:
strbnd.o:
strtor.o:
suffix.o: ascii.o
superpose.o: align.o atomid.o atoms.o bound.o files.o inform.o iounit.o sizes.o titles.o
surface.o: atoms.o inform.o iounit.o math.o sizes.o usage.o
surfatom.o: atoms.o iounit.o math.o sizes.o
switch.o: limits.o nonpol.o shunt.o
sybylxyz.o: files.o iounit.o titles.o
syntrn.o:
tarray.o:
temper.o: atomid.o atoms.o bath.o group.o mdstuf.o molcul.o moldyn.o rgddyn.o sizes.o units.o usage.o
testgrad.o: atoms.o deriv.o energi.o files.o inform.o inter.o iounit.o sizes.o solute.o usage.o
testhess.o: atoms.o files.o hescut.o inform.o iounit.o sizes.o usage.o
testpair.o: atoms.o deriv.o energi.o inform.o iounit.o light.o limits.o neigh.o polpot.o potent.o sizes.o tarray.o vdwpot.o
testpol.o: atoms.o bound.o inform.o iounit.o limits.o minima.o polar.o polpot.o potent.o rigid.o sizes.o units.o usage.o
testrot.o: domega.o energi.o inform.o iounit.o math.o omega.o sizes.o zcoord.o
timer.o: atoms.o hescut.o inform.o iounit.o limits.o sizes.o
timerot.o: iounit.o limits.o omega.o sizes.o
titles.o:
tncg.o: atoms.o hescut.o inform.o iounit.o keys.o linmin.o math.o minima.o output.o piorbs.o potent.o sizes.o
torphase.o:
torpot.o:
torque.o: atoms.o deriv.o mpole.o sizes.o
tors.o:
torsfit.o: atomid.o atoms.o files.o inform.o iounit.o keys.o ktorsn.o math.o output.o potent.o qmstuf.o restrn.o scales.o sizes.o tors.o usage.o
torsions.o: atoms.o bndstr.o couple.o iounit.o sizes.o tors.o
tortor.o:
tree.o:
trimtext.o:
unitcell.o: bound.o boxes.o iounit.o keys.o
units.o:
uprior.o:
urey.o:
urypot.o:
usage.o:
usolve.o:
valence.o: angbnd.o angpot.o atomid.o atoms.o bndpot.o bndstr.o couple.o files.o hescut.o inform.o iounit.o kangs.o kbonds.o keys.o kopbnd.o kstbnd.o ktorsn.o kurybr.o kvdws.o linmin.o math.o minima.o opbend.o output.o potent.o qmstuf.o scales.o sizes.o strbnd.o torpot.o tors.o units.o urey.o urypot.o usage.o valfit.o vdwpot.o
valfit.o:
vdw.o:
vdwpot.o:
verlet.o: atomid.o atoms.o freeze.o moldyn.o sizes.o units.o usage.o
version.o: iounit.o output.o
vibbig.o: atomid.o atoms.o bound.o couple.o files.o hescut.o hessn.o inform.o iounit.o keys.o limits.o mpole.o potent.o rigid.o sizes.o units.o usage.o vdw.o vdwpot.o vibs.o
vibrate.o: atomid.o atoms.o files.o hescut.o iounit.o math.o sizes.o units.o usage.o
vibrot.o: iounit.o omega.o sizes.o
vibs.o:
virial.o:
volume.o: atoms.o iounit.o math.o sizes.o
warp.o:
xtalfit.o: atomid.o atoms.o bound.o boxes.o charge.o dipole.o energi.o files.o fracs.o inform.o iounit.o kvdws.o limits.o math.o molcul.o mpole.o polar.o potent.o sizes.o vdw.o vdwpot.o xtals.o
xtalmin.o: atoms.o boxes.o files.o inform.o iounit.o keys.o math.o scales.o sizes.o
xtals.o:
xyzatm.o: atoms.o inform.o iounit.o math.o sizes.o
xyzedit.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o iounit.o limits.o math.o molcul.o ptable.o refer.o sizes.o titles.o units.o usage.o
xyzint.o: files.o iounit.o titles.o
xyzpdb.o: atomid.o atoms.o couple.o fields.o files.o inform.o molcul.o pdb.o resdue.o sequen.o sizes.o
xyzsybyl.o: files.o iounit.o sizes.o titles.o
zatom.o: angbnd.o atomid.o atoms.o bndstr.o fields.o iounit.o kangs.o katoms.o kbonds.o sizes.o zclose.o zcoord.o
zclose.o: sizes.o
zcoord.o: sizes.o

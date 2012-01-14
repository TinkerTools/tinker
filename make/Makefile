##
###################################################################
##                                                               ##
##  Makefile for Building the Tinker Molecular Modeling Package  ##
##                                                               ##
###################################################################
##
##  Invocation Options:
##
##   1. make all              Build all the TINKER executables
##   2. make rename           Move the executables to BINDIR
##   3. make remove_links     Create soft links in LINKDIR
##   4. make create_links     Remove soft links from LINKDIR
##   6. make listing          Concatenate source to tinker.txt
##   5. make clean            Delete objects and executables
##
##  Original version of this file is due to Peter Happersberger
##  and Jochen Buehler of the University of Konstanz, Jan. 1998,
##  Modifications by Reece Hart & Jay Ponder, Washington University
##
###################################################################

###################################################################
##  Master Directory Locations; Change as Needed for Local Site  ##
###################################################################

##
##  TINKERDIR    TINKER Distribution Directory
##  BINDIR       Hard Copies of TINKER Executables
##  LINKDIR      Linked Copies of TINKER Executables
##

TINKERDIR = /Users/ponder/tinker
BINDIR = $(TINKERDIR)/bin
LINKDIR = /usr/local/bin

####################################################################
##  Known Machine Types; Uncomment One of the Following Sections  ##
##  May Need Editing to Match Your Desired OS & Compiler Version  ##
####################################################################

##
##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible
##  Oper Sys: Fedora Core
##  Compiler: GNU gfortran
##

#F77 = /usr/bin/gfortran
#LIBS =
#F77FLAGS = -c
#OPTFLAGS = -O
#LIBFLAGS = -crusv
#LINKFLAGS =

##
##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible (also AMD)
##  Oper Sys: Fedora Core
##  Compiler: Intel Fortran for Linux 12.0
##

#F77 = /opt/intel/bin/ifort
#LIBS =
#F77FLAGS = -c -xHost -vec-report0
#OPTFLAGS = -O3 -no-ipo -no-prec-div
#LIBFLAGS = -crusv
#LINKFLAGS = -static-intel

##
##  Machine:  Generic Linux
##  CPU Type: Intel x86 Compatible (also AMD)
##  Oper Sys: Fedora Core
##  Compiler: Intel Fortran for Linux 12.0
##  Parallel: OpenMP
##

#F77 = /opt/intel/bin/ifort
#LIBS = -L$(TINKERDIR)/fftw/lib -lfftw3_omp -lfftw3
#F77FLAGS = -c -xHost -assume cc_omp
#OPTFLAGS = -O3 -no-ipo -no-prec-div -openmp
#LIBFLAGS = -crusv
#LINKFLAGS = $(OPTFLAGS) -static-intel

##
##  Machine:  Macintosh
##  CPU Type: Intel Xeon
##  Oper Sys: OS X 10.6 (Snow Leopard)
##  Compiler: Intel Fortran for Mac 12.0
##

#F77 = /opt/intel/bin/ifort
#LIBS =
#F77FLAGS = -c -axSSSE3 -vec-report0
#OPTFLAGS = -O3 -no-ipo -no-prec-div
#LIBFLAGS = -crusv
#LINKFLAGS = -static-intel -mmacosx-version-min=10.4

##
##  Machine:  Macintosh
##  CPU Type: Intel Xeon
##  Oper Sys: OS X 10.6 (Snow Leopard)
##  Compiler: Intel Fortran for Mac 12.0
##  Parallel: OpenMP
##

F77 = /opt/intel/bin/ifort
LIBS = -L$(TINKERDIR)/fftw/lib -lfftw3_omp -lfftw3
F77FLAGS = -c -axSSSE3 -assume cc_omp
OPTFLAGS = -O3 -no-ipo -no-prec-div -openmp
LIBFLAGS = -crusv
LINKFLAGS = $(OPTFLAGS) -static-intel -mmacosx-version-min=10.4

##
##  Machine:  Macintosh
##  CPU Type: PowerPC G4/G5
##  Oper Sys: OS X 10.4 (Tiger)
##  Compiler: IBM XLF 8.1
##

#F77 = /opt/ibmcmp/xlf/8.1/bin/xlf
#F77FLAGS = -c -qextname
#OPTFLAGS = -O2 -qmaxmem=-1
#LIBFLAGS = -r
#LINKFLAGS =

##
##  Machine:  HP Alpha
##  CPU Type: Alpha 21264
##  Oper Sys: HP Tru64 Unix
##  Compiler: HP Fortran
##

#F77 = /usr/bin/f77
#LIBS =
#F77FLAGS = -c
#OPTFLAGS = -fast -arch host -tune host
#LIBFLAGS = -rclvs
#LINKFLAGS = -fast -non_shared -om -WL,-om_no_inst_sched

##
##  Machine:  Silicon Graphics
##  CPU Type: MIPS R10000
##  Oper Sys: SGI Irix 6.5
##  Compiler: MIPSPro Fortran
##

#F77 = /bin/f77
#LIBS =
#F77FLAGS = -c
#OPTFLAGS = -O -mips4
#LIBFLAGS = -rclvs
#LINKFLAGS = -O -mips4

##
##  Machine:  SUN Workstation
##  CPU Type: UltraSPARC
##  Oper Sys: Solaris 4.0
##  Compiler: SUN Fortran
##

#F77 = /bin/f77
#LIBS =
#F77FLAGS = -c
#OPTFLAGS = -fast -temp=.
#LIBFLAGS = rcv
#LINKFLAGS = -fast

#################################################################
##  Should not be Necessary to Change Things Below this Point  ##
#################################################################

OBJS = active.o \
       alchemy.o \
       analysis.o \
       analyze.o \
       angles.o \
       anneal.o \
       archive.o \
       attach.o \
       basefile.o \
       beeman.o \
       bicubic.o \
       bitors.o \
       bonds.o \
       born.o \
       bounds.o \
       bussi.o \
       calendar.o \
       center.o \
       chkpole.o \
       chkring.o \
       chkxyz.o \
       cholesky.o \
       clock.o \
       cluster.o \
       column.o \
       command.o \
       connect.o \
       connolly.o \
       control.o \
       correlate.o \
       crystal.o \
       cspline.o \
       cutoffs.o \
       deflate.o \
       delete.o \
       diagq.o \
       diffeq.o \
       diffuse.o \
       distgeom.o \
       document.o \
       dynamic.o \
       eangang.o \
       eangang1.o \
       eangang2.o \
       eangang3.o \
       eangle.o \
       eangle1.o \
       eangle2.o \
       eangle3.o \
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
       extra.o \
       extra1.o \
       extra2.o \
       extra3.o \
       fatal.o \
       fft3d.o \
       fftpack.o \
       field.o \
       final.o \
       flatten.o \
       freeunit.o \
       gda.o \
       geometry.o \
       getint.o \
       getkey.o \
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
       gradient.o \
       gradrgd.o \
       gradrot.o \
       groups.o \
       grpline.o \
       gyrate.o \
       hessian.o \
       hessrgd.o \
       hessrot.o \
       hybrid.o \
       image.o \
       impose.o \
       induce.o \
       inertia.o \
       initial.o \
       initprm.o \
       initres.o \
       initrot.o \
       insert.o \
       intedit.o \
       intxyz.o \
       invbeta.o \
       invert.o \
       jacobi.o \
       kangang.o \
       kangle.o \
       katom.o \
       kbond.o \
       kcharge.o \
       kdipole.o \
       kewald.o \
       kgeom.o \
       kimprop.o \
       kimptor.o \
       kinetic.o \
       kmetal.o \
       kmpole.o \
       kopbend.o \
       kopdist.o \
       korbit.o \
       kpitors.o \
       kpolar.o \
       ksolv.o \
       kstrbnd.o \
       kstrtor.o \
       ktors.o \
       ktortor.o \
       kurey.o \
       kvdw.o \
       lattice.o \
       lbfgs.o \
       lights.o \
       makeint.o \
       makeref.o \
       makexyz.o \
       maxwell.o \
       mdinit.o \
       mdrest.o \
       mdsave.o \
       mdstat.o \
       mechanic.o \
       merge.o \
       minimize.o \
       minirot.o \
       minrigid.o \
       molecule.o \
       moments.o \
       monte.o \
       mutate.o \
       nblist.o \
       newton.o \
       newtrot.o \
       nextarg.o \
       nexttext.o \
       nose.o \
       nspline.o \
       nucleic.o \
       number.o \
       numeral.o \
       numgrad.o \
       ocvm.o \
       openend.o \
       optimize.o \
       optirot.o \
       optrigid.o \
       optsave.o \
       orbital.o \
       orient.o \
       orthog.o \
       overlap.o \
       path.o \
       pdbxyz.o \
       piscf.o \
       pmestuff.o \
       pmpb.o \
       polarize.o \
       poledit.o \
       polymer.o \
       potential.o \
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
       quatfit.o \
       radial.o \
       random.o \
       rattle.o \
       readdyn.o \
       readgau.o \
       readint.o \
       readmol2.o \
       readpdb.o \
       readprm.o \
       readseq.o \
       readxyz.o \
       replica.o \
       respa.o \
       rgdstep.o \
       rings.o \
       rmsfit.o \
       rotlist.o \
       rotpole.o \
       saddle.o \
       scan.o \
       sdstep.o \
       search.o \
       server.o \
       shakeup.o \
       sigmoid.o \
       sktstuff.o \
       sniffer.o \
       sort.o \
       spacefill.o \
       spectrum.o \
       square.o \
       suffix.o \
       superpose.o \
       surface.o \
       surfatom.o \
       switch.o \
       sybylxyz.o \
       temper.o \
       testgrad.o \
       testhess.o \
       testpair.o \
       testrot.o \
       timer.o \
       timerot.o \
       tncg.o \
       torphase.o \
       torque.o \
       torsfit.o \
       torsions.o \
       trimtext.o \
       unitcell.o \
       valence.o \
       verlet.o \
       version.o \
       vibbig.o \
       vibrate.o \
       vibrot.o \
       volume.o \
       xtalfit.o \
       xtalmin.o \
       xyzatm.o \
       xyzedit.o \
       xyzint.o \
       xyzpdb.o \
       xyzsybyl.o \
       zatom.o

EXEFILES = alchemy.x \
           analyze.x \
           anneal.x \
           archive.x \
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

.f.o:
	${F77} ${F77FLAGS} ${OPTFLAGS} $*.f 

pmpb.o:
	${F77} ${F77FLAGS} ${OPTFLAGS} pmpb.f

server.o:
	${F77} ${F77FLAGS} ${OPTFLAGS} server.f

%.x: %.o libtinker.a
	${F77} ${LINKFLAGS} -o $@ $^ ${LIBS}; strip $@

all:	${EXEFILES}

clean:
	rm -f *.o *.a *.x

listing:
	cat *.i *.f *.c > tinker.txt

rename:
	mv  alchemy.x    $(BINDIR)/alchemy
	mv  analyze.x    $(BINDIR)/analyze
	mv  anneal.x     $(BINDIR)/anneal
	mv  archive.x    $(BINDIR)/archive
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
        active.o \
        analysis.o \
        angles.o \
        attach.o \
        basefile.o \
        beeman.o \
        bicubic.o \
        bitors.o \
        bonds.o \
        born.o \
        bounds.o \
        bussi.o \
        calendar.o \
        center.o \
        chkpole.o \
        chkring.o \
        chkxyz.o \
        cholesky.o \
        clock.o \
        cluster.o \
        column.o \
        command.o \
        connect.o \
        connolly.o \
        control.o \
        cspline.o \
        cutoffs.o \
        deflate.o \
        delete.o \
        diagq.o \
        diffeq.o \
        eangang.o \
        eangang1.o \
        eangang2.o \
        eangang3.o \
        eangle.o \
        eangle1.o \
        eangle2.o \
        eangle3.o \
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
        extra.o \
        extra1.o \
        extra2.o \
        extra3.o \
        fatal.o \
        fft3d.o \
        fftpack.o \
        field.o \
        final.o \
        flatten.o \
        freeunit.o \
        geometry.o \
        getint.o \
        getkey.o \
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
        gradient.o \
        gradrgd.o \
        gradrot.o \
        groups.o \
        grpline.o \
        gyrate.o \
        hessian.o \
        hessrgd.o \
        hessrot.o \
        hybrid.o \
        image.o \
        impose.o \
        induce.o \
        inertia.o \
        initial.o \
        initprm.o \
        initres.o \
        initrot.o \
        insert.o \
        invbeta.o \
        invert.o \
        jacobi.o \
        kangang.o \
        kangle.o \
        katom.o \
        kbond.o \
        kcharge.o \
        kdipole.o \
        kewald.o \
        kgeom.o \
        kimprop.o \
        kimptor.o \
        kinetic.o \
        kmetal.o \
        kmpole.o \
        kopbend.o \
        kopdist.o \
        korbit.o \
        kpitors.o \
        kpolar.o \
        ksolv.o \
        kstrbnd.o \
        kstrtor.o \
        ktors.o \
        ktortor.o \
        kurey.o \
        kvdw.o \
        lattice.o \
        lbfgs.o \
        lights.o \
        makeint.o \
        makeref.o \
        makexyz.o \
        maxwell.o \
        mdinit.o \
        mdrest.o \
        mdsave.o \
        mdstat.o \
        mechanic.o \
        merge.o \
        molecule.o \
        moments.o \
        mutate.o \
        nblist.o \
        nextarg.o \
        nexttext.o \
        nose.o \
        nspline.o \
        number.o \
        numeral.o \
        numgrad.o \
        ocvm.o \
        openend.o \
        optsave.o \
        orbital.o \
        orient.o \
        orthog.o \
        overlap.o \
        piscf.o \
        pmestuff.o \
        pmpb.o \
        polymer.o \
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
        quatfit.o \
        random.o \
        rattle.o \
        readdyn.o \
        readgau.o \
        readint.o \
        readmol2.o \
        readpdb.o \
        readprm.o \
        readseq.o \
        readxyz.o \
        replica.o \
        respa.o \
        rgdstep.o \
        rings.o \
        rmsfit.o \
        rotlist.o \
        rotpole.o \
        sdstep.o \
        search.o \
        server.o \
        shakeup.o \
        sigmoid.o \
        sktstuff.o \
        sort.o \
        square.o \
        suffix.o \
        surface.o \
        surfatom.o \
        switch.o \
        temper.o \
        tncg.o \
        torphase.o \
        torque.o \
        torsions.o \
        trimtext.o \
        unitcell.o \
        verlet.o \
        version.o \
        volume.o \
        xyzatm.o \
        zatom.o

###############################################################
##  Next Section has Explicit Dependencies on Include Files  ##
###############################################################

active.o: atoms.i inform.i iounit.i keys.i sizes.i usage.i
alchemy.o: analyz.i atoms.i energi.i files.i inform.i iounit.i katoms.i mutant.i potent.i sizes.i units.i usage.i
analysis.o: analyz.i atoms.i bound.i cutoff.i energi.i group.i inter.i iounit.i potent.i sizes.i vdwpot.i
analyze.o: action.i analyz.i angang.i angle.i angpot.i atmtyp.i atoms.i bitor.i bond.i bound.i boxes.i charge.i chgpot.i cutoff.i dipole.i energi.i ewald.i fields.i files.i improp.i imptor.i inform.i inter.i iounit.i korbs.i ktrtor.i kvdws.i math.i molcul.i moment.i mpole.i opbend.i opdist.i piorbs.i pistuf.i pitors.i pme.i polar.i polgrp.i potent.i sizes.i solute.i strbnd.i strtor.i tors.i tortor.i units.i urey.i vdw.i vdwpot.i virial.i
angles.o: angle.i atmlst.i atoms.i couple.i iounit.i sizes.i
anneal.o: atmtyp.i atoms.i bath.i bond.i bound.i inform.i iounit.i mdstuf.i potent.i sizes.i solute.i usage.i warp.i
archive.o: atmtyp.i atoms.i boxes.i couple.i files.i inform.i iounit.i sizes.i titles.i usage.i
attach.o: atoms.i couple.i iounit.i sizes.i
basefile.o: ascii.i files.i
beeman.o: atmtyp.i atoms.i freeze.i mdstuf.i moldyn.i sizes.i units.i usage.i
bicubic.o:
bitors.o: angle.i bitor.i couple.i iounit.i sizes.i
bonds.o: atmlst.i atoms.i bond.i couple.i iounit.i sizes.i
born.o: atmtyp.i atoms.i bath.i chgpot.i couple.i deriv.i inform.i iounit.i math.i pbstuf.i potent.i sizes.i solute.i virial.i
bounds.o: atmtyp.i atoms.i boxes.i molcul.i sizes.i
bussi.o: atmtyp.i atoms.i bath.i boxes.i freeze.i mdstuf.i moldyn.i sizes.i units.i usage.i
calendar.o:
center.o: align.i sizes.i
chkpole.o: atoms.i mpole.i sizes.i
chkring.o: couple.i sizes.i
chkxyz.o: atoms.i iounit.i sizes.i
cholesky.o:
clock.o: chrono.i
cluster.o: atmtyp.i atoms.i bound.i cutoff.i group.i inform.i iounit.i keys.i molcul.i sizes.i
column.o: sizes.i
command.o: argue.i
connect.o: atoms.i couple.i sizes.i zclose.i zcoord.i
connolly.o: atoms.i faces.i inform.i iounit.i math.i sizes.i
control.o: argue.i inform.i keys.i output.i sizes.i
correlate.o: ascii.i atmtyp.i atoms.i files.i inform.i iounit.i sizes.i
crystal.o: atmtyp.i atoms.i bound.i boxes.i couple.i files.i iounit.i math.i molcul.i sizes.i
cspline.o: iounit.i
cutoffs.o: bound.i cutoff.i hescut.i keys.i neigh.i sizes.i
deflate.o: iounit.i
delete.o: atmtyp.i atoms.i couple.i inform.i iounit.i sizes.i
diagq.o:
diffeq.o: atoms.i iounit.i math.i sizes.i warp.i
diffuse.o: atmtyp.i atoms.i boxes.i iounit.i molcul.i sizes.i
distgeom.o: angle.i atmtyp.i atoms.i bond.i couple.i disgeo.i files.i inform.i iounit.i kgeoms.i math.i refer.i sizes.i tors.i
document.o: iounit.i
dynamic.o: atoms.i bath.i bond.i bound.i inform.i iounit.i keys.i mdstuf.i potent.i sizes.i solute.i stodyn.i usage.i
eangang.o: angang.i angle.i angpot.i atoms.i bound.i energi.i group.i math.i sizes.i usage.i
eangang1.o: angang.i angle.i angpot.i atoms.i bound.i deriv.i energi.i group.i math.i sizes.i usage.i virial.i
eangang2.o: angang.i angle.i angpot.i atoms.i bound.i deriv.i group.i hessn.i math.i sizes.i
eangang3.o: action.i analyz.i angang.i angle.i angpot.i atmtyp.i atoms.i bound.i energi.i group.i inform.i iounit.i math.i sizes.i usage.i
eangle.o: angle.i angpot.i atoms.i bound.i energi.i group.i math.i sizes.i usage.i
eangle1.o: angle.i angpot.i atoms.i bound.i deriv.i energi.i group.i math.i sizes.i usage.i virial.i
eangle2.o: angle.i angpot.i atoms.i bound.i deriv.i group.i hessn.i math.i sizes.i
eangle3.o: action.i analyz.i angle.i angpot.i atmtyp.i atoms.i bound.i energi.i group.i inform.i iounit.i math.i sizes.i usage.i
ebond.o: atoms.i bndpot.i bond.i bound.i energi.i group.i sizes.i usage.i
ebond1.o: atoms.i bndpot.i bond.i bound.i deriv.i energi.i group.i sizes.i usage.i virial.i
ebond2.o: atmlst.i atoms.i bndpot.i bond.i bound.i couple.i group.i hessn.i sizes.i
ebond3.o: action.i analyz.i atmtyp.i atoms.i bndpot.i bond.i bound.i energi.i group.i inform.i iounit.i sizes.i usage.i
ebuck.o: atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i energi.i group.i iounit.i light.i math.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i warp.i
ebuck1.o: atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i deriv.i energi.i group.i inter.i iounit.i light.i math.i molcul.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i virial.i warp.i
ebuck2.o: atmtyp.i atoms.i bound.i cell.i couple.i group.i hessn.i iounit.i math.i shunt.i sizes.i vdw.i vdwpot.i warp.i
ebuck3.o: action.i analyz.i atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i energi.i group.i inform.i inter.i iounit.i light.i math.i molcul.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i warp.i
echarge.o: atoms.i bound.i boxes.i cell.i charge.i chgpot.i couple.i cutoff.i energi.i ewald.i group.i iounit.i light.i math.i neigh.i pme.i shunt.i sizes.i usage.i warp.i
echarge1.o: atoms.i bound.i boxes.i cell.i charge.i chgpot.i couple.i cutoff.i deriv.i energi.i ewald.i group.i inter.i light.i math.i molcul.i neigh.i pme.i shunt.i sizes.i usage.i virial.i warp.i
echarge2.o: atoms.i bound.i cell.i charge.i chgpot.i couple.i cutoff.i ewald.i group.i hessn.i math.i shunt.i sizes.i warp.i
echarge3.o: action.i analyz.i atmtyp.i atoms.i bound.i boxes.i cell.i charge.i chgpot.i couple.i cutoff.i energi.i ewald.i group.i inform.i inter.i iounit.i light.i math.i molcul.i neigh.i shunt.i sizes.i usage.i warp.i
echgdpl.o: atoms.i bound.i cell.i charge.i chgpot.i couple.i dipole.i energi.i group.i shunt.i sizes.i units.i usage.i
echgdpl1.o: atoms.i bound.i cell.i charge.i chgpot.i couple.i deriv.i dipole.i energi.i group.i inter.i molcul.i shunt.i sizes.i units.i usage.i virial.i
echgdpl2.o: atoms.i bound.i cell.i charge.i chgpot.i couple.i dipole.i group.i hessn.i shunt.i sizes.i units.i
echgdpl3.o: action.i analyz.i atmtyp.i atoms.i bound.i cell.i charge.i chgpot.i couple.i dipole.i energi.i group.i inform.i inter.i iounit.i molcul.i shunt.i sizes.i units.i usage.i
edipole.o: atoms.i bound.i cell.i chgpot.i dipole.i energi.i group.i shunt.i sizes.i units.i usage.i
edipole1.o: atoms.i bound.i cell.i chgpot.i deriv.i dipole.i energi.i group.i inter.i molcul.i shunt.i sizes.i units.i usage.i virial.i
edipole2.o: atoms.i bound.i cell.i chgpot.i dipole.i group.i hessn.i shunt.i sizes.i units.i
edipole3.o: action.i analyz.i atmtyp.i atoms.i bound.i cell.i chgpot.i dipole.i energi.i group.i inform.i inter.i iounit.i molcul.i shunt.i sizes.i units.i usage.i
egauss.o: atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i energi.i group.i light.i math.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i warp.i
egauss1.o: atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i deriv.i energi.i group.i inter.i iounit.i light.i math.i molcul.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i virial.i warp.i
egauss2.o: atmtyp.i atoms.i bound.i cell.i couple.i group.i hessn.i shunt.i sizes.i vdw.i vdwpot.i warp.i
egauss3.o: action.i analyz.i atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i energi.i group.i inform.i inter.i iounit.i light.i math.i molcul.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i warp.i
egeom.o: atmtyp.i atoms.i bound.i energi.i group.i kgeoms.i math.i molcul.i sizes.i usage.i
egeom1.o: atmtyp.i atoms.i bound.i deriv.i energi.i group.i inter.i kgeoms.i math.i molcul.i sizes.i usage.i virial.i
egeom2.o: atmtyp.i atoms.i bound.i deriv.i group.i hessn.i kgeoms.i math.i molcul.i sizes.i
egeom3.o: action.i analyz.i atmtyp.i atoms.i bound.i energi.i group.i inform.i inter.i iounit.i kgeoms.i math.i molcul.i sizes.i usage.i
ehal.o: atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i energi.i group.i light.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i
ehal1.o: atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i deriv.i energi.i group.i inter.i iounit.i light.i molcul.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i virial.i
ehal2.o: atmtyp.i atoms.i bound.i cell.i couple.i group.i hessn.i shunt.i sizes.i vdw.i vdwpot.i
ehal3.o: action.i analyz.i atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i energi.i group.i inform.i inter.i iounit.i light.i molcul.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i
eimprop.o: atoms.i bound.i energi.i group.i improp.i math.i sizes.i torpot.i usage.i
eimprop1.o: atoms.i bound.i deriv.i energi.i group.i improp.i math.i sizes.i torpot.i usage.i virial.i
eimprop2.o: atoms.i bound.i group.i hessn.i improp.i math.i sizes.i torpot.i
eimprop3.o: action.i analyz.i atmtyp.i atoms.i bound.i energi.i group.i improp.i inform.i iounit.i math.i sizes.i torpot.i usage.i
eimptor.o: atoms.i bound.i energi.i group.i imptor.i sizes.i torpot.i usage.i
eimptor1.o: atoms.i bound.i deriv.i energi.i group.i imptor.i sizes.i torpot.i usage.i virial.i
eimptor2.o: atoms.i bound.i group.i hessn.i imptor.i sizes.i torpot.i
eimptor3.o: action.i analyz.i atmtyp.i atoms.i bound.i energi.i group.i imptor.i inform.i iounit.i math.i sizes.i torpot.i usage.i
elj.o: atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i energi.i group.i light.i math.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i warp.i
elj1.o: atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i deriv.i energi.i group.i inter.i light.i math.i molcul.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i virial.i warp.i
elj2.o: atmtyp.i atoms.i bound.i cell.i couple.i group.i hessn.i math.i shunt.i sizes.i vdw.i vdwpot.i warp.i
elj3.o: action.i analyz.i atmtyp.i atoms.i bound.i boxes.i cell.i couple.i cutoff.i energi.i group.i inform.i inter.i iounit.i light.i math.i molcul.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i warp.i
embed.o: angle.i atoms.i bond.i couple.i disgeo.i files.i inform.i iounit.i keys.i kgeoms.i light.i math.i minima.i output.i refer.i sizes.i tors.i units.i
emetal.o: atmtyp.i atoms.i couple.i energi.i kchrge.i sizes.i
emetal1.o: atmtyp.i atoms.i couple.i deriv.i energi.i kchrge.i sizes.i
emetal2.o:
emetal3.o: action.i analyz.i atmtyp.i atoms.i energi.i kchrge.i sizes.i
emm3hb.o: atmlst.i atmtyp.i atoms.i bond.i bound.i boxes.i cell.i chgpot.i couple.i cutoff.i energi.i group.i light.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i
emm3hb1.o: atmlst.i atmtyp.i atoms.i bond.i bound.i boxes.i cell.i chgpot.i couple.i cutoff.i deriv.i energi.i group.i inter.i iounit.i light.i molcul.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i virial.i
emm3hb2.o: atmlst.i atmtyp.i atoms.i bond.i bound.i cell.i chgpot.i couple.i group.i hessn.i shunt.i sizes.i vdw.i vdwpot.i
emm3hb3.o: action.i analyz.i atmlst.i atmtyp.i atoms.i bond.i bound.i boxes.i cell.i chgpot.i couple.i cutoff.i energi.i group.i inform.i inter.i iounit.i light.i molcul.i neigh.i shunt.i sizes.i usage.i vdw.i vdwpot.i
empole.o: atoms.i bound.i boxes.i cell.i chgpot.i couple.i cutoff.i energi.i ewald.i group.i math.i mplpot.i mpole.i neigh.i pme.i polar.i polgrp.i polpot.i potent.i shunt.i sizes.i usage.i
empole1.o: atoms.i bound.i boxes.i cell.i chgpot.i couple.i cutoff.i deriv.i energi.i ewald.i group.i inter.i math.i molcul.i mplpot.i mpole.i neigh.i pme.i polar.i polgrp.i polpot.i potent.i shunt.i sizes.i usage.i virial.i
empole2.o: atoms.i bound.i boxes.i chgpot.i couple.i cutoff.i deriv.i group.i hessn.i molcul.i mplpot.i mpole.i polar.i polgrp.i polpot.i potent.i shunt.i sizes.i usage.i
empole3.o: action.i analyz.i atmtyp.i atoms.i bound.i boxes.i cell.i chgpot.i couple.i cutoff.i energi.i ewald.i group.i inform.i inter.i iounit.i math.i molcul.i mplpot.i mpole.i neigh.i polar.i polgrp.i polpot.i potent.i shunt.i sizes.i usage.i
energy.o: bound.i cutoff.i energi.i iounit.i potent.i rigid.i sizes.i vdwpot.i
eopbend.o: angle.i angpot.i atoms.i bound.i energi.i fields.i group.i math.i opbend.i sizes.i usage.i
eopbend1.o: angle.i angpot.i atoms.i bound.i deriv.i energi.i group.i math.i opbend.i sizes.i usage.i virial.i
eopbend2.o: angle.i angpot.i atoms.i bound.i deriv.i group.i hessn.i math.i opbend.i sizes.i
eopbend3.o: action.i analyz.i angle.i angpot.i atmtyp.i atoms.i bound.i energi.i group.i inform.i iounit.i math.i opbend.i sizes.i usage.i
eopdist.o: angpot.i atoms.i bound.i energi.i group.i opdist.i sizes.i usage.i
eopdist1.o: angpot.i atoms.i bound.i deriv.i energi.i group.i opdist.i sizes.i usage.i virial.i
eopdist2.o: angpot.i atoms.i bound.i group.i hessn.i opdist.i sizes.i usage.i
eopdist3.o: action.i analyz.i angpot.i atmtyp.i atoms.i bound.i energi.i group.i inform.i iounit.i opdist.i sizes.i usage.i
epitors.o: atoms.i bound.i energi.i group.i pitors.i sizes.i torpot.i usage.i
epitors1.o: atoms.i bound.i deriv.i energi.i group.i pitors.i sizes.i torpot.i usage.i virial.i
epitors2.o: angle.i atoms.i bound.i deriv.i group.i hessn.i pitors.i sizes.i torpot.i usage.i
epitors3.o: action.i analyz.i atmtyp.i atoms.i bound.i energi.i group.i inform.i iounit.i math.i pitors.i sizes.i torpot.i usage.i
erf.o: iounit.i math.i
erxnfld.o: atoms.i chgpot.i energi.i mpole.i rxnfld.i rxnpot.i shunt.i sizes.i usage.i
erxnfld1.o: atoms.i deriv.i energi.i sizes.i
erxnfld2.o:
erxnfld3.o: action.i analyz.i atmtyp.i atoms.i chgpot.i energi.i inform.i iounit.i mpole.i shunt.i sizes.i usage.i
esolv.o: atmtyp.i atoms.i bound.i charge.i chgpot.i couple.i deriv.i energi.i gkstuf.i group.i hpmf.i kvdws.i math.i mpole.i npolar.i pbstuf.i polar.i polgrp.i polpot.i potent.i shunt.i sizes.i solute.i usage.i vdw.i warp.i
esolv1.o: atmtyp.i atoms.i bound.i boxes.i charge.i chgpot.i couple.i cutoff.i deriv.i energi.i gkstuf.i group.i hpmf.i inter.i iounit.i kvdws.i math.i molcul.i mplpot.i mpole.i npolar.i pbstuf.i polar.i polgrp.i polpot.i potent.i shunt.i sizes.i solute.i usage.i vdw.i virial.i warp.i
esolv2.o: atoms.i charge.i chgpot.i hessn.i math.i potent.i shunt.i sizes.i solute.i warp.i
esolv3.o: action.i analyz.i atmtyp.i atoms.i bound.i charge.i chgpot.i couple.i deriv.i energi.i gkstuf.i group.i hpmf.i inform.i inter.i iounit.i kvdws.i math.i molcul.i mpole.i npolar.i pbstuf.i polar.i polgrp.i polpot.i potent.i shunt.i sizes.i solute.i usage.i vdw.i warp.i
estrbnd.o: angle.i angpot.i atoms.i bond.i bound.i energi.i group.i math.i sizes.i strbnd.i usage.i
estrbnd1.o: angle.i angpot.i atoms.i bond.i bound.i deriv.i energi.i group.i math.i sizes.i strbnd.i usage.i virial.i
estrbnd2.o: angle.i angpot.i atoms.i bond.i bound.i group.i hessn.i math.i sizes.i strbnd.i
estrbnd3.o: action.i analyz.i angle.i angpot.i atmtyp.i atoms.i bond.i bound.i energi.i group.i inform.i iounit.i math.i sizes.i strbnd.i usage.i
estrtor.o: atoms.i bond.i bound.i energi.i group.i sizes.i strtor.i torpot.i tors.i usage.i
estrtor1.o: atoms.i bond.i bound.i deriv.i energi.i group.i sizes.i strtor.i torpot.i tors.i usage.i virial.i
estrtor2.o: atoms.i bond.i bound.i group.i hessn.i sizes.i strtor.i torpot.i tors.i
estrtor3.o: action.i analyz.i atmtyp.i atoms.i bond.i bound.i energi.i group.i inform.i iounit.i math.i sizes.i strtor.i torpot.i tors.i usage.i
etors.o: atoms.i bound.i energi.i group.i math.i sizes.i torpot.i tors.i usage.i warp.i
etors1.o: atoms.i bound.i deriv.i energi.i group.i math.i sizes.i torpot.i tors.i usage.i virial.i warp.i
etors2.o: atoms.i bound.i group.i hessn.i math.i sizes.i torpot.i tors.i warp.i
etors3.o: action.i analyz.i atmtyp.i atoms.i bound.i energi.i group.i inform.i iounit.i math.i sizes.i torpot.i tors.i usage.i warp.i
etortor.o: atmtyp.i atoms.i bitor.i bound.i couple.i energi.i group.i ktrtor.i math.i sizes.i torpot.i tortor.i usage.i
etortor1.o: atoms.i bitor.i bound.i deriv.i energi.i group.i ktrtor.i math.i sizes.i torpot.i tortor.i usage.i virial.i
etortor2.o: atoms.i bitor.i bound.i group.i hessn.i ktrtor.i math.i sizes.i torpot.i tortor.i units.i
etortor3.o: action.i analyz.i atoms.i bitor.i bound.i energi.i group.i inform.i iounit.i ktrtor.i math.i sizes.i torpot.i tortor.i usage.i
eurey.o: atoms.i bound.i energi.i group.i sizes.i urey.i urypot.i usage.i
eurey1.o: atoms.i bound.i deriv.i energi.i group.i sizes.i urey.i urypot.i usage.i virial.i
eurey2.o: atoms.i bound.i couple.i group.i hessn.i sizes.i urey.i urypot.i
eurey3.o: action.i analyz.i atmtyp.i atoms.i bound.i energi.i group.i inform.i iounit.i sizes.i urey.i urypot.i usage.i
evcorr.o: bound.i boxes.i cutoff.i math.i shunt.i sizes.i vdw.i vdwpot.i
extra.o: energi.i
extra1.o: atoms.i deriv.i energi.i sizes.i
extra2.o: atoms.i hessn.i sizes.i
extra3.o: action.i analyz.i atoms.i energi.i sizes.i
fatal.o: iounit.i
fft3d.o: fft.i openmp.i pme.i sizes.i
fftpack.o: math.i
field.o: keys.i potent.i sizes.i
final.o: chunks.i inform.i iounit.i pme.i sizes.i socket.i solute.i uprior.i usage.i vibs.i
flatten.o: atoms.i fields.i inform.i iounit.i keys.i sizes.i warp.i
freeunit.o: iounit.i
gda.o: atoms.i files.i iounit.i minima.i potent.i sizes.i vdwpot.i warp.i
geometry.o: atoms.i math.i sizes.i
getint.o: inform.i iounit.i output.i
getkey.o: argue.i files.i iounit.i keys.i openmp.i sizes.i
getmol2.o: files.i iounit.i
getnumb.o: ascii.i
getpdb.o: iounit.i
getprm.o: files.i iounit.i keys.i params.i sizes.i
getref.o: atmtyp.i atoms.i couple.i files.i refer.i sizes.i titles.i
getstring.o: ascii.i
gettext.o: ascii.i
getword.o: ascii.i
getxyz.o: inform.i iounit.i output.i
ghmcstep.o: atmtyp.i atoms.i bath.i freeze.i iounit.i mdstuf.i moldyn.i sizes.i stodyn.i units.i usage.i virial.i
gradient.o: atoms.i bound.i cutoff.i deriv.i energi.i inter.i iounit.i potent.i rigid.i sizes.i vdwpot.i virial.i
gradrgd.o: atoms.i group.i rigid.i sizes.i
gradrot.o: atoms.i deriv.i domega.i omega.i potent.i rotate.i sizes.i
groups.o: group.i sizes.i
grpline.o: atmtyp.i atoms.i group.i rgddyn.i sizes.i
gyrate.o: atoms.i sizes.i usage.i
hessian.o: atoms.i bound.i couple.i cutoff.i hescut.i hessn.i inform.i iounit.i mpole.i potent.i rigid.i sizes.i usage.i vdw.i vdwpot.i
hessrgd.o: atoms.i group.i rigid.i sizes.i
hessrot.o: math.i omega.i sizes.i zcoord.i
hybrid.o: angle.i atmlst.i atmtyp.i atoms.i bond.i charge.i couple.i dipole.i imptor.i inform.i iounit.i kangs.i katoms.i kbonds.i kchrge.i kdipol.i kitors.i kpolr.i kstbnd.i ksttor.i ktorsn.i kvdws.i math.i mpole.i mutant.i polar.i sizes.i strbnd.i strtor.i tors.i vdw.i vdwpot.i
image.o: boxes.i cell.i sizes.i
impose.o: align.i inform.i iounit.i sizes.i
induce.o: atoms.i bound.i boxes.i cell.i couple.i cutoff.i ewald.i gkstuf.i group.i inform.i iounit.i math.i mpole.i neigh.i pbstuf.i pme.i polar.i polgrp.i polpot.i potent.i shunt.i sizes.i solute.i units.i uprior.i
inertia.o: atmtyp.i atoms.i iounit.i math.i sizes.i
initial.o: align.i atoms.i bath.i bound.i cell.i files.i group.i inform.i iounit.i keys.i linmin.i minima.i molcul.i mutant.i neigh.i openmp.i output.i params.i pdb.i precis.i rigid.i scales.i sequen.i sizes.i socket.i warp.i zclose.i
initprm.o: angpot.i bndpot.i chgpot.i fields.i kanang.i kangs.i katoms.i kbonds.i kchrge.i kdipol.i khbond.i kiprop.i kitors.i kmulti.i kopbnd.i kopdst.i korbs.i kpitor.i kpolr.i kstbnd.i ksttor.i ktorsn.i ktrtor.i kurybr.i kvdwpr.i kvdws.i math.i mplpot.i polpot.i rxnpot.i sizes.i solute.i torpot.i units.i urypot.i vdwpot.i
initres.o: resdue.i sizes.i
initrot.o: atoms.i couple.i group.i inform.i iounit.i kgeoms.i math.i omega.i rotate.i sizes.i usage.i zcoord.i
insert.o: atmtyp.i atoms.i couple.i inform.i iounit.i sizes.i
intedit.o: atmtyp.i atoms.i files.i iounit.i katoms.i sizes.i zcoord.i
intxyz.o: files.i iounit.i titles.i
invbeta.o:
invert.o: iounit.i
jacobi.o: iounit.i
kangang.o: angang.i angle.i atmlst.i atmtyp.i atoms.i couple.i iounit.i kanang.i keys.i potent.i sizes.i
kangle.o: angle.i angpot.i atmtyp.i atoms.i bond.i couple.i fields.i inform.i iounit.i kangs.i keys.i merck.i potent.i ring.i sizes.i usage.i
katom.o: atmtyp.i atoms.i couple.i inform.i iounit.i katoms.i keys.i sizes.i
kbond.o: angle.i atmlst.i atmtyp.i atoms.i bond.i couple.i fields.i inform.i iounit.i kbonds.i keys.i merck.i potent.i sizes.i tors.i usage.i
kcharge.o: atmtyp.i atoms.i charge.i chgpot.i couple.i fields.i inform.i iounit.i kchrge.i keys.i merck.i potent.i sizes.i
kdipole.o: atmlst.i atoms.i bond.i couple.i dipole.i inform.i iounit.i kdipol.i keys.i potent.i sizes.i
kewald.o: atoms.i bound.i boxes.i chunks.i cutoff.i ewald.i fft.i inform.i iounit.i keys.i math.i openmp.i pme.i sizes.i
kgeom.o: atmtyp.i atoms.i bound.i couple.i group.i iounit.i keys.i kgeoms.i molcul.i potent.i sizes.i
kimprop.o: atmtyp.i atoms.i couple.i improp.i inform.i iounit.i keys.i kiprop.i potent.i sizes.i
kimptor.o: atmtyp.i atoms.i couple.i imptor.i inform.i iounit.i keys.i kitors.i math.i potent.i sizes.i
kinetic.o: atmtyp.i atoms.i bath.i group.i mdstuf.i moldyn.i rgddyn.i sizes.i units.i usage.i
kmetal.o:
kmpole.o: atoms.i couple.i inform.i iounit.i keys.i kmulti.i mpole.i polar.i polgrp.i potent.i sizes.i units.i
kopbend.o: angle.i angpot.i atmtyp.i atoms.i couple.i fields.i inform.i iounit.i keys.i kopbnd.i merck.i opbend.i potent.i sizes.i usage.i
kopdist.o: angle.i angpot.i atmlst.i atmtyp.i atoms.i couple.i inform.i iounit.i keys.i kopdst.i opdist.i potent.i sizes.i
korbit.o: atmtyp.i atoms.i bond.i inform.i iounit.i keys.i korbs.i orbits.i piorbs.i pistuf.i sizes.i tors.i units.i
kpitors.o: atmtyp.i atoms.i bond.i couple.i inform.i iounit.i keys.i kpitor.i pitors.i potent.i sizes.i
kpolar.o: atoms.i couple.i inform.i iounit.i keys.i kpolr.i mpole.i polar.i polgrp.i polpot.i potent.i sizes.i
ksolv.o: angle.i atmlst.i atmtyp.i atoms.i bath.i bond.i chgpot.i couple.i gkstuf.i hpmf.i inform.i iounit.i keys.i kvdws.i math.i npolar.i pbstuf.i potent.i sizes.i solute.i
kstrbnd.o: angle.i angpot.i atmlst.i atmtyp.i atoms.i couple.i fields.i inform.i iounit.i keys.i kstbnd.i merck.i potent.i ring.i sizes.i strbnd.i
kstrtor.o: atmlst.i atmtyp.i atoms.i couple.i inform.i iounit.i keys.i ksttor.i potent.i sizes.i strtor.i tors.i
ktors.o: atmtyp.i atoms.i couple.i fields.i inform.i iounit.i keys.i ktorsn.i math.i merck.i potent.i ring.i sizes.i tors.i usage.i
ktortor.o: atmtyp.i atoms.i bitor.i inform.i iounit.i keys.i ktrtor.i potent.i sizes.i tortor.i
kurey.o: angle.i atmtyp.i atoms.i inform.i iounit.i keys.i kurybr.i potent.i sizes.i urey.i
kvdw.o: atmtyp.i atoms.i couple.i fields.i inform.i iounit.i keys.i khbond.i kvdwpr.i kvdws.i math.i merck.i potent.i sizes.i vdw.i vdwpot.i
lattice.o: boxes.i cell.i inform.i iounit.i math.i sizes.i
lbfgs.o: inform.i iounit.i keys.i linmin.i math.i minima.i output.i scales.i sizes.i
lights.o: bound.i boxes.i cell.i iounit.i light.i sizes.i
makeint.o: atoms.i couple.i inform.i iounit.i math.i sizes.i zclose.i zcoord.i
makeref.o: atmtyp.i atoms.i couple.i files.i refer.i sizes.i titles.i
makexyz.o: atoms.i sizes.i zcoord.i
maxwell.o: units.i
mdinit.o: atmtyp.i atoms.i bath.i bound.i files.i freeze.i group.i inform.i iounit.i keys.i mdstuf.i molcul.i moldyn.i mpole.i rgddyn.i rigid.i sizes.i stodyn.i units.i uprior.i usage.i
mdrest.o: atmtyp.i atoms.i bound.i group.i inform.i iounit.i mdstuf.i moldyn.i rgddyn.i sizes.i units.i
mdsave.o: atmtyp.i atoms.i bound.i boxes.i files.i group.i inform.i iounit.i mdstuf.i moldyn.i mpole.i output.i polar.i potent.i rgddyn.i sizes.i socket.i titles.i units.i
mdstat.o: atoms.i bath.i bound.i boxes.i cutoff.i inform.i inter.i iounit.i mdstuf.i molcul.i sizes.i units.i usage.i warp.i
mechanic.o: cutoff.i inform.i iounit.i potent.i vdwpot.i
merge.o: atmtyp.i atoms.i couple.i iounit.i refer.i sizes.i
minimize.o: atoms.i files.i inform.i iounit.i keys.i scales.i sizes.i usage.i
minirot.o: files.i inform.i iounit.i keys.i math.i omega.i scales.i sizes.i zcoord.i
minrigid.o: files.i group.i inform.i iounit.i keys.i math.i output.i rigid.i sizes.i
molecule.o: atmtyp.i atoms.i couple.i molcul.i sizes.i
moments.o: atmtyp.i atoms.i charge.i dipole.i moment.i mpole.i polar.i potent.i sizes.i solute.i units.i usage.i
monte.o: atoms.i files.i inform.i iounit.i omega.i output.i sizes.i units.i usage.i zcoord.i
mutate.o: atmtyp.i atoms.i iounit.i katoms.i keys.i mutant.i sizes.i
nblist.o: atoms.i bound.i boxes.i charge.i cutoff.i iounit.i light.i mpole.i neigh.i potent.i sizes.i vdw.i
newton.o: atoms.i files.i inform.i iounit.i keys.i sizes.i usage.i
newtrot.o: files.i hescut.i inform.i iounit.i keys.i math.i omega.i sizes.i zcoord.i
nextarg.o: argue.i
nexttext.o:
nose.o: atmtyp.i atoms.i bath.i boxes.i freeze.i mdstuf.i moldyn.i sizes.i units.i usage.i virial.i
nspline.o:
nucleic.o: atoms.i couple.i files.i group.i inform.i iounit.i katoms.i kgeoms.i math.i molcul.i nucleo.i output.i potent.i resdue.i rigid.i sequen.i sizes.i titles.i usage.i
number.o: inform.i iounit.i
numeral.o:
numgrad.o: atoms.i sizes.i
ocvm.o: inform.i iounit.i keys.i linmin.i math.i minima.i output.i potent.i scales.i sizes.i
openend.o:
optimize.o: atoms.i files.i inform.i iounit.i keys.i scales.i sizes.i usage.i
optirot.o: files.i inform.i iounit.i keys.i math.i omega.i scales.i sizes.i zcoord.i
optrigid.o: files.i group.i inform.i iounit.i keys.i math.i output.i rigid.i sizes.i
optsave.o: atoms.i files.i iounit.i math.i omega.i output.i scales.i sizes.i socket.i usage.i zcoord.i
orbital.o: atmtyp.i atoms.i bond.i couple.i iounit.i keys.i piorbs.i potent.i sizes.i tors.i
orient.o: atmtyp.i atoms.i group.i math.i rigid.i sizes.i
orthog.o:
overlap.o: units.i
path.o: align.i atmtyp.i atoms.i files.i inform.i iounit.i linmin.i minima.i output.i paths.i sizes.i
pdbxyz.o: atmtyp.i atoms.i couple.i fields.i files.i inform.i iounit.i katoms.i pdb.i resdue.i sequen.i sizes.i
piscf.o: atmtyp.i atoms.i bond.i border.i couple.i inform.i iounit.i orbits.i piorbs.i pistuf.i sizes.i tors.i units.i
pmestuff.o: atoms.i boxes.i charge.i chunks.i mpole.i pme.i potent.i sizes.i
pmpb.o: iounit.i
polarize.o: atoms.i inform.i iounit.i molcul.i mpole.i polar.i polgrp.i polpot.i potent.i sizes.i units.i
poledit.o: atmtyp.i atoms.i couple.i dma.i files.i iounit.i keys.i kpolr.i mpole.i polar.i polgrp.i polpot.i potent.i sizes.i units.i
polymer.o: atoms.i bond.i bound.i boxes.i iounit.i keys.i sizes.i
potential.o: atmtyp.i atoms.i bond.i charge.i chgpot.i couple.i dipole.i files.i inform.i iounit.i katoms.i kchrge.i kdipol.i keys.i kmulti.i kpolr.i math.i minima.i moment.i mpole.i output.i polar.i potent.i potfit.i refer.i sizes.i titles.i units.i
precise.o:
pressure.o: atmtyp.i atoms.i bath.i bound.i boxes.i group.i iounit.i math.i mdstuf.i molcul.i sizes.i units.i usage.i virial.i
prmedit.o: angpot.i bndpot.i iounit.i math.i params.i sizes.i urypot.i vdwpot.i
prmkey.o: angpot.i bndpot.i chgpot.i fields.i mplpot.i polpot.i potent.i rxnpot.i sizes.i torpot.i urypot.i vdwpot.i
promo.o: iounit.i
protein.o: atmtyp.i atoms.i couple.i files.i group.i inform.i iounit.i katoms.i kgeoms.i math.i molcul.i output.i phipsi.i potent.i resdue.i rigid.i sequen.i sizes.i titles.i usage.i
prtdyn.o: atoms.i boxes.i files.i group.i mdstuf.i moldyn.i rgddyn.i sizes.i titles.i
prterr.o: files.i output.i
prtint.o: atmtyp.i atoms.i files.i inform.i sizes.i titles.i zclose.i zcoord.i
prtmol2.o: atmtyp.i atoms.i bond.i couple.i files.i iounit.i sizes.i titles.i
prtpdb.o: files.i pdb.i sequen.i sizes.i titles.i
prtprm.o: angpot.i bndpot.i chgpot.i fields.i kanang.i kangs.i katoms.i kbonds.i kchrge.i kdipol.i khbond.i kiprop.i kitors.i kmulti.i kopbnd.i kopdst.i korbs.i kpitor.i kpolr.i kstbnd.i ksttor.i ktorsn.i ktrtor.i kurybr.i kvdwpr.i kvdws.i mplpot.i polpot.i sizes.i urypot.i vdwpot.i
prtseq.o: files.i sequen.i sizes.i
prtxyz.o: atmtyp.i atoms.i couple.i files.i inform.i sizes.i titles.i
pss.o: atoms.i files.i hescut.i inform.i iounit.i math.i omega.i refer.i sizes.i tree.i warp.i zcoord.i
pssrigid.o: atoms.i files.i group.i inform.i iounit.i math.i minima.i molcul.i refer.i rigid.i sizes.i warp.i
pssrot.o: atoms.i files.i inform.i iounit.i math.i minima.i omega.i refer.i sizes.i warp.i zcoord.i
quatfit.o: align.i sizes.i
radial.o: argue.i atmtyp.i atoms.i bound.i boxes.i cutoff.i files.i inform.i iounit.i math.i molcul.i potent.i sizes.i
random.o: inform.i iounit.i keys.i math.i sizes.i
rattle.o: atmtyp.i atoms.i freeze.i group.i inform.i iounit.i moldyn.i sizes.i units.i usage.i virial.i
readdyn.o: atoms.i boxes.i files.i group.i iounit.i mdstuf.i moldyn.i rgddyn.i sizes.i
readgau.o: ascii.i iounit.i qmstuf.i sizes.i units.i
readint.o: atmtyp.i atoms.i files.i inform.i iounit.i sizes.i titles.i zclose.i zcoord.i
readmol2.o: atmtyp.i atoms.i couple.i files.i iounit.i sizes.i titles.i
readpdb.o: files.i inform.i iounit.i pdb.i resdue.i sequen.i sizes.i titles.i
readprm.o: fields.i iounit.i kanang.i kangs.i katoms.i kbonds.i kchrge.i kdipol.i khbond.i kiprop.i kitors.i kmulti.i kopbnd.i kopdst.i korbs.i kpitor.i kpolr.i kstbnd.i ksttor.i ktorsn.i ktrtor.i kurybr.i kvdwpr.i kvdws.i merck.i params.i sizes.i
readseq.o: files.i iounit.i resdue.i sequen.i sizes.i
readxyz.o: atmtyp.i atoms.i couple.i files.i inform.i iounit.i sizes.i titles.i
replica.o: bound.i boxes.i cell.i inform.i iounit.i sizes.i
respa.o: atmtyp.i atoms.i cutoff.i freeze.i moldyn.i potent.i sizes.i units.i usage.i virial.i
rgdstep.o: atmtyp.i atoms.i bound.i group.i iounit.i rgddyn.i sizes.i units.i virial.i
rings.o: angle.i atoms.i bitor.i bond.i couple.i inform.i iounit.i ring.i sizes.i tors.i
rmsfit.o: align.i sizes.i
rotlist.o: atoms.i couple.i iounit.i molcul.i rotate.i sizes.i zclose.i
rotpole.o: atoms.i mpole.i sizes.i
saddle.o: atoms.i inform.i iounit.i keys.i linmin.i minima.i sizes.i syntrn.i titles.i zcoord.i
scan.o: atoms.i files.i inform.i iounit.i math.i minima.i omega.i output.i sizes.i zcoord.i
sdstep.o: atmtyp.i atoms.i bath.i couple.i freeze.i kvdws.i math.i mdstuf.i moldyn.i sizes.i stodyn.i units.i usage.i virial.i
search.o: linmin.i math.i sizes.i
server.o:
shakeup.o: angle.i atmlst.i atmtyp.i atoms.i bond.i bound.i couple.i freeze.i keys.i math.i molcul.i ring.i sizes.i usage.i
sigmoid.o:
sktstuff.o: atmtyp.i atoms.i charge.i couple.i deriv.i fields.i files.i inform.i iounit.i keys.i moldyn.i mpole.i polar.i potent.i sizes.i socket.i
sniffer.o: atoms.i files.i inform.i iounit.i linmin.i math.i minima.i output.i scales.i sizes.i usage.i
sort.o:
spacefill.o: atmtyp.i atoms.i files.i inform.i iounit.i kvdws.i math.i sizes.i usage.i
spectrum.o: files.i iounit.i math.i units.i
square.o: inform.i iounit.i keys.i minima.i sizes.i
suffix.o: ascii.i
superpose.o: align.i atmtyp.i atoms.i bound.i files.i inform.i iounit.i sizes.i titles.i
surface.o: atoms.i inform.i iounit.i math.i sizes.i usage.i
surfatom.o: atoms.i iounit.i math.i sizes.i
switch.o: cutoff.i npolar.i shunt.i sizes.i
sybylxyz.o: files.i iounit.i titles.i
temper.o: atmtyp.i atoms.i bath.i group.i mdstuf.i molcul.i moldyn.i rgddyn.i sizes.i units.i usage.i
testgrad.o: atoms.i deriv.i energi.i inform.i inter.i iounit.i sizes.i solute.i usage.i
testhess.o: atoms.i files.i hescut.i inform.i iounit.i sizes.i usage.i
testpair.o: atoms.i cutoff.i deriv.i energi.i inform.i iounit.i light.i neigh.i potent.i sizes.i vdwpot.i
testrot.o: domega.i energi.i inform.i iounit.i math.i omega.i sizes.i zcoord.i
timer.o: atoms.i cutoff.i hescut.i inform.i iounit.i sizes.i
timerot.o: cutoff.i iounit.i sizes.i
tncg.o: hescut.i inform.i iounit.i keys.i linmin.i math.i minima.i output.i piorbs.i potent.i sizes.i
torphase.o:
torque.o: atoms.i deriv.i mpole.i sizes.i
torsfit.o: atmtyp.i atoms.i files.i inform.i iounit.i keys.i kgeoms.i ktorsn.i math.i output.i potent.i qmstuf.i scales.i sizes.i tors.i usage.i
torsions.o: bond.i couple.i iounit.i sizes.i tors.i
trimtext.o:
unitcell.o: bound.i boxes.i iounit.i keys.i sizes.i
valence.o: angle.i angpot.i atmtyp.i atoms.i bndpot.i bond.i couple.i files.i hescut.i inform.i iounit.i kangs.i kbonds.i keys.i kopbnd.i kstbnd.i ktorsn.i kurybr.i kvdws.i linmin.i math.i minima.i opbend.i output.i potent.i qmstuf.i scales.i sizes.i strbnd.i torpot.i tors.i units.i urey.i urypot.i usage.i valfit.i vdwpot.i
verlet.o: atmtyp.i atoms.i freeze.i moldyn.i sizes.i units.i usage.i
version.o: iounit.i output.i
vibbig.o: atmtyp.i atoms.i bound.i couple.i cutoff.i files.i hescut.i hessn.i inform.i iounit.i keys.i mpole.i potent.i rigid.i sizes.i units.i usage.i vdw.i vdwpot.i vibs.i
vibrate.o: atmtyp.i atoms.i files.i hescut.i iounit.i math.i sizes.i units.i usage.i
vibrot.o: iounit.i omega.i sizes.i
volume.o: atoms.i iounit.i math.i sizes.i
xtalfit.o: atmtyp.i atoms.i bound.i boxes.i charge.i couple.i dipole.i files.i fracs.i iounit.i kvdws.i math.i molcul.i potent.i sizes.i vdw.i xtals.i
xtalmin.o: atoms.i boxes.i files.i inform.i iounit.i keys.i math.i scales.i sizes.i
xyzatm.o: atoms.i inform.i iounit.i math.i sizes.i
xyzedit.o: atmtyp.i atoms.i bound.i couple.i cutoff.i files.i iounit.i math.i molcul.i refer.i sizes.i titles.i units.i usage.i
xyzint.o: files.i iounit.i titles.i
xyzpdb.o: atmtyp.i atoms.i couple.i fields.i files.i inform.i molcul.i pdb.i resdue.i sequen.i sizes.i
xyzsybyl.o: files.i iounit.i sizes.i titles.i
zatom.o: atmtyp.i atoms.i fields.i iounit.i katoms.i sizes.i zclose.i zcoord.i

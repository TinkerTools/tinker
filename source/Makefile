##
###################################################################
##                                                               ##
##  Makefile for Building the Tinker Molecular Modeling Package  ##
##                                                               ##
###################################################################
##

# source directory
src := $(PWD)
# release / debug / profile / etc.
opt := release
# gfortran / gfortran8 / gfortran-8 / ifort / ifort14 / ifort-14 / etc.
F77 := gfortran
# FFTW directory
fftw := /home/kchung25/tinker/fftw

###################################################################
##  Master Directory Locations; Change as Needed for Local Site  ##
###################################################################

##  TINKERDIR       Tinker Distribution Directory
##  TINKER_LIBDIR   Libraries needed to build Tinker
##  BINDIR          Directory with Tinker Executables
##  LINKDIR         Linked Copies of Tinker Executables
##  FFTWDIR         FFTW Top-Level Directory
##  FFTW_LIB_DIR    Directory with FFTW Libraries
##  FFTW_LIBS       FFTW Libraries needed to build Tinker
##  APBSDIR         APBS Top-Level Directory
##  APBS_INC_DIR    Directory with APBS Include Files
##  APBS_LIB_DIR    Directory with APBS Libraries
##  APBS_LIBS       APBS Libraries needed to build Tinker

TINKERDIR := /home/kchung25/tinker_rel
TINKER_LIBDIR := $(TINKERDIR)/lib
BINDIR := $(TINKERDIR)/bin
LINKDIR := /usr/local/bin

ifeq ($(fftw), default__)
  FFTWDIR := $(TINKERDIR)/fftw
else
  FFTWDIR := $(fftw)
endif
FFTW_LIBDIR := -L$(FFTWDIR)/lib
FFTW_LIBS := -lfftw3_threads -lfftw3

APBSDIR := $(TINKERDIR)/apbs
APBS_INCDIR := -I$(APBSDIR)/include
APBS_LIBDIR := -L$(APBSDIR)/lib
APBS_LIBS := -lapbsmainroutines -lapbs -lmaloc -lapbsblas

######################
## Operating System ##
######################
# Darwin (macOS) / Linux / CYGWIN_NT
os__ := $(shell uname -s)

ifeq ($(os__), Darwin)
else ifeq ($(os__), Linux)
else ifeq ($(shell echo $(os__) | cut -c 1-9), CYGWIN_NT)
else
$(error Unknown OS -- $(os__); Please help with us)
endif

##################
## F77 Compiler ##
##################
use_gfortran__ := false
use_ifort := false
found__ := false

f77__ := $(shell echo $(F77) | cut -c 1-8)
ifeq ($(f77__), gfortran)
  use_gfortran__ := true
  found__ := true
endif
f77__ := $(shell echo $(F77) | cut -c 1-5)
ifeq ($(f77__), ifort)
  use_ifort__ := true
  found__ := true
endif
ifneq ($(found__), true)
$(error Unknown fortran compiler -- $(F77); Please help with us)
endif

###################
## Compile Flags ##
###################

ifeq ($(use_gfortran__), true)
  ifeq ($(os__), Linux)
    ifeq ($(opt), release)
      OPTFLAGS := -Ofast -msse3 -fopenmp
    else ifeq ($(opt), profile)
      OPTFLAGS :=
    else
      # debug
      OPTFLAGS :=
    endif
    F77FLAGS := -c
    LIBDIR := -L. -L$(TINKER_LIBDIR)/linux
    LIBS :=
    LIBFLAGS := -crusv
    RANLIB := ranlib
    LINKFLAGS := $(OPTFLAGS) -static-libgcc
    RENAME := rename_bin
  endif

  ifeq ($(os__), Darwin)
    ifeq ($(opt), release)
      OPTFLAGS := -Ofast -mssse3 -fopenmp
    else ifeq ($(opt), profile)
      OPTFLAGS :=
    else
      # debug
      OPTFLAGS :=
    endif
    F77FLAGS := -c
    LIBDIR := -L. -L$(TINKER_LIBDIR)/macos
    LIBS :=
    LIBFLAGS := -crusv
    RANLIB := ranlib -c
    LINKFLAGS := $(OPTFLAGS) -static-libgcc
    RENAME := rename_bin
  endif

  ifeq ($(shell echo $(os__) | cut -c 1-9), CYGWIN_NT)
    ifeq ($(opt), release)
      OPTFLAGS := -Ofast -msse3 -fopenmp
    else ifeq ($(opt), profile)
      OPTFLAGS :=
    else
      # debug
      OPTFLAGS :=
    endif
    F77FLAGS := -c
    LIBDIR := -L. -L$(TINKER_LIBDIR)/windows
    LIBS :=
    LIBFLAGS := -crusv
    RANLIB := ranlib
    LINKFLAGS := $(OPTFLAGS) -static
    RENAME := rename_exe
  endif
endif

ifeq ($(use_ifort__), true)
  ifeq ($(os__), Linux)
    ifeq ($(opt), release)
      OPTFLAGS := -O3 -no-ipo -no-prec-div -recursive -qopenmp
    else ifeq ($(opt), profile)
      OPTFLAGS :=
    else
      # debug
      OPTFLAGS :=
    endif
    F77FLAGS := -c -xHost
    LIBDIR := -L. -L$(TINKER_LIBDIR)/linux
    LIBS :=
    LIBFLAGS := -crusv
    RANLIB := echo
    LINKFLAGS := $(OPTFLAGS) -static-libgcc -static-intel
    RENAME := rename_bin
  endif

  ifeq ($(os__), Darwin)
    ifeq ($(opt), release)
      OPTFLAGS := -O3 -no-ipo -no-prec-div -mdynamic-no-pic -qopenmp
    else ifeq ($(opt), profile)
      OPTFLAGS :=
    else
      # debug
      OPTFLAGS :=
    endif
    F77FLAGS := -c -axSSSE3
    LIBDIR := -L. -L$(TINKER_LIBDIR)/macos
    LIBS :=
    LIBFLAGS := -crusv
    RANLIB := ranlib -c
    LINKFLAGS := $(OPTFLAGS) -static-intel -Wl,-stack_size,0x10000000
    RENAME := rename_bin
  endif
endif

#################################################################
##  Should not be Necessary to Change Things Below this Point  ##
#################################################################

LIBOBJS := action.o active.o align.o alterchg.o alterpol.o analysis.o \
 analyz.o angang.o angbnd.o angles.o angpot.o angtor.o argue.o \
 ascii.o atmlst.o atomid.o atoms.o attach.o baoab.o basefile.o \
 bath.o beeman.o bicubic.o bitor.o bitors.o bndpot.o bndstr.o \
 bonds.o born.o bound.o bounds.o boxes.o bussi.o calendar.o cell.o \
 center.o cflux.o charge.o chgpen.o chgpot.o chgtrn.o chkpole.o \
 chkring.o chkxyz.o cholesky.o chrono.o chunks.o clock.o cluster.o \
 column.o command.o connect.o connolly.o control.o couple.o cspline.o \
 ctrpot.o cutoffs.o damping.o dcflux.o deflate.o delete.o deriv.o \
 dexpol.o diagq.o diffeq.o dipole.o disgeo.o disp.o dma.o domega.o \
 dsppot.o eangang.o eangang1.o eangang2.o eangang3.o eangle.o \
 eangle1.o eangle2.o eangle3.o eangtor.o eangtor1.o eangtor2.o \
 eangtor3.o ebond.o ebond1.o ebond2.o ebond3.o ebuck.o ebuck1.o \
 ebuck2.o ebuck3.o echarge.o echarge1.o echarge2.o echarge3.o \
 echgdpl.o echgdpl1.o echgdpl2.o echgdpl3.o echgtrn.o echgtrn1.o \
 echgtrn2.o echgtrn3.o edipole.o edipole1.o edipole2.o edipole3.o \
 edisp.o edisp1.o edisp2.o edisp3.o egauss.o egauss1.o egauss2.o \
 egauss3.o egeom.o egeom1.o egeom2.o egeom3.o ehal.o ehal1.o ehal2.o \
 ehal3.o eimprop.o eimprop1.o eimprop2.o eimprop3.o eimptor.o \
 eimptor1.o eimptor2.o eimptor3.o elj.o elj1.o elj2.o elj3.o embed.o \
 emetal.o emetal1.o emetal2.o emetal3.o emm3hb.o emm3hb1.o emm3hb2.o \
 emm3hb3.o empole.o empole1.o empole2.o empole3.o energi.o energy.o \
 eopbend.o eopbend1.o eopbend2.o eopbend3.o eopdist.o eopdist1.o \
 eopdist2.o eopdist3.o epitors.o epitors1.o epitors2.o epitors3.o \
 epolar.o epolar1.o epolar2.o epolar3.o erepel.o erepel1.o erepel2.o \
 erepel3.o erf.o erxnfld.o erxnfld1.o erxnfld2.o erxnfld3.o esolv.o \
 esolv1.o esolv2.o esolv3.o estrbnd.o estrbnd1.o estrbnd2.o estrbnd3.o \
 estrtor.o estrtor1.o estrtor2.o estrtor3.o etors.o etors1.o etors2.o \
 etors3.o etortor.o etortor1.o etortor2.o etortor3.o eurey.o eurey1.o \
 eurey2.o eurey3.o evcorr.o ewald.o exfield.o expol.o exrepel1.o \
 exrepel3.o extfld.o extra.o extra1.o extra2.o extra3.o faces.o \
 fatal.o fft.o fft3d.o fftpack.o field.o fields.o files.o final.o \
 flatten.o fracs.o freeunit.o freeze.o geometry.o getarc.o getcart.o \
 getdcd.o getint.o getkey.o getmol.o getmol2.o getnumb.o getpdb.o \
 getprm.o getref.o getstring.o gettext.o getword.o getxyz.o ghmcstep.o \
 gkstuf.o gradient.o gradrgd.o gradrot.o group.o groups.o grpline.o \
 gyrate.o hescut.o hessian.o hessn.o hessrgd.o hessrot.o hpmf.o \
 hybrid.o ielscf.o image.o impose.o improp.o imptor.o induce.o \
 inertia.o inform.o initatom.o initial.o initprm.o initres.o initrot.o \
 insert.o inter.o invbeta.o invert.o iounit.o jacobi.o kanang.o \
 kangang.o kangle.o kangs.o kangtor.o kantor.o katom.o katoms.o \
 kbond.o kbonds.o kcflux.o kcharge.o kchgflx.o kchgtrn.o kchrge.o \
 kcpen.o kctrn.o kdipol.o kdipole.o kdisp.o kdsp.o kewald.o kexpl.o \
 kexpol.o kextra.o keys.o kgeom.o khbond.o kimprop.o kimptor.o \
 kinetic.o kiprop.o kitors.o kmetal.o kmpole.o kmulti.o kopbend.o \
 kopbnd.o kopdist.o kopdst.o korbit.o korbs.o kpitor.o kpitors.o \
 kpolar.o kpolpr.o kpolr.o krepel.o krepl.o ksolut.o ksolv.o kstbnd.o \
 kstrbnd.o kstrtor.o ksttor.o ktors.o ktorsn.o ktortor.o ktrtor.o \
 kurey.o kurybr.o kvdw.o kvdwpr.o kvdws.o kxrepel.o kxrepl.o lattice.o \
 lbfgs.o light.o lights.o limits.o linmin.o lusolve.o makeint.o \
 makeref.o makexyz.o math.o maxwell.o mdinit.o mdrest.o mdsave.o \
 mdstat.o mdstuf.o mechanic.o merck.o merge.o minima.o molcul.o \
 moldyn.o molecule.o moment.o moments.o mplpot.o mpole.o mrecip.o \
 mutant.o mutate.o nblist.o neigh.o nextarg.o nexttext.o nonpol.o \
 nose.o nspline.o nucleo.o number.o numeral.o numgrad.o ocvm.o \
 omega.o opbend.o opdist.o openend.o openmp.o optinit.o optsave.o \
 orbital.o orbits.o orient.o orthog.o output.o overlap.o params.o \
 paths.o pbstuf.o pdb.o phipsi.o picalc.o piorbs.o pistuf.o pitors.o \
 pme.o pmestuf.o pmpb.o polar.o polgrp.o polopt.o polpcg.o polpot.o \
 poltcg.o polymer.o potent.o potfit.o predict.o pressure.o prmkey.o \
 promo.o prtarc.o prtdcd.o prtdyn.o prterr.o prtint.o prtmol2.o \
 prtpdb.o prtprm.o prtseq.o prtxyz.o ptable.o qmstuf.o qrsolve.o \
 quatfit.o random.o rattle.o readcart.o readdcd.o readdyn.o readgau.o \
 readgdma.o readint.o readmol.o readmol2.o readpdb.o readprm.o \
 readseq.o readxyz.o refer.o repel.o replica.o reppot.o resdue.o \
 respa.o restrn.o rgddyn.o rgdstep.o rigid.o ring.o rings.o rmsfit.o \
 rotbnd.o rotlist.o rotpole.o rxnfld.o rxnpot.o scales.o sdstep.o \
 search.o sequen.o server.o setprm.o shakeup.o shunt.o sigmoid.o \
 simplex.o sizes.o sktstuf.o socket.o solpot.o solute.o sort.o \
 square.o stodyn.o strbnd.o strtor.o suffix.o surface.o surfatom.o \
 switch.o syntrn.o tarray.o tcgstuf.o temper.o titles.o tncg.o \
 torphase.o torpot.o torque.o tors.o torsions.o tortor.o tree.o \
 trimtext.o unitcell.o units.o uprior.o urey.o urypot.o usage.o \
 valfit.o vdw.o vdwpot.o verlet.o version.o vibs.o virial.o volume.o \
 warp.o xrepel.o xtals.o xyzatm.o zatom.o zclose.o zcoord.o

EXEOBJS := alchemy.o analyze.o anneal.o arcedit.o bar.o correlate.o \
 critical.o crystal.o diffuse.o distgeom.o document.o dynamic.o \
 freefix.o gda.o intedit.o intxyz.o minimize.o minirot.o minrigid.o \
 mol2xyz.o molxyz.o monte.o newton.o newtrot.o nucleic.o optimize.o \
 optirot.o optrigid.o path.o pdbxyz.o polarize.o poledit.o potential.o \
 prmedit.o protein.o pss.o pssrigid.o pssrot.o radial.o saddle.o \
 scan.o sniffer.o spacefill.o spectrum.o superpose.o testgrad.o \
 testhess.o testpair.o testpol.o testrot.o testvir.o timer.o timerot.o \
 torsfit.o valence.o vibbig.o vibrate.o vibrot.o xtalfit.o xtalmin.o \
 xyzedit.o xyzint.o xyzmol2.o xyzpdb.o

OBJS := $(LIBOBJS) $(EXEOBJS)

EXEFILES := alchemy.x analyze.x anneal.x arcedit.x bar.x correlate.x \
 critical.x crystal.x diffuse.x distgeom.x document.x dynamic.x \
 freefix.x gda.x intedit.x intxyz.x minimize.x minirot.x minrigid.x \
 mol2xyz.x molxyz.x monte.x newton.x newtrot.x nucleic.x optimize.x \
 optirot.x optrigid.x path.x pdbxyz.x polarize.x poledit.x potential.x \
 prmedit.x protein.x pss.x pssrigid.x pssrot.x radial.x saddle.x \
 scan.x sniffer.x spacefill.x spectrum.x superpose.x testgrad.x \
 testhess.x testpair.x testpol.x testrot.x testvir.x timer.x timerot.x \
 torsfit.x valence.x vibbig.x vibrate.x vibrot.x xtalfit.x xtalmin.x \
 xyzedit.x xyzint.x xyzmol2.x xyzpdb.x

%.o: $(src)/%.f
	$(F77) $(F77FLAGS) $(OPTFLAGS) $< -o $@

%.x: %.o libtinker.a
	$(F77) $(LINKFLAGS) -o $@ $(LIBDIR) $(FFTW_LIBDIR) $^ $(LIBS) $(FFTW_LIBS); strip $@

all: $(EXEFILES)

install: $(RENAME)

clean:
	rm -f *.o *.mod *.a *.x

listing:
	cat $(src)/*.f $(src)/*.c > tinker.txt

rename_bin:
	mv alchemy.x            $(BINDIR)/alchemy
	mv analyze.x            $(BINDIR)/analyze
	mv anneal.x             $(BINDIR)/anneal
	mv arcedit.x            $(BINDIR)/arcedit
	mv bar.x                $(BINDIR)/bar
	mv correlate.x          $(BINDIR)/correlate
	mv critical.x           $(BINDIR)/critical
	mv crystal.x            $(BINDIR)/crystal
	mv diffuse.x            $(BINDIR)/diffuse
	mv distgeom.x           $(BINDIR)/distgeom
	mv document.x           $(BINDIR)/document
	mv dynamic.x            $(BINDIR)/dynamic
	mv freefix.x            $(BINDIR)/freefix
	mv gda.x                $(BINDIR)/gda
	mv intedit.x            $(BINDIR)/intedit
	mv intxyz.x             $(BINDIR)/intxyz
	mv minimize.x           $(BINDIR)/minimize
	mv minirot.x            $(BINDIR)/minirot
	mv minrigid.x           $(BINDIR)/minrigid
	mv mol2xyz.x            $(BINDIR)/mol2xyz
	mv molxyz.x             $(BINDIR)/molxyz
	mv monte.x              $(BINDIR)/monte
	mv newton.x             $(BINDIR)/newton
	mv newtrot.x            $(BINDIR)/newtrot
	mv nucleic.x            $(BINDIR)/nucleic
	mv optimize.x           $(BINDIR)/optimize
	mv optirot.x            $(BINDIR)/optirot
	mv optrigid.x           $(BINDIR)/optrigid
	mv path.x               $(BINDIR)/path
	mv pdbxyz.x             $(BINDIR)/pdbxyz
	mv polarize.x           $(BINDIR)/polarize
	mv poledit.x            $(BINDIR)/poledit
	mv potential.x          $(BINDIR)/potential
	mv prmedit.x            $(BINDIR)/prmedit
	mv protein.x            $(BINDIR)/protein
	mv pss.x                $(BINDIR)/pss
	mv pssrigid.x           $(BINDIR)/pssrigid
	mv pssrot.x             $(BINDIR)/pssrot
	mv radial.x             $(BINDIR)/radial
	mv saddle.x             $(BINDIR)/saddle
	mv scan.x               $(BINDIR)/scan
	mv sniffer.x            $(BINDIR)/sniffer
	mv spacefill.x          $(BINDIR)/spacefill
	mv spectrum.x           $(BINDIR)/spectrum
	mv superpose.x          $(BINDIR)/superpose
	mv testgrad.x           $(BINDIR)/testgrad
	mv testhess.x           $(BINDIR)/testhess
	mv testpair.x           $(BINDIR)/testpair
	mv testpol.x            $(BINDIR)/testpol
	mv testrot.x            $(BINDIR)/testrot
	mv testvir.x            $(BINDIR)/testvir
	mv timer.x              $(BINDIR)/timer
	mv timerot.x            $(BINDIR)/timerot
	mv torsfit.x            $(BINDIR)/torsfit
	mv valence.x            $(BINDIR)/valence
	mv vibbig.x             $(BINDIR)/vibbig
	mv vibrate.x            $(BINDIR)/vibrate
	mv vibrot.x             $(BINDIR)/vibrot
	mv xtalfit.x            $(BINDIR)/xtalfit
	mv xtalmin.x            $(BINDIR)/xtalmin
	mv xyzedit.x            $(BINDIR)/xyzedit
	mv xyzint.x             $(BINDIR)/xyzint
	mv xyzmol2.x            $(BINDIR)/xyzmol2
	mv xyzpdb.x             $(BINDIR)/xyzpdb

rename_exe:
	mv alchemy.x            $(BINDIR)/alchemy.exe
	mv analyze.x            $(BINDIR)/analyze.exe
	mv anneal.x             $(BINDIR)/anneal.exe
	mv arcedit.x            $(BINDIR)/arcedit.exe
	mv bar.x                $(BINDIR)/bar.exe
	mv correlate.x          $(BINDIR)/correlate.exe
	mv critical.x           $(BINDIR)/critical.exe
	mv crystal.x            $(BINDIR)/crystal.exe
	mv diffuse.x            $(BINDIR)/diffuse.exe
	mv distgeom.x           $(BINDIR)/distgeom.exe
	mv document.x           $(BINDIR)/document.exe
	mv dynamic.x            $(BINDIR)/dynamic.exe
	mv freefix.x            $(BINDIR)/freefix.exe
	mv gda.x                $(BINDIR)/gda.exe
	mv intedit.x            $(BINDIR)/intedit.exe
	mv intxyz.x             $(BINDIR)/intxyz.exe
	mv minimize.x           $(BINDIR)/minimize.exe
	mv minirot.x            $(BINDIR)/minirot.exe
	mv minrigid.x           $(BINDIR)/minrigid.exe
	mv mol2xyz.x            $(BINDIR)/mol2xyz.exe
	mv molxyz.x             $(BINDIR)/molxyz.exe
	mv monte.x              $(BINDIR)/monte.exe
	mv newton.x             $(BINDIR)/newton.exe
	mv newtrot.x            $(BINDIR)/newtrot.exe
	mv nucleic.x            $(BINDIR)/nucleic.exe
	mv optimize.x           $(BINDIR)/optimize.exe
	mv optirot.x            $(BINDIR)/optirot.exe
	mv optrigid.x           $(BINDIR)/optrigid.exe
	mv path.x               $(BINDIR)/path.exe
	mv pdbxyz.x             $(BINDIR)/pdbxyz.exe
	mv polarize.x           $(BINDIR)/polarize.exe
	mv poledit.x            $(BINDIR)/poledit.exe
	mv potential.x          $(BINDIR)/potential.exe
	mv prmedit.x            $(BINDIR)/prmedit.exe
	mv protein.x            $(BINDIR)/protein.exe
	mv pss.x                $(BINDIR)/pss.exe
	mv pssrigid.x           $(BINDIR)/pssrigid.exe
	mv pssrot.x             $(BINDIR)/pssrot.exe
	mv radial.x             $(BINDIR)/radial.exe
	mv saddle.x             $(BINDIR)/saddle.exe
	mv scan.x               $(BINDIR)/scan.exe
	mv sniffer.x            $(BINDIR)/sniffer.exe
	mv spacefill.x          $(BINDIR)/spacefill.exe
	mv spectrum.x           $(BINDIR)/spectrum.exe
	mv superpose.x          $(BINDIR)/superpose.exe
	mv testgrad.x           $(BINDIR)/testgrad.exe
	mv testhess.x           $(BINDIR)/testhess.exe
	mv testpair.x           $(BINDIR)/testpair.exe
	mv testpol.x            $(BINDIR)/testpol.exe
	mv testrot.x            $(BINDIR)/testrot.exe
	mv testvir.x            $(BINDIR)/testvir.exe
	mv timer.x              $(BINDIR)/timer.exe
	mv timerot.x            $(BINDIR)/timerot.exe
	mv torsfit.x            $(BINDIR)/torsfit.exe
	mv valence.x            $(BINDIR)/valence.exe
	mv vibbig.x             $(BINDIR)/vibbig.exe
	mv vibrate.x            $(BINDIR)/vibrate.exe
	mv vibrot.x             $(BINDIR)/vibrot.exe
	mv xtalfit.x            $(BINDIR)/xtalfit.exe
	mv xtalmin.x            $(BINDIR)/xtalmin.exe
	mv xyzedit.x            $(BINDIR)/xyzedit.exe
	mv xyzint.x             $(BINDIR)/xyzint.exe
	mv xyzmol2.x            $(BINDIR)/xyzmol2.exe
	mv xyzpdb.x             $(BINDIR)/xyzpdb.exe

remove_links:
	rm -f $(LINKDIR)/alchemy
	rm -f $(LINKDIR)/analyze
	rm -f $(LINKDIR)/anneal
	rm -f $(LINKDIR)/arcedit
	rm -f $(LINKDIR)/bar
	rm -f $(LINKDIR)/correlate
	rm -f $(LINKDIR)/critical
	rm -f $(LINKDIR)/crystal
	rm -f $(LINKDIR)/diffuse
	rm -f $(LINKDIR)/distgeom
	rm -f $(LINKDIR)/document
	rm -f $(LINKDIR)/dynamic
	rm -f $(LINKDIR)/freefix
	rm -f $(LINKDIR)/gda
	rm -f $(LINKDIR)/intedit
	rm -f $(LINKDIR)/intxyz
	rm -f $(LINKDIR)/minimize
	rm -f $(LINKDIR)/minirot
	rm -f $(LINKDIR)/minrigid
	rm -f $(LINKDIR)/mol2xyz
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
	rm -f $(LINKDIR)/testgrad
	rm -f $(LINKDIR)/testhess
	rm -f $(LINKDIR)/testpair
	rm -f $(LINKDIR)/testpol
	rm -f $(LINKDIR)/testrot
	rm -f $(LINKDIR)/testvir
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
	rm -f $(LINKDIR)/xyzmol2
	rm -f $(LINKDIR)/xyzpdb

create_links:
	ln -s $(BINDIR)/alchemy         $(LINKDIR)/alchemy
	ln -s $(BINDIR)/analyze         $(LINKDIR)/analyze
	ln -s $(BINDIR)/anneal          $(LINKDIR)/anneal
	ln -s $(BINDIR)/arcedit         $(LINKDIR)/arcedit
	ln -s $(BINDIR)/bar             $(LINKDIR)/bar
	ln -s $(BINDIR)/correlate       $(LINKDIR)/correlate
	ln -s $(BINDIR)/critical        $(LINKDIR)/critical
	ln -s $(BINDIR)/crystal         $(LINKDIR)/crystal
	ln -s $(BINDIR)/diffuse         $(LINKDIR)/diffuse
	ln -s $(BINDIR)/distgeom        $(LINKDIR)/distgeom
	ln -s $(BINDIR)/document        $(LINKDIR)/document
	ln -s $(BINDIR)/dynamic         $(LINKDIR)/dynamic
	ln -s $(BINDIR)/freefix         $(LINKDIR)/freefix
	ln -s $(BINDIR)/gda             $(LINKDIR)/gda
	ln -s $(BINDIR)/intedit         $(LINKDIR)/intedit
	ln -s $(BINDIR)/intxyz          $(LINKDIR)/intxyz
	ln -s $(BINDIR)/minimize        $(LINKDIR)/minimize
	ln -s $(BINDIR)/minirot         $(LINKDIR)/minirot
	ln -s $(BINDIR)/minrigid        $(LINKDIR)/minrigid
	ln -s $(BINDIR)/mol2xyz         $(LINKDIR)/mol2xyz
	ln -s $(BINDIR)/molxyz          $(LINKDIR)/molxyz
	ln -s $(BINDIR)/monte           $(LINKDIR)/monte
	ln -s $(BINDIR)/newton          $(LINKDIR)/newton
	ln -s $(BINDIR)/newtrot         $(LINKDIR)/newtrot
	ln -s $(BINDIR)/nucleic         $(LINKDIR)/nucleic
	ln -s $(BINDIR)/optimize        $(LINKDIR)/optimize
	ln -s $(BINDIR)/optirot         $(LINKDIR)/optirot
	ln -s $(BINDIR)/optrigid        $(LINKDIR)/optrigid
	ln -s $(BINDIR)/path            $(LINKDIR)/path
	ln -s $(BINDIR)/pdbxyz          $(LINKDIR)/pdbxyz
	ln -s $(BINDIR)/polarize        $(LINKDIR)/polarize
	ln -s $(BINDIR)/poledit         $(LINKDIR)/poledit
	ln -s $(BINDIR)/potential       $(LINKDIR)/potential
	ln -s $(BINDIR)/prmedit         $(LINKDIR)/prmedit
	ln -s $(BINDIR)/protein         $(LINKDIR)/protein
	ln -s $(BINDIR)/pss             $(LINKDIR)/pss
	ln -s $(BINDIR)/pssrigid        $(LINKDIR)/pssrigid
	ln -s $(BINDIR)/pssrot          $(LINKDIR)/pssrot
	ln -s $(BINDIR)/radial          $(LINKDIR)/radial
	ln -s $(BINDIR)/saddle          $(LINKDIR)/saddle
	ln -s $(BINDIR)/scan            $(LINKDIR)/scan
	ln -s $(BINDIR)/sniffer         $(LINKDIR)/sniffer
	ln -s $(BINDIR)/spacefill       $(LINKDIR)/spacefill
	ln -s $(BINDIR)/spectrum        $(LINKDIR)/spectrum
	ln -s $(BINDIR)/superpose       $(LINKDIR)/superpose
	ln -s $(BINDIR)/testgrad        $(LINKDIR)/testgrad
	ln -s $(BINDIR)/testhess        $(LINKDIR)/testhess
	ln -s $(BINDIR)/testpair        $(LINKDIR)/testpair
	ln -s $(BINDIR)/testpol         $(LINKDIR)/testpol
	ln -s $(BINDIR)/testrot         $(LINKDIR)/testrot
	ln -s $(BINDIR)/testvir         $(LINKDIR)/testvir
	ln -s $(BINDIR)/timer           $(LINKDIR)/timer
	ln -s $(BINDIR)/timerot         $(LINKDIR)/timerot
	ln -s $(BINDIR)/torsfit         $(LINKDIR)/torsfit
	ln -s $(BINDIR)/valence         $(LINKDIR)/valence
	ln -s $(BINDIR)/vibbig          $(LINKDIR)/vibbig
	ln -s $(BINDIR)/vibrate         $(LINKDIR)/vibrate
	ln -s $(BINDIR)/vibrot          $(LINKDIR)/vibrot
	ln -s $(BINDIR)/xtalfit         $(LINKDIR)/xtalfit
	ln -s $(BINDIR)/xtalmin         $(LINKDIR)/xtalmin
	ln -s $(BINDIR)/xyzedit         $(LINKDIR)/xyzedit
	ln -s $(BINDIR)/xyzint          $(LINKDIR)/xyzint
	ln -s $(BINDIR)/xyzmol2         $(LINKDIR)/xyzmol2
	ln -s $(BINDIR)/xyzpdb          $(LINKDIR)/xyzpdb

libtinker.a: $(LIBOBJS)
	ar $(LIBFLAGS) libtinker.a $(LIBOBJS)
	$(RANLIB) libtinker.a

###############################################################
##  Next Section has Explicit Dependencies on Include Files  ##
###############################################################

action.o:
active.o: atoms.o inform.o iounit.o keys.o usage.o
alchemy.o: analyz.o atoms.o energi.o files.o inform.o iounit.o katoms.o mutant.o potent.o units.o usage.o
align.o:
alterchg.o: angbnd.o atmlst.o atoms.o bndstr.o bound.o cflux.o charge.o chgpen.o inform.o iounit.o math.o mplpot.o mpole.o sizes.o
alterpol.o: atoms.o bound.o cell.o couple.o expol.o limits.o mpole.o neigh.o polgrp.o polpot.o shunt.o
analysis.o: analyz.o atoms.o energi.o group.o inter.o iounit.o limits.o potent.o reppot.o vdwpot.o
analyz.o:
analyze.o: action.o analyz.o angang.o angbnd.o angpot.o angtor.o atomid.o atoms.o bath.o bitor.o bndstr.o bound.o boxes.o cflux.o charge.o chgpen.o chgpot.o chgtrn.o couple.o deriv.o dipole.o disp.o energi.o ewald.o fields.o files.o improp.o imptor.o inform.o inter.o iounit.o korbs.o ktrtor.o kvdws.o limits.o math.o molcul.o moment.o mplpot.o mpole.o opbend.o opdist.o output.o piorbs.o pistuf.o pitors.o pme.o polar.o polgrp.o polpot.o potent.o repel.o solute.o strbnd.o strtor.o titles.o tors.o tortor.o units.o urey.o vdw.o vdwpot.o virial.o
angang.o:
angbnd.o:
angles.o: angbnd.o atmlst.o atoms.o couple.o iounit.o
angpot.o:
angtor.o:
anneal.o: atomid.o atoms.o bath.o bndstr.o bound.o inform.o iounit.o mdstuf.o potent.o solute.o usage.o warp.o
arcedit.o: atoms.o bound.o files.o inform.o iounit.o output.o usage.o
argue.o:
ascii.o:
atmlst.o:
atomid.o: sizes.o
atoms.o: sizes.o
attach.o: atoms.o couple.o iounit.o
baoab.o: atomid.o atoms.o bath.o freeze.o limits.o mdstuf.o moldyn.o potent.o stodyn.o units.o usage.o virial.o
bar.o: boxes.o files.o inform.o iounit.o keys.o output.o titles.o units.o
basefile.o: ascii.o files.o
bath.o:
beeman.o: atomid.o atoms.o freeze.o ielscf.o mdstuf.o moldyn.o polar.o units.o usage.o
bicubic.o:
bitor.o:
bitors.o: angbnd.o atoms.o bitor.o couple.o iounit.o
bndpot.o:
bndstr.o:
bonds.o: atmlst.o atoms.o bndstr.o couple.o iounit.o
born.o: atomid.o atoms.o bath.o chgpot.o couple.o deriv.o inform.o iounit.o math.o mpole.o pbstuf.o solpot.o solute.o virial.o
bound.o:
bounds.o: atomid.o atoms.o boxes.o math.o molcul.o
boxes.o:
bussi.o: atomid.o atoms.o bath.o boxes.o freeze.o ielscf.o mdstuf.o moldyn.o polar.o units.o usage.o
calendar.o:
cell.o:
center.o: align.o
cflux.o:
charge.o:
chgpen.o:
chgpot.o:
chgtrn.o:
chkpole.o: atoms.o mpole.o repel.o
chkring.o: couple.o
chkxyz.o: atoms.o iounit.o
cholesky.o:
chrono.o:
chunks.o:
clock.o: chrono.o
cluster.o: atomid.o atoms.o bound.o group.o inform.o iounit.o keys.o limits.o molcul.o
column.o:
command.o: argue.o
connect.o: atoms.o couple.o zclose.o zcoord.o
connolly.o: atoms.o faces.o inform.o iounit.o math.o
control.o: argue.o inform.o keys.o output.o
correlate.o: ascii.o atomid.o atoms.o files.o inform.o iounit.o
couple.o: sizes.o
critical.o: atoms.o files.o inform.o iounit.o keys.o minima.o output.o sizes.o usage.o
crystal.o: atomid.o atoms.o bound.o boxes.o couple.o files.o iounit.o math.o molcul.o
cspline.o: iounit.o
ctrpot.o:
cutoffs.o: atoms.o bound.o hescut.o keys.o limits.o neigh.o polpot.o tarray.o
damping.o: atoms.o ewald.o math.o mplpot.o mpole.o polar.o polpot.o
dcflux.o: angbnd.o atoms.o bndstr.o bound.o cflux.o sizes.o
deflate.o: iounit.o
delete.o: atomid.o atoms.o couple.o inform.o iounit.o
deriv.o:
dexpol.o: atoms.o bound.o cell.o chgpot.o couple.o deriv.o expol.o limits.o math.o mpole.o neigh.o polar.o polgrp.o polpot.o shunt.o units.o virial.o
diagq.o:
diffeq.o: atoms.o iounit.o math.o warp.o
diffuse.o: atomid.o atoms.o bound.o inform.o iounit.o molcul.o usage.o
dipole.o:
disgeo.o:
disp.o:
distgeom.o: angbnd.o atomid.o atoms.o bndstr.o couple.o disgeo.o files.o inform.o iounit.o kvdws.o math.o refer.o restrn.o tors.o
dma.o:
document.o: iounit.o
domega.o:
dsppot.o:
dynamic.o: atoms.o bath.o bndstr.o bound.o inform.o iounit.o keys.o mdstuf.o potent.o stodyn.o usage.o
eangang.o: angang.o angbnd.o angpot.o atoms.o bound.o energi.o group.o math.o usage.o
eangang1.o: angang.o angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o usage.o virial.o
eangang2.o: angang.o angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o
eangang3.o: action.o analyz.o angang.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o usage.o
eangle.o: angbnd.o angpot.o atoms.o bound.o energi.o group.o math.o usage.o
eangle1.o: angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o usage.o virial.o
eangle2.o: angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o
eangle3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o usage.o
eangtor.o: angbnd.o angtor.o atoms.o bound.o energi.o group.o math.o torpot.o tors.o usage.o
eangtor1.o: angbnd.o angtor.o atoms.o bound.o deriv.o energi.o group.o math.o torpot.o tors.o usage.o virial.o
eangtor2.o: angbnd.o angtor.o atoms.o bound.o group.o hessn.o math.o torpot.o tors.o
eangtor3.o: action.o analyz.o angbnd.o angtor.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o torpot.o tors.o usage.o
ebond.o: atoms.o bndpot.o bndstr.o bound.o energi.o group.o usage.o
ebond1.o: atoms.o bndpot.o bndstr.o bound.o deriv.o energi.o group.o usage.o virial.o
ebond2.o: atmlst.o atoms.o bndpot.o bndstr.o bound.o couple.o group.o hessn.o
ebond3.o: action.o analyz.o atomid.o atoms.o bndpot.o bndstr.o bound.o energi.o group.o inform.o iounit.o usage.o
ebuck.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o iounit.o light.o limits.o math.o neigh.o shunt.o usage.o vdw.o vdwpot.o warp.o
ebuck1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o iounit.o light.o limits.o math.o neigh.o shunt.o usage.o vdw.o vdwpot.o virial.o warp.o
ebuck2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o iounit.o math.o shunt.o vdw.o vdwpot.o warp.o
ebuck3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o usage.o vdw.o vdwpot.o warp.o
echarge.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o energi.o ewald.o extfld.o group.o iounit.o light.o limits.o math.o neigh.o pme.o shunt.o usage.o warp.o
echarge1.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o deriv.o energi.o ewald.o extfld.o group.o light.o limits.o math.o neigh.o pme.o shunt.o usage.o virial.o warp.o
echarge2.o: atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o deriv.o ewald.o group.o hessn.o limits.o math.o neigh.o pme.o shunt.o warp.o
echarge3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o charge.o chgpot.o couple.o energi.o ewald.o extfld.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o pme.o shunt.o usage.o warp.o
echgdpl.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o energi.o group.o shunt.o units.o usage.o
echgdpl1.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o deriv.o dipole.o energi.o group.o shunt.o units.o usage.o virial.o
echgdpl2.o: atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o group.o hessn.o shunt.o units.o
echgdpl3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o charge.o chgpot.o couple.o dipole.o energi.o group.o inform.o inter.o iounit.o molcul.o shunt.o units.o usage.o
echgtrn.o: atoms.o bound.o boxes.o cell.o chgpot.o chgtrn.o couple.o ctrpot.o energi.o group.o light.o limits.o mplpot.o mpole.o mutant.o neigh.o shunt.o usage.o
echgtrn1.o: atoms.o bound.o cell.o chgpot.o chgtrn.o couple.o ctrpot.o deriv.o energi.o group.o limits.o mplpot.o mpole.o mutant.o neigh.o shunt.o usage.o virial.o
echgtrn2.o: atoms.o bound.o cell.o chgpot.o chgtrn.o couple.o ctrpot.o group.o hessn.o mplpot.o mpole.o mutant.o shunt.o usage.o
echgtrn3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o chgpot.o chgtrn.o couple.o ctrpot.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o mplpot.o mpole.o mutant.o neigh.o shunt.o usage.o
edipole.o: atoms.o bound.o cell.o chgpot.o dipole.o energi.o group.o shunt.o units.o usage.o
edipole1.o: atoms.o bound.o cell.o chgpot.o deriv.o dipole.o energi.o group.o shunt.o units.o usage.o virial.o
edipole2.o: atoms.o bound.o cell.o chgpot.o dipole.o group.o hessn.o shunt.o units.o
edipole3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o chgpot.o dipole.o energi.o group.o inform.o inter.o iounit.o molcul.o shunt.o units.o usage.o
edisp.o: atoms.o bound.o boxes.o cell.o couple.o disp.o dsppot.o energi.o ewald.o group.o limits.o math.o mutant.o neigh.o pme.o shunt.o usage.o
edisp1.o: atoms.o bound.o boxes.o cell.o couple.o deriv.o disp.o dsppot.o energi.o ewald.o group.o limits.o math.o mutant.o neigh.o pme.o shunt.o usage.o virial.o
edisp2.o: atoms.o bound.o cell.o couple.o disp.o dsppot.o group.o hessn.o shunt.o usage.o
edisp3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o disp.o dsppot.o energi.o ewald.o group.o inform.o inter.o iounit.o limits.o molcul.o mutant.o neigh.o pme.o shunt.o usage.o
egauss.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o math.o neigh.o shunt.o usage.o vdw.o vdwpot.o warp.o
egauss1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o light.o limits.o math.o neigh.o shunt.o usage.o vdw.o vdwpot.o virial.o warp.o
egauss2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o shunt.o vdw.o vdwpot.o warp.o
egauss3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o neigh.o shunt.o usage.o vdw.o vdwpot.o warp.o
egeom.o: atomid.o atoms.o bound.o boxes.o cell.o energi.o group.o math.o molcul.o restrn.o usage.o
egeom1.o: atomid.o atoms.o bound.o boxes.o cell.o deriv.o energi.o group.o math.o molcul.o restrn.o usage.o virial.o
egeom2.o: atomid.o atoms.o bound.o boxes.o cell.o deriv.o group.o hessn.o math.o molcul.o restrn.o
egeom3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o energi.o group.o inform.o inter.o iounit.o math.o molcul.o restrn.o usage.o
ehal.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o mutant.o neigh.o shunt.o usage.o vdw.o vdwpot.o
ehal1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o light.o limits.o mutant.o neigh.o shunt.o usage.o vdw.o vdwpot.o virial.o
ehal2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o shunt.o vdw.o vdwpot.o
ehal3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o mutant.o neigh.o shunt.o usage.o vdw.o vdwpot.o
eimprop.o: atoms.o bound.o energi.o group.o improp.o math.o torpot.o usage.o
eimprop1.o: atoms.o bound.o deriv.o energi.o group.o improp.o math.o torpot.o usage.o virial.o
eimprop2.o: atoms.o bound.o group.o hessn.o improp.o math.o torpot.o
eimprop3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o improp.o inform.o iounit.o math.o torpot.o usage.o
eimptor.o: atoms.o bound.o energi.o group.o imptor.o torpot.o usage.o
eimptor1.o: atoms.o bound.o deriv.o energi.o group.o imptor.o torpot.o usage.o virial.o
eimptor2.o: atoms.o bound.o group.o hessn.o imptor.o torpot.o
eimptor3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o imptor.o inform.o iounit.o math.o torpot.o usage.o
elj.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o light.o limits.o math.o mutant.o neigh.o shunt.o usage.o vdw.o vdwpot.o warp.o
elj1.o: atomid.o atoms.o bound.o boxes.o cell.o couple.o deriv.o energi.o group.o light.o limits.o math.o mutant.o neigh.o shunt.o usage.o vdw.o vdwpot.o virial.o warp.o
elj2.o: atomid.o atoms.o bound.o cell.o couple.o group.o hessn.o math.o shunt.o vdw.o vdwpot.o warp.o
elj3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o math.o molcul.o mutant.o neigh.o shunt.o usage.o vdw.o vdwpot.o warp.o
embed.o: angbnd.o atoms.o bndstr.o couple.o disgeo.o files.o inform.o iounit.o keys.o light.o math.o minima.o output.o refer.o restrn.o tors.o units.o
emetal.o: atomid.o atoms.o couple.o energi.o kchrge.o
emetal1.o: atomid.o atoms.o couple.o deriv.o energi.o kchrge.o
emetal2.o:
emetal3.o: action.o analyz.o atomid.o atoms.o energi.o kchrge.o
emm3hb.o: atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o energi.o group.o light.o limits.o neigh.o shunt.o usage.o vdw.o vdwpot.o
emm3hb1.o: atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o deriv.o energi.o group.o light.o limits.o neigh.o shunt.o usage.o vdw.o vdwpot.o virial.o
emm3hb2.o: atmlst.o atomid.o atoms.o bndstr.o bound.o cell.o chgpot.o couple.o group.o hessn.o shunt.o vdw.o vdwpot.o
emm3hb3.o: action.o analyz.o atmlst.o atomid.o atoms.o bndstr.o bound.o boxes.o cell.o chgpot.o couple.o energi.o group.o inform.o inter.o iounit.o light.o limits.o molcul.o neigh.o shunt.o usage.o vdw.o vdwpot.o
empole.o: atoms.o bound.o boxes.o cell.o chgpen.o chgpot.o couple.o energi.o ewald.o extfld.o group.o limits.o math.o mplpot.o mpole.o mrecip.o neigh.o pme.o polpot.o shunt.o usage.o
empole1.o: atoms.o bound.o boxes.o cell.o chgpen.o chgpot.o couple.o deriv.o energi.o ewald.o extfld.o group.o limits.o math.o mplpot.o mpole.o mrecip.o neigh.o pme.o potent.o shunt.o usage.o virial.o
empole2.o: atoms.o bound.o boxes.o cell.o chgpen.o chgpot.o couple.o deriv.o energi.o group.o hessn.o limits.o molcul.o mplpot.o mpole.o potent.o shunt.o usage.o
empole3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o chgpen.o chgpot.o couple.o energi.o ewald.o extfld.o group.o inform.o inter.o iounit.o limits.o math.o molcul.o mplpot.o mpole.o mrecip.o neigh.o pme.o potent.o shunt.o usage.o
energi.o:
energy.o: energi.o iounit.o limits.o potent.o rigid.o vdwpot.o
eopbend.o: angbnd.o angpot.o atoms.o bound.o energi.o fields.o group.o math.o opbend.o usage.o
eopbend1.o: angbnd.o angpot.o atoms.o bound.o deriv.o energi.o group.o math.o opbend.o usage.o virial.o
eopbend2.o: angbnd.o angpot.o atoms.o bound.o group.o hessn.o math.o opbend.o
eopbend3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o opbend.o usage.o
eopdist.o: angpot.o atoms.o bound.o energi.o group.o opdist.o usage.o
eopdist1.o: angpot.o atoms.o bound.o deriv.o energi.o group.o opdist.o usage.o virial.o
eopdist2.o: angpot.o atoms.o bound.o group.o hessn.o opdist.o usage.o
eopdist3.o: action.o analyz.o angpot.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o opdist.o usage.o
epitors.o: atoms.o bound.o energi.o group.o pitors.o torpot.o usage.o
epitors1.o: atoms.o bound.o deriv.o energi.o group.o pitors.o torpot.o usage.o virial.o
epitors2.o: angbnd.o atoms.o bound.o deriv.o group.o hessn.o pitors.o torpot.o usage.o
epitors3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o pitors.o torpot.o usage.o
epolar.o: atoms.o bound.o boxes.o cell.o chgpen.o chgpot.o couple.o energi.o ewald.o extfld.o limits.o math.o mplpot.o mpole.o mrecip.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o
epolar1.o: atoms.o bound.o boxes.o cell.o chgpen.o chgpot.o couple.o deriv.o energi.o ewald.o iounit.o limits.o math.o molcul.o mplpot.o mpole.o mrecip.o neigh.o pme.o polar.o polgrp.o polopt.o polpot.o poltcg.o potent.o shunt.o virial.o
epolar2.o: atoms.o bound.o cell.o chgpen.o chgpot.o couple.o deriv.o hessn.o limits.o mplpot.o mpole.o polar.o polgrp.o polopt.o polpot.o poltcg.o potent.o shunt.o
epolar3.o: action.o analyz.o atomid.o atoms.o bound.o boxes.o cell.o chgpen.o chgpot.o couple.o energi.o ewald.o extfld.o inform.o inter.o iounit.o limits.o math.o molcul.o mplpot.o mpole.o mrecip.o neigh.o pme.o polar.o polgrp.o polpot.o potent.o shunt.o units.o
erepel.o: atoms.o bound.o cell.o couple.o energi.o group.o inform.o inter.o limits.o mpole.o mutant.o neigh.o polar.o repel.o reppot.o shunt.o usage.o
erepel1.o: atoms.o bound.o cell.o couple.o deriv.o energi.o group.o inform.o limits.o mpole.o mutant.o neigh.o repel.o reppot.o shunt.o usage.o virial.o
erepel2.o: atoms.o bound.o cell.o couple.o deriv.o group.o hessn.o mpole.o potent.o repel.o reppot.o shunt.o usage.o
erepel3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o limits.o molcul.o mpole.o mutant.o neigh.o polar.o repel.o reppot.o shunt.o usage.o
erf.o: iounit.o math.o
erxnfld.o: atoms.o chgpot.o energi.o mpole.o rxnfld.o rxnpot.o shunt.o usage.o
erxnfld1.o: atoms.o deriv.o energi.o
erxnfld2.o:
erxnfld3.o: action.o analyz.o atomid.o atoms.o chgpot.o energi.o inform.o iounit.o mpole.o shunt.o usage.o
esolv.o: atomid.o atoms.o bound.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o kvdws.o limits.o math.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o solpot.o solute.o usage.o vdw.o warp.o
esolv1.o: atomid.o atoms.o bound.o boxes.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o kvdws.o limits.o math.o mplpot.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o solpot.o solute.o usage.o vdw.o virial.o warp.o
esolv2.o: atoms.o charge.o chgpot.o deriv.o hessn.o math.o mpole.o potent.o shunt.o solpot.o solute.o warp.o
esolv3.o: action.o analyz.o atomid.o atoms.o bound.o charge.o chgpot.o couple.o deriv.o energi.o gkstuf.o group.o hpmf.o inform.o inter.o iounit.o kvdws.o limits.o math.o molcul.o mpole.o neigh.o nonpol.o pbstuf.o polar.o polgrp.o polpot.o potent.o shunt.o solpot.o solute.o usage.o vdw.o warp.o
estrbnd.o: angbnd.o angpot.o atoms.o bndstr.o bound.o energi.o group.o math.o strbnd.o usage.o
estrbnd1.o: angbnd.o angpot.o atoms.o bndstr.o bound.o deriv.o energi.o group.o math.o strbnd.o usage.o virial.o
estrbnd2.o: angbnd.o angpot.o atoms.o bndstr.o bound.o group.o hessn.o math.o strbnd.o
estrbnd3.o: action.o analyz.o angbnd.o angpot.o atomid.o atoms.o bndstr.o bound.o energi.o group.o inform.o iounit.o math.o strbnd.o usage.o
estrtor.o: atoms.o bndstr.o bound.o energi.o group.o strtor.o torpot.o tors.o usage.o
estrtor1.o: atoms.o bndstr.o bound.o deriv.o energi.o group.o strtor.o torpot.o tors.o usage.o virial.o
estrtor2.o: atoms.o bndstr.o bound.o group.o hessn.o strtor.o torpot.o tors.o
estrtor3.o: action.o analyz.o atomid.o atoms.o bndstr.o bound.o energi.o group.o inform.o iounit.o math.o strtor.o torpot.o tors.o usage.o
etors.o: atoms.o bound.o energi.o group.o math.o torpot.o tors.o usage.o warp.o
etors1.o: atoms.o bound.o deriv.o energi.o group.o math.o torpot.o tors.o usage.o virial.o warp.o
etors2.o: atoms.o bound.o group.o hessn.o math.o torpot.o tors.o warp.o
etors3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o math.o torpot.o tors.o usage.o warp.o
etortor.o: atomid.o atoms.o bitor.o bound.o couple.o energi.o group.o ktrtor.o math.o torpot.o tortor.o usage.o
etortor1.o: atoms.o bitor.o bound.o deriv.o energi.o group.o ktrtor.o math.o torpot.o tortor.o usage.o virial.o
etortor2.o: atoms.o bitor.o bound.o group.o hessn.o ktrtor.o math.o torpot.o tortor.o units.o
etortor3.o: action.o analyz.o atoms.o bitor.o bound.o energi.o group.o inform.o iounit.o ktrtor.o math.o torpot.o tortor.o usage.o
eurey.o: atoms.o bound.o energi.o group.o urey.o urypot.o usage.o
eurey1.o: atoms.o bound.o deriv.o energi.o group.o urey.o urypot.o usage.o virial.o
eurey2.o: atoms.o bound.o couple.o group.o hessn.o urey.o urypot.o
eurey3.o: action.o analyz.o atomid.o atoms.o bound.o energi.o group.o inform.o iounit.o urey.o urypot.o usage.o
evcorr.o: atomid.o atoms.o bound.o boxes.o kdsp.o limits.o math.o mutant.o potent.o shunt.o vdw.o vdwpot.o
ewald.o:
exfield.o: action.o analyz.o atoms.o charge.o chgpot.o deriv.o energi.o extfld.o mpole.o usage.o virial.o
expol.o:
exrepel1.o: atoms.o bound.o cell.o couple.o deriv.o energi.o group.o limits.o mpole.o mutant.o repel.o reppot.o shunt.o units.o usage.o virial.o xrepel.o
exrepel3.o: action.o analyz.o atomid.o atoms.o bound.o cell.o couple.o energi.o group.o inform.o inter.o iounit.o limits.o molcul.o mutant.o repel.o reppot.o shunt.o units.o usage.o xrepel.o
extfld.o:
extra.o: energi.o
extra1.o: atoms.o deriv.o energi.o
extra2.o: atoms.o hessn.o
extra3.o: action.o analyz.o atoms.o energi.o
faces.o:
fatal.o: iounit.o
fft.o:
fft3d.o: fft.o openmp.o pme.o
fftpack.o: math.o
field.o: fields.o inform.o iounit.o keys.o potent.o sizes.o
fields.o:
files.o:
final.o: align.o analyz.o angang.o angbnd.o angtor.o atmlst.o bitor.o bndstr.o cell.o cflux.o charge.o chgpen.o chgtrn.o chunks.o couple.o deriv.o dipole.o disgeo.o domega.o expol.o faces.o fft.o fields.o fracs.o freeze.o group.o hessn.o hpmf.o ielscf.o improp.o imptor.o inform.o iounit.o kanang.o kangs.o kantor.o katoms.o kbonds.o kcflux.o kchrge.o kcpen.o kctrn.o kdipol.o kdsp.o kexpl.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolpr.o kpolr.o krepl.o ksolut.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o kxrepl.o light.o limits.o merck.o molcul.o moldyn.o mpole.o mrecip.o mutant.o neigh.o nonpol.o omega.o opbend.o opdist.o orbits.o paths.o pbstuf.o pdb.o piorbs.o pistuf.o pitors.o pme.o polar.o polgrp.o polopt.o polpcg.o poltcg.o potfit.o qmstuf.o refer.o repel.o restrn.o rgddyn.o rigid.o ring.o rotbnd.o socket.o solpot.o solute.o stodyn.o strbnd.o strtor.o syntrn.o tarray.o tors.o tortor.o uprior.o urey.o usage.o vdw.o vibs.o warp.o xrepel.o
flatten.o: atoms.o fields.o inform.o iounit.o keys.o warp.o
fracs.o:
freefix.o: iounit.o math.o units.o
freeunit.o: iounit.o
freeze.o:
gda.o: atoms.o files.o iounit.o minima.o potent.o vdwpot.o warp.o
geometry.o: atoms.o math.o
getarc.o: files.o inform.o iounit.o output.o
getcart.o: files.o inform.o iounit.o output.o
getdcd.o: files.o inform.o iounit.o output.o
getint.o: atoms.o files.o inform.o iounit.o output.o
getkey.o: argue.o files.o iounit.o keys.o openmp.o
getmol.o: files.o inform.o iounit.o
getmol2.o: files.o inform.o iounit.o
getnumb.o: ascii.o
getpdb.o: files.o inform.o iounit.o
getprm.o: files.o inform.o iounit.o keys.o params.o
getref.o: atomid.o atoms.o boxes.o couple.o files.o refer.o titles.o
getstring.o: ascii.o
gettext.o: ascii.o
getword.o: ascii.o
getxyz.o: files.o inform.o iounit.o output.o
ghmcstep.o: atomid.o atoms.o bath.o freeze.o iounit.o moldyn.o stodyn.o units.o usage.o virial.o
gkstuf.o: sizes.o
gradient.o: atoms.o couple.o deriv.o energi.o inform.o inter.o iounit.o limits.o potent.o rigid.o vdwpot.o virial.o
gradrgd.o: atoms.o group.o rigid.o
gradrot.o: atoms.o deriv.o domega.o omega.o potent.o rotbnd.o
group.o:
groups.o: group.o
grpline.o: atomid.o atoms.o group.o rgddyn.o
gyrate.o: atoms.o usage.o
hescut.o:
hessian.o: atoms.o couple.o hescut.o hessn.o inform.o iounit.o limits.o mpole.o potent.o rigid.o usage.o vdw.o vdwpot.o
hessn.o:
hessrgd.o: atoms.o group.o rigid.o
hessrot.o: math.o omega.o zcoord.o
hpmf.o:
hybrid.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o charge.o couple.o dipole.o imptor.o inform.o iounit.o kangs.o katoms.o kbonds.o kchrge.o kdipol.o kitors.o kstbnd.o ksttor.o ktorsn.o kvdws.o math.o mutant.o strbnd.o strtor.o tors.o vdw.o vdwpot.o
ielscf.o:
image.o: boxes.o cell.o math.o
impose.o: align.o inform.o iounit.o
improp.o:
imptor.o:
induce.o: atoms.o bound.o boxes.o cell.o chgpen.o couple.o ewald.o expol.o extfld.o gkstuf.o group.o ielscf.o inform.o iounit.o limits.o math.o mplpot.o mpole.o neigh.o openmp.o pbstuf.o pme.o polar.o polgrp.o polopt.o polpcg.o polpot.o potent.o shunt.o solpot.o solute.o tarray.o units.o uprior.o
inertia.o: atomid.o atoms.o iounit.o math.o
inform.o:
initatom.o: ptable.o
initial.o: align.o atoms.o bath.o bound.o boxes.o cell.o fft.o files.o group.o inform.o iounit.o keys.o linmin.o minima.o molcul.o mutant.o neigh.o openmp.o output.o params.o pdb.o rigid.o scales.o sequen.o socket.o virial.o warp.o zclose.o
initprm.o: angpot.o bndpot.o chgpot.o ctrpot.o dsppot.o expol.o extfld.o fields.o ielscf.o kanang.o kangs.o kantor.o katoms.o kbonds.o kcflux.o kchrge.o kcpen.o kctrn.o kdipol.o kdsp.o kexpl.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolpr.o kpolr.o krepl.o ksolut.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o kxrepl.o math.o merck.o mplpot.o polpot.o reppot.o rxnpot.o sizes.o solpot.o solute.o torpot.o units.o uprior.o urypot.o vdwpot.o xrepel.o
initres.o: resdue.o
initrot.o: atoms.o couple.o group.o inform.o iounit.o math.o omega.o potent.o restrn.o rotbnd.o usage.o zcoord.o
insert.o: atomid.o atoms.o couple.o inform.o iounit.o
intedit.o: atomid.o atoms.o files.o iounit.o katoms.o zcoord.o
inter.o:
intxyz.o: files.o iounit.o titles.o
invbeta.o:
invert.o: iounit.o
iounit.o:
jacobi.o: iounit.o
kanang.o:
kangang.o: angang.o angbnd.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o kanang.o keys.o potent.o tors.o
kangle.o: angbnd.o angpot.o atomid.o atoms.o bndstr.o couple.o fields.o inform.o iounit.o kangs.o keys.o merck.o potent.o ring.o usage.o
kangs.o:
kangtor.o: angtor.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o kantor.o keys.o potent.o tors.o
kantor.o:
katom.o: atomid.o atoms.o couple.o inform.o iounit.o katoms.o keys.o
katoms.o:
kbond.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o couple.o fields.o inform.o iounit.o kbonds.o keys.o merck.o potent.o tors.o usage.o
kbonds.o:
kcflux.o:
kcharge.o: atomid.o atoms.o charge.o chgpot.o couple.o fields.o inform.o iounit.o kchrge.o keys.o merck.o potent.o
kchgflx.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o cflux.o couple.o inform.o iounit.o kangs.o kbonds.o kcflux.o keys.o potent.o sizes.o usage.o
kchgtrn.o: atomid.o atoms.o chgpen.o chgtrn.o expol.o inform.o iounit.o kctrn.o keys.o mplpot.o mpole.o polar.o polpot.o potent.o sizes.o
kchrge.o:
kcpen.o:
kctrn.o:
kdipol.o:
kdipole.o: atmlst.o atoms.o bndstr.o couple.o dipole.o inform.o iounit.o kdipol.o keys.o potent.o
kdisp.o: atomid.o atoms.o disp.o dsppot.o inform.o iounit.o kdsp.o keys.o limits.o potent.o sizes.o
kdsp.o:
kewald.o: atoms.o bound.o boxes.o chunks.o ewald.o fft.o inform.o iounit.o keys.o limits.o openmp.o pme.o potent.o
kexpl.o:
kexpol.o: atomid.o atoms.o expol.o inform.o iounit.o kexpl.o keys.o sizes.o
kextra.o:
keys.o:
kgeom.o: atomid.o atoms.o bound.o couple.o group.o iounit.o keys.o molcul.o potent.o restrn.o
khbond.o:
kimprop.o: atomid.o atoms.o couple.o improp.o inform.o iounit.o keys.o kiprop.o potent.o tors.o
kimptor.o: atomid.o atoms.o couple.o imptor.o inform.o iounit.o keys.o kitors.o math.o potent.o tors.o
kinetic.o: atomid.o atoms.o bath.o group.o ielscf.o mdstuf.o moldyn.o rgddyn.o units.o usage.o
kiprop.o:
kitors.o:
kmetal.o:
kmpole.o: atomid.o atoms.o chgpen.o couple.o inform.o iounit.o kcpen.o keys.o kmulti.o math.o mplpot.o mpole.o polar.o polgrp.o potent.o units.o
kmulti.o:
kopbend.o: angbnd.o atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o kopbnd.o merck.o opbend.o potent.o usage.o
kopbnd.o:
kopdist.o: angbnd.o atmlst.o atomid.o atoms.o couple.o inform.o iounit.o keys.o kopdst.o opdist.o potent.o
kopdst.o:
korbit.o: atomid.o atoms.o bndstr.o inform.o iounit.o keys.o korbs.o orbits.o piorbs.o pistuf.o tors.o units.o
korbs.o:
kpitor.o:
kpitors.o: atomid.o atoms.o bndstr.o couple.o inform.o iounit.o keys.o kpitor.o pitors.o potent.o tors.o
kpolar.o: atoms.o chgpen.o couple.o expol.o inform.o iounit.o keys.o kpolpr.o kpolr.o mplpot.o mpole.o polar.o polgrp.o polopt.o polpcg.o polpot.o poltcg.o potent.o
kpolpr.o:
kpolr.o:
krepel.o: atomid.o atoms.o inform.o iounit.o keys.o krepl.o mpole.o potent.o repel.o reppot.o sizes.o
krepl.o:
ksolut.o:
ksolv.o: angbnd.o atmlst.o atomid.o atoms.o bath.o bndstr.o chgpot.o couple.o gkstuf.o hpmf.o inform.o iounit.o keys.o ksolut.o kvdws.o math.o nonpol.o pbstuf.o polar.o polopt.o polpot.o potent.o ptable.o sizes.o solpot.o solute.o vdw.o
kstbnd.o:
kstrbnd.o: angbnd.o angpot.o atmlst.o atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o kstbnd.o merck.o potent.o ring.o strbnd.o
kstrtor.o: atmlst.o atomid.o atoms.o couple.o inform.o iounit.o keys.o ksttor.o potent.o strtor.o tors.o
ksttor.o:
ktors.o: atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o ktorsn.o math.o merck.o potent.o ring.o tors.o usage.o
ktorsn.o:
ktortor.o: atomid.o atoms.o bitor.o inform.o iounit.o keys.o ktrtor.o potent.o tortor.o
ktrtor.o:
kurey.o: angbnd.o atomid.o atoms.o inform.o iounit.o keys.o kurybr.o potent.o urey.o
kurybr.o:
kvdw.o: atomid.o atoms.o couple.o fields.o inform.o iounit.o keys.o khbond.o kvdwpr.o kvdws.o math.o merck.o potent.o vdw.o vdwpot.o
kvdwpr.o:
kvdws.o:
kxrepel.o: atomid.o atoms.o inform.o iounit.o keys.o kxrepl.o mpole.o potent.o repel.o reppot.o sizes.o xrepel.o
kxrepl.o:
lattice.o: bound.o boxes.o cell.o inform.o iounit.o math.o
lbfgs.o: inform.o iounit.o keys.o linmin.o math.o minima.o output.o scales.o
light.o:
lights.o: bound.o boxes.o cell.o iounit.o light.o
limits.o:
linmin.o:
lusolve.o: iounit.o
makeint.o: atoms.o couple.o inform.o iounit.o math.o sizes.o zclose.o zcoord.o
makeref.o: atomid.o atoms.o boxes.o couple.o files.o refer.o titles.o
makexyz.o: atoms.o zcoord.o
math.o:
maxwell.o: units.o
mdinit.o: atomid.o atoms.o bath.o bound.o couple.o extfld.o files.o freeze.o group.o inform.o iounit.o keys.o mdstuf.o molcul.o moldyn.o mpole.o output.o potent.o rgddyn.o rigid.o stodyn.o units.o usage.o
mdrest.o: atomid.o atoms.o bound.o group.o inform.o iounit.o mdstuf.o moldyn.o rgddyn.o units.o
mdsave.o: atomid.o atoms.o bound.o boxes.o deriv.o files.o group.o inform.o iounit.o mdstuf.o moldyn.o mpole.o output.o polar.o potent.o rgddyn.o socket.o titles.o units.o
mdstat.o: atoms.o bath.o bound.o boxes.o inform.o iounit.o limits.o mdstuf.o molcul.o units.o usage.o warp.o
mdstuf.o:
mechanic.o: inform.o iounit.o
merck.o: sizes.o
merge.o: atomid.o atoms.o couple.o iounit.o refer.o
minima.o:
minimize.o: atoms.o bound.o files.o freeze.o inform.o iounit.o scales.o usage.o
minirot.o: files.o inform.o iounit.o math.o omega.o scales.o zcoord.o
minrigid.o: files.o group.o inform.o iounit.o math.o output.o rigid.o
mol2xyz.o: files.o iounit.o titles.o
molcul.o:
moldyn.o:
molecule.o: atomid.o atoms.o couple.o molcul.o
molxyz.o: files.o iounit.o titles.o
moment.o:
moments.o: atomid.o atoms.o bound.o charge.o dipole.o limits.o moment.o mpole.o polar.o potent.o rigid.o solpot.o units.o usage.o
monte.o: atoms.o bound.o files.o inform.o iounit.o omega.o output.o potent.o units.o usage.o zcoord.o
mplpot.o:
mpole.o:
mrecip.o:
mutant.o:
mutate.o: angbnd.o atomid.o atoms.o bndstr.o cflux.o charge.o chgpen.o dipole.o inform.o iounit.o katoms.o keys.o mplpot.o mpole.o mutant.o polar.o potent.o tors.o
nblist.o: atoms.o bound.o boxes.o cell.o charge.o disp.o iounit.o light.o limits.o mpole.o neigh.o potent.o vdw.o
neigh.o:
newton.o: atoms.o bound.o files.o inform.o iounit.o usage.o
newtrot.o: files.o hescut.o inform.o iounit.o math.o omega.o zcoord.o
nextarg.o: argue.o
nexttext.o:
nonpol.o:
nose.o: atomid.o atoms.o bath.o boxes.o freeze.o mdstuf.o moldyn.o units.o usage.o virial.o
nspline.o:
nucleic.o: atoms.o couple.o files.o group.o inform.o iounit.o katoms.o math.o molcul.o nucleo.o output.o potent.o resdue.o restrn.o rigid.o sequen.o titles.o usage.o
nucleo.o: sizes.o
number.o: inform.o iounit.o
numeral.o:
numgrad.o: atoms.o
ocvm.o: inform.o iounit.o keys.o linmin.o math.o minima.o output.o potent.o scales.o
omega.o:
opbend.o:
opdist.o:
openend.o:
openmp.o:
optimize.o: atoms.o bound.o files.o inform.o iounit.o scales.o usage.o
optinit.o: inform.o keys.o output.o potent.o
optirot.o: files.o inform.o iounit.o math.o omega.o scales.o zcoord.o
optrigid.o: files.o group.o inform.o iounit.o math.o output.o rigid.o
optsave.o: atomid.o atoms.o bound.o deriv.o files.o group.o iounit.o math.o mpole.o omega.o output.o polar.o potent.o scales.o socket.o titles.o units.o usage.o zcoord.o
orbital.o: atomid.o atoms.o bndstr.o couple.o iounit.o keys.o piorbs.o potent.o tors.o
orbits.o:
orient.o: atomid.o atoms.o group.o math.o rigid.o
orthog.o:
output.o:
overlap.o: units.o
params.o:
path.o: align.o atomid.o atoms.o files.o inform.o iounit.o linmin.o minima.o output.o paths.o
paths.o:
pbstuf.o:
pdb.o:
pdbxyz.o: atomid.o atoms.o couple.o fields.o files.o inform.o iounit.o katoms.o pdb.o resdue.o sequen.o titles.o
phipsi.o: sizes.o
picalc.o: atomid.o atoms.o bndstr.o couple.o inform.o iounit.o orbits.o piorbs.o pistuf.o tors.o units.o
piorbs.o:
pistuf.o:
pitors.o:
pme.o:
pmestuf.o: atoms.o boxes.o charge.o chunks.o disp.o math.o mpole.o openmp.o pme.o potent.o
pmpb.o: iounit.o
polar.o:
polarize.o: atoms.o inform.o iounit.o molcul.o mpole.o polar.o polopt.o polpcg.o polpot.o potent.o units.o
poledit.o: atomid.o atoms.o bndstr.o chgpen.o couple.o fields.o files.o inform.o iounit.o keys.o kpolr.o mplpot.o mpole.o polar.o polgrp.o polpot.o potent.o ptable.o ring.o sizes.o units.o
polgrp.o:
polopt.o:
polpcg.o:
polpot.o:
poltcg.o:
polymer.o: atoms.o bndstr.o bound.o boxes.o iounit.o keys.o
potent.o:
potential.o: atomid.o atoms.o charge.o chgpen.o chgpot.o dipole.o files.o inform.o iounit.o katoms.o keys.o math.o moment.o mplpot.o mpole.o neigh.o polar.o potent.o potfit.o ptable.o refer.o titles.o units.o
potfit.o: sizes.o
predict.o: atomid.o atoms.o ielscf.o keys.o polar.o uprior.o
pressure.o: atomid.o atoms.o bath.o bound.o boxes.o group.o math.o mdstuf.o molcul.o moldyn.o units.o usage.o virial.o
prmedit.o: angpot.o bndpot.o iounit.o math.o params.o sizes.o urypot.o vdwpot.o
prmkey.o: angpot.o bndpot.o chgpot.o ctrpot.o dsppot.o expol.o extfld.o fields.o mplpot.o polpot.o potent.o reppot.o rxnpot.o torpot.o units.o urypot.o vdwpot.o xrepel.o
promo.o: iounit.o
protein.o: atomid.o atoms.o couple.o files.o group.o inform.o iounit.o katoms.o math.o molcul.o output.o phipsi.o potent.o resdue.o restrn.o rigid.o sequen.o titles.o usage.o
prtarc.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o output.o titles.o usage.o
prtdcd.o: atoms.o bound.o boxes.o files.o titles.o
prtdyn.o: atoms.o boxes.o files.o group.o mdstuf.o moldyn.o rgddyn.o titles.o
prterr.o: files.o output.o
prtint.o: atomid.o atoms.o files.o inform.o titles.o zclose.o zcoord.o
prtmol2.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o couple.o files.o iounit.o ptable.o ring.o titles.o tors.o
prtpdb.o: bound.o boxes.o files.o pdb.o sequen.o titles.o
prtprm.o: angpot.o bndpot.o chgpot.o fields.o kanang.o kangs.o kantor.o katoms.o kbonds.o kcflux.o kchrge.o kcpen.o kctrn.o kdipol.o kdsp.o kexpl.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolpr.o kpolr.o krepl.o ksolut.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o mplpot.o polpot.o sizes.o urypot.o vdwpot.o
prtseq.o: files.o sequen.o
prtxyz.o: atomid.o atoms.o bound.o boxes.o couple.o files.o inform.o titles.o
pss.o: atoms.o files.o hescut.o inform.o iounit.o math.o omega.o refer.o tree.o warp.o zcoord.o
pssrigid.o: atoms.o files.o group.o inform.o iounit.o math.o minima.o molcul.o refer.o rigid.o sizes.o warp.o
pssrot.o: atoms.o files.o inform.o iounit.o math.o minima.o omega.o refer.o warp.o zcoord.o
ptable.o:
qmstuf.o:
qrsolve.o:
quatfit.o: align.o
radial.o: argue.o atomid.o atoms.o bound.o boxes.o files.o inform.o iounit.o limits.o math.o molcul.o potent.o
random.o: inform.o iounit.o keys.o math.o
rattle.o: atomid.o atoms.o freeze.o group.o inform.o iounit.o moldyn.o units.o usage.o virial.o
readcart.o: output.o
readdcd.o: atoms.o bound.o boxes.o files.o inform.o iounit.o math.o titles.o
readdyn.o: atoms.o boxes.o files.o group.o iounit.o mdstuf.o moldyn.o rgddyn.o
readgau.o: ascii.o iounit.o qmstuf.o units.o
readgdma.o: atomid.o atoms.o dma.o files.o iounit.o mpole.o units.o
readint.o: atomid.o atoms.o files.o inform.o iounit.o titles.o zclose.o zcoord.o
readmol.o: atomid.o atoms.o couple.o files.o iounit.o ptable.o titles.o
readmol2.o: atomid.o atoms.o couple.o files.o iounit.o ptable.o titles.o
readpdb.o: files.o inform.o iounit.o pdb.o resdue.o sequen.o titles.o
readprm.o: fields.o iounit.o kanang.o kangs.o kantor.o katoms.o kbonds.o kcflux.o kchrge.o kcpen.o kctrn.o kdipol.o kdsp.o kexpl.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolpr.o kpolr.o krepl.o ksolut.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o kvdws.o kxrepl.o merck.o params.o solute.o
readseq.o: files.o iounit.o resdue.o sequen.o
readxyz.o: atomid.o atoms.o boxes.o couple.o files.o inform.o iounit.o titles.o
refer.o: sizes.o
repel.o:
replica.o: bound.o boxes.o cell.o inform.o iounit.o math.o
reppot.o:
resdue.o:
respa.o: atomid.o atoms.o freeze.o ielscf.o limits.o mdstuf.o moldyn.o polar.o potent.o units.o usage.o virial.o
restrn.o:
rgddyn.o:
rgdstep.o: atomid.o atoms.o bound.o group.o iounit.o rgddyn.o units.o virial.o
rigid.o:
ring.o:
rings.o: angbnd.o atoms.o bitor.o bndstr.o couple.o inform.o iounit.o ring.o tors.o
rmsfit.o: align.o
rotbnd.o:
rotlist.o: atoms.o couple.o iounit.o molcul.o rotbnd.o zclose.o
rotpole.o: atoms.o math.o mpole.o repel.o
rxnfld.o:
rxnpot.o:
saddle.o: atoms.o inform.o iounit.o keys.o linmin.o minima.o syntrn.o titles.o zcoord.o
scales.o:
scan.o: atoms.o files.o inform.o iounit.o math.o minima.o omega.o output.o potent.o scales.o zcoord.o
sdstep.o: atomid.o atoms.o bath.o couple.o freeze.o kvdws.o math.o moldyn.o stodyn.o units.o usage.o virial.o
search.o: linmin.o math.o
sequen.o: sizes.o
server.o:
setprm.o: atoms.o couple.o kangs.o kantor.o kbonds.o kcflux.o kdipol.o keys.o khbond.o kiprop.o kitors.o kmulti.o kopbnd.o kopdst.o korbs.o kpitor.o kpolpr.o kstbnd.o ksttor.o ktorsn.o ktrtor.o kurybr.o kvdwpr.o params.o restrn.o
shakeup.o: angbnd.o atmlst.o atomid.o atoms.o bndstr.o bound.o couple.o freeze.o keys.o math.o molcul.o ring.o usage.o
shunt.o:
sigmoid.o:
simplex.o: inform.o iounit.o keys.o minima.o
sizes.o:
sktstuf.o: atomid.o atoms.o charge.o couple.o deriv.o fields.o files.o inform.o iounit.o keys.o moldyn.o mpole.o polar.o potent.o socket.o
sniffer.o: atoms.o files.o inform.o iounit.o linmin.o math.o minima.o output.o scales.o usage.o
socket.o:
solpot.o:
solute.o:
sort.o:
spacefill.o: atomid.o atoms.o files.o inform.o iounit.o kvdws.o math.o ptable.o usage.o
spectrum.o: files.o iounit.o math.o units.o
square.o: inform.o iounit.o keys.o minima.o
stodyn.o:
strbnd.o:
strtor.o:
suffix.o: ascii.o
superpose.o: align.o atomid.o atoms.o bound.o files.o inform.o iounit.o titles.o
surface.o: atoms.o inform.o iounit.o math.o usage.o
surfatom.o: atoms.o iounit.o math.o
switch.o: limits.o nonpol.o shunt.o
syntrn.o:
tarray.o:
tcgstuf.o: atoms.o iounit.o limits.o mpole.o polar.o poltcg.o potent.o
temper.o: atomid.o atoms.o bath.o group.o ielscf.o mdstuf.o molcul.o moldyn.o rgddyn.o units.o usage.o
testgrad.o: atoms.o deriv.o energi.o files.o inform.o inter.o iounit.o solute.o usage.o
testhess.o: atoms.o files.o hescut.o inform.o iounit.o usage.o
testpair.o: atoms.o deriv.o energi.o inform.o iounit.o light.o limits.o neigh.o polpot.o potent.o tarray.o vdwpot.o
testpol.o: atoms.o bound.o inform.o iounit.o limits.o minima.o mpole.o polar.o polopt.o polpot.o poltcg.o potent.o rigid.o units.o usage.o
testrot.o: domega.o energi.o inform.o iounit.o math.o omega.o zcoord.o
testvir.o: atoms.o bath.o bound.o boxes.o inform.o iounit.o math.o units.o virial.o
timer.o: atoms.o hescut.o inform.o iounit.o limits.o polpot.o
timerot.o: hescut.o inform.o iounit.o limits.o omega.o polpot.o
titles.o:
tncg.o: atoms.o hescut.o inform.o iounit.o keys.o linmin.o math.o minima.o output.o piorbs.o potent.o
torphase.o:
torpot.o:
torque.o: atoms.o deriv.o mpole.o
tors.o:
torsfit.o: atomid.o atoms.o files.o inform.o iounit.o keys.o ktorsn.o math.o output.o potent.o qmstuf.o restrn.o scales.o tors.o usage.o
torsions.o: atoms.o bndstr.o couple.o iounit.o tors.o
tortor.o:
tree.o:
trimtext.o:
unitcell.o: bound.o boxes.o keys.o math.o
units.o:
uprior.o:
urey.o:
urypot.o:
usage.o:
valence.o: angbnd.o angpot.o atomid.o atoms.o bndpot.o bndstr.o couple.o files.o hescut.o inform.o iounit.o kangs.o kbonds.o keys.o kopbnd.o kstbnd.o ktorsn.o kurybr.o kvdws.o linmin.o math.o minima.o opbend.o output.o potent.o qmstuf.o scales.o strbnd.o torpot.o tors.o units.o urey.o urypot.o usage.o valfit.o vdwpot.o
valfit.o:
vdw.o:
vdwpot.o:
verlet.o: atomid.o atoms.o freeze.o ielscf.o moldyn.o polar.o units.o usage.o
version.o: iounit.o output.o
vibbig.o: atomid.o atoms.o bound.o couple.o files.o hescut.o hessn.o inform.o iounit.o keys.o limits.o mpole.o potent.o rigid.o units.o usage.o vdw.o vdwpot.o vibs.o
vibrate.o: atomid.o atoms.o files.o hescut.o iounit.o math.o units.o usage.o
vibrot.o: iounit.o omega.o
vibs.o:
virial.o:
volume.o: atoms.o iounit.o math.o sizes.o
warp.o:
xrepel.o:
xtalfit.o: atomid.o atoms.o bound.o boxes.o charge.o dipole.o energi.o files.o fracs.o inform.o iounit.o kvdws.o limits.o math.o molcul.o mpole.o polar.o potent.o sizes.o vdw.o vdwpot.o xtals.o
xtalmin.o: atoms.o boxes.o files.o inform.o iounit.o keys.o math.o scales.o
xtals.o:
xyzatm.o: atoms.o inform.o iounit.o math.o
xyzedit.o: atomid.o atoms.o bound.o boxes.o couple.o energi.o files.o inform.o iounit.o katoms.o limits.o linmin.o math.o minima.o molcul.o neigh.o output.o potent.o ptable.o refer.o repel.o scales.o titles.o units.o usage.o vdw.o vdwpot.o
xyzint.o: files.o iounit.o titles.o
xyzmol2.o: files.o iounit.o titles.o
xyzpdb.o: atomid.o atoms.o couple.o fields.o files.o inform.o molcul.o pdb.o resdue.o sequen.o
zatom.o: angbnd.o atomid.o atoms.o bndstr.o fields.o iounit.o kangs.o katoms.o kbonds.o sizes.o zclose.o zcoord.o
zclose.o: sizes.o
zcoord.o: sizes.o

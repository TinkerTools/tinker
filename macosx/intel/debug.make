#
#
#  ################################################################
#  ##                                                            ##
#  ##  debug.make  --  debug compile each of the TINKER modules  ##
#  ##        (Intel Fortran Compiler for Mac OSX Version)        ##
#  ##                                                            ##
#  ################################################################
#
#
ifort -c -g -warn unused -check uninit -check bounds active.f
ifort -c -g -warn unused -check uninit -check bounds alchemy.f
ifort -c -g -warn unused -check uninit -check bounds analysis.f
ifort -c -g -warn unused -check uninit -check bounds analyze.f
ifort -c -g -warn unused -check uninit -check bounds angles.f
ifort -c -g -warn unused -check uninit -check bounds anneal.f
ifort -c -g -warn unused -check uninit -check bounds archive.f
ifort -c -g -warn unused -check uninit -check bounds attach.f
ifort -c -g -warn unused -check uninit -check bounds basefile.f
ifort -c -g -warn unused -check uninit -check bounds beeman.f
ifort -c -g -warn unused -check uninit -check bounds bicubic.f
ifort -c -g -warn unused -check uninit -check bounds bitors.f
ifort -c -g -warn unused -check uninit -check bounds bonds.f
ifort -c -g -warn unused -check uninit -check bounds born.f
ifort -c -g -warn unused -check uninit -check bounds bounds.f
ifort -c -g -warn unused -check uninit -check bounds bussi.f
ifort -c -g -warn unused -check uninit -check bounds calendar.f
ifort -c -g -warn unused -check uninit -check bounds center.f
ifort -c -g -warn unused -check uninit -check bounds chkpole.f
ifort -c -g -warn unused -check uninit -check bounds chkring.f
ifort -c -g -warn unused -check uninit -check bounds chkxyz.f
ifort -c -g -warn unused -check uninit -check bounds cholesky.f
ifort -c -g -warn unused -check uninit -check bounds clock.f
ifort -c -g -warn unused -check uninit -check bounds cluster.f
ifort -c -g -warn unused -check uninit -check bounds column.f
ifort -c -g -warn unused -check uninit -check bounds command.f
ifort -c -g -warn unused -check uninit -check bounds connect.f
ifort -c -g -warn unused -check uninit -check bounds connolly.f
ifort -c -g -warn unused -check uninit -check bounds control.f
ifort -c -g -warn unused -check uninit -check bounds correlate.f
ifort -c -g -warn unused -check uninit -check bounds crystal.f
ifort -c -g -warn unused -check uninit -check bounds cspline.f
ifort -c -g -warn unused -check uninit -check bounds cutoffs.f
ifort -c -g -warn unused -check uninit -check bounds deflate.f
ifort -c -g -warn unused -check uninit -check bounds delete.f
ifort -c -g -warn unused -check uninit -check bounds diagq.f
ifort -c -g -warn unused -check uninit -check bounds diffeq.f
ifort -c -g -warn unused -check uninit -check bounds diffuse.f
ifort -c -g -warn unused -check uninit -check bounds distgeom.f
ifort -c -g -warn unused -check uninit -check bounds document.f
ifort -c -g -warn unused -check uninit -check bounds dynamic.f
ifort -c -g -warn unused -check uninit -check bounds eangang.f
ifort -c -g -warn unused -check uninit -check bounds eangang1.f
ifort -c -g -warn unused -check uninit -check bounds eangang2.f
ifort -c -g -warn unused -check uninit -check bounds eangang3.f
ifort -c -g -warn unused -check uninit -check bounds eangle.f
ifort -c -g -warn unused -check uninit -check bounds eangle1.f
ifort -c -g -warn unused -check uninit -check bounds eangle2.f
ifort -c -g -warn unused -check uninit -check bounds eangle3.f
ifort -c -g -warn unused -check uninit -check bounds ebond.f
ifort -c -g -warn unused -check uninit -check bounds ebond1.f
ifort -c -g -warn unused -check uninit -check bounds ebond2.f
ifort -c -g -warn unused -check uninit -check bounds ebond3.f
ifort -c -g -warn unused -check uninit -check bounds ebuck.f
ifort -c -g -warn unused -check uninit -check bounds ebuck1.f
ifort -c -g -warn unused -check uninit -check bounds ebuck2.f
ifort -c -g -warn unused -check uninit -check bounds ebuck3.f
ifort -c -g -warn unused -check uninit -check bounds echarge.f
ifort -c -g -warn unused -check uninit -check bounds echarge1.f
ifort -c -g -warn unused -check uninit -check bounds echarge2.f
ifort -c -g -warn unused -check uninit -check bounds echarge3.f
ifort -c -g -warn unused -check uninit -check bounds echgdpl.f
ifort -c -g -warn unused -check uninit -check bounds echgdpl1.f
ifort -c -g -warn unused -check uninit -check bounds echgdpl2.f
ifort -c -g -warn unused -check uninit -check bounds echgdpl3.f
ifort -c -g -warn unused -check uninit -check bounds edipole.f
ifort -c -g -warn unused -check uninit -check bounds edipole1.f
ifort -c -g -warn unused -check uninit -check bounds edipole2.f
ifort -c -g -warn unused -check uninit -check bounds edipole3.f
ifort -c -g -warn unused -check uninit -check bounds egauss.f
ifort -c -g -warn unused -check uninit -check bounds egauss1.f
ifort -c -g -warn unused -check uninit -check bounds egauss2.f
ifort -c -g -warn unused -check uninit -check bounds egauss3.f
ifort -c -g -warn unused -check uninit -check bounds egeom.f
ifort -c -g -warn unused -check uninit -check bounds egeom1.f
ifort -c -g -warn unused -check uninit -check bounds egeom2.f
ifort -c -g -warn unused -check uninit -check bounds egeom3.f
ifort -c -g -warn unused -check uninit -check bounds ehal.f
ifort -c -g -warn unused -check uninit -check bounds ehal1.f
ifort -c -g -warn unused -check uninit -check bounds ehal2.f
ifort -c -g -warn unused -check uninit -check bounds ehal3.f
ifort -c -g -warn unused -check uninit -check bounds eimprop.f
ifort -c -g -warn unused -check uninit -check bounds eimprop1.f
ifort -c -g -warn unused -check uninit -check bounds eimprop2.f
ifort -c -g -warn unused -check uninit -check bounds eimprop3.f
ifort -c -g -warn unused -check uninit -check bounds eimptor.f
ifort -c -g -warn unused -check uninit -check bounds eimptor1.f
ifort -c -g -warn unused -check uninit -check bounds eimptor2.f
ifort -c -g -warn unused -check uninit -check bounds eimptor3.f
ifort -c -g -warn unused -check uninit -check bounds elj.f
ifort -c -g -warn unused -check uninit -check bounds elj1.f
ifort -c -g -warn unused -check uninit -check bounds elj2.f
ifort -c -g -warn unused -check uninit -check bounds elj3.f
ifort -c -g -warn unused -check uninit -check bounds embed.f
ifort -c -g -warn unused -check uninit -check bounds emetal.f
ifort -c -g -warn unused -check uninit -check bounds emetal1.f
ifort -c -g -warn unused -check uninit -check bounds emetal2.f
ifort -c -g -warn unused -check uninit -check bounds emetal3.f
ifort -c -g -warn unused -check uninit -check bounds emm3hb.f
ifort -c -g -warn unused -check uninit -check bounds emm3hb1.f
ifort -c -g -warn unused -check uninit -check bounds emm3hb2.f
ifort -c -g -warn unused -check uninit -check bounds emm3hb3.f
ifort -c -g -warn unused -check uninit -check bounds empole.f
ifort -c -g -warn unused -check uninit -check bounds empole1.f
ifort -c -g -warn unused -check uninit -check bounds empole2.f
ifort -c -g -warn unused -check uninit -check bounds empole3.f
ifort -c -g -warn unused -check uninit -check bounds energy.f
ifort -c -g -warn unused -check uninit -check bounds eopbend.f
ifort -c -g -warn unused -check uninit -check bounds eopbend1.f
ifort -c -g -warn unused -check uninit -check bounds eopbend2.f
ifort -c -g -warn unused -check uninit -check bounds eopbend3.f
ifort -c -g -warn unused -check uninit -check bounds eopdist.f
ifort -c -g -warn unused -check uninit -check bounds eopdist1.f
ifort -c -g -warn unused -check uninit -check bounds eopdist2.f
ifort -c -g -warn unused -check uninit -check bounds eopdist3.f
ifort -c -g -warn unused -check uninit -check bounds epitors.f
ifort -c -g -warn unused -check uninit -check bounds epitors1.f
ifort -c -g -warn unused -check uninit -check bounds epitors2.f
ifort -c -g -warn unused -check uninit -check bounds epitors3.f
ifort -c -g -warn unused -check uninit -check bounds erf.f
ifort -c -g -warn unused -check uninit -check bounds erxnfld.f
ifort -c -g -warn unused -check uninit -check bounds erxnfld1.f
ifort -c -g -warn unused -check uninit -check bounds erxnfld2.f
ifort -c -g -warn unused -check uninit -check bounds erxnfld3.f
ifort -c -g -warn unused -check uninit -check bounds esolv.f
ifort -c -g -warn unused -check uninit -check bounds esolv1.f
ifort -c -g -warn unused -check uninit -check bounds esolv2.f
ifort -c -g -warn unused -check uninit -check bounds esolv3.f
ifort -c -g -warn unused -check uninit -check bounds estrbnd.f
ifort -c -g -warn unused -check uninit -check bounds estrbnd1.f
ifort -c -g -warn unused -check uninit -check bounds estrbnd2.f
ifort -c -g -warn unused -check uninit -check bounds estrbnd3.f
ifort -c -g -warn unused -check uninit -check bounds estrtor.f
ifort -c -g -warn unused -check uninit -check bounds estrtor1.f
ifort -c -g -warn unused -check uninit -check bounds estrtor2.f
ifort -c -g -warn unused -check uninit -check bounds estrtor3.f
ifort -c -g -warn unused -check uninit -check bounds etors.f
ifort -c -g -warn unused -check uninit -check bounds etors1.f
ifort -c -g -warn unused -check uninit -check bounds etors2.f
ifort -c -g -warn unused -check uninit -check bounds etors3.f
ifort -c -g -warn unused -check uninit -check bounds etortor.f
ifort -c -g -warn unused -check uninit -check bounds etortor1.f
ifort -c -g -warn unused -check uninit -check bounds etortor2.f
ifort -c -g -warn unused -check uninit -check bounds etortor3.f
ifort -c -g -warn unused -check uninit -check bounds eurey.f
ifort -c -g -warn unused -check uninit -check bounds eurey1.f
ifort -c -g -warn unused -check uninit -check bounds eurey2.f
ifort -c -g -warn unused -check uninit -check bounds eurey3.f
ifort -c -g -warn unused -check uninit -check bounds evcorr.f
ifort -c -g -warn unused -check uninit -check bounds extra.f
ifort -c -g -warn unused -check uninit -check bounds extra1.f
ifort -c -g -warn unused -check uninit -check bounds extra2.f
ifort -c -g -warn unused -check uninit -check bounds extra3.f
ifort -c -g -warn unused -check uninit -check bounds fatal.f
ifort -c -g -warn unused -check uninit -check bounds fft3d.f
ifort -c -g -warn unused -check uninit -check bounds fftpack.f
ifort -c -g -warn unused -check uninit -check bounds field.f
ifort -c -g -warn unused -check uninit -check bounds final.f
ifort -c -g -warn unused -check uninit -check bounds flatten.f
ifort -c -g -warn unused -check uninit -check bounds freeunit.f
ifort -c -g -warn unused -check uninit -check bounds gda.f
ifort -c -g -warn unused -check uninit -check bounds geometry.f
ifort -c -g -warn unused -check uninit -check bounds getint.f
ifort -c -g -warn unused -check uninit -check bounds getkey.f
ifort -c -g -warn unused -check uninit -check bounds getmol2.f
ifort -c -g -warn unused -check uninit -check bounds getnumb.f
ifort -c -g -warn unused -check uninit -check bounds getpdb.f
ifort -c -g -warn unused -check uninit -check bounds getprm.f
ifort -c -g -warn unused -check uninit -check bounds getref.f
ifort -c -g -warn unused -check uninit -check bounds getstring.f
ifort -c -g -warn unused -check uninit -check bounds gettext.f
ifort -c -g -warn unused -check uninit -check bounds getword.f
ifort -c -g -warn unused -check uninit -check bounds getxyz.f
ifort -c -g -warn unused -check uninit -check bounds ghmcstep.f
ifort -c -g -warn unused -check uninit -check bounds gradient.f
ifort -c -g -warn unused -check uninit -check bounds gradrgd.f
ifort -c -g -warn unused -check uninit -check bounds gradrot.f
ifort -c -g -warn unused -check uninit -check bounds groups.f
ifort -c -g -warn unused -check uninit -check bounds grpline.f
ifort -c -g -warn unused -check uninit -check bounds gyrate.f
ifort -c -g -warn unused -check uninit -check bounds hessian.f
ifort -c -g -warn unused -check uninit -check bounds hessrgd.f
ifort -c -g -warn unused -check uninit -check bounds hessrot.f
ifort -c -g -warn unused -check uninit -check bounds hybrid.f
ifort -c -g -warn unused -check uninit -check bounds image.f
ifort -c -g -warn unused -check uninit -check bounds impose.f
ifort -c -g -warn unused -check uninit -check bounds induce.f
ifort -c -g -warn unused -check uninit -check bounds inertia.f
ifort -c -g -warn unused -check uninit -check bounds initial.f
ifort -c -g -warn unused -check uninit -check bounds initprm.f
ifort -c -g -warn unused -check uninit -check bounds initres.f
ifort -c -g -warn unused -check uninit -check bounds initrot.f
ifort -c -g -warn unused -check uninit -check bounds insert.f
ifort -c -g -warn unused -check uninit -check bounds intedit.f
ifort -c -g -warn unused -check uninit -check bounds intxyz.f
ifort -c -g -warn unused -check uninit -check bounds invbeta.f
ifort -c -g -warn unused -check uninit -check bounds invert.f
ifort -c -g -warn unused -check uninit -check bounds jacobi.f
ifort -c -g -warn unused -check uninit -check bounds kangang.f
ifort -c -g -warn unused -check uninit -check bounds kangle.f
ifort -c -g -warn unused -check uninit -check bounds katom.f
ifort -c -g -warn unused -check uninit -check bounds kbond.f
ifort -c -g -warn unused -check uninit -check bounds kcharge.f
ifort -c -g -warn unused -check uninit -check bounds kdipole.f
ifort -c -g -warn unused -check uninit -check bounds kewald.f
ifort -c -g -warn unused -check uninit -check bounds kgeom.f
ifort -c -g -warn unused -check uninit -check bounds kimprop.f
ifort -c -g -warn unused -check uninit -check bounds kimptor.f
ifort -c -g -warn unused -check uninit -check bounds kinetic.f
ifort -c -g -warn unused -check uninit -check bounds kmetal.f
ifort -c -g -warn unused -check uninit -check bounds kmpole.f
ifort -c -g -warn unused -check uninit -check bounds kopbend.f
ifort -c -g -warn unused -check uninit -check bounds kopdist.f
ifort -c -g -warn unused -check uninit -check bounds korbit.f
ifort -c -g -warn unused -check uninit -check bounds kpitors.f
ifort -c -g -warn unused -check uninit -check bounds kpolar.f
ifort -c -g -warn unused -check uninit -check bounds ksolv.f
ifort -c -g -warn unused -check uninit -check bounds kstrbnd.f
ifort -c -g -warn unused -check uninit -check bounds kstrtor.f
ifort -c -g -warn unused -check uninit -check bounds ktors.f
ifort -c -g -warn unused -check uninit -check bounds ktortor.f
ifort -c -g -warn unused -check uninit -check bounds kurey.f
ifort -c -g -warn unused -check uninit -check bounds kvdw.f
ifort -c -g -warn unused -check uninit -check bounds lattice.f
ifort -c -g -warn unused -check uninit -check bounds lbfgs.f
ifort -c -g -warn unused -check uninit -check bounds lights.f
ifort -c -g -warn unused -check uninit -check bounds makeint.f
ifort -c -g -warn unused -check uninit -check bounds makeref.f
ifort -c -g -warn unused -check uninit -check bounds makexyz.f
ifort -c -g -warn unused -check uninit -check bounds maxwell.f
ifort -c -g -warn unused -check uninit -check bounds mdinit.f
ifort -c -g -warn unused -check uninit -check bounds mdrest.f
ifort -c -g -warn unused -check uninit -check bounds mdsave.f
ifort -c -g -warn unused -check uninit -check bounds mdstat.f
ifort -c -g -warn unused -check uninit -check bounds mechanic.f
ifort -c -g -warn unused -check uninit -check bounds merge.f
ifort -c -g -warn unused -check uninit -check bounds minimize.f
ifort -c -g -warn unused -check uninit -check bounds minirot.f
ifort -c -g -warn unused -check uninit -check bounds minrigid.f
ifort -c -g -warn unused -check uninit -check bounds molecule.f
ifort -c -g -warn unused -check uninit -check bounds moments.f
ifort -c -g -warn unused -check uninit -check bounds monte.f
ifort -c -g -warn unused -check uninit -check bounds mutate.f
ifort -c -g -warn unused -check uninit -check bounds nblist.f
ifort -c -g -warn unused -check uninit -check bounds newton.f
ifort -c -g -warn unused -check uninit -check bounds newtrot.f
ifort -c -g -warn unused -check uninit -check bounds nextarg.f
ifort -c -g -warn unused -check uninit -check bounds nexttext.f
ifort -c -g -warn unused -check uninit -check bounds nose.f
ifort -c -g -warn unused -check uninit -check bounds nspline.f
ifort -c -g -warn unused -check uninit -check bounds nucleic.f
ifort -c -g -warn unused -check uninit -check bounds number.f
ifort -c -g -warn unused -check uninit -check bounds numeral.f
ifort -c -g -warn unused -check uninit -check bounds numgrad.f
ifort -c -g -warn unused -check uninit -check bounds ocvm.f
ifort -c -g -warn unused -check uninit -check bounds openend.f
ifort -c -g -warn unused -check uninit -check bounds optimize.f
ifort -c -g -warn unused -check uninit -check bounds optirot.f
ifort -c -g -warn unused -check uninit -check bounds optrigid.f
ifort -c -g -warn unused -check uninit -check bounds optsave.f
ifort -c -g -warn unused -check uninit -check bounds orbital.f
ifort -c -g -warn unused -check uninit -check bounds orient.f
ifort -c -g -warn unused -check uninit -check bounds orthog.f
ifort -c -g -warn unused -check uninit -check bounds overlap.f
ifort -c -g -warn unused -check uninit -check bounds path.f
ifort -c -g -warn unused -check uninit -check bounds pdbxyz.f
ifort -c -g -warn unused -check uninit -check bounds piscf.f
ifort -c -g -warn unused -check uninit -check bounds pmestuff.f
ifort -c -g -warn unused -check uninit -check bounds pmpb.f
ifort -c -g -warn unused -check uninit -check bounds polarize.f
ifort -c -g -warn unused -check uninit -check bounds poledit.f
ifort -c -g -warn unused -check uninit -check bounds polymer.f
ifort -c -g -warn unused -check uninit -check bounds potential.f
ifort -c -g -warn unused -check uninit -check bounds precise.f
ifort -c -g -warn unused -check uninit -check bounds pressure.f
ifort -c -g -warn unused -check uninit -check bounds prmedit.f
ifort -c -g -warn unused -check uninit -check bounds prmkey.f
ifort -c -g -warn unused -check uninit -check bounds promo.f
ifort -c -g -warn unused -check uninit -check bounds protein.f
ifort -c -g -warn unused -check uninit -check bounds prtdyn.f
ifort -c -g -warn unused -check uninit -check bounds prterr.f
ifort -c -g -warn unused -check uninit -check bounds prtint.f
ifort -c -g -warn unused -check uninit -check bounds prtmol2.f
ifort -c -g -warn unused -check uninit -check bounds prtpdb.f
ifort -c -g -warn unused -check uninit -check bounds prtprm.f
ifort -c -g -warn unused -check uninit -check bounds prtseq.f
ifort -c -g -warn unused -check uninit -check bounds prtxyz.f
ifort -c -g -warn unused -check uninit -check bounds pss.f
ifort -c -g -warn unused -check uninit -check bounds pssrigid.f
ifort -c -g -warn unused -check uninit -check bounds pssrot.f
ifort -c -g -warn unused -check uninit -check bounds quatfit.f
ifort -c -g -warn unused -check uninit -check bounds radial.f
ifort -c -g -warn unused -check uninit -check bounds random.f
ifort -c -g -warn unused -check uninit -check bounds rattle.f
ifort -c -g -warn unused -check uninit -check bounds readdyn.f
ifort -c -g -warn unused -check uninit -check bounds readgau.f
ifort -c -g -warn unused -check uninit -check bounds readint.f
ifort -c -g -warn unused -check uninit -check bounds readmol2.f
ifort -c -g -warn unused -check uninit -check bounds readpdb.f
ifort -c -g -warn unused -check uninit -check bounds readprm.f
ifort -c -g -warn unused -check uninit -check bounds readseq.f
ifort -c -g -warn unused -check uninit -check bounds readxyz.f
ifort -c -g -warn unused -check uninit -check bounds replica.f
ifort -c -g -warn unused -check uninit -check bounds respa.f
ifort -c -g -warn unused -check uninit -check bounds rgdstep.f
ifort -c -g -warn unused -check uninit -check bounds rings.f
ifort -c -g -warn unused -check uninit -check bounds rmsfit.f
ifort -c -g -warn unused -check uninit -check bounds rotlist.f
ifort -c -g -warn unused -check uninit -check bounds rotpole.f
ifort -c -g -warn unused -check uninit -check bounds saddle.f
ifort -c -g -warn unused -check uninit -check bounds scan.f
ifort -c -g -warn unused -check uninit -check bounds sdstep.f
ifort -c -g -warn unused -check uninit -check bounds search.f
ifort -c -g -warn unused -check uninit -check bounds server.f
ifort -c -g -warn unused -check uninit -check bounds shakeup.f
ifort -c -g -warn unused -check uninit -check bounds sigmoid.f
ifort -c -g -warn unused -check uninit -check bounds sktstuff.f
ifort -c -g -warn unused -check uninit -check bounds sniffer.f
ifort -c -g -warn unused -check uninit -check bounds sort.f
ifort -c -g -warn unused -check uninit -check bounds spacefill.f
ifort -c -g -warn unused -check uninit -check bounds spectrum.f
ifort -c -g -warn unused -check uninit -check bounds square.f
ifort -c -g -warn unused -check uninit -check bounds suffix.f
ifort -c -g -warn unused -check uninit -check bounds superpose.f
ifort -c -g -warn unused -check uninit -check bounds surface.f
ifort -c -g -warn unused -check uninit -check bounds surfatom.f
ifort -c -g -warn unused -check uninit -check bounds switch.f
ifort -c -g -warn unused -check uninit -check bounds sybylxyz.f
ifort -c -g -warn unused -check uninit -check bounds temper.f
ifort -c -g -warn unused -check uninit -check bounds testgrad.f
ifort -c -g -warn unused -check uninit -check bounds testhess.f
ifort -c -g -warn unused -check uninit -check bounds testpair.f
ifort -c -g -warn unused -check uninit -check bounds testrot.f
ifort -c -g -warn unused -check uninit -check bounds timer.f
ifort -c -g -warn unused -check uninit -check bounds timerot.f
ifort -c -g -warn unused -check uninit -check bounds tncg.f
ifort -c -g -warn unused -check uninit -check bounds torphase.f
ifort -c -g -warn unused -check uninit -check bounds torque.f
ifort -c -g -warn unused -check uninit -check bounds torsfit.f
ifort -c -g -warn unused -check uninit -check bounds torsions.f
ifort -c -g -warn unused -check uninit -check bounds trimtext.f
ifort -c -g -warn unused -check uninit -check bounds unitcell.f
ifort -c -g -warn unused -check uninit -check bounds valence.f
ifort -c -g -warn unused -check uninit -check bounds verlet.f
ifort -c -g -warn unused -check uninit -check bounds version.f
ifort -c -g -warn unused -check uninit -check bounds vibbig.f
ifort -c -g -warn unused -check uninit -check bounds vibrate.f
ifort -c -g -warn unused -check uninit -check bounds vibrot.f
ifort -c -g -warn unused -check uninit -check bounds volume.f
ifort -c -g -warn unused -check uninit -check bounds xtalfit.f
ifort -c -g -warn unused -check uninit -check bounds xtalmin.f
ifort -c -g -warn unused -check uninit -check bounds xyzatm.f
ifort -c -g -warn unused -check uninit -check bounds xyzedit.f
ifort -c -g -warn unused -check uninit -check bounds xyzint.f
ifort -c -g -warn unused -check uninit -check bounds xyzpdb.f
ifort -c -g -warn unused -check uninit -check bounds xyzsybyl.f
ifort -c -g -warn unused -check uninit -check bounds zatom.f

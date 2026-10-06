c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_fep  --  trial lambda energy tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_fep" checks the energies that "efeptrial" returns at
c     trial values of the main lambda, and that getting them leaves
c     the native state of the system as it was found
c
c
      subroutine test_fep
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_fep')
c
c
      if (skiptest(tname,'fep'))  return
      call test_fep_trial
      call test_fep_terms
      call test_fep_state
      call test_fep_predict
      call test_fep_dynamic
      call test_fep_sample
      call test_fep_run
      call test_fep_stop
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_fep_trial  --  trial energy values  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_fep_trial" compares the trial energies against separate
c     setups of the same system made at each trial lambda, for the
c     absolute and staged relative modes with single and dual
c     topology polarization
c
c
      subroutine test_fep_trial
      implicit none
      real*8 ltrial(7)
      character*11 dm
      parameter (dm='file/mutate')
c
c
c     absolute single topology with all three sublambdas mapped
c
      ltrial(1) = 0.0d0
      ltrial(2) = 0.4d0
      ltrial(3) = 0.6d0
      ltrial(4) = 1.0d0
      call fepcase (dm,'water2','164_water_ast_ne_mcut_d1_l05.key',
     &              4,ltrial,'ast')
      call fepcase (dm,'ionwat','179_ionwat_ast_l05.key',
     &              4,ltrial,'ast ewald ion')
c
c     van der Waals annihilation with its long range correction
c
      call fepcase (dm,'water2','170_water_ast_vcorr_annih_d1_l05.key',
     &              4,ltrial,'ast vcorr')
c
c     the quintic taper map of the electrostatic sublambda
c
      call fepcase (dm,'water2','093_water_lmda_qnt_l05.key',
     &              4,ltrial,'ast qnt')
c
c     absolute dual topology polarization, with and without Ewald
c
      ltrial(1) = 0.0d0
      ltrial(2) = 0.5d0
      ltrial(3) = 0.7d0
      ltrial(4) = 1.0d0
      call fepcase (dm,'water2','169_water_adt_d1_x2_l06.key',
     &              4,ltrial,'adt ewald')
      call fepcase (dm,'water2','178_water_adt_d1_ne_l06.key',
     &              4,ltrial,'adt')
c
c     staged relative, native in the last leg with trials in every
c     leg, on a stage boundary and at both ends
c
      ltrial(1) = 0.0d0
      ltrial(2) = 0.15d0
      ltrial(3) = 0.5d0
      ltrial(4) = 0.7d0
      ltrial(5) = 0.8d0
      ltrial(6) = 0.9d0
      ltrial(7) = 1.0d0
      call fepcase (dm,'water2','190_water_rels3_l085.key',
     &              7,ltrial,'rels dt')
      call fepcase (dm,'water2','196_water_rels3_st_l085.key',
     &              7,ltrial,'rels st')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_fep_terms  --  trial energy of each term  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_fep_terms" repeats the trial energy and native state
c     checks for the kinds of potential term and polarization model
c     that the other cases do not use, since the method relies only
c     on each term being correct at whatever lambda is in effect
c
c
      subroutine test_fep_terms
      implicit none
      real*8 ltrial(4)
      character*8 df
      parameter (df='file/fep')
c
c
      ltrial(1) = 0.0d0
      ltrial(2) = 0.4d0
      ltrial(3) = 0.6d0
      ltrial(4) = 1.0d0
c
c     extrapolated and truncated conjugate gradient polarization
c
      call fepcase (df,'../mutate/water2','opt.key',4,ltrial,'opt')
      call fepstate (df,'../mutate/water2','opt.key',4,ltrial,'opt')
      call fepcase (df,'../mutate/water2','tcg.key',4,ltrial,'tcg')
      call fepstate (df,'../mutate/water2','tcg.key',4,ltrial,'tcg')
c
c     HIPPO repulsion, dispersion, charge transfer and penetration
c
      call fepcase (df,'../ephippo/dmso','hippo.key',4,ltrial,'hippo')
      call fepstate (df,'../ephippo/dmso','hippo.key',4,ltrial,'hippo')
c
c     AMOEBA+ with charge flux, charge transfer and penetration
c
      call fepcase (df,'../aplusliquid/tetramer','aplus.key',
     &              4,ltrial,'aplus')
      call fepstate (df,'../aplusliquid/tetramer','aplus.key',
     &               4,ltrial,'aplus')
c
c     partial charges in a periodic box
c
      call fepcase (df,'../chglj/trp_charmm','charge.key',
     &              4,ltrial,'charge')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine fepcase  --  trial energies of one system  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "fepcase" gets the trial energies of one fixture and compares
c     each of them against the energy of that same fixture set up
c     afresh with the trial value as its main lambda; it also checks
c     the native energy and that a second call repeats the first
c
c
      subroutine fepcase (dir,base,key,ntrial,ltrial,tag)
      use energi
      use mutant
      implicit none
      integer i
      integer ntrial
      real*8 energy
      real*8 eps
      real*8 eref
      real*8 enative
      real*8 enative2
      real*8 lnative
      real*8 ltrial(*)
      real*8 etrial(7)
      real*8 etrial2(7)
      character*(*) dir,base,key,tag
      character*40 addkey
      character*60 label
c
c
c     get the trial energies twice from the native setup
c
      eps = 1.0d-8
      call pushdir (dir)
      call loadfix (base,key)
      lnative = lambda
      call efeptrial (ntrial,ltrial,etrial,enative)
      call assert_real (lambda,lnative,0.0d0,
     &                  'fep '//tag//' lambda restored')
      call assert_real (esum,enative,0.0d0,
     &                  'fep '//tag//' native energy left')
      call efeptrial (ntrial,ltrial,etrial2,enative2)
      call assert_real (enative2,enative,eps,
     &                  'fep '//tag//' native repeat')
      do i = 1, ntrial
         write (label,10)  tag,ltrial(i)
   10    format ('fep ',a,' repeat at',f6.2)
         call assert_real (etrial2(i),etrial(i),eps,label)
      end do
c
c     the native energy is that of an ordinary energy call
c
      eref = energy ()
      call assert_real (enative,eref,eps,'fep '//tag//' native energy')
      call final
c
c     set the system up again at each trial lambda for a reference
c
      do i = 1, ntrial
         write (addkey,20)  ltrial(i)
   20    format ('lambda ',f12.8)
         call loadfix_keyadd (base,key,addkey)
         eref = energy ()
         write (label,30)  tag,ltrial(i)
   30    format ('fep ',a,' trial at',f6.2)
         call assert_real (etrial(i),eref,eps,label)
         call final
      end do
      call popdir
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_fep_state  --  native state restored  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_fep_state" checks that the lambda values, the lambda
c     scaled parameters, the induced dipoles and the structure are
c     all exactly as they were before the trial energies were found
c
c
      subroutine test_fep_state
      implicit none
      real*8 ltrial(4)
      character*11 dm
      parameter (dm='file/mutate')
c
c
      ltrial(1) = 0.0d0
      ltrial(2) = 0.3d0
      ltrial(3) = 0.75d0
      ltrial(4) = 1.0d0
      call fepstate (dm,'ionwat','179_ionwat_ast_l05.key',
     &               4,ltrial,'ast')
      call fepstate (dm,'water2','169_water_adt_d1_x2_l06.key',
     &               4,ltrial,'adt')
      call fepstate (dm,'water2','190_water_rels3_l085.key',
     &               4,ltrial,'rels dt')
      call fepstate (dm,'water2','196_water_rels3_st_l085.key',
     &               4,ltrial,'rels st')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine fepstate  --  state of one system restored  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "fepstate" records the native state of one fixture after an
c     energy evaluation, gets a set of trial energies, and finds the
c     largest change in each part of that state, which must be zero
c
c
      subroutine fepstate (dir,base,key,ntrial,ltrial,tag)
      use atoms
      use boxes
      use dlmda
      use mpole
      use mutant
      use polar
      implicit none
      integer i,j
      integer ntrial
      integer nflip
      real*8 energy
      real*8 e
      real*8 enative
      real*8 dmax
      real*8 xboxold
      real*8 lold(4)
      real*8 dlold(3)
      real*8 ltrial(*)
      real*8 etrial(7)
      real*8, allocatable :: xold(:,:)
      real*8, allocatable :: poleold(:,:)
      real*8, allocatable :: polarold(:)
      real*8, allocatable :: uold(:,:,:)
      logical, allocatable :: douold(:)
      character*4 stageold
      character*(*) dir,base,key,tag
c
c
c     evaluate the energy once, as the dynamics would have done
c
      call pushdir (dir)
      call loadfix (base,key)
      e = energy ()
c
c     record the native state
c
      allocate (xold(3,n))
      allocate (poleold(13,n))
      allocate (polarold(n))
      allocate (uold(3,n,4))
      allocate (douold(n))
      lold(1) = lambda
      lold(2) = elambda
      lold(3) = plambda
      lold(4) = vlambda
      dlold(1) = deldlmda
      dlold(2) = dpldlmda
      dlold(3) = dvldlmda
      stageold = relstage
      xboxold = xbox
      do i = 1, n
         xold(1,i) = x(i)
         xold(2,i) = y(i)
         xold(3,i) = z(i)
         polarold(i) = polarity(i)
         douold(i) = douind(i)
         do j = 1, 13
            poleold(j,i) = pole(j,i)
         end do
         do j = 1, 3
            uold(j,i,1) = uind(j,i)
            uold(j,i,2) = uinp(j,i)
            uold(j,i,3) = udir(j,i)
            uold(j,i,4) = udirp(j,i)
         end do
      end do
c
c     get the trial energies, then compare each part of the state
c
      call efeptrial (ntrial,ltrial,etrial,enative)
      dmax = abs(lambda-lold(1))
      dmax = max(dmax,abs(elambda-lold(2)))
      dmax = max(dmax,abs(plambda-lold(3)))
      dmax = max(dmax,abs(vlambda-lold(4)))
      call assert_real (dmax,0.0d0,0.0d0,'fep '//tag//' state lambdas')
      dmax = abs(deldlmda-dlold(1))
      dmax = max(dmax,abs(dpldlmda-dlold(2)))
      dmax = max(dmax,abs(dvldlmda-dlold(3)))
      call assert_real (dmax,0.0d0,0.0d0,'fep '//tag//' state maps')
      call assert_logical (relstage.eq.stageold,.true.,
     &                     'fep '//tag//' state stage')
      dmax = abs(xbox-xboxold)
      do i = 1, n
         dmax = max(dmax,abs(x(i)-xold(1,i)))
         dmax = max(dmax,abs(y(i)-xold(2,i)))
         dmax = max(dmax,abs(z(i)-xold(3,i)))
      end do
      call assert_real (dmax,0.0d0,0.0d0,'fep '//tag//' state coords')
      dmax = 0.0d0
      nflip = 0
      do i = 1, n
         dmax = max(dmax,abs(polarity(i)-polarold(i)))
         if (douind(i) .neqv. douold(i))  nflip = nflip + 1
         do j = 1, 13
            dmax = max(dmax,abs(pole(j,i)-poleold(j,i)))
         end do
      end do
      call assert_real (dmax,0.0d0,0.0d0,'fep '//tag//' state params')
      call assert_int (nflip,0,'fep '//tag//' state douind')
      dmax = 0.0d0
      do i = 1, n
         do j = 1, 3
            dmax = max(dmax,abs(uind(j,i)-uold(j,i,1)))
            dmax = max(dmax,abs(uinp(j,i)-uold(j,i,2)))
            dmax = max(dmax,abs(udir(j,i)-uold(j,i,3)))
            dmax = max(dmax,abs(udirp(j,i)-uold(j,i,4)))
         end do
      end do
      call assert_real (dmax,0.0d0,0.0d0,'fep '//tag//' state dipoles')
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (poleold)
      deallocate (polarold)
      deallocate (uold)
      deallocate (douold)
      call final
      call popdir
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_fep_predict  --  prediction forced off  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_fep_predict" checks that a request for induced dipole
c     prediction is honored in general, but is turned off when free
c     energy perturbation over lambda windows is in use
c
c
      subroutine test_fep_predict
      use dlmda
      use ielscf
      use uprior
      implicit none
      real*8 energy
      real*8 e
c
c
c     a polynomial predictor is kept unless the method is in use
c
      call pushdir ('file/mutate')
      call loadfix_keyadd ('water2','164_water_ast_ne_mcut_d1_l05.key',
     &                     'polar-predict aspc')
      call predict
      call assert_logical (use_pred,.true.,'fep predict aspc kept')
      call assert_int (maxualt,17,'fep predict aspc history')
      use_fep = .true.
      call predict
      call assert_logical (use_pred,.false.,'fep predict aspc off')
      call assert_int (maxualt,0,'fep predict aspc no history')
      use_fep = .false.
      call final
c
c     the extended Lagrangian predictor is treated the same way;
c     its setup solves for the induced dipoles, so an energy is
c     found first to put the multipoles in the global frame
c
      call loadfix_keyadd ('water2','164_water_ast_ne_mcut_d1_l05.key',
     &                     'polar-predict iel')
      e = energy ()
      call predict
      call assert_logical (use_ielscf,.true.,'fep predict iel kept')
      use_fep = .true.
      call predict
      call assert_logical (use_ielscf,.false.,'fep predict iel off')
      call assert_logical (use_pred,.false.,'fep predict iel no pred')
      use_fep = .false.
      call final
      call popdir
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_fep_dynamic  --  dynamics undisturbed  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_fep_dynamic" checks that a short trajectory is the same
c     whether or not trial energies are found between its steps, with
c     and without neighbor lists, and when the induced dipole solver
c     starts from zero instead of from the direct field dipoles
c
c
      subroutine test_fep_dynamic
      implicit none
c
c
      call fepdyn ('polar-eps 0.00001','default')
      call fepdyn ('neighbor-list','lists')
      call fepdyn ('pcg-noguess','noguess')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine fepdyn  --  one pair of short trajectories  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "fepdyn" runs the same short trajectory twice from rest, the
c     second time finding trial energies after every step, and then
c     compares the final coordinates and velocities of the two runs
c
c
      subroutine fepdyn (addkey,tag)
      use atoms
      use moldyn
      implicit none
      integer i,j
      integer irun
      integer nat
      real*8 dxmax,dvmax
      real*8, allocatable :: xref(:,:)
      real*8, allocatable :: vref(:,:)
      character*(*) addkey,tag
c
c
      call pushdir ('file/mutate')
      do irun = 1, 2
         call loadfix_keyadd ('ionwat','179_ionwat_ast_l05.key',addkey)
         call fepdynrun (irun.eq.2)
         if (irun .eq. 1) then
            nat = n
            allocate (xref(3,nat))
            allocate (vref(3,nat))
            do i = 1, nat
               xref(1,i) = x(i)
               xref(2,i) = y(i)
               xref(3,i) = z(i)
               do j = 1, 3
                  vref(j,i) = v(j,i)
               end do
            end do
         else
            dxmax = 0.0d0
            dvmax = 0.0d0
            do i = 1, nat
               dxmax = max(dxmax,abs(x(i)-xref(1,i)))
               dxmax = max(dxmax,abs(y(i)-xref(2,i)))
               dxmax = max(dxmax,abs(z(i)-xref(3,i)))
               do j = 1, 3
                  dvmax = max(dvmax,abs(v(j,i)-vref(j,i)))
               end do
            end do
            call assert_real (dxmax,0.0d0,1.0d-10,
     &                        'fep dynamic '//tag//' coordinates')
            call assert_real (dvmax,0.0d0,1.0d-8,
     &                        'fep dynamic '//tag//' velocities')
         end if
         call sd_final
      end do
      deallocate (xref)
      deallocate (vref)
      call popdir
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine fepdynrun  --  short trajectory in process  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "fepdynrun" takes twenty Beeman steps from rest for the system
c     that is loaded, with the output of the dynamics silenced, and
c     finds trial energies after each step when "observe" is set
c
c
      subroutine fepdynrun (observe)
      use inform
      use iounit
      implicit none
      integer istep
      integer isave,freeunit
      real*8 dt
      real*8 enative
      real*8 ltrial(2)
      real*8 etrial(2)
      logical observe
c
c
c     set up the dynamics, and keep it from saving or printing
c
      call sd_mdinit (0.0d0)
      iwrite = 1000000
      iprint = 1000000
      dt = 0.001d0
      ltrial(1) = 0.4d0
      ltrial(2) = 0.6d0
      isave = iout
      iout = freeunit ()
      open (unit=iout,status='scratch')
      do istep = 1, 20
         call beeman (istep,dt)
         if (observe)  call efeptrial (2,ltrial,etrial,enative)
      end do
      close (unit=iout)
      iout = isave
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_fep_sample  --  sampler bookkeeping  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_fep_sample" checks the settings made by the free energy
c     perturbation mode, the dynamics steps at which trial energies
c     are found, the order and the lambda values of the records, and
c     that each record is written to the output file exactly once
c
c
      subroutine test_fep_sample
      use dlmda
      use fep
      use files
      use mutant
      use uprior
      implicit none
      integer i
      integer istep
      integer irow
      integer nrow
      integer oldleng
      real*8 eps
      real*8 simlmda,trial,ene
      character*240 oldname
c
c
c     the mode owns the main lambda without the lambda derivative
c
      eps = 1.0d-8
      call pushdir ('file/fep')
      call loadfix ('../mutate/water2','window.key')
      call assert_logical (use_fep,.true.,'fep sample use_fep')
      call assert_logical (use_mainlmda,.true.,'fep sample mainlmda')
      call assert_logical (use_dlmda,.false.,'fep sample no deriv')
      call assert_int (nlmdawin,3,'fep sample window count')
      call assert_int (fepintv,4,'fep sample interval')
      call assert_real (lmdawinratio,0.5d0,0.0d0,'fep sample ratio')
      call assert_real (lambda,1.0d0,0.0d0,'fep sample first lambda')
c
c     sixty steps give three windows of twenty steps, each with ten
c     production steps and so two samples; the end windows have one
c     neighbor and the middle window has two
c
      call initfepdyn (60)
      call assert_int (nfeptot,14,'fep sample capacity')
      do istep = 1, 60
         call lmdawindyn (istep)
         call lmdawinstep (istep)
      end do
      call assert_int (nfep,14,'fep sample records')
c
c     the first window samples at steps 14 and 18, native first
c
      call assert_int (fepstep(1),14,'fep sample first step')
      call assert_real (feplmda(1),1.0d0,0.0d0,'fep sample 1 sim')
      call assert_real (feptrial(1),1.0d0,0.0d0,'fep sample 1 trial')
      call assert_int (fepstep(2),14,'fep sample second row step')
      call assert_real (feplmda(2),1.0d0,0.0d0,'fep sample 2 sim')
      call assert_real (feptrial(2),0.5d0,0.0d0,'fep sample 2 trial')
      call assert_int (fepstep(3),18,'fep sample second step')
c
c     the middle window lists the previous window, then the next
c
      call assert_int (fepstep(5),34,'fep sample middle step')
      call assert_real (feplmda(5),0.5d0,0.0d0,'fep sample 5 sim')
      call assert_real (feptrial(5),0.5d0,0.0d0,'fep sample 5 trial')
      call assert_real (feptrial(6),1.0d0,0.0d0,'fep sample 6 trial')
      call assert_real (feptrial(7),0.0d0,0.0d0,'fep sample 7 trial')
      call assert_real (feplmda(7),0.5d0,0.0d0,'fep sample 7 sim')
c
c     the last window samples at steps 54 and 58
c
      call assert_int (fepstep(13),58,'fep sample last step')
      call assert_real (feplmda(13),0.0d0,0.0d0,'fep sample 13 sim')
      call assert_real (feptrial(14),0.5d0,0.0d0,'fep sample 14 trial')
c
c     the structure never moved, so a trial energy found from one
c     window is the native energy found in the neighboring window
c
      call assert_real (fepene(3),fepene(1),eps,'fep sample repeat')
      call assert_real (fepene(2),fepene(5),eps,'fep sample down')
      call assert_real (fepene(6),fepene(1),eps,'fep sample up')
      call assert_real (fepene(7),fepene(11),eps,'fep sample last')
c
c     write the records into a file of our own in two parts, and
c     check that a call with nothing new adds no rows
c
      oldname = filename
      oldleng = leng
      filename = 'fepsave_tmp'
      leng = 11
      call tiwipe ('fepsave_tmp.fep')
      nfep = 5
      call prtfephead
      call fepcount ('fepsave_tmp.fep',nrow)
      call assert_int (nrow,0,'fep sample header only')
      call savefep
      call fepcount ('fepsave_tmp.fep',nrow)
      call assert_int (nrow,5,'fep sample first rows')
      nfep = 14
      call savefep
      call savefep
      call fepcount ('fepsave_tmp.fep',nrow)
      call assert_int (nrow,14,'fep sample all rows once')
      irow = 6
      call feprow ('fepsave_tmp.fep',irow,i,simlmda,trial,ene)
      call assert_int (i,34,'fep sample row step')
      call assert_real (simlmda,0.5d0,eps,'fep sample row sim')
      call assert_real (trial,1.0d0,eps,'fep sample row trial')
      call assert_real (ene,fepene(6),eps,'fep sample row energy')
      call tiwipe ('fepsave_tmp.fep')
      filename = oldname
      leng = oldleng
c
c     an ascending schedule lists its neighbors the same way
c
      lmdawinlist(1) = 0.0d0
      lmdawinlist(3) = 1.0d0
      lambda = lmdawinlist(1)
      call initfepdyn (60)
      do istep = 1, 60
         call lmdawindyn (istep)
         call lmdawinstep (istep)
      end do
      call assert_int (nfep,14,'fep sample ascending records')
      call assert_real (feplmda(1),0.0d0,0.0d0,'fep sample up 1 sim')
      call assert_real (feptrial(2),0.5d0,0.0d0,'fep sample up 2 trial')
      call assert_real (feptrial(6),0.0d0,0.0d0,'fep sample up 6 trial')
      call assert_real (feptrial(7),1.0d0,0.0d0,'fep sample up 7 trial')
      call assert_real (feplmda(14),1.0d0,0.0d0,'fep sample up 14 sim')
c
c     the dynamics setup turns off the predictor the key asks for
c
      call sd_mdinit (0.0d0)
      call assert_logical (use_pred,.false.,'fep sample no predictor')
      call sd_final
      call popdir
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_fep_run  --  energies of a dynamics run  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_fep_run" runs dynamics over three lambda windows with
c     free energy perturbation in use, and checks each energy in the
c     output file against the saved frame it was found from, set up
c     afresh at the trial lambda, at constant volume and pressure
c
c
      subroutine test_fep_run
      implicit none
c
c
      call feprun ('water2 60 0.1 0.001 2 298','digits 8','nvt')
      call feprun ('water2 60 0.1 0.001 4 298 1.0',
     &             'barostat montecarlo','npt')
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine feprun  --  check one dynamics run  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "feprun" makes one dynamics run in a scratch directory, with
c     windows of twenty steps that each take a single sample on
c     their final step, where a frame is also saved, and compares
c     the seven energies written against those of the saved frames
c
c
      subroutine feprun (args,addkey,tag)
      implicit none
      integer i
      integer ist
      integer istep
      integer nrow
      integer iwin(7)
      real*8 energy
      real*8 eref
      real*8 simlmda,trial,ene
      real*8 sim(7)
      real*8 tri(7)
      character*1 quote
      character*(*) args,addkey,tag
      character*40 lmdakey
      character*60 label
      character*240 cmd
c
c
c     the window, and both lambda values, that each row should hold
c
      data iwin  / 1, 1, 2, 2, 2, 3, 3 /
      data sim  / 1.0d0, 1.0d0, 0.5d0, 0.5d0, 0.5d0, 0.0d0, 0.0d0 /
      data tri  / 1.0d0, 0.5d0, 0.5d0, 1.0d0, 0.0d0, 0.0d0, 0.5d0 /
c
c     run the dynamics with a sample at the end of every window
c
      quote = char(39)
      call fepprep ('fep_run')
      call pushdir ('file/fep_run')
      call tnist_append ('water2.key','lambda-eqratio 0.5')
      call tnist_append ('water2.key','lambda-interval 10')
      call tnist_append ('water2.key','digits 8')
      call tnist_append ('water2.key',addkey)
      call run_prog ('dynamic',args,'out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'fep run '//tag//' status')
         call fepcount ('water2.fep',nrow)
         call assert_int (nrow,7,'fep run '//tag//' rows')
c
c     the header carries the temperature given to the dynamics
c
         cmd = 'grep -q "Simulation Temperature  *298.00" water2.fep'
         call execute_command_line (cmd,exitstat=ist)
         call assert_int (ist,0,'fep run '//tag//' temperature')
c
c     set each saved frame up afresh at the trial lambda of the row
c
         do i = 1, min(nrow,7)
            call feprow ('water2.fep',i,istep,simlmda,trial,ene)
            write (label,10)  tag,i
   10       format ('fep run ',a,' row ',i0)
            call assert_int (istep,20*iwin(i),trim(label)//' step')
            call assert_real (simlmda,sim(i),1.0d-8,
     &                        trim(label)//' sim lambda')
            call assert_real (trial,tri(i),1.0d-8,
     &                        trim(label)//' trial lambda')
            write (cmd,20)  2*iwin(i),quote,quote
   20       format ('awk -v k=',i0,1x,a1,'/Water/{f++} f==k',a1,
     &                 ' water2.arc > frame.xyz')
            call execute_command_line (cmd)
            write (lmdakey,30)  trial
   30       format ('lambda ',f12.8)
            call loadfix_keyadd ('frame','ref.key',lmdakey)
            eref = energy ()
            call assert_real (ene,eref,1.0d-5,trim(label)//' energy')
            call final
         end do
      end if
      call popdir
      call tnist_clean ('fep_run')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_fep_stop  --  requested stop records  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_fep_stop" runs dynamics that is asked to stop at its
c     first trajectory save, and checks that the samples taken up to
c     and on that final step have all reached the output file
c
c
      subroutine test_fep_stop
      implicit none
      integer ist
      integer nrow
      integer iend
      integer freeunit
c
c
c     with no equilibration and samples every five steps, the stop
c     at step ten follows two samples of two energies each
c
      call fepprep ('fep_stop')
      call pushdir ('file/fep_stop')
      call tnist_append ('water2.key','lambda-eqratio 0.0')
      call tnist_append ('water2.key','lambda-interval 5')
      iend = freeunit ()
      open (unit=iend,file='water2.end',status='new')
      close (unit=iend)
      call run_prog ('dynamic','water2 60 0.1 0.001 2 298',
     &               'out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'fep stop status')
         call fepcount ('water2.fep',nrow)
         call assert_int (nrow,4,'fep stop rows at stop step')
      end if
      call popdir
      call tnist_clean ('fep_stop')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine fepprep  --  fixture for an FEP dynamics run  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "fepprep" creates a scratch directory holding a small water
c     box and two keyfiles; "ref.key" has the force field and the
c     sublambda maps and takes its box from each structure, while
c     "water2.key" adds the box and three lambda windows for a
c     free energy perturbation run
c
c
      subroutine fepprep (work)
      implicit none
      character*1 quote
      character*(*) work
      character*512 cmd
c
c
c     copy the structure, and its keyfile less the box dimension
c
      quote = char(39)
      call pushdir ('file/mutate')
      cmd = 'rm -rf ../'//trim(work)//' ; mkdir -p ../'//trim(work)//
     &      ' ; cp water2.xyz ../'//trim(work)//'/ ; sed '//quote//
     &      '/^lambda-deriv/d;/^a-axis/d'//quote//
     &      ' 164_water_ast_ne_mcut_d1_l05.key > ../'//trim(work)//
     &      '/ref.key ; cp ../'//trim(work)//'/ref.key ../'//
     &      trim(work)//'/water2.key'
      call execute_command_line (cmd)
      call popdir
c
c     add the box and the lambda window schedule for the dynamics
c
      call pushdir ('file/'//trim(work))
      call tnist_append ('water2.key','a-axis 18.643')
      call tnist_append ('water2.key','integrator verlet')
      call tnist_append ('water2.key','lambda-mode fep')
      call tnist_append ('water2.key','lambda-window 1.0')
      call tnist_append ('water2.key','lambda-window 0.5')
      call tnist_append ('water2.key','lambda-window 0.0')
      call popdir
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine fepcount  --  count FEP file data records  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "fepcount" returns the number of data records in a trial
c     energy file, counting only the nonblank records that follow
c     the line labeling the columns of the energies
c
c
      subroutine fepcount (fname,nrow)
      implicit none
      integer nrow
      integer ifep
      integer freeunit
      logical header
      character*(*) fname
      character*240 record
c
c
      nrow = 0
      header = .false.
      ifep = freeunit ()
      open (unit=ifep,file=fname,status='old')
   10 continue
      read (ifep,20,end=30)  record
   20 format (a240)
      if (.not. header) then
         if (index(record,'Sim Lambda') .ne. 0)  header = .true.
      else if (record .ne. ' ') then
         nrow = nrow + 1
      end if
      goto 10
   30 continue
      close (unit=ifep)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine feprow  --  values of one FEP file record  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "feprow" returns the step, the simulation and trial lambda
c     values and the energy held in the requested data record of a
c     trial energy file; the step is negative if there is no such
c     record
c
c
      subroutine feprow (fname,irow,istep,simlmda,trial,ene)
      implicit none
      integer irow
      integer istep
      integer nrow
      integer ifep
      integer freeunit
      real*8 simlmda,trial,ene
      logical header
      character*(*) fname
      character*240 record
c
c
      istep = -1
      simlmda = -1.0d0
      trial = -1.0d0
      ene = 0.0d0
      nrow = 0
      header = .false.
      ifep = freeunit ()
      open (unit=ifep,file=fname,status='old')
   10 continue
      read (ifep,20,end=30)  record
   20 format (a240)
      if (.not. header) then
         if (index(record,'Sim Lambda') .ne. 0)  header = .true.
      else if (record .ne. ' ') then
         nrow = nrow + 1
         if (nrow .eq. irow) then
            read (record,*,err=30,end=30)  istep,simlmda,trial,ene
            goto 30
         end if
      end if
      goto 10
   30 continue
      close (unit=ifep)
      return
      end

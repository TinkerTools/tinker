c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################
c     ##                                              ##
c     ##  subroutine test_mutate  --  mutation tests  ##
c     ##                                              ##
c     ##################################################
c
c
c     "test_mutate" checks AMOEBA mutation energy, gradient, virial
c     and named energy component regressions
c
c
      subroutine test_mutate
      implicit none
c
c
      call test_mutate_refresh
      call test_mutate_lmda
      call test_mutate_lmdamode
      call test_mutate_mv
      call test_mutate_mp
      call test_mutate_ast
      call test_mutate_adt
      call test_mutate_qnt
      call test_mutate_exp
      call test_mutate_inv
      call test_mutate_exf
      call test_mutate_emplar
      call test_mutate_qntrng
      call test_mutate_rels
      call test_mutate_lmdafix
      call test_mutate_lmdadrv
      call test_mutate_lmdapin
      call test_mutate_relsmap
      call test_mutate_apm
      call test_mutate_vsoft
      call test_mutate_gate
      call test_mutate_chiral
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_refresh  --  lambda state refresh  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_refresh" advances a TI window while leaving the old
c     mapped parameter arrays installed, then checks that energy,
c     analysis and gradient refresh the complete lambda state before use
c
c
      subroutine test_mutate_refresh
      use atoms
      use dlmda
      use energi
      use mutant
      use thrmint
      implicit none
      real*8 energy
      real*8 e,eref,earef
      real*8 eoneref,emoneref
      real*8 emref,epref,evref
      real*8 emaref,eparef,evaref
      real*8 dedlref,d2edlref
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: derivsref(:,:)
      logical skiptest
c
c
      if (skiptest('test_mutate_refresh','mutate'))  return
      call pushdir ('file/mutate')
      call loadfix ('water2','085_water_exp_ast_l05.key')
      allocate (derivs(3,n))
      allocate (derivsref(3,n))
c
c     replace OST ownership with a deterministic three-window TI schedule
c
      use_ost = .false.
      use_ostdyn = .false.
      use_meta = .false.
      use_metadyn = .false.
      use_ti = .true.
      use_mainlmda = .true.
      tinbin = 3
      tibin = 1
      tieqratio = 0.0d0
      tinstepavg = 1
      if (allocated(tilmdalist))  deallocate (tilmdalist)
      if (allocated(tiwinend))  deallocate (tiwinend)
      allocate (tilmdalist(tinbin))
      allocate (tiwinend(tinbin))
      tilmdalist(1) = 0.50d0
      tilmdalist(2) = 0.25d0
      tilmdalist(3) = 1.00d0
      tiwinend(1) = 10
      tiwinend(2) = 20
      tiwinend(3) = 30
      lambda = tilmdalist(2)
      call mapsublmda (lambda)
      call altelec
      eref = energy ()
      emref = em
      epref = ep
      evref = ev
      call analysis (earef)
      emaref = em
      eparef = ep
      evaref = ev
      call gradient (eref,derivsref)
      dedlref = dedl
      d2edlref = d2edl2
c
c     install the first-window arrays, advance only the authoritative TI
c     lambda, and require energy to rebuild all dependent parameter state
c
      tibin = 1
      lambda = tilmdalist(1)
      call mapsublmda (lambda)
      call altelec
      call tischedule
      call assert_real (lambda,0.25d0,0.0d0,
     &                  'energy refresh TI window')
      call assert_real (elambda,0.50d0,0.0d0,
     &                  'energy refresh stale elambda')
      e = energy ()
      call assert_real (elambda,0.25d0,0.0d0,
     &                  'energy refresh elambda')
      call assert_real (e,eref,1.0d-10,
     &                  'energy refresh total')
      call assert_real (em,emref,1.0d-10,
     &                  'energy refresh multipole')
      call assert_real (ep,epref,1.0d-10,
     &                  'energy refresh polarization')
      call assert_real (ev,evref,1.0d-10,
     &                  'energy refresh van der Waals')
c
c     returning from a fractional state to one must restore originals
c
      lambda = tilmdalist(3)
      call mapsublmda (lambda)
      call altelec
      eoneref = energy ()
      emoneref = em
      tibin = 2
      lambda = tilmdalist(2)
      call mapsublmda (lambda)
      call altelec
      call tischedule
      e = energy ()
      call assert_real (elambda,1.0d0,0.0d0,
     &                  'energy refresh endpoint elambda')
      call assert_real (e,eoneref,1.0d-10,
     &                  'energy refresh endpoint total')
      call assert_real (em,emoneref,1.0d-10,
     &                  'energy refresh endpoint multipole')
c
c     repeat with stale first-window arrays at the analysis boundary
c
      tibin = 1
      lambda = tilmdalist(1)
      call mapsublmda (lambda)
      call altelec
      call tischedule
      call analysis (e)
      call assert_real (e,earef,1.0d-10,
     &                  'analysis refresh total')
      call assert_real (em,emaref,1.0d-10,
     &                  'analysis refresh multipole')
      call assert_real (ep,eparef,1.0d-10,
     &                  'analysis refresh polarization')
      call assert_real (ev,evaref,1.0d-10,
     &                  'analysis refresh van der Waals')
c
c     poison the mapped scalars as well as leaving first-window arrays
c     installed, then require gradient to refresh values and derivatives
c
      tibin = 1
      lambda = tilmdalist(1)
      call mapsublmda (lambda)
      call altelec
      call tischedule
      elambda = -1.0d0
      plambda = -1.0d0
      vlambda = -1.0d0
      deldlmda = -1.0d0
      dpldlmda = -1.0d0
      dvldlmda = -1.0d0
      call gradient (e,derivs)
      call assert_real (elambda,0.25d0,0.0d0,
     &                  'gradient refresh elambda')
      call assert_real (plambda,0.0625d0,1.0d-15,
     &                  'gradient refresh plambda')
      call assert_real (vlambda,0.015625d0,1.0d-15,
     &                  'gradient refresh vlambda')
      call assert_real (e,eref,1.0d-10,
     &                  'gradient refresh total')
      call assert_real (dedl,dedlref,1.0d-10,
     &                  'gradient refresh dE/dL')
      call assert_real (d2edl2,d2edlref,1.0d-10,
     &                  'gradient refresh d2E/dL2')
      call assert_grad (derivs,derivsref,n,1.0d-10,
     &                  'gradient refresh Cartesian gradient')
c
c     clean up the molecular fixture and schedule storage
c
      deallocate (derivs)
      deallocate (derivsref)
      call popdir
      call final
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_mutate_lmda  --  main lambda keyword  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_mutate_lmda" checks that the "LAMBDA" keyword drives the
c     sublambdas that name a map; with no map keyword at all it drives
c     all three linearly, while a map keyword restricts it to the terms
c     named and leaves the rest fully coupled at one; each mapped
c     fixture is compared against the explicit sublambda fixture it
c     should match
c
c
      subroutine test_mutate_lmda
      use atoms
      use dlmda
      use mutant
      implicit none
      real*8 energy
      real*8 e,eref
      logical skiptest
c
c
      if (skiptest('test_mutate_lmda','mutate'))  return
      call pushdir ('file/mutate')
c
c     a main lambda with no map keyword drives all three sublambdas
c
      call loadfix ('water2','146_water_lmda_ast_l05.key')
      call assert_logical (use_mainlmda,.true.,
     &                  'lmda linear use_mainlmda')
      call assert_logical (use_elmdamap,.true.,
     &                  'lmda linear use_elmdamap')
      call assert_logical (use_plmdamap,.true.,
     &                  'lmda linear use_plmdamap')
      call assert_logical (use_vlmdamap,.true.,
     &                  'lmda linear use_vlmdamap')
      call assert_real (elambda,0.5d0,0.0d0,'lmda linear elambda')
      call assert_real (plambda,0.5d0,0.0d0,'lmda linear plambda')
      call assert_real (vlambda,0.5d0,0.0d0,'lmda linear vlambda')
      call assert_real (deldlmda,1.0d0,0.0d0,'lmda linear deldlmda')
      call assert_real (dpldlmda,1.0d0,0.0d0,'lmda linear dpldlmda')
      call assert_real (dvldlmda,1.0d0,0.0d0,'lmda linear dvldlmda')
      e = energy ()
      call final
      call loadfix ('water2','147_water_lmda_ast_e05.key')
      eref = energy ()
      call assert_real (e,eref,1.0d-10,'lmda linear energy')
      call final
c
c     the fully coupled endpoint matches a fixture with no lambda
c
      call loadfix ('water2','148_water_lmda_ast_l10.key')
      call assert_real (elambda,1.0d0,0.0d0,'lmda endpoint elambda')
      call assert_real (plambda,1.0d0,0.0d0,'lmda endpoint plambda')
      call assert_real (vlambda,1.0d0,0.0d0,'lmda endpoint vlambda')
      call assert_real (scalphav,0.0d0,0.0d0,
     &                  'lmda endpoint scalphav')
      call assert_real (scexp,2.0d0,0.0d0,'lmda endpoint scexp')
      e = energy ()
      call final
      call loadfix ('water2','149_water_lmda_ast_none.key')
      eref = energy ()
      call assert_real (e,eref,1.0d-10,'lmda endpoint energy')
      call final
c
c     an explicit electrostatic map drives electrostatics alone, with
c     polarization no longer following it
c
      call loadfix ('water2','150_water_lmda_qnt_l05.key')
      call assert_logical (use_elmdamap,.true.,
     &                  'lmda ele only use_elmdamap')
      call assert_logical (use_plmdamap,.false.,
     &                  'lmda ele only use_plmdamap')
      call assert_logical (use_vlmdamap,.false.,
     &                  'lmda ele only use_vlmdamap')
      call assert_real (elambda,0.5d0,1.0d-15,'lmda ele only elambda')
      call assert_real (plambda,1.0d0,0.0d0,'lmda ele only plambda')
      call assert_real (vlambda,1.0d0,0.0d0,'lmda ele only vlambda')
      call assert_real (deldlmda,1.875d0,1.0d-15,
     &                  'lmda ele only deldlmda')
      call assert_real (dpldlmda,0.0d0,0.0d0,'lmda ele only dpldlmda')
      call assert_real (dvldlmda,0.0d0,0.0d0,'lmda ele only dvldlmda')
      call final
c
c     each sublambda follows the map named for it, so the three reach
c     different values with different chain factors; the polarization
c     window starts at the main lambda, leaving it on the flat side of
c     its taper and fully decoupled
c
      call loadfix ('water2','151_water_lmda_vexp_l05.key')
      call assert_logical (vlmdamap.eq.'EXP',.true.,
     &                  'lmda mixed vlmdamap')
      call assert_real (elambda,0.31744d0,1.0d-12,
     &                  'lmda mixed elambda')
      call assert_real (plambda,0.0d0,0.0d0,'lmda mixed plambda')
      call assert_real (vlambda,0.125d0,1.0d-15,'lmda mixed vlambda')
      call assert_real (deldlmda,3.456d0,1.0d-12,
     &                  'lmda mixed deldlmda')
      call assert_real (dpldlmda,0.0d0,0.0d0,'lmda mixed dpldlmda')
      call assert_real (dvldlmda,0.75d0,1.0d-15,'lmda mixed dvldlmda')
      call final
c
c     naming two maps drives those terms and leaves the third fully
c     coupled, matching the same state set by explicit sublambdas
c
      call loadfix ('water2','152_water_lmda_mp05.key')
      call assert_logical (use_elmdamap,.true.,
     &                  'lmda pair use_elmdamap')
      call assert_logical (use_plmdamap,.true.,
     &                  'lmda pair use_plmdamap')
      call assert_logical (use_vlmdamap,.false.,
     &                  'lmda pair use_vlmdamap')
      call assert_real (elambda,0.5d0,1.0d-15,'lmda pair elambda')
      call assert_real (plambda,0.5d0,1.0d-15,'lmda pair plambda')
      call assert_real (vlambda,1.0d0,0.0d0,'lmda pair vlambda')
      call assert_real (deldlmda,1.0d0,1.0d-15,'lmda pair deldlmda')
      call assert_real (dpldlmda,1.0d0,1.0d-15,'lmda pair dpldlmda')
      call assert_real (dvldlmda,0.0d0,0.0d0,'lmda pair dvldlmda')
      e = energy ()
      call final
      call loadfix ('water2','153_water_lmda_mp05_expl.key')
      call assert_real (vlambda,1.0d0,0.0d0,'lmda pair expl vlambda')
      eref = energy ()
      call assert_real (e,eref,1.0d-10,'lmda pair energy')
      call final
      call popdir
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_lmdamode  --  lambda sampling mode  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_lmdamode" checks that the "LAMBDA-MODE" keyword
c     turns on the sampling method it names along with the lambda
c     derivative and main lambda, and that an unknown mode leaves
c     every sampling method off
c
c
      subroutine test_mutate_lmdamode
      use dlmda
      use mutant
      use ost
      implicit none
      logical skiptest
c
c
      if (skiptest('test_mutate_lmdamode','mutate'))  return
      call pushdir ('file/mutate')
c
c     an abf mode turns on abf alone and drives the main lambda
c
      call loadfix_keyadd ('water2','150_water_lmda_qnt_l05.key',
     &                     'LAMBDA-MODE abf')
      call assert_logical (lmdasampmode.eq.'ABF',.true.,
     &                     'lmdamode abf lmdasampmode')
      call assert_logical (use_abf,.true.,'lmdamode abf use_abf')
      call assert_logical (use_ost,.false.,'lmdamode abf use_ost')
      call assert_logical (use_meta,.false.,'lmdamode abf use_meta')
      call assert_logical (use_ti,.false.,'lmdamode abf use_ti')
      call assert_logical (use_dlmda,.true.,'lmdamode abf use_dlmda')
      call assert_logical (use_mainlmda,.true.,
     &                     'lmdamode abf use_mainlmda')
c
c     abf allocates only its samples, interval lists and lambda bins
c
      call assert_logical (allocated(lmdaihist),.true.,
     &                     'lmdamode abf allocates lmdaihist')
      call assert_logical (allocated(lmdalhist),.true.,
     &                     'lmdamode abf allocates lmdalhist')
      call assert_logical (allocated(lmdafhist),.true.,
     &                     'lmdamode abf allocates lmdafhist')
      call assert_logical (allocated(lmdallist),.true.,
     &                     'lmdamode abf allocates lmdallist')
      call assert_logical (allocated(lmdaflist),.true.,
     &                     'lmdamode abf allocates lmdaflist')
      call assert_logical (allocated(lmdafmean),.true.,
     &                     'lmdamode abf allocates lmdafmean')
      call assert_logical (allocated(lmdafsum),.true.,
     &                     'lmdamode abf allocates lmdafsum')
      call assert_logical (allocated(lmdafwt),.true.,
     &                     'lmdamode abf allocates lmdafwt')
      call assert_logical (allocated(osthist),.false.,
     &                     'lmdamode abf allocates no osthist')
      call assert_logical (allocated(osthead),.false.,
     &                     'lmdamode abf allocates no osthead')
      call assert_logical (allocated(osthhist),.false.,
     &                     'lmdamode abf allocates no osthhist')
      call assert_logical (allocated(gkernel),.false.,
     &                     'lmdamode abf allocates no gkernel')
      call assert_logical (allocated(vkernelmax),.false.,
     &                     'lmdamode abf allocates no vkernelmax')
      call final
c
c     ost leaves pinned polarization in absolute single topology
c
      call loadfix_keyadd ('water2','170_water_lmda_ast_epin_l05.key',
     &                     'LAMBDA-MODE ost')
      call assert_logical (use_ost,.true.,'lmdamode ost use_ost')
      call assert_logical (use_pdlmda,.false.,
     &                     'lmdamode ost pinned use_pdlmda')
      call assert_logical (use_epdt,.false.,
     &                     'lmdamode ost pinned use_epdt')
      call assert_logical (use_prst,.true.,
     &                     'lmdamode ost pinned use_prst')
      call final
c
c     ost uses dual topology when polarization follows the main lambda
c
      call loadfix_keyadd ('water2','152_water_lmda_mp05.key',
     &                     'LAMBDA-MODE ost')
      call assert_logical (use_pdlmda,.true.,
     &                     'lmdamode ost mapped use_pdlmda')
      call assert_logical (use_epdt,.true.,
     &                     'lmdamode ost mapped use_epdt')
      call assert_logical (use_prst,.false.,
     &                     'lmdamode ost mapped use_prst')
      call final
c
c     the first lambda derivative alone keeps polarization on one state
c
      call loadfix_keyadd ('water2','152_water_lmda_mp05.key',
     &                     'LAMBDA-DERIV')
      call assert_logical (use_dlmda,.true.,'lmdaderiv use_dlmda')
      call assert_logical (use_d2lmda,.false.,'lmdaderiv use_d2lmda')
      call assert_logical (use_epdt,.false.,'lmdaderiv use_epdt')
      call assert_logical (use_prst,.true.,'lmdaderiv use_prst')
      call final
c
c     second lambda derivatives need dual topology polarization
c
      call loadfix_keyadd ('water2','152_water_lmda_mp05.key',
     &                     'LAMBDA-DERIV2')
      call assert_logical (use_dlmda,.true.,'lmdaderiv2 use_dlmda')
      call assert_logical (use_d2lmda,.true.,'lmdaderiv2 use_d2lmda')
      call assert_logical (use_epdt,.true.,'lmdaderiv2 use_epdt')
      call assert_logical (use_prst,.false.,'lmdaderiv2 use_prst')
      call final
c
c     the staged relative legs choose the polarization path the same way
c
      call loadfix ('water2','203_water_rels_st_l085.key')
      call assert_logical (use_rel,.true.,'rels deriv use_rel')
      call assert_logical (use_epdt,.false.,'rels deriv use_epdt')
      call assert_logical (use_prst,.true.,'rels deriv use_prst')
      call final
      call loadfix ('water2','136_water_rels_ye_l085.key')
      call assert_logical (use_rel,.true.,'rels deriv2 use_rel')
      call assert_logical (use_epdt,.true.,'rels deriv2 use_epdt')
      call assert_logical (use_prst,.false.,'rels deriv2 use_prst')
      call final
c
c     an unknown mode leaves every sampling method off
c
      call loadfix_keyadd ('water2','150_water_lmda_qnt_l05.key',
     &                     'LAMBDA-MODE bogus')
      call assert_logical (use_abf,.false.,'lmdamode bogus use_abf')
      call assert_logical (use_ost,.false.,'lmdamode bogus use_ost')
      call assert_logical (use_meta,.false.,'lmdamode bogus use_meta')
      call assert_logical (use_ti,.false.,'lmdamode bogus use_ti')
      call assert_logical (use_dlmda,.false.,
     &                     'lmdamode bogus use_dlmda')
      call assert_logical (use_lmdacv,.false.,
     &                     'lmdamode convcri gate off by default')
      call final
c
c     a meta mode turns on metadynamics alone
c
      call loadfix_keyadd ('water2','150_water_lmda_qnt_l05.key',
     &                     'LAMBDA-MODE meta')
      call assert_logical (use_meta,.true.,'lmdamode meta use_meta')
      call assert_logical (use_ost,.false.,'lmdamode meta use_ost')
      call final
c
c     the convergence gate reads the deviation then the ratio
c
      call loadfix_keyadd ('water2','150_water_lmda_qnt_l05.key',
     &                     'LAMBDA-CONVCRI 5.0 0.2')
      call assert_logical (use_lmdacv,.true.,
     &                     'lmdamode convcri use_lmdacv')
      call assert_real (lmdacvstd,5.0d0,0.0d0,
     &                  'lmdamode convcri lmdacvstd')
      call assert_real (lmdacvrat,0.2d0,0.0d0,
     &                  'lmdamode convcri lmdacvrat')
      call final
c
c     a bare convergence keyword keeps the default tolerances
c
      call loadfix_keyadd ('water2','150_water_lmda_qnt_l05.key',
     &                     'LAMBDA-CONVCRI')
      call assert_logical (use_lmdacv,.true.,
     &                     'lmdamode bare convcri use_lmdacv')
      call assert_real (lmdacvstd,50.0d0,0.0d0,
     &                  'lmdamode bare convcri lmdacvstd')
      call assert_real (lmdacvrat,0.5d0,0.0d0,
     &                  'lmdamode bare convcri lmdacvrat')
      call final
      call popdir
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_mv  --  mutation lambda-scan cases  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_mv" runs the nine water mutation fixtures; cases
c     001-006 isolate electrostatics at three ele-lambda values with
c     Ewald on and off, while cases 007-009 isolate van der Waals at
c     three vdw-lambda values
c
c
      subroutine test_mutate_mv
      implicit none
c
c
      call test_mutate_fixed ('001_water_ye_m10.key',
     &   '001_water_ye_m10.txt','001_water_ye_m10',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('002_water_ne_m10.key',
     &   '002_water_ne_m10.txt','002_water_ne_m10',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('003_water_ye_m05.key',
     &   '003_water_ye_m05.txt','003_water_ye_m05',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('004_water_ne_m05.key',
     &   '004_water_ne_m05.txt','004_water_ne_m05',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('005_water_ye_m00.key',
     &   '005_water_ye_m00.txt','005_water_ye_m00',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('006_water_ne_m00.key',
     &   '006_water_ne_m00.txt','006_water_ne_m00',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('007_water_v10.key',
     &   '007_water_v10.txt','007_water_v10',
     &   .false., .false., .true.,  .true.)
      call test_mutate_fixed ('008_water_v05.key',
     &   '008_water_v05.txt','008_water_v05',
     &   .false., .false., .true.,  .true.)
      call test_mutate_fixed ('009_water_v00.key',
     &   '009_water_v00.txt','009_water_v00',
     &   .false., .false., .true.,  .true.)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_mp  --  electrostatic lambda cases  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_mp" runs the twenty water mutation fixtures that
c     scan the electrostatic lambda values; cases 010-015 keep only the
c     multipole term at three ele-lambda values with Ewald on and off,
c     cases 016-021 keep only the polarization term at three pol-lambda
c     values, and cases 022-029 leave both terms active while varying
c     ele-lambda and pol-lambda together
c
c
      subroutine test_mutate_mp
      implicit none
c
c
      call test_mutate_fixed ('010_water_ye_m10.key',
     &   '010_water_ye_m10.txt','010_water_ye_m10',
     &   .true.,  .false., .false., .true.)
      call test_mutate_fixed ('011_water_ne_m10.key',
     &   '011_water_ne_m10.txt','011_water_ne_m10',
     &   .true.,  .false., .false., .false.)
      call test_mutate_fixed ('012_water_ye_m05.key',
     &   '012_water_ye_m05.txt','012_water_ye_m05',
     &   .true.,  .false., .false., .true.)
      call test_mutate_fixed ('013_water_ne_m05.key',
     &   '013_water_ne_m05.txt','013_water_ne_m05',
     &   .true.,  .false., .false., .false.)
      call test_mutate_fixed ('014_water_ye_m00.key',
     &   '014_water_ye_m00.txt','014_water_ye_m00',
     &   .true.,  .false., .false., .true.)
      call test_mutate_fixed ('015_water_ne_m00.key',
     &   '015_water_ne_m00.txt','015_water_ne_m00',
     &   .true.,  .false., .false., .false.)
      call test_mutate_fixed ('016_water_ye_p10.key',
     &   '016_water_ye_p10.txt','016_water_ye_p10',
     &   .false., .true.,  .false., .true.)
      call test_mutate_fixed ('017_water_ne_p10.key',
     &   '017_water_ne_p10.txt','017_water_ne_p10',
     &   .false., .true.,  .false., .false.)
      call test_mutate_fixed ('018_water_ye_p05.key',
     &   '018_water_ye_p05.txt','018_water_ye_p05',
     &   .false., .true.,  .false., .true.)
      call test_mutate_fixed ('019_water_ne_p05.key',
     &   '019_water_ne_p05.txt','019_water_ne_p05',
     &   .false., .true.,  .false., .false.)
      call test_mutate_fixed ('020_water_ye_p00.key',
     &   '020_water_ye_p00.txt','020_water_ye_p00',
     &   .false., .true.,  .false., .true.)
      call test_mutate_fixed ('021_water_ne_p00.key',
     &   '021_water_ne_p00.txt','021_water_ne_p00',
     &   .false., .true.,  .false., .false.)
      call test_mutate_fixed ('022_water_ye_m10p05.key',
     &   '022_water_ye_m10p05.txt','022_water_ye_m10p05',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('023_water_ne_m10p05.key',
     &   '023_water_ne_m10p05.txt','023_water_ne_m10p05',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('024_water_ye_m05p10.key',
     &   '024_water_ye_m05p10.txt','024_water_ye_m05p10',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('025_water_ne_m05p10.key',
     &   '025_water_ne_m05p10.txt','025_water_ne_m05p10',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('026_water_ye_m05p00.key',
     &   '026_water_ye_m05p00.txt','026_water_ye_m05p00',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('027_water_ne_m05p00.key',
     &   '027_water_ne_m05p00.txt','027_water_ne_m05p00',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('028_water_ye_m00p05.key',
     &   '028_water_ye_m00p05.txt','028_water_ye_m00p05',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('029_water_ne_m00p05.key',
     &   '029_water_ne_m00p05.txt','029_water_ne_m00p05',
     &   .true.,  .true.,  .false., .false.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_mutate_ast  --  single topology lambda  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_mutate_ast" runs the eighteen absolute single topology water
c     fixtures 030-040, 182 and 187-192, each carrying "lambda-deriv2"
c     (187-192 carry "lambda-deriv", first derivatives only), and drives
c     them through "test_mutate_calc" with the level 4 lambda derivative
c     checks enabled; cases 030-035 keep only the
c     multipole term at three ele-lambda values with Ewald on and off;
c     cases 036-038 keep only the van der Waals term at three vdw-lambda
c     values, and cases 039-040 leave the multipole and polarization
c     terms active with Ewald on and off; case 182 leaves all three
c     nonbonded terms active and annihilates van der Waals interactions;
c     cases 187-192 leave all three nonbonded terms active at matched
c     vdw-lambda, ele-lambda and pol-lambda values 1.0, 0.5 and 0.0 with
c     Ewald on and off; the no-Ewald cases cannot use a pairwise
c     neighbor list
c
c
      subroutine test_mutate_ast
      implicit none
c
c
      call test_mutate_calc ('water2','030_water_ast_ye_m10.key',
     &   '030_water_ast_ye_m10.txt','030_water_ast_ye_m10',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','031_water_ast_ne_m10.key',
     &   '031_water_ast_ne_m10.txt','031_water_ast_ne_m10',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','032_water_ast_ye_m05.key',
     &   '032_water_ast_ye_m05.txt','032_water_ast_ye_m05',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','033_water_ast_ne_m05.key',
     &   '033_water_ast_ne_m05.txt','033_water_ast_ne_m05',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','034_water_ast_ye_m00.key',
     &   '034_water_ast_ye_m00.txt','034_water_ast_ye_m00',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','035_water_ast_ne_m00.key',
     &   '035_water_ast_ne_m00.txt','035_water_ast_ne_m00',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','036_water_ast_v10.key',
     &   '036_water_ast_v10.txt','036_water_ast_v10',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','037_water_ast_v05.key',
     &   '037_water_ast_v05.txt','037_water_ast_v05',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','038_water_ast_v00.key',
     &   '038_water_ast_v00.txt','038_water_ast_v00',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','039_water_ast_ye_mp05.key',
     &   '039_water_ast_ye_mp05.txt','039_water_ast_ye_mp05',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','040_water_ast_ne_mp05.key',
     &   '040_water_ast_ne_mp05.txt','040_water_ast_ne_mp05',
     &   .true.,  .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2',
     &   '182_water_ast_v05_annihilate.key',
     &   '182_water_ast_v05_annihilate.txt',
     &   '182_water_ast_v05_annihilate',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','187_water_ast_ye_l10.key',
     &   '187_water_ast_ye_l10.txt','187_water_ast_ye_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','188_water_ast_ne_l10.key',
     &   '188_water_ast_ne_l10.txt','188_water_ast_ne_l10',
     &   .true.,  .true.,  .true.,  .false., .true.)
      call test_mutate_calc ('water2','189_water_ast_ye_l05.key',
     &   '189_water_ast_ye_l05.txt','189_water_ast_ye_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','190_water_ast_ne_l05.key',
     &   '190_water_ast_ne_l05.txt','190_water_ast_ne_l05',
     &   .true.,  .true.,  .true.,  .false., .true.)
      call test_mutate_calc ('water2','191_water_ast_ye_l00.key',
     &   '191_water_ast_ye_l00.txt','191_water_ast_ye_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','192_water_ast_ne_l00.key',
     &   '192_water_ast_ne_l00.txt','192_water_ast_ne_l00',
     &   .true.,  .true.,  .true.,  .false., .true.)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_mutate_adt  --  abs dual topo lambda  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_mutate_adt" runs the eight absolute polarization dual
c     topology water fixtures 047-052 and 056-057, each carrying the
c     "lambda-deriv2" keyword and the POL-DUALTOPO keyword, and drives
c     them through "test_mutate_calc" with the level 4 lambda derivative
c     checks enabled; cases 047-052 keep only the polarization term at
c     three pol-lambda values, and cases 056-057 also keep the multipole
c     term, on its single topology path, with Ewald on and off; the
c     no-Ewald cases cannot use a pairwise neighbor list
c
c
      subroutine test_mutate_adt
      implicit none
c
c
      call test_mutate_calc ('water2','047_water_adt_ye_p10.key',
     &   '047_water_adt_ye_p10.txt','047_water_adt_ye_p10',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','048_water_adt_ne_p10.key',
     &   '048_water_adt_ne_p10.txt','048_water_adt_ne_p10',
     &   .false., .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','049_water_adt_ye_p05.key',
     &   '049_water_adt_ye_p05.txt','049_water_adt_ye_p05',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','050_water_adt_ne_p05.key',
     &   '050_water_adt_ne_p05.txt','050_water_adt_ne_p05',
     &   .false., .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','051_water_adt_ye_p00.key',
     &   '051_water_adt_ye_p00.txt','051_water_adt_ye_p00',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','052_water_adt_ne_p00.key',
     &   '052_water_adt_ne_p00.txt','052_water_adt_ne_p00',
     &   .false., .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','056_water_adt_ye_mp05.key',
     &   '056_water_adt_ye_mp05.txt','056_water_adt_ye_mp05',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','057_water_adt_ne_mp05.key',
     &   '057_water_adt_ne_mp05.txt','057_water_adt_ne_mp05',
     &   .true.,  .true.,  .false., .false., .true.)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_mutate_qnt  --  quintic lambda map  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_mutate_qnt" runs the six water fixtures 075-080 that map the
c     main lambda to the electrostatics, polarization and van der Waals
c     sub-lambdas with the quintic "qnt" scheme; each carries the "ost"
c     keyword and drives "test_mutate_calc" with the level 4 lambda
c     derivative checks enabled; cases 075-077 use single topology and
c     078-080 use polarization dual topology with single topology
c     multipoles and van der Waals switched off, each at main lambda
c     values 1.0, 0.5 and 0.0; all fixtures use Ewald and support a
c     pairwise neighbor list
c
c
      subroutine test_mutate_qnt
      implicit none
c
c
      call test_mutate_calc ('water2','075_water_qnt_ast_l10.key',
     &   '075_water_qnt_ast_l10.txt','075_water_qnt_ast_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','076_water_qnt_ast_l05.key',
     &   '076_water_qnt_ast_l05.txt','076_water_qnt_ast_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','077_water_qnt_ast_l00.key',
     &   '077_water_qnt_ast_l00.txt','077_water_qnt_ast_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','078_water_qnt_adt_l10.key',
     &   '078_water_qnt_adt_l10.txt','078_water_qnt_adt_l10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','079_water_qnt_adt_l05.key',
     &   '079_water_qnt_adt_l05.txt','079_water_qnt_adt_l05',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','080_water_qnt_adt_l00.key',
     &   '080_water_qnt_adt_l00.txt','080_water_qnt_adt_l00',
     &   .true.,  .true.,  .false., .true.,  .true.)
      return
      end
c
c
c     ######################################################
c     ##                                                  ##
c     ##  subroutine test_mutate_exp  --  exp lambda map  ##
c     ##                                                  ##
c     ######################################################
c
c
c     "test_mutate_exp" runs the six water fixtures 084-089 that map the
c     main lambda to the electrostatics, polarization and van der Waals
c     sub-lambdas with the exponential "exp" scheme; each carries the
c     "ost" keyword and drives "test_mutate_calc" with the level 4
c     lambda derivative checks enabled; cases 084-086 use single
c     topology and 087-089 use polarization dual topology with single
c     topology multipoles and van der Waals switched off, each at main
c     lambda values 1.0, 0.5 and 0.0; all fixtures use Ewald and
c     support a pairwise neighbor list
c
c
      subroutine test_mutate_exp
      implicit none
c
c
      call test_mutate_calc ('water2','084_water_exp_ast_l10.key',
     &   '084_water_exp_ast_l10.txt','084_water_exp_ast_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','085_water_exp_ast_l05.key',
     &   '085_water_exp_ast_l05.txt','085_water_exp_ast_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','086_water_exp_ast_l00.key',
     &   '086_water_exp_ast_l00.txt','086_water_exp_ast_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','087_water_exp_adt_l10.key',
     &   '087_water_exp_adt_l10.txt','087_water_exp_adt_l10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','088_water_exp_adt_l05.key',
     &   '088_water_exp_adt_l05.txt','088_water_exp_adt_l05',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','089_water_exp_adt_l00.key',
     &   '089_water_exp_adt_l00.txt','089_water_exp_adt_l00',
     &   .true.,  .true.,  .false., .true.,  .true.)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_mutate_inv  --  inverse lambda map  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_mutate_inv" runs the six water fixtures 093-098 that map the
c     main lambda to the electrostatics, polarization and van der Waals
c     sub-lambdas with the inverse "inv" scheme; each carries the "ost"
c     keyword and drives "test_mutate_calc" with the level 4 lambda
c     derivative checks enabled; cases 093-095 use single topology and
c     096-098 use polarization dual topology with single topology
c     multipoles and van der Waals switched off, each at main lambda
c     values 1.0, 0.5 and 0.0; all fixtures use Ewald and support a
c     pairwise neighbor list
c
c
      subroutine test_mutate_inv
      implicit none
c
c
      call test_mutate_calc ('water2','093_water_inv_ast_l10.key',
     &   '093_water_inv_ast_l10.txt','093_water_inv_ast_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','094_water_inv_ast_l05.key',
     &   '094_water_inv_ast_l05.txt','094_water_inv_ast_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','095_water_inv_ast_l00.key',
     &   '095_water_inv_ast_l00.txt','095_water_inv_ast_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','096_water_inv_adt_l10.key',
     &   '096_water_inv_adt_l10.txt','096_water_inv_adt_l10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','097_water_inv_adt_l05.key',
     &   '097_water_inv_adt_l05.txt','097_water_inv_adt_l05',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','098_water_inv_adt_l00.key',
     &   '098_water_inv_adt_l00.txt','098_water_inv_adt_l00',
     &   .true.,  .true.,  .false., .true.,  .true.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_mutate_exf  --  applied external field  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_mutate_exf" runs the thirteen water fixtures 102-104,
c     111-113, 117, 184-186 and 193-195 that apply an external electric
c     field to a mutated system, each carrying "lambda-deriv2" (193-195
c     carry "lambda-deriv") so the level 4 checks are enabled;
c     cases 102-104 keep only the multipole term under single topology
c     at ele-lambda values 1.0, 0.5 and 0.0; cases 111-113 keep only the
c     polarization term under absolute dual topology at pol-lambda
c     values 1.0, 0.5 and 0.0, since the
c     induced dipoles respond to the applied field through the direct
c     field; case 117 leaves the multipole and polarization terms active
c     together under absolute polarization dual topology with single
c     topology multipoles; cases 184-186 leave all three nonbonded terms
c     active under the same topology at matched ele-lambda and
c     pol-lambda values 1.0, 0.5 and 0.0, so the unscaled van der Waals
c     term is carried alongside the field-driven multipole and
c     polarization terms; cases 193-195 repeat that all-term case under
c     single topology at matched vdw-lambda, ele-lambda and pol-lambda
c     values 1.0, 0.5 and 0.0; the dual topology cases use an exponent
c     of two, since with the default exponent of one the dual topology
c     weighting of the external field term cannot be distinguished from
c     a linear scaling; all fixtures use Ewald and support a pairwise
c     neighbor list
c
c
      subroutine test_mutate_exf
      implicit none
c
c
      call test_mutate_calc ('water2','102_water_exf_ast_m10.key',
     &   '102_water_exf_ast_m10.txt','102_water_exf_ast_m10',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','103_water_exf_ast_m05.key',
     &   '103_water_exf_ast_m05.txt','103_water_exf_ast_m05',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','104_water_exf_ast_m00.key',
     &   '104_water_exf_ast_m00.txt','104_water_exf_ast_m00',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','111_water_exf_adt_p10.key',
     &   '111_water_exf_adt_p10.txt','111_water_exf_adt_p10',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','112_water_exf_adt_p05.key',
     &   '112_water_exf_adt_p05.txt','112_water_exf_adt_p05',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','113_water_exf_adt_p00.key',
     &   '113_water_exf_adt_p00.txt','113_water_exf_adt_p00',
     &   .false., .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','117_water_exf_adt_mp05.key',
     &   '117_water_exf_adt_mp05.txt','117_water_exf_adt_mp05',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','184_water_exf_adt_l10.key',
     &   '184_water_exf_adt_l10.txt','184_water_exf_adt_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','185_water_exf_adt_l05.key',
     &   '185_water_exf_adt_l05.txt','185_water_exf_adt_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','186_water_exf_adt_l00.key',
     &   '186_water_exf_adt_l00.txt','186_water_exf_adt_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','193_water_exf_ast_l10.key',
     &   '193_water_exf_ast_l10.txt','193_water_exf_ast_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','194_water_exf_ast_l05.key',
     &   '194_water_exf_ast_l05.txt','194_water_exf_ast_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','195_water_exf_ast_l00.key',
     &   '195_water_exf_ast_l00.txt','195_water_exf_ast_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_mutate_emplar  --  mpole plus polar dt  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_mutate_emplar" runs the six water fixtures 119-124 that
c     keep the multipole and polarization terms active together at
c     matched ele-lambda and pol-lambda values, each carrying the
c     "lambda-deriv2" keyword and the POL-DUALTOPO keyword, and drives
c     them through "test_mutate_calc" with the level 4 lambda derivative
c     checks enabled; the multipoles use single topology, and each case
c     sits at ele-/pol-lambda values 1.0, 0.5 and 0.0 with Ewald on and
c     off; the polarization term uses an interpolation exponent of
c     three; the no-Ewald cases cannot use a pairwise neighbor list
c
c
      subroutine test_mutate_emplar
      implicit none
c
c
      call test_mutate_calc ('water2','119_water_adt_ye_l10.key',
     &   '119_water_adt_ye_l10.txt','119_water_adt_ye_l10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','120_water_adt_ne_l10.key',
     &   '120_water_adt_ne_l10.txt','120_water_adt_ne_l10',
     &   .true.,  .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','121_water_adt_ye_l05.key',
     &   '121_water_adt_ye_l05.txt','121_water_adt_ye_l05',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','122_water_adt_ne_l05.key',
     &   '122_water_adt_ne_l05.txt','122_water_adt_ne_l05',
     &   .true.,  .true.,  .false., .false., .true.)
      call test_mutate_calc ('water2','123_water_adt_ye_l00.key',
     &   '123_water_adt_ye_l00.txt','123_water_adt_ye_l00',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','124_water_adt_ne_l00.key',
     &   '124_water_adt_ne_l00.txt','124_water_adt_ne_l00',
     &   .true.,  .true.,  .false., .false., .true.)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_qntrng  --  narrowed quintic range  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_qntrng" runs the two water fixtures 131-132 that
c     narrow the quintic sub-lambda windows to 0.1-0.9 with the
c     "ele-lmda-range", "pol-lmda-range" and "vdw-lmda-range" keywords,
c     then sample the main lambda at the 1.0 and 0.0 endpoints which
c     now sit outside every window; there the quintic taper is on its
c     flat plateau, so each sub-lambda is pinned fully coupled at 1.0
c     and fully decoupled at 0.0 while every first and second lambda
c     derivative vanishes, which the level 4 checks confirm; both use
c     polarization dual topology with single topology multipoles and
c     van der Waals switched off, and both use Ewald and support a
c     pairwise neighbor list
c
c
      subroutine test_mutate_qntrng
      implicit none
c
c
      call test_mutate_calc ('water2','131_water_qnt_adt_l10.key',
     &   '131_water_qnt_adt_l10.txt','131_water_qnt_adt_l10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','132_water_qnt_adt_l00.key',
     &   '132_water_qnt_adt_l00.txt','132_water_qnt_adt_l00',
     &   .true.,  .true.,  .false., .true.,  .true.)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_mutate_rels  --  staged rel dual topo  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_mutate_rels" runs the twelve water fixtures 135-141, 168-169
c     and 203-205 that drive the staged relative free energy schedule,
c     each of them
c     naming the one leg it walks, so ligand 2 is discharged on the LIG2
c     leg, van der Waals morphs between the ligands on the VDWM leg, and
c     ligand 1 is charged on the LIG1 leg; the cases sit one per regime,
c     at main lambda values 1.0 and 0.0 where a leg is flat and only the
c     coupled endpoint is built, 0.85 and 0.15 inside the two mixing
c     legs, and 0.5 in the middle of the morph window where both ligands
c     are annihilated and every electrostatic lambda derivative vanishes
c
c     the leg boundaries at 0.7 and 0.3 are each run twice, once from
c     either side, since the quintic taper is flat at both ends of every
c     window and leaves the boundary state derivative free; 139 and 168
c     are the same state reached from the LIG2 and VDWM legs, and 137 and
c     169 the same state reached from the VDWM and LIG1 legs, so each
c     pair shares one set of reference values
c
c     the first nine carry the "lambda-deriv2" keyword, so polarization
c     takes the dual topology path between the annihilated endpoints;
c     203 and 204 repeat 136 and 177 with the "lambda-deriv" keyword, so
c     polarization takes the single topology path and only the first
c     lambda derivative is checked; 205 drops the "lambda" keyword from
c     135, so REL-STAGE must default the main lambda to one and match
c     the 135 reference; all run the level 4 checks
c
c
      subroutine test_mutate_rels
      implicit none
c
c
      call test_mutate_calc ('water2','135_water_rels_ye_l100.key',
     &   '135_water_rels_ye_l100.txt','135_water_rels_ye_l100',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','136_water_rels_ye_l085.key',
     &   '136_water_rels_ye_l085.txt','136_water_rels_ye_l085',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','137_water_rels_ye_l070.key',
     &   '137_water_rels_ye_l070.txt','137_water_rels_ye_l070',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','169_water_rels_ye_lig1_l070.key',
     &   '169_water_rels_ye_lig1_l070.txt',
     &   '169_water_rels_ye_lig1_l070',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','138_water_rels_ye_l050.key',
     &   '138_water_rels_ye_l050.txt','138_water_rels_ye_l050',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','139_water_rels_ye_l030.key',
     &   '139_water_rels_ye_l030.txt','139_water_rels_ye_l030',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','168_water_rels_ye_vdwm_l030.key',
     &   '168_water_rels_ye_vdwm_l030.txt',
     &   '168_water_rels_ye_vdwm_l030',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','140_water_rels_ye_l015.key',
     &   '140_water_rels_ye_l015.txt','140_water_rels_ye_l015',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','141_water_rels_ye_l000.key',
     &   '141_water_rels_ye_l000.txt','141_water_rels_ye_l000',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','203_water_rels_st_l085.key',
     &   '203_water_rels_st_l085.txt','203_water_rels_st_l085',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2',
     &   '204_water_rels_st_lig2_exp_l030.key',
     &   '204_water_rels_st_lig2_exp_l030.txt',
     &   '204_water_rels_st_lig2_exp_l030',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','205_water_rels_nolmda.key',
     &   '135_water_rels_ye_l100.txt','205_water_rels_nolmda',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_mutate_lmdafix  --  main lambda fixtures  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_mutate_lmdafix" runs the eight water fixtures 146-153 that
c     drive the sublambdas from the "lambda" keyword at a fixed lambda
c     value; no fixture carries a lambda derivative keyword, so the level
c     4 checks stay off and only the energy, gradient, virial and named
c     components are verified; cases 146 and 147, 148 and 149, and 152
c     and 153 are pairs that must agree, respectively matching the
c     linear map against explicit sublambdas, the fully coupled endpoint
c     against a fixture with no lambda at all, and two named maps
c     against the same state set by explicit sublambdas; case 150 maps
c     only electrostatics while 151 maps only van der Waals, each
c     leaving the unnamed terms fully coupled; all fixtures use Ewald
c     and support a neighbor list
c
c
      subroutine test_mutate_lmdafix
      implicit none
c
c
      call test_mutate_calc ('water2','146_water_lmda_ast_l05.key',
     &   '146_water_lmda_ast_l05.txt','146_water_lmda_ast_l05',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','147_water_lmda_ast_e05.key',
     &   '147_water_lmda_ast_e05.txt','147_water_lmda_ast_e05',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','148_water_lmda_ast_l10.key',
     &   '148_water_lmda_ast_l10.txt','148_water_lmda_ast_l10',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','149_water_lmda_ast_none.key',
     &   '149_water_lmda_ast_none.txt','149_water_lmda_ast_none',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','150_water_lmda_qnt_l05.key',
     &   '150_water_lmda_qnt_l05.txt','150_water_lmda_qnt_l05',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','151_water_lmda_vexp_l05.key',
     &   '151_water_lmda_vexp_l05.txt','151_water_lmda_vexp_l05',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','152_water_lmda_mp05.key',
     &   '152_water_lmda_mp05.txt','152_water_lmda_mp05',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','153_water_lmda_mp05_expl.key',
     &   '153_water_lmda_mp05_expl.txt','153_water_lmda_mp05_expl',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_lmdadrv  --  driven sublambda sets  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_lmdadrv" runs the fourteen water fixtures 154-167
c     that name each combination of the electrostatic, polarization and
c     van der Waals maps at a main lambda of 0.6, so every fixture
c     drives one subset of the sublambdas on the identity map and holds
c     the rest fully coupled at one; all three nonbonded terms stay
c     active throughout, so an undriven term is present in the energy
c     rather than switched off; cases 154-160 carry no "lambda-deriv2"
c     keyword while 161-167 repeat them with it, and the two halves
c     share their energy, gradient and virial values, which holds the
c     plain and the lambda derivative energy routines to the same
c     result; the level 4 checks run only for the derivative half
c
c
      subroutine test_mutate_lmdadrv
      implicit none
c
c
      call test_mutate_calc ('water2','154_water_lmda_e_l06.key',
     &   '154_water_lmda_e_l06.txt','154_water_lmda_e_l06',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','155_water_lmda_p_l06.key',
     &   '155_water_lmda_p_l06.txt','155_water_lmda_p_l06',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','156_water_lmda_v_l06.key',
     &   '156_water_lmda_v_l06.txt','156_water_lmda_v_l06',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','157_water_lmda_ep_l06.key',
     &   '157_water_lmda_ep_l06.txt','157_water_lmda_ep_l06',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','158_water_lmda_ev_l06.key',
     &   '158_water_lmda_ev_l06.txt','158_water_lmda_ev_l06',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','159_water_lmda_pv_l06.key',
     &   '159_water_lmda_pv_l06.txt','159_water_lmda_pv_l06',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','160_water_lmda_epv_l06.key',
     &   '160_water_lmda_epv_l06.txt','160_water_lmda_epv_l06',
     &   .true.,  .true.,  .true.,  .true.,  .false.)
      call test_mutate_calc ('water2','161_water_dlmda_e_l06.key',
     &   '161_water_dlmda_e_l06.txt','161_water_dlmda_e_l06',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','162_water_dlmda_p_l06.key',
     &   '162_water_dlmda_p_l06.txt','162_water_dlmda_p_l06',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','163_water_dlmda_v_l06.key',
     &   '163_water_dlmda_v_l06.txt','163_water_dlmda_v_l06',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','164_water_dlmda_ep_l06.key',
     &   '164_water_dlmda_ep_l06.txt','164_water_dlmda_ep_l06',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','165_water_dlmda_ev_l06.key',
     &   '165_water_dlmda_ev_l06.txt','165_water_dlmda_ev_l06',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','166_water_dlmda_pv_l06.key',
     &   '166_water_dlmda_pv_l06.txt','166_water_dlmda_pv_l06',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','167_water_dlmda_epv_l06.key',
     &   '167_water_dlmda_epv_l06.txt','167_water_dlmda_epv_l06',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_lmdapin  --  pinned sublambda legs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_lmdapin" runs the three water fixtures 170, 171 and
c     173 that hold one term at a fixed coupling state while the main
c     lambda morphs another, the way a leg of a two simulation free
c     energy split is set up; the "epin" fixtures pin electrostatics and
c     polarization decoupled by their own keywords and drive van der
c     Waals from the "lambda" keyword, while the "vpin" fixtures pin
c     van der Waals decoupled and drive the other two
c
c     the fixtures walk one topology and one map form each, 170 and
c     171 absolute single topology on the exponential map at a main
c     lambda of 0.5, and 173 absolute polarization dual topology on the
c     quintic map at 0.6; 171 carries POL-DUALTOPO even though
c     the rest of that pair is single topology, since "epolar4"
c     supplies the polarization lambda derivative only through dual
c     topology
c
c     a pinned sublambda is held at its own keyword value and leaves the
c     chain rule, so every fixture names the "lambda-deriv2" keyword and
c     runs the level 4 checks, where the pinned terms report exact zeros
c     and the driven ones carry the whole derivative; all three
c     nonbonded terms stay active, so a pinned term is still present in
c     the energy rather than switched off, and all three fixtures use
c     Ewald and support a neighbor list
c
c
      subroutine test_mutate_lmdapin
      implicit none
c
c
      call test_mutate_calc ('water2','170_water_lmda_ast_epin_l05.key',
     &   '170_water_lmda_ast_epin_l05.txt',
     &   '170_water_lmda_ast_epin_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','171_water_lmda_ast_vpin_l05.key',
     &   '171_water_lmda_ast_vpin_l05.txt',
     &   '171_water_lmda_ast_vpin_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','173_water_lmda_adt_vpin_l06.key',
     &   '173_water_lmda_adt_vpin_l06.txt',
     &   '173_water_lmda_adt_vpin_l06',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_mutate_relsmap  --  staged leg map forms  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_mutate_relsmap" runs the six water fixtures 176-181 that
c     drive a staged relative leg with a map form and a dual topology
c     exponent that the staged schedule once refused, every leg now
c     walking the map named for the one term it drives and honoring the
c     exponent of that term as any other relative leg does
c
c     176 to 178 name a map other than the quintic, taking the VDWM leg
c     on the exponential, the LIG2 leg on the exponential through the
c     complement that discharges ligand 2, and the LIG1 leg on the
c     inverse power; the exponential and inverse power maps carry no
c     window, so these three walk the whole main lambda range rather
c     than a slice of it and are each a run of their own
c
c     179 to 181 exercise the dual topology exponents, 179 and 180 on
c     the quintic windows of 138 and 136 and 181 on the inverse power,
c     so each is the state of an existing fixture reached with a
c     different endpoint weight; the exponent is placed on the term the
c     leg drives, since a ligand leg pins the van der Waals endpoint and
c     the VDWM leg holds both electrostatic endpoints at the same
c     decoupled state, leaving the other exponents with nothing to scale
c
c     every fixture names the "lambda-deriv2" keyword and runs the level
c     4 checks, all three nonbonded terms stay active, and all six use
c     Ewald and support a neighbor list
c
c
      subroutine test_mutate_relsmap
      implicit none
c
c
      call test_mutate_calc ('water2',
     &   '176_water_rels_ye_vdwm_exp_l050.key',
     &   '176_water_rels_ye_vdwm_exp_l050.txt',
     &   '176_water_rels_ye_vdwm_exp_l050',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2',
     &   '177_water_rels_ye_lig2_exp_l030.key',
     &   '177_water_rels_ye_lig2_exp_l030.txt',
     &   '177_water_rels_ye_lig2_exp_l030',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2',
     &   '178_water_rels_ye_lig1_inv_l070.key',
     &   '178_water_rels_ye_lig1_inv_l070.txt',
     &   '178_water_rels_ye_lig1_inv_l070',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2',
     &   '179_water_rels_ye_vdwm_vx3_l050.key',
     &   '179_water_rels_ye_vdwm_vx3_l050.txt',
     &   '179_water_rels_ye_vdwm_vx3_l050',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2',
     &   '180_water_rels_ye_lig1_ex3_l085.key',
     &   '180_water_rels_ye_lig1_ex3_l085.txt',
     &   '180_water_rels_ye_lig1_ex3_l085',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2',
     &   '181_water_rels_ye_lig2_ix2_l015.key',
     &   '181_water_rels_ye_lig2_ix2_l015.txt',
     &   '181_water_rels_ye_lig2_ix2_l015',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_apm  --  asymmetric power map legs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_apm" runs the four water fixtures 196-199 that map
c     the main lambda onto a sublambda with the asymmetric power "apm"
c     scheme, whose shape is compiled into "mutate_dlmda" rather than
c     read from the keyfile, so a fixture names the map and nothing
c     else; each carries the "lambda-deriv" keyword and drives
c     "test_mutate_calc" with the level 4 lambda derivative checks
c     enabled, which is where the map earns its keep, since its chain
c     rule factor is neither the unit slope of the linear map nor the
c     vanishing endpoint slope of the quintic taper
c
c     the fixtures pin one set of terms and drive the rest the way the
c     170-175 pair does; 196 copies 170, pinning electrostatics and
c     polarization decoupled by their own keywords and driving van der
c     Waals across the map at a main lambda of 0.5, while 197-199 copy
c     171, pinning van der Waals decoupled and driving electrostatics
c     and polarization at main lambda values of 0.0, 0.5 and 1.0; the
c     three walk the whole interval, so they cover both endpoints of
c     the map, where the slope ratio it normalizes to is defined, as
c     well as the interior; 197-199 carry POL-DUALTOPO for the same
c     reason 171 does, "epolar4" supplying the polarization lambda
c     derivative only through dual topology
c
c     all four use Ewald and support a pairwise neighbor list, and all
c     three nonbonded terms stay active, so a pinned term is still
c     present in the energy rather than switched off
c
c
      subroutine test_mutate_apm
      implicit none
c
c
      call test_mutate_calc ('water2','196_water_apm_ast_vpin_l05.key',
     &   '196_water_apm_ast_vpin_l05.txt',
     &   '196_water_apm_ast_vpin_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','197_water_apm_ast_epin_l00.key',
     &   '197_water_apm_ast_epin_l00.txt',
     &   '197_water_apm_ast_epin_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','198_water_apm_ast_epin_l05.key',
     &   '198_water_apm_ast_epin_l05.txt',
     &   '198_water_apm_ast_epin_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','199_water_apm_ast_epin_l10.key',
     &   '199_water_apm_ast_epin_l10.txt',
     &   '199_water_apm_ast_epin_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_mutate_vsoft  --  softcore vdw settings  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_mutate_vsoft" runs the three water fixtures 200-202 that
c     scan an exponential van der Waals lambda map using a small
c     softcore alpha and a cubic softcore exponent; electrostatics and
c     polarization are pinned decoupled while van der Waals is checked
c     at both endpoints and the midpoint; each fixture carries the
c     "lambda-deriv" keyword and runs the level 4 lambda derivative
c     checks, and all three use Ewald and support a neighbor list
c
c
      subroutine test_mutate_vsoft
      implicit none
c
c
      call test_mutate_calc ('water2','200_water_vsoft_l10.key',
     &   '200_water_vsoft_l10.txt','200_water_vsoft_l10',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','201_water_vsoft_l05.key',
     &   '201_water_vsoft_l05.txt','201_water_vsoft_l05',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','202_water_vsoft_l00.key',
     &   '202_water_vsoft_l00.txt','202_water_vsoft_l00',
     &   .true.,  .true.,  .true.,  .true.,  .true.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_mutate_chiral  --  chiral frame refresh  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_mutate_chiral" mutates the first residues of trp-cage, whose
c     alpha carbons carry chiral multipole frames, with the multipole
c     term alone at a fixed electrostatic lambda; the shipped alpha
c     carbon multipoles have no y components, so a y dipole is added at
c     every chiral site to make the inversion visible; after the
c     coordinates are mirrored, "chkpole" inverts the chiral multipoles,
c     and a later reinstall of the scaled parameters by "altelec" must
c     keep that inversion, so the mirrored energy equals the original
c
c
      subroutine test_mutate_chiral
      use atoms
      use dlmda
      use mpole
      implicit none
      integer i
      real*8 e0,e1,e2
      real*8 energy
      logical skiptest
c
c
      if (skiptest('test_mutate_chiral','mutate'))  return
      call pushdir ('file/mutate')
      call loadfix ('../angle/trpcage','206_trpcage_chiral_m05.key')
      do i = 1, n
         if (polaxe(i).eq.'Z-then-X' .and. yaxis(i).ne.0) then
            poleorig(3,i) = 0.1d0
         end if
      end do
      call altelec
      e0 = energy ()
c
c     mirror the structure, which inverts every chiral frame
c
      do i = 1, n
         x(i) = -x(i)
      end do
      e1 = energy ()
      call assert_real (e1,e0,1.0d-8,'test_mutate_chiral mirrored')
c
c     reinstall the scaled parameters from their original values
c
      call altelec
      e2 = energy ()
      call assert_real (e2,e0,1.0d-8,'test_mutate_chiral reinstall')
      call popdir
      call final
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_mutate_gate  --  lambda derivative gate  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_mutate_gate" runs the staged ligand 1 charging fixture 136,
c     whose multipole second lambda derivative is nonzero, with and
c     without the second, force and virial lambda derivatives; turning
c     them off must leave the energy, gradient and first lambda
c     derivatives unchanged and clear the multipole second, force and
c     virial lambda derivatives
c
c
      subroutine test_mutate_gate
      use atoms
      use dlmda
      implicit none
      integer i,j
      real*8 e,e1,dedl1,demdl1
      real*8 dfm,dvm
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: g1(:,:)
      logical skiptest
c
c
      if (skiptest('test_mutate_gate','mutate'))  return
      call pushdir ('file/mutate')
      call loadfix ('water2','136_water_rels_ye_l085.key')
      allocate (derivs(3,n))
      allocate (g1(3,n))
c
c     full lambda derivatives as the LAMBDA-DERIV keyword requests
c
      call gradient (e,derivs)
      e1 = e
      dedl1 = dedl
      demdl1 = demdl
      do i = 1, n
         do j = 1, 3
            g1(j,i) = derivs(j,i)
         end do
      end do
      call assert_logical (d2emdl2.ne.0.0d0,.true.,
     &                     'test_mutate_gate d2EM/dL2 built')
c
c     only the first lambda derivative, as TI, META and ABF request
c
      use_d2lmda = .false.
      call gradient (e,derivs)
      call assert_real (e,e1,1.0d-10,'test_mutate_gate energy')
      call assert_grad (derivs,g1,n,1.0d-10,'test_mutate_gate gradient')
      call assert_real (dedl,dedl1,1.0d-10,'test_mutate_gate dE/dL')
      call assert_real (demdl,demdl1,1.0d-10,'test_mutate_gate dEM/dL')
      call assert_real (d2emdl2,0.0d0,0.0d0,
     &                  'test_mutate_gate d2EM/dL2 cleared')
      dfm = 0.0d0
      do i = 1, n
         do j = 1, 3
            dfm = dfm + abs(dfmdl(j,i))
         end do
      end do
      dvm = 0.0d0
      do i = 1, 3
         do j = 1, 3
            dvm = dvm + abs(demvirdl(j,i))
         end do
      end do
      call assert_real (dfm,0.0d0,0.0d0,
     &                  'test_mutate_gate dFM/dL cleared')
      call assert_real (dvm,0.0d0,0.0d0,
     &                  'test_mutate_gate dVM/dL cleared')
      use_d2lmda = .true.
      deallocate (derivs)
      deallocate (g1)
      call popdir
      call final
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_mutate_fixed  --  one mutation case  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_mutate_fixed" runs a single mutation fixture with and
c     without the neighbor-list keyword, checking the level 0/1/3
c     energy, gradient, virial and named component regressions; it is
c     a thin wrapper around "test_mutate_calc" with the level 4 lambda
c     derivative checks disabled, and backs the "test_mutate_mv" and
c     "test_mutate_mp" case lists
c
c
      subroutine test_mutate_fixed
     &   (key,ref,cname,checkm,checkp,checkv,canlist)
      implicit none
      logical checkm,checkp,checkv,canlist
      character*(*) key,ref,cname
c
c
      call test_mutate_calc ('water',key,ref,cname,checkm,checkp,
     &                       checkv,canlist,.false.)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_mutate_calc  --  one mutation case  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_mutate_calc" runs a single mutation fixture with and
c     without the neighbor-list keyword; for each neighbor-list
c     variant the force field is built once, then the checks are
c     repeated twice before teardown, giving four passes per fixture;
c     the "checkm", "checkp" and "checkv" flags select which level 3
c     energy components are verified (Atomic Multipoles, Polarization
c     and Van der Waals); the "canlist" flag records whether the
c     neighbor-list variant is compatible with this fixture (false for
c     the no-Ewald cases); when "dolmda" is set, the level 4 checks
c     also verify the analytical lambda derivatives, second lambda
c     derivatives, per-atom lambda gradient and dV/dL tensor against
c     the reference read by "load_lmdaref", reusing the values already
c     produced by the level 1 gradient call
c
c
      subroutine test_mutate_calc
     &   (base,key,ref,cname,checkm,checkp,checkv,canlist,dolmda)
      use action
      use atoms
      use dlmda
      use energi
      use virial
      implicit none
      integer nat,natl,irun,ilist,nlist
      real*8 energy,e,ref_e,ref_ei
      real*8 eps_e,eps_g,eps_v,eps_l,refv(3,3)
      real*8 ref_dedl(4),ref_d2edl2(4),ref_dvdl(3,3)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: refg(:,:)
      real*8, allocatable :: ref_lg(:,:)
      logical skiptest,checkm,checkp,checkv,canlist,uselist,dolmda
      character*(*) base,key,ref,cname
      character*240 rpath,pre
      character*8 rtag
c
c
      if (skiptest('test_'//trim(cname),'mutate'))  return
c
c     run each fixture without the neighbor-list keyword, and with it
c     when the fixture supports a pairwise neighbor list
c
      nlist = 2
      if (.not. canlist)  nlist = 1
      do ilist = 1, nlist
         uselist = (ilist .eq. 2)
c
c     set up the force field once for this neighbor-list variant
c
         call pushdir ('file/mutate')
         if (uselist) then
            call loadfix_keyadd (base,key,'neighbor-list')
            pre = 'test_'//trim(cname)//' list'
         else
            call loadfix (base,key)
            pre = 'test_'//trim(cname)//' nolist'
         end if
         allocate (derivs(3,n))
         allocate (refg(3,n))
         call refpath ('mutate',ref,rpath)
         call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
         eps_e = 1.0d-4
         eps_g = 1.0d-4
         eps_v = 1.0d-3
         if (index(cname,'_vcorr_') .ne. 0)  eps_v = 1.0d-2
c
c     read the reference lambda derivatives for the level 4 checks
c
         if (dolmda) then
            allocate (ref_lg(3,n))
            call load_lmdaref (rpath,n,ref_dedl,ref_d2edl2,ref_lg,
     &                         ref_dvdl,natl)
            eps_l = 1.0d-4
         end if
c
c     repeat the checks twice against the built system
c
         do irun = 1, 2
            if (irun .eq. 1) then
               rtag = ' run1'
            else
               rtag = ' run2'
            end if
c
c     level 0  --  total potential energy
c
            e = energy ()
            call assert_real (esum,ref_e,eps_e,
     &                        trim(pre)//' energy (v0)'//trim(rtag))
c
c     level 1  --  total energy, Cartesian gradient and virial
c
            call gradient (e,derivs)
            call assert_real (esum,ref_e,eps_e,
     &                        trim(pre)//' grad-e (v1)'//trim(rtag))
            call assert_grad (derivs,refg,n,eps_g,
     &                        trim(pre)//' grad (v1)'//trim(rtag))
            call assert_grad (vir,refv,3,eps_v,
     &                        trim(pre)//' virial (v1)'//trim(rtag))
c
c     level 4  --  lambda derivatives from the level 1 gradient call
c
            if (dolmda) then
               call assert_real (dedl,ref_dedl(1),eps_l,
     &                        trim(pre)//' dE/dL (v4)'//trim(rtag))
               call assert_real (devdl,ref_dedl(2),eps_l,
     &                        trim(pre)//' dEV/dL (v4)'//trim(rtag))
               call assert_real (demdl,ref_dedl(3),eps_l,
     &                        trim(pre)//' dEM/dL (v4)'//trim(rtag))
               call assert_real (depdl,ref_dedl(4),eps_l,
     &                        trim(pre)//' dEP/dL (v4)'//trim(rtag))
               call assert_real (d2edl2,ref_d2edl2(1),eps_l,
     &                        trim(pre)//' d2E/dL2 (v4)'//trim(rtag))
               call assert_real (d2evdl2,ref_d2edl2(2),eps_l,
     &                        trim(pre)//' d2EV/dL2 (v4)'//trim(rtag))
               call assert_real (d2emdl2,ref_d2edl2(3),eps_l,
     &                        trim(pre)//' d2EM/dL2 (v4)'//trim(rtag))
               call assert_real (d2epdl2,ref_d2edl2(4),eps_l,
     &                        trim(pre)//' d2EP/dL2 (v4)'//trim(rtag))
               call assert_grad (dfsumdl,ref_lg,n,eps_g,
     &                        trim(pre)//' lgrad (v4)'//trim(rtag))
               call assert_grad (dvirdl,ref_dvdl,3,eps_v,
     &                        trim(pre)//' dV/dL (v4)'//trim(rtag))
            end if
c
c     level 3  --  total and named AMOEBA energy components
c
            call analysis (e)
            call assert_real (esum,ref_e,eps_e,
     &                        trim(pre)//' analysis (v3)'//trim(rtag))
            if (checkm) then
               call check_engcnt (rpath,'Atomic Multipoles',em,nem,
     &                       eps_e,trim(pre)//' mpole (v3)'//trim(rtag))
            end if
            if (checkp) then
               call check_engcnt (rpath,'Polarization',ep,nep,
     &                       eps_e,trim(pre)//' polar (v3)'//trim(rtag))
            end if
            if (checkv) then
               call check_engcnt (rpath,'Van der Waals',ev,nev,
     &                       eps_e,trim(pre)//' vdw (v3)'//trim(rtag))
            end if
         end do
c
c     clean up this neighbor-list variant
c
         deallocate (derivs)
         deallocate (refg)
         if (dolmda)  deallocate (ref_lg)
         call popdir
         call final
      end do
      return
      end

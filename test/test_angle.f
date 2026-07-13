c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_angle  --  angle bending FF tests  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_angle" checks the same angle bending cases exercised by
c     the tinker-gpu angle.cpp tests: a protein AMOEBA angle-only
c     system and a Fourier angle bending system
c
c
      subroutine test_angle
      implicit none
c
c
      call test_angle_trpcage
      call test_angle_fourier
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_angle_trpcage  --  AMOEBA protein angles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_angle_trpcage" checks harmonic AMOEBA angle bending for
c     trpcage against the upstream tinker-gpu angle.1 reference
c
c
      subroutine test_angle_trpcage
      implicit none
c
c
      call test_angle_case ('trpcage','angle.1.txt','amoeba',
     &                      'angle_trpcage',1.0d-4,1.0d-3,3.0d-3)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_angle_fourier  --  Fourier angle bending  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_angle_fourier" checks the Fourier angle bending case
c     against the upstream tinker-gpu angle.2 reference
c
c
      subroutine test_angle_fourier
      implicit none
c
c
      call test_angle_case ('anglef','angle.2.txt','amoeba',
     &                      'angle_fourier',1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_angle_case  --  run one angle test case  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_angle_case" loads one angle fixture, evaluates the angle
c     term through energy, gradient and analysis paths, and compares
c     energy, count, virial and gradient to the stored reference
c
c
      subroutine test_angle_case (base,reffile,tags,tname,
     &                            eps_e,eps_g,eps_v)
      use action
      use angbnd
      use atoms
      use energi
      use virial
      implicit none
      integer nat,refcnt
      real*8 energy,e,ref_e,ref_ei,refeng
      real*8 eps_e,eps_g,eps_v
      real*8, allocatable :: derivs(:,:)
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      logical skiptest
      character*(*) base,reffile,tags,tname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,tags))  return
c
c     load the angle-only fixture and the reference values
c
      call pushdir ('file/angle')
      call loadfix (base,base//'.key')
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('angle',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Angle Bending',refeng,refcnt)
c
c     setup and level 0  --  total angle bending energy
c
      call assert_int (nangle,refcnt,tpre//tname//' count')
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
      call assert_real (ea,refeng,eps_e,tpre//tname//' angle (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,tpre//tname//' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,tpre//tname//' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,tpre//tname//' virial (v1)')
c
c     level 3  --  energy and per-term analysis
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,
     &                  tpre//tname//' analysis (v3)')
      call check_engcnt (rpath,'Angle Bending',ea,nea,eps_e,
     &                   tpre//tname//' angle (v3)')
c
c     clean up
c
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

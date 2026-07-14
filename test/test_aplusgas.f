c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine test_aplusgas  --  AMOEBA+ gas tests  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "test_aplusgas" checks the AMOEBA+ gas-phase case exercised
c     by the tinker-gpu aplusgas.cpp tests
c
c
      subroutine test_aplusgas
      use action
      use atoms
      use energi
      use virial
      implicit none
      integer nat
      real*8 energy,e,ref_e,ref_ei
      real*8 eps_e,eps_g,eps_v
      real*8, allocatable :: derivs(:,:)
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      logical skiptest
      character*240 rpath
      character*(*) tname
      character*(*) tpre
      parameter (tname='aplusgas')
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'amoebaplus'))  return
      call pushdir ('file/aplusgas')
      call loadfix ('aplusgas','gas.key')
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('aplusgas','aplusgas.1.txt',rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total AMOEBA+ gas energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,tpre//tname//' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,tpre//tname//' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,tpre//tname//' virial (v1)')
c
c     level 3  --  total energy via analysis path
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,
     &                  tpre//tname//' analysis (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

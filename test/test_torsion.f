c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_torsion  --  torsional angle tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_torsion" checks the torsional angle case exercised by
c     the tinker-gpu torsion.cpp test
c
c
      subroutine test_torsion
      use action
      use atoms
      use energi
      use tors
      use virial
      implicit none
      integer nat,refcnt
      real*8 energy,e,ref_e,ref_ei,refeng
      real*8 eps_e,eps_g,eps_v
      real*8, allocatable :: derivs(:,:)
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      logical skiptest
      character*240 rpath
      character*(*) tname
      character*(*) tpre
      parameter (tname='torsion_trpcage')
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'amoeba'))  return
      call pushdir ('file/torsion')
      call loadfix ('trpcage','trpcage.key')
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('torsion','torsion.txt',rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Torsional Angle',refeng,refcnt)
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     setup and level 0  --  total torsional angle energy
c
      call assert_int (ntors,refcnt,tpre//tname//' count')
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
      call assert_real (et,refeng,eps_e,tpre//tname//' torsion (v0)')
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
      call check_engcnt (rpath,'Torsional Angle',et,net,eps_e,
     &                   tpre//tname//' torsion (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

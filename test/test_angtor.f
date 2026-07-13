c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_angtor  --  angle-torsion FF tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_angtor" checks the angle-torsion case exercised by the
c     tinker-gpu angtor.cpp test
c
c
      subroutine test_angtor
      implicit none
c
c
      call test_angtor_a4
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_angtor_a4  --  A4 angle-torsion test  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_angtor_a4" checks the AMOEBA nucleic acid angle-torsion
c     interaction for A4 against the upstream tinker-gpu reference
c
c
      subroutine test_angtor_a4
      use action
      use angtor
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
      character*240 rpath
      character*(*) tname
      character*(*) tpre
      parameter (tname='angtor_a4')
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'amoeba'))  return
c
c     load the angle-torsion-only fixture and reference values
c
      call pushdir ('file/angtor')
      call loadfix ('a4')
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('angtor','angtor.txt',rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Angle-Torsion',refeng,refcnt)
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     setup and level 0  --  total angle-torsion energy
c
      call assert_int (nangtor,refcnt,tpre//tname//' count')
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
      call assert_real (eat,refeng,eps_e,tpre//tname//' angtor (v0)')
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
      call check_engcnt (rpath,'Angle-Torsion',eat,neat,eps_e,
     &                   tpre//tname//' angtor (v3)')
c
c     clean up
c
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_improp  --  improper dihedral tests  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_improp" checks the improper dihedral cases exercised by
c     the tinker-gpu improp.cpp tests
c
c
      subroutine test_improp
      implicit none
c
c
      call test_improp_trpcage
      call test_improp_bdna3
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_improp_trpcage  --  CHARMM19 impropers  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_improp_trpcage
      implicit none
c
c
      call test_improp_case ('trp_charmm','improp.1.txt',
     &                       'improp_trpcage',1.0d-4,3.0d-4,1.0d-3)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_improp_bdna3  --  CHARMM36 impropers  ##
c     ##                                                        ##
c     ############################################################
c
c
      subroutine test_improp_bdna3
      implicit none
c
c
      call test_improp_case ('bdna3','improp.2.txt',
     &                       'improp_bdna3',1.0d-4,1.0d-3,1.0d-3)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_improp_case  --  run improper case  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine test_improp_case (base,reffile,tname,
     &                             eps_e,eps_g,eps_v)
      use action
      use atoms
      use energi
      use improp
      use virial
      implicit none
      integer nat,refcnt
      real*8 energy,e,ref_e,ref_ei,refeng
      real*8 eps_e,eps_g,eps_v
      real*8, allocatable :: derivs(:,:)
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      logical skiptest
      character*(*) base,reffile,tname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'charmm'))  return
      call pushdir ('file/improp')
      call loadfix (base,base//'.key')
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('improp',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Improper Dihedral',refeng,refcnt)
c
c     setup and level 0  --  total improper dihedral energy
c
      call assert_int (niprop,refcnt,tpre//tname//' count')
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
      call assert_real (eid,refeng,eps_e,tpre//tname//' improp (v0)')
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
      call check_engcnt (rpath,'Improper Dihedral',eid,neid,eps_e,
     &                   tpre//tname//' improp (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

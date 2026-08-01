c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ################################################
c     ##                                            ##
c     ##  subroutine test_vdw14  --  1-4 vdW tests  ##
c     ##                                            ##
c     ################################################
c
c
c     "test_vdw14" checks the 1-4 Lennard-Jones cases exercised by
c     the tinker-gpu vdw14.cpp tests
c
c
      subroutine test_vdw14
      implicit none
c
c
      call test_vdw14_nopbc
      call test_vdw14_cutoff
      return
      end
c
c
c     ##################################################
c     ##                                              ##
c     ##  subroutine test_vdw14_nopbc  --  no cutoff  ##
c     ##                                              ##
c     ##################################################
c
c
      subroutine test_vdw14_nopbc
      implicit none
c
c
      call test_vdw14_case ('vdw14.1','vdw14.1.txt',
     &                      'vdw14_nopbc',.false.)
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine test_vdw14_cutoff  --  pbc cutoff  ##
c     ##                                                ##
c     ####################################################
c
c
      subroutine test_vdw14_cutoff
      implicit none
c
c
      call test_vdw14_case ('vdw14.2','vdw14.2.txt',
     &                      'vdw14_cutoff_list',.true.)
      call test_vdw14_case ('vdw14.2','vdw14.2.txt',
     &                      'vdw14_cutoff_nolist',.false.)
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  subroutine test_vdw14_case  --  run 1-4 vdW  ##
c     ##                                               ##
c     ###################################################
c
c
      subroutine test_vdw14_case (key,reffile,tname,uselist)
      use action
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
      logical skiptest,uselist
      character*(*) key,reffile,tname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'charmm'))  return
      call pushdir ('file/vdw14')
      if (uselist) then
         call loadfix_keyadd ('trp_charmm',key//'.key',
     &                        'neighbor-list')
      else
         call loadfix ('trp_charmm',key//'.key')
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('vdw14',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Van der Waals',refeng,refcnt)
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total 1-4 van der Waals energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
      call assert_real (ev,refeng,eps_e,tpre//tname//' vdw14 (v0)')
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
      call check_engcnt (rpath,'Van der Waals',ev,nev,eps_e,
     &                   tpre//tname//' vdw14 (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

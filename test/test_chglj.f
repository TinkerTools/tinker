c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_chglj  --  charge plus Lennard-Jones  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_chglj" checks the charge plus Lennard-Jones cases
c     exercised by the tinker-gpu chglj.cpp tests
c
c
      subroutine test_chglj
      implicit none
c
c
      call test_chglj_nopbc
      call test_chglj_cutoff
      return
      end
c
c
c     ##################################################
c     ##                                              ##
c     ##  subroutine test_chglj_nopbc  --  no cutoff  ##
c     ##                                              ##
c     ##################################################
c
c
      subroutine test_chglj_nopbc
      implicit none
c
c
      call test_chglj_case ('chglj.1','chglj.1.txt',
     &                      'chglj_nopbc',.false.)
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine test_chglj_cutoff  --  pbc cutoff  ##
c     ##                                                ##
c     ####################################################
c
c
      subroutine test_chglj_cutoff
      implicit none
c
c
      call test_chglj_case ('chglj.2','chglj.2.txt',
     &                      'chglj_cutoff_list',.true.)
      call test_chglj_case ('chglj.2','chglj.2.txt',
     &                      'chglj_cutoff_nolist',.false.)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_chglj_case  --  run charge plus LJ  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine test_chglj_case (key,reffile,tname,uselist)
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
      logical skiptest,uselist
      character*(*) key,reffile,tname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'charmm'))  return
      call pushdir ('file/chglj')
      if (uselist) then
         call loadfix_keyadd ('trp_charmm',key//'.key',
     &                        'neighbor-list')
      else
         call loadfix ('trp_charmm',key//'.key')
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('chglj',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      eps_e = 1.0d-4
      eps_g = 5.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total charge plus Lennard-Jones energy
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
c     level 3  --  energy analysis
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

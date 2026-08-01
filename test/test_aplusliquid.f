c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_aplusliquid  --  AMOEBA+ liquid tests  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_aplusliquid" checks the AMOEBA+ liquid cases exercised
c     by the tinker-gpu aplusliquid.cpp tests
c
c
      subroutine test_aplusliquid
      implicit none
c
c
      call test_aplusliquid_ewald
      call test_aplusliquid_nonewald
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_aplusliquid_ewald  --  Ewald liquid  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_aplusliquid_ewald
      implicit none
c
c
      call test_aplusliquid_case ('liquid_ewald','aplusliquid.1.txt',
     &                           'aplusliquid_ewald_list',.true.,
     &                           1.0d-4,1.0d-4,1.0d-3)
      call test_aplusliquid_case ('liquid_ewald','aplusliquid.1.txt',
     &                           'aplusliquid_ewald_nolist',.false.,
     &                           1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_aplusliquid_nonewald  --  non-Ewald liquid  ##
c     ##                                                              ##
c     ##################################################################
c
c
      subroutine test_aplusliquid_nonewald
      implicit none
c
c
      call test_aplusliquid_case ('liquid','aplusliquid.2.txt',
     &                           'aplusliquid_nonewald_list',.true.,
     &                           1.0d-4,1.0d-4,1.0d-3)
      call test_aplusliquid_case ('liquid','aplusliquid.2.txt',
     &                           'aplusliquid_nonewald_nolist',
     &                           .false.,1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_aplusliquid_case  --  run liquid case  ##
c     ##                                                         ##
c     #############################################################
c
c
      subroutine test_aplusliquid_case (key,reffile,tname,uselist,
     &                                  eps_e,eps_g,eps_v)
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
      if (skiptest(tpre//tname,'amoebaplus'))  return
      call pushdir ('file/aplusliquid')
      if (uselist) then
         call loadfix_keyadd ('tetramer',key//'.key','neighbor-list')
      else
         call loadfix ('tetramer',key//'.key')
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('aplusliquid',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
c
c     level 0  --  total AMOEBA+ liquid energy
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
c     level 3  --  energy and per-term analysis
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,
     &                  tpre//tname//' analysis (v3)')
      call check_engcnt (rpath,'Van der Waals',ev,nev,eps_e,
     &                   tpre//tname//' vdw (v3)')
      call check_engcnt (rpath,'Atomic Multipoles',em,nem,eps_e,
     &                   tpre//tname//' mpole (v3)')
      call check_engcnt (rpath,'Polarization',ep,nep,eps_e,
     &                   tpre//tname//' polar (v3)')
      call check_engcnt (rpath,'Charge Transfer',ect,nect,eps_e,
     &                   tpre//tname//' chgtrn (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_rephippo  --  HIPPO repulsion tests  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_rephippo" checks the HIPPO Pauli repulsion cases
c     exercised by the tinker-gpu rephippo.cpp tests
c
c
      subroutine test_rephippo
      implicit none
c
c
      call test_rephippo_c5h12
      call test_rephippo_water
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine test_rephippo_c5h12  --  HIPPO repuls  ##
c     ##                                                    ##
c     ########################################################
c
c
      subroutine test_rephippo_c5h12
      implicit none
c
c
      call test_rephippo_case ('c5h12acnh2','rep','rephippo.1.txt',
     &                         'rephippo_c5h12',.true.)
      call test_rephippo_case ('c5h12acnh2','rep','rephippo.1.txt',
     &                         'rephippo_c5h12',.false.)
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine test_rephippo_water  --  water repuls  ##
c     ##                                                    ##
c     ########################################################
c
c
      subroutine test_rephippo_water
      implicit none
c
c
      call test_rephippo_case ('h2o10','repwater','rephippo.2.txt',
     &                         'rephippo_water',.true.)
      call test_rephippo_case ('h2o10','repwater','rephippo.2.txt',
     &                         'rephippo_water',.false.)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_rephippo_case  --  run HIPPO repuls  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_rephippo_case (base,key,reffile,tname,uselist)
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
      character*(*) base,key,reffile,tname
      character*48 cname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (uselist) then
         cname = tpre//trim(tname)//'_list'
      else
         cname = tpre//trim(tname)//'_nolist'
      end if
      if (skiptest(cname,'hippo'))  return
      call pushdir ('file/rephippo')
      if (uselist) then
         call loadfix_keyadd (base,key//'.key','neighbor-list')
      else
         call loadfix (base,key//'.key')
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('rephippo',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Repulsion',refeng,refcnt)
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total HIPPO repulsion energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,trim(cname)//' energy (v0)')
      call assert_real (er,refeng,eps_e,trim(cname)//' repuls (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,trim(cname)//' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,trim(cname)//' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,trim(cname)//' virial (v1)')
c
c     level 3  --  energy and per-term analysis
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,
     &                  trim(cname)//' analysis (v3)')
      call check_engcnt (rpath,'Repulsion',er,ner,eps_e,
     &                   trim(cname)//' repuls (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

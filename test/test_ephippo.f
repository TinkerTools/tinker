c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_ephippo  --  HIPPO polarization tests  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_ephippo" checks the HIPPO polarization cases exercised
c     by the tinker-gpu ephippo.cpp tests
c
c
      subroutine test_ephippo
      implicit none
c
c
      call test_ephippo_ewald
      call test_ephippo_nonewald
      return
      end
c
c
c     ######################################################
c     ##                                                  ##
c     ##  subroutine test_ephippo_ewald  --  Ewald polar  ##
c     ##                                                  ##
c     ######################################################
c
c
      subroutine test_ephippo_ewald
      implicit none
c
c
      call test_ephippo_case ('ewald','ephippo.1.txt',
     &                       'ephippo_ewald',.true.,
     &                       1.0d-4,1.0d-4,1.0d-3)
      call test_ephippo_case ('ewald','ephippo.1.txt',
     &                       'ephippo_ewald',.false.,
     &                       1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_ephippo_nonewald  --  non-Ewald polar  ##
c     ##                                                         ##
c     #############################################################
c
c
      subroutine test_ephippo_nonewald
      implicit none
c
c
      call test_ephippo_case ('newald','ephippo.2.txt',
     &                       'ephippo_nonewald',.true.,1.0d-4,1.0d-4,
     &                       1.0d-3)
      call test_ephippo_case ('newald','ephippo.2.txt',
     &                       'ephippo_nonewald',.false.,1.0d-4,1.0d-4,
     &                       1.0d-3)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_ephippo_case  --  run HIPPO polar  ##
c     ##                                                     ##
c     #########################################################
c
c
      subroutine test_ephippo_case (key,reffile,tname,uselist,
     &                              eps_e,eps_g,eps_v)
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
      call pushdir ('file/ephippo')
      if (uselist) then
         call loadfix_keyadd ('dmso',key//'.key','neighbor-list')
      else
         call loadfix ('dmso',key//'.key')
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('ephippo',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Polarization',refeng,refcnt)
c
c     level 0  --  total HIPPO polarization energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,trim(cname)//' energy (v0)')
      call assert_real (ep,refeng,eps_e,trim(cname)//' polar (v0)')
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
      call check_engcnt (rpath,'Polarization',ep,nep,eps_e,
     &                   trim(cname)//' polar (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

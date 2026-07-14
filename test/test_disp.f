c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################
c     ##                                              ##
c     ##  subroutine test_disp  --  dispersion tests  ##
c     ##                                              ##
c     ##################################################
c
c
c     "test_disp" checks the HIPPO dispersion cases exercised by
c     the tinker-gpu disp.cpp tests
c
c
      subroutine test_disp
      implicit none
c
c
      call test_disp_nonewald
      call test_disp_dewald
      return
      end
c
c
c     ######################################################
c     ##                                                  ##
c     ##  subroutine test_disp_nonewald  --  non-DEwald  ##
c     ##                                                  ##
c     ######################################################
c
c
      subroutine test_disp_nonewald
      implicit none
c
c
      call test_disp_case ('ndewald','disp.1.txt','disp_nonewald',
     &                     1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     ###############################################
c     ##                                           ##
c     ##  subroutine test_disp_dewald  --  DEwald  ##
c     ##                                           ##
c     ###############################################
c
c
      subroutine test_disp_dewald
      implicit none
c
c
      call test_disp_case ('dewald','disp.2.txt','disp_dewald',
     &                     1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine test_disp_case  --  run dispersion  ##
c     ##                                                 ##
c     #####################################################
c
c
      subroutine test_disp_case (key,reffile,tname,
     &                           eps_e,eps_g,eps_v)
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
      logical skiptest
      character*(*) key,reffile,tname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'hippo'))  return
      call pushdir ('file/disp')
      call loadfix ('c5h12acnh2',key//'.key')
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('disp',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Dispersion',refeng,refcnt)
c
c     level 0  --  total dispersion energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
      call assert_real (edsp,refeng,eps_e,tpre//tname//' disp (v0)')
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
      call check_engcnt (rpath,'Dispersion',edsp,nedsp,eps_e,
     &                   tpre//tname//' disp (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_expol  --  exchange polarization tests  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_expol" checks the exchange polarization cases exercised
c     by the tinker-gpu expol.cpp tests
c
c
      subroutine test_expol
      implicit none
c
c
      call test_expol_case ('NaCl','expol','expol.1.txt',
     &                      'expol_nacl',1.0d-4,1.0d-4,1.0d-3)
      call test_expol_case ('Nawater','expol','expol.2.txt',
     &                      'expol_nawater',1.0d-4,1.0d-4,1.0d-3)
      call test_expol_case ('Clwater','expol','expol.3.txt',
     &                      'expol_clwater',1.0d-4,2.0d-4,1.0d-3)
      call test_expol_case ('water2Na2Clbox','expol','expol.4.txt',
     &                      'expol_box',1.0d-4,1.0d-4,1.0d-3)
      call test_expol_case ('crys','expols2','expol.5.txt',
     &                      'expol_s2',1.0d-4,1.0d-4,1.0d-3)
      call test_expol_case ('crys','expolg','expol.6.txt',
     &                      'expol_g',1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_expol_case  --  run exchange polar  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine test_expol_case (base,key,reffile,tname,
     &                            eps_e,eps_g,eps_v)
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
      character*(*) base,key,reffile,tname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'hippo'))  return
      call pushdir ('file/expol')
      call loadfix (base,key//'.key')
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('expol',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Polarization',refeng,refcnt)
c
c     level 0  --  total exchange polarization energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
      call assert_real (ep,refeng,eps_e,tpre//tname//' expol (v0)')
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
      call check_engcnt (rpath,'Polarization',ep,nep,eps_e,
     &                   tpre//tname//' expol (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

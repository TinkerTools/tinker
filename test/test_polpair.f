c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_polpair  --  polarization pair tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_polpair" checks the polarization pair cases exercised by
c     the tinker-gpu polpair.cpp tests
c
c
      subroutine test_polpair
      implicit none
c
c
      call test_polpair_ewald
      call test_polpair_nonewald
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine test_polpair_ewald  --  Ewald pair  ##
c     ##                                                 ##
c     #####################################################
c
c
      subroutine test_polpair_ewald
      implicit none
c
c
      call test_polpair_case ('nacl_ewald','polpair.1.txt',
     &                        'polpair_ewald')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_polpair_nonewald  --  non-Ewald pair  ##
c     ##                                                        ##
c     ############################################################
c
c
      subroutine test_polpair_nonewald
      implicit none
c
c
      call test_polpair_case ('nacl','polpair.2.txt',
     &                        'polpair_nonewald')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_polpair_case  --  run polarization pair  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine test_polpair_case (key,reffile,tname)
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
      logical skiptest
      character*(*) key,reffile,tname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'amoeba'))  return
      call pushdir ('file/polpair')
      call loadfix ('nacl',key//'.key')
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('polpair',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total polarization pair energy
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
c     level 3  --  total energy via analysis path
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

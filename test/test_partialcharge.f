c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_partialcharge  --  partial charge tests  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_partialcharge" checks the partial charge cases exercised
c     by the tinker-gpu partialcharge.cpp tests
c
c
      subroutine test_partialcharge
      implicit none
c
c
      call test_partialcharge_nonewald
      call test_partialcharge_ewald
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_partialcharge_nonewald  --  no cutoff  ##
c     ##                                                         ##
c     #############################################################
c
c
      subroutine test_partialcharge_nonewald
      implicit none
c
c
      call test_partialcharge_case ('partialcharge.1',
     &                              'partialcharge.1.txt',
     &                              'partialcharge_nonewald',.false.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_partialcharge_ewald  --  pbc and ewald  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_partialcharge_ewald
      implicit none
c
c
      call test_partialcharge_case ('partialcharge.2',
     &                              'partialcharge.2.txt',
     &                              'partialcharge_ewald_list',.true.)
      call test_partialcharge_case ('partialcharge.2',
     &                              'partialcharge.2.txt',
     &                              'partialcharge_ewald_nolist',
     &                              .false.)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_partialcharge_case  --  run charges  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_partialcharge_case (key,reffile,tname,uselist)
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
      call pushdir ('file/partialcharge')
      if (uselist) then
         call loadfix_keyadd ('trp_charmm',key//'.key',
     &                        'neighbor-list')
      else
         call loadfix ('trp_charmm',key//'.key')
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('partialcharge',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Charge-Charge',refeng,refcnt)
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total charge-charge energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
      call assert_real (ec,refeng,eps_e,tpre//tname//' charge (v0)')
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
      call check_engcnt (rpath,'Charge-Charge',ec,nec,eps_e,
     &                   tpre//tname//' charge (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

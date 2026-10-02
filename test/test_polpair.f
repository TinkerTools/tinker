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
      call test_polpair_thole
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
      call test_polpair_case ('nacl','nacl_ewald','polpair.1.txt',
     &                        'polpair_ewald',.true.,.false.)
      call test_polpair_case ('nacl','nacl_ewald','polpair.1.txt',
     &                        'polpair_ewald',.false.,.false.)
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
      call test_polpair_case ('nacl','nacl','polpair.2.txt',
     &                        'polpair_nonewald',.true.,.false.)
      call test_polpair_case ('nacl','nacl','polpair.2.txt',
     &                        'polpair_nonewald',.false.,.false.)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_polpair_thole  --  zero Thole damping  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_polpair_thole" checks systems whose Thole pair values are
c     zero or vanishingly small; a zero pair value means the pair is
c     undamped, while 1.0d-30 damps it fully, and such fully damped
c     pairs must not count as interactions
c
c
      subroutine test_polpair_thole
      implicit none
c
c
c     AMOEBA water with an H-H Thole pair value of zero, set either
c     by "polpair" or by a zero H Thole that falls back to the O value
c     in the O-H pairs and leaves the H-H pairs at zero
c
      call test_polpair_case ('../mutate/water','water_zeropair',
     &                        'polpair.3.txt','polpair_zeropair',
     &                        .true.,.true.)
      call test_polpair_case ('../mutate/water','water_zeropair',
     &                        'polpair.3.txt','polpair_zeropair',
     &                        .false.,.true.)
      call test_polpair_case ('../mutate/water','water_zerothole',
     &                        'polpair.3.txt','polpair_zerothole',
     &                        .true.,.true.)
      call test_polpair_case ('../mutate/water','water_zerothole',
     &                        'polpair.3.txt','polpair_zerothole',
     &                        .false.,.true.)
c
c     AMOEBA water with fully damped H-H pairs
c
      call test_polpair_case ('../mutate/water','water_tinypair',
     &                        'polpair.4.txt','polpair_tinypair',
     &                        .true.,.true.)
      call test_polpair_case ('../mutate/water','water_tinypair',
     &                        'polpair.4.txt','polpair_tinypair',
     &                        .false.,.true.)
c
c     Dang-Chang water dimer with undamped point dipoles; its O sites
c     carry neither multipoles nor polarizability, so they are dropped
c     from the multipole sites and atom and site numbers differ
c
      call test_polpair_case ('dang','dang','polpair.5.txt',
     &                        'polpair_dang',.true.,.true.)
      call test_polpair_case ('dang','dang','polpair.5.txt',
     &                        'polpair_dang',.false.,.true.)
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
      subroutine test_polpair_case (base,key,reffile,tname,uselist,
     &                              chkcnt)
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
      logical skiptest,uselist,chkcnt
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
      if (skiptest(cname,'amoeba'))  return
      call pushdir ('file/polpair')
      if (uselist) then
         call loadfix_keyadd (base,key//'.key','neighbor-list')
      else
         call loadfix (base,key//'.key')
      end if
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
      call assert_real (esum,ref_e,eps_e,trim(cname)//' energy (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,trim(cname)//' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,trim(cname)//' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,trim(cname)//' virial (v1)')
c
c     level 3  --  total energy via analysis path
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,
     &                  trim(cname)//' analysis (v3)')
      if (chkcnt) then
         call check_engcnt (rpath,'Polarization',ep,nep,eps_e,
     &                      trim(cname)//' polar (v3)')
      end if
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

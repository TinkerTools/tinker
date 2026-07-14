c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################
c     ##                                              ##
c     ##  subroutine test_nacl  --  NaCl force tests  ##
c     ##                                              ##
c     ##################################################
c
c
c     "test_nacl" checks the NaCl subtypes exercised by the
c     tinker-gpu nacl.cpp tests
c
c
      subroutine test_nacl
      implicit none
c
c
      call test_nacl_ehal_noswitch
      call test_nacl_ehal_switch_near_cut
      call test_nacl_ehal_switch_near_off
      call test_nacl_ehal_evcorr_full
      call test_nacl_ehal_evcorr_lambda
      call test_nacl_empole_nonewald
      call test_nacl_empole_ewald
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_nacl_ehal_noswitch  --  ehal no switch  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_nacl_ehal_noswitch
      implicit none
c
c
      call test_nacl_vdw_case ('nacl1','nacl_vdw.key',
     &   'nacl_ehal_noswitch',51.4242d0,1,184.4899d0,
     &   -405.878d0,0.0d0,0.0d0)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_nacl_ehal_switch_near_cut  --  near cut  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine test_nacl_ehal_switch_near_cut
      implicit none
c
c
      call test_nacl_vdw_case ('nacl2','nacl_vdw.key',
     &   'nacl_ehal_switch_near_cut',25.8420d0,1,149.0904d0,
     &   -354.835d0,0.0d0,0.0d0)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_nacl_ehal_switch_near_off  --  near off  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine test_nacl_ehal_switch_near_off
      implicit none
c
c
      call test_nacl_vdw_case ('nacl3','nacl_vdw.key',
     &   'nacl_ehal_switch_near_off',4.8849d0,1,127.6639d0,
     &   -319.160d0,0.0d0,0.0d0)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_nacl_ehal_evcorr_full  --  full evcorr  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_nacl_ehal_evcorr_full
      implicit none
c
c
      call test_nacl_vdw_case ('nacl1','nacl_vdw_corr.key',
     &   'nacl_ehal_evcorr_full',52.0802d0,1,184.4899d0,
     &   -408.398d0,-2.521d0,-2.521d0)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_nacl_ehal_evcorr_lambda  --  lambda evcorr  ##
c     ##                                                              ##
c     ##################################################################
c
c
      subroutine test_nacl_ehal_evcorr_lambda
      implicit none
c
c
      call test_nacl_vdw_case ('nacl1','nacl_vdw_lambda.key',
     &   'nacl_ehal_evcorr_lambda',25.8947d0,1,81.9570d0,
     &   -182.721d0,-2.415d0,-2.415d0)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_nacl_empole_nonewald  --  non-Ewald mpole  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine test_nacl_empole_nonewald
      implicit none
c
c
      call test_nacl_mpole_case ('nacl1','nacl_mpole.key',
     &   'nacl_empole_nonewald',-150.9381d0,1,1.0d-4)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_nacl_empole_ewald  --  Ewald multipole  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_nacl_empole_ewald
      implicit none
c
c
      call test_nacl_mpole_case ('nacl4','nacl_mpole_ewald.key',
     &   'nacl_empole_ewald',-0.3146d0,2,5.0d-4)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine test_nacl_vdw_case  --  run NaCl vdW  ##
c     ##                                                   ##
c     #######################################################
c
c
      subroutine test_nacl_vdw_case (base,key,tname,refeng,refcnt,
     &                               gx,vxx,vyy,vzz)
      use action
      use atoms
      use energi
      use virial
      implicit none
      integer refcnt
      real*8 energy,e,refeng,gx,vxx,vyy,vzz
      real*8 eps,refv(3,3)
      real*8, allocatable :: derivs(:,:),refg(:,:)
      logical skiptest
      character*(*) base,key,tname
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'amoeba'))
     &   return
      call pushdir ('file/nacl')
      call loadfix (base,key)
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call zero_ref2 (refg,refv,n)
      refg(1,1) = gx
      refg(1,2) = -gx
      refv(1,1) = vxx
      refv(2,2) = vyy
      refv(3,3) = vzz
      eps = 1.0d-3
c
c     level 0  --  buffered 14-7 van der Waals energy
c
      e = energy ()
      call assert_real (ev,refeng,eps,tpre//tname//' vdw (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (ev,refeng,eps,tpre//tname//' vdw (v1)')
      call assert_grad (derivs,refg,n,eps,tpre//tname//' grad (v1)')
      call assert_grad (vir,refv,3,eps,tpre//tname//' virial (v1)')
c
c     level 3  --  energy and per-term analysis
c
      call analysis (e)
      call assert_real (ev,refeng,eps,tpre//tname//' vdw (v3)')
      call assert_int (nev,refcnt,tpre//tname//' vdw count (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_nacl_mpole_case  --  run NaCl mpole  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_nacl_mpole_case (base,key,tname,refeng,refcnt,
     &                                 eps)
      use action
      use atoms
      use energi
      use virial
      implicit none
      integer refcnt
      real*8 energy,e,refeng,eps,refv(3,3)
      real*8, allocatable :: derivs(:,:),refg(:,:)
      logical skiptest,ewald
      character*(*) base,key,tname
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'amoeba'))
     &   return
      call pushdir ('file/nacl')
      call loadfix (base,key)
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call zero_ref2 (refg,refv,n)
      ewald = index(tname,'ewald') .gt. 0
     &        .and. index(tname,'nonewald') .eq. 0
      if (ewald) then
         refg(1,1) = 0.1731d0
         refg(2,1) = 0.1921d0
         refg(3,1) = 0.2103d0
         refg(1,2) = -0.1668d0
         refg(2,2) = -0.1917d0
         refg(3,2) = -0.2081d0
         refv(1,1) = 0.059d0
         refv(2,1) = -0.284d0
         refv(3,1) = -0.311d0
         refv(1,2) = -0.284d0
         refv(2,2) = 0.105d0
         refv(3,2) = -0.342d0
         refv(1,3) = -0.311d0
         refv(2,3) = -0.342d0
         refv(3,3) = 0.155d0
      else
         refg(1,1) = -68.6082d0
         refg(1,2) = 68.6082d0
         refv(1,1) = 150.938d0
      end if
c
c     level 0  --  atomic multipole energy
c
      e = energy ()
      call assert_real (em,refeng,eps,tpre//tname//' mpole (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (em,refeng,eps,tpre//tname//' mpole (v1)')
      call assert_grad (derivs,refg,n,eps,tpre//tname//' grad (v1)')
      call assert_grad (vir,refv,3,eps,tpre//tname//' virial (v1)')
c
c     level 3  --  energy and per-term analysis
c
      call analysis (e)
      call assert_real (em,refeng,eps,tpre//tname//' mpole (v3)')
      call assert_int (nem,refcnt,tpre//tname//' mpole count (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  subroutine zero_ref2  --  clear 2-atom refs  ##
c     ##                                               ##
c     ###################################################
c
c
      subroutine zero_ref2 (refg,refv,n)
      implicit none
      integer n,i,j
      real*8 refg(3,n),refv(3,3)
c
c
      do i = 1, n
         do j = 1, 3
            refg(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            refv(j,i) = 0.0d0
         end do
      end do
      return
      end

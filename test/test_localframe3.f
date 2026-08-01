c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_localframe3  --  mpole plus polar  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_localframe3" checks the combined multipole/polarization
c     local-frame cases exercised by tinker-gpu localframe3.cpp
c
c
      subroutine test_localframe3
      implicit none
c
c
      call test_localframe3_nonewald
      call test_localframe3_ewald
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_localframe3_nonewald  --  non-Ewald  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_localframe3_nonewald
      implicit none
c
c
      call test_localframe3_case ('localframe3_nonewald.key',
     &   'localframe3_nonewald','localframe3.1.txt',120,.true.)
      call test_localframe3_case ('localframe3_nonewald.key',
     &   'localframe3_nonewald','localframe3.1.txt',120,.false.)
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine test_localframe3_ewald  --  Ewald PME  ##
c     ##                                                    ##
c     ########################################################
c
c
      subroutine test_localframe3_ewald
      implicit none
c
c
      call test_localframe3_case ('localframe3_ewald.key',
     &   'localframe3_ewald','localframe3.2.txt',138,.true.)
      call test_localframe3_case ('localframe3_ewald.key',
     &   'localframe3_ewald','localframe3.2.txt',138,.false.)
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine test_localframe3_case  --  run emplar  ##
c     ##                                                    ##
c     ########################################################
c
c
      subroutine test_localframe3_case (key,tname,rfile,refcnt,
     &                                  uselist)
      use action
      use atoms
      use energi
      use virial
      implicit none
      integer nat,refcnt
      integer refcntm,refcntp
      real*8 energy,e,ref_e,ref_ei,refem,refep
      real*8 eps_e,eps_g,eps_v
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: refg(:,:)
      real*8 refv(3,3)
      logical skiptest,uselist
      character*(*) key,tname,rfile
      character*72 cname
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
      call pushdir ('file/localframe')
      if (uselist) then
         call loadfix_keyadd ('localframe2',key,'neighbor-list')
      else
         call loadfix ('localframe2',key)
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('localframe',rfile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Atomic Multipoles',refem,refcntm)
      call load_engcnt (rpath,'Polarization',refep,refcntp)
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total multipole plus polarization energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,trim(cname)//
     &                  ' energy (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,trim(cname)//
     &                  ' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,trim(cname)//
     &                  ' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,trim(cname)//
     &                  ' virial (v1)')
c
c     level 3  --  energy and per-term analysis
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,trim(cname)//
     &                  ' energy (v3)')
      call assert_real (em,refem,eps_e,trim(cname)//
     &                  ' mpole (v3)')
      call assert_real (ep,refep,eps_e,trim(cname)//
     &                  ' polar (v3)')
      call assert_int (nem,refcnt,trim(cname)//
     &                 ' mpole count (v3)')
      call assert_int (nep,refcnt,trim(cname)//
     &                 ' polar count (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

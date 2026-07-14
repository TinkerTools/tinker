c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_extfield  --  external field tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_extfield" checks the external field cases exercised by
c     the tinker-gpu extfield.cpp tests, except the dynamic frequency
c     molecular dynamics case
c
c
      subroutine test_extfield
      implicit none
c
c
      call test_extfield_mpole
      call test_extfield_polar
      call test_extfield_mpolar
      call test_extfield_vdwpchg
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_extfield_mpole  --  multipole field  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_extfield_mpole
      implicit none
c
c
      call test_extfield_case ('water4_Na_shifted',
     &   'water4_Na_shifted_mpole.key','extfield.1.txt',
     &   'extfield_mpole','Atomic Multipoles','amoeba')
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_extfield_polar  --  polarization field  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_extfield_polar
      implicit none
c
c
      call test_extfield_case ('water4_Na_shifted',
     &   'water4_Na_shifted_polar.key','extfield.2.txt',
     &   'extfield_polar','Polarization','amoeba')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_extfield_mpolar  --  mpole + polar field  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine test_extfield_mpolar
      implicit none
c
c
      call test_extfield_case ('water4_Na_shifted',
     &   'water4_Na_shifted_mpolar.key','extfield.3.txt',
     &   'extfield_mpolar',' ','amoeba')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_extfield_vdwpchg  --  vdw + charge field  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine test_extfield_vdwpchg
      implicit none
c
c
      call test_extfield_case ('amber','amber.key','extfield.4.txt',
     &   'extfield_vdwpchg',' ','amber')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_extfield_case  --  run external field  ##
c     ##                                                         ##
c     #############################################################
c
c
      subroutine test_extfield_case (base,key,reffile,tname,term,fftag)
      use action
      use atoms
      use energi
      use extfld
      use units
      use virial
      implicit none
      integer nat,refcnt,i
      real*8 energy,e,ref_e,ref_ei,refeng
      real*8 eps_e,eps_g,eps_v
      real*8, allocatable :: derivs(:,:)
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      logical skiptest,hasterm
      character*(*) base,key,reffile,tname,term,fftag
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,fftag))  return
      call pushdir ('file/extfield')
      call loadfix (base,key)
      use_exfld = .true.
      exfld(1) = 150.0d0 / elefield
      exfld(2) = -300.0d0 / elefield
      exfld(3) = 450.0d0 / elefield
      do i = 1, 3
         texfld(i) = exfld(i)
      end do
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('extfield',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      hasterm = len_trim(term) .ne. 0
      if (hasterm) then
         call load_engcnt (rpath,term,refeng,refcnt)
         ref_e = refeng
      end if
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  external field energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
      if (trim(term) .eq. 'Atomic Multipoles')
     &   call assert_real (em,refeng,eps_e,
     &                     tpre//tname//' mpole (v0)')
      if (trim(term) .eq. 'Polarization')
     &   call assert_real (ep,refeng,eps_e,
     &                     tpre//tname//' polar (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,tpre//tname//' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,tpre//tname//' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,tpre//tname//' virial (v1)')
c
c     level 3  --  total energy and optional per-term analysis
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,tpre//tname//' analysis (v3)')
      if (trim(term) .eq. 'Atomic Multipoles')
     &   call check_engcnt (rpath,'Atomic Multipoles',em,nem,
     &                      eps_e,tpre//tname//' mpole (v3)')
      if (trim(term) .eq. 'Polarization')
     &   call check_engcnt (rpath,'Polarization',ep,nep,
     &                      eps_e,tpre//tname//' polar (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

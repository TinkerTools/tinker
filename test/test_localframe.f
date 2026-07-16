c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_localframe  --  local frame mpole  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_localframe" checks the local-frame multipole and
c     polarization cases exercised by tinker-gpu localframe.cpp
c
c
      subroutine test_localframe
      implicit none
c
c
      call test_localframe_nonewald
      call test_localframe_ewald
      call test_localframe_polar_nonewald
      call test_localframe_polar_ewald
      call test_localframe_dfield_nonewald
      call test_localframe_ufield_nonewald
      call test_localframe_induce_nonewald
      call test_localframe_dfield_ewald
      call test_localframe_ufield_ewald
      call test_localframe_induce_ewald
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_localframe_nonewald  --  gas mpole  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine test_localframe_nonewald
      implicit none
c
c
      call test_localframe_mpole_case ('localframe_mpole.key',
     &   'localframe_nonewald','localframe.1.txt',.false.,18,
     &   1.0d-4,1.0d-4,1.0d-3)
      call test_localframe_mpole_case ('localframe_mpole.key',
     &   'localframe_nonewald','localframe.1.txt',.true.,18,
     &   1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine test_localframe_ewald  --  PME mpole  ##
c     ##                                                   ##
c     #######################################################
c
c
      subroutine test_localframe_ewald
      implicit none
c
c
      call test_localframe_mpole_case ('localframe_mpole_ewald.key',
     &   'localframe_ewald','localframe.2.txt',.true.,22,
     &   1.0d-4,1.0d-4,1.0d-3)
      call test_localframe_mpole_case ('localframe_mpole_ewald.key',
     &   'localframe_ewald','localframe.2.txt',.false.,22,
     &   1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_localframe_polar_nonewald  --  gas polar  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine test_localframe_polar_nonewald
      implicit none
c
c
      call test_localframe_polar_case ('localframe_polar.key',
     &   'localframe_polar_nonewald','localframe.3.txt',.false.,18,
     &   1.0d-4,1.0d-4,1.0d-3)
      call test_localframe_polar_case ('localframe_polar.key',
     &   'localframe_polar_nonewald','localframe.3.txt',.true.,18,
     &   1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_localframe_polar_ewald  --  PME polar  ##
c     ##                                                         ##
c     #############################################################
c
c
      subroutine test_localframe_polar_ewald
      implicit none
c
c
      call test_localframe_polar_case ('localframe_polar_ewald.key',
     &   'localframe_polar_ewald','localframe.4.txt',.true.,18,
     &   1.0d-4,1.0d-4,1.0d-3)
      call test_localframe_polar_case ('localframe_polar_ewald.key',
     &   'localframe_polar_ewald','localframe.4.txt',.false.,18,
     &   1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_localframe_dfield_nonewald  --  gas field  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine test_localframe_dfield_nonewald
      implicit none
c
c
      call test_localframe_dfield_case ('localframe_polar.key',
     &   'localframe_dfield_nonewald','localframe_dfield.1.txt',
     &   'localframe_dfield.2.txt',.false.,.false.)
      call test_localframe_dfield_case ('localframe_polar.key',
     &   'localframe_dfield_nonewald','localframe_dfield.1.txt',
     &   'localframe_dfield.2.txt',.true.,.false.)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_localframe_ufield_nonewald  --  gas field  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine test_localframe_ufield_nonewald
      implicit none
c
c
      call test_localframe_ufield_case ('localframe_polar.key',
     &   'localframe_ufield_nonewald','localframe_ufield.1.txt',
     &   'localframe_ufield.2.txt',.false.,.false.)
      call test_localframe_ufield_case ('localframe_polar.key',
     &   'localframe_ufield_nonewald','localframe_ufield.1.txt',
     &   'localframe_ufield.2.txt',.true.,.false.)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_localframe_induce_nonewald  --  gas dipole  ##
c     ##                                                              ##
c     ##################################################################
c
c
      subroutine test_localframe_induce_nonewald
      implicit none
c
c
      call test_localframe_induce_case ('localframe_polar.key',
     &   'localframe_induce_nonewald','localframe_induce.1.txt',
     &   'localframe_induce.2.txt',.false.)
      call test_localframe_induce_case ('localframe_polar.key',
     &   'localframe_induce_nonewald','localframe_induce.1.txt',
     &   'localframe_induce.2.txt',.true.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_localframe_dfield_ewald  --  PME field  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_localframe_dfield_ewald
      implicit none
c
c
      call test_localframe_dfield_case ('localframe_polar_ewald.key',
     &   'localframe_dfield_ewald','localframe_dfield.3.txt',
     &   'localframe_dfield.4.txt',.true.,.true.)
      call test_localframe_dfield_case ('localframe_polar_ewald.key',
     &   'localframe_dfield_ewald','localframe_dfield.3.txt',
     &   'localframe_dfield.4.txt',.false.,.true.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_localframe_ufield_ewald  --  PME field  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_localframe_ufield_ewald
      implicit none
c
c
      call test_localframe_ufield_case ('localframe_polar_ewald.key',
     &   'localframe_ufield_ewald','localframe_ufield.3.txt',
     &   'localframe_ufield.4.txt',.true.,.true.)
      call test_localframe_ufield_case ('localframe_polar_ewald.key',
     &   'localframe_ufield_ewald','localframe_ufield.3.txt',
     &   'localframe_ufield.4.txt',.false.,.true.)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_localframe_induce_ewald  --  PME dipole  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine test_localframe_induce_ewald
      implicit none
c
c
      call test_localframe_induce_case ('localframe_polar_ewald.key',
     &   'localframe_induce_ewald','localframe_induce.3.txt',
     &   'localframe_induce.4.txt',.true.)
      call test_localframe_induce_case ('localframe_polar_ewald.key',
     &   'localframe_induce_ewald','localframe_induce.3.txt',
     &   'localframe_induce.4.txt',.false.)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_localframe_mpole_case  --  run mpole  ##
c     ##                                                        ##
c     ############################################################
c
c
      subroutine test_localframe_mpole_case (key,tname,rfile,uselist,
     &                                      ngrad,eps_e,eps_g,eps_v)
      use action
      use atoms
      use energi
      use virial
      implicit none
      integer nat,ngrad,refcnt
      real*8 energy,e,ref_e,ref_ei,refeng
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
         call loadfix_keyadd ('localframe',key,'neighbor-list')
      else
         call loadfix ('localframe',key)
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('localframe',rfile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Atomic Multipoles',refeng,refcnt)
c
c     level 0  --  multipole energy
c
      e = energy ()
      call assert_real (em,refeng,eps_e,trim(cname)//
     &                  ' mpole (v0)')
c
c     level 1  --  multipole energy, gradient and virial
c
      call gradient (e,derivs)
      call assert_real (em,refeng,eps_e,trim(cname)//
     &                  ' mpole (v1)')
      call assert_grad (derivs,refg,ngrad,eps_g,trim(cname)//
     &                  ' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,trim(cname)//
     &                  ' virial (v1)')
c
c     level 3  --  multipole analysis
c
      call analysis (e)
      call check_engcnt (rpath,'Atomic Multipoles',em,nem,eps_e,
     &                   trim(cname)//' mpole (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_localframe_polar_case  --  run polar  ##
c     ##                                                        ##
c     ############################################################
c
c
      subroutine test_localframe_polar_case (key,tname,rfile,uselist,
     &                                      ngrad,eps_e,eps_g,eps_v)
      use action
      use atoms
      use energi
      use virial
      implicit none
      integer nat,ngrad,refcnt
      real*8 energy,e,ref_e,ref_ei,refeng
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
         call loadfix_keyadd ('localframe',key,'neighbor-list')
      else
         call loadfix ('localframe',key)
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('localframe',rfile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Polarization',refeng,refcnt)
c
c     direct polarization energy, matching the upstream dot-product case
c
      if (uselist)  call nblist
      call epolar
      call assert_real (ep,refeng,eps_e,trim(cname)//
     &                  ' epolar direct')
c
c     level 0  --  polarization energy
c
      e = energy ()
      call assert_real (ep,refeng,eps_e,trim(cname)//
     &                  ' polar (v0)')
c
c     level 1  --  polarization energy, gradient and virial
c
      call gradient (e,derivs)
      call assert_real (ep,refeng,eps_e,trim(cname)//
     &                  ' polar (v1)')
      call assert_grad (derivs,refg,ngrad,eps_g,trim(cname)//
     &                  ' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,trim(cname)//
     &                  ' virial (v1)')
c
c     level 3  --  polarization analysis
c
      call analysis (e)
      call check_engcnt (rpath,'Polarization',ep,nep,eps_e,
     &                   trim(cname)//' polar (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_localframe_dfield_case  --  dfield  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine test_localframe_dfield_case (key,tname,rfiled,rfilep,
     &                                       uselist,useewald)
      use atoms
      use limits
      use polar
      implicit none
      integer nat
      logical skiptest,uselist,useewald
      real*8, allocatable :: refd(:,:)
      real*8, allocatable :: refp(:,:)
      character*(*) key,tname,rfiled,rfilep
      character*72 cname
      character*240 rpathd,rpathp
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
         call loadfix_keyadd ('localframe',key,'neighbor-list')
      else
         call loadfix ('localframe',key)
      end if
      allocate (refd(3,n))
      allocate (refp(3,n))
      call refpath ('localframe',rfiled,rpathd)
      call refpath ('localframe',rfilep,rpathp)
      call read_gradient (rpathd,n,refd,nat)
      call read_gradient (rpathp,n,refp,nat)
      call chkpole
      call rotpole ('MPOLE')
      if (use_list)  call nblist
      if (useewald) then
         call dfield0c (udir,udirp)
      else
         call dfield0a (udir,udirp)
      end if
      call assert_grad (udir,refd,n,1.0d-4,trim(cname)//
     &                  ' udir')
      call assert_grad (udirp,refp,n,1.0d-4,trim(cname)//
     &                  ' udirp')
      deallocate (refd)
      deallocate (refp)
      call popdir
      call final
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_localframe_ufield_case  --  ufield  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine test_localframe_ufield_case (key,tname,rfiled,rfilep,
     &                                       uselist,useewald)
      use atoms
      use limits
      use polar
      implicit none
      integer i,j,nat
      logical skiptest,uselist,useewald
      real*8, allocatable :: refd(:,:)
      real*8, allocatable :: refp(:,:)
      character*(*) key,tname,rfiled,rfilep
      character*72 cname
      character*240 rpathd,rpathp
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
         call loadfix_keyadd ('localframe',key,'neighbor-list')
      else
         call loadfix ('localframe',key)
      end if
      allocate (refd(3,n))
      allocate (refp(3,n))
      call refpath ('localframe',rfiled,rpathd)
      call refpath ('localframe',rfilep,rpathp)
      call read_gradient (rpathd,n,refd,nat)
      call read_gradient (rpathp,n,refp,nat)
      call chkpole
      call rotpole ('MPOLE')
      if (use_list)  call nblist
      if (useewald)  call dfield0c (udir,udirp)
      do i = 1, n
         do j = 1, 3
            uind(j,i) = 0.1d0*dble(i) + 0.03d0*dble(j)
            uinp(j,i) = 0.1d0*dble(i) - 0.03d0*dble(j)
         end do
      end do
      if (useewald) then
         call ufield0c (udir,udirp)
      else
         call ufield0a (udir,udirp)
      end if
      call assert_grad (udir,refd,n,1.0d-4,trim(cname)//
     &                  ' udir')
      call assert_grad (udirp,refp,n,1.0d-4,trim(cname)//
     &                  ' udirp')
      deallocate (refd)
      deallocate (refp)
      call popdir
      call final
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_localframe_induce_case  --  induced  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_localframe_induce_case (key,tname,rfiled,rfilep,
     &                                       uselist)
      use atoms
      use limits
      use polar
      use units
      implicit none
      integer i,j,nat
      logical skiptest,uselist
      real*8, allocatable :: refd(:,:)
      real*8, allocatable :: refp(:,:)
      real*8, allocatable :: dipd(:,:)
      real*8, allocatable :: dipp(:,:)
      character*(*) key,tname,rfiled,rfilep
      character*72 cname
      character*240 rpathd,rpathp
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
         call loadfix_keyadd ('localframe',key,'neighbor-list')
      else
         call loadfix ('localframe',key)
      end if
      allocate (refd(3,n))
      allocate (refp(3,n))
      allocate (dipd(3,n))
      allocate (dipp(3,n))
      call refpath ('localframe',rfiled,rpathd)
      call refpath ('localframe',rfilep,rpathp)
      call read_gradient (rpathd,n,refd,nat)
      call read_gradient (rpathp,n,refp,nat)
      call chkpole
      call rotpole ('MPOLE')
      if (use_list)  call nblist
      call induce
      do i = 1, n
         do j = 1, 3
            dipd(j,i) = debye * uind(j,i)
            dipp(j,i) = debye * uinp(j,i)
         end do
      end do
      call assert_grad (dipd,refd,n,1.0d-4,trim(cname)//
     &                  ' uind')
      call assert_grad (dipp,refp,n,1.0d-4,trim(cname)//
     &                  ' uinp')
      deallocate (refd)
      deallocate (refp)
      deallocate (dipd)
      deallocate (dipp)
      call popdir
      call final
      return
      end

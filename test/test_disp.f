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
c     the tinker-gpu disp.cpp tests, then the replicate cell loops
c
c
      subroutine test_disp
      implicit none
c
c
      call test_disp_nonewald
      call test_disp_dewald
      call test_disp_replica
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
     &                     .true.,1.0d-4,1.0d-4,1.0d-3)
      call test_disp_case ('ndewald','disp.1.txt','disp_nonewald',
     &                     .false.,1.0d-4,1.0d-4,1.0d-3)
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
     &                     .true.,1.0d-4,1.0d-4,1.0d-3)
      call test_disp_case ('dewald','disp.2.txt','disp_dewald',
     &                     .false.,1.0d-4,1.0d-4,1.0d-3)
c
c     no dispersion on the methyl H, so its excluded pairs have no
c     energy and must not be taken away from the interaction count
c
      call test_disp_case ('dewald_zeroc6','disp.3.txt',
     &                     'disp_zeroc6',.true.,1.0d-4,1.0d-4,1.0d-3)
      call test_disp_case ('dewald_zeroc6','disp.3.txt',
     &                     'disp_zeroc6',.false.,1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_disp_replica  --  replicate cell paths  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_disp_replica" checks that the energy, gradient and analysis
c     routines give the same dispersion energy for two DMSO molecules
c     in a box small enough for the cutoff to reach replicate cells,
c     where images of bonded partners lie inside the cutoff, with and
c     without Ewald summation and with one molecule decoupled under
c     the soft core
c
c
      subroutine test_disp_replica
      implicit none
c
c
      call test_disp_replcase ('replica_cut','disp_replica_cut')
      call test_disp_replcase ('replica_ewald','disp_replica_ewald')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_disp_replcase  --  one replicate system  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine test_disp_replcase (key,tname)
      use action
      use atoms
      use energi
      implicit none
      integer nsave
      real*8 energy,e,e0,e1
      real*8, allocatable :: derivs(:,:)
      logical skiptest
      character*(*) key,tname
      character*48 cname
c
c
      cname = 'test_'//trim(tname)
      if (skiptest(cname,'hippo'))  return
      call pushdir ('file/disp')
      call loadfix ('../testgrad/dmso2',key//'.key')
      allocate (derivs(3,n))
c
c     compare the energy from each of the three routines
c
      e = energy ()
      e0 = edsp
      call gradient (e,derivs)
      e1 = edsp
      call assert_real (e1,e0,1.0d-8,trim(cname)//' gradient')
      call analysis (e)
      nsave = nedsp
      call assert_real (edsp,e0,1.0d-8,trim(cname)//' analysis')
      call assert_logical (nsave.gt.0,.true.,
     &                     trim(cname)//' analysis count')
      deallocate (derivs)
      call popdir
      call final
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
      subroutine test_disp_case (key,reffile,tname,uselist,
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
      call pushdir ('file/disp')
      if (uselist) then
         call loadfix_keyadd ('c5h12acnh2',key//'.key',
     &                        'neighbor-list')
      else
         call loadfix ('c5h12acnh2',key//'.key')
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('disp',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Dispersion',refeng,refcnt)
c
c     level 0  --  total dispersion energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,trim(cname)//' energy (v0)')
      call assert_real (edsp,refeng,eps_e,trim(cname)//' disp (v0)')
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
      call check_engcnt (rpath,'Dispersion',edsp,nedsp,eps_e,
     &                   trim(cname)//' disp (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

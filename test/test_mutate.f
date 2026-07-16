c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################
c     ##                                              ##
c     ##  subroutine test_mutate  --  mutation tests  ##
c     ##                                              ##
c     ##################################################
c
c
c     "test_mutate" checks AMOEBA mutation energy, gradient, virial
c     and named energy component regressions
c
c
      subroutine test_mutate
      implicit none
c
c
      call test_mutate_mv
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_mv  --  mutation lambda-scan cases  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_mv" runs the nine water mutation fixtures; cases
c     001-006 isolate electrostatics at three ele-lambda values with
c     Ewald on and off, while cases 007-009 isolate van der Waals at
c     three vdw-lambda values
c
c
      subroutine test_mutate_mv
      implicit none
c
c
c
c     the checkm/checkp/checkv flags select which level 3 energy
c     components are verified (Atomic Multipoles, Polarization, and
c     Van der Waals); cases 001-006 check multipole and polarization,
c     cases 007-009 check van der Waals.  the final canlist flag
c     enables the neighbor-list variant.  it is off for the no-Ewald
c     (ne) cases 002/004/006, whose full-cutoff electrostatics require
c     periodic replicas that Tinker forbids combining with a pairwise
c     neighbor list; every other fixture keeps its cutoffs plus the
c     list buffer under half the box and so runs both variants
c
c                                              checkm  checkp  checkv  canlist
      call test_mutate_mv_case ('001_water_ye_m10.key',
     &   '001_water_ye_m10.txt','001_water_ye_m10',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_mv_case ('002_water_ne_m10.key',
     &   '002_water_ne_m10.txt','002_water_ne_m10',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_mv_case ('003_water_ye_m05.key',
     &   '003_water_ye_m05.txt','003_water_ye_m05',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_mv_case ('004_water_ne_m05.key',
     &   '004_water_ne_m05.txt','004_water_ne_m05',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_mv_case ('005_water_ye_m00.key',
     &   '005_water_ye_m00.txt','005_water_ye_m00',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_mv_case ('006_water_ne_m00.key',
     &   '006_water_ne_m00.txt','006_water_ne_m00',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_mv_case ('007_water_v10.key',
     &   '007_water_v10.txt','007_water_v10',
     &   .false., .false., .true.,  .true.)
      call test_mutate_mv_case ('008_water_v05.key',
     &   '008_water_v05.txt','008_water_v05',
     &   .false., .false., .true.,  .true.)
      call test_mutate_mv_case ('009_water_v00.key',
     &   '009_water_v00.txt','009_water_v00',
     &   .false., .false., .true.,  .true.)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_mutate_mv_case  --  one mutation case  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_mutate_mv_case" runs a single mutation fixture with and
c     without the neighbor-list keyword; for each neighbor-list
c     variant the force field is built once, then the level 0/1/3
c     checks are repeated twice before teardown, giving four passes
c     per fixture; the "checkm", "checkp" and "checkv" flags select
c     which level 3 energy components are verified (Atomic Multipoles,
c     Polarization and Van der Waals); the "canlist" flag records
c     whether the neighbor-list variant is compatible with this
c     fixture (false for the no-Ewald cases)
c
c
      subroutine test_mutate_mv_case
     &   (key,ref,cname,checkm,checkp,checkv,canlist)
      use action
      use atoms
      use energi
      use virial
      implicit none
      integer nat,irun,ilist,nlist
      real*8 energy,e,ref_e,ref_ei
      real*8 eps_e,eps_g,eps_v,refv(3,3)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: refg(:,:)
      logical skiptest,checkm,checkp,checkv,canlist,uselist
      character*(*) key,ref,cname
      character*240 rpath,pre
      character*8 rtag
c
c
      if (skiptest('test_'//trim(cname),'mutate'))  return
c
c     run each fixture without the neighbor-list keyword, and with it
c     when the fixture supports a pairwise neighbor list
c
      nlist = 2
      if (.not. canlist)  nlist = 1
      do ilist = 1, nlist
         uselist = (ilist .eq. 2)
c
c     set up the force field once for this neighbor-list variant
c
         call pushdir ('file/mutate')
         if (uselist) then
            call loadfix_keyadd ('water',key,'neighbor-list')
            pre = 'test_'//trim(cname)//' list'
         else
            call loadfix ('water',key)
            pre = 'test_'//trim(cname)//' nolist'
         end if
         allocate (derivs(3,n))
         allocate (refg(3,n))
         call refpath ('mutate',ref,rpath)
         call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
         eps_e = 1.0d-4
         eps_g = 1.0d-4
         eps_v = 1.0d-3
c
c     repeat the level 0/1/3 checks twice against the built system
c
         do irun = 1, 2
            if (irun .eq. 1) then
               rtag = ' run1'
            else
               rtag = ' run2'
            end if
c
c     level 0  --  total potential energy
c
            e = energy ()
            call assert_real (esum,ref_e,eps_e,
     &                        trim(pre)//' energy (v0)'//trim(rtag))
c
c     level 1  --  total energy, Cartesian gradient and virial
c
            call gradient (e,derivs)
            call assert_real (esum,ref_e,eps_e,
     &                        trim(pre)//' grad-e (v1)'//trim(rtag))
            call assert_grad (derivs,refg,n,eps_g,
     &                        trim(pre)//' grad (v1)'//trim(rtag))
            call assert_grad (vir,refv,3,eps_v,
     &                        trim(pre)//' virial (v1)'//trim(rtag))
c
c     level 3  --  total and named AMOEBA energy components
c
            call analysis (e)
            call assert_real (esum,ref_e,eps_e,
     &                        trim(pre)//' analysis (v3)'//trim(rtag))
            if (checkm) then
               call check_engcnt (rpath,'Atomic Multipoles',em,nem,
     &                       eps_e,trim(pre)//' mpole (v3)'//trim(rtag))
            end if
            if (checkp) then
               call check_engcnt (rpath,'Polarization',ep,nep,
     &                       eps_e,trim(pre)//' polar (v3)'//trim(rtag))
            end if
            if (checkv) then
               call check_engcnt (rpath,'Van der Waals',ev,nev,
     &                       eps_e,trim(pre)//' vdw (v3)'//trim(rtag))
            end if
         end do
c
c     clean up this neighbor-list variant
c
         deallocate (derivs)
         deallocate (refg)
         call popdir
         call final
      end do
      return
      end

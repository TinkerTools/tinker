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
      call test_mutate_mp
      call test_mutate_ast
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
      call test_mutate_fixed ('001_water_ye_m10.key',
     &   '001_water_ye_m10.txt','001_water_ye_m10',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('002_water_ne_m10.key',
     &   '002_water_ne_m10.txt','002_water_ne_m10',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('003_water_ye_m05.key',
     &   '003_water_ye_m05.txt','003_water_ye_m05',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('004_water_ne_m05.key',
     &   '004_water_ne_m05.txt','004_water_ne_m05',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('005_water_ye_m00.key',
     &   '005_water_ye_m00.txt','005_water_ye_m00',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('006_water_ne_m00.key',
     &   '006_water_ne_m00.txt','006_water_ne_m00',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('007_water_v10.key',
     &   '007_water_v10.txt','007_water_v10',
     &   .false., .false., .true.,  .true.)
      call test_mutate_fixed ('008_water_v05.key',
     &   '008_water_v05.txt','008_water_v05',
     &   .false., .false., .true.,  .true.)
      call test_mutate_fixed ('009_water_v00.key',
     &   '009_water_v00.txt','009_water_v00',
     &   .false., .false., .true.,  .true.)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_mutate_mp  --  electrostatic lambda cases  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_mutate_mp" runs the twenty water mutation fixtures that
c     scan the electrostatic lambda values; cases 010-015 keep only the
c     multipole term at three ele-lambda values with Ewald on and off,
c     cases 016-021 keep only the polarization term at three pol-lambda
c     values, and cases 022-029 leave both terms active while varying
c     ele-lambda and pol-lambda together
c
c
      subroutine test_mutate_mp
      implicit none
c
c
      call test_mutate_fixed ('010_water_ye_m10.key',
     &   '010_water_ye_m10.txt','010_water_ye_m10',
     &   .true.,  .false., .false., .true.)
      call test_mutate_fixed ('011_water_ne_m10.key',
     &   '011_water_ne_m10.txt','011_water_ne_m10',
     &   .true.,  .false., .false., .false.)
      call test_mutate_fixed ('012_water_ye_m05.key',
     &   '012_water_ye_m05.txt','012_water_ye_m05',
     &   .true.,  .false., .false., .true.)
      call test_mutate_fixed ('013_water_ne_m05.key',
     &   '013_water_ne_m05.txt','013_water_ne_m05',
     &   .true.,  .false., .false., .false.)
      call test_mutate_fixed ('014_water_ye_m00.key',
     &   '014_water_ye_m00.txt','014_water_ye_m00',
     &   .true.,  .false., .false., .true.)
      call test_mutate_fixed ('015_water_ne_m00.key',
     &   '015_water_ne_m00.txt','015_water_ne_m00',
     &   .true.,  .false., .false., .false.)
      call test_mutate_fixed ('016_water_ye_p10.key',
     &   '016_water_ye_p10.txt','016_water_ye_p10',
     &   .false., .true.,  .false., .true.)
      call test_mutate_fixed ('017_water_ne_p10.key',
     &   '017_water_ne_p10.txt','017_water_ne_p10',
     &   .false., .true.,  .false., .false.)
      call test_mutate_fixed ('018_water_ye_p05.key',
     &   '018_water_ye_p05.txt','018_water_ye_p05',
     &   .false., .true.,  .false., .true.)
      call test_mutate_fixed ('019_water_ne_p05.key',
     &   '019_water_ne_p05.txt','019_water_ne_p05',
     &   .false., .true.,  .false., .false.)
      call test_mutate_fixed ('020_water_ye_p00.key',
     &   '020_water_ye_p00.txt','020_water_ye_p00',
     &   .false., .true.,  .false., .true.)
      call test_mutate_fixed ('021_water_ne_p00.key',
     &   '021_water_ne_p00.txt','021_water_ne_p00',
     &   .false., .true.,  .false., .false.)
      call test_mutate_fixed ('022_water_ye_m10p05.key',
     &   '022_water_ye_m10p05.txt','022_water_ye_m10p05',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('023_water_ne_m10p05.key',
     &   '023_water_ne_m10p05.txt','023_water_ne_m10p05',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('024_water_ye_m05p10.key',
     &   '024_water_ye_m05p10.txt','024_water_ye_m05p10',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('025_water_ne_m05p10.key',
     &   '025_water_ne_m05p10.txt','025_water_ne_m05p10',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('026_water_ye_m05p00.key',
     &   '026_water_ye_m05p00.txt','026_water_ye_m05p00',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('027_water_ne_m05p00.key',
     &   '027_water_ne_m05p00.txt','027_water_ne_m05p00',
     &   .true.,  .true.,  .false., .false.)
      call test_mutate_fixed ('028_water_ye_m00p05.key',
     &   '028_water_ye_m00p05.txt','028_water_ye_m00p05',
     &   .true.,  .true.,  .false., .true.)
      call test_mutate_fixed ('029_water_ne_m00p05.key',
     &   '029_water_ne_m00p05.txt','029_water_ne_m00p05',
     &   .true.,  .true.,  .false., .false.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_mutate_ast  --  single topology lambda  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_mutate_ast" runs the eleven absolute single topology water
c     fixtures 030-040, each carrying the "lambda-deriv" keyword, and
c     drives them through "test_mutate_calc" with the level 4 lambda
c     derivative checks enabled; cases 030-035 keep only the multipole
c     term at three ele-lambda values with Ewald on and off, cases
c     036-038 keep only the van der Waals term at three vdw-lambda
c     values, and cases 039-040 leave the multipole and polarization
c     terms active with Ewald on and off; the no-Ewald cases cannot use
c     a pairwise neighbor list
c
c
      subroutine test_mutate_ast
      implicit none
c
c
      call test_mutate_calc ('water2','030_water_ast_ye_m10.key',
     &   '030_water_ast_ye_m10.txt','030_water_ast_ye_m10',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','031_water_ast_ne_m10.key',
     &   '031_water_ast_ne_m10.txt','031_water_ast_ne_m10',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','032_water_ast_ye_m05.key',
     &   '032_water_ast_ye_m05.txt','032_water_ast_ye_m05',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','033_water_ast_ne_m05.key',
     &   '033_water_ast_ne_m05.txt','033_water_ast_ne_m05',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','034_water_ast_ye_m00.key',
     &   '034_water_ast_ye_m00.txt','034_water_ast_ye_m00',
     &   .true.,  .false., .false., .true.,  .true.)
      call test_mutate_calc ('water2','035_water_ast_ne_m00.key',
     &   '035_water_ast_ne_m00.txt','035_water_ast_ne_m00',
     &   .true.,  .false., .false., .false., .true.)
      call test_mutate_calc ('water2','036_water_ast_v10.key',
     &   '036_water_ast_v10.txt','036_water_ast_v10',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','037_water_ast_v05.key',
     &   '037_water_ast_v05.txt','037_water_ast_v05',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','038_water_ast_v00.key',
     &   '038_water_ast_v00.txt','038_water_ast_v00',
     &   .false., .false., .true.,  .true.,  .true.)
      call test_mutate_calc ('water2','039_water_ast_ye_m05p10.key',
     &   '039_water_ast_ye_m05p10.txt','039_water_ast_ye_m05p10',
     &   .true.,  .true.,  .false., .true.,  .true.)
      call test_mutate_calc ('water2','040_water_ast_ne_m05p10.key',
     &   '040_water_ast_ne_m05p10.txt','040_water_ast_ne_m05p10',
     &   .true.,  .true.,  .false., .false., .true.)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_mutate_fixed  --  one mutation case  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_mutate_fixed" runs a single mutation fixture with and
c     without the neighbor-list keyword, checking the level 0/1/3
c     energy, gradient, virial and named component regressions; it is
c     a thin wrapper around "test_mutate_calc" with the level 4 lambda
c     derivative checks disabled, and backs the "test_mutate_mv" and
c     "test_mutate_mp" case lists
c
c
      subroutine test_mutate_fixed
     &   (key,ref,cname,checkm,checkp,checkv,canlist)
      implicit none
      logical checkm,checkp,checkv,canlist
      character*(*) key,ref,cname
c
c
      call test_mutate_calc ('water',key,ref,cname,checkm,checkp,
     &                       checkv,canlist,.false.)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_mutate_calc  --  one mutation case  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_mutate_calc" runs a single mutation fixture with and
c     without the neighbor-list keyword; for each neighbor-list
c     variant the force field is built once, then the checks are
c     repeated twice before teardown, giving four passes per fixture;
c     the "checkm", "checkp" and "checkv" flags select which level 3
c     energy components are verified (Atomic Multipoles, Polarization
c     and Van der Waals); the "canlist" flag records whether the
c     neighbor-list variant is compatible with this fixture (false for
c     the no-Ewald cases); when "dolmda" is set, the level 4 checks
c     also verify the analytical lambda derivatives, second lambda
c     derivatives, per-atom lambda gradient and dV/dL tensor against
c     the reference read by "load_lmdaref", reusing the values already
c     produced by the level 1 gradient call
c
c
      subroutine test_mutate_calc
     &   (base,key,ref,cname,checkm,checkp,checkv,canlist,dolmda)
      use action
      use atoms
      use dlmda
      use energi
      use virial
      implicit none
      integer nat,natl,irun,ilist,nlist
      real*8 energy,e,ref_e,ref_ei
      real*8 eps_e,eps_g,eps_v,eps_l,refv(3,3)
      real*8 ref_dedl(4),ref_d2edl2(4),ref_dvdl(3,3)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: refg(:,:)
      real*8, allocatable :: ref_lg(:,:)
      logical skiptest,checkm,checkp,checkv,canlist,uselist,dolmda
      character*(*) base,key,ref,cname
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
            call loadfix_keyadd (base,key,'neighbor-list')
            pre = 'test_'//trim(cname)//' list'
         else
            call loadfix (base,key)
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
c     read the reference lambda derivatives for the level 4 checks
c
         if (dolmda) then
            allocate (ref_lg(3,n))
            call load_lmdaref (rpath,n,ref_dedl,ref_d2edl2,ref_lg,
     &                         ref_dvdl,natl)
            eps_l = 1.0d-4
         end if
c
c     repeat the checks twice against the built system
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
c     level 4  --  lambda derivatives from the level 1 gradient call
c
            if (dolmda) then
               call assert_real (dedl,ref_dedl(1),eps_l,
     &                        trim(pre)//' dE/dL (v4)'//trim(rtag))
               call assert_real (devdl,ref_dedl(2),eps_l,
     &                        trim(pre)//' dEV/dL (v4)'//trim(rtag))
               call assert_real (demdl,ref_dedl(3),eps_l,
     &                        trim(pre)//' dEM/dL (v4)'//trim(rtag))
               call assert_real (depdl,ref_dedl(4),eps_l,
     &                        trim(pre)//' dEP/dL (v4)'//trim(rtag))
               call assert_real (d2edl2,ref_d2edl2(1),eps_l,
     &                        trim(pre)//' d2E/dL2 (v4)'//trim(rtag))
               call assert_real (d2evdl2,ref_d2edl2(2),eps_l,
     &                        trim(pre)//' d2EV/dL2 (v4)'//trim(rtag))
               call assert_real (d2emdl2,ref_d2edl2(3),eps_l,
     &                        trim(pre)//' d2EM/dL2 (v4)'//trim(rtag))
               call assert_real (d2epdl2,ref_d2edl2(4),eps_l,
     &                        trim(pre)//' d2EP/dL2 (v4)'//trim(rtag))
               call assert_grad (dfsumdl,ref_lg,n,eps_g,
     &                        trim(pre)//' lgrad (v4)'//trim(rtag))
               call assert_grad (dvirdl,ref_dvdl,3,eps_v,
     &                        trim(pre)//' dV/dL (v4)'//trim(rtag))
            end if
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
         if (dolmda)  deallocate (ref_lg)
         call popdir
         call final
      end do
      return
      end

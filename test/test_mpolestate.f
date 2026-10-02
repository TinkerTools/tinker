c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_mpolestate  --  multipole state tests  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_mpolestate" checks the multipole and induced dipole state
c     the energy routines leave behind, as exercised by the tinker-gpu
c     mpolestate.cpp tests
c
c
      subroutine test_mpolestate
      implicit none
c
c
      call test_mpolestate_physical
      call test_mpolestate_induced
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_mpolestate_physical  --  physical frames  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_mpolestate_physical" checks that energy, gradient and
c     analysis calls leave "pole" at the electrostatic lambda, with
c     charge flux, since "moments" and later energy terms rotate it;
c     each fixture leaves a different last pass behind: polarization
c     at a lambda other than the electrostatic one, its lambda
c     derivative, the dual topology end states, multipoles alone, a
c     pinned polarization lambda, and charge flux
c
c     "rpole" is not held to the physical state, since every reader
c     rotates "pole" first; "empole4" rotates the unscaled "poleorig"
c     into it and scales on the fly, and leaves it that way when no
c     polarization pass follows
c
c
      subroutine test_mpolestate_physical
      implicit none
      integer i,j
      integer nfix
      parameter (nfix=8)
      character*10 dir(nfix)
      character*18 base(nfix)
      character*36 key(nfix)
      character*14 name(nfix)
      character*10 tags(nfix)
      character*1 lvl(3)
      data dir   / 'mpolestate', 'testlmda', 'mutate', 'mutate',
     &             'mpolestate', 'mpolestate', 'mpolestate',
     &             'mpolestate' /
      data base  / '../mutate/water2', 'water2', 'water2', 'water2',
     &             '../mutate/water2', '../mutate/water2',
     &             '../rephippo/h2o10', '../mutate/water2' /
      data key   / 'rels_prng_nodl.key',
     &             '27_water_rels_lig1_st_prng_l088.key',
     &             '203_water_rels_st_l085.key',
     &             '136_water_rels_ye_l085.key',
     &             'emast_mponly.key', 'polpinned.key',
     &             'hippo_cflux.key', 'rels_prng_l100.key' /
      data name  / 'rels_prng_nodl', 'rels_prng_dl', 'rels_st',
     &             'epdt', 'emast_mponly', 'polpinned', 'hippo_cflux',
     &             'rels_prng_l100' /
      data tags  / 'mutate', 'mutate', 'mutate', 'mutate', 'mutate',
     &             'mutate', 'hippo', 'mutate' /
      data lvl   / '0', '1', '3' /
c
c
      do i = 1, nfix
         do j = 1, 3
            call test_mpolestate_case (trim(dir(i)),trim(base(i)),
     &                                 trim(key(i)),
     &                                 'test_mpolestate_'//
     &                                 trim(name(i))//'_v'//lvl(j),
     &                                 'mpolestate,'//trim(tags(i)),
     &                                 lvl(j))
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_mpolestate_case  --  one physical check  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_mpolestate_case" evaluates one fixture at the level named
c     by "lvl" and compares "pole" with the state rebuilt from scratch
c     the way "altemdt" restores the electrostatic lambda
c
c
      subroutine test_mpolestate_case (dir,base,key,tname,tags,lvl)
      use atoms
      use mpole
      use potent
      implicit none
      integer i,j
      real*8 energy,e
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: pole0(:,:)
      logical skiptest
      character*(*) dir,base,key
      character*(*) tname,tags,lvl
c
c
      if (skiptest(tname,tags))  return
      call pushdir ('file/'//dir)
      call loadfix (base,key)
      allocate (derivs(3,n))
      allocate (pole0(maxpole,n))
c
c     evaluate the energy at the requested level
c
      if (lvl .eq. '0') then
         e = energy ()
      else if (lvl .eq. '1') then
         call gradient (e,derivs)
      else
         call analysis (e)
      end if
      do i = 1, n
         do j = 1, maxpole
            pole0(j,i) = pole(j,i)
         end do
      end do
c
c     rebuild the electrostatic lambda state at these coordinates
c
      if (use_mutate)  call altelec
      if (use_chgflx)  call alterchg
      call chkpole
      call assert_array2 (pole0,pole,maxpole,n,1.0d-10,
     &                    tname//' pole')
      deallocate (derivs)
      deallocate (pole0)
      call popdir
      call final
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_mpolestate_induced  --  reported dipoles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_mpolestate_induced" checks that single topology leaves the
c     induced dipoles solved at the polarization lambda in "uind" after
c     a gradient call, even though the lambda derivative passes run
c     after that solve; two solves of one state agree to a few 1.0d-9
c     here, far inside the polarization convergence, while a wrong
c     state or weighting is off by 1.0d-2 or more
c
c
      subroutine test_mpolestate_induced
      use atoms
      use dlmda
      use mutant
      use polar
      implicit none
      integer i,j
      real*8 e
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: uind0(:,:)
      logical skiptest,same
      character*(*) tname
      parameter (tname='test_mpolestate_induced_prst')
c
c
      if (skiptest(tname,'mpolestate,mutate'))  return
      call pushdir ('file/testlmda')
      call loadfix ('water2','27_water_rels_lig1_st_prng_l088.key')
      allocate (derivs(3,n))
      allocate (uind0(3,n))
      call assert_logical (use_prst,.true.,tname//' single topology')
      call assert_logical (abs(plambda-elambda).gt.1.0d-12,.true.,
     &                     tname//' lambdas differ')
      call gradient (e,derivs)
      do i = 1, n
         do j = 1, 3
            uind0(j,i) = uind(j,i)
         end do
      end do
c
c     solve the dipoles again at the polarization lambda
c
      call altepset (same)
      call induce
      if (.not. same)  call alteprst
      call assert_array2 (uind,uind0,3,n,1.0d-6,tname//' uind')
      deallocate (derivs)
      deallocate (uind0)
      call popdir
      call final
      return
      end

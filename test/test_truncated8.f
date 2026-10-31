c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_truncated8  --  truncated octahedron PBC  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_truncated8" checks the truncated octahedron periodic box
c     cases exercised by the tinker-gpu truncated8.cpp test: an argon
c     box with van der Waals only, and a water box with AMOEBA, HIPPO
c     with Ewald dispersion, AMOEBA+ and TIP3P
c
c
      subroutine test_truncated8
      implicit none
c
c
      call test_truncated8_argon
      call test_truncated8_water ('test_truncated8_amoeba',
     &                            'octamoeba.key','truncated8.2.txt')
      call test_truncated8_water ('test_truncated8_hippo',
     &                            'octhippo.key','truncated8.3.txt')
      call test_truncated8_water ('test_truncated8_aplus',
     &                            'octaplus.key','truncated8.4.txt')
      call test_truncated8_water ('test_truncated8_tip3p',
     &                            'octtip3p.key','truncated8.5.txt')
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_truncated8_argon  --  argon box vdw  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_truncated8_argon" checks van der Waals alone with a
c     neighbor list on the octahedral argon box
c
c
      subroutine test_truncated8_argon
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
      logical skiptest
      character*240 rpath
      character*(*) tname
      parameter (tname='test_truncated8')
c
c
      if (skiptest(tname,'amoeba'))  return
c
c     load the octahedral argon box and the reference values
c
      call pushdir ('file/truncated8')
      call loadfix ('arbox','arbox.key')
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('truncated8','truncated8.1.txt',rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Van der Waals',refeng,refcnt)
      eps_e = 1.0d-4
      eps_g = 2.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total van der Waals energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tname//' energy (v0)')
      call assert_real (ev,refeng,eps_e,tname//' vdw (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,tname//' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,tname//' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,tname//' virial (v1)')
c
c     level 3  --  energy and per-term analysis
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,tname//' analysis (v3)')
      call check_engcnt (rpath,'Van der Waals',ev,nev,eps_e,
     &                   tname//' vdw (v3)')
c
c     clean up
c
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
c     ##  subroutine test_truncated8_water  --  water box PBC  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_truncated8_water" checks the energy, gradient, virial and
c     every nonbonded term and count on the octahedral water box with
c     neighbor lists, polarization and Ewald sums under "key"
c
c     the references hold half the energy, interaction count and
c     virial of the cube of edge a containing the box twice, the second
c     copy shifted by (a/2,a/2,a/2), and the gradient rows of the first
c     copy; that cube uses only the orthogonal image code
c
c
      subroutine test_truncated8_water (tname,key,ref)
      use action
      use atoms
      use energi
      use potent
      use virial
      implicit none
      integer nat
      real*8 energy,e,ref_e,ref_ei
      real*8 eps_e,eps_g,eps_v
      real*8, allocatable :: derivs(:,:)
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      logical skiptest
      character*(*) tname,key,ref
      character*240 rpath
c
c
      if (skiptest(tname,'amoeba'))  return
c
c     load the octahedral water box and the reference values
c
      call pushdir ('file/truncated8')
      call loadfix ('octwater',key)
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('truncated8',ref,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      eps_e = 1.0d-4
      eps_g = 2.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tname//' energy (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,tname//' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,tname//' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,tname//' virial (v1)')
c
c     level 3  --  energy and each nonbonded term with its count
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,tname//' analysis (v3)')
      if (use_vdw)  call check_engcnt (rpath,'Van der Waals',ev,nev,
     &                                 eps_e,tname//' vdw (v3)')
      if (use_repel)  call check_engcnt (rpath,'Repulsion',er,ner,
     &                                   eps_e,tname//' repel (v3)')
      if (use_disp)  call check_engcnt (rpath,'Dispersion',edsp,nedsp,
     &                                  eps_e,tname//' disp (v3)')
      if (use_charge)  call check_engcnt (rpath,'Charge-Charge',ec,nec,
     &                                    eps_e,tname//' charge (v3)')
      if (use_mpole)  call check_engcnt (rpath,'Atomic Multipoles',
     &                                   em,nem,eps_e,
     &                                   tname//' mpole (v3)')
      if (use_polar)  call check_engcnt (rpath,'Polarization',ep,nep,
     &                                   eps_e,tname//' polar (v3)')
      if (use_chgtrn)  call check_engcnt (rpath,'Charge Transfer',
     &                                    ect,nect,eps_e,
     &                                    tname//' chgtrn (v3)')
c
c     clean up
c
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

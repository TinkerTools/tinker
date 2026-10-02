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
c     case exercised by the tinker-gpu truncated8.cpp test: an argon
c     box with van der Waals only and a neighbor list
c
c
      subroutine test_truncated8
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

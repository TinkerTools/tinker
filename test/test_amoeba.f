c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_amoeba  --  AMOEBA force field tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_amoeba" calls the Tinker unit tests that exercise the AMOEBA
c     polarizable force field; each routine loads a fixture from the
c     test tree and checks the energy, gradient and analysis levels
c     against a stored reference
c
c
      subroutine test_amoeba
      implicit none
c
c
      call test_vasopressin
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_vasopressin  --  AMOEBA energy and gradient ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "test_vasopressin" checks the AMOEBA potential energy, Cartesian
c     gradient and analysis of vasopressin against the reference output,
c     confirming that energy, gradient and analysis all agree
c
c
      subroutine test_vasopressin
      use action
      use atoms
      use energi
      use virial
      implicit none
      integer nat
      real*8 energy,e,ref_e,ref_ei
      real*8 eps_e,eps_g,eps_v
      real*8, allocatable :: derivs(:,:)
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      logical skiptest
      character*240 rpath
      character*(*) tname
      character*(*) tpre
      parameter (tname='vasopressin')
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'amoeba'))  return
c
c     load the structure and force field, and the reference values
c
      call pushdir ('file/'//tname)
      call loadfix (tname)
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath (tname,tname//'.txt',rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total potential energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
c
c     level 1  --  energy and Cartesian gradient
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,tpre//tname//' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,tpre//tname//' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,tpre//tname//' virial (v1)')
c
c     level 3  --  energy and per-term analysis
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,tpre//tname//' analysis (v3)')
      call check_engcnt (rpath,'Bond Stretching',eb,neb,eps_e,
     &                   tpre//tname//' bond (v3)')
      call check_engcnt (rpath,'Angle Bending',ea,nea,eps_e,
     &                   tpre//tname//' angle (v3)')
      call check_engcnt (rpath,'Stretch-Bend',eba,neba,eps_e,
     &                   tpre//tname//' strbnd (v3)')
      call check_engcnt (rpath,'Out-of-Plane Bend',eopb,neopb,eps_e,
     &                   tpre//tname//' opbend (v3)')
      call check_engcnt (rpath,'Torsional Angle',et,net,eps_e,
     &                   tpre//tname//' torsion (v3)')
      call check_engcnt (rpath,'Pi-Orbital Torsion',ept,nept,eps_e,
     &                   tpre//tname//' pitors (v3)')
      call check_engcnt (rpath,'Torsion-Torsion',ett,nett,eps_e,
     &                   tpre//tname//' tortor (v3)')
      call check_engcnt (rpath,'Van der Waals',ev,nev,eps_e,
     &                   tpre//tname//' vdw (v3)')
      call check_engcnt (rpath,'Atomic Multipoles',em,nem,eps_e,
     &                   tpre//tname//' mpole (v3)')
      call check_engcnt (rpath,'Polarization',ep,nep,eps_e,
     &                   tpre//tname//' polar (v3)')
c
c     clean up
c
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end

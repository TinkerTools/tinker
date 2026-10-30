c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_dodecahedron  --  rhombic dodecahedron PBC  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "test_dodecahedron" checks the rhombic dodecahedron periodic box
c     on the example/waterdodec.xyz water box, for van der Waals with
c     and without a neighbor list and for Ewald charges, plus the
c     wrapping of molecules into the box by "bounds"
c
c     the references hold half the energy, interaction count and
c     virial of an a x a x sqrt(2)a orthorhombic cell containing the
c     box twice, the second copy shifted by (a/2,a/2,a/sqrt(2)), and
c     the gradient rows of the first copy; that cell uses only the
c     orthogonal image code
c
c
      subroutine test_dodecahedron
      implicit none
c
c
      call test_dodecahedron_calc ('test_dodecahedron_vdw','vdw.key',
     &                             'vdw.1.txt','Van der Waals')
      call test_dodecahedron_calc ('test_dodecahedron_vdwnolist',
     &                             'vdwnolist.key','vdw.1.txt',
     &                             'Van der Waals')
      call test_dodecahedron_calc ('test_dodecahedron_ewald',
     &                             'ewald.key','ewald.1.txt',
     &                             'Charge-Charge')
      call test_dodecahedron_bounds
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_dodecahedron_calc  --  energy and derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine test_dodecahedron_calc (tname,key,ref,label)
      use action
      use atoms
      use energi
      use virial
      implicit none
      integer nat,refcnt
      real*8 energy,e,ref_e,ref_ei,refeng
      real*8 eps_e,eps_g,eps_v,eterm
      real*8, allocatable :: derivs(:,:)
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      logical skiptest
      character*(*) tname,key,ref,label
      character*240 rpath
c
c
      if (skiptest(tname,'amoeba'))  return
c
c     load the dodecahedral water box and the reference values
c
      call pushdir ('file/dodecahedron')
      call loadfix ('waterdodec',key)
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('dodecahedron',ref,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,label,refeng,refcnt)
      eps_e = 1.0d-4
      eps_g = 2.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  total energy and the single active term
c
      e = energy ()
      eterm = ev + ec
      call assert_real (esum,ref_e,eps_e,tname//' energy (v0)')
      call assert_real (eterm,refeng,eps_e,tname//' term (v0)')
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
      if (label .eq. 'Van der Waals') then
         call check_engcnt (rpath,label,ev,nev,eps_e,tname//' vdw (v3)')
      else
         call check_engcnt (rpath,label,ec,nec,eps_e,
     &                      tname//' charge (v3)')
      end if
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_dodecahedron_bounds  --  wrap molecules  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_dodecahedron_bounds" translates the whole box by a vector
c     that is not a lattice vector, wraps it with "bounds", and checks
c     that the energy is unchanged and every molecular center of mass
c     lies inside the dodecahedron
c
c
      subroutine test_dodecahedron_bounds
      use atomid
      use atoms
      use boxes
      use energi
      use molcul
      implicit none
      integer i,j,k,nout
      real*8 energy,e,ref_e,ref_ei
      real*8 xmid,ymid,zmid,wt
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      integer nat
      logical skiptest
      character*240 rpath
      character*(*) tname
      parameter (tname='test_dodecahedron_bounds')
c
c
      if (skiptest(tname,'amoeba'))  return
      call pushdir ('file/dodecahedron')
      call loadfix ('waterdodec','vdwnolist.key')
      allocate (refg(3,n))
      call refpath ('dodecahedron','vdw.1.txt',rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
c
c     translate, wrap and recompute the energy
c
      do i = 1, n
         x(i) = x(i) + 13.1d0
         y(i) = y(i) - 5.7d0
         z(i) = z(i) + 20.3d0
      end do
      call bounds
      e = energy ()
      call assert_real (esum,ref_e,1.0d-4,tname//' energy')
c
c     count molecular centers outside the dodecahedron
c
      nout = 0
      do i = 1, nmol
         xmid = 0.0d0
         ymid = 0.0d0
         zmid = 0.0d0
         do j = imol(1,i), imol(2,i)
            k = kmol(j)
            wt = mass(k)
            xmid = xmid + x(k)*wt
            ymid = ymid + y(k)*wt
            zmid = zmid + z(k)*wt
         end do
         xmid = xmid / molmass(i)
         ymid = ymid / molmass(i)
         zmid = zmid / molmass(i)
         if (abs(xmid/xbox)+abs(ymid/ybox)+abs(zmid/zbox2)
     &          .gt. 1.0d0+1.0d-8)  nout = nout + 1
      end do
      call assert_int (nout,0,tname//' centers inside')
      deallocate (refg)
      call popdir
      call final
      return
      end

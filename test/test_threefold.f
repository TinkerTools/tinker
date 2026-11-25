c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_threefold  --  3-Fold frame torque test  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_threefold" checks the force distribution of torques on a
c     3-Fold local frame, whose z-axis is the trisector of the three
c     bond vectors and whose x-axis comes from the x-axis atom alone
c
c     the fixture is triethylamine with an asymmetric 3-Fold N
c     multipole (in-plane dipole and off-diagonal quadrupole) and a
c     sodium ion 4 Ang from the N, 60 degrees off its lone-pair axis,
c     so the torque about the trisector is large; the older symmetric
c     distribution gave gradient errors near 1 kcal/mole/Ang here,
c     while the exact frame derivative matches central finite
c     differences of the energy; geometries are pyramidal N and N
c     0.08 Ang from the plane of its neighbors
c
c
      subroutine test_threefold
      implicit none
c
c
      call test_threefold_case ('pyramid')
      call test_threefold_case ('flat')
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_threefold_case  --  one 3-Fold geometry  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_threefold_case (tag)
      use atoms
      use virial
      implicit none
      integer i,j,k
      real*8 energy,e
      real*8 eplus,eminus
      real*8 eps,save
      real*8 vnum(3,3)
      real*8 pos(3)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: numg(:,:)
      logical skiptest
      character*(*) tag
      character*72 cname
c
c
      cname = 'test_threefold_'//trim(tag)
      if (skiptest(cname,'threefold amoeba',1))  return
      call pushdir ('file/threefold')
      call loadfix ('threefold_'//trim(tag),'threefold.key')
c
c     analytic gradient and virial
c
      allocate (derivs(3,n))
      allocate (numg(3,n))
      call gradient (e,derivs)
c
c     numerical gradient by central differences of the energy
c
      eps = 1.0d-5
      do i = 1, n
         do j = 1, 3
            if (j .eq. 1)  save = x(i)
            if (j .eq. 2)  save = y(i)
            if (j .eq. 3)  save = z(i)
            if (j .eq. 1)  x(i) = save + eps
            if (j .eq. 2)  y(i) = save + eps
            if (j .eq. 3)  z(i) = save + eps
            eplus = energy ()
            if (j .eq. 1)  x(i) = save - eps
            if (j .eq. 2)  y(i) = save - eps
            if (j .eq. 3)  z(i) = save - eps
            eminus = energy ()
            if (j .eq. 1)  x(i) = save
            if (j .eq. 2)  y(i) = save
            if (j .eq. 3)  z(i) = save
            numg(j,i) = (eplus-eminus) / (2.0d0*eps)
         end do
      end do
      call assert_grad (derivs,numg,n,1.0d-5,trim(cname)//
     &                  ' gradient')
c
c     virial of a free molecule from positions and gradient
c
      do j = 1, 3
         do k = 1, 3
            vnum(k,j) = 0.0d0
         end do
      end do
      do i = 1, n
         pos(1) = x(i)
         pos(2) = y(i)
         pos(3) = z(i)
         do j = 1, 3
            do k = 1, 3
               vnum(k,j) = vnum(k,j) + 0.5d0*(pos(k)*derivs(j,i)
     &                                       +pos(j)*derivs(k,i))
            end do
         end do
      end do
      call assert_grad (vir,vnum,3,1.0d-6,trim(cname)//' virial')
      deallocate (derivs)
      deallocate (numg)
      call popdir
      call final
      return
      end

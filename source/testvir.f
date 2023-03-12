c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2018  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program testvir  --  check analytical & numerical virial  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "testvir" computes the analytical internal virial and compares
c     it to a numerical virial derived from the finite difference
c     derivative of the energy with respect to lattice vectors
c
c
      program testvir
      use atoms
      use inform
      use iounit
      use virial
      implicit none
      integer i
      real*8 energy
      real*8, allocatable :: derivs(:,:)
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     set option control flags based desired analysis types
c
      debug = .false.
      allocate (derivs(3,n))
      call gradient (energy,derivs)
      deallocate (derivs)
c
c     print the components of the analytical internal virial
c
      write (iout,10)  (vir(1,i),vir(2,i),vir(3,i),i=1,3)
   10 format (/,' Analytical Virial Tensor :',9x,3f13.3,
     &           /,36x,3f13.3,/,36x,3f13.3)
c
c     get the numerical dE/dV value and a pressure estimate
c
      call ptest
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ptest  --  find pressure via finite differences  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ptest" determines the numerical virial tensor, and compares
c     analytical to numerical values for dE/dV and isotropic pressure
c
c     original version written by John D. Chodera, University of
c     California, Berkeley, December 2010
c
c     modified for off-diagonal numerical virial by Jay W. Ponder,
c     Saint Louis, August 2018
c
c
      subroutine ptest
      use atoms
      use bath
      use bound
      use boxes
      use iounit
      use math
      use units
      use virial
      implicit none
      integer i,j,k
      real*8 energy,epos,eneg
      real*8 eps,temp,lorig
      real*8 dedv_vir,dedv_num
      real*8 pres_vir,pres_num
      real*8 dedl(3,3)
      real*8 virn(3,3)
      real*8, allocatable :: xf(:)
      real*8, allocatable :: yf(:)
      real*8, allocatable :: zf(:)
c
c
c     set relative volume change for finite-differences
c
      if (.not. use_bounds)  return
      eps = 0.00003d0
c
c     set prism lattice type to the general triclinic case
c
      if (.not. nonprism) then
         orthogonal = .false.
         monoclinic = .false.
         triclinic = .true.
      end if
c
c     print out the lattice vectors as matrix rows
c
      write (iout,10)  (lvec(1,i),lvec(2,i),lvec(3,i),i=1,3)
   10 format (/,' Lattice Vectors (Lvec) :',11x,3f13.3,
     &           /,36x,3f13.3,/,36x,3f13.3)
c
c     perform dynamic allocation of some local arrays
c
      allocate (xf(n))
      allocate (yf(n))
      allocate (zf(n))
c
c     store the original fractional coordinate values
c
      do i = 1, n
         xf(i) = x(i)*recip(1,1) + y(i)*recip(2,1) + z(i)*recip(3,1)
         yf(i) = x(i)*recip(1,2) + y(i)*recip(2,2) + z(i)*recip(3,2)
         zf(i) = x(i)*recip(1,3) + y(i)*recip(2,3) + z(i)*recip(3,3)
      end do
c
c     get energy derivatives with respect to lattice vectors
c
      do i = 1, 3
         do j = i, 3
            dedl(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j = i, 3
            lorig = lvec(j,i)
            lvec(j,i) = lorig - eps
            call cellang (xf,yf,zf)
            eneg = energy ()
            lvec(j,i) = lorig + eps
            call cellang (xf,yf,zf)
            epos = energy ()
            lvec(j,i) = lorig
            call cellang (xf,yf,zf)
            dedl(j,i) = 0.5d0 * (epos-eneg) / eps
         end do
      end do
c
c     print out the partial derivatives of the energy
c
      write (iout,20)  (dedl(1,i),dedl(2,i),dedl(3,i),i=1,3)
   20 format (/,' dE/dLvec Derivatives :',13x,3f13.3,
     &           /,36x,3f13.3,/,36x,3f13.3)
c
c     perform deallocation of some local arrays
c
      deallocate (xf)
      deallocate (yf)
      deallocate (zf)
c
c     compute and print numerical virial tensor components
c
      do i = 1, 3
         do j = 1, i
            virn(j,i) = 0.0d0
            do k = 1, 3
               virn(j,i) = virn(j,i) + dedl(k,j)*lvec(k,i)
            end do
            virn(i,j) = virn(j,i)
         end do
      end do
      if (dodecadron) then
         write (iout,30)  (virn(1,1)+virn(2,2)+virn(3,3))/3.0d0
   30    format (/,' Numerical Mean Diagonal :',10x,f13.3)
      else if (octahedron) then
         write (iout,40)  virn(1,1),virn(2,2),virn(3,3)
   40    format (/,' Numerical Virial Diagonal :',8x,3f13.3)
      else
         write (iout,50)  (virn(1,i),virn(2,i),virn(3,i),i=1,3)
   50    format (/,' Numerical Virial Tensor :',10x,3f13.3,
     &              /,36x,3f13.3,/,36x,3f13.3)
      end if
c
c     find the analytical and numerical values of dE/dV
c
      dedv_vir = (vir(1,1)+vir(2,2)+vir(3,3)) / (3.0d0*volbox)
      dedv_num = (virn(1,1)+virn(2,2)+virn(3,3)) / (3.0d0*volbox)
c
c     get analytical and numerical isotropic pressure values
c
      temp = kelvin
      if (temp .eq. 0.0d0)  temp = 298.0d0
      pres_vir = prescon * (dble(n)*gasconst*temp/volbox-dedv_vir)
      pres_num = prescon * (dble(n)*gasconst*temp/volbox-dedv_num)
      write (iout,60)  nint(temp),pres_vir
   60 format (/,' Pressure (Analytical,',i4,' K) :',5x,f13.3,
     &           ' Atmospheres')
      write (iout,70)  nint(temp),pres_num
   70 format (' Pressure (Numerical,',i4,' K) :',6x,f13.3,
     &           ' Atmospheres')
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine cellang  --  lattice vectors to cell parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "cellang" computes atomic coordinates and unit cell parameters
c     from fractional coordinates and lattice vectors
c
c
      subroutine cellang (xf,yf,zf)
      use atoms
      use boxes
      use math
      implicit none
      integer i
      real*8 amag,bmag,cmag
      real*8 abdot,acdot,bcdot
      real*8 xf(*),yf(*),zf(*)
c
c
c     update coordinates via fractionals and lattice vectors
c
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
c
c     compute unit cell lengths and angles from lattice vectors
c
      amag = sqrt(lvec(1,1)**2+lvec(1,2)**2+lvec(1,3)**2)
      bmag = sqrt(lvec(2,1)**2+lvec(2,2)**2+lvec(2,3)**2)
      cmag = sqrt(lvec(3,1)**2+lvec(3,2)**2+lvec(3,3)**2)
      abdot = lvec(1,1)*lvec(2,1) + lvec(1,2)*lvec(2,2)
     &           + lvec(1,3)*lvec(2,3)
      acdot = lvec(1,1)*lvec(3,1) + lvec(1,2)*lvec(3,2)
     &           + lvec(1,3)*lvec(3,3)
      bcdot = lvec(2,1)*lvec(3,1) + lvec(2,2)*lvec(3,2)
     &           + lvec(2,3)*lvec(3,3)
      xbox = amag
      ybox = bmag
      zbox = cmag
      alpha = radian * acos(bcdot/(bmag*cmag))
      beta = radian * acos(acdot/(amag*cmag))
      gamma = radian * acos(abdot/(amag*bmag))
c
c     reset lattice parameters, box dimensions and volume
c
      call lattice
      return
      end

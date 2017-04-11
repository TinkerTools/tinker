c
c
c     ################################################################
c     ##  COPYRIGHT (C)  2009  by Chuanjie Wu & Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program testpres  --  compare virial to dU/dV pressure  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "testpres" computes the internal virial for a periodic system
c     and compares it to a numerical evaluation of dU/dV
c
c
      program testpres
      use sizes
      use atoms
      use boxes
      use bound
      use iounit
      use potent
      use units
      use usage
      use virial
      implicit none
      integer i,j,mode
      real*8 epot,epot1,epot2
      real*8 eps,factor
      real*8 dudv,virial
      real*8 stress(3,3)
      real*8, allocatable :: derivs(:,:)
      logical exist,query
      character*240 string
c
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     only appropriate if periodic boundaries are in use
c
      if (.not. use_bounds)  stop
c
c     get the method to use for scaling the coordinates
c
      mode = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         query = .false.
      end if
   10 continue
      if (query) then
         dowhile (mode.lt.1 .or. mode.gt.3)
            mode = 0
            write (iout,20)
   20       format (/,' Enter Coordinate Scaling Mode',
     &                 ' (1=Atom, 2=Molecule, 3=Group) :  ',$)
            read (input,30,err=40,end=40)  mode
   30       format (i10)
   40       continue
         end do
      end if
c
c     get the magnitude of the volume change to be used
c
      eps = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  eps
   50 continue
      if (eps .lt. 0.0d0) then
         write (iout,60)
   60    format (/,' Enter the Volume Perturbation in Cubic Angs',
     &              ' [0.1] :  ',$)
         read (input,70)  eps
   70    format (f20.0)
      end if
      if (eps .le. 0.0d0)  eps = 0.1d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     get the potential energy, forces and internal virial
c
      call gradient (epot,derivs)
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     calculate the full stress tensor for the periodic system
c
      factor = prescon / volbox
      virial = 0.0d0
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = -factor * vir(j,i)
         end do
      end do
c
c     set isotropic virial to the average of tensor diagonal
c
      virial = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
c
c     find dU/dV numerically by varying the periodic box size
c
      call vscale (epot1,mode,-0.5d0*eps)
      call vscale (epot2,mode,0.5d0*eps)
      dudv = -prescon * (epot2-epot1)/eps
      write (iout,80)
   80 format (/,' Comparison of Internal Virial and Numerical',
     &           ' dU/dV Term :')
      write (iout,90)  eps,virial,dudv,dudv-virial
   90 format (/,6x,'Volume Change',8x,'Virial',10x,'dU/dV',11x,'Delta',
     &        /,1x,f15.4,3x,3f15.4)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine vscale  --  energy due to small volume change  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "vscale" finds the potential energy for a small perturbation
c     of the periodic box volume
c
c
      subroutine vscale (epot,mode,eps)
      use sizes
      use atomid
      use atoms
      use boxes
      use group
      use molcul
      use usage
      implicit none
      integer i,j,k,mode
      integer start,stop
      real*8 epot,eps,vold
      real*8 energy,term
      real*8 scale,diff
      real*8 third,weigh
      real*8 xcm,ycm,zcm
      real*8 xmove,xboxold
      real*8 ymove,yboxold
      real*8 zmove,zboxold
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
c
c     save the original box size and atomic coordinates
c
      xboxold = xbox
      yboxold = ybox
      zboxold = zbox
      vold = volbox
      do i = 1, n
         xold(i) = x(i)
         yold(i) = y(i)
         zold(i) = z(i)
      end do
c
c     get scale factor that reflects the chosen volume change
c
      volbox = volbox + eps
      third = 1.0d0 / 3.0d0
      scale = (volbox/vold)**third
      if (monoclinic) then
         term = 1.0d0 / beta_sin
         scale = 1.0d0 + (scale-1.0d0)*term**third
      else if (triclinic) then
         term = 1.0d0 / (gamma_sin*gamma_term)
         scale = 1.0d0 + (scale-1.0d0)*term**third
      else if (octahedron) then
         term = 2.0d0
         scale = 1.0d0 + (scale-1.0d0)*term**third
      end if
c
c     set the new box dimensions and other lattice values
c
      xbox = xbox * scale
      ybox = ybox * scale
      zbox = zbox * scale
      call lattice
c
c     scale the coordinates by groups, molecules or atoms
c
      if (mode .eq. 1) then
         do i = 1, n
            if (use(i)) then
               x(i) = scale * x(i)
               y(i) = scale * y(i)
               z(i) = scale * z(i)
            end if
         end do
      else if (mode .eq. 2) then
         diff = scale - 1.0
         do i = 1, nmol
            xcm = 0.0d0
            ycm = 0.0d0
            zcm = 0.0d0
            start = imol(1,i)
            stop = imol(2,i)
            do j = start, stop
               k = kmol(j)
               weigh = mass(k)
               xcm = xcm + x(k)*weigh
               ycm = ycm + y(k)*weigh
               zcm = zcm + z(k)*weigh
            end do
            xmove = diff * xcm/molmass(i)
            ymove = diff * ycm/molmass(i)
            zmove = diff * zcm/molmass(i)
            do j = start, stop
               k = kmol(j)
               if (use(k)) then
                  x(k) = x(k) + xmove
                  y(k) = y(k) + ymove
                  z(k) = z(k) + zmove
               end if
            end do
         end do
      else if (mode .eq. 3) then
         diff = scale - 1.0
         do i = 1, ngrp
            xcm = 0.0d0
            ycm = 0.0d0
            zcm = 0.0d0
            start = igrp(1,i)
            stop = igrp(2,i)
            do j = start, stop
               k = kgrp(j)
               weigh = mass(k)
               xcm = xcm + x(k)*weigh
               ycm = ycm + y(k)*weigh
               zcm = zcm + z(k)*weigh
            end do
            xmove = diff * xcm/grpmass(i)
            ymove = diff * ycm/grpmass(i)
            zmove = diff * zcm/grpmass(i)
            do j = start, stop
               k = kgrp(j)
               x(k) = x(k) + xmove
               y(k) = y(k) + ymove
               z(k) = z(k) + zmove
            end do
         end do
      end if
c
c     find the potential energy of the system after the move
c
      epot = energy ()
c
c     reverse the move, restore old box size and coordinates
c
      xbox = xboxold
      ybox = yboxold
      zbox = zboxold
      call lattice
      volbox = vold
      do i = 1, n
         x(i) = xold(i)
         y(i) = yold(i)
         z(i) = zold(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      return
      end

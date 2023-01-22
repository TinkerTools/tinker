c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program minrigid  --  low store BFGS rigid body optimizer  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "minrigid" performs an energy minimization of rigid body atom
c     groups using a low storage BFGS nonlinear optimization
c
c
      program minrigid
      use files
      use group
      use inform
      use iounit
      use output
      use rigid
      implicit none
      integer i,j
      integer imin,nvar
      integer freeunit
      real*8 minimum,minrigid1
      real*8 grdmin,grms,gnorm
      real*8, allocatable :: xx(:)
      real*8, allocatable :: derivs(:,:)
      logical exist
      character*240 minfile
      character*240 string
      external minrigid1
      external optsave
c
c
c     set up the molecular mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     set up the use of rigid body coordinate system
c
      coordtype = 'RIGIDBODY'
      use_rigid = .true.
      call orient
c
c     perform the setup functions needed for optimization
c
      call optinit
c
c     get termination criterion as RMS rigid body gradient
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  grdmin
   10 continue
      if (grdmin .le. 0.0d0) then
         write (iout,20)
   20    format (/,' Enter RMS Gradient per Rigid Body Criterion',
     &              ' [0.01] :  ',$)
         read (input,30)  grdmin
   30    format (f20.0)
      end if
      if (grdmin .eq. 0.0d0)  grdmin = 0.01d0
c
c     write out a copy of coordinates for later update
c
      imin = freeunit ()
      minfile = basename(1:leng)//'.xyz'
      call version (minfile,'new')
      open (unit=imin,file=minfile,status='new')
      call prtxyz (imin)
      close (unit=imin)
      outfile = minfile
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(6*ngrp))
c
c     convert rigid body coordinates to optimization parameters
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            xx(nvar) = rbc(j,i)
         end do
      end do
c
c     make the call to the optimization routine
c
      call lbfgs (nvar,xx,minimum,grdmin,minrigid1,optsave)
c
c     convert optimization parameters to rigid body coordinates
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            rbc(j,i) = xx(nvar)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(6,ngrp))
c
c     compute the final function and RMS gradient values
c
      call gradrgd (minimum,derivs)
      gnorm = 0.0d0
      do i = 1, ngrp
         do j = 1, 6
            gnorm = gnorm + derivs(j,i)**2
         end do
      end do
      gnorm = sqrt(gnorm)
      grms = gnorm / sqrt(dble(ngrp))
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     write out the final function and gradient values
c
      if (digits .ge. 8) then
         if (grms .gt. 1.0d-8) then
            write (iout,40)  minimum,grms,gnorm
   40       format (/,' Final Function Value :',2x,f20.8,
     &              /,' Final RMS Gradient :',4x,f20.8,
     &              /,' Final Gradient Norm :',3x,f20.8)
         else
            write (iout,50)  minimum,grms,gnorm
   50       format (/,' Final Function Value :',2x,f20.8,
     &              /,' Final RMS Gradient :',4x,d20.8,
     &              /,' Final Gradient Norm :',3x,d20.8)
         end if
      else if (digits .ge. 6) then
         if (grms .gt. 1.0d-6) then
            write (iout,60)  minimum,grms,gnorm
   60       format (/,' Final Function Value :',2x,f18.6,
     &              /,' Final RMS Gradient :',4x,f18.6,
     &              /,' Final Gradient Norm :',3x,f18.6)
         else
            write (iout,70)  minimum,grms,gnorm
   70       format (/,' Final Function Value :',2x,f18.6,
     &              /,' Final RMS Gradient :',4x,d18.6,
     &              /,' Final Gradient Norm :',3x,d18.6)
         end if
      else
         if (grms .gt. 1.0d-4) then
            write (iout,80)  minimum,grms,gnorm
   80       format (/,' Final Function Value :',2x,f16.4,
     &              /,' Final RMS Gradient :',4x,f16.4,
     &              /,' Final Gradient Norm :',3x,f16.4)
         else
            write (iout,90)  minimum,grms,gnorm
   90       format (/,' Final Function Value :',2x,f16.4,
     &              /,' Final RMS Gradient :',4x,d16.4,
     &              /,' Final Gradient Norm :',3x,d16.4)
         end if
      end if
c
c     write the final coordinates into a file
c
      imin = freeunit ()
      open (unit=imin,file=minfile,status='old')
      rewind (unit=imin)
      call prtxyz (imin)
      close (unit=imin)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function minrigid1  --  energy and gradient for minrigid  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "minrigid1" is a service routine that computes the energy
c     and gradient for a low storage BFGS nonlinear optimization
c     of rigid bodies
c
c
      function minrigid1 (xx,g)
      use group
      use math
      use rigid
      implicit none
      integer i,j,nvar
      real*8 minrigid1,e
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     convert optimization parameters to rigid body coordinates
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            rbc(j,i) = xx(nvar)
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(6,ngrp))
c
c     compute and store the energy and gradient
c
      call rigidxyz
      call gradrgd (e,derivs)
      minrigid1 = e
c
c     convert rigid body gradient to optimization parameters
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            g(nvar) = derivs(j,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end

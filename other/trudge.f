c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2006  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program trudge  --  Nelder-Mead derivative-free optimizer  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "trudge" performs energy minimization in Cartesian coordinate
c     space using a derivative-free Nelder-Mead simplex optimization
c
c
      program trudge
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'scales.i'
      include 'usage.i'
      integer i,j,imin,nvar
      integer next,freeunit
      real*8 minimiz1,minimum
      real*8 grdmin,gnorm,grms
      real*8 xx(maxvar)
      real*8 derivs(3,maxatm)
      logical exist
      character*20 keyword
      character*240 minfile
      character*240 record
      character*240 string
      external minimiz1
      external optsave
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     search the keywords for output frequency parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (string,*,err=10,end=10)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (string,*,err=10,end=10)  iwrite
         end if
   10    continue
      end do
c
c     get termination criterion as RMS gradient per atom
c
      ftol = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=20,end=20)  ftol
   20 continue
      if (ftol .le. 0.0d0) then
         write (iout,30)
   30    format (/,' Enter Energy Tolerance Convergence Criterion',
     &              ' [0.01] :  ',$)
         read (input,40)  ftol
   40    format (f20.0)
      end if
      if (ftol .le. 0.0d0)  ftol = 0.01d0
c
c     write out a copy of coordinates for later update
c
      imin = freeunit ()
      minfile = filename(1:leng)//'.xyz'
      call version (minfile,'new')
      open (unit=imin,file=minfile,status='new')
      call prtxyz (imin)
      close (unit=imin)
      outfile = minfile
c
c     set active atom coordinates as optimization variables
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            xx(1,nvar) = x(i)
            nvar = nvar + 1
            xx(1,nvar) = y(i)
            nvar = nvar + 1
            xx(1,nvar) = z(i)
         end if
      end do
      y(1) = trudge0 (xx(1,1))
c
c     build initial simplex by incrementing each variable
c
      do k = 2, nvar+1
         do i = 1, nvar
            xx(j,
      do i = 1, nvar
         xold = xx(i)
         xx(i) = xx(i) + eps
         y(i+1) = trudge
      end do
c
c     make the call to the optimization routine
c
      call simplex (nvar,xx,y,ftol,trudge0,iter)
c
c     unpack the final coordinates for active atoms
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
c
c     compute the final function and RMS gradient values
c
      minimum = energy ()
c
c     write out the final function and gradient values
c
      if (digits .ge. 8) then
         if (grms .gt. 1.0d-8) then
            write (iout,50)  minimum
   50       format (/,' Final Function Value :',2x,f20.8)
         else
            write (iout,60)  minimum
   60       format (/,' Final Function Value :',2x,f20.8)
         end if
      else if (digits .ge. 6) then
         if (grms .gt. 1.0d-6) then
            write (iout,70)  minimum
   70       format (/,' Final Function Value :',2x,f18.6)
         else
            write (iout,80)  minimum
   80       format (/,' Final Function Value :',2x,f18.6)
         end if
      else
         if (grms .gt. 1.0d-4) then
            write (iout,90)  minimum
   90       format (/,' Final Function Value :',2x,f16.4)
         else
            write (iout,100)  minimum
  100       format (/,' Final Function Value :',2x,f16.4)
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
c     ##############################################################
c     ##                                                          ##
c     ##  function trudge0  --  energy for Nelder-Mead optimizer  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "trudge0" is a service routine that computes the energy for
c     a derivative-free Nelder-Mead simplex optimization in Cartesian
c     coordinate space
c
c
      function trudge0 (xx)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'usage.i'
      integer i,nvar
      real*8 trudge0,e
      real*8 xx(maxvar)
c
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar)
         end if
      end do
c
c     compute and store the energy and gradient
c
      e = energy ()
      trudge0 = e
      return
      end

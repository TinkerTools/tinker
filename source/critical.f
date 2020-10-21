c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2020  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program critical  --  stationary point by least squares  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "critical" finds a stationary point for a molecular system via
c     least squares minimization of the atomic gradient components
c
c
      program critical
      use sizes
      use atoms
      use files
      use inform
      use iounit
      use keys
      use minima
      use usage
      implicit none
      integer i,j,k,next
      integer nvar,nrsd
      integer imin,freeunit
      real*8 epot,grdmin
      real*8 gnorm,grms
      real*8, allocatable :: xx(:)
      real*8, allocatable :: xlo(:)
      real*8, allocatable :: xhi(:)
      real*8, allocatable :: rsd(:)
      real*8, allocatable :: grd(:)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: fjac(:,:)
      logical exist
      character*20 keyword
      character*240 minfile
      character*240 record
      character*240 string
      external critical1
      external critsave
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
         string = record(next:120)
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
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=20,end=20)  grdmin
   20 continue
      if (grdmin .le. 0.0d0) then
         write (iout,30)
   30    format (/,' Enter Residual Gradient Convergence Criterion',
     &              ' [0.01] :  ',$)
         read (input,40)  grdmin
   40    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.01d0
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
c     perform dynamic allocation of some local arrays
c
      allocate (xx(3*n))
      allocate (derivs(3,n))
c
c     set active atom coordinates as optimization variables
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            xx(nvar) = x(i)
            nvar = nvar + 1
            xx(nvar) = y(i)
            nvar = nvar + 1
            xx(nvar) = z(i)
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xlo(3*n))
      allocate (xhi(3*n))
      allocate (rsd(3*n))
      allocate (grd(3*n))
      allocate (fjac(3*n,3*n))
c
c     make the call to the least squares optimization routine
c
      maxiter = 10000
      nrsd = nvar
      do i = 1, nvar
         xlo(i) = -1000000.0d0
         xhi(i) = 1000000.0d0
      end do
      call square (nvar,nrsd,xlo,xhi,xx,rsd,grd,fjac,
     &                grdmin,critical1,critsave)
c
c     perform deallocation of some local arrays
c
      deallocate (xlo)
      deallocate (xhi)
      deallocate (rsd)
      deallocate (grd)
      deallocate (fjac)
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
      call gradient (epot,derivs)
      gnorm = 0.0d0
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            gnorm = gnorm + derivs(j,k)**2
         end do
      end do
      gnorm = sqrt(gnorm)
      grms = gnorm / sqrt(dble(nvar/3))
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      deallocate (derivs)
c
c     write out the final function and gradient values
c
      if (digits .ge. 8) then
         if (grms .gt. 1.0d-8) then
            write (iout,50)  epot,grms,gnorm
   50       format (/,' Final Function Value :',2x,f20.8,
     &              /,' Final RMS Gradient :',4x,f20.8,
     &              /,' Final Gradient Norm :',3x,f20.8)
         else
            write (iout,60)  epot,grms,gnorm
   60       format (/,' Final Function Value :',2x,f20.8,
     &              /,' Final RMS Gradient :',4x,d20.8,
     &              /,' Final Gradient Norm :',3x,d20.8)
         end if
      else if (digits .ge. 6) then
         if (grms .gt. 1.0d-6) then
            write (iout,70)  epot,grms,gnorm
   70       format (/,' Final Function Value :',2x,f18.6,
     &              /,' Final RMS Gradient :',4x,f18.6,
     &              /,' Final Gradient Norm :',3x,f18.6)
         else
            write (iout,80)  epot,grms,gnorm
   80       format (/,' Final Function Value :',2x,f18.6,
     &              /,' Final RMS Gradient :',4x,d18.6,
     &              /,' Final Gradient Norm :',3x,d18.6)
         end if
      else
         if (grms .gt. 1.0d-4) then
            write (iout,90)  epot,grms,gnorm
   90       format (/,' Final Function Value :',2x,f16.4,
     &              /,' Final RMS Gradient :',4x,f16.4,
     &              /,' Final Gradient Norm :',3x,f16.4)
         else
            write (iout,100)  epot,grms,gnorm
  100       format (/,' Final Function Value :',2x,f16.4,
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine critical1  --  least squares gradient residual  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "trudge1" is a service routine to compute gradient components
c     for a least squares minimization to a stationary point
c
c
      subroutine critical1 (nvar,nrsd,xx,rsd)
      use sizes
      use atoms
      use usage
      implicit none
      integer i
      integer nvar
      integer nrsd
      real*8 epot 
      real*8 xx(*)
      real*8 rsd(*)
      real*8, allocatable :: derivs(:,:)
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
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     compute energy and gradient for the current structure
c
      call gradient (epot,derivs)
c
c     store the gradient components as the residual vector
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            rsd(nvar) = derivs(1,i)
            nvar = nvar + 1
            rsd(nvar) = derivs(2,i)
            nvar = nvar + 1
            rsd(nvar) = derivs(3,i)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine critsave  --  critical point output routine  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine critsave (niter,nrsd,xx,gs,rsd)
      use files
      use iounit
      use output
      implicit none
      integer niter,nrsd
      integer iopt,iend
      integer lext
      integer freeunit
      real*8 xx(*)
      real*8 gs(*)
      real*8 rsd(*)
      logical exist
      character*7 ext
      character*240 optfile
      character*240 endfile
c
c
c     get name of archive or intermediate coordinates file
c
      iopt = freeunit ()
      if (cyclesave) then
         if (archive) then
            optfile = filename(1:leng)
            call suffix (optfile,'arc','old')
            inquire (file=optfile,exist=exist)
            if (exist) then
               call openend (iopt,optfile)
            else
               open (unit=iopt,file=optfile,status='new')
            end if
         else
            lext = 3
            call numeral (niter,ext,lext)
            optfile = filename(1:leng)//'.'//ext(1:lext)
            call version (optfile,'new')
            open (unit=iopt,file=optfile,status='new')
         end if
      else
         optfile = outfile
         call version (optfile,'old')
         open (unit=iopt,file=optfile,status='old')
         rewind (unit=iopt)
      end if
c
c     update intermediate file with desired coordinate type
c
      call prtxyz (iopt)
      close (unit=iopt)
c
c     test for requested termination of the optimization
c
      endfile = 'tinker.end'
      inquire (file=endfile,exist=exist)
      if (.not. exist) then
         endfile = filename(1:leng)//'.end'
         inquire (file=endfile,exist=exist)
         if (exist) then
            iend = freeunit ()
            open (unit=iend,file=endfile,status='old')
            close (unit=iend,status='delete')
         end if
      end if
      if (exist) then
         write (iout,10)
   10    format (/,' CRITSAVE  --  Optimization Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
      return
      end

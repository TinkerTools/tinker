c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1998 by Rohit Pappu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program scan  --  maps minima on potential energy surface  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "scan" attempts to find all the local minima on a potential
c     energy surface via an iterative series of local searches along
c     normal mode directions
c
c     literature reference:
c
c     I. Kolossvary and W. C. Guida, "Low-Mode Conformational Search
c     Elucidated: Application to C39H80 and Flexible Docking of
c     9-Deazaguanine Inhibitors into PNP, Journal of Computational
c     Chemistry, 20, 1671-1684 (1999)
c
c
      program scan
      use files
      use inform
      use iounit
      use omega
      use output
      implicit none
      integer maxmap
      parameter (maxmap=100000)
      integer i,ixyz
      integer lext,freeunit
      integer nmap,niter
      integer nvec,neigen
      real*8 minimum,grdmin,range
      real*8 emap(maxmap)
      logical exist
      character*7 ext
      character*240 xyzfile
      character*240 string
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     initialize the number of minima and coordinate type
c
      nmap = 0
      coordtype = 'CARTESIAN'
c
c     get the rotatable bonds for torsional local search
c
      call makeint (0)
      call initrot
      call active
c
c     get the number of eigenvectors to use for the local search
c
      neigen = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  neigen
   10 continue
      nvec = min(nomega,5)
      if (neigen .le. 0) then
         write (iout,20)  nvec
   20    format(/,' Enter the Number of Eigenvectors for Local',
     &             ' Search [',i1,'] :  ',$)
         read (input,30)  neigen
   30    format (i10)
         if (neigen .le. 0)  neigen = nvec
      end if
      neigen = min(neigen,nvec)
c
c     get the energy threshold criterion for map membership
c
      range = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  range
   40 continue
      if (range .le. 0.0d0) then
         write (iout,50)
   50    format (/,' Enter the Energy Threshold for Local Minima',
     &              ' [100.0] :  ',$)
         read (input,60)  range
   60    format (f20.0)
      end if
      if (range .le. 0.0d0)  range = 100.0d0
c
c     get the termination criterion as RMS gradient per atom
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  grdmin
   70 continue
      if (grdmin .le. 0.0d0) then
         write (iout,80)
   80    format (/,' Enter RMS Gradient per Atom Criterion',
     &              ' [0.0001] :  ',$)
         read (input,90)  grdmin
   90    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.0001d0
c
c     set the energy output precision via convergence criterion
c
      if (grdmin .le. 0.000001d0)  digits = 6
      if (grdmin .le. 0.00000001d0)  digits = 8
c
c     create and open an output file if using archive mode
c
      if (archive) then
         ixyz = freeunit ()
         xyzfile = filename(1:leng)
         call suffix (xyzfile,'arc','new')
         open (unit=ixyz,file=xyzfile,status='new')
         close (unit=ixyz)
      end if
c
c     find the first map point from the input structure
c
      write (iout,100)
  100 format (/,' Generating Seed Point for Potential Energy',
     &           ' Surface Scan',/)
      call localmin (minimum,grdmin)
      call mapcheck (nmap,emap,range,minimum,grdmin)
c
c     use normal mode local search to explore adjacent minima
c
      niter = 0
      do while (niter .lt. nmap)
         niter = niter + 1
         write (iout,110)  niter
  110    format (/,' Normal Mode Local Search',7x,'Minimum',i7,/)
         ixyz = freeunit ()
         if (archive) then
            xyzfile = filename(1:leng)
            call suffix (xyzfile,'arc','old')
            open (unit=ixyz,file=xyzfile,status='old')
            do i = 1, niter-1
               call readxyz (ixyz)
            end do
         else
            lext = 3
            call numeral (niter,ext,lext)
            xyzfile = filename(1:leng)//'.'//ext(1:lext)
            call version (xyzfile,'old')
            open (unit=ixyz,file=xyzfile,status='old')
         end if
         call readxyz (ixyz)
         close (unit=ixyz)
         call modesrch (nmap,emap,range,neigen,grdmin)
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mapcheck  --  addition to local minimum list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mapcheck" checks the current minimum energy structure
c     for possible addition to the master list of local minima
c
c
      subroutine mapcheck (nmap,emap,range,minimum,grdmin)
      use files
      use inform
      use iounit
      use output
      implicit none
      integer i,ixyz,lext
      integer nmap,freeunit
      real*8 minimum,grdmin
      real*8 delta,eps,range
      real*8 emap(*)
      logical unique,exist
      character*7 ext
      character*240 xyzfile
c
c
c     check to see if the current minimum was previously found
c
      eps = grdmin
      unique = .true.
      do i = 1, nmap
         delta = minimum - emap(i)
         if (abs(delta) .lt. eps)  unique = .false.
         if (delta .gt. range)  unique = .false.
      end do
c
c     add minimum to master list if it was not previously known
c
      if (unique) then
         nmap = nmap + 1
         emap(nmap) = minimum
         if (digits .ge. 8) then
            write (iout,10)  nmap,minimum
   10       format (/,4x,'Potential Surface Map',7x,'Minimum',
     &                 i7,6x,f20.8,/)
         else if (digits .ge. 6) then
            write (iout,20)  nmap,minimum
   20       format (/,4x,'Potential Surface Map',7x,'Minimum',
     &                 i7,6x,f18.6,/)
         else
            write (iout,30)  nmap,minimum
   30       format (/,4x,'Potential Surface Map',7x,'Minimum',
     &                 i7,6x,f16.4,/)
         end if
c
c     write the coordinates of the new minimum to a file
c
         ixyz = freeunit ()
         if (archive) then
            xyzfile = filename(1:leng)
            call suffix (xyzfile,'arc','old')
            inquire (file=xyzfile,exist=exist)
            if (exist) then
               call openend (ixyz,xyzfile)
            else
               open (unit=ixyz,file=xyzfile,status='new')
            end if
         else
            lext = 3
            call numeral (nmap,ext,lext)
            xyzfile = filename(1:leng)//'.'//ext(1:lext)
            call version (xyzfile,'new')
            open (unit=ixyz,file=xyzfile,status='new')
         end if
         call prtxyz (ixyz)
         close (unit=ixyz)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function scan1  --  energy and gradient values for scan  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "scan1" is a service routine that computes the energy and
c     gradient during exploration of a potential energy surface
c     via iterative local search
c
c
      function scan1 (xx,g)
      use atoms
      implicit none
      integer i,nvar
      real*8 scan1,e
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     convert optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     compute and store the energy and gradient
c
      call gradient (e,derivs)
      scan1 = e
c
c     convert gradient components to optimization parameters
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         g(nvar) = derivs(1,i)
         nvar = nvar + 1
         g(nvar) = derivs(2,i)
         nvar = nvar + 1
         g(nvar) = derivs(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine scan2  --  Hessian matrix values for scan  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "scan2" is a service routine that computes the sparse matrix
c     Hessian elements during exploration of a potential energy
c     surface via iterative local search
c
c
      subroutine scan2 (mode,xx,h,hinit,hstop,hindex,hdiag)
      use atoms
      implicit none
      integer i,nvar
      integer hinit(*)
      integer hstop(*)
      integer hindex(*)
      real*8 xx(*)
      real*8 hdiag(*)
      real*8 h(*)
      character*4 mode
c
c
c     translate optimization parameters to atomic coordinates
c
      if (mode .eq. 'NONE')  return
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
c
c     compute and store the Hessian elements
c
      call hessian (h,hinit,hstop,hindex,hdiag)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine modesrch  --  normal mode local search  ##
c     ##                                                     ##
c     #########################################################
c
c
      subroutine modesrch (nmap,emap,range,neigen,grdmin)
      use iounit
      use omega
      implicit none
      integer i,k,nsearch
      integer nmap,neigen
      real*8 minimum,grdmin,range
      real*8 emap(*)
      real*8, allocatable :: step(:)
      real*8, allocatable :: eigen(:)
      real*8, allocatable :: vects(:,:)
c
c
c     store the current coordinates as the reference set
c
      call makeref (1)
c
c     perform dynamic allocation of some local arrays
c
      allocate (step(nomega))
      allocate (eigen(nomega))
      allocate (vects(nomega,nomega))
c
c     convert to internal coordinates and find torsional modes
c
      call makeint (0)
      call eigenrot (eigen,vects)
c
c     search both directions along each torsional eigenvector
c
      nsearch = 0
      do i = 1, neigen
         do k = 1, nomega
            step(k) = vects(k,nomega-i+1)
         end do
         nsearch = nsearch + 1
         call climber (nsearch,minimum,step,grdmin)
         call mapcheck (nmap,emap,range,minimum,grdmin)
         do k = 1, nomega
            step(k) = -vects(k,nomega-i+1)
         end do
         nsearch = nsearch + 1
         call climber (nsearch,minimum,step,grdmin)
         call mapcheck (nmap,emap,range,minimum,grdmin)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (step)
      deallocate (eigen)
      deallocate (vects)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eigenrot  --  torsional Hessian eigenvectors  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine eigenrot (eigen,vects)
      use atoms
      use omega
      implicit none
      integer i,j,ihess
      real*8 vnorm
      real*8 eigen(*)
      real*8, allocatable :: matrix(:)
      real*8 vects(nomega,*)
      real*8, allocatable :: hrot(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (matrix(nomega*(nomega+1)/2))
      allocate (hrot(nomega,nomega))
c
c     compute the Hessian in torsional space
c
      call hessrot ('FULL',hrot)
c
c     place Hessian elements into triangular form
c
      ihess = 0
      do i = 1, nomega
         do j = i, nomega
            ihess = ihess + 1
            matrix(ihess) = hrot(i,j)
         end do
      end do
c
c     diagonalize the Hessian to obtain eigenvalues
c
      call diagq (nomega,nomega,matrix,eigen,vects)
c
c     perform deallocation of some local arrays
c
      deallocate (matrix)
      deallocate (hrot)
c
c     normalize the torsional Hessian eigenvectors
c
      do i = 1, nomega
         vnorm = 0.0d0
         do j = 1, nomega
            vnorm = vnorm + vects(j,i)**2
         end do
         vnorm = sqrt(vnorm)
         do j = 1, nomega
            vects(j,i) = vects(j,i) / vnorm
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine climber  --  explore single search direction  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine climber (nsearch,minimum,step,grdmin)
      use inform
      use iounit
      use math
      use omega
      use potent
      use zcoord
      implicit none
      integer maxstep
      parameter (maxstep=500)
      integer i,kstep
      integer nstep,nsearch
      real*8 minimum,grdmin
      real*8 big,energy,size
      real*8 estep(0:maxstep)
      real*8 step(*)
      logical done
      logical oldpolar
c
c
c     convert current reference coordinates to a Z-matrix
c
      call getref (1)
      call makeint (0)
c
c     set the maximum number of steps and the step size
c
      done = .false.
      big = 100000.0d0
      minimum = big
      kstep = 0
      nstep = 65
      size = 0.1d0 * radian
      do i = 1, nomega
         step(i) = size * step(i)
      end do
c
c     scan the search direction for a minimization candidate
c
      do while (.not. done)
         if (kstep .ne. 0) then
            do i = 1, nomega
               ztors(zline(i)) = ztors(zline(i)) + step(i)
            end do
         end if
         call makexyz
         oldpolar = use_polar
         use_polar = .false.
         estep(kstep) = energy ()
         use_polar = oldpolar
         if (kstep .ge. 2) then
            if (estep(kstep) .lt. estep(kstep-2) .and.
     &          estep(kstep-1) .lt. estep(kstep-2)) then
               done = .true.
               do i = 1, nomega
                  ztors(zline(i)) = ztors(zline(i)) - step(i)
               end do
               call makexyz
               call localmin (minimum,grdmin)
               if (minimum .le. -big) then
                  minimum = big
                  write (iout,10)  nsearch
   10             format (4x,'Search Direction',i4,38x,'<<<<<<')
               else if (minimum .ge. big) then
                  minimum = big
                  write (iout,20)  nsearch
   20             format (4x,'Search Direction',i4,38x,'>>>>>>')
               else
                  if (digits .ge. 8) then
                     write (iout,30)  nsearch,kstep-1,minimum
   30                format (4x,'Search Direction',i4,11x,'Step',
     &                          i7,6x,f20.8)
                  else if (digits .ge. 6) then
                     write (iout,40)  nsearch,kstep-1,minimum
   40                format (4x,'Search Direction',i4,11x,'Step',
     &                          i7,6x,f18.6)
                  else
                     write (iout,50)  nsearch,kstep-1,minimum
   50                format (4x,'Search Direction',i4,11x,'Step',
     &                          i7,6x,f16.4)
                  end if
               end if
            end if
         end if
         if (kstep.ge.nstep .and. .not.done) then
            done = .true.
            write (iout,60)  nsearch
   60       format (4x,'Search Direction',i4,38x,'------')
         end if
         kstep = kstep + 1
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine localmin  --  optimize local search candidate  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "localmin" is used during normal mode local search to
c     perform a Cartesian coordinate energy minimization
c
c
      subroutine localmin (minimum,grdmin)
      use atoms
      use inform
      use minima
      use output
      use potent
      use scales
      implicit none
      integer i,j,nvar
      real*8 minimum,scan1
      real*8 grdmin,oldgrd
      real*8 gnorm,grms,big
      real*8, allocatable :: xx(:)
      real*8, allocatable :: derivs(:,:)
      logical oldverb,oldpolar
      character*6 mode,method
      external scan1,scan2
      external optsave
c
c
c     initialize optimization output and maximum energy
c
      iwrite = 0
      iprint = 0
      big = 100000.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(3*n))
c
c     convert atomic coordinates to optimization parameters
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         xx(nvar) = x(i)
         nvar = nvar + 1
         xx(nvar) = y(i)
         nvar = nvar + 1
         xx(nvar) = z(i)
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (.not. set_scale) then
         if (.not. allocated(scale))  allocate (scale(nvar))
c
c     set scaling parameters to unity due to mixed optimization
c
         set_scale = .true.
         do i = 1, nvar
            scale(i) = 1.0d0
         end do
      end if
c
c     adjust polarization and set initial optimization values
c
      oldverb = verbose
      oldpolar = use_polar
      oldgrd = grdmin
      verbose = .false.
      use_polar = .false.
      grdmin = 3.0
c
c     initial optimizaton to get close to approximate minimum
c
      call lbfgs (nvar,xx,minimum,grdmin,scan1,optsave)
c
c     secondary optimization to reach the exact local minimum
c
      use_polar = oldpolar
      grdmin = oldgrd
      mode = 'AUTO'
      method = 'AUTO'
      call tncg (mode,method,nvar,xx,minimum,
     &           grdmin,scan1,scan2,optsave)
      verbose = oldverb
c
c     convert optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     independently check the gradient convergence criterion
c
      call gradient (minimum,derivs)
      gnorm = 0.0d0
      do i = 1, n
         do j = 1, 3
            gnorm = gnorm + derivs(j,i)**2
         end do
      end do
      gnorm = sqrt(gnorm)
      grms = gnorm / sqrt(dble(n))
      if (grms .gt. grdmin)  minimum = big
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program testpol  --  check convergence of induced dipoles  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "testpol" computes the induced dipole moments for direct
c     polarization and increasing numbers of SCF iterations in
c     order to monitor the degree and rate of convergence
c
c
      program testpol
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'inform.i'
      include 'iounit.i'
      include 'polar.i'
      include 'polpot.i'
      include 'rigid.i'
      include 'units.i'
      integer i,j,k
      integer maxiter
      real*8 sum
      real*8, allocatable :: rms(:)
      real*8, allocatable :: drms(:)
      real*8, allocatable :: uexact(:,:)
      real*8, allocatable :: ustore(:,:,:)
c
c
c     get the coordinates and required force field parameters
c
      call initial
      call getxyz
      call mechanic
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c
c     perform dynamic allocation of some local arrays
c
      maxiter = 100
      allocate (rms(0:maxiter))
      allocate (drms(maxiter))
      allocate (uexact(3,n))
      allocate (ustore(3,n,0:maxiter))
c
c     set tolerance and rotate multipoles to global frame
c
      debug = .false.
      poleps = 0.0000000001d0
      call chkpole
      call rotpole
c
c     get induced dipoles for noniterative direct polarization
c
      poltyp = 'DIRECT'
      call induce
      do i = 1, n
         do j = 1, 3
            ustore(j,i,0) = debye * uind(j,i)
         end do
      end do
c
c     compute induced dipoles for increasing iteration counts
c
      poltyp = 'MUTUAL'
      do k = 1, maxiter
         politer = k
         call induce
         do i = 1, n
            do j = 1, 3
               ustore(j,i,k) = debye * uind(j,i)
            end do
         end do
      end do
c
c     find the RMS difference per atom between iterations
c
      do k = 1, maxiter
         sum = 0.0d0
         do i = 1, n
            do j = 1, 3
               sum = sum + (ustore(j,i,k)-ustore(j,i,k-1))**2
            end do
         end do
         drms(k) = sqrt(sum/dble(n)) 
      end do
c
c     find the RMS of each iteration from converged dipoles
c
      do i = 1, n
         do j = 1, 3
            uexact(j,i) = ustore(j,i,maxiter)
         end do
      end do
      do k = 0, maxiter
         sum = 0.0d0
         do i = 1, n
            do j = 1, 3
               sum = sum + (ustore(j,i,k)-uexact(j,i))**2
            end do
         end do
         rms(k) = sqrt(sum/dble(n)) 
      end do
c
c     print the RMS between iterations and to converged dipoles
c
      write (iout,10)
   10 format (/,' Convergence of Induced Dipole Moments :',
     &        //,2x,'Iteration',10x,'RMS Change',11x,'RMS vs Final')
      write (iout,20)  0,rms(0)
   20 format (/,i8,16x,'----',6x,f20.10)
      do k = 1, maxiter
         write (iout,30)  k,drms(k),rms(k)
   30    format (i8,3x,f20.10,3x,f20.10)
         if (rms(k) .lt. 0.5d0*poleps)  goto 40
      end do
   40 continue
c
c     perform deallocation of some local arrays
c
      deallocate (rms)
      deallocate (drms)
      deallocate (uexact)
      deallocate (ustore)
      end

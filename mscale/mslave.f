c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program mslave  --  client for CHARMM MSCALE communication  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "mslave" is a simple Tinker client/slave program for use with
c     the CHARMM MSCALE master
c
c     note use of the "mpif.h" header file is deprecated in MPI and
c     may be removed in a future version; also it may be necessary
c     to use the C-style pre-processor directive #include <mpif.h>
c     and the -cpp command line option, instead of include 'mpif.h' 
c
c
      program mslave
      implicit none
      include 'mpif.h'
      use sizes
      use atoms
      use energi
      use hescut
      use iounit
      use output
      integer length
      parameter (length=128)
      integer i,j,k,m
      integer n3,n6
      integer ixyz
      integer nhess
      integer freeunit
      integer meroot
      integer mpidp,mpilog
      integer ncontrol
      integer irank,ierr
      integer parent
      integer, allocatable :: hindex(:)
      integer, allocatable :: hinit(:::)
      integer, allocatable :: hstop(:,:)
      real*8 e
      real*8 cheterm(length)
      real*8, allocatable :: dx(:)
      real*8, allocatable :: dy(:)
      real*8, allocatable :: dz(:)
      real*8, allocatable :: ddx2(:)
      real*8, allocatable :: h(:)
      real*8, allocatable :: g(:,:)
      real*8, allocatable :: hdiag(:,:)
      logical exist,qsecd
      character*240 xyzfile
c
c
c     set up the structure and mechanics calculation
c
      call initial
      xyzfile = 'mslave'
      call basefile (xyzfile)
      call suffix (xyzfile,'xyz')
      call version (xyzfile,'old')
      inquire (file=xyzfile,exist=exist)
      if (exist) then
         coordtype = 'CARTESIAN'
         ixyz = freeunit ()
         open (unit=ixyz,file=xyzfile,status='old')
         rewind (unit=ixyz)
         call readxyz (ixyz)
         close (unit=ixyz)
      else
         call getxyz
      end if
      call mechanic
c
c     perform dynamic allocation of some local arrays
c
      nhess = 3*n*(3*n-1) / 2
      allocate (hindex(nhess))
      allocate (hinit(3,n))
      allocate (hstop(3,n))
      allocate (dx(n))
      allocate (dy(n))
      allocate (dz(n))
      allocate (g(3,n))
      allocate (hdiag(3,n))
      allocate (h(nhess))
c
c     initialization of various MPI stuff
c
      mpidp = MPI_DOUBLE_PRECISION
      mpiint = MPI_INTEGER
      mpilog = MPI_LOGICAL
      call mpi_init (ierr)
      call mpi_comm_get_parent (parent,ierr)
      call mpi_comm_rank (parent,irank,ierr)
      if (ierr .eq. 0) then
         write (iout,'(a,i3)')  ' Value of PARENT = ',parent
      else
         write (iout,'(a,i3)')  ' Call to set PARENT returned ',ierr
      end if

c     receive the coordinates from the master process
c
      dowhile (.true.)
         meroot = 0
         call mpi_bcast (ncontrol,i,mpiint,meroot,parent,ierr)
         if (ncontrol .lt. 0) then
            call final
            stop
         end if
         call mpi_bcast (x,n,mpidp,meroot,parent,ierr)
         call mpi_bcast (y,n,mpidp,meroot,parent,ierr)
         call mpi_bcast (z,n,mpidp,meroot,parent,ierr)
         call mpi_bcast (qsecd,1,mpilog,meroot,parent,ierr)
c
c     compute the TINKER energy and gradient values
c
         call gradient (e,g)
c
c     copy the gradient into 1-D arrays for CHARMM
c
         do i = 1, n
            dx(i) = g(1,i)
            dy(i) = g(2,i)
            dz(i) = g(3,i)
         end do
c
c     transfer the gradient values to the MSCALE master
c
         if (irank .eq. 0) then
            meroot = MPI_ROOT
         else
            meroot = MPI_PROC_NULL
         end if
         call mpi_bcast (dx,n,mpidp,meroot,parent,ierr)
         call mpi_bcast (dy,n,mpidp,meroot,parent,ierr)
         call mpi_bcast (dz,n,mpidp,meroot,parent,ierr)
c
c     copy the energy terms into a 1-D array for CHARMM
c
         do i = 1, length
            cheterm(i) = 0.0d0
         end do
         cheterm(1) = eb
         cheterm(2) = ea + eba + eub + eaa + eopd + eopb + eid + eit
         cheterm(4) = et + ept + ebt + ett
         cheterm(6) = ev
         cheterm(7) = ec + ecd + ed + em + ep
         cheterm(60) = er + es
c
c     transfer the energy values to the MSCALE master
c
         call mpi_bcast (cheterm,length,mpidp,meroot,parent,ierr)
c
c     print the energy terms as a check of passed values
c
c        write (iout,*)  'TINKER> Total Energy = ',e
c        write (iout,*)  'TINKER> Bond Energy = ',cheterm(1)
c        write (iout,*)  'TINKER> Angl Energy = ',cheterm(2)
c        write (iout,*)  'TINKER> Dihe Energy = ',cheterm(4)
c        write (iout,*)  'TINKER> vDW  Energy = ',cheterm(6)
c        write (iout,*)  'TINKER> Elec Energy = ',cheterm(7)
c        write (iout,*)  'TINKER> rdp  Energy = ',cheterm(60)
c
c     compute the TINKER Hessian matrix values
c
         if (qsecd) then
            hesscut = 0.0d0
            call hessian (h,hinit,hstop,hindex,hdiag)
            n3 = 3 * n
            n6 = (n3*(n3+1)) / 2
            allocate (ddx2(n6))
c
c     copy the Hessian into a 1-D array for CHARMM
c
            m = 0
            do i = 1, n
               do j = 1, 3
                  m = m + 1
                  ddx2(m) = hdiag(j,i)
                  do k = hinit(j,i), hstop(j,i)
                     m = m + 1
                     ddx2(m) = h(k)
                  end do
               end do
            end do
c
c     transfer the Hessian values to the MSCALE master
c
c           write(*,'(A,I4)') 'TINKER> BCAST DDX2, N6 = ', n6
            call mpi_bcast (ddx2,n6,mpidp,meroot,parent,ierr)
            deallocate (ddx2)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dx)
      deallocate (dy)
      deallocate (dz)
      deallocate (g)
      deallocate (hdiag)
      deallocate (h)
c
c     perform any final tasks before program exit
c
      call final
      end

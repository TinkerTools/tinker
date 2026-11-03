c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_rattle  --  RATTLE bond constraints  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_rattle" takes one Verlet step from rest on a solvated
c     DHFR box with every bond and the waters frozen, and checks the
c     positions and velocities of RATTLE and the virial of RATTLE2,
c     as in the tinker-gpu rattle.cpp test
c
c     the references hold the same step from "rattle_step" and are
c     saved as gradient rows; the system has no potential terms, so
c     the step only moves atoms onto their constraints
c
c
      subroutine test_rattle
      use atoms
      use moldyn
      use virial
      implicit none
      integer i,j,nat
!$    integer nsave
!$    integer omp_get_max_threads
      real*8 ref_e,ref_ei
      real*8 refv(3,3)
      real*8, allocatable :: cur(:,:)
      real*8, allocatable :: refq(:,:)
      real*8, allocatable :: refp(:,:)
      logical skiptest
      character*240 rpath
      character*(*) tname
      parameter (tname='test_rattle')
c
c
      if (skiptest(tname,'rattle freeze'))  return
c
c     the keyfile asks for one OpenMP thread, so restore the count
c
!$    nsave = omp_get_max_threads ()
      call pushdir ('file/rattle')
      call rattle_step
      allocate (cur(3,n))
      allocate (refq(3,n))
      allocate (refp(3,n))
      call refpath ('rattle','rattle.q.txt',rpath)
      call read_gradient (rpath,n,refq,nat)
      call assert_int (nat,n,tname//' position count')
      call refpath ('rattle','rattle.p.txt',rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refp,nat)
      call assert_int (nat,n,tname//' velocity count')
      do i = 1, n
         cur(1,i) = x(i)
         cur(2,i) = y(i)
         cur(3,i) = z(i)
      end do
      call assert_grad (cur,refq,n,1.0d-6,tname//' positions')
      call assert_grad (v,refp,n,1.0d-5,tname//' velocities')
      call assert_grad (vir,refv,3,1.0d-4,tname//' virial')
      deallocate (cur)
      deallocate (refq)
      deallocate (refp)
      call popdir
      call sd_final
!$    call omp_set_num_threads (nsave)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine rattle_step  --  one constrained MD step  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "rattle_step" sets up the DHFR box at rest and takes one Verlet
c     step of 2 fs with the dynamics output silenced and nothing saved
c
c
      subroutine rattle_step
      use inform
      use iounit
      implicit none
      integer isave,freeunit
      real*8 dt
c
c
      call loadfix ('dhfr','dhfr.key')
      call sd_mdinit (0.0d0)
      iwrite = 1000000
      iprint = 1000000
      dt = 0.002d0
      isave = iout
      iout = freeunit ()
      open (unit=iout,status='scratch')
      call verlet (1,dt)
      close (unit=iout)
      iout = isave
      return
      end

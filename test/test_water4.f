c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_water4  --  rigid four-site water  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_water4" checks the extra sites of rigid four-site TIP4P
c     waters under FREEZE WATER, which are placed on their O-H-H
c     frames by "watfour" and pass their gradient to the O-H-H atoms
c     in "watfour2"; it covers the cases of the tinker-gpu water4.cpp
c     tests
c
c
      subroutine test_water4
      implicit none
c
c
      call test_water4_verlet
      call test_water4_respa
      call test_water4_grad
      call test_water4_minimize
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_water4_verlet  --  one step from rest  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_water4_verlet" takes one Verlet step from rest and checks
c     the positions and velocities against "dynamic" for the same
c     step, then checks that the waters are rigid, with the extra
c     sites on their frames and at rest
c
c     the references hold one step of "dynamic tip4p 1 1.0 0.001 1"
c     with "integrator verlet" and the extra sites placed on their
c     frames at load, shared with the tinker-gpu water4.cpp test,
c     which saves positions and velocities as gradient rows
c
c
      subroutine test_water4_verlet
      use freeze
      use inform
      use iounit
      implicit none
      integer isave,freeunit
      real*8 dt
      logical skiptest
      character*(*) tname
      parameter (tname='test_water4_verlet')
c
c
      if (skiptest(tname,'water4 freeze'))  return
      call pushdir ('file/water4')
      call loadfix ('tip4p','tip4p.key')
      call assert_int (nwat4,1000,tname//' four-site waters')
c
c     take one step with the dynamics output silenced and nothing
c     saved, from the velocities "mdinit" sets at zero kelvin
c
      call sd_mdinit (0.0d0)
      iwrite = 1000000
      iprint = 1000000
      dt = 0.001d0
      isave = iout
      iout = freeunit ()
      open (unit=iout,status='scratch')
      call verlet (1,dt)
      close (unit=iout)
      iout = isave
      call w4check (tname,'water4.q.txt','water4.p.txt')
      call popdir
      call sd_final
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_water4_respa  --  RESPA steps from rest  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_water4_respa" takes three VRESPA steps from rest with four
c     inner steps each, which constrain the waters and place the extra
c     sites at every inner step, and checks them as for one Verlet step
c
c     the references come from "vrespa" with the extra sites placed on
c     their frames at load, shared with the tinker-gpu water4.cpp test
c
c
      subroutine test_water4_respa
      use freeze
      use inform
      use iounit
      use mdstuf
      implicit none
      integer k,isave,freeunit
      real*8 dt
      logical skiptest
      character*(*) tname
      parameter (tname='test_water4_respa')
c
c
      if (skiptest(tname,'water4 freeze respa'))  return
      call pushdir ('file/water4')
      call loadfix_keyadd ('tip4p','tip4p.key','integrator vrespa')
      call assert_int (nwat4,1000,tname//' four-site waters')
c
c     the integrator must be known to "mdinit", which sets the fast
c     and slow accelerations; the inner steps only matter to "vrespa"
c
      call sd_mdinit (0.0d0)
      nrespa = 4
      iwrite = 1000000
      iprint = 1000000
      dt = 0.001d0
      isave = iout
      iout = freeunit ()
      open (unit=iout,status='scratch')
      do k = 1, 3
         call vrespa (k,dt)
      end do
      close (unit=iout)
      iout = isave
      call w4check (tname,'water4respa.q.txt','water4respa.p.txt')
      call popdir
      call sd_final
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine w4check  --  rigid four-site water checks  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "w4check" compares the positions and velocities with references
c     saved as gradient rows, then checks that the waters are rigid,
c     with the extra sites on their frames and at rest
c
c
      subroutine w4check (tname,qname,pname)
      use atoms
      use freeze
      use math
      use moldyn
      implicit none
      integer i,j,nat
      integer ia,ib,ic,id
      real*8 eps
      real*8 roh,rhh
      real*8 dmax,vmax
      real*8 rmax,wmax
      real*8 xr(3)
      real*8, allocatable :: cur(:,:)
      real*8, allocatable :: refq(:,:)
      real*8, allocatable :: refp(:,:)
      character*240 rpath
      character*(*) tname,qname,pname
c
c
c     compare positions under the minimum image, in case a molecule
c     is wrapped into the box at a different time
c
      allocate (cur(3,n))
      allocate (refq(3,n))
      allocate (refp(3,n))
      call refpath ('water4',qname,rpath)
      call read_gradient (rpath,n,refq,nat)
      call refpath ('water4',pname,rpath)
      call read_gradient (rpath,n,refp,nat)
      do i = 1, n
         cur(1,i) = x(i)
         cur(2,i) = y(i)
         cur(3,i) = z(i)
         do j = 1, 3
            xr(j) = cur(j,i) - refq(j,i)
         end do
         call image (xr(1),xr(2),xr(3))
         do j = 1, 3
            refq(j,i) = cur(j,i) - xr(j)
         end do
      end do
      eps = 1.0d-8
      call assert_grad (cur,refq,n,eps,tname//' positions')
      call assert_grad (v,refp,n,eps,tname//' velocities')
c
c     each extra site sits on its O-H-H frame with no velocity
c
      dmax = 0.0d0
      vmax = 0.0d0
      do i = 1, nwat4
         ia = iwat4(1,i)
         ib = iwat4(2,i)
         ic = iwat4(3,i)
         id = iwat4(4,i)
         xr(1) = x(id) - kwat4(1,i)*x(ia) - kwat4(2,i)*(x(ib)+x(ic))
         xr(2) = y(id) - kwat4(1,i)*y(ia) - kwat4(2,i)*(y(ib)+y(ic))
         xr(3) = z(id) - kwat4(1,i)*z(ia) - kwat4(2,i)*(z(ib)+z(ic))
         do j = 1, 3
            dmax = max(dmax,abs(xr(j)))
            vmax = max(vmax,abs(v(j,id)))
         end do
      end do
      call assert_real (dmax,0.0d0,1.0d-12,tname//' extra site frame')
      call assert_real (vmax,0.0d0,0.0d0,tname//' extra site at rest')
c
c     O-H and H-H distances keep the TIP4P geometry of tip4p.prm,
c     and no atom pair moves along the line joining it
c
      roh = 0.9572d0
      rhh = 2.0d0 * roh * sin(0.5d0*104.52d0/radian)
      rmax = 0.0d0
      wmax = 0.0d0
      do i = 1, nwat
         ia = iwat(1,i)
         ib = iwat(2,i)
         ic = iwat(3,i)
         call w4pair (ia,ib,roh,rmax,wmax)
         call w4pair (ia,ic,roh,rmax,wmax)
         call w4pair (ib,ic,rhh,rmax,wmax)
      end do
      call assert_real (rmax,0.0d0,1.0d-8,tname//' rigid distances')
      call assert_real (wmax,0.0d0,1.0d-8,tname//' rigid velocities')
      deallocate (cur)
      deallocate (refq)
      deallocate (refp)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine w4pair  --  constraint error of an atom pair  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "w4pair" updates the largest deviation from an ideal distance
c     and the largest velocity along the line joining an atom pair
c
c
      subroutine w4pair (ia,ib,r0,rmax,wmax)
      use atoms
      use moldyn
      implicit none
      integer ia,ib
      real*8 r0,rmax,wmax
      real*8 xr,yr,zr,r
      real*8 vxr,vyr,vzr
c
c
      xr = x(ia) - x(ib)
      yr = y(ia) - y(ib)
      zr = z(ia) - z(ib)
      vxr = v(1,ia) - v(1,ib)
      vyr = v(2,ia) - v(2,ib)
      vzr = v(3,ia) - v(3,ib)
      r = sqrt(xr*xr + yr*yr + zr*zr)
      rmax = max(rmax,abs(r-r0))
      wmax = max(wmax,abs(xr*vxr+yr*vyr+zr*vzr)/r)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_water4_grad  --  extra site gradient  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_water4_grad" checks that the extra sites are put on their
c     frames when the coordinates are read, that "gradient" leaves no
c     gradient on them, in "derivs" and in "desum", and that what it
c     gives the O-H-H atoms of a water is the derivative of the energy
c     when the extra site moves with them
c
c
      subroutine test_water4_grad
      use atoms
      use deriv
      use freeze
      implicit none
      integer i,j,k,ia
      integer ib,ic,id
      real*8 e,ep,em,h
      real*8 dmax,energy
      real*8 an(3,3),fd(3,3)
      real*8, allocatable :: derivs(:,:)
      logical skiptest
      character*(*) tname
      parameter (tname='test_water4_grad')
c
c
      if (skiptest(tname,'water4 freeze'))  return
      call pushdir ('file/water4')
      call loadfix ('tip4p','tip4p.key')
c
c     the extra sites of tip4p.xyz are off their frames in the file
c
      dmax = 0.0d0
      do i = 1, nwat4
         ia = iwat4(1,i)
         ib = iwat4(2,i)
         ic = iwat4(3,i)
         id = iwat4(4,i)
         dmax = max(dmax,abs(x(id)-kwat4(1,i)*x(ia)
     &                          -kwat4(2,i)*(x(ib)+x(ic))))
         dmax = max(dmax,abs(y(id)-kwat4(1,i)*y(ia)
     &                          -kwat4(2,i)*(y(ib)+y(ic))))
         dmax = max(dmax,abs(z(id)-kwat4(1,i)*z(ia)
     &                          -kwat4(2,i)*(z(ib)+z(ic))))
      end do
      call assert_real (dmax,0.0d0,1.0d-12,tname//' frames at load')
      allocate (derivs(3,n))
      call gradient (e,derivs)
      dmax = 0.0d0
      do i = 1, nwat4
         k = iwat4(4,i)
         do j = 1, 3
            dmax = max(dmax,abs(derivs(j,k)))
         end do
      end do
      call assert_real (dmax,0.0d0,0.0d0,tname//' extra site rows')
      call assert_grad (desum,derivs,n,0.0d0,tname//' desum')
c
c     central differences over the O-H-H atoms of the first water,
c     placing its extra site again after each displacement
c
      h = 1.0d-4
      do k = 1, 3
         ia = iwat4(k,1)
         do j = 1, 3
            an(j,k) = derivs(j,ia)
            call w4shift (ia,j,h)
            ep = energy ()
            call w4shift (ia,j,-2.0d0*h)
            em = energy ()
            call w4shift (ia,j,h)
            fd(j,k) = (ep-em) / (2.0d0*h)
         end do
      end do
      call assert_grad (an,fd,3,1.0d-4,tname//' frame atoms')
      deallocate (derivs)
      call popdir
      call final
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine w4shift  --  move one coordinate of an atom  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "w4shift" displaces one Cartesian coordinate of an atom and puts
c     the extra sites of four-site waters back on their frames
c
c
      subroutine w4shift (i,j,d)
      use atoms
      implicit none
      integer i,j
      real*8 d
c
c
      if (j .eq. 1)  x(i) = x(i) + d
      if (j .eq. 2)  y(i) = y(i) + d
      if (j .eq. 3)  z(i) = z(i) + d
      call watfour (0.0d0)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_water4_minimize  --  first iteration  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_water4_minimize" runs "minimize" with a gradient criterion
c     the starting structure meets, so it stops after the function it
c     optimizes is found once; "minimiz1" first moves the waters onto
c     their constraints and then projects the gradient off them
c
c     the references are iteration 0 of testmanual/tip4p.log
c
c
      subroutine test_water4_minimize
      implicit none
      integer ist
      real*8 e,grms,err
      logical found,skiptest
      character*(*) tname
      parameter (tname='test_water4_minimize')
c
c
      if (skiptest(tname,'water4 freeze'))  return
      call pushdir ('file/water4')
      call execute_command_line ('rm -f tip4p.xyz_* out.txt')
      call run_prog ('minimize','tip4p 2.0','out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call barfind ('out.txt','SmallGrad',found,e,err)
      call assert_logical (found,.true.,tname//' stops at iteration 0')
      call barfind ('out.txt','Final Function Value :',found,e,err)
      call assert_real (e,-12336.9614d0,1.0d-4,tname//' energy')
      call barfind ('out.txt','Final RMS Gradient :',found,grms,err)
      call assert_real (grms,1.1078d0,1.0d-4,tname//' rms gradient')
      call execute_command_line ('rm -f tip4p.xyz_* out.txt')
      call popdir
      return
      end

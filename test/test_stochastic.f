c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_stochastic  --  stochastic dynamics tests  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_stochastic" runs the stochastic dynamics tests mirrored
c     from the tinker-gpu stochastic.cpp test cases; each dynamics
c     test uses the 20-atom TEMOA guest G3 with AMOEBA polarization
c     in up to four cells, with or without a periodic box and with
c     or without removal of inertia every three steps
c
c     the periodic cells use a 30 Angstrom box with 6 Angstrom
c     cutoffs and no Ewald, so the box changes the forces; they
c     restart from g3pbc.dyn, which is g3.dyn with the box set and
c     accelerations consistent with the periodic keywords
c
c
      subroutine test_stochastic
      implicit none
c
c
      call test_sd_coef
      call test_sd_setup
      call test_sd_mdrest
      call test_sd_zerofric
      call test_sd_momentum
      call test_sd_trajectory
      call test_sd_fricsign
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_sd_coef  --  friction coefficients  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_sd_coef" checks the frictional terms from "sdterm"
c     against values computed in 50-digit arithmetic for a time
c     step of 1 fs; the friction is set separately for each atom,
c     so one call spans the series expansion used below a value
c     of 0.05 for friction times time step, the closed form used
c     above it, the seam between them and the zero friction limit
c
c
      subroutine test_sd_coef
      use atoms
      use stodyn
      implicit none
      integer nref
      parameter (nref=8)
      integer i,j
      real*8 dt,rmax
      real*8 gref(nref),tol(nref)
      real*8 pref(nref)
      real*8 vref(nref)
      real*8 aref(nref)
      real*8, allocatable :: pfric(:)
      real*8, allocatable :: vfric(:)
      real*8, allocatable :: afric(:)
      real*8, allocatable :: prand(:,:)
      real*8, allocatable :: vrand(:,:)
      logical skiptest
      character*1 tag
      character*(*) tname
      parameter (tname='test_sd_coef')
      data gref  / 0.0001d0, 0.001d0, 0.005d0, 0.02d0, 0.049d0,
     &             0.05d0, 0.5d0, 5.0d0 /
      data tol   / 1.0d-13, 1.0d-13, 1.0d-13, 1.0d-13, 1.0d-12,
     &             1.0d-11, 1.0d-13, 1.0d-13 /
      data pref  / 9.9990000499983d-01, 9.9900049983338d-01,
     &             9.9501247919268d-01, 9.8019867330676d-01,
     &             9.5218112969850d-01, 9.5122942450071d-01,
     &             6.0653065971263d-01, 6.7379469990855d-03 /
      data vref  / 9.9995000166663d-04, 9.9950016662501d-04,
     &             9.9750416146354d-04, 9.9006633466223d-04,
     &             9.7589531227541d-04, 9.7541150998572d-04,
     &             7.8693868057473d-04, 1.9865241060018d-04 /
      data aref  / 4.9998333374999d-07, 4.9983337499167d-07,
     &             4.9916770729253d-07, 4.9668326688826d-07,
     &             4.9193240254263d-07, 4.9176980028560d-07,
     &             4.2612263885053d-07, 1.6026951787996d-07 /
c
c
      if (skiptest(tname,'stochastic'))  return
      call sd_prep ('stochastic_coef',.false.,.false.)
      call loadfix ('g3','g3.key')
      call sd_mdinit (298.0d0)
      dt = 0.001d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (pfric(n))
      allocate (vfric(n))
      allocate (afric(n))
      allocate (prand(3,n))
      allocate (vrand(3,n))
c
c     the first call in a process resets the atomic friction,
c     so make one call before setting the values to be tested
c
      if (.not. allocated(fgamma))  allocate (fgamma(n))
      call sdterm (1,dt,pfric,vfric,afric,prand,vrand)
c
c     reference rows, then both sides of the seam, then the
c     zero friction limit and a value just above that limit
c
      do i = 1, nref
         fgamma(i) = gref(i) / dt
      end do
      fgamma(9) = 0.05d0 * (1.0d0+1.0d-12) / dt
      fgamma(10) = 0.05d0 * (1.0d0-1.0d-12) / dt
      fgamma(11) = 0.0d0
      fgamma(12) = 1.0d-8 / dt
      call sdterm (2,dt,pfric,vfric,afric,prand,vrand)
c
c     compare with the reference via the relative error; the
c     tolerance is looser near the seam, where the truncated
c     series and the closed form are each at their worst
c
      do i = 1, nref
         write (tag,10)  i
   10    format (i1)
         call assert_real (pfric(i)/pref(i),1.0d0,tol(i),
     &                     tname//' pfric '//tag)
         call assert_real (vfric(i)/vref(i),1.0d0,tol(i),
     &                     tname//' vfric '//tag)
         call assert_real (afric(i)/aref(i),1.0d0,tol(i),
     &                     tname//' afric '//tag)
      end do
c
c     the series and the closed form must agree across the seam
c
      call assert_real (pfric(10)/pfric(9),1.0d0,1.0d-9,
     &                  tname//' seam pfric')
      call assert_real (vfric(10)/vfric(9),1.0d0,1.0d-9,
     &                  tname//' seam vfric')
      call assert_real (afric(10)/afric(9),1.0d0,1.0d-9,
     &                  tname//' seam afric')
c
c     zero friction is velocity Verlet without any random terms
c
      rmax = 0.0d0
      do j = 1, 3
         rmax = max(rmax,abs(prand(j,11)),abs(vrand(j,11)))
      end do
      call assert_real (pfric(11),1.0d0,0.0d0,tname//' zero pfric')
      call assert_real (vfric(11),dt,0.0d0,tname//' zero vfric')
      call assert_real (afric(11),0.5d0*dt*dt,0.0d0,
     &                  tname//' zero afric')
      call assert_real (rmax,0.0d0,0.0d0,tname//' zero random')
c
c     and the series must approach that limit continuously
c
      call assert_real (pfric(12),1.0d0,1.0d-7,tname//' tiny pfric')
      call assert_real (vfric(12)/dt,1.0d0,1.0d-7,
     &                  tname//' tiny vfric')
      call assert_real (afric(12)/(0.5d0*dt*dt),1.0d0,1.0d-7,
     &                  tname//' tiny afric')
c
c     perform deallocation of some local arrays
c
      deallocate (pfric)
      deallocate (vfric)
      deallocate (afric)
      deallocate (prand)
      deallocate (vrand)
      call sd_final
      call sd_clean ('stochastic_coef')
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_sd_setup  --  inertia removal setup  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_sd_setup" checks what "mdinit" decides for stochastic
c     dynamics; REMOVE-INERTIA N turns removal on every N steps and
c     takes the removed modes out of the degrees of freedom, while
c     without the keyword nothing is removed and all modes count;
c     an explicit DEGREES-FREEDOM value is taken as given
c
c
      subroutine test_sd_setup
      use mdstuf
      implicit none
      integer ic,idf
      integer wantrest,wantfree
      logical pbc,rest
      logical skiptest
      character*12 cell,sd_cell
      character*40 label
      character*(*) tname
      parameter (tname='test_sd_setup')
c
c
      if (skiptest(tname,'stochastic'))  return
      do ic = 1, 4
         pbc = (ic .ge. 3)
         rest = (mod(ic,2) .eq. 0)
         cell = sd_cell (pbc,rest)
         do idf = 1, 2
            call sd_prep ('stochastic_setup',pbc,rest)
            if (idf .eq. 2)  call sd_addkey ('DEGREES-FREEDOM 50')
            call loadfix ('g3','g3.key')
            call sd_mdinit (298.0d0)
            wantrest = 0
            wantfree = 60
            if (rest) then
               wantrest = 3
               wantfree = 54
               if (pbc)  wantfree = 57
            end if
            label = tname//' '//trim(cell)
            if (idf .eq. 2) then
               wantfree = 50
               label = trim(label)//' dof'
            end if
            call assert_logical (dorest,rest,trim(label)//' dorest')
            call assert_int (irest,wantrest,trim(label)//' irest')
            call assert_int (nfree,wantfree,trim(label)//' nfree')
            call sd_final
            call sd_clean ('stochastic_setup')
         end do
      end do
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_sd_mdrest  --  removal of inertia  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_sd_mdrest" checks that "mdrest" leaves the velocities
c     alone unless removal is on and due at the current step, and
c     that it then zeroes the linear momentum and, without a box,
c     the angular momentum about the center of mass; saved frames
c     are written before removal, so this is checked in memory
c
c
      subroutine test_sd_mdrest
      use atomid
      use atoms
      use mdstuf
      use moldyn
      implicit none
      integer i,j,ib
      real*8 ps,ls,ps0,ls0
      real*8 diff,sd_norm
      real*8 p(3),l(3)
      real*8 p0(3),l0(3)
      real*8, allocatable :: xyz(:,:)
      real*8, allocatable :: vel(:,:)
      logical pbc
      logical skiptest
      character*12 cell,sd_cell
      character*40 label
      character*(*) tname
      parameter (tname='test_sd_mdrest')
c
c
      if (skiptest(tname,'stochastic'))  return
      do ib = 1, 2
         pbc = (ib .eq. 2)
         cell = sd_cell (pbc,.true.)
         label = tname//' '//trim(cell)
         call sd_prep ('stochastic_mdrest',pbc,.true.)
         call loadfix ('g3','g3.key')
         call sd_mdinit (298.0d0)
         allocate (xyz(3,n))
         allocate (vel(3,n))
         do i = 1, n
            xyz(1,i) = x(i)
            xyz(2,i) = y(i)
            xyz(3,i) = z(i)
            do j = 1, 3
               vel(j,i) = v(j,i)
            end do
         end do
         call sd_momenta (n,mass,xyz,vel,p0,l0,ps0,ls0)
c
c     removal is real only if there is something to remove
c
         call assert_logical (sd_norm(p0).gt.1.0d-3*ps0,.true.,
     &                        trim(label)//' linear present')
         call assert_logical (sd_norm(l0).gt.1.0d-3*ls0,.true.,
     &                        trim(label)//' angular present')
c
c     nothing happens if removal is off or not due at this step
c
         dorest = .false.
         call mdrest (3)
         dorest = .true.
         call mdrest (2)
         diff = 0.0d0
         do i = 1, n
            do j = 1, 3
               diff = max(diff,abs(v(j,i)-vel(j,i)))
            end do
         end do
         call assert_real (diff,0.0d0,0.0d0,trim(label)//' gate')
c
c     removal zeroes the linear momentum and, without a box, the
c     angular momentum; in a box the angular momentum is kept
c
         call mdrest (3)
         do i = 1, n
            do j = 1, 3
               vel(j,i) = v(j,i)
            end do
         end do
         call sd_momenta (n,mass,xyz,vel,p,l,ps,ls)
         call assert_real (sd_norm(p)/ps0,0.0d0,1.0d-12,
     &                     trim(label)//' linear removed')
         if (pbc) then
            do j = 1, 3
               l(j) = l(j) - l0(j)
            end do
            call assert_real (sd_norm(l)/ls0,0.0d0,1.0d-12,
     &                        trim(label)//' angular kept')
         else
            call assert_real (sd_norm(l)/ls0,0.0d0,1.0d-12,
     &                        trim(label)//' angular removed')
         end if
         deallocate (xyz)
         deallocate (vel)
         call sd_final
         call sd_clean ('stochastic_mdrest')
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_sd_zerofric  --  zero friction limit  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_sd_zerofric" checks that stochastic dynamics with zero
c     friction reproduces velocity Verlet step for step, including
c     the removal steps; the linear momentum and, without a box,
c     the angular momentum must then stay at their starting values,
c     or at zero in frames saved after the first removal
c
c
      subroutine test_sd_zerofric
      implicit none
      integer nat,nstep
      parameter (nat=20,nstep=10)
      integer i,j,ic
      integer ist1,ist2
      integer nv,ns,nx,nw
      real*8 emax,pmax,lmax
      real*8 ps,ls,ref
      real*8 m(nat)
      real*8 p(3),l(3)
      real*8 p1(3),l1(3)
      real*8 vpot(nstep),vkin(nstep)
      real*8 spot(nstep),skin(nstep)
      real*8 xyz(3,nat,nstep)
      real*8 vel(3,nat,nstep)
      logical pbc,rest
      logical skiptest
      character*12 cell,sd_cell
      character*40 label
      character*(*) tname,work,args
      parameter (tname='test_sd_zerofric')
      parameter (work='stochastic_zerofric')
      parameter (args='g3 10 1.0 0.001 1')
c
c
      if (skiptest(tname,'stochastic'))  return
      call sd_masses (nat,m)
      do ic = 1, 4
         pbc = (ic .ge. 3)
         rest = (mod(ic,2) .eq. 0)
         cell = sd_cell (pbc,rest)
         label = tname//' '//trim(cell)
c
c     velocity Verlet removes inertia every 100 steps by default,
c     so turn it off explicitly to match the stochastic setting
c
         call sd_prep (work,pbc,rest)
         call sd_addkey ('INTEGRATOR VERLET')
         if (.not. rest)  call sd_addkey ('REMOVE-INERTIA 0')
         call sd_addkey ('DIGITS 8')
         call run_prog ('dynamic',args,'out.txt',ist1)
         call sd_energies ('out.txt',nstep,nv,vpot,vkin)
         call sd_clean (work)
         call sd_prep (work,pbc,rest)
         call sd_addkey ('FRICTION 0.0')
         call sd_addkey ('DIGITS 8')
         call sd_addkey ('SAVE-VELOCITY')
         call run_prog ('dynamic',args,'out.txt',ist2)
         call sd_energies ('out.txt',nstep,ns,spot,skin)
         call sd_frames ('g3.arc',nat,nstep,nx,xyz)
         call sd_frames ('g3.vel',nat,nstep,nw,vel)
         call sd_clean (work)
         if (ist1.ne.-1 .and. ist2.ne.-1) then
            call assert_int (ist1,0,trim(label)//' verlet status')
            call assert_int (ist2,0,trim(label)//' status')
            call assert_int (nv,nstep,trim(label)//' verlet steps')
            call assert_int (ns,nstep,trim(label)//' steps')
            call assert_int (nx,nstep,trim(label)//' arc frames')
            call assert_int (nw,nstep,trim(label)//' vel frames')
         end if
         if (min(nv,ns,nx,nw) .eq. nstep) then
            emax = 0.0d0
            do i = 1, nstep
               emax = max(emax,abs(spot(i)-vpot(i)),
     &                       abs(skin(i)-vkin(i)))
            end do
            call assert_real (emax,0.0d0,1.0d-6,
     &                        trim(label)//' energy')
c
c     frames up to the first removal step are saved before any
c     removal, and the later ones after it; the tolerance on the
c     angular momentum allows for the net torque that is left by
c     induced dipoles converged to a finite POLAR-EPS
c
            pmax = 0.0d0
            lmax = 0.0d0
            call sd_momenta (nat,m,xyz(1,1,1),vel(1,1,1),
     &                       p1,l1,ps,ls)
            do i = 2, nstep
               call sd_momenta (nat,m,xyz(1,1,i),vel(1,1,i),
     &                          p,l,ps,ls)
               do j = 1, 3
                  ref = p1(j)
                  if (rest .and. i.gt.3)  ref = 0.0d0
                  pmax = max(pmax,abs(p(j)-ref)/ps)
                  ref = l1(j)
                  if (rest .and. i.gt.3)  ref = 0.0d0
                  lmax = max(lmax,abs(l(j)-ref)/ls)
               end do
            end do
            call assert_real (pmax,0.0d0,1.0d-6,
     &                        trim(label)//' linear')
            if (.not. pbc) then
               call assert_real (lmax,0.0d0,2.0d-6,
     &                           trim(label)//' angular')
            end if
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_sd_momentum  --  momentum with friction  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_sd_momentum" checks the momentum of a stochastic run
c     at 20 ps-1 friction; without removal, the linear momentum and,
c     without a box, the angular momentum of each saved frame are
c     compared with recorded values; these were recorded from this
c     code at eight decimals and agree with Tinker 8.10.5 to the six
c     digits its velocity files carry
c
c     with uniform friction the internal forces cancel in the sum
c     over atoms, so the linear momentum obeys p' = pfric*p + noise;
c     the noise is the same with and without removal, so the frames
c     of a removal run must differ from those without removal by
c     the decayed momentum that was last taken out
c
c
      subroutine test_sd_momentum
      implicit none
      integer nat,nstep
      parameter (nat=20,nstep=10)
      integer i,j,k,ib,ir
      integer ist,nx,nw
      real*8 pfric,pred
      real*8 pmax,lmax,rmax,pmin
      real*8 sd_norm
      real*8 m(nat)
      real*8 ps(nstep,2),ls(nstep,2)
      real*8 p(3,nstep,2),l(3,nstep,2)
      real*8 pref(3,nstep,2),lref(3,nstep)
      real*8 xyz(3,nat,nstep)
      real*8 vel(3,nat,nstep)
      logical pbc,rest,ok
      logical skiptest
      character*12 cell,sd_cell
      character*40 label
      character*(*) tname,work,args
      parameter (tname='test_sd_momentum')
      parameter (work='stochastic_momentum')
      parameter (args='g3 10 1.0 0.001 2 298')
c
c
      if (skiptest(tname,'stochastic'))  return
      call sd_masses (nat,m)
      call sd_ref_momentum (pref,lref)
      pfric = exp(-20.0d0*0.001d0)
      do ib = 1, 2
         pbc = (ib .eq. 2)
         ok = .true.
         do ir = 1, 2
            rest = (ir .eq. 2)
            cell = sd_cell (pbc,rest)
            label = tname//' '//trim(cell)
            call sd_prep (work,pbc,rest)
            call sd_addkey ('FRICTION 20.0')
            call sd_addkey ('DIGITS 8')
            call sd_addkey ('SAVE-VELOCITY')
            call run_prog ('dynamic',args,'out.txt',ist)
            call sd_frames ('g3.arc',nat,nstep,nx,xyz)
            call sd_frames ('g3.vel',nat,nstep,nw,vel)
            call sd_clean (work)
            if (ist .ne. -1) then
               call assert_int (ist,0,trim(label)//' status')
               call assert_int (nx,nstep,trim(label)//' arc frames')
               call assert_int (nw,nstep,trim(label)//' vel frames')
            end if
            if (min(nx,nw) .eq. nstep) then
               do i = 1, nstep
                  call sd_momenta (nat,m,xyz(1,1,i),vel(1,1,i),
     &                             p(1,i,ir),l(1,i,ir),
     &                             ps(i,ir),ls(i,ir))
               end do
            else
               ok = .false.
            end if
         end do
         if (ok) then
c
c     compare the run without removal against recorded values
c
            cell = sd_cell (pbc,.false.)
            label = tname//' '//trim(cell)
            pmax = 0.0d0
            lmax = 0.0d0
            do i = 1, nstep
               do j = 1, 3
                  pmax = max(pmax,abs(p(j,i,1)-pref(j,i,ib))
     &                               /ps(i,1))
                  lmax = max(lmax,abs(l(j,i,1)-lref(j,i))/ls(i,1))
               end do
            end do
            call assert_real (pmax,0.0d0,1.0d-6,
     &                        trim(label)//' linear')
            if (.not. pbc) then
               call assert_real (lmax,0.0d0,1.0d-6,
     &                           trim(label)//' angular')
            end if
c
c     predict the removal run from the run without removal, where
c     k is the last removal step before the current frame
c
            cell = sd_cell (pbc,.true.)
            label = tname//' '//trim(cell)
            rmax = 0.0d0
            pmin = 1.0d0
            do i = 1, nstep
               k = 3 * ((i-1)/3)
               if (k .ge. 3) then
                  pmin = min(pmin,sd_norm(p(1,k,1))/ps(k,1))
               end if
               do j = 1, 3
                  pred = p(j,i,1)
                  if (k .ge. 3) then
                     pred = pred - p(j,k,1)*pfric**(i-k)
                  end if
                  rmax = max(rmax,abs(p(j,i,2)-pred)/ps(i,1))
               end do
            end do
            call assert_logical (pmin.gt.1.0d-3,.true.,
     &                           trim(label)//' linear present')
            call assert_real (rmax,0.0d0,1.0d-7,
     &                        trim(label)//' linear removed')
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_sd_trajectory  --  match Tinker 8.10.5  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_sd_trajectory" compares a stochastic trajectory step by
c     step against Fortran Tinker 8.10.5, restarted from a dynamics
c     file so that no random numbers are drawn for the velocities;
c     the friction is set explicitly since 8.10.5 defaults to 91 ps-1
c
c     the reference never removes inertia; removal does not change
c     the energy of the step it happens on, so removal cells match
c     up to the first removal; in a box removing translation only
c     shifts every atom equally, so the potential energy matches at
c     every step and the final frame up to a common offset
c
c
      subroutine test_sd_trajectory
      implicit none
      integer nat,nstep
      parameter (nat=20,nstep=10)
      integer i,j,ic,ib
      integer ist,ns,nx,nsame
      real*8 umax,kmax,xmax
      real*8 off(3)
      real*8 pot(nstep),kin(nstep)
      real*8 gpot(nstep,2),gkin(nstep,2)
      real*8 gxyz(3,nat,2)
      real*8 xyz(3,nat,nstep)
      logical pbc,rest
      logical skiptest
      character*12 cell,sd_cell
      character*40 label
      character*(*) tname,work,args
      parameter (tname='test_sd_trajectory')
      parameter (work='stochastic_trajectory')
      parameter (args='g3 10 0.1 0.0001 2 298')
c
c
      if (skiptest(tname,'stochastic'))  return
      call sd_ref_trajectory (gpot,gkin,gxyz)
      do ic = 1, 4
         pbc = (ic .ge. 3)
         rest = (mod(ic,2) .eq. 0)
         cell = sd_cell (pbc,rest)
         label = tname//' '//trim(cell)
         ib = 1
         if (pbc)  ib = 2
         call sd_prep (work,pbc,rest)
         call sd_addkey ('FRICTION 91.0')
         call run_prog ('dynamic',args,'out.txt',ist)
         call sd_energies ('out.txt',nstep,ns,pot,kin)
         call sd_frames ('g3.arc',nat,nstep,nx,xyz)
         call sd_clean (work)
         if (ist .ne. -1) then
            call assert_int (ist,0,trim(label)//' status')
            call assert_int (ns,nstep,trim(label)//' steps')
            call assert_int (nx,nstep,trim(label)//' arc frames')
         end if
         if (min(ns,nx) .eq. nstep) then
c
c     the reference energies are printed to four decimals
c
            nsame = nstep
            if (rest)  nsame = 3
            umax = 0.0d0
            kmax = 0.0d0
            do i = 1, nstep
               if (i.le.nsame .or. pbc) then
                  umax = max(umax,abs(pot(i)-gpot(i,ib)))
               end if
               if (i .le. nsame) then
                  kmax = max(kmax,abs(kin(i)-gkin(i,ib)))
               end if
            end do
            call assert_real (umax,0.0d0,5.0d-4,
     &                        trim(label)//' potential')
            call assert_real (kmax,0.0d0,5.0d-4,
     &                        trim(label)//' kinetic')
c
c     energies can hide a single bad atom, so also check the final
c     frame atom by atom, which the reference has to six decimals
c
            if (pbc .or. .not.rest) then
               do j = 1, 3
                  off(j) = 0.0d0
                  if (rest) then
                     do i = 1, nat
                        off(j) = off(j) + (xyz(j,i,nstep)
     &                              -gxyz(j,i,ib))/dble(nat)
                     end do
                  end if
               end do
               xmax = 0.0d0
               do i = 1, nat
                  do j = 1, 3
                     xmax = max(xmax,abs(xyz(j,i,nstep)-off(j)
     &                                      -gxyz(j,i,ib)))
                  end do
               end do
               call assert_real (xmax,0.0d0,1.0d-5,
     &                           trim(label)//' final frame')
            end if
         end if
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_sd_fricsign  --  negative friction  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_sd_fricsign" checks that a negative FRICTION value has
c     its sign flipped, by reproducing the reference trajectory
c
c
      subroutine test_sd_fricsign
      implicit none
      integer nat,nstep
      parameter (nat=20,nstep=10)
      integer i,ist,ns
      real*8 emax
      real*8 pot(nstep),kin(nstep)
      real*8 gpot(nstep,2),gkin(nstep,2)
      real*8 gxyz(3,nat,2)
      logical skiptest
      character*(*) tname,work
      parameter (tname='test_sd_fricsign')
      parameter (work='stochastic_fricsign')
c
c
      if (skiptest(tname,'stochastic'))  return
      call sd_ref_trajectory (gpot,gkin,gxyz)
      call sd_prep (work,.false.,.false.)
      call sd_addkey ('FRICTION -91.0')
      call run_prog ('dynamic','g3 10 0.1 0.0001 2 298','out.txt',ist)
      call sd_energies ('out.txt',nstep,ns,pot,kin)
      call sd_clean (work)
      if (ist .ne. -1) then
         call assert_int (ist,0,tname//' status')
         call assert_int (ns,nstep,tname//' steps')
      end if
      if (ns .eq. nstep) then
         emax = 0.0d0
         do i = 1, nstep
            emax = max(emax,abs(pot(i)-gpot(i,1)),
     &                    abs(kin(i)-gkin(i,1)))
         end do
         call assert_real (emax,0.0d0,5.0d-4,tname//' energy')
      end if
      return
      end
c
c
c     #################################################
c     ##                                             ##
c     ##  function sd_cell  --  name of a test cell  ##
c     ##                                             ##
c     #################################################
c
c
c     "sd_cell" returns a short name for one of the four test cells
c
c
      character*12 function sd_cell (pbc,rest)
      implicit none
      logical pbc,rest
c
c
      if (pbc) then
         sd_cell = 'pbc_norest'
         if (rest)  sd_cell = 'pbc_rest'
      else
         sd_cell = 'nopbc_norest'
         if (rest)  sd_cell = 'nopbc_rest'
      end if
      return
      end
c
c
c     ######################################################
c     ##                                                  ##
c     ##  subroutine sd_prep  --  create scratch fixture  ##
c     ##                                                  ##
c     ######################################################
c
c
c     "sd_prep" copies the G3 fixture into a scratch directory with
c     the restart file and keywords of the requested cell, and then
c     enters it; "dynamic" overwrites its restart file, so each run
c     needs a fresh copy; "sd_clean" leaves and removes the directory
c
c
      subroutine sd_prep (work,pbc,rest)
      implicit none
      logical pbc,rest
      character*(*) work
      character*9 dyn
      character*512 cmd
c
c
      dyn = 'g3.dyn'
      if (pbc)  dyn = 'g3pbc.dyn'
      call pushdir ('file/stochastic')
      cmd = 'rm -rf ../'//trim(work)//' ; mkdir -p ../'//trim(work)//
     &      ' ; cp g3.xyz g3.key hostsG3.prm ../'//trim(work)//
     &      '/ ; cp '//trim(dyn)//' ../'//trim(work)//'/g3.dyn'
      call execute_command_line (cmd)
      call popdir
      call pushdir ('file/'//trim(work))
      if (pbc) then
         call sd_addkey ('A-AXIS 30.0')
         call sd_addkey ('VDW-CUTOFF 6.0')
         call sd_addkey ('MPOLE-CUTOFF 6.0')
      end if
      if (rest)  call sd_addkey ('REMOVE-INERTIA 3')
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine sd_clean  --  remove scratch fixture  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "sd_clean" leaves the scratch directory and removes it
c
c
      subroutine sd_clean (work)
      implicit none
      character*(*) work
      character*512 cmd
c
c
      call popdir
      call pushdir ('file')
      cmd = 'rm -rf '//trim(work)
      call execute_command_line (cmd)
      call popdir
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine sd_addkey  --  append a keyword line  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "sd_addkey" appends one keyword line to the scratch keyfile
c
c
      subroutine sd_addkey (line)
      implicit none
      integer unit,freeunit
      character*(*) line
c
c
      unit = freeunit ()
      open (unit=unit,file='g3.key',status='old',position='append')
      write (unit,10)  trim(line)
   10 format (a)
      close (unit=unit)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine sd_mdinit  --  set up dynamics in process  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "sd_mdinit" sets the bath values that the "dynamic" program
c     would provide and calls "mdinit" with its messages silenced
c
c
      subroutine sd_mdinit (temp)
      use bath
      use iounit
      implicit none
      integer isave,freeunit
      real*8 temp,dt
c
c
      kelvin = temp
      atmsph = 0.0d0
      isothermal = (temp .gt. 0.0d0)
      isobaric = .false.
      dt = 0.001d0
      isave = iout
      iout = freeunit ()
      open (unit=iout,status='scratch')
      call mdinit (dt)
      close (unit=iout)
      iout = isave
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine sd_final  --  undo an in-process setup  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "sd_final" resets the bath values and releases the system
c
c
      subroutine sd_final
      use bath
      implicit none
c
c
      kelvin = 0.0d0
      isothermal = .false.
      call final
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sd_masses  --  atomic masses of the fixture  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sd_masses" returns the atomic masses of the G3 fixture, which
c     the program-run tests need to turn saved frames into momenta
c
c
      subroutine sd_masses (nat,m)
      use atomid
      use atoms
      implicit none
      integer i,nat
      real*8 m(*)
c
c
      call pushdir ('file/stochastic')
      call loadfix ('g3','g3.key')
      do i = 1, min(n,nat)
         m(i) = mass(i)
      end do
      call popdir
      call final
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine sd_energies  --  read energies from a log  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "sd_energies" reads the potential and kinetic energy printed
c     for each saved frame of a "dynamic" run
c
c
      subroutine sd_energies (file,mx,nstep,pot,kin)
      implicit none
      integer k,mx,nstep
      integer unit,freeunit
      real*8 pot(*),kin(*)
      character*(*) file
      character*240 record
c
c
      nstep = 0
      unit = freeunit ()
      open (unit=unit,file=file,status='old',err=30)
      do while (.true.)
         read (unit,10,end=20)  record
   10    format (a240)
         if (nstep .lt. mx) then
            k = index(record,'Current Potential')
            if (k .ne. 0)  read (record(k+17:),*)  pot(nstep+1)
            k = index(record,'Current Kinetic')
            if (k .ne. 0) then
               nstep = nstep + 1
               read (record(k+15:),*)  kin(nstep)
            end if
         end if
      end do
   20 continue
      close (unit=unit)
   30 continue
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine sd_frames  --  read saved vector frames  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "sd_frames" reads the frames of a coordinate archive or a
c     velocity file; each frame is a title line followed by one line
c     per atom, and any line that is not an atom line, such as the
c     periodic box line, is skipped
c
c
      subroutine sd_frames (file,nat,mx,nfrm,vec)
      implicit none
      integer i,k,ios
      integer nat,mx,nfrm
      integer unit,freeunit
      real*8 xx,yy,zz
      real*8 vec(3,nat,*)
      character*3 name
      character*(*) file
      character*240 record
c
c
      nfrm = 0
      unit = freeunit ()
      open (unit=unit,file=file,status='old',err=30)
      do while (nfrm .lt. mx)
         read (unit,10,end=20)  record
   10    format (a240)
         i = 0
         do while (i .lt. nat)
            read (unit,10,end=20)  record
            read (record,*,iostat=ios)  k,name,xx,yy,zz
            if (ios .eq. 0) then
               i = i + 1
               vec(1,i,nfrm+1) = xx
               vec(2,i,nfrm+1) = yy
               vec(3,i,nfrm+1) = zz
            end if
         end do
         nfrm = nfrm + 1
      end do
   20 continue
      close (unit=unit)
   30 continue
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sd_momenta  --  linear and angular momentum  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sd_momenta" finds the overall linear momentum and the angular
c     momentum about the center of mass, each with a scale to judge
c     zero against; the sum of m|v| for the linear momentum, and the
c     sum of m|r-rcom||v-vcom| for the angular momentum
c
c
      subroutine sd_momenta (nat,m,xyz,vel,p,l,pscale,lscale)
      implicit none
      integer i,j,nat
      real*8 mtot,pscale,lscale
      real*8 sd_norm
      real*8 m(*),p(3),l(3)
      real*8 r(3),d(3),u(3)
      real*8 xyz(3,*),vel(3,*)
c
c
      mtot = 0.0d0
      pscale = 0.0d0
      lscale = 0.0d0
      do j = 1, 3
         r(j) = 0.0d0
         p(j) = 0.0d0
         l(j) = 0.0d0
      end do
      do i = 1, nat
         mtot = mtot + m(i)
         do j = 1, 3
            r(j) = r(j) + m(i)*xyz(j,i)
            p(j) = p(j) + m(i)*vel(j,i)
         end do
         pscale = pscale + m(i)*sd_norm(vel(1,i))
      end do
      do i = 1, nat
         do j = 1, 3
            d(j) = xyz(j,i) - r(j)/mtot
            u(j) = vel(j,i) - p(j)/mtot
         end do
         l(1) = l(1) + m(i)*(d(2)*u(3)-d(3)*u(2))
         l(2) = l(2) + m(i)*(d(3)*u(1)-d(1)*u(3))
         l(3) = l(3) + m(i)*(d(1)*u(2)-d(2)*u(1))
         lscale = lscale + m(i)*sd_norm(d)*sd_norm(u)
      end do
      return
      end
c
c
c     ################################################
c     ##                                            ##
c     ##  function sd_norm  --  length of a vector  ##
c     ##                                            ##
c     ################################################
c
c
c     "sd_norm" returns the length of a three-component vector
c
c
      real*8 function sd_norm (a)
      implicit none
      real*8 a(3)
c
c
      sd_norm = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sd_ref_trajectory  --  Tinker 8.10.5 values  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sd_ref_trajectory" returns the energies of each step and the
c     final coordinates from Fortran Tinker 8.10.5, produced by
c
c        dynamic.x g3.xyz 10 0.1 0.0001 2 298
c
c     with FRICTION 91.0 added to the keyfile; the second index is
c     one without a box and two for the periodic cell
c
c
      subroutine sd_ref_trajectory (pot,kin,xyz)
      implicit none
      integer i,j,k
      real*8 pot(10,2),kin(10,2)
      real*8 xyz(3,20,2)
      real*8 gpot(10,2),gkin(10,2)
      real*8 gxyz(3,20,2)
      data (gpot(i,1),i=1,10)  /
     &   17.6071d0, 17.4366d0, 17.2742d0, 17.1252d0, 16.9497d0,
     &   16.7513d0, 16.5434d0, 16.3543d0, 16.1759d0, 16.0185d0 /
      data (gpot(i,2),i=1,10)  /
     &   28.2011d0, 28.0317d0, 27.8706d0, 27.7226d0, 27.5481d0,
     &   27.3507d0, 27.1438d0, 26.9557d0, 26.7781d0, 26.6214d0 /
      data (gkin(i,1),i=1,10)  /
     &   15.0109d0, 15.3748d0, 16.2659d0, 16.2372d0, 16.6306d0,
     &   17.4837d0, 15.8362d0, 15.7196d0, 16.0082d0, 15.2176d0 /
      data (gkin(i,2),i=1,10)  /
     &   15.0099d0, 15.3722d0, 16.2621d0, 16.2334d0, 16.6251d0,
     &   17.4758d0, 15.8277d0, 15.7110d0, 15.9995d0, 15.2080d0 /
      data ((gxyz(j,i,1),j=1,3),i=1,20)  /
     &     2.251125d0,  -8.674469d0,  -9.320083d0,
     &     3.056947d0,  -7.238430d0, -13.687479d0,
     &     4.472119d0,  -7.221294d0, -13.394276d0,
     &     2.353982d0,  -7.439405d0, -12.359208d0,
     &     4.695021d0,  -8.336345d0, -12.398479d0,
     &     3.368089d0,  -8.412901d0, -11.617271d0,
     &     3.416223d0,  -8.135841d0, -10.116180d0,
     &     2.449207d0,  -8.793252d0,  -8.091180d0,
     &     1.297022d0,  -8.949186d0, -10.075506d0,
     &     2.776613d0,  -6.290624d0, -14.100571d0,
     &     2.726457d0,  -7.983473d0, -14.406847d0,
     &     5.060731d0,  -7.467839d0, -14.299131d0,
     &     4.810927d0,  -6.263950d0, -13.015115d0,
     &     2.247034d0,  -6.507958d0, -11.805922d0,
     &     1.388113d0,  -7.928338d0, -12.532347d0,
     &     4.875142d0,  -9.382288d0, -12.836548d0,
     &     5.596862d0,  -8.114126d0, -11.730311d0,
     &     2.999723d0,  -9.408701d0, -11.739652d0,
     &     3.443032d0,  -7.015822d0, -10.056027d0,
     &     4.379849d0,  -8.462221d0,  -9.669014d0 /
      data ((gxyz(j,i,2),j=1,3),i=1,20)  /
     &     2.251129d0,  -8.674468d0,  -9.320085d0,
     &     3.056945d0,  -7.238431d0, -13.687476d0,
     &     4.472120d0,  -7.221295d0, -13.394274d0,
     &     2.353981d0,  -7.439404d0, -12.359208d0,
     &     4.695022d0,  -8.336346d0, -12.398480d0,
     &     3.368088d0,  -8.412899d0, -11.617274d0,
     &     3.416224d0,  -8.135841d0, -10.116179d0,
     &     2.449201d0,  -8.793259d0,  -8.091161d0,
     &     1.297020d0,  -8.949187d0, -10.075505d0,
     &     2.776631d0,  -6.290613d0, -14.100654d0,
     &     2.726485d0,  -7.983442d0, -14.406907d0,
     &     5.060735d0,  -7.467815d0, -14.299181d0,
     &     4.810960d0,  -6.263924d0, -13.015227d0,
     &     2.247038d0,  -6.507959d0, -11.805928d0,
     &     1.388118d0,  -7.928335d0, -12.532348d0,
     &     4.875140d0,  -9.382286d0, -12.836548d0,
     &     5.596858d0,  -8.114121d0, -11.730314d0,
     &     2.999729d0,  -9.408699d0, -11.739657d0,
     &     3.443032d0,  -7.015821d0, -10.056028d0,
     &     4.379852d0,  -8.462220d0,  -9.669016d0 /
c
c
      do k = 1, 2
         do i = 1, 10
            pot(i,k) = gpot(i,k)
            kin(i,k) = gkin(i,k)
         end do
         do i = 1, 20
            do j = 1, 3
               xyz(j,i,k) = gxyz(j,i,k)
            end do
         end do
      end do
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine sd_ref_momentum  --  recorded momenta  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "sd_ref_momentum" returns the linear momentum and, without a
c     box, the angular momentum for each frame of the runs without
c     removal in "test_sd_momentum"; the last index of the linear
c     momentum is one without a box and two for the periodic cell
c
c
      subroutine sd_ref_momentum (p,l)
      implicit none
      integer i,j,k
      real*8 p(3,10,2),l(3,10)
      real*8 pref(3,10,2)
      real*8 lref(3,10)
      data ((pref(j,i,1),j=1,3),i=1,10)  /
     &   -21.44323141d0,  48.24403649d0, -65.67896471d0,
     &   -24.89926107d0, 104.18007632d0, -75.28894082d0,
     &   -35.82187310d0, 136.76726519d0, -54.65471061d0,
     &   -55.86937582d0, 104.34955659d0, -74.20144516d0,
     &   -86.54538138d0, 130.24691122d0, -68.49943829d0,
     &   -62.14976749d0, 102.59555184d0, -35.44686772d0,
     &   -26.23660340d0,  68.58038281d0, -21.95300828d0,
     &   -36.88670361d0,  18.44197010d0, -36.52759348d0,
     &   -10.30290253d0,   5.27234216d0, -68.30669078d0,
     &    -3.28957809d0,   0.08904002d0, -87.13163771d0 /
      data ((pref(j,i,2),j=1,3),i=1,10)  /
     &   -21.44323119d0,  48.24403659d0, -65.67896459d0,
     &   -24.89926112d0, 104.18007592d0, -75.28894095d0,
     &   -35.82187309d0, 136.76726490d0, -54.65471059d0,
     &   -55.86937579d0, 104.34955663d0, -74.20144509d0,
     &   -86.54538128d0, 130.24691122d0, -68.49943828d0,
     &   -62.14976757d0, 102.59555214d0, -35.44686799d0,
     &   -26.23660328d0,  68.58038285d0, -21.95300822d0,
     &   -36.88670357d0,  18.44197026d0, -36.52759315d0,
     &   -10.30290254d0,   5.27234231d0, -68.30669076d0,
     &    -3.28957799d0,   0.08903982d0, -87.13163787d0 /
      data ((lref(j,i),j=1,3),i=1,10)  /
     &   -340.76409588d0, -264.67332328d0, -536.57739809d0,
     &   -368.68848921d0, -292.60531404d0, -574.38037878d0,
     &   -409.61384425d0, -197.75646543d0, -550.89895397d0,
     &   -426.47365802d0, -128.02351598d0, -520.21156239d0,
     &   -337.93628604d0,  -51.85512254d0, -446.96784008d0,
     &   -327.27123923d0, -113.30665058d0, -529.98162725d0,
     &   -380.67773698d0,  -72.78805072d0, -472.17965623d0,
     &   -346.24271259d0, -149.33862985d0, -449.77628712d0,
     &   -410.05437902d0, -145.85463888d0, -485.39992627d0,
     &   -362.28080636d0,  -33.62269233d0, -394.12229391d0 /
c
c
      do i = 1, 10
         do j = 1, 3
            do k = 1, 2
               p(j,i,k) = pref(j,i,k)
            end do
            l(j,i) = lref(j,i)
         end do
      end do
      return
      end

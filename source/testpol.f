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
c     polarization, perturbation theory extrapolation, and for
c     SCF iterations in order to monitor convergence
c
c
      program testpol
      use sizes
      use atoms
      use bound
      use inform
      use iounit
      use limits
      use minima
      use polar
      use polpot
      use potent
      use rigid
      use units
      use usage
      implicit none
      integer i,j,k
      integer next,kpart
      integer nvar,iter
      integer miny
      real*8 sum,epscut
      real*8 ux,uy,uz,u2
      real*8 rdirect
      real*8 rexpt,rpart
      real*8 eps,delta
      real*8 extrap0
      real*8, allocatable :: var(:)
      real*8, allocatable :: yval(:)
      real*8, allocatable :: p(:,:)
      real*8, allocatable :: rms(:)
      real*8, allocatable :: drms(:)
      real*8, allocatable :: tdirect(:)
      real*8, allocatable :: texpt(:)
      real*8, allocatable :: tpart(:)
      real*8, allocatable :: ddirect(:,:)
      real*8, allocatable :: dexpt(:,:)
      real*8, allocatable :: dpart(:,:)
      real*8, allocatable :: udirect(:,:)
      real*8, allocatable :: uexpt(:,:)
      real*8, allocatable :: upart(:,:)
      real*8, allocatable :: ustore(:,:,:)
      logical exist,dofull
      character*1 answer
      character*120 record
      external extrap0
c
c
c     get the coordinates and required force field parameters
c
      call initial
      call getxyz
      call mechanic
c
c     check to make sure mutual polarization is being used
c
      if (.not. use_polar) then
         write (iout,10)
   10    format (/,' TESTPOL  --  Induced Dipole Polarization Model',
     &              ' is not in Use')
         call fatal
      end if
c
c     decide whether to output results by gradient component
c
      dofull = .true.
      if (n .gt. 100) then
         dofull = .false.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,20)
   20       format (/,' Output Induced Dipole Components by Atom',
     &                 ' [N] :  ',$)
            read (input,30)  record
   30       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  dofull = .true.
      end if
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     generate neighbor lists for iterative SCF solver
c
      poltyp = 'MUTUAL'
      call cutoffs
      if (use_list)  call nblist
c
c     set tolerance and rotate multipoles to global frame
c
      maxiter = 500
      epscut = poleps
      poleps = 0.0000000001d0
      debug = .false.
      call chkpole
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (rms(0:maxiter))
      allocate (drms(maxiter))
      allocate (tdirect(n))
      allocate (texpt(n))
      allocate (tpart(n))
      allocate (ddirect(3,n))
      allocate (dexpt(3,n))
      allocate (dpart(3,n))
      allocate (udirect(3,n))
      allocate (uexpt(3,n))
      allocate (upart(3,n))
      allocate (ustore(3,n,0:maxiter))
c
c     perform dynamic allocation of some global arrays
c
      allocate (uexact(3,n))
c
c     get induced dipoles for direct polarization only
c
      poltyp = 'DIRECT'
      call induce
      do i = 1, n
         do j = 1, 3
            udirect(j,i) = debye * uind(j,i)
            ustore(j,i,0) = udirect(j,i)
         end do
      end do
c
c     print the direct polarization induced dipole moments
c
      if (dofull) then
         write (iout,40)
   40    format (/,' Direct Induced Dipole Moments :',
     &           //,4x,'Atom',15x,'X',13x,'Y',13x,'Z',12x,'Norm',/)
         do i = 1, n
            if (use(i)) then
               ux = udirect(1,i)
               uy = udirect(2,i)
               uz = udirect(3,i)
               u2 = sqrt(ux*ux+uy*uy+uz*uz)
               write (iout,50)  i,ux,uy,uz,u2
   50          format (i8,4x,4f14.6)
            end if
         end do
      end if
c
c     get induced dipoles from ExPT extrapolation method
c
      poltyp = 'EXTRAP'
      call induce
      do i = 1, n
         do j = 1, 3
            uexpt(j,i) = debye * uind(j,i)
         end do
      end do
c
c     print the ExPT extrapolation induced dipole moments
c
      if (dofull) then
         write (iout,60)  cxtr(0),cxtr(1),cxtr(2),cxtr(3)
   60    format (/,' Extrapolated Induced Dipole Moments :',
     &              4x,'(',4f7.3,')',
     &           //,4x,'Atom',15x,'X',13x,'Y',13x,'Z',12x,'Norm',/)
         do i = 1, n
            if (use(i)) then
               ux = uexpt(1,i)
               uy = uexpt(2,i)
               uz = uexpt(3,i)
               u2 = sqrt(ux*ux+uy*uy+uz*uz)
               write (iout,70)  i,ux,uy,uz,u2
   70          format (i8,4x,4f14.6)
            end if
         end do
      end if
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
         sum = 0.0d0
         do i = 1, n
            do j = 1, 3
               sum = sum + (ustore(j,i,k)-ustore(j,i,k-1))**2
            end do
         end do
         drms(k) = sqrt(sum/dble(npolar))
         if (drms(k) .lt. epscut) then
            kpart = k
            epscut = -epscut
            do i = 1, n
               do j = 1, 3
                  upart(j,i) = ustore(j,i,k)
               end do
            end do
         end if
         if (drms(k) .lt. 0.5d0*poleps)  goto 80
      end do
   80 continue
      maxiter = politer
      epscut = -epscut
      do i = 1, n
         do j = 1, 3
            uexact(j,i) = ustore(j,i,maxiter)
         end do
      end do
c
c     print the partially converged and exact induced dipoles
c
      if (dofull) then
         write (iout,90)  epscut,kpart
   90    format (/,' Partially Converged SCF Induced Dipoles :',
     &              4x,'(',d9.2,' at',i3,' Iter)',
     &           //,4x,'Atom',15x,'X',13x,'Y',13x,'Z',12x,'Norm',/)
         do i = 1, n
            if (use(i)) then
               ux = upart(1,i)
               uy = upart(2,i)
               uz = upart(3,i)
               u2 = sqrt(ux*ux+uy*uy+uz*uz)
               write (iout,100)  i,ux,uy,uz,u2
  100          format (i8,4x,4f14.6)
            end if
         end do
         write (iout,110)
  110    format (/,' Exact SCF Induced Dipole Moments :',
     &           //,4x,'Atom',14x,'X',13x,'Y',13x,'Z',12x,'Norm',/)
         do i = 1, n
            if (use(i)) then
               ux = uexact(1,i)
               uy = uexact(2,i)
               uz = uexact(3,i)
               u2 = sqrt(ux*ux+uy*uy+uz*uz)
               write (iout,120)  i,ux,uy,uz,u2
  120          format (i8,4x,4f14.6)
            end if
         end do
      end if
c
c     find differences between approximate and exact dipoles
c
      rdirect = 0.0d0
      rexpt = 0.0d0
      rpart = 0.0d0
      do i = 1, n
         do j = 1, 3
            ddirect(j,i) = udirect(j,i) - uexact(j,i)
            dexpt(j,i) = uexpt(j,i) - uexact(j,i)
            dpart(j,i) = upart(j,i) - uexact(j,i)
         end do
         tdirect(i) = sqrt(ddirect(1,i)**2+ddirect(2,i)**2
     &                           +ddirect(3,i)**2)
         texpt(i) = sqrt(dexpt(1,i)**2+dexpt(2,i)**2+dexpt(3,i)**2)
         tpart(i) = sqrt(dpart(1,i)**2+dpart(2,i)**2+dpart(3,i)**2)
         rdirect = rdirect + tdirect(i)**2
         rexpt = rexpt + texpt(i)**2
         rpart = rpart + tpart(i)**2
      end do
      rdirect = sqrt(rdirect/dble(n))
      rexpt = sqrt(rexpt/dble(n))
      rpart = sqrt(rpart/dble(n))
c
c     print the RMS between approximate and exact dipoles
c
      write (iout,130)
  130 format (/,' Approximate vs. Exact Induced Dipoles :',
     &        //,4x,'Atom',14x,'Direct',13x,'ExPT',10x,'Partial SCF')
      if (dofull) then
         write (iout,140)
  140    format ()
         do i = 1, n
            if (use(i)) then
               write (iout,150)  i,tdirect(i),texpt(i),tpart(i)
  150          format (i8,4x,3f18.10)
            end if
         end do
      end if
      write (iout,160)  rdirect,rexpt,rpart
  160 format (/,5x,'RMS',4x,3f18.10)
c
c     find the RMS of each iteration from the exact dipoles
c
      do k = 0, maxiter
         sum = 0.0d0
         do i = 1, n
            do j = 1, 3
               sum = sum + (ustore(j,i,k)-uexact(j,i))**2
            end do
         end do
         rms(k) = sqrt(sum/dble(npolar))
      end do
c
c     print the RMS between iterations and versus exact dipoles
c
      write (iout,170)
  170 format (/,' Iterative Induced Dipole Convergence :',
     &        //,4x,'Iter',12x,'RMS Change',11x,'RMS vs Exact')
      write (iout,180)  0,rms(0)
  180 format (/,i8,15x,'----',6x,f20.10)
      do k = 1, maxiter
         write (iout,190)  k,drms(k),rms(k)
  190    format (i8,2x,f20.10,3x,f20.10)
         if (rms(k) .lt. 0.5d0*poleps)  goto 200
      end do
  200 continue
c
c     refine the extrapolation coefficients via optimization
c
      maxiter = 10000
      eps = 0.033d0
      delta = 0.0001d0
      write (iout,210)
  210 format (/,' Extrapolation Coefficient Refinement :',
     &        //,4x,'Iter',5x,'C0',5x,'C1',5x,'C2',5x,'C3',5x,
     &           'C4',5x,'C5',5x,'C6',5x,'C7',5x,'Residual',/)
c
c     count number of variables and define the initial simplex
c
      nvar = 0
      do i = 0, cxmax
         if (cxtr(i) .ne. 0.0d0)  nvar = nvar + 1
      end do
      allocate (var(nvar))
      allocate (yval(nvar+1))
      allocate (p(nvar+1,nvar))
      nvar = 0
      do i = 0, cxmax
         if (cxtr(i) .ne. 0.0d0) then
            nvar = nvar + 1
            var(nvar) = cxtr(i)
         end if
      end do
      do i = 1, nvar
         p(1,i) = var(i)
      end do
      yval(1) = extrap0 (var)
      do k = 1, nvar
         var(k) = var(k) + eps
         do i = 1, nvar
            p(k+1,i) = var(i)
         end do
         yval(k+1) = extrap0 (var)
         var(k) = var(k) - eps
      end do
c
c     optimize coefficients, then find and print refined values
c
      call simplex (nvar,p,yval,delta,extrap0,iter)
      iter = iter + nvar + 1
      rexpt = 1000000.0d0
      do i = 1, nvar+1
         if (yval(i) .lt. rexpt) then
            miny = i
            rexpt = yval(i)
         end if
      end do
      nvar = 0
      do i = 0, cxmax
         if (cxtr(i) .ne. 0.0d0) then
            nvar = nvar + 1
            cxtr(i) = p(miny,nvar)
         end if
      end do
      write (iout,220)  iter,(cxtr(i),i=0,maxxtr),rexpt
  220 format (i8,1x,8f7.3,f14.10)
c
c     perform deallocation of some local arrays
c
      deallocate (var)
      deallocate (yval)
      deallocate (p)
      deallocate (rms)
      deallocate (drms)
      deallocate (tdirect)
      deallocate (texpt)
      deallocate (tpart)
      deallocate (ddirect)
      deallocate (dexpt)
      deallocate (dpart)
      deallocate (udirect)
      deallocate (uexpt)
      deallocate (upart)
      deallocate (ustore)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function extrap0  --  extrapolation coefficient refinement  ##
c     ##                                                              ##
c     ##################################################################
c
c
      function extrap0 (var)
      use sizes
      use atoms
      use iounit
      use polar
      use polpot
      use units
      implicit none
      integer i,j
      integer iter,nvar
      real*8 extrap0
      real*8 rexpt
      real*8 var(*)
      real*8, allocatable :: uexpt(:,:)
      logical first
      save first,iter
      data first  / .true. /
c
c
c     count the number of times the function has been called
c
      if (first) then
         first = .false.
         iter = -1
      end if
      iter = iter + 1
c
c     copy optimization variables into extrapolation coefficients
c
      nvar = 0
      do i = 0, maxxtr
         if (cxtr(i) .ne. 0.0d0) then
            nvar = nvar + 1
            cxtr(i) = var(nvar)
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (uexpt(3,n))
c
c     compute RMS error between ExPT and exact dipoles
c
      poltyp = 'EXTRAP'
      call induce
      do i = 1, n
         do j = 1, 3
            uexpt(j,i) = debye * uind(j,i)
         end do
      end do
      rexpt = 0.0d0
      do i = 1, n
         do j = 1, 3
c           rexpt = rexpt + (uexpt(j,i)-uexact(j,i))**2
            rexpt = rexpt + (uexpt(j,i)-uexact(j,i))**6
         end do
      end do
      rexpt = sqrt(rexpt/dble(n))
      if (mod(iter,10) .eq. 0) then
         write (iout,10)  iter,(cxtr(i),i=0,maxxtr),rexpt
   10    format (i8,1x,8f7.3,f14.10)
      end if
c
c     set the return value equal to the RMS error
c
      extrap0 = rexpt
c
c     perform deallocation of some local arrays
c
      deallocate (uexpt)
      return
      end

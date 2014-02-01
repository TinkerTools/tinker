c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1998 by Rohit Pappu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program agda  --  adiabatic gaussian density annealing  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "agda" implements the Adiabatic Gaussian Density Annealing
c     method (AGDA) for global optimization using a conjugate
c     gradient optimization on differently annealed potential
c     surfaces and a numerical integrator to control the widths
c     of the Gaussian densities
c
c     literature references :
c
c     J. Ma and J. E. Straub, "Simulated Annealing using the
c     Classical Density Distribution", Journal of Chemical Physics,
c     101, 533-541 (1994)
c
c     C. Tsoo and C. L. Brooks III "Cluster Structure Determintation
c     using Gaussian Density Distribution Global Minimization
c     Methods", Journal of Chemical Physics, 101, 6405-6411 (1994)
c
c     J. Kostrowicki and H. A. Scheraga, "Application of the Diffusion
c     Equation Method for Global Optimization to Oligopeptides",
c     Journal of Physical Chemistry, 96, 7442-7449 (1992)
c
c
      program agda
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'iounit.i'
      include 'minima.i'
      include 'potent.i'
      include 'warp.i'
      integer maxgda
      parameter (maxgda=4*maxatm)
      integer i,igda,lext,freeunit
      integer itrial,ntrial,nstep
      integer nvar,nok,nbad
      integer im2
      real*8 bstart,bstop
      real*8 random,boxsize
      real*8 eps,h1,hmin
      real*8 minimum,grdmin
      real*8 xcm,ycm,zcm
      real*8 gda1,gda2
      real*8 m2init(maxatm)
      real*8 xx(maxgda)
      logical exist,randomize
      character*1 answer
      character*6 mode,method
      character*7 ext,status
      character*60 gdafile
      character*60 m2file
      character*80 record
      character*80 string
      external gda1,gda2,gda3,optsave
c
c
c     initialize and get the structure to be optimized
c
      call initial
      call getxyz
c
c     set flag for GDA method, then set up molecular mechanics
c
      use_smooth = .true.
      use_gda = .true.
      call mechanic
c
c     store the initial values of the squared mean Gaussian width
c
      do i = 1, n
         m2init(i) = m2(1)
      end do
c
c     get the number of optimized structures to be constructed
c
      ntrial = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  ntrial
   10 continue
      if (ntrial .le. 0) then
         write (iout,20)
   20    format (/,' Enter Number of Annealing Trials [1] :  ',$)
         read (input,30)  ntrial
   30    format (i10)
      end if
      if (ntrial .le. 0)  ntrial = 1
c
c     see if random coordinates are desired as starting structures
c
      randomize = .true.
      if (ntrial .eq. 1) then
         randomize = .false.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,40)
   40       format (/,' Use Randomized Initial Coordinates [N] :  ',$)
            read (input,50)  answer
   50       format (a1)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  randomize = .true.
      end if
      if (randomize)  boxsize = 10.0d0 * (dble(n))**(1.0d0/3.0d0)
c
c     get initial and final values of inverse temperature
c
      bstart = -1.0d0
      bstop = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  bstart
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  bstop
   60 continue
      if (bstart.le.0.0d0 .or. bstop.le.0.0d0) then
         write (iout,70)
   70    format (/,' Enter Initial and Final Beta [0.01, 10**10] :  ',$)
         read (input,80)  record
   80    format (a80)
         read (record,*,err=90,end=90)  bstart,bstop
   90    continue
      end if
      if (bstart .le. 0.0d0)  bstart = 0.01d0
      if (bstop .le. 0.0d0)  bstop = 1.0d10
c
c     write out a copy of coordinates for later update
c
      do itrial = 1, ntrial
         lext = 3
         call numeral (itrial,ext,lext)
         gdafile = filename(1:leng)//'.'//ext(1:lext)
         call version (gdafile,'new')
         igda = freeunit ()
         open (unit=igda,file=gdafile,status='new')
         call prtxyz (igda)
         close (unit=igda)
         outfile = gdafile
c
c     set an initial box size and generate random coordinates
c
         if (randomize) then
            do i = 1, n
               x(i) = boxsize * random ()
               y(i) = boxsize * random ()
               z(i) = boxsize * random ()
            end do
         end if
c
c     translate M2's to parameters for numerical integration
c
         nvar = 0
         do i = 1, n
            nvar = nvar + 1
            xx(nvar) = m2init(i)
         end do
c
c        make changes to the potential to use potential smoothing
c
         use_smooth = .true.
         use_geom = .true.
         use_gda = .true.
c
c        make the call to the Bulirsch-Stoer integration routine
c
         nstep = 0
         status = '       '
         eps = 1.0d-8
         h1 = 0.1d0
         hmin = 0.0d0
         write (iout,100)
c
         m2file = 'm2.'//ext(1:lext)
         call version (m2file,'new')
         im2 = freeunit ()
         open (unit=im2,file=m2file,status='new')
  100    format (//,' Adiabatic Gaussian Density Annealing',
     &              ' Optimization :',
     &           //,' BS Step',5x,'Log(Beta)',6x,'Energy',
     &              9x,'Rg',8x,'Log(M2)',7x,'Status',/)
c
         call gdastat (im2,nstep,bstart,xx,status)
c        call diffeq (im2,nvar,xx,bstart,bstop,eps,h1,hmin,nok,
c    &                nbad,gda1)
         call diffeq (nvar,xx,bstart,bstop,eps,h1,hmin,nok,nbad,gda1)
         nstep = nok + nbad
c
c        make changes to the potential for standard optimization
c
         use_smooth = .false.
         use_gda = .false.
         use_geom = .false.
c
c        make the call to the energy minimization routine
c
         nvar = 0
         do i = 1,n
            nvar = nvar + 1
            xx(nvar) = x(i)
            nvar = nvar + 1
            xx(nvar) = y(i)
            nvar = nvar + 1
            xx(nvar) = z(i)
         end do
         mode = 'DTNCG'
         method = 'AUTO'
         grdmin = 0.0001d0
         nextiter = nstep + 1
         call tncg (mode,method,nvar,xx,minimum,grdmin,
     &                     gda2,gda3,optsave)
         write (iout,110)  itrial,minimum
  110    format (/,' Global Energy Minimum for Trial',i4,' :',f15.4)
c
c     translate optimization parameters into atomic coordinates
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
c     move the center of mass to the origin
c
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         do i = 1, n
            xcm = xcm + x(i)
            ycm = ycm + y(i)
            zcm = zcm + z(i)
         end do
         xcm = xcm / dble(n)
         ycm = ycm / dble(n)
         zcm = zcm / dble(n)
         do i = 1, n
            x(i) = x(i) - xcm
            y(i) = y(i) - xcm
            z(i) = z(i) - xcm
         end do
c
c     write the final coordinates into a file
c
         close(unit=im2)
         igda = freeunit ()
         open (unit=igda,file=gdafile,status='old')
         rewind (unit=igda)
         call prtxyz (igda)
         close (unit=igda)
      end do
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function gda1  --  gradient for gaussian density annealing  ##
c     ##                                                              ##
c     ##################################################################
c
c
      function gda1 (beta,xx,g)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'warp.i'
      integer maxgda
      parameter (maxgda=4*maxatm)
      integer i,nvar
      integer hinit(maxvar),hstop(maxvar)
      integer hindex(maxhess)
      real*8 energy
      real*8 gda1,beta,sum
      real*8 hdiag(maxvar),h(maxhess)
      real*8 xx(maxgda),g(maxgda)
c
c
c     translate optimization parameters to M2's
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         m2(i) = xx(nvar)
         if (m2(i) .lt. 0.0d0) then
            write (iout,10)  i,m2(i)
   10       format (' GDA1  --  Warning, Negative M2 at Atom',i6,
     &                  ' with Value',f14.4)
            m2(i) = -m2(i)
         end if
      end do
c
c     compute and store the energy and gradient
c
      gda1 = energy ()
c
c
c     compute and store the Hessian elements
c
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     translate the Hessian diagonal into a dM2/dbeta vector
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         sum = hdiag(3*i-2) + hdiag(3*i-1) + hdiag(3*i)
         g(nvar) = -(m2(i)/3.0d0)**2 * sum
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function gda2  --  energy/gradient for TNCG optimization  ##
c     ##                                                            ##
c     ################################################################
c
c
      function gda2 (xx,g)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,nvar
      real*8 gda2,e,derivs(3,maxatm)
      real*8 xx(maxvar),g(maxvar)
c
c
c     translate optimization parameters to atomic coordinates
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
c     compute and store the energy and gradient
c
      call gradient (e,derivs)
      gda2 = e
c
c     store atom gradients as optimization gradient, also
c     translate the coordinates of each active atom; the
c     latter may be needed when using periodic boundaries
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         xx(nvar) = x(i)
         g(nvar) = derivs(1,i)
         nvar = nvar + 1
         xx(nvar) = y(i)
         g(nvar) = derivs(2,i)
         nvar = nvar + 1
         xx(nvar) = z(i)
         g(nvar) = derivs(3,i)
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine gda3  --  Hessian values for TNCG optimization  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine gda3 (mode,xx,h,hinit,hstop,hindex,hdiag)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,nvar
      integer hinit(maxvar),hstop(maxvar)
      integer hindex(maxhess)
      real*8 xx(maxvar),hdiag(maxvar),h(maxhess)
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
c
c     translate the coordinates of each active atom;
c     this may be needed when using periodic boundaries
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
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine gdastat  --  compute GDA trajectory averages  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine gdastat (munit,nstep,beta,xx,status)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'math.i'
      include 'warp.i'
      include 'minima.i'
      include 'inform.i'
      integer maxgda,munit
      parameter (maxgda=4*maxatm)
      integer i,nstep,nvar,nmin,oldprt
      integer neigen,iter
      real*8 beta,xx(maxgda),xc(maxvar)
      real*8 e,energy,rg,m2ave
      real*8 grdmin,minimum
      real*8 mean,std,dev
      logical oldverb
      logical done,use_local
      character*6 mode,method
      character*7 status
      external gda2,gda3,optsave
c
c
c     translate optimization parameters to M2's
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         m2(i) = abs(xx(nvar))
      end do
c
c     translate coordinates to optimization parameters
c
      nmin = 0
      do i = 1, n
         nmin = nmin + 1
         xc(nmin) = x(i)
         nmin = nmin + 1
         xc(nmin) = y(i)
         nmin = nmin + 1
         xc(nmin) = z(i)
      end do
c
c     make the call to the energy minimization routine
c
      mode = 'DTNCG'
      method = 'AUTO'
      grdmin = 0.0001d0
      oldverb = verbose
      oldprt = iprint
c     verbose = .false.
      verbose = .true.
      iprint = 0
c     call lbfgs (nmin,xc,minimum,grdmin,gda2,optsave)
      call tncg (mode,method,nmin,xc,minimum,grdmin,
     &                     gda2,gda3,optsave)
c
c     translate optimization parameters to coordinates
c
      nmin = 0
      do i = 1, n
         nmin = nmin + 1
         x(i) = xc(nmin)
         nmin = nmin + 1
         y(i) = xc(nmin)
         nmin = nmin + 1
         z(i) = xc(nmin)
      end do
c
      e = energy ()
c
      neigen = 5
      call stat (n,xx,mean,std)
      use_local = .false.
c     if (mean .gt. 0.5d0) then
c        done = .false.
c        iter = 0
c        dowhile (.not. done)
c           iter = iter + 1
c           if (iter .gt. n) then
c              use_local = .false.
c              done = .true.
c           else
c              dev = abs(m2(iter) - mean)
c              if (dev .gt. 2.0d0*std) then
c                 use_local = .true.
c                 done = .true.
c              else
c                 use_local = .false.
c                 done = .false.
c              end if
c           end if
c        end do
c     end if
c     if (use_local)  call loclsrch (neigen,minimum)
      call stat (n,xx,mean,std)
      write (munit,10)  mean,std,log(beta)/logten,(m2(i),i=1,n)
   10 format(22f12.4)
      call gyrate (rg)
      m2ave = 0.0d0
      do i = 1, n
         m2ave = m2ave + m2(i)
      end do
      m2ave = m2ave / dble(n)
c
c
      write (iout,20)  nstep,log(beta)/logten,e,rg,
     &                 log(m2ave)/logten,status
   20 format (i6,2x,4f13.4,6x,a7)
c
c     save the current coordinates to a disk file
c
      nmin = 0
      do i = 1, n
         nmin = nmin + 1
         xc(nmin)=x(i)
         nmin = nmin + 1
         xc(nmin)=y(i)
         nmin = nmin + 1
         xc(nmin)=z(i)
      end do
      call optsave (nstep,e,xc)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine loclsrch  --  local search on smoothed surface  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine loclsrch (neigen,minimum)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'omega.i'
      include 'refer.i'
      integer i,j,k,neigen,ndoi,nsearch
      real*8 minimum,minref,minbest
      real*8 eps,rms,step(3,maxvib)
      real*8 eigen(maxvib),vects(maxvib,maxvib)
      real*8 xbest(maxatm),ybest(maxatm),zbest(maxatm)
      logical done
c
c
c     store the current coordinates as the reference set
c
      call makeref
c
c     set parameters related to the local search procedure
c
      done = .false.
      eps = 1.0d-4
      minref = minimum
      minbest = minimum
      ndoi = 0
c
c     find local minimum along each of the steepest directions
c
      dowhile (.not. done)
         ndoi = ndoi + 1
         write (iout,10)  ndoi,minref
   10    format (/,' Local Search :',10x,'Iteration',i4,
     &              5x,'Energy',f12.4,/)
         call eigencart (eigen,vects)
c
c     search both directions along each eigenvector in turn
c
         nsearch = 0
         do i = 1, neigen
            do k = 1, n
               j = 3*(k-1)
               step(1,k) = vects(j+1,n-i+1)
               step(2,k) = vects(j+2,n-i+1)
               step(3,k) = vects(j+3,n-i+1)
            end do
            nsearch = nsearch + 1
            call getref
            call climber (nsearch,minimum,step)
            if (minimum .lt. minbest) then
               minbest = minimum
               do k = 1, n
                  xbest(k) = x(k)
                  ybest(k) = y(k)
                  zbest(k) = z(k)
               end do
            end if
            do k = 1, n
               j = 3*(k-1)
               step(1,k) = -vects(j+1,n-i+1)
               step(2,k) = -vects(j+2,n-i+1)
               step(3,k) = -vects(j+3,n-i+1)
            end do
            nsearch = nsearch + 1
            call getref
            call climber (nsearch,minimum,step)
            if (minimum .lt. minbest) then
               minbest = minimum
               do k = 1, n
                  xbest(k) = x(k)
                  ybest(k) = y(k)
                  zbest(k) = z(k)
               end do
            end if
         end do
c
c     check for convergence of the local search procedure
c
         if (minbest .lt. minref-eps) then
            done = .false.
            minref = minbest
            call impose (n,xref,yref,zref,n,xbest,ybest,zbest,rms)
            do k = 1, n
               x(k) = xbest(k)
               y(k) = ybest(k)
               z(k) = zbest(k)
            end do
            call makeref
         else
            done = .true.
            minimum = minref
            call getref
         end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eigencart  --  Cartesian Hessian eigenvectors  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine eigencart (eigen,vects)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'hescut.i'
      include 'iounit.i'
      include 'math.i'
      integer i,j,k,ihess
      integer nfreq,ndummy,nvib
      integer hinit(3,maxatm),hstop(3,maxatm),hindex(maxhess)
      real*8 h(maxhess),hdiag(3,maxatm)
      real*8 matrix((maxvib+1)*maxvib/2)
      real*8 eigen(maxvib),vects(maxvib,maxvib)
      real*8 a(maxvib+1),b(maxvib+1),p(maxvib+1),w(maxvib+1)
      real*8 ta(maxvib+1),tb(maxvib+1),ty(maxvib+1)
c
      nfreq = 3 * n
      ndummy = 0
      nvib = nfreq - 3*ndummy
c     calculate the Hessian matrix of second derivatives
c
      hesscut = 0.0d0
      call hessian (h,hinit,hstop,hindex,hdiag)
c
c     store upper triangle of the Hessian in "matrix"
c
      ihess = 0
      do i = 1, n
         do j = 1, 3
            ihess = ihess + 1
            matrix(ihess) = hdiag(j,i)
            do k = hinit(j,i), hstop(j,i)
               ihess = ihess + 1
               matrix(ihess) = h(k)
            end do
         end do
      end do
c
c     perform diagonalization to get Hessian eigenvalues
c
      call diagq (nfreq,maxvib,nfreq,matrix,eigen,vects,
     &                     a,b,p,w,ta,tb,ty)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine climber  --  find minimum along torsional mode  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine climber (nsearch,minimum,step)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'refer.i'
      include 'iounit.i'
      include 'math.i'
      integer maxstep
      parameter (maxstep=500)
      integer i,kstep,nstep,nsearch
      real*8 minimum,big,energy
      real*8 step(3,maxvib),estep(0:maxstep)
      logical done
c
c
c     convert current reference coordinates to a Z-matrix
c
      call getref
c
c     set the maximum number of steps and the step size
c
      done = .false.
      big = 100000.0d0
      minimum = big
      kstep = 0
      nstep = 65
c
c     scan the search direction for a minimization candidate
c
      dowhile (.not. done)
         if (kstep .ne. 0) then
            do i = 1, n
               x(i) = x(i) + step(1,i)
               y(i) = y(i) + step(2,i)
               z(i) = z(i) + step(3,i)
            end do
         end if
         estep(kstep) = energy ()
         if (kstep .ge. 2) then
            if (estep(kstep) .lt. estep(kstep-2) .and.
     &          estep(kstep-1) .lt. estep(kstep-2)) then
               done = .true.
               do i = 1, n
                  x(i) = x(i) - step(1,i)
                  y(i) = y(i) - step(2,i)
                  z(i) = z(i) - step(3,i)
               end do
               call loclmin (minimum)
               if (minimum .ge. -big) then
                  write (iout,10)  nsearch,kstep+1,minimum
   10             format (4x,'Search Direction',i4,10x,'Step',
     &                       i4,11x,f12.4)
               else
                  minimum = big
                  write (iout,20)  nsearch
   20             format (4x,'Search Direction',i4,35x,'------')
               end if
            end if
         end if
         if (kstep.ge.nstep .and. .not.done) then
            done = .true.
            write (iout,30)  nsearch
   30       format (4x,'Search Direction',i4,35x,'------')
         end if
         kstep = kstep + 1
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine loclmin  --  optimization for Doi local search  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "loclmin" performs an energy minimization in Cartesian
c     coordinate space using a truncated Newton method
c
c
      subroutine loclmin (minimum)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'minima.i'
      integer i,nvar,oldprt
      real*8 minimum,grdmin
      real*8 gda2,xx(maxvar)
      logical oldverb
      character*6 mode,method
      external gda2,gda3,optsave
c
c
c     translate the coordinates of each atom
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
c     make the call to the optimization routine
c
      grdmin = 0.0001d0
      oldverb = verbose
      oldprt = iprint
      verbose = .false.
      iprint = 0
      call lbfgs (nvar,xx,minimum,grdmin,gda2,optsave)
c
c     untranslate the final coordinates for each atom
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
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine stat  --  mean and standard deviation  ##
c     ##                                                    ##
c     ########################################################
c
c
      subroutine stat (nz,xx,mean,std)
      implicit none
      include 'sizes.i'
      integer maxgda
      parameter (maxgda=maxatm)
      integer i,nz
      real*8 mean,std
      real*8 var,dev
      real*8 xx(maxgda)
c
c
      mean = 0.0d0
      do i = 1, nz
        mean = mean + xx(i)
      end do
      mean = mean / nz
      var = 0.0d0
      do i = 1, nz
        dev = xx(i) - mean
        var = var + dev*dev
      end do
      var = var / dble(nz-1)
      std = sqrt(var)
      return
      end

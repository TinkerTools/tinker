c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2013  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program bar  --  free energy perturbation via BAR method  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "bar" computes the free energy difference between two states
c     via the Bennett acceptance ratio (BAR) method
c
c     current version takes as input the coordinate archive and key
c     file for state A, and similar for state B; it then computes the
c     total potential energy for all frames of each state under the
c     control of both key files; finally the BAR equation is iterated
c     to compute the free energy difference for state A --> state B
c     and a bootstrap estimate of the standard deviation
c
c     literature reference:
c
c     C. H. Bennett, "Efficient Estimation of Free Energy Differences
c     from Monte Carlo Data", Journal of Computational Physics, 22,
c     245-268 (1976)
c
c
      program bar
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'units.i'
      integer i,j,k
      integer ixyz,nbst
      integer iter,maxiter
      integer nkey1,nfrm1
      integer nkey2,nfrm2
      integer maxframe
      integer freeunit
      integer, allocatable :: bst1(:)
      integer, allocatable :: bst2(:)
      real*8 energy
      real*8 temp,rt
      real*8 delta,eps
      real*8 cold,cnew
      real*8 top,bot,rfrm
      real*8 sum,sum2
      real*8 mean,stdev
      real*8 random,ratio
      real*8, allocatable :: ua0(:)
      real*8, allocatable :: ua1(:)
      real*8, allocatable :: ub0(:)
      real*8, allocatable :: ub1(:)
      logical done
      character*120 fname1
      character*120 fname2
      character*120 xyzfile
      character*120, allocatable :: keys1(:)
      character*120, allocatable :: keys2(:)
c
c     set up trajectory file A and the mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys1(nkey))
c
c     store the filename and keyword file for trajectory A
c
      fname1 = filename
      nkey1 = nkey
      do i = 1, nkey1
         keys1(i) = keyline(i)
      end do
c
c     set up trajectory file B and the mechanics calculation
c
      call getxyz
      call mechanic
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys2(nkey))
c
c     store the filename and keyword file for trajectory B
c
      fname2 = filename
      nkey2 = nkey
      do i = 1, nkey2
         keys2(i) = keyline(i)
      end do
c
c     perform dynamic allocation of some local arrays
c
      maxframe = 20000
      allocate (ua0(maxframe))
      allocate (ua1(maxframe))
      allocate (ub0(maxframe))
      allocate (ub1(maxframe))
c
c     reopen trajectory file A and read the initial structure
c
      ixyz = freeunit ()
      xyzfile = fname1
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     find the energies for each successive trajectory frame
c
      write (iout,10)
   10 format (/,' Potential Energy Values for Trajectory A :',
     &        //,7x,'Frame',11x,'Key A',11x,'Key B',12x,'Delta',/)
      k = 0
      do while (.not. abort)
         k = k + 1
         nkey = nkey1
         do i = 1, nkey
            keyline(i) = keys1(i)
         end do
c        call mechanic
         call kmpole
         call kpolar
         call mutate
         ua0(k) = energy ()
         nkey = nkey2
         do i = 1, nkey
            keyline(i) = keys2(i)
         end do
c        call mechanic
         call kmpole
         call kpolar
         call mutate
         ua1(k) = energy ()
         write (iout,20)  k,ua0(k),ua1(k),ua1(k)-ua0(k)
   20    format (i11,2x,3f16.4)
         call readxyz (ixyz)
      end do
      nfrm1 = k
      close (unit=ixyz)
c
c     reopen trajectory file B and read the initial structure
c
      ixyz = freeunit ()
      xyzfile = fname2
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     find the energies for each successive trajectory frame
c
      write (iout,30)
   30 format (/,' Potential Energy Values for Trajectory B :',
     &        //,7x,'Frame',11x,'Key A',11x,'Key B',12x,'Delta',/)
      k = 0
      do while (.not. abort)
         k = k + 1
         nkey = nkey1
         do i = 1, nkey
            keyline(i) = keys1(i)
         end do
c        call mechanic
         call kmpole
         call kpolar
         call mutate
         ub0(k) = energy ()
         nkey = nkey2
         do i = 1, nkey
            keyline(i) = keys2(i)
         end do
c        call mechanic
         call kmpole
         call kpolar
         call mutate
         ub1(k) = energy ()
         write (iout,40)  k,ub0(k),ub1(k),ub0(k)-ub1(k)
   40    format (i11,2x,3f16.4)
         call readxyz (ixyz)
      end do
      nfrm2 = k
      close (unit=ixyz)
c
c     perform deallocation of some local arrays
c
      deallocate (keys1)
      deallocate (keys2)
c
c     set various constants and parameter values for BAR
c
      write (iout,50)
   50 format (/,' Estimation of Free Energy Difference',
     &           ' via BAR Method :',/)
      maxiter = 100
      eps = 0.0001d0
      rfrm = dble(nfrm1) / dble(nfrm2)
      temp = 298.0d0
      rt = gasconst * temp
c
c     compute the initial free energy via the BAR equation
c
      done = .false.
      iter = 0
      cold = 0.0d0
      top = 0.0d0
      do i = 1, nfrm2
         top = top + 1.0d0/(1.0d0+exp((ub0(i)-ub1(i)+cold)/rt)) 
      end do
      bot = 0.0d0
      do i = 1, nfrm1
         bot = bot + 1.0d0/(1.0d0+exp((ua1(i)-ua0(i)-cold)/rt)) 
      end do
      cnew = rt*log(rfrm*top/bot) + cold
      delta = abs(cnew-cold)
      if (delta .lt. eps) then
         done = .true.
         write (iout,60)  cnew,iter
   60    format (' BAR Free Energy Estimate',1x,f12.4,
     &              ' Kcal/mol at',i4,' Iterations')  
      else
         write (iout,70)  iter,cnew
   70    format (' BAR Iteration',i4,8x,f12.4,' Kcal/mol')  
      end if
c
c     iterate the BAR equation to converge the free energy
c
      do while (.not. done)
         iter = iter + 1
         cold = cnew
         top = 0.0d0
         do i = 1, nfrm2
            top = top + 1.0d0/(1.0d0+exp((ub0(i)-ub1(i)+cold)/rt)) 
         end do
         bot = 0.0d0
         do i = 1, nfrm1
            bot = bot + 1.0d0/(1.0d0+exp((ua1(i)-ua0(i)-cold)/rt)) 
         end do
         cnew = rt*log(rfrm*top/bot) + cold
         delta = abs(cnew-cold)
         if (delta .lt. eps) then
            done = .true.
            write (iout,80)  cnew,iter
   80       format (/,' BAR Free Energy Estimate',1x,f12.4,
     &                 ' Kcal/mol at',i4,' Iterations')  
         else
            write (iout,90)  iter,cnew
   90       format (' BAR Iteration',i4,8x,f12.4,' Kcal/mol')  
         end if
         if (iter.ge.maxiter .and. .not.done) then
            done = .true.
            write (iout,100)  maxiter
  100       format (/,' BAR Free Energy Estimate not Converged',
     &                 ' after',i4,' Iterations')
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      nbst = 50000
      allocate (bst1(nbst))
      allocate (bst2(nbst))
c
c     use bootstrapping analysis to estimate statistical error
c
      write (iout,110)
  110 format (/,' Bootstrapping Error Analysis of BAR',
     &           ' Free Energy :',/)
      sum = 0.0d0
      sum2 = 0.0d0
      do k = 1, nbst
         done = .false.
         iter = 0
         cold = 0.0d0
         top = 0.0d0
         do i = 1, nfrm2
            bst2(i) = int(dble(nfrm2)*random()) + 1
            j = bst2(i)
            top = top + 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)+cold)/rt)) 
         end do
         bot = 0.0d0
         do i = 1, nfrm1
            bst1(i) = int(dble(nfrm1)*random()) + 1
            j = bst1(i)
            bot = bot + 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)-cold)/rt)) 
         end do
         cnew = rt*log(rfrm*top/bot) + cold
         delta = abs(cnew-cold)
         do while (.not. done)
            iter = iter + 1
            cold = cnew
            top = 0.0d0
            do i = 1, nfrm2
               j = bst2(i)
               top = top + 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)+cold)/rt)) 
            end do
            bot = 0.0d0
            do i = 1, nfrm1
               j = bst1(i)
               bot = bot + 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)-cold)/rt)) 
            end do
            cnew = rt*log(rfrm*top/bot) + cold
            delta = abs(cnew-cold)
            if (delta .lt. eps) then
               done = .true.
               sum = sum + cnew
               sum2 = sum2 + cnew*cnew
               if (debug) then
                  write (iout,120)  k,cnew,iter
  120             format (' Bootstrap Estimate',i6,1x,f12.4,
     &                       ' Kcal/mol at',i4,' Iterations')
               end if  
            end if
         end do
      end do
      mean = sum / dble(nbst)
      ratio = dble(nbst/(nbst-1))
      stdev = sqrt(ratio*(sum2/dble(nbst)-mean*mean))
      if (verbose) then
         write (iout,130)
  130    format ()
      end if
      write (iout,140)  mean,stdev
  140 format (' Bootstrap Free Energy',4x,f12.4,
     &           ' +/-',f8.4,' Kcal/mol')  
c
c     perform deallocation of some local arrays
c
      deallocate (ua0)
      deallocate (ua1)
      deallocate (ub0)
      deallocate (ub1)
      deallocate (bst1)
      deallocate (bst2)
c
c     perform any final tasks before program exit
c
      call final
      end

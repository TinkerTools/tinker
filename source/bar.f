c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2013  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program bar  --  free energy differences via FEP and BAR  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "bar" computes the free energy difference between two states
c     via the Zwanzig free energy perturbation (FEP) and Bennett
c     acceptance ratio (BAR) methods, while also decomposing free
c     energy into its enthalpy and entropy components
c
c     note current version takes as input the trajectory archives A
c     and B, originally generated for states 0 and 1, respectively;
c     then finds the total potential energy for all frames of each
c     trajectory under control of key files for both states 0 and 1;
c     finally the FEP and BAR algorithms are used to compute the
c     free energy for state 0 --> state 1, and similar estimators
c     provide the enthalpy and entropy for state 0 --> state 1
c
c     modifications to handle NPT simulation data by Chengwen Liu,
c     University of Texas at Austin during October 2015; addition of
c     enthalpy and entropy estimates by Aaron Gordon, December 2016
c
c     literature references:
c
c     C. H. Bennett, "Efficient Estimation of Free Energy Differences
c     from Monte Carlo Data", Journal of Computational Physics, 22,
c     245-268 (1976)
c
c     K. B. Daly, J. B. Benziger, P. G. Debenedetti and
c     A. Z. Panagiotopoulos, "Massively Parallel Chemical Potential
c     Calculation on Graphics Processing Units", Computer Physics
c     Communications, 183, 2054-2062 (2012)  [modification for NPT]
c
c     M. A. Wyczalkowski, A. Vitalis and R. V. Pappu, "New Estimators
c     for Calculating Solvation Entropy and Enthalpy and Comparative
c     Assessments of Their Accuracy and Precision, Journal of Physical
c     Chemistry, 114, 8166-8180 (2010)  [entropy and enthalpy]
c
c
      program bar
      use sizes
      use atoms
      use boxes
      use energi
      use files
      use inform
      use iounit
      use keys
      use titles
      use units
      implicit none
      integer i,j,k
      integer n0,n1
      integer ixyz,ibar,nbst
      integer iter,maxiter
      integer leng0,leng1
      integer ltitle0,ltitle1
      integer nkey0,nkey1
      integer nfrma,nfrmb
      integer starta,startb
      integer stopa,stopb
      integer stepa,stepb
      integer maxframe
      integer freeunit
      integer trimtext
      integer, allocatable :: bsta(:)
      integer, allocatable :: bstb(:)
      real*8 energy,term
      real*8 rt,rta,rtb
      real*8 delta,eps
      real*8 frma,frmb
      real*8 tempa,tempb
      real*8 cold,cnew
      real*8 top,top2
      real*8 bot,bot2
      real*8 fterm,rfrm
      real*8 sum,sum2,vave
      real*8 mean,stdev
      real*8 random,ratio
      real*8 cfore,cback
      real*8 uave0,uave1
      real*8 hfore,hback
      real*8 fore,back
      real*8 hdir,hbar
      real*8 sbar,tsbar
      real*8 numer,denom
      real*8 fsum,bsum
      real*8 fvsum,bvsum
      real*8 fbvsum,vsum
      real*8 fbsum0,fbsum1
      real*8 alpha0,alpha1
      real*8, allocatable :: ua0(:)
      real*8, allocatable :: ua1(:)
      real*8, allocatable :: ub0(:)
      real*8, allocatable :: ub1(:)
      real*8, allocatable :: vola(:)
      real*8, allocatable :: volb(:)
      logical exist,query,done
      character*120 record
      character*120 string
      character*120 fname0
      character*120 fname1
      character*120 title0
      character*120 title1
      character*120 xyzfile
      character*120 barfile
      character*120, allocatable :: keys0(:)
      character*120, allocatable :: keys1(:)
c
c
c     perform dynamic allocation of some local arrays
c
      maxframe = 1000000
      allocate (ua0(maxframe))
      allocate (ua1(maxframe))
      allocate (ub0(maxframe))
      allocate (ub1(maxframe))
      allocate (vola(maxframe))
      allocate (volb(maxframe))
c
c     get trajectory A archive and setup mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     set beginning and ending frame for trajectory A
c
      starta = 0
      stopa = 0
      stepa = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  starta
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  stopa
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  stepa
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Numbers of First & Last File and Step',
     &              ' Increment :  ',$)
         read (input,30)  record
   30    format (a120)
         read (record,*,err=40,end=40)  starta,stopa,stepa
   40    continue
      end if
      if (starta .eq. 0)  starta = 1
      if (stopa .eq. 0)  stopa = maxframe
      if (stepa .eq. 0)  stepa = 1
c
c     find the original temperature value for trajectory A
c
      tempa = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=50,end=50)  tempa
   50 continue
      do while (tempa .lt. 0.0d0)
         write (iout,60)
   60    format (/,' Enter the Trajectory Temperature in Degrees',
     &              ' K [298] :  ',$)
         read (input,70,err=80)  tempa
   70    format (f20.0)
         if (tempa .le. 0.0d0)  tempa = 298.0d0
   80    continue
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys1(nkey))
c
c     store the filename and the keyword file for state 0
c
      n0 = n
      fname0 = filename
      leng0 = leng
      title0 = title
      ltitle0 = ltitle
      nkey0 = nkey
      do i = 1, nkey0
         keys0(i) = keyline(i)
      end do
c
c     get trajectory B archive and setup mechanics calculation
c
      call getxyz
      call mechanic
      silent = .true.
c
c     set beginning and ending frame for trajectory B
c
      startb = 0
      stopb = 0
      stepb = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=90,end=90)  startb
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  stopb
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  stepb
   90 continue
      if (query) then
         write (iout,100)
  100    format (/,' Numbers of First & Last File and Step',
     &              ' Increment :  ',$)
         read (input,110)  record
  110    format (a120)
         read (record,*,err=120,end=120)  startb,stopb,stepb
  120    continue
      end if
      if (startb .eq. 0)  startb = 1
      if (stopb .eq. 0)  stopb = maxframe
      if (stepb .eq. 0)  stepb = 1
c
c     find the original temperature value for trajectory B
c
      tempb = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=130,end=130)  tempb
  130 continue
      do while (tempb .lt. 0.0d0)
         write (iout,140)
  140    format (/,' Enter the Trajectory Temperature in Degrees',
     &              ' K [298] :  ',$)
         read (input,150,err=160)  tempb
  150    format (f20.0)
         if (tempb .le. 0.0d0)  tempb = 298.0d0
  160    continue
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (keys1(nkey))
c
c     store the filename and the keyword file for state 1
c
      n1 = n
      fname1 = filename
      leng1 = leng
      title1 = title
      ltitle1 = ltitle
      nkey1 = nkey
      do i = 1, nkey1
         keys1(i) = keyline(i)
      end do
c
c     reopen trajectory A and process the initial structure
c
      ixyz = freeunit ()
      xyzfile = fname0
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      do i = 1, starta-1
         call readxyz (ixyz)
      end do
      call readxyz (ixyz)
      nkey = nkey0
      do i = 1, nkey
         keyline(i) = keys0(i)
      end do
      if (abort) then
         write (iout,170)
  170    format (/,' BAR  --  No Coordinate Frames Available',
     &              ' from First Input File')
         call fatal
      end if
      call mechanic
c
c     find potential energies for trajectory A in state 0
c
      write (iout,180)
  180 format (/,' Initial Processing for Trajectory A :',/)
      j = 0
      k = starta - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         call cutoffs
         ua0(j) = energy ()
         vola(j) = volbox
         do i = 1, stepa-1
            call readxyz (ixyz)
         end do
         call readxyz (ixyz)
         k = k + stepa - 1
         if (k .ge. stopa)  abort = .true.
         if (mod(j,100).eq.0 .or. abort) then
            write (iout,190)  j
  190       format (7x,'Completed',i8,' Coordinate Frames')
            flush (iout)
         end if
      end do
c
c     reset trajectory A and process the initial structure
c
      rewind (unit=ixyz)
      do i = 1, starta-1
         call readxyz (ixyz)
      end do
      call readxyz (ixyz)
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
      call mechanic
c
c     find potential energies for trajectory A in state 1
c
      if (verbose) then
         write (iout,200)
  200    format (/,' Potential Energy Values for Trajectory A :',
     &           //,7x,'Frame',9x,'State 0',9x,'State 1',12x,'Delta',/)
      end if
      j = 0
      k = starta - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         call cutoffs
         ua1(j) = energy ()
         if (verbose) then
            write (iout,210)  k,ua0(j),ua1(j),ua1(j)-ua0(j)
  210       format (i11,2x,3f16.4)
         end if
         do i = 1, stepa-1
            call readxyz (ixyz)
         end do
         call readxyz (ixyz)
         k = k + stepa - 1
         if (k .ge. stopa)  abort = .true.
      end do
      nfrma = j
      frma = dble(nfrma)
      close (unit=ixyz)
c
c     reopen trajectory B and process the initial structure
c
      ixyz = freeunit ()
      xyzfile = fname1
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      do i = 1, startb-1
         call readxyz (ixyz)
      end do
      call readxyz (ixyz)
      nkey = nkey0
      do i = 1, nkey
         keyline(i) = keys0(i)
      end do
      if (abort) then
         write (iout,220)
  220    format (/,' BAR  --  No Coordinate Frames Available',
     &              ' from Second Input File')
         call fatal
      end if
      call mechanic
c
c     find potential energies for trajectory B in state 0
c
      write (iout,230)
  230 format (/,' Initial Processing for Trajectory B :',/)
      j = 0
      k = startb - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         call cutoffs
         ub0(j) = energy ()
         volb(j) = volbox
         do i = 1, stepb-1
            call readxyz (ixyz)
         end do
         call readxyz (ixyz)
         k = k + stepb - 1
         if (k .ge. stopb)  abort = .true.
         if (mod(j,100).eq.0 .or. abort) then
            write (iout,240)  j
  240       format (7x,'Completed',i8,' Coordinate Frames')
            flush (iout)
         end if
      end do
c
c     reset trajectory B and process the initial structure
c
      rewind (unit=ixyz)
      do i = 1, startb-1
         call readxyz (ixyz)
      end do
      call readxyz (ixyz)
      nkey = nkey1
      do i = 1, nkey
         keyline(i) = keys1(i)
      end do
      call mechanic
c
c     find potential energies for trajectory B in state 1
c
      if (verbose) then
         write (iout,250)
  250    format (/,' Potential Energy Values for Trajectory B :',
     &           //,7x,'Frame',9x,'State 0',9x,'State 1',12x,'Delta',/)
      end if
      j = 0
      k = startb - 1
      do while (.not. abort)
         j = j + 1
         k = k + 1
         call cutoffs
         ub1(j) = energy ()
         if (verbose) then
            write (iout,260)  k,ub0(j),ub1(j),ub0(j)-ub1(j)
  260       format (i11,2x,3f16.4)
         end if
         do i = 1, stepb-1
            call readxyz (ixyz)
         end do
         call readxyz (ixyz)
         k = k + stepb - 1
         if (k .ge. stopb)  abort = .true.
      end do
      nfrmb = j
      frmb = dble(nfrmb)
      close (unit=ixyz)
c
c     perform deallocation of some local arrays
c
      deallocate (keys0)
      deallocate (keys1)
c
c     set the frame ratio, temperature and Boltzmann factor
c
      rfrm = frma / frmb
      rta = gasconst * tempa
      rtb = gasconst * tempb
      rt = 0.5d0 * (rta+rtb)
c
c     compute the volume corrections for trajectories A and B
c
      vave = 0.0d0
      do i = 1, nfrma
         vave = vave + vola(i)
      end do
      vave = vave / frma
      if (vave .ne. 0.0d0) then
         do i = 1, nfrma
            if (vola(i) .ne. 0.0d0)  vola(i) = -rta * log(vola(i)/vave)
         end do
      end if
      vave = 0.0d0
      do i = 1, nfrmb
         vave = vave + volb(i)
      end do
      vave = vave / frmb
      if (vave .ne. 0.0d0) then
         do i = 1, nfrmb
            if (volb(i) .ne. 0.0d0)  volb(i) = -rtb * log(volb(i)/vave)
         end do
      end if
c
c     get the free energy change via thermodynamic perturbation
c
      write (iout,270)
  270 format (/,' Estimation of Free Energy Difference',
     &           ' via FEP Method :',/)
      sum = 0.0d0
      do i = 1, nfrma
         sum = sum + exp((ua0(i)-ua1(i)+vola(i))/rta)
      end do
      cfore = -rta * log(sum/frma)
      sum = 0.0d0
      do i = 1, nfrmb
         sum = sum + exp((ub1(i)-ub0(i)+volb(i))/rtb)
      end do
      cback = -rtb * log(sum/frmb)
      write (iout,280)  cfore
  280 format (' FEP Forward Free Energy',4x,f12.4,' Kcal/mol')
      write (iout,290)  cback
  290 format (' FEP Backward Free Energy',3x,f12.4,' Kcal/mol')
c
c     determine the initial free energy via the BAR method
c
      write (iout,300)
  300 format (/,' Estimation of Free Energy Difference',
     &           ' via BAR Method :',/)
      maxiter = 100
      eps = 0.0001d0
      done = .false.
      iter = 0
      cold = 0.0d0
      top = 0.0d0
      top2 = 0.0d0
      do i = 1, nfrmb
         fterm = 1.0d0 / (1.0d0+exp((ub0(i)-ub1(i)+volb(i)+cold)/rtb))
         top = top + fterm
         top2 = top2 + fterm*fterm
      end do
      bot = 0.0d0
      bot2 = 0.0d0
      do i = 1, nfrma
         fterm = 1.0d0 / (1.0d0+exp((ua1(i)-ua0(i)+vola(i)-cold)/rta))
         bot = bot + fterm
         bot2 = bot2 + fterm*fterm
      end do
      cnew = rt*log(rfrm*top/bot) + cold
      stdev = sqrt((bot2-bot*bot/frma)/(bot*bot)
     &                + (top2-top*top/frmb)/(top*top))
      delta = abs(cnew-cold)
      write (iout,310)  iter,cnew
  310 format (' BAR Iteration',i4,10x,f12.4,' Kcal/mol')
      if (delta .lt. eps) then
         done = .true.
         write (iout,320)  cnew,stdev
  320    format (' BAR Free Energy Estimate',3x,f12.4,
     &              ' +/-',f8.4,' Kcal/mol')
      end if
c
c     iterate the BAR equation to converge the free energy
c
      do while (.not. done)
         iter = iter + 1
         cold = cnew
         top = 0.0d0
         top2 = 0.0d0
         do i = 1, nfrmb
            fterm = 1.0d0 / (1.0d0+exp((ub0(i)-ub1(i)+volb(i)
     &                                     +cold)/rtb))
            top = top + fterm
            top2 = top2 + fterm*fterm
         end do
         bot = 0.0d0
         bot2 = 0.0d0
         do i = 1, nfrma
            fterm = 1.0d0 / (1.0d0+exp((ua1(i)-ua0(i)+vola(i)
     &                                     -cold)/rta))
            bot = bot + fterm
            bot2 = bot2 + fterm*fterm
         end do
         cnew = rt*log(rfrm*top/bot) + cold
         stdev = sqrt((bot2-bot*bot/frma)/(bot*bot)
     &                   + (top2-top*top/frmb)/(top*top))
         delta = abs(cnew-cold)
         write (iout,330)  iter,cnew
  330    format (' BAR Iteration',i4,10x,f12.4,' Kcal/mol')
         if (delta .lt. eps) then
            done = .true.
            write (iout,340)  cnew,stdev
  340       format (/,' BAR Free Energy Estimate',3x,f12.4,
     &                 ' +/-',f8.4,' Kcal/mol')
         end if
         if (iter.ge.maxiter .and. .not.done) then
            done = .true.
            write (iout,350)  maxiter
  350       format (/,' BAR Free Energy Estimate not Converged',
     &                 ' after',i4,' Iterations')
            call fatal
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      nbst = 100000
      allocate (bsta(nbst))
      allocate (bstb(nbst))
c
c     use bootstrap analysis to estimate statistical error
c
      write (iout,360)
  360 format (/,' Bootstrap Error Analysis for BAR',
     &           ' Free Energy :')
      if (debug) then
         write (iout,370)
  370    format ()
      end if
      sum = 0.0d0
      sum2 = 0.0d0
      do k = 1, nbst
         done = .false.
         iter = 0
         cold = 0.0d0
         top = 0.0d0
         do i = 1, nfrmb
            bstb(i) = int(frmb*random()) + 1
            j = bstb(i)
            top = top + 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)+volb(i)
     &                                       +cold)/rtb))
         end do
         bot = 0.0d0
         do i = 1, nfrma
            bsta(i) = int(frma*random()) + 1
            j = bsta(i)
            bot = bot + 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)+vola(i)
     &                                       -cold)/rta))
         end do
         cnew = rt*log(rfrm*top/bot) + cold
         delta = abs(cnew-cold)
         do while (.not. done)
            iter = iter + 1
            cold = cnew
            top = 0.0d0
            do i = 1, nfrmb
               j = bstb(i)
               top = top + 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)+volb(i)
     &                                          +cold)/rtb))
            end do
            bot = 0.0d0
            do i = 1, nfrma
               j = bsta(i)
               bot = bot + 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)+vola(i)
     &                                          -cold)/rta))
            end do
            cnew = rt*log(rfrm*top/bot) + cold
            delta = abs(cnew-cold)
            if (delta .lt. eps) then
               done = .true.
               sum = sum + cnew
               sum2 = sum2 + cnew*cnew
               if (debug) then
                  write (iout,380)  k,cnew,iter
  380             format (' Bootstrap Estimate',i7,2x,f12.4,
     &                       ' Kcal/mol at',i4,' Resamples')
               end if
            end if
         end do
      end do
      mean = sum / dble(nbst)
      ratio = dble(nbst/(nbst-1))
      stdev = sqrt(ratio*(sum2/dble(nbst)-mean*mean))
      write (iout,390)  mean,stdev
  390 format (/,' BAR Bootstrap Free Energy',2x,f12.4,
     &           ' +/-',f8.4,' Kcal/mol')
c
c     perform deallocation of some local arrays
c
      deallocate (bsta)
      deallocate (bstb)
c
c     find the enthalpy directly from average potential energies
c
      sum = 0.0d0
      do i = 1, nfrma
         sum = sum + ua0(i)
      end do
      uave0 = sum / frma
      sum = 0.0d0
      do i = 1, nfrmb
         sum = sum + ub1(i)
      end do
      uave1 = sum / frmb
      hdir = uave1 - uave0
      write (iout,400)  hdir
  400 format (/,' Direct Enthalpy Estimate',3x,f12.4,' Kcal/mol')
c
c     get the enthalpy and entropy via thermodynamic perturbation
c
      write (iout,410)
  410 format (/,' Estimation of Enthalpy and Entropy',
     &           ' via FEP Method :',/)
      numer = 0.0d0
      denom = 0.0d0
      do i = 1, nfrma
         term = exp((ua0(i)-ua1(i)+vola(i))/rta)
         numer = numer + ua1(i)*term
         denom = denom + term
      end do
      hfore = (numer/denom) - uave0
      numer = 0.0d0
      denom = 0.0d0
      do i = 1, nfrmb
         term = exp((ub1(i)-ub0(i)+volb(i))/rtb)
         numer = numer + ub0(i)*term
         denom = denom + term
      end do
      hback = -(numer/denom) + uave1
      write (iout,420)  hfore
  420 format (' FEP Forward Enthalpy',7x,f12.4,' Kcal/mol')
      write (iout,430)  hback
  430 format (' FEP Backward Enthalpy',6x,f12.4,' Kcal/mol')
c
c     determine the enthalpy and entropy via the BAR method
c
      write (iout,440)
  440 format (/,' Estimation of Enthalpy and Entropy',
     &           ' via BAR Method :',/)
      fsum = 0.0d0
      fvsum = 0.0d0
      fbvsum = 0.0d0
      vsum = 0.0d0
      fbsum0 = 0.0d0
      do i = 1, nfrma
         fore = 1.0d0 / (1.0d0+exp((ua1(i)-ua0(i)+vola(i)-cnew)/rta))
         back = 1.0d0 / (1.0d0+exp((ua0(i)-ua1(i)-vola(i)+cnew)/rta))
         fsum = fsum + fore
         fvsum = fvsum + fore*ua0(i)
         fbvsum = fbvsum + fore*back*(ua1(i)-ua0(i)+vola(i))
         vsum = vsum + ua0(i)
         fbsum0 = fbsum0 + fore*back
      end do
      alpha0 = fvsum - fsum*(vsum/nfrma) + fbvsum
      bsum = 0.0d0
      bvsum = 0.0d0
      fbvsum = 0.0d0
      vsum = 0.0d0
      fbsum1 = 0.0d0
      do i = 1, nfrmb
         fore = 1.0d0 / (1.0d0+exp((ub1(i)-ub0(i)+volb(i)-cnew)/rtb))
         back = 1.0d0 / (1.0d0+exp((ub0(i)-ub1(i)-volb(i)+cnew)/rtb))
         bsum = bsum + back
         bvsum = bvsum + back*ub1(i)
         fbvsum = fbvsum + fore*back*(ub1(i)-ub0(i)+volb(i))
         vsum = vsum + ub1(i)
         fbsum1 = fbsum1 + fore*back
      end do
      alpha1 = bvsum - bsum*(vsum/nfrmb) - fbvsum
      hbar = (alpha0-alpha1) / (fbsum0+fbsum1)
      tsbar = hbar - cnew
      sbar = tsbar / (0.5d0*(tempa+tempb))
      write (iout,450)  hbar
  450 format (' BAR Enthalpy Estimate',6x,f12.4)
      write (iout,460)  sbar
  460 format (' BAR Entropy Estimate',7x,f12.6)
      write (iout,470)  tsbar
  470 format (' BAR T*dS Estimate',10x,f12.4)
c
c     save the energies and volume corrections to a file
c
      ibar = freeunit ()
      filename = fname0
      barfile = filename(1:leng1)//'.bar'
      call version (barfile,'new')
      open (unit=ibar,file=barfile,status ='new')
      write (ibar,480)  nfrma,title0(1:ltitle0)
  480 format (i6,2x,a)
      do i = 1, nfrma
         write (ibar,490)  i,ua0(i),ua1(i),vola(i)
  490    format (i6,2x,3f16.4)
      end do
      write (ibar,500)  nfrmb,title1(1:ltitle1)
  500 format (i6,2x,a)
      do i = 1, nfrmb
         write (ibar,510)  i,ub0(i),ub1(i),volb(i)
  510    format (i6,2x,3f16.4)
      end do
      close (unit=ibar)
c
c     perform deallocation of some local arrays
c
      deallocate (ua0)
      deallocate (ua1)
      deallocate (ub0)
      deallocate (ub1)
      deallocate (vola)
      deallocate (volb)
c
c     perform any final tasks before program exit
c
      call final
      end

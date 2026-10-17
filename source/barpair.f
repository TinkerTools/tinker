c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine barpair  --  FEP and BAR for a pair of states  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "barpair" computes the free energy, enthalpy and entropy change
c     from state 0 to state 1 via the FEP and BAR methods, given the
c     potential energies of frames sampled in each state evaluated
c     under both states, with bootstrap estimates of the error
c
c     variables and parameters:
c
c     nfrma    number of frames sampled in state 0
c     nfrmb    number of frames sampled in state 1
c     ua0      energies of state 0 frames under state 0
c     ua1      energies of state 0 frames under state 1
c     ub0      energies of state 1 frames under state 0
c     ub1      energies of state 1 frames under state 1
c     vola     volume of state 0 frames for the PV enthalpy term,
c                or zero if not used
c     volb     volume of state 1 frames for the PV enthalpy term,
c                or zero if not used
c     tempa    temperature of the state 0 simulation
c     tempb    temperature of the state 1 simulation
c     patm     pressure in atmospheres for the PV enthalpy term
c     nbst     number of bootstrap trials, at least two
c     verbose  logical flag to print the results as they are found
c     dgfep    forward FEP value and error, backward value and error
c     dgbar    BAR iteration value and error, bootstrap value and error
c     nboot    number of BAR bootstrap trials converged and failed
c     sdu      standard deviation of the forward and backward energy
c                differences over the frames
c     uave     average state 0 energy and error, state 1 energy
c                and error
c     hbar     BAR enthalpy and error, BAR -T*dS value
c     istat    0 if all values are found, 1 if the BAR iteration
c                fails to converge
c
c     each value is computed once from all of the frames, and the BAR
c     enthalpy and entropy use the BAR iteration free energy; the BAR
c     iteration error is the analytical estimate of Bennett, the other
c     errors are standard deviations over bootstrap trials that apply
c     the same estimators to frames resampled with replacement; the
c     BAR bootstrap value is the mean over converged trials only, and
c     the BAR bootstrap and BAR enthalpy errors are found only if at
c     least two trials converge, as given by "nboot", else are zero
c
c
      subroutine barpair (nfrma,nfrmb,ua0,ua1,ub0,ub1,vola,volb,
     &                    tempa,tempb,patm,nbst,verbose,dgfep,dgbar,
     &                    nboot,sdu,uave,hbar,istat)
      use iounit
      use units
      implicit none
      integer i,k,m
      integer nfrma,nfrmb
      integer nfrm,nbst
      integer niter,ntrial
      integer maxiter
      integer nconv,istat
      integer nboot(2)
      integer, allocatable :: ia(:)
      integer, allocatable :: ib(:)
      real*8 rta,rtb
      real*8 frma,frmb
      real*8 tempa,tempb
      real*8 patm,pres
      real*8 bst,ratio
      real*8 random,dave
      real*8 sum2,rconv
      real*8 cfore,cback
      real*8 cinit,cbar
      real*8 stdev,ctrial
      real*8 etrial,htrial
      real*8 vtrla,vtrlb
      real*8 csum,csum2
      real*8 hsum,hsum2
      real*8 uave0,uave1
      real*8 vavea,vaveb
      real*8 hfore,hback
      real*8 sfore,sback
      real*8 epv,hdir
      real*8 hbr,stdhbr
      real*8 sbar,tsbar
      real*8 val(8),sd(8)
      real*8 vsum(8),vsum2(8)
      real*8 ua0(*),ua1(*)
      real*8 ub0(*),ub1(*)
      real*8 vola(*),volb(*)
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
      real*8, allocatable :: hist(:)
      real*8, allocatable :: work(:)
      logical verbose
      logical conv,tconv
      external random
c
c
c     zero out the values returned to the calling routine
c
      istat = 0
      do i = 1, 4
         dgfep(i) = 0.0d0
         dgbar(i) = 0.0d0
         uave(i) = 0.0d0
      end do
      do i = 1, 3
         hbar(i) = 0.0d0
      end do
      do i = 1, 2
         nboot(i) = 0
         sdu(i) = 0.0d0
      end do
c
c     set the frame counts, Boltzmann factors and pressure
c
      frma = dble(nfrma)
      frmb = dble(nfrmb)
      rta = gasconst * tempa
      rtb = gasconst * tempb
      pres = patm / prescon
      nfrm = max(nfrma,nfrmb)
      bst = dble(nbst)
      ratio = bst / (bst-1.0d0)
      maxiter = 100
c
c     perform dynamic allocation of some local arrays
c
      allocate (ia(nfrm))
      allocate (ib(nfrm))
      allocate (hist(0:maxiter))
      allocate (work(0:maxiter))
c
c     find the spread of the energy differences over the frames
c
      if (nfrma .gt. 1) then
         dave = 0.0d0
         do i = 1, nfrma
            dave = dave + ua1(i) - ua0(i)
         end do
         dave = dave / frma
         sum2 = 0.0d0
         do i = 1, nfrma
            sum2 = sum2 + (ua1(i)-ua0(i)-dave)**2
         end do
         sdu(1) = sqrt(sum2/(frma-1.0d0))
      end if
      if (nfrmb .gt. 1) then
         dave = 0.0d0
         do i = 1, nfrmb
            dave = dave + ub0(i) - ub1(i)
         end do
         dave = dave / frmb
         sum2 = 0.0d0
         do i = 1, nfrmb
            sum2 = sum2 + (ub0(i)-ub1(i)-dave)**2
         end do
         sdu(2) = sqrt(sum2/(frmb-1.0d0))
      end if
c
c     find the FEP free energies, enthalpies and averages from all
c     of the frames, where the enthalpies include the PV term
c
      do i = 1, nfrma
         ia(i) = i
      end do
      do i = 1, nfrmb
         ib(i) = i
      end do
      call barfep (nfrma,ia,ua0,ua1,vola,rta,pres,cfore,hfore,
     &             uave0,vavea)
      call barfep (nfrmb,ib,ub1,ub0,volb,rtb,pres,cback,hback,
     &             uave1,vaveb)
      cback = -cback
      hback = -hback
      epv = pres * (vaveb-vavea)
      hdir = uave1 - uave0 + epv
c
c     find the BAR free energy from all of the frames, starting from
c     the mean of the FEP values, then the BAR enthalpy at that value
c
      cinit = 0.5d0 * (cfore+cback)
      if (.not. (abs(cinit) .lt. 1.0d100))  cinit = 0.0d0
      call bariter (nfrma,ia,nfrmb,ib,ua0,ua1,ub0,ub1,rta,rtb,cinit,
     &              maxiter,cbar,stdev,niter,hist,conv)
      hbr = 0.0d0
      if (conv) then
         call barhbar (nfrma,ia,nfrmb,ib,ua0,ua1,ub0,ub1,vola,volb,
     &                 rta,rtb,pres,cbar,hbr)
      end if
c
c     apply the same estimators to bootstrap resamples of the frames,
c     with each BAR trial started from the converged BAR value and
c     dropped if it fails to converge
c
      do m = 1, 8
         vsum(m) = 0.0d0
         vsum2(m) = 0.0d0
      end do
      nconv = 0
      csum = 0.0d0
      csum2 = 0.0d0
      hsum = 0.0d0
      hsum2 = 0.0d0
      do k = 1, nbst
         do i = 1, nfrma
            ia(i) = int(frma*random()) + 1
         end do
         do i = 1, nfrmb
            ib(i) = int(frmb*random()) + 1
         end do
         call barfep (nfrma,ia,ua0,ua1,vola,rta,pres,val(1),val(7),
     &                val(3),vtrla)
         call barfep (nfrmb,ib,ub1,ub0,volb,rtb,pres,val(2),val(8),
     &                val(4),vtrlb)
         val(2) = -val(2)
         val(8) = -val(8)
         val(5) = pres * (vtrlb-vtrla)
         val(6) = val(4) - val(3) + val(5)
         do m = 1, 8
            vsum(m) = vsum(m) + val(m)
            vsum2(m) = vsum2(m) + val(m)*val(m)
         end do
         if (conv) then
            call bariter (nfrma,ia,nfrmb,ib,ua0,ua1,ub0,ub1,rta,rtb,
     &                    cbar,maxiter,ctrial,etrial,ntrial,work,tconv)
            if (tconv) then
               nconv = nconv + 1
               csum = csum + ctrial
               csum2 = csum2 + ctrial*ctrial
               call barhbar (nfrma,ia,nfrmb,ib,ua0,ua1,ub0,ub1,vola,
     &                       volb,rta,rtb,pres,ctrial,htrial)
               hsum = hsum + htrial
               hsum2 = hsum2 + htrial*htrial
            end if
         end if
      end do
      do m = 1, 8
         sd(m) = vsum(m) / bst
         sd(m) = sqrt(max(0.0d0,ratio*(vsum2(m)/bst-sd(m)*sd(m))))
      end do
c
c     report the free energy change via thermodynamic perturbation
c
      dgfep(1) = cfore
      dgfep(2) = sd(1)
      dgfep(3) = cback
      dgfep(4) = sd(2)
      if (verbose) then
         write (iout,10)
   10    format (/,' Free Energy Difference via FEP Method :',/)
         write (iout,20)  cfore,sd(1)
   20    format (' Free Energy via Forward FEP',9x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         write (iout,30)  cback,sd(2)
   30    format (' Free Energy via Backward FEP',8x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
      end if
c
c     report the free energy change via the BAR iteration
c
      if (verbose) then
         write (iout,40)
   40    format (/,' Free Energy Difference via BAR Method :',/)
         do i = 0, niter
            write (iout,50)  i,hist(i)
   50       format (' BAR Iteration',i4,19x,f12.4,' Kcal/mol')
         end do
         if (conv .and. niter.eq.0) then
            write (iout,60)  cbar,stdev
   60       format (' BAR Free Energy Estimate',8x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         else if (conv) then
            write (iout,70)  cbar,stdev
   70       format (/,' Free Energy via BAR Estimate',8x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         end if
      end if
c
c     report failure of the BAR iteration to the calling routine
c
      if (.not. conv) then
         istat = 1
         if (verbose) then
            write (iout,80)  niter
   80       format (/,' BAR Free Energy Estimate not Converged',
     &                 ' after',i4,' Iterations')
         end if
         deallocate (ia)
         deallocate (ib)
         deallocate (hist)
         deallocate (work)
         return
      end if
      dgbar(1) = cbar
      dgbar(2) = stdev
c
c     report the mean and spread of the BAR bootstrap trials
c
      nboot(1) = nconv
      nboot(2) = nbst - nconv
      rconv = dble(nconv)
      stdhbr = 0.0d0
      if (nconv .ge. 1) then
         dgbar(3) = csum / rconv
      end if
      if (nconv .ge. 2) then
         dgbar(4) = sqrt(max(0.0d0,(rconv/(rconv-1.0d0))
     &                       *(csum2/rconv-dgbar(3)*dgbar(3))))
         stdhbr = hsum / rconv
         stdhbr = sqrt(max(0.0d0,(rconv/(rconv-1.0d0))
     &                       *(hsum2/rconv-stdhbr*stdhbr)))
         if (verbose) then
            write (iout,90)  dgbar(3),dgbar(4)
   90       format (' Free Energy via BAR Bootstrap',7x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         end if
      else if (verbose) then
         write (iout,100)
  100    format (' Free Energy via BAR Bootstrap',13x,
     &              'not Available')
      end if
      if (verbose .and. nconv.lt.nbst) then
         write (iout,110)  nbst-nconv,nbst
  110    format (' BAR Bootstrap :',i8,' of',i8,
     &              ' Trials did not Converge')
      end if
c
c     report the enthalpy directly via average potential energy
c
      uave(1) = uave0
      uave(2) = sd(3)
      uave(3) = uave1
      uave(4) = sd(4)
      if (verbose) then
         write (iout,120)
  120    format (/,' Enthalpy from Potential Energy Averages :',/)
         write (iout,130)  uave0,sd(3)
  130    format (' Average Energy for State 0',10x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         write (iout,140)  uave1,sd(4)
  140    format (' Average Energy for State 1',10x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         if (epv .ne. 0.0d0) then
            write (iout,150)  patm,epv,sd(5)
  150       format (' PdV Work Term at',f9.2,' Atm',7x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         end if
         write (iout,160)  hdir,sd(6)
  160    format (' Enthalpy via Direct Estimate',8x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
      end if
c
c     report the enthalpy and entropy via thermodynamic perturbation
c
      sfore = (hfore-cfore) / tempa
      sback = (hback-cback) / tempb
      if (verbose) then
         write (iout,170)
  170    format (/,' Enthalpy and Entropy via FEP Method :',/)
         write (iout,180)  hfore,sd(7)
  180    format (' Enthalpy via Forward FEP',12x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         write (iout,190)  sfore
  190    format (' Entropy via Forward FEP',13x,f12.6,' Kcal/mol/K')
         write (iout,200)  -tempa*sfore
  200    format (' Forward FEP -T*dS Value',13x,f12.4,' Kcal/mol')
         write (iout,210)  hback,sd(8)
  210    format (/,' Enthalpy via Backward FEP',11x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         write (iout,220)  sback
  220    format (' Entropy via Backward FEP',12x,f12.6,' Kcal/mol/K')
         write (iout,230)  -tempb*sback
  230    format (' Backward FEP -T*dS Value',12x,f12.4,' Kcal/mol')
      end if
c
c     report the enthalpy and entropy via the BAR method, where the
c     enthalpy error needs at least two converged bootstrap trials
c
      tsbar = hbr - cbar
      sbar = tsbar / (0.5d0*(tempa+tempb))
      hbar(1) = hbr
      hbar(2) = stdhbr
      hbar(3) = -tsbar
      if (verbose) then
         write (iout,240)
  240    format (/,' Enthalpy and Entropy via BAR Method :',/)
         if (nconv .ge. 2) then
            write (iout,250)  hbr,stdhbr
  250       format (' Enthalpy via BAR Estimate',11x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         else
            write (iout,260)  hbr
  260       format (' Enthalpy via BAR Estimate',11x,f12.4,
     &                 ' Kcal/mol, Error not Available')
         end if
         write (iout,270)  sbar
  270    format (' Entropy via BAR Estimate',12x,f12.6,' Kcal/mol/K')
         write (iout,280)  -tsbar
  280    format (' BAR Estimate -T*dS Value',12x,f12.4,' Kcal/mol')
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (ia)
      deallocate (ib)
      deallocate (hist)
      deallocate (work)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine barfep  --  one-direction FEP over a sample  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "barfep" finds the free energy and enthalpy change from the
c     sampled state to the other state via thermodynamic perturbation
c     over a list of frames, which may repeat, along with the average
c     energy and volume of those frames; the exponentials are shifted
c     by their largest exponent so the sum cannot overflow, and the
c     enthalpy includes the PV term
c
c     variables and parameters:
c
c     n      number of frames in the list
c     idx    frame indices to use
c     uo     energies of the frames under the sampled state
c     ut     energies of the frames under the other state
c     vol    volumes of the frames, or zero if not used
c     rt     gas constant times the sampling temperature
c     pres   pressure in Kcal/mole/Ang**3
c     dg     free energy change from sampled to other state
c     dh     enthalpy change from sampled to other state
c     uave   average energy of the frames under the sampled state
c     vave   average volume of the frames
c
c
      subroutine barfep (n,idx,uo,ut,vol,rt,pres,dg,dh,uave,vave)
      implicit none
      integer i,j,n
      integer idx(*)
      real*8 rt,pres
      real*8 dg,dh
      real*8 uave,vave
      real*8 x,xmax,w
      real*8 wsum,hwsum
      real*8 usum,vsum
      real*8 uo(*),ut(*)
      real*8 vol(*)
c
c
c     find the largest exponent over the frames
c
      xmax = 0.0d0
      do i = 1, n
         j = idx(i)
         x = (uo(j)-ut(j)) / rt
         if (i.eq.1 .or. x.gt.xmax)  xmax = x
      end do
c
c     accumulate the shifted weights and weighted enthalpies
c
      wsum = 0.0d0
      hwsum = 0.0d0
      usum = 0.0d0
      vsum = 0.0d0
      do i = 1, n
         j = idx(i)
         w = exp((uo(j)-ut(j))/rt-xmax)
         wsum = wsum + w
         hwsum = hwsum + (ut(j)+pres*vol(j))*w
         usum = usum + uo(j)
         vsum = vsum + vol(j)
      end do
      uave = usum / dble(n)
      vave = vsum / dble(n)
      dg = -rt * (xmax+log(wsum/dble(n)))
      dh = hwsum/wsum - (uave+pres*vave)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine bariter  --  BAR free energy by iteration  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "bariter" solves the Bennett acceptance ratio equation for the
c     free energy change from state 0 to state 1 by iteration over
c     lists of frames sampled in each state, starting from an initial
c     value; also finds the analytical error of Bennett
c
c     variables and parameters:
c
c     na       number of state 0 frames in the list
c     ia       state 0 frame indices to use
c     nb       number of state 1 frames in the list
c     ib       state 1 frame indices to use
c     ua0      energies of state 0 frames under state 0
c     ua1      energies of state 0 frames under state 1
c     ub0      energies of state 1 frames under state 0
c     ub1      energies of state 1 frames under state 1
c     rta      gas constant times the state 0 temperature
c     rtb      gas constant times the state 1 temperature
c     cinit    initial value of the free energy change
c     maxiter  maximum number of iterations
c     cnew     final value of the free energy change
c     stdev    analytical error of the final value
c     niter    number of the last iteration, counting from zero
c     hist     value of the free energy change at each iteration
c     conv     logical flag set if the iteration converged
c
c
      subroutine bariter (na,ia,nb,ib,ua0,ua1,ub0,ub1,rta,rtb,cinit,
     &                    maxiter,cnew,stdev,niter,hist,conv)
      implicit none
      integer i,j
      integer na,nb
      integer niter,maxiter
      integer ia(*),ib(*)
      real*8 rta,rtb,rt
      real*8 cinit,cold,cnew
      real*8 stdev,eps
      real*8 frma,frmb,rfrm
      real*8 top,top2
      real*8 bot,bot2
      real*8 fterm
      real*8 ua0(*),ua1(*)
      real*8 ub0(*),ub1(*)
      real*8 hist(0:*)
      logical conv,fail
c
c
c     set the convergence criterion and frame ratio
c
      eps = 0.0001d0
      rt = 0.5d0 * (rta+rtb)
      frma = dble(na)
      frmb = dble(nb)
      rfrm = frma / frmb
c
c     iterate the BAR equation until the value stops changing
c
      cnew = cinit
      stdev = 0.0d0
      niter = -1
      conv = .false.
      fail = .false.
      do while (.not.conv .and. .not.fail)
         niter = niter + 1
         cold = cnew
         top = 0.0d0
         top2 = 0.0d0
         do i = 1, nb
            j = ib(i)
            fterm = 1.0d0 / (1.0d0+exp((ub0(j)-ub1(j)+cold)/rtb))
            top = top + fterm
            top2 = top2 + fterm*fterm
         end do
         bot = 0.0d0
         bot2 = 0.0d0
         do i = 1, na
            j = ia(i)
            fterm = 1.0d0 / (1.0d0+exp((ua1(j)-ua0(j)-cold)/rta))
            bot = bot + fterm
            bot2 = bot2 + fterm*fterm
         end do
         cnew = rt*log(rfrm*top/bot) + cold
         stdev = rt * sqrt(max(0.0d0,(bot2-bot*bot/frma)/(bot*bot)
     &                        +(top2-top*top/frmb)/(top*top)))
         hist(niter) = cnew
         if (abs(cnew-cold) .lt. eps) then
            conv = .true.
         else if (niter .ge. maxiter) then
            fail = .true.
         else if (.not. (abs(cnew) .lt. 1.0d100)) then
            fail = .true.
         end if
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine barhbar  --  BAR enthalpy at a free energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "barhbar" finds the enthalpy change from state 0 to state 1 via
c     the BAR estimator of Wyczalkowski et al. over lists of frames
c     sampled in each state, given the BAR free energy change for the
c     same frames; the state energies include the PV term, which
c     cancels in the energy differences
c
c     variables and parameters:
c
c     na       number of state 0 frames in the list
c     ia       state 0 frame indices to use
c     nb       number of state 1 frames in the list
c     ib       state 1 frame indices to use
c     ua0      energies of state 0 frames under state 0
c     ua1      energies of state 0 frames under state 1
c     ub0      energies of state 1 frames under state 0
c     ub1      energies of state 1 frames under state 1
c     vola     volume of state 0 frames, or zero if not used
c     volb     volume of state 1 frames, or zero if not used
c     rta      gas constant times the state 0 temperature
c     rtb      gas constant times the state 1 temperature
c     pres     pressure in Kcal/mole/Ang**3
c     cbar     BAR free energy change for these frames
c     hbr      BAR enthalpy change
c
c
      subroutine barhbar (na,ia,nb,ib,ua0,ua1,ub0,ub1,vola,volb,
     &                    rta,rtb,pres,cbar,hbr)
      implicit none
      integer i,j
      integer na,nb
      integer ia(*),ib(*)
      real*8 rta,rtb
      real*8 pres,cbar,hbr
      real*8 fore,back,h
      real*8 fsum,bsum
      real*8 fvsum,bvsum
      real*8 fbvsum,vsum
      real*8 fbsum0,fbsum1
      real*8 alpha0,alpha1
      real*8 ua0(*),ua1(*)
      real*8 ub0(*),ub1(*)
      real*8 vola(*),volb(*)
c
c
c     find the terms from the state 0 frames
c
      fsum = 0.0d0
      fvsum = 0.0d0
      fbvsum = 0.0d0
      vsum = 0.0d0
      fbsum0 = 0.0d0
      do i = 1, na
         j = ia(i)
         fore = 1.0d0/(1.0d0+exp((ua1(j)-ua0(j)-cbar)/rta))
         back = 1.0d0/(1.0d0+exp((ua0(j)-ua1(j)+cbar)/rta))
         h = ua0(j) + pres*vola(j)
         fsum = fsum + fore
         fvsum = fvsum + fore*h
         fbvsum = fbvsum + fore*back*(ua1(j)-ua0(j))
         vsum = vsum + h
         fbsum0 = fbsum0 + fore*back
      end do
      alpha0 = fvsum - fsum*(vsum/dble(na)) + fbvsum
c
c     find the terms from the state 1 frames
c
      bsum = 0.0d0
      bvsum = 0.0d0
      fbvsum = 0.0d0
      vsum = 0.0d0
      fbsum1 = 0.0d0
      do i = 1, nb
         j = ib(i)
         fore = 1.0d0/(1.0d0+exp((ub1(j)-ub0(j)-cbar)/rtb))
         back = 1.0d0/(1.0d0+exp((ub0(j)-ub1(j)+cbar)/rtb))
         h = ub1(j) + pres*volb(j)
         bsum = bsum + back
         bvsum = bvsum + back*h
         fbvsum = fbvsum + fore*back*(ub1(j)-ub0(j))
         vsum = vsum + h
         fbsum1 = fbsum1 + fore*back
      end do
      alpha1 = bvsum - bsum*(vsum/dble(nb)) - fbvsum
      hbr = (alpha0-alpha1) / (fbsum0+fbsum1)
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine lmqn  --  conjugate gradient optimization  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "lmqn" is a general conjugate gradient optimization routine
c     which implements steepest descent, Fletcher-Reeves CG,
c     Polak-Ribiere CG, Hestenes-Stiefel CG, Powell-Beale CG with
c     restarts, and a memoryless BFGS quasi-Newton method
c
c     literature references:
c
c     D. G. Luenberger, "Linear and Nonlinear Programming", 2nd
c     Edition, Addison-Wesley, Reading, MA, 1984, Section 9-7
c
c     S. G. Nash and A. Sofer, "Linear and Nonlinear Programming",
c     McGraw-Hill, New York, 1996, Chapter 12
c
c     parameters used in the main iteration:
c
c     nvar     number of parameters in the objective function
c     method   steepest descent, FRCG, PRCG, HSCG, Powell or LMQN
c     fctmin   normal exit if function gets less than
c     grdmin   normal exit if gradient norm gets less than
c     maxiter  error return if number of iterations exceeds
c     period   restart at least every period iterations
c     iprint   print iteration results every iprint iterations
c     iwrite   call user-supplied output every iwrite iterations
c     fast     steepest descent while function decrease exceeds
c     slow     restart if relative function decrease drops below
c     epsln    error if total move is less than
c     scale    factor by which actual function has been multiplied
c     rms      factor to convert grad norm and movement to rms
c     minimum  final value of the objective function
c     ncalls   number of function/gradient (fgvalue) calls
c     niter    number of conjugate gradient iterations performed
c     status   string containing informative termination message
c
c     parameters used in the line search:
c
c     cappa    reduction in projected gradient for termination
c     stpmin   minimum allowed line search step size
c     stpmax   maximum allowed line search step size
c     angmax   maximum angle between search and -grad directions
c     intmax   maximum number of interpolations in line search
c
c     vectors stored by the routine:
c
c     x        current parameter values
c     x_old    previous parameter values
c     g        current gradient vector
c     g_old    previous gradient vector
c     s        current minus previous parameter values
c     d        current minus previous gradient values
c     p        conjugate direction search vector
c
c     requried external routines:
c
c     fgvalue    function to evaluate function and gradient values
c     optsave    subroutine to write out info about current status
c
c
      subroutine lmqn (nvar,x,minimum,grdmin,fgvalue,optsave)
      implicit none
      include 'sizes.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'linmin.i'
      include 'math.i'
      include 'minima.i'
      include 'output.i'
      include 'potent.i'
      include 'scales.i'
      integer i,nvar,next
      integer ncalls,nerror
      integer ipow,period
      integer niter,nstart
      real*8 grdmin,fast
      real*8 slow,epsln
      real*8 f,f_old
      real*8 f_new,f_move
      real*8 rms,beta,gamma
      real*8 angle,ratio
      real*8 fgvalue,minimum
      real*8 x_move,g_rms
      real*8 gg,gg_old
      real*8 sg,dg,dp,sd,dd
      real*8 dtg,dtpt
      real*8 x(maxvar)
      real*8 g(maxvar)
      real*8 x_old(maxvar)
      real*8 g_old(maxvar)
      real*8 p(maxvar)
      real*8 s(maxvar)
      real*8 d(maxvar)
      real*8 pt(maxvar)
      real*8 dt(maxvar)
      logical restart,done
      character*3 method
      character*9 blank,status
      character*20 keyword
      character*80 record
      external fgvalue,optsave
c
c
c     initialize some values to be used below
c
      if (nvar .gt. maxvar) then
         write (iout,10)
   10    format (/,' LMQN  --  Too many Parameters,',
     &              ' Increase the Value of MAXVAR')
         return
      end if
      ncalls = 0
      rms = sqrt(dble(nvar))
      if (coordtype .eq. 'CARTESIAN') then
         rms = rms / sqrt(3.0d0)
      else if (coordtype .eq. 'RIGIDBODY') then
         rms = rms / sqrt(6.0d0)
      end if
      blank = '         '
      status = blank
      method = 'MQN'
      restart = .true.
      done = .false.
c
c     set default values for variable scale factors
c
      if (.not. set_scale) then
         do i = 1, nvar
            if (scale(i) .eq. 0.0d0)  scale(i) = 1.0d0
         end do
      end if
c
c     set default parameters for the optimization
c
      if (fctmin .eq. 0.0d0)  fctmin = -1000000.0d0
      if (maxiter .eq. 0)  maxiter = 1000000
      if (nextiter .eq. 0)  nextiter = 1
      if (iprint .lt. 0)  iprint = 1
      if (iwrite .lt. 0)  iwrite = 1
      fast = 0.5d0
      slow = 0.0d0
      epsln = 1.0d-16
      period = max(200,nvar)
c
c     set default parameters for the line search
c
      if (stpmin .eq. 0.0d0)  stpmin = 1.0d-16
      if (stpmax .eq. 0.0d0)  stpmax = 5.0d0
      if (cappa .eq. 0.0d0)  cappa = 0.1d0
      slpmax = 10000.0d0
      angmax = 88.0d0
      intmax = 5
c
c     search the keywords for optimization parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:17) .eq. 'STEEPEST-DESCENT ') then
            method = 'SD'
         else if (keyword(1:16) .eq. 'FLETCHER-REEVES ') then
            method = 'FR'
         else if (keyword(1:14) .eq. 'POLAK-RIBIERE ') then
            method = 'PR'
         else if (keyword(1:17) .eq. 'HESTENES-STIEFEL ') then
            method = 'HS'
         else if (keyword(1:13) .eq. 'POWELL-BEALE ') then
            method = 'POW'
         else if (keyword(1:7) .eq. 'FCTMIN ') then
            read (record(next:80),*,err=20,end=20)  fctmin
         else if (keyword(1:8) .eq. 'MAXITER ') then
            read (record(next:80),*,err=20,end=20)  maxiter
         else if (keyword(1:9) .eq. 'NEXTITER ') then
            read (record(next:80),*,err=20,end=20)  nextiter
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (record(next:80),*,err=20,end=20)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (record(next:80),*,err=20,end=20)  iwrite
         else if (keyword(1:8) .eq. 'STEPMAX ') then
            read (record(next:80),*,err=20,end=20)  stpmax
         else if (keyword(1:8) .eq. 'STEPMIN ') then
            read (record(next:80),*,err=20,end=20)  stpmin
         else if (keyword(1:6) .eq. 'CAPPA ') then
            read (record(next:80),*,err=20,end=20)  cappa
         else if (keyword(1:9) .eq. 'SLOPEMAX ') then
            read (record(next:80),*,err=20,end=20)  slpmax
         else if (keyword(1:7) .eq. 'ANGMAX ') then
            read (record(next:80),*,err=20,end=20)  angmax
         else if (keyword(1:6) .eq. 'EPSLN ') then
            read (record(next:80),*,err=20,end=20)  epsln
         else if (keyword(1:5) .eq. 'FAST ') then
            read (record(next:80),*,err=20,end=20)  fast
         else if (keyword(1:5) .eq. 'SLOW ') then
            read (record(next:80),*,err=20,end=20)  slow
         else if (keyword(1:7) .eq. 'PERIOD ') then
            read (record(next:80),*,err=20,end=20)  period
         else if (keyword(1:7) .eq. 'INTMAX ') then
            read (record(next:80),*,err=20,end=20)  intmax
         end if
   20    continue
      end do
c
c     print header information prior to iterations
c
      if (iprint .gt. 0) then
         if (method .eq. 'SD') then
            write (iout,30)
   30       format (/,' Steepest Descent Gradient Optimization :')
         else if (method .eq. 'FR') then
            write (iout,40)
   40       format (/,' Fletcher-Reeves Conjugate Gradient',
     &                 ' Optimization :')
         else if (method .eq. 'PR') then
            write (iout,50)
   50       format (/,' Polak-Ribiere Conjugate Gradient',
     &                 ' Optimization :')
         else if (method .eq. 'HS') then
            write (iout,60)
   60       format (/,' Hestenes-Stiefel Conjugate Gradient',
     &                 ' Optimization :')
         else if (method .eq. 'POW') then
            write (iout,70)
   70       format (/,' Powell-Beale Conjugate Gradient Optimization :')
         else if (method .eq. 'MQN') then
            write (iout,80)
   80       format (/,' Memoryless BFGS Quasi-Newton Optimization :')
         end if
         if (method .eq. 'SD') then
            write (iout,90)
   90       format (/,' SD Iter    F Value      G RMS     F Move',
     &                 '    X Move    Angle  FG Call  Comment',/)
         else
            write (iout,100)
  100       format (/,' CG Iter    F Value      G RMS     F Move',
     &                 '    X Move    Angle  FG Call  Comment',/)
         end if
      end if
c
c     get initial function and gradient values
c
      niter = nextiter - 1
      maxiter = niter + maxiter
      ncalls = ncalls + 1
      f = fgvalue (x,g)
      g_rms = 0.0d0
      f_move = 0.0d0
      do i = 1, nvar
         x_old(i) = x(i)
         g_old(i) = g(i)
         g_rms = g_rms + (g(i)*scale(i))**2
         f_move = f_move + g(i)**2
      end do
      g_rms = sqrt(g_rms) / rms
      f_move = 0.5d0 * stpmax * sqrt(f_move)
c
c     print initial information prior to first iteration
c
      if (iprint .gt. 0) then
         if (f.lt.1.0d7 .and. f.gt.-1.0d6 .and. g_rms.lt.1.0d5) then
            write (iout,110)  niter,f,g_rms,ncalls
  110       format (i6,f13.4,f11.4,30x,i7)
         else
            write (iout,120)  niter,f,g_rms,ncalls
  120       format (i6,d13.4,d11.4,30x,i7)
         end if
      end if
c
c     write initial intermediate prior to first iteration
c
      if (iwrite .gt. 0)  call optsave (niter,f,x)
c
c     tests of the various termination criteria
c
      if (niter .ge. maxiter) then
         status = 'IterLimit'
         done = .true.
      end if
      if (f .le. fctmin) then
         status = 'SmallFct '
         done = .true.
      end if
      if (g_rms .le. grdmin) then
         status = 'SmallGrad'
         done = .true.
      end if
c
c     start of a new conjugate gradient iteration
c
      dowhile (.not. done)
         niter = niter + 1
         if (status .eq. blank)  nerror = 0
c
c     compute the next search direction using Steepest Descent,
c     Fletcher-Reeves CG, Polak-Ribiere CG or Memoryless BFGS
c
  130    continue
         status = blank
         if (method.eq.'SD' .or. restart) then
            do i = 1, nvar
               p(i) = -g(i)
            end do
            nstart = niter
            restart = .false.
         else if (method .eq. 'FR') then
            gg = 0.0d0
            gg_old = 0.0d0
            do i = 1, nvar
               gg = gg + g(i)*g(i)
               gg_old = gg_old + g_old(i)*g_old(i)
            end do
            beta = gg / gg_old
            do i = 1, nvar
               p(i) = -g(i) + beta*p(i)
            end do
         else if (method .eq. 'PR') then
            dg = 0.0d0
            gg_old = 0.0d0
            do i = 1, nvar
               dg = dg + d(i)*g(i)
               gg_old = gg_old + g_old(i)*g_old(i)
            end do
            beta = dg / gg_old
            do i = 1, nvar
               p(i) = -g(i) + beta*p(i)
            end do
         else if (method .eq. 'HS') then
            dg = 0.0d0
            dp = 0.0d0
            do i = 1, nvar
               dg = dg + d(i)*g(i)
               dp = dp + d(i)*p(i)
            end do
            beta = dg / dp
            do i = 1, nvar
               p(i) = -g(i) + beta*p(i)
            end do
         else if (method .eq. 'POW') then
            dg = 0.0d0
            dp = 0.0d0
            do i = 1, nvar
               dg = dg + d(i)*g(i)
               dp = dp + d(i)*p(i)
            end do
            beta = dg / dp
            if (niter .eq. ipow) then
               gamma = 0.0d0
            else
               dtg = 0.0d0
               dtpt = 0.0d0
               do i = 1, nvar
                  dtg = dtg + dt(i)*g(i)
                  dtpt = dtpt + dt(i)*pt(i)
               end do
               gamma = dtg / dtpt
            end if
            do i = 1, nvar
               p(i) = -g(i) + beta*p(i) + gamma*pt(i)
            end do
         else if (method .eq. 'MQN') then
            sg = 0.0d0
            dg = 0.0d0
            dd = 0.0d0
            sd = 0.0d0
            do i = 1, nvar
               sg = sg + s(i)*g(i)
               dg = dg + d(i)*g(i)
               dd = dd + d(i)*d(i)
               sd = sd + s(i)*d(i)
            end do
            do i = 1, nvar
               p(i) = -g(i) + (d(i)*sg+s(i)*dg)/sd
     &                   - (1.0d0+dd/sd)*(s(i)*sg/sd)
            end do
         end if
c
c     perform line search along the new conjugate direction
c
         f_old = f
         call search (nvar,f,g,x,p,f_move,angle,ncalls,fgvalue,status)
         f_new = f
c
c     if angle between the search direction and the negative
c     gradient vector was too large, use steepest descent
c
         if (status .eq. 'WideAngle') then
            restart = .true.
            goto 130
         end if
c
c     special test for reset of the Powell method
c
         if (method .eq. 'POW') then
            restart = .false.
            gg = 0.0d0
            gg_old = 0.0d0
            do i = 1, nvar
               gg = gg + g(i)*g(i)
               gg_old = gg_old + g(i)*g_old(i)
            end do
            ratio = gg_old / gg
            if (niter.eq.1 .or. ratio.ge.0.2d0) then
               ipow = niter + 1
               do i = 1, nvar
                  pt(i) = p(i)
                  dt(i) = g(i) - g_old(i)
               end do
               nstart = niter
               status = 'Resetting'
            end if
         end if
c
c     update variables based on results of this iteration
c
         f_move = f_old - f_new
         x_move = 0.0d0
         g_rms = 0.0d0
         do i = 1, nvar
            s(i) = x(i) - x_old(i)
            d(i) = g(i) - g_old(i)
            x_move = x_move + (s(i)/scale(i))**2
            g_rms = g_rms + (g(i)*scale(i))**2
            x_old(i) = x(i)
            g_old(i) = g(i)
         end do
         x_move = sqrt(x_move) / rms
         if (coordtype .eq. 'INTERNAL') then
            x_move = radian * x_move
         end if
         g_rms = sqrt(g_rms) / rms
c
c     test for error/restart due to line search problems
c
         if (status .eq. 'BadIntpln') then
            nerror = nerror + 1
            if (nerror .eq. 3) then
               done = .true.
            else
               restart = .true.
            end if
         end if
         if (status .eq. 'IntplnErr') then
            nerror = nerror + 1
            if (nerror .eq. 3) then
               done = .true.
            else
               restart = .true.
            end if
            do i = 1, nvar
               x(i) = x_old(i)
               g(i) = g_old(i)
            end do
         end if
c
c     test for error/restart due to lack of movement
c
         if (x_move .lt. epsln) then
            status = 'SmallMove'
            nerror = nerror + 1
            if (nerror .eq. 3) then
               done = .true.
            else
               restart = .true.
            end if
         end if
c
c     test for error/restart due to function increase
c
         if (f_move .lt. 0.0d0) then
            status = 'Increase '
            nerror = nerror + 1
            if (nerror .eq. 3) then
               done = .true.
            else
               restart = .true.
            end if
            do i = 1, nvar
               x(i) = x_old(i)
               g(i) = g_old(i)
            end do
         end if
c
c     test for error/restart due to slow progress
c
         ratio = f_move / abs(f_new-fctmin)
         if (abs(ratio) .lt. slow) then
            status = 'SlowDecr '
            nerror = nerror + 1
            if (nerror .eq. 3) then
               done = .true.
            else
               restart = .true.
            end if
         end if
c
c     test for restart due to fast progress
c
         if (ratio .gt. fast) then
            status = 'FastDecr '
            nerror = 0
            restart = .true.
         end if
c
c     test for restart based on iterations since last restart
c
         if (niter-nstart .ge. period) then
            status = 'Periodic '
            restart = .true.
         end if
c
c     test for too many total iterations
c
         if (niter .ge. maxiter) then
            status = 'IterLimit'
            done = .true.
         end if
c
c     test the normal termination criteria
c
         if (f .le. fctmin) then
            status = 'SmallFct '
            done = .true.
         end if
         if (g_rms .le. grdmin) then
            status = 'SmallGrad'
            done = .true.
         end if
c
c     print intermediate results for the current iteration
c
         if (iprint .gt. 0) then
            if (done .or. mod(niter,iprint).eq.0) then
               if (f.lt.1.0d7 .and. f.gt.-1.0d6 .and.
     &             g_rms.lt.1.0d5 .and. f_move.lt.1.0d5) then
                  write (iout,140)  niter,f,g_rms,f_move,
     &                              x_move,angle,ncalls,status
  140             format (i6,f13.4,f11.4,f11.4,f10.4,f9.2,i7,3x,a9)
               else
                  write (iout,150)  niter,f,g_rms,f_move,
     &                              x_move,angle,ncalls,status
  150             format (i6,d13.4,d11.4,d11.4,f10.4,f9.2,i7,3x,a9)
               end if
            end if
         end if
c
c     write intermediate results for the current iteration
c
         if (iwrite .gt. 0) then
            if (done .or. mod(niter,iwrite).eq.0) then
               call optsave (niter,f,x)
            end if
         end if
      end do
c
c     set final value of the objective function
c
      minimum = f
      if (iprint .gt. 0) then
         if (status.eq.'SmallGrad' .or. status.eq.'SmallFct ') then
            write (iout,160)  status
  160       format (/,' LMQN  --  Normal Termination due to ',a9)
         else
            write (iout,170)  status
  170       format (/,' LMQN  --  Incomplete Convergence due to ',a9)
         end if
      end if
      return
      end

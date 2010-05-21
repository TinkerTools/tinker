c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine square  --  nonlinear least squares with bounds  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "square" is a nonlinear least squares routine derived from
c     the IMSL routine BCLSF and More's Minpack routine LMDER; the
c     Jacobian is estimated by finite differences and bounds can
c     be specified for the variables to be refined
c
c     arguments and variables :
c
c     n        number of variables to optimize
c     m        number of residual functions
c     xlo      vector of length n containing the lower bounds for
c                the variables
c     xhi      vector of length n containing the upper bounds for
c                the variables
c     xscale   vector of length n containing the diagonal scaling
c                matrix for the variables
c     xc       vector of length n containing the variable values
c                at the approximate solution
c     fc       vector of length m containing the residuals at the
c                approximate solution
c     fp       vector of length m containing the updated residual
c     xp       vector of length n containing the updated point
c     sc       vector of length n containing the last step taken
c     gc       vector of length n containing an estimate of
c                the gradient at the approximate solution
c     fjac     real m by n matrix containing an estimate of
c                the Jacobian at the approximate solution
c     mdim     leading dimension of fjac exactly as specified in
c                the dimension statement of the calling program
c     iactive  integer vector of length n indicating if xc(i) had
c                to be moved to an upper or lower bound
c     ipvt     vector of length n containing the permutation matrix
c                used in the QR factorization of the Jacobian at
c                the approximate solution
c     stpmax   real scalar containing the maximum allowed step size
c     delta    real scalar containing the trust region radius
c
c     required external routines :
c
c     rsdvalue   subroutine to evaluate residual function values
c     lsqwrite   subroutine to write out info about current status
c
c
      subroutine square (m,n,xlo,xhi,xc,fc,gc,fjac,mdim,
     &                      grdmin,rsdvalue,lsqwrite)
      implicit none
      include 'sizes.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'minima.i'
      integer maxlsq,maxrsd
      parameter (maxlsq=50)
      parameter (maxrsd=100)
      integer i,j,k
      integer icode,next
      integer niter,ncalls
      integer nactive
      integer m,n,mdim
      integer nbigstp,ndigit
      integer iactive(maxlsq)
      integer ipvt(maxlsq)
      real*8 amu,delta,epsfcn
      real*8 fcnorm,fpnorm
      real*8 gcnorm,ganorm
      real*8 precise,eps
      real*8 grdmin,stpnorm
      real*8 stpmax,stpmin
      real*8 rftol,faketol
      real*8 xtemp,stepsz
      real*8 sum,temp
      real*8 xc(maxlsq)
      real*8 xlo(maxlsq)
      real*8 xhi(maxlsq)
      real*8 fc(maxrsd),gc(maxlsq)
      real*8 fjac(mdim,maxlsq)
      real*8 xp(maxlsq),fp(maxrsd)
      real*8 ga(maxlsq),gs(maxlsq)
      real*8 sc(maxlsq),sa(maxlsq)
      real*8 ftemp(maxrsd)
      real*8 xscale(maxlsq)
      real*8 xsa(maxlsq)
      real*8 rdiag(maxlsq)
      real*8 qtf(maxlsq)
      real*8 work(maxlsq)
      logical done,first
      logical gauss,bigstp
      character*20 keyword
      character*120 record
      character*120 string
      external rsdvalue
      external lsqwrite
      external precise
c
c
c     check for too many variables or residuals
c
      if (n .gt. maxlsq) then
         write (iout,10)
   10    format (/,' SQUARE  --  Too many Parameters,',
     &              ' Increase the Value of MAXLSQ')
         return
      else if (m .gt. maxrsd) then
         write (iout,20)
   20    format (/,' SQUARE  --  Too many Residuals,',
     &              ' Increase the Value of MAXRSD')
         return
      end if
c
c     initialize various counters and status code
c
      niter = 0
      ncalls = 0
      nbigstp = 0
      done = .false.
c
c     setup the default tolerances and parameter values
c
      ndigit = 10
      eps = precise (2)
      eps = max(eps,10.0d0**(-ndigit))
      if (maxiter .eq. 0)  maxiter = 100
      if (iprint .lt. 0)  iprint = 1
      if (iwrite .lt. 0)  iwrite = 1
      if (fctmin .eq. 0.0d0)  fctmin = eps
      if (grdmin .eq. 0.0d0)  grdmin = eps**(1.0d0/3.0d0)
      epsfcn = sqrt(eps)
      delta = 0.0d0
      stpmax = 1000.0d0 * sqrt(dble(n))
      stpmin = eps**(2.0d0/3.0d0)
      rftol = eps**(2.0d0/3.0d0)
      faketol = 100.0d0 * eps
c
c     search each line of the keyword file for options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:7) .eq. 'FCTMIN ') then
            read (string,*,err=30,end=30)  fctmin
         else if (keyword(1:8) .eq. 'MAXITER ') then
            read (string,*,err=30,end=30)  maxiter
         else if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (string,*,err=30,end=30)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (string,*,err=30,end=30)  iwrite
         end if
   30    continue
      end do
c
c     check feasibility of variables and use bounds if needed
c
      nactive = 0
      do j = 1, n
         if (xc(j) .lt. xlo(j)) then
            xc(j) = xlo(j)
            iactive(j) = -1
         else if (xc(j) .gt. xhi(j)) then
            xc(j) = xhi(j)
            iactive(j) = 1
         else
            nactive = nactive + 1
            iactive(j) = 0
         end if
      end do
c
c     evaluate the function at the initial point
c
      ncalls = ncalls + 1
      call rsdvalue (m,n,xc,fc)
      fcnorm = 0.0d0
      do i = 1, m
         fcnorm = fcnorm + fc(i)**2
      end do
      fcnorm = 0.5d0 * fcnorm
c
c     evaluate the Jacobian at the initial point by finite
c     differences; replace loop with user routine if desired
c
      do j = 1, n
         stepsz = epsfcn * abs(xc(j))
         if (stepsz .lt. epsfcn)  stepsz = epsfcn
         if (xc(j) .lt. 0.0d0)  stepsz = -stepsz
         xtemp = xc(j)
         xc(j) = xtemp + stepsz
         ncalls = ncalls + 1
         call rsdvalue (m,n,xc,ftemp)
         xc(j) = xtemp
         do i = 1, m
            fjac(i,j) = (ftemp(i)-fc(i)) / stepsz
         end do
      end do
c
c     compute More's adaptive variable scale factors
c
      do j = 1, n
         temp = 0.0d0
         do i = 1, m
            temp = temp + fjac(i,j)**2
         end do
         xscale(j) = sqrt(temp)
         if (xscale(j) .eq. 0.0d0)  xscale(j) = 1.0d0
      end do
c
c     compute the total gradient vector for all variables
c
      do j = 1, n
         gc(j) = 0.0d0
         do i = 1, m
            gc(j) = gc(j) + fjac(i,j)*fc(i)
         end do
      end do
c
c     compute the norm of the scaled total gradient
c     and the scaled gradient for active variables
c
      gcnorm = 0.0d0
      ganorm = 0.0d0
      do j = 1, n
         gs(j) = gc(j) * max(abs(xc(j)),1.0d0/xscale(j))
         gcnorm = gcnorm + gs(j)**2
         if (iactive(j) .eq. 0) then
            ganorm = ganorm + gs(j)**2
         end if
      end do
      gcnorm = sqrt(gcnorm/n)
      if (nactive .ne. 0)  ganorm = sqrt(ganorm/nactive)
c
c     print out information about initial conditions
c
      if (iprint .gt. 0) then
         write (iout,40)
   40    format (/,' Levenberg-Marquardt Nonlinear Least Squares :')
         write (iout,50)
   50    format (/,' LS Iter    F Value      Total G     Active G',
     &              '    N Active   F Calls',/)
         if (max(fcnorm,gcnorm) .lt. 10000000.0d0) then
            write (iout,60)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
   60       format (i6,3f13.4,2i10)
         else
            write (iout,70)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
   70       format (i6,3d13.4,2i10)
         end if
      end if
c
c     write out the parameters, derivatives and residuals
c
      if (iwrite .ne. 0)  call lsqwrite (niter,xc,gs,m,fc)
c
c     check stopping criteria at the initial point; test the
c     absolute function value and gradient norm for termination
c
      if (fcnorm .le. fctmin)  return
      if (ganorm .le. grdmin)  return
c
c     start of the main body of least squares iteration
c
   80 continue
      niter = niter + 1
c
c     repack the Jacobian to include only active variables
c
      if (nactive .ne. n) then
         k = 0
         do j = 1, n
            if (iactive(j) .ne. 0) then
               if (k .eq. 0)  k = j
            else
               if (k .ne. 0) then
                  do i = 1, m
                     fjac(i,k) = fjac(i,j)
                  end do
                  k = k + 1
               end if
            end if
         end do
      end if
c
c     repack scale factors and gradient for active variables
c
      k = 0
      do j = 1, n
         if (iactive(j) .eq. 0) then
            k = k + 1
            xsa(k) = xscale(j)
            ga(k) = gc(j)
         end if
      end do
c
c     compute the QR factorization of the Jacobian
c
      call qrfact (m,nactive,fjac,mdim,.true.,ipvt,rdiag,work)
c
c     compute the vector Q(transpose) * residuals
c
      do i = 1, m
         qtf(i) = fc(i)
      end do
      do j = 1, nactive
         if (fjac(j,j) .ne. 0.0d0) then
            sum = 0.0d0
            do i = j, m
               sum = sum + fjac(i,j)*qtf(i)
            end do
            temp = -sum / fjac(j,j)
            do i = j, m
               qtf(i) = qtf(i) + fjac(i,j)*temp
            end do
         end if
         fjac(j,j) = rdiag(j)
      end do
c
c     compute the Levenberg-Marquardt step
c
      icode = 6
      first = .true.
      dowhile (icode .ge. 4)
         call lmstep (nactive,ga,fjac,mdim,ipvt,xsa,qtf,stpmax,
     &                       delta,amu,first,sa,gauss)
c
c     unpack the step vector to include all variables
c
         k = 0
         do i = 1, n
            if (iactive(i) .ne. 0) then
               sc(i) = 0.0d0
            else
               k = k + 1
               sc(i) = sa(k)
            end if
         end do
c
c     check new point and update the trust region
c
         call trust (rsdvalue,m,n,xc,fcnorm,gc,fjac,mdim,ipvt,sc,sa,
     &               xscale,gauss,stpmax,delta,icode,xp,fc,fp,fpnorm,
     &               bigstp,ncalls,xlo,xhi,nactive,stpmin,rftol,faketol)
      end do
      if (icode .eq. 1)  done = .true.
c
c     update to the new variables and residuals
c
      do j = 1, n
         xc(j) = xp(j)
      end do
      do i = 1, m
         fc(i) = fp(i)
      end do
      fcnorm = fpnorm
c
c     update the active vs inactive status of the variables;
c     in a true active set strategy, at most one constraint
c     is added to the active set per iteration
c
      do j = 1, n
         if (iactive(j) .eq. 0) then
            if (abs(xc(j)-xlo(j)) .le. eps) then
               nactive = nactive - 1
               iactive(j) = -1
c              goto 110
            else if (abs(xc(j)-xhi(j)) .le. eps) then
               nactive = nactive - 1
               iactive(j) = 1
c              goto 110
            end if
         end if
      end do
   90 continue
c
c     evaluate the Jacobian at the new point using finite
c     differences; replace loop with user routine if desired
c
      do j = 1, n
         stepsz = epsfcn * max(abs(xc(j)),1.0d0/xscale(j))
         if (xc(j) .lt. 0.0d0)  stepsz = -stepsz
         xtemp = xc(j)
         xc(j) = xtemp + stepsz
         ncalls = ncalls + 1
         call rsdvalue (m,n,xc,ftemp)
         xc(j) = xtemp
         do i = 1, m
            fjac(i,j) = (ftemp(i)-fc(i)) / stepsz
         end do
      end do
c
c     compute More's adaptive variable scale factors
c
      do j = 1, n
         temp = 0.0d0
         do i = 1, m
            temp = temp + fjac(i,j)**2
         end do
         xscale(j) = max(xscale(j),sqrt(temp))
      end do
c
c     compute the total gradient vector for all variables
c
      do j = 1, n
         gc(j) = 0.0d0
         do i = 1, m
            gc(j) = gc(j) + fjac(i,j)*fc(i)
         end do
      end do
c
c     compute the norm of the scaled total gradient
c     and the scaled gradient for active variables
c
      gcnorm = 0.0d0
      ganorm = 0.0d0
      do j = 1, n
         gs(j) = gc(j) * max(abs(xc(j)),1.0d0/xscale(j))
         gcnorm = gcnorm + gs(j)**2
         if (iactive(j) .eq. 0) then
            ganorm = ganorm + gs(j)**2
         end if
      end do
      gcnorm = sqrt(gcnorm/n)
      if (nactive .ne. 0)  ganorm = sqrt(ganorm/nactive)
c
c     print out information about current iteration
c
      if (iprint.ne.0 .and. mod(niter,iprint).eq.0) then
         if (max(fcnorm,gcnorm) .lt. 10000000.0d0) then
            write (iout,100)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
  100       format (i6,3f13.4,2i10)
         else
            write (iout,110)  niter,fcnorm,gcnorm,ganorm,nactive,ncalls
  110       format (i6,3d13.4,2i10)
         end if
      end if
c
c     check stopping criteria at the new point; test the absolute
c     function value, the gradient norm and step for termination
c
      if (fcnorm .le. fctmin)  done = .true.
      if (ganorm .le. grdmin)  done = .true.
      stpnorm = 0.0d0
      do j = 1, n
         temp = max(abs(xc(j)),1.0d0/xscale(j))
         stpnorm = stpnorm + (sc(j)/temp)**2
      end do
      stpnorm = sqrt(stpnorm/n)
      if (stpnorm .le. stpmin)  done = .true.
c
c     check for inactive variables that can be made active;
c     in a true active set strategy, variables are released
c     one at a time at a minimum of the current active set
c
c     if (done) then
      if (nactive .ne. n) then
         do j = 1, n
            if (iactive(j).eq.-1 .and. gc(j).lt.0.0d0) then
               nactive = nactive + 1
               iactive(j) = 0
               done = .false.
c              goto 140
            else if (iactive(j).eq.1 .and. gc(j).gt.0.0d0) then
               nactive = nactive + 1
               iactive(j) = 0
               done = .false.
c              goto 140
            end if
         end do
  120    continue
      end if
c     end if
c
c     if still done, then normal termination has been achieved
c
      if (done) then
         write (iout,130)
  130    format (/,' SQUARE  --  Normal Termination of Least Squares')
c
c     check the limit on the number of iterations
c
      else if (niter .ge. maxiter) then
         done = .true.
         write (iout,140)
  140    format (/,' SQUARE  --  Maximum Number of Allowed Iterations')
c
c     check for termination due to relative function convergence
c
      else if (icode .eq. 2) then
         done = .true.
         write (iout,150)
  150    format (/,' SQUARE  --  Relative Function Convergence',
     &           //,' Both the scaled actual and predicted',
     &              ' reductions in the function',
     &           /,' are less than or equal to the relative',
     &              ' convergence tolerance')
c
c     check for termination due to false convergence
c
      else if (icode .eq. 3) then
         done = .true.
         write (iout,160)
  160    format (/,' SQUARE  --  Possible False Convergence',
     &           //,' The iterates appear to be converging to',
     &              ' a noncritical point due',
     &           /,' to bad gradient information, discontinuous',
     &              ' function, or stopping',
     &           /,' tolerances being too tight')
c
c     check for several consecutive maximum steps taken
c
      else if (bigstp) then
         nbigstp = nbigstp + 1
         if (nbigstp .eq. 5) then
            done = .true.
            write (iout,170)
  170       format (/,' SQUARE  --  Five Consecutive Maximum',
     &                 ' Length Steps',
     &              //,' Either the function is unbounded below,',
     &                 ' or has a finite',
     &              /,' asymptote in some direction, or STEPMAX',
     &                 ' is too small')
         end if
c
c     no reason to quit, so prepare to take another step
c
      else
         nbigstp = 0
      end if
c
c     write out the parameters, derivatives and residuals
c
      if (iwrite.ne.0 .and. mod(niter,iwrite).eq.0) then
         if (.not. done)  call lsqwrite (niter,xc,gs,m,fc)
      end if
c
c     continue with the next iteration if not finished
c
      if (.not. done)  goto 80
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine lmstep  --  finds Levenberg-Marquardt step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "lmstep" computes the Levenberg-Marquardt step during a
c     nonlinear least squares calculation; this version is based
c     upon ideas from the Minpack routine LMPAR together with
c     with the internal doubling strategy of Dennis and Schnabel
c
c     arguments and variables :
c
c     n        dimension of the problem
c     ga       real vector of length n containing the gradient
c                of the residual vector
c     a        real n by n array which on input contains in the full
c                upper triangle the upper triangle of the matrix r
c                resulting from the QR factorization of the Jacobian;
c                on output the full upper triangle is unaltered, and
c                the strict lower triangle contains the strict lower
c                triangle of the lower triangular matrix l which is
c                equal to the Cholesky factor of (j**t)*j + amu*xscale
c     lda      leading dimension of a exactly as specified in the
c                dimension statement of the calling program
c     ipvt     integer array of length n containing the pivoting
c                infomation from the QR factorization routine
c     xscale   real vector of length n containing the diagonal
c                scaling matrix for the variables
c     qtf      real vector containing the first n elements of
c                Q(transpose) * (the scaled residual vector)
c     amu      scalar containing an initial estimate of the
c                Levenberg-Marquardt parameter on input; on
c                output, amu contains the final estimate of
c                the Levenberg-Marquardt parameter
c     first    logical variable set true only if this is the first
c                call to this routine in this iteration
c     sa       real vector of length n containing the
c                Levenberg-Marquardt step
c     gnstep   real vector of length n containing the
c                Gauss-Newton step
c     gauss    logical variable which is true if the Gauss-Newton
c                step is acceptable, and is false otherwise
c     diag     vector of length n containing the diagonal elements
c                of the Cholesky factor of (j**t)*j + amu*xscale
c
c
      subroutine lmstep (n,ga,a,lda,ipvt,xscale,qtf,stpmax,
     &                       delta,amu,first,sa,gauss)
      implicit none
      integer maxlsq,maxrsd
      parameter (maxlsq=50)
      parameter (maxrsd=100)
      integer i,j,k,n
      integer lda,nsing
      integer ipvt(maxlsq)
      real*8 stpmax,delta,amu
      real*8 alow,alpha
      real*8 amulow,amuhi
      real*8 beta,small,sum
      real*8 deltap,gnleng
      real*8 phi,phip,phipi
      real*8 sgnorm,precise
      real*8 high,tiny
      real*8 stplen,temp
      real*8 ga(maxlsq)
      real*8 a(lda,maxlsq)
      real*8 xscale(maxlsq)
      real*8 qtf(maxlsq)
      real*8 sa(maxlsq)
      real*8 gnstep(maxlsq)
      real*8 diag(maxlsq)
      real*8 work1(maxlsq)
      real*8 work2(maxlsq)
      logical first,gauss
      logical done
      save deltap,nsing
      save phi,phip
      external precise
c
c
c     set smallest floating point magnitude and spacing
c
      tiny = precise (1)
      small = precise (2)
c
c     if initial trust region is not provided by the user,
c     compute and use the length of the Cauchy step given
c     by beta = norm2(r*trans(p)*d**(-2)*g)**2
c
      if (delta .eq. 0.0d0) then
         amu = 0.0d0
         do i = 1, n
            work1(i) = ga(i) / xscale(i)
         end do
         alpha = 0.0d0
         do i = 1, n
            alpha = alpha + work1(i)**2
         end do
         beta = 0.0d0
         do i = 1, n
            temp = 0.0d0
            do j = i, n
               k = ipvt(j)
               temp = temp + a(i,j)*ga(k)/xscale(k)**2
            end do
            beta = beta + temp**2
         end do
         if (beta .le. tiny) then
            delta = alpha * sqrt(alpha)
         else
            delta = alpha * sqrt(alpha)/beta
         end if
         delta = min(delta,stpmax)
      end if
c
c     the following is done only on the first time through
c     this iteration: (1) compute the Gauss-Newton step;
c     if the Jacobian is rank-deficient, obtain a least
c     squares solution, (2) compute the length of the scaled
c     Gauss-Newton step, (3) compute the norm of the scaled
c     gradient used in computing an upper bound for "amu"
c
      if (first) then
         nsing = n
         do j = 1, n
            if (a(j,j).eq.0.0d0 .and. nsing.eq.n)  nsing = j - 1
            if (nsing .lt. n)  work1(j) = 0.0d0
         end do
         work1(nsing) = qtf(nsing) / a(nsing,nsing)
         do j = nsing-1, 1, -1
            sum = 0.0d0
            do i = j+1, nsing
               sum = sum + a(j,i)*work1(i)
            end do
            work1(j) = (qtf(j)-sum) / a(j,j)
         end do
         do j = 1, n
            gnstep(ipvt(j)) = -work1(j)
         end do
c
c     find the length of scaled Gauss-Newton step
c
         do j = 1, n
            work1(j) = xscale(j) * gnstep(j)
         end do
         gnleng = 0.0d0
         do j = 1, n
            gnleng = gnleng + work1(j)**2
         end do
         gnleng = sqrt(gnleng)
c
c     find the length of the scaled gradient
c
         do j = 1, n
            work1(j) = ga(j) / xscale(j)
         end do
         sgnorm = 0.0d0
         do j = 1, n
            sgnorm = sgnorm + work1(j)**2
         end do
         sgnorm = sqrt(sgnorm)
      end if
c
c     set the bounds on the computed step
c
      high = 1.5d0
      alow = 0.75d0
c
c     check to see if the Gauss-Newton step is acceptable
c
      if (gnleng .le. high*delta) then
         gauss = .true.
         do j = 1, n
            sa(j) = gnstep(j)
         end do
         amu = 0.0d0
         delta = min(delta,gnleng)
c
c     the Gauss-Newton step is rejected, find a nontrivial step;
c     first compute a starting value of "amu" if previous step
c     was not a Gauss-Newton step
c
      else
         gauss = .false.
         if (amu .gt. 0.0d0)
     &      amu = amu - ((phi+deltap)/delta)*(((deltap-delta)+phi)/phip)
         phi = gnleng - delta
c
c     if the Jacobian is not rank deficient, the Newton step
c     provides a lower bound for "amu"; else set bound to zero
c
         if (nsing .eq. n) then
            if (first) then
               first = .false.
               do j = 1, n
                  k = ipvt(j)
                  work1(j) = gnstep(k) * xscale(k)**2
               end do
c
c     obtain trans(r**-1)*(trans(p)*s) by solving the
c     system of equations trans(r)*work1 = work1
c
               work1(n) = work1(n) / a(n,n)
               do j = n-1, 1, -1
                  sum = 0.0d0
                  do i = j+1, n
                     sum = sum + a(j,i)*work1(i)
                  end do
                  work1(j) = (work1(j)-sum) / a(j,j)
               end do
               phipi = 0.0d0
               do j = 1, n
                  phipi = phipi - work1(j)**2
               end do
               phipi = phipi / gnleng
            end if
            amulow = -phi / phipi
         else
            first = .false.
            amulow = 0.0d0
         end if
         amuhi = sgnorm / delta
c
c     iterate until a satisfactory "amu" is generated
c
         done = .false.
         dowhile (.not. done)
            if (amu.lt.amulow .or. amu.gt.amuhi) then
               amu = max(sqrt(amulow*amuhi),0.001d0*amuhi)
            end if
            temp = sqrt(amu)
            do j = 1, n
               work1(j) = temp * xscale(j)
            end do
c
c     solve the damped least squares system for the value of the
c     Levenberg-Marquardt step using More's Minpack technique
c
            call qrsolve (n,a,lda,ipvt,work1,qtf,sa,diag,work2)
            do j = 1, n
               sa(j) = -sa(j)
            end do
            do j = 1, n
               work2(j) = xscale(j) * sa(j)
            end do
            stplen = 0.0d0
            do j = 1, n
               stplen = stplen + work2(j)**2
            end do
            stplen = sqrt(stplen)
            phi = stplen - delta
            do j = 1, n
               k = ipvt(j)
               work1(j) = xscale(k) * work2(k)
            end do
            do j = 1, n
               if (abs(diag(j)) .ge. tiny) then
                  work1(j) = work1(j) / diag(j)
               end if
               if (j .lt. n) then
                  do i = j+1, n
                     work1(i) = work1(i) - work1(j)*a(i,j)
                  end do
               end if
            end do
            phip = 0.0d0
            do j = 1, n
               phip = phip - work1(j)**2
            end do
            phip = phip / stplen
c
c     check to see if the step is acceptable; if not,
c     update amulow, amuhi and amu for next iteration
c
            if ((stplen.ge.alow*delta.and.stplen.le.high*delta)
     &            .or. (amuhi-amulow).le.small) then
               done = .true.
            else
               amulow = max(amulow,amu-(phi/phip))
               if (phi .lt. 0.0d0)  amuhi = amu
               amu = amu - (stplen/delta)*(phi/phip)
            end if
         end do
      end if
      deltap = delta
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine qrfact  --  QR factorization of a matrix  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "qrfact" performs Householder transformations with column
c     pivoting (optional) to compute a QR factorization of the
c     m by n matrix a; the routine determines an orthogonal
c     matrix q, a permutation matrix p, and an upper trapezoidal
c     matrix r with diagonal elements of nonincreasing magnitude,
c     such that a*p = q*r; the Householder transformation for
c     column k, k = 1,2,...,min(m,n), is of the form
c
c               i - (1/u(k))*u*u(transpose)
c
c     where u has zeros in the first k-1 positions
c
c     arguments and variables :
c
c     m        positive integer input variable set to
c                the number of rows of a
c     n        positive integer input variable set to
c                the number of columns of a
c     a        m by n array; on input a contains the matrix for
c                which the QR factorization is to be computed; on
c                output the strict upper trapezoidal part contains
c                the strict upper trapezoidal part of r, the lower
c                trapezoidal part contains a factored form of q
c                (non-trivial elements of u vectors described above)
c     lda      positive integer input variable not less than m
c                which specifies the leading dimension of array a
c     pivot    logical input variable; if pivot is set true,
c                then column pivoting is enforced; if pivot is
c                set false, then no column pivoting is done
c     ipvt     integer output array which defines the permutation
c                matrix p such that a*p = q*r; column j of p is
c                column ipvt(j) of the identity matrix
c     rdiag    an output array of length n which contains
c                the diagonal elements of r
c
c
      subroutine qrfact (m,n,a,lda,pivot,ipvt,rdiag,work)
      implicit none
      integer i,j,k,jmax
      integer m,n,lda
      integer minmn,itemp
      integer ipvt(n)
      real*8 aknorm,temp
      real*8 rdiag(n),work(n)
      real*8 a(lda,n)
      logical pivot
c
c
c     find the initial column norms and initialize some arrays
c
      do j = 1, n
         temp = 0.0d0
         do i = 1, m
            temp = temp + a(i,j)**2
         end do
         rdiag(j) = sqrt(temp)
         work(j) = rdiag(j)
         if (pivot)  ipvt(j) = j
      end do
c
c     reduce the matrix with Householder transformations
c
      minmn = min(m,n)
      do k = 1, minmn
c
c     bring the column of largest norm into the pivot position
c
         if (pivot) then
            jmax = k
            do j = k, n
               if (rdiag(j) .gt. rdiag(jmax))  jmax = j
            end do
            if (jmax .ne. k) then
               do i = 1, m
                  temp = a(i,k)
                  a(i,k) = a(i,jmax)
                  a(i,jmax) = temp
               end do
               rdiag(jmax) = rdiag(k)
               work(jmax) = work(k)
               itemp = ipvt(k)
               ipvt(k) = ipvt(jmax)
               ipvt(jmax) = itemp
            end if
         end if
c
c     compute the Householder transformation to reduce the
c     k-th column of a to a multiple of the k-th unit vector
c
         aknorm = 0.0d0
         do i = k, m
            aknorm = aknorm + a(i,k)**2
         end do
         aknorm = sqrt(aknorm)
         if (aknorm .ne. 0.0d0) then
            if (a(k,k) .lt. 0.0d0)  aknorm = -aknorm
            do i = k, m
               a(i,k) = a(i,k) / aknorm
            end do
            a(k,k) = a(k,k) + 1.0d0
c
c     apply the transformation to the remaining columns
c     and update the column norms
c
            if (n .ge. k+1) then
               do j = k+1, n
                  temp = 0.0d0
                  do i = k, m
                     temp = temp + a(i,k)*a(i,j)
                  end do
                  temp = temp / a(k,k)
                  do i = k, m
                     a(i,j) = a(i,j) - temp*a(i,k)
                  end do
                  if (pivot .and. rdiag(j).ne.0.0d0) then
                     temp = a(k,j) / rdiag(j)
                     if (abs(temp) .lt. 1.0d0) then
                        rdiag(j) = rdiag(j) * sqrt(1.0d0-temp**2)
                     else
                        temp = 0.0d0
                        do i = k+1, m
                           temp = temp + a(i,j)**2
                        end do
                        rdiag(j) = sqrt(temp)
                        work(j) = rdiag(j)
                     end if
                  end if
               end do
            end if
         end if
         rdiag(k) = -aknorm
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine qrsolve  --  get least squares via QR factors  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "qrsolve" solves a*x=b and d*x=0 in the least squares sense;
c     normally used in combination with routine "qrfact" to solve
c     least squares problems
c
c     arguments and variables :
c
c     n        number of rows and columns in the matrix r
c     r        an n by n array containing the upper triangular
c                matrix r; on output the full triangle is unaltered,
c                and the strict lower triangle contains the transpose
c                of the strict upper triangular matrix s
c     ldr      leading dimension of r exactly as specified in
c                the dimension statement of the calling program
c     ipvt     vector of length n which defines the permutation
c                matrix p such that a*p = q*r; column j of p is
c                column ipvt(j) of the identity matrix
c     diag     vector of length n containing the diagonal elements
c                of the matrix d
c     qtb      vector of length n containing the first n elements
c                of the vector q(transpose)*b
c     x        vector of length n containing the least squares
c                solution of the systems a*x = b, d*x = 0
c     sdiag    vector of length n containing the diagonal elements
c                of the upper triangular matrix s
c
c
      subroutine qrsolve (n,r,ldr,ipvt,diag,qtb,x,sdiag,work)
      implicit none
      integer i,j,k,jj
      integer nsing,n,ldr
      integer ipvt(n)
      real*8 sine,cosine
      real*8 tangent
      real*8 cotangent
      real*8 qtbpj,temp
      real*8 diag(n),qtb(n)
      real*8 x(n),sdiag(n)
      real*8 work(n)
      real*8 r(ldr,n)
c
c
c     copy r and (q transpose)*b to preserve input and
c     initialize s; in particular, save the diagonal
c     elements of r in x
c
      do j = 1, n-1
         do k = j+1, n
            r(k,j) = r(j,k)
         end do
      end do
      do j = 1, n
         x(j) = r(j,j)
         work(j) = qtb(j)
      end do
c
c     eliminate the diagonal matrix d using a Givens rotation
c
      do j = 1, n
c
c     prepare the row of d to be eliminated, locating
c     the diagonal element using p from the QR factorization
c
         jj = ipvt(j)
         if (diag(jj) .ne. 0.0d0) then
            do k = j, n
               sdiag(k) = 0.0d0
            end do
            sdiag(j) = diag(jj)
c
c     the transformations to eliminate the row of d modify
c     only a single element of (q transpose)*b beyond the
c     first n, which is initially zero
c
            qtbpj = 0.0d0
            do k = j, n
c
c     determine a Givens rotation which eliminates the
c     appropriate element in the current row of d
c
               if (sdiag(k) .ne. 0.0d0) then
                  if (abs(r(k,k)) .lt. abs(sdiag(k))) then
                     cotangent = r(k,k) / sdiag(k)
                     sine = 0.5d0 / sqrt(0.25d0+0.25d0*cotangent**2)
                     cosine = sine * cotangent
                  else
                     tangent = sdiag(k) / r(k,k)
                     cosine = 0.5d0 / sqrt(0.25d0+0.25d0*tangent**2)
                     sine = cosine * tangent
                  end if
c
c     compute the modified diagonal element of r
c     and the modified element of ((q transpose)*b,0)
c
                  r(k,k) = cosine*r(k,k) + sine*sdiag(k)
                  temp = cosine*work(k) + sine*qtbpj
                  qtbpj = -sine*work(k) + cosine*qtbpj
                  work(k) = temp
c
c     accumulate the tranformation in the row of s
c
                  if (n .ge. k+1) then
                     do i = k+1, n
                        temp = cosine*r(i,k) + sine*sdiag(i)
                        sdiag(i) = -sine*r(i,k) + cosine*sdiag(i)
                        r(i,k) = temp
                     end do
                  end if
               end if
            end do
         end if
c
c     store the diagonal element of s and restore
c     the corresponding diagonal element of r
c
         sdiag(j) = r(j,j)
         r(j,j) = x(j)
      end do
c
c     solve the triangular system for z; if the system
c     is singular, then obtain a least squares solution
c
      nsing = n
      do j = 1, n
         if (sdiag(j).eq.0.0d0 .and. nsing.eq.n)  nsing = j - 1
         if (nsing .lt. n)  work(j) = 0.0d0
      end do
      if (nsing .ge. 1) then
         do k = 1, nsing
            j = nsing - k + 1
            temp = 0.0d0
            if (nsing .ge. j+1) then
               do i = j+1, nsing
                  temp = temp + r(i,j)*work(i)
               end do
            end if
            work(j) = (work(j)-temp) / sdiag(j)
         end do
      end if
c
c     permute the components of z back to components of x
c
      do j = 1, n
         k = ipvt(j)
         x(k) = work(j)
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine trust  --  update the model trust region  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "trust" updates the model trust region for a nonlinear
c     least squares calculation; this version is based on the
c     ideas found in NL2SOL and in Dennis and Schnabel's book
c
c     arguments and variables :
c
c     m        number of functions
c     n        number of variables
c     xc       vector of length n containing the current iterate
c     fcnorm   real scalar containing the norm of f(xc)
c     gc       vector of length n containing the gradient at xc
c     a        real m by n matrix containing the upper triangular
c                matrix r from the QR factorization of the current
c                Jacobian in the upper triangle
c     lda      leading dimension of A exactly as specified in
c                the dimension statement of the calling program
c     ipvt     integer vector of length n containing the permutation
c                matrix from QR factorization of the Jacobian
c     sc       vector of length n containing the Newton step
c     sa       vector of length n containing current step
c     xscale   vector of length n containing the diagonal
c                scaling matrix for x
c     gauss    logical variable equal is true when the Gauss-Newton
c                step is taken
c     stpmax   maximum allowable step size
c     delta    trust region radius with value retained between calls
c     icode    return code set upon exit
c                0  means xp accepted as next iterate, delta
c                     is trust region for next iteration
c                1  means the algorithm was unable to find a
c                     satisfactory xp sufficiently distinct from xc
c                2  means both the scaled actual and predicted
c                     function reductions are smaller than rftol
c                3  means that false convergence is detected
c                4  means fpnorm is too large, current iteration is
c                     continued with a new, reduced trust region
c                5  means fpnorm is sufficiently small, but the
c                     chance of taking a longer successful step
c                     seems good that the current iteration is to
c                     be continued with a new, doubled trust region
c     xpprev   vector of length n containing the value of xp
c                at the previous call within this iteration
c     fpprev   vector of length m containing f(xpprev)
c     xp       vector of length n containing the new iterate
c     fp       vector of length m containing the functions at xp
c     fpnorm   scalar containing the norm of f(xp)
c     bigstp   logical variable; true, if maximum step length was taken,
c                false  otherwise
c     ncalls   number of function evaluations used
c     xlo      vector of length n containing the lower bounds
c     xhi      vector of length n containing the upper bounds
c     nactive  number of columns in the active Jacobian
c
c     required external routines :
c
c     rsdvalue   subroutine to evaluate residual function values
c
c
      subroutine trust (rsdvalue,m,n,xc,fcnorm,gc,a,lda,ipvt,sc,sa,
     &                  xscale,gauss,stpmax,delta,icode,xp,fc,fp,
     &                  fpnorm,bigstp,ncalls,xlo,xhi,nactive,stpmin,
     &                  rftol,faketol)
      implicit none
      integer maxlsq,maxrsd
      parameter (maxlsq=50)
      parameter (maxrsd=100)
      integer i,j,k,m,n
      integer lda,icode
      integer ncalls,nactive
      integer ipvt(maxlsq)
      real*8 fcnorm,stpmax
      real*8 delta,fpnorm,fpnrmp
      real*8 reduce,predict
      real*8 rellen,slope,tiny
      real*8 stplen,stpmin
      real*8 rftol,faketol
      real*8 alpha,temp,precise
      real*8 xc(maxlsq),gc(maxlsq)
      real*8 a(lda,maxlsq)
      real*8 sc(maxlsq),sa(maxlsq)
      real*8 xp(maxlsq)
      real*8 fc(maxrsd),fp(maxrsd)
      real*8 xlo(maxlsq),xhi(maxlsq)
      real*8 xscale(maxlsq)
      real*8 xpprev(maxlsq)
      real*8 fpprev(maxrsd)
      logical gauss,bigstp
      logical feas,ltemp
      save xpprev,fpprev
      save fpnrmp
      external rsdvalue
      external precise
c
c
c     set value of alpha, logicals and step length
c
      alpha = 0.0001d0
      bigstp = .false.
      feas = .true.
      stplen = 0.0d0
      do i = 1, n
         stplen = stplen + (xscale(i)*sc(i))**2
      end do
      stplen = sqrt(stplen)
c
c     compute new trial point and new function values
c
      do i = 1, n
         xp(i) = xc(i) + sc(i)
         if (xp(i) .gt. xhi(i)) then
            sc(i) = xhi(i) - xc(i)
            xp(i) = xhi(i)
            feas = .false.
         else if (xp(i) .lt. xlo(i)) then
            sc(i) = xlo(i) - xc(i)
            xp(i) = xlo(i)
            feas = .false.
         end if
      end do
      ncalls = ncalls + 1
      call rsdvalue (m,n,xp,fp)
      fpnorm = 0.0d0
      do i = 1, m
         fpnorm = fpnorm + fp(i)**2
      end do
      fpnorm = 0.5d0 * fpnorm
      reduce = fpnorm - fcnorm
      slope = 0.0d0
      do i = 1, n
         slope = slope + gc(i)*sc(i)
      end do
      if (icode .ne. 5)  fpnrmp = 0.0d0
c
c     internal doubling no good; reset to previous and quit
c
      if (icode.eq.5 .and.
     &     ((fpnorm.ge.fpnrmp).or.(reduce.gt.alpha*slope))) then
         icode = 0
         do i = 1, n
            xp(i) = xpprev(i)
         end do
         do i = 1, m
            fp(i) = fpprev(i)
         end do
         fpnorm = fpnrmp
         delta = 0.5d0 * delta
c
c     fpnorm is too large; the step is unacceptable
c
      else if (reduce .ge. alpha*slope) then
         rellen = 0.0d0
         do i = 1, n
            temp = abs(sc(i))/max(abs(xp(i)),1.0d0/xscale(i))
            rellen = max(rellen,temp)
         end do
c
c     magnitude of (xp-xc) is too small, end the global step
c
         if (rellen .lt. stpmin) then
            icode = 1
            do i = 1, n
               xp(i) = xc(i)
            end do
            do i = 1, m
               fp(i) = fc(i)
            end do
c
c     quadratic interpolation step; reduce delta and continue
c
         else
            icode = 4
            tiny = precise (1)
            if (abs(reduce-slope) .gt. tiny) then
               temp = -slope*stplen / (2.0d0*(reduce-slope))
            else
               temp = -slope*stplen / 2.0d0
            end if
            if (temp .lt. 0.1d0*delta) then
               delta = 0.1d0 * delta
            else if (temp .gt. 0.5d0*delta) then
               delta = 0.5d0 * delta
            else
               delta = temp
            end if
         end if
c
c     fpnorm is sufficiently small; the step is acceptable compute
c     the predicted reduction as predict = g(T)*s + (1/2)*s(T)*h*s
c     with h = p * r**t * r * p**t
c
      else
         predict = slope
         do i = 1, nactive
            k = ipvt(i)
            temp = 0.0d0
            do j = i, nactive
               temp = temp + sa(k)*a(i,j)
            end do
            predict = predict + 0.5d0*temp**2
         end do
         ltemp = (abs(predict-reduce) .le. 0.1d0*abs(reduce))
c
c     if reduce and predict agree to within relative error of 0.1
c     or if negative curvature is indicated, and a longer step is
c     possible and delta has not been decreased this iteration,
c     then double trust region and continue global step
c
         if (icode.ne.4 .and. (ltemp.or.(reduce.le.slope)) .and. feas
     &        .and. .not.gauss .and. (delta.le.0.99d0*stpmax)) then
            icode = 5
            do i = 1, n
               xpprev(i) = xp(i)
            end do
            do i = 1, m
               fpprev(i) = fp(i)
            end do
            fpnrmp = fpnorm
            delta = min(2.0d0*delta,stpmax)
c
c     accept the point; choose new trust region for next iteration
c
         else
            icode = 0
            if (stplen .gt. 0.99d0*stpmax)  bigstp = .true.
            if (reduce .ge. 0.1d0*predict) then
               delta = 0.5d0 * delta
            else if (reduce .le. 0.75d0*predict) then
               delta = min(2.0d0*delta,stpmax)
            end if
         end if
c
c     check relative function convergence and false convergence
c
         if (reduce .le. 2.0d0*predict) then
            if (abs(reduce).le.rftol*abs(fcnorm) .and.
     &          abs(predict).le.rftol*abs(fcnorm)) then
               icode = 2
            end if
         else
            rellen = 0.0d0
            do i = 1, n
               temp = abs(sc(i))/max(abs(xp(i)),1.0d0/xscale(i))
               rellen = max(rellen,temp)
            end do
            if (rellen .lt. faketol)  icode = 3
         end if
      end if
      return
      end

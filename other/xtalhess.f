c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2004 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program xtalhess  --  crystal lattice crystal Hessian  ##
c     ##                                                         ##
c     #############################################################
c
c
      program xtalhess
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'math.i'
      include 'scales.i'
      integer i,j,nvar
      real*8 minimum,mean
      real*8 eps,xxold
      real*8 xtalmin1
      real*8 xx(maxvar)
      real*8 glat(maxvar)
      real*8 glat0(maxvar)
      real*8 glat1(maxvar)
      real*8 derivs(3,maxatm)
      real*8 hess(6,6)
      real*8 eval(6)
      real*8 evec(6,6)
      real*8 temp1(6),temp2(6)
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     write out the initial values of the lattice parameters
c
      write (iout,10)  xbox,ybox,zbox,alpha,beta,gamma
   10 format (/,' Initial Lattice Dimensions :    a   ',f12.4,
     &        /,'                                 b   ',f12.4,
     &        /,'                                 c   ',f12.4,
     &        /,'                                Alpha',f12.4,
     &        /,'                                Beta ',f12.4,
     &        /,'                                Gamma',f12.4)
c
c     set scale factors to apply to optimization variables
c
      set_scale = .true.
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         scale(nvar) = 12.0d0 * xbox
         nvar = nvar + 1
         scale(nvar) = 12.0d0 * ybox
         nvar = nvar + 1
         scale(nvar) = 12.0d0 * zbox
      end do
      scale(nvar+1) = 1.0d0
      scale(nvar+2) = 1.0d0
      scale(nvar+3) = 1.0d0
      scale(nvar+4) = 1.0d0
      scale(nvar+5) = 1.0d0
      scale(nvar+6) = 1.0d0
      scale(nvar+1) = 4.0d0 * sqrt(xbox)
      scale(nvar+2) = 4.0d0 * sqrt(ybox)
      scale(nvar+3) = 4.0d0 * sqrt(zbox)
      scale(nvar+4) = 0.02d0 * sqrt(volbox)
      scale(nvar+5) = 0.02d0 * sqrt(volbox)
      scale(nvar+6) = 0.02d0 * sqrt(volbox)
c
c     compute the fractional coordinates for each atom
c
      call lattice
      do i = 1, n
         j = 3*i - 3
         xx(j+1) = recip(1,1)*x(i) + recip(2,1)*y(i) + recip(3,1)*z(i)
         xx(j+2) = recip(1,2)*x(i) + recip(2,2)*y(i) + recip(3,2)*z(i)
         xx(j+3) = recip(1,3)*x(i) + recip(2,3)*y(i) + recip(3,3)*z(i)
      end do
c
c     scale the fractional coordinates and lattice parameters
c
      nvar = 3 * n
      do i = 1, nvar
         xx(i) = xx(i) * scale(i)
      end do
      xx(nvar+1) = xbox * scale(nvar+1)
      xx(nvar+2) = ybox * scale(nvar+2)
      xx(nvar+3) = zbox * scale(nvar+3)
      xx(nvar+4) = alpha * scale(nvar+4)
      xx(nvar+5) = beta * scale(nvar+5)
      xx(nvar+6) = gamma * scale(nvar+6)
c
c     get lattice gradient at initial lattice parameters
c
      eps = 0.0001d0
      minimum = xtalmin1 (xx,glat)
      write (iout,20)  0,(glat(i),i=nvar+1,nvar+6)
   20 format (/,i6,6f12.8)
c
c     compute the numerical crystal lattice Hessian elements
c
      do j = 1, 6
         eps = 0.0001d0
         if (j .ge. 4)  eps = eps * radian
         xxold = xx(nvar+j)
         xx(nvar+j) = xx(nvar+j) - 0.5d0*eps
         xbox = xbox - 0.5d0*eps
         minimum = xtalmin1 (xx,glat0)
         xx(nvar+j) = xxold
         write (iout,30)  1,(glat0(i),i=nvar+1,nvar+6)
   30    format (/,i6,6f12.8)
         xx(nvar+j) = xx(nvar+j) + 0.5d0*eps
         xbox = xbox - 0.5d0*eps
         minimum = xtalmin1 (xx,glat1)
         xx(nvar+j) = xxold
         write (iout,40)  2,(glat1(i),i=nvar+1,nvar+6)
   40    format (i6,6f12.8)
         do i = 1, 6
            hess(i,j) = (glat1(nvar+i)-glat0(nvar+i)) / eps
         end do
      end do
c
c     output the raw Hessian elements and a symmetrized version
c
      write (iout,50)
   50 format ()
      do j = 1, 6
         write (iout,60)  (hess(i,j),i=1,6)
   60    format (6f12.6)
      end do
      do i = 1, 5
         do j = i+1, 6
            mean = 0.5d0 * (hess(i,j) + hess(j,i))
            hess(i,j) = mean
            hess(j,i) = mean
         end do
      end do
      write (iout,70)
   70 format ()
      do j = 1, 6
         write (iout,80)  (hess(i,j),i=1,6)
   80    format (6f12.6)
      end do
c
c     diagonalize the matrix of lattice Hessian elements
c
      call jacobi (6,6,hess,eval,evec,temp1,temp2)
      write (iout,90)  (eval(i),i=1,6)
   90 format (/,6f12.6)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function xtalmin1  --  energy and gradient for lattice  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "xtalmin1" is a service routine that computes the energy and
c     gradient with respect to fractional coordinates and lattice
c     dimensions for a crystal energy minimization
c
c
      function xtalmin1 (xx,g)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'math.i'
      include 'scales.i'
      integer i,j
      real*8 xtalmin1,energy
      real*8 e,e0,old,eps
      real*8 xx(maxvar)
      real*8 g(maxvar)
      real*8 xf(maxatm)
      real*8 yf(maxatm)
      real*8 zf(maxatm)
      real*8 derivs(3,maxatm)
c
c
c     translate optimization variables to fractional coordinates
c
      do i = 1, n
         j = 3*i - 3
         xf(i) = xx(j+1) / scale(j+1)
         yf(i) = xx(j+2) / scale(j+2)
         zf(i) = xx(j+3) / scale(j+3)
      end do
c
c     translate optimization variables to lattice parameters
c
      xbox = xx(3*n+1) / scale(3*n+1)
      ybox = xx(3*n+2) / scale(3*n+2)
      zbox = xx(3*n+3) / scale(3*n+3)
      alpha = xx(3*n+4) / scale(3*n+4)
      beta = xx(3*n+5) / scale(3*n+5)
      gamma = xx(3*n+6) / scale(3*n+6)
c
c     update current atomic coordinates based on optimization values
c
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
c
c     find energy and fractional coordinates deriviatives
c
      call gradient (e,derivs)
      xtalmin1 = e
      do i = 1, n
         j = 3*i - 3
         g(j+1) = derivs(1,i)*lvec(1,1) + derivs(2,i)*lvec(1,2)
     &               + derivs(3,i)*lvec(1,3)
         g(j+2) = derivs(1,i)*lvec(2,1) + derivs(2,i)*lvec(2,2)
     &               + derivs(3,i)*lvec(2,3)
         g(j+3) = derivs(1,i)*lvec(3,1) + derivs(2,i)*lvec(3,2)
     &               + derivs(3,i)*lvec(3,3)
      end do
c
c     find derivative with respect to lattice a-axis length
c
      eps = 0.0001d0
      old = xbox
      xbox = xbox - 0.5d0*eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e0 = energy ()
      xbox = xbox + eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e = energy ()
      g(3*n+1) = (e - e0) / eps
      xbox = old
c
c     find derivative with respect to lattice b-axis length
c
      old = ybox
      ybox = ybox - 0.5d0*eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e0 = energy ()
      ybox = ybox + eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e = energy ()
      g(3*n+2) = (e - e0) / eps
      ybox = old
c
c     find derivative with respect to lattice c-axis length
c
      old = zbox
      zbox = zbox - 0.5d0*eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e0 = energy ()
      zbox = zbox + eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e = energy ()
      g(3*n+3) = (e - e0) / eps
      zbox = old
c
c     find derivative with respect to lattice alpha angle
c
      eps = eps * radian
      old = alpha
      alpha = alpha - 0.5d0*eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e0 = energy ()
      alpha = alpha + eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e = energy ()
      g(3*n+4) = (e - e0) / eps
      alpha = old
c
c     find derivative with respect to lattice beta angle
c
      old = beta
      beta = beta - 0.5d0*eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e0 = energy ()
      beta = beta + eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e = energy ()
      g(3*n+5) = (e - e0) / eps
      beta = old
c
c     find derivative with respect to lattice gamma angle
c
      old = gamma
      gamma = gamma - 0.5d0*eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e0 = energy ()
      gamma = gamma + eps
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
      e = energy ()
      g(3*n+6) = (e - e0) / eps
      gamma = old
c
c     revert to the original atomic coordinate values
c
      call lattice
      do i = 1, n
         x(i) = xf(i)*lvec(1,1) + yf(i)*lvec(2,1) + zf(i)*lvec(3,1)
         y(i) = xf(i)*lvec(1,2) + yf(i)*lvec(2,2) + zf(i)*lvec(3,2)
         z(i) = xf(i)*lvec(1,3) + yf(i)*lvec(2,3) + zf(i)*lvec(3,3)
      end do
c
c     apply scale factors to the coordinate and lattice gradient
c
      do i = 1, 3*n+6
         g(i) = g(i) / scale(i)
      end do
      return
      end

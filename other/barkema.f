c
c
c     ###########################################################
c     ##                                                       ##
c     ##  program barkema  --  ART code to find saddle points  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "barkema" is a program to find the transition states on the
c     Muller-Brown test potential surface via the Barkema method
c
c
      program barkema
      implicit none
      include 'sizes.i'
      include 'minima.i'
      include 'output.i'
      integer nvar
      real*8 art1,fmuller
      real*8 minimum,grdmin
      real*8 xx(maxopt),random
      real*8 value,grad(maxopt)
      real*8 xref(2)
      common /art/ xref
      external art1,optsave
c
c
c     set the coordinates of Muller-Brown minimum to be used
c
      call initial
      nvar = 2
c     xx(1) = -0.558d0
c     xx(2) = 1.442d0
      xx(1) = 0.623d0
      xx(2) = 0.028d0
c     xx(1) = -0.050d0
c     xx(2) = 0.467d0
c
c     write the starting point and its function value
c
      value = fmuller(xx(1),xx(2))
      call gmuller (xx(1),xx(2),grad(1),grad(2))
      write (*,10)  xx(1),xx(2),value,grad(1),grad(2)
   10 format (/,' Start Point :     ',2f16.4,
     &        //,' Function Value :  ',f16.4,
     &        //,' Gradients :       ',2f16.4)
c
c     store the coordinates of the original point
c
      xref(1) = xx(1)
      xref(2) = xx(2)
c
c     generate a random perturbation in original values
c
      xx(1) = xref(1) + 0.1d0 * random ()
      xx(2) = xref(2) + 0.1d0 * random ()
c
c     termination criterion on the RMS gradient
c
      iprint = 1
      coordtype = 'NONE'
c     fctmin = 0.0001d0
      grdmin = 0.0001d0
c
c     make the call to the optimization routine
c
      call vmetric (nvar,xx,minimum,grdmin,art1,optsave)
c     call lmqn (nvar,xx,minimum,grdmin,art1,optsave)
c
c     write out final function value and gradient
c
      write (*,20)  minimum
   20 format (/,' Final Barkema Norm Value :  ',f16.4)
c
c     write the function value and saddle point coordinates
c
      value = fmuller(xx(1),xx(2))
      call gmuller (xx(1),xx(2),grad(1),grad(2))
      write (*,30)  xx(1),xx(2),value,grad(1),grad(2)
   30 format (/,' Saddle Point :    ',2f16.4,
     &        //,' Function Value :  ',f16.4,
     &        //,' Gradients :       ',2f16.4)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  function art1  --  function & gradient for vmetric  ##
c     ##                                                      ##
c     ##########################################################
c
c
      function art1 (xx,g)
      implicit none
      include 'sizes.i'
      real*8 art1,barkema
      real*8 eps,old,fore,back
      real*8 xx(maxopt),g(maxopt)
c
c
c     value of the Barkema G norm
c
      eps = 0.0001d0
      art1 = barkema (xx)
c
c     gradient of Barkema G norm in x direction
c
      old = xx(1)
      xx(1) = xx(1) + 0.5d0*eps
      fore = barkema (xx)
      xx(1) = xx(1) - eps
      back = barkema (xx)
      xx(1) = old
      g(1) = (fore - back) / eps
c
c     gradient of Barkema G norm in y direction
c
      old = xx(2)
      xx(2) = xx(2) + 0.5d0*eps
      fore = barkema (xx)
      xx(2) = xx(2) - eps
      back = barkema (xx)
      xx(2) = old
      g(2) = (fore - back) / eps
c
c     write out the Barkema G norm and its gradient
c
c     write (*,10)  art1,g(1),g(2)
c  10 format (' ART1  --  ',5x,3f16.4)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  function barkema  --  norm of the Barkema G vector  ##
c     ##                                                      ##
c     ##########################################################
c
c
      function barkema (xx)
      implicit none
      include 'sizes.i'
      real*8 barkema,xx(*)
      real*8 grad(2),gbark(2)
      real*8 force(2),dx(2)
      real*8 norm,project,alpha
      real*8 xref(2)
      common /art/ xref
c
c
c     compute the normalized displacement from local minimum
c
      dx(1) = xx(1) - xref(1)
      dx(2) = xx(2) - xref(2)
      norm = sqrt(dx(1)**2 + dx(2)**2)
      dx(1) = dx(1) / norm
      dx(2) = dx(2) / norm
c
c     get the gradient on the Muller-Brown surface
c
      call gmuller (xx(1),xx(2),grad(1),grad(2))
      force(1) = -grad(1)
      force(2) = -grad(2)
c
c     set alpha and find its product with the projected force
c
      alpha = 0.25d0
      project = (1.0d0 + alpha) * (force(1)*dx(1) + force(2)*dx(2))
c
c     calculate the Barkema force components
c
      gbark(1) = force(1) - project * dx(1)
      gbark(2) = force(2) - project * dx(2)
c
c     calculate object function as the norm of Barkema force
c
      barkema = sqrt(gbark(1)**2 + gbark(2)**2)
      return
      end
c
c
c     ########################
c     ##                    ##
c     ##  function fmuller  ##
c     ##                    ##
c     ########################
c
c
      function fmuller (x,y)
      implicit none
      real*8 x,y,fmuller
      real*8 t1,t2,t3,t4
c
c
c     compute the Muller-Brown saddle point test surface
c
      t1 = -200.0d0 * exp(-(x-1.0d0)**2-10.0d0*y**2)
      t2 = -100.0d0 * exp(-x**2-10.0d0*(y-0.5d0)**2)
      t3 = -170.0d0 * exp(-6.5d0*(x+0.5d0)**2
     &                    +11.0d0*(x+0.5d0)*(y-1.5d0)
     &                    -6.5d0*(y-1.5d0)**2)
      t4 = 15.0d0 * exp(0.7d0*(x+1.0d0)**2
     &                  +0.6d0*(x+1.0d0)*(y-1.0d0)
     &                  +0.7d0*(y-1.0d0)**2)
      fmuller = t1 + t2 + t3 + t4
c
c     write out the Muller-Brown function value
c
c     write (*,10)  fmuller
c  10 format (' FMULLER  --  ',5x,f16.4)
      return
      end
c
c
c     ##########################
c     ##                      ##
c     ##  subroutine gmuller  ##
c     ##                      ##
c     ##########################
c
c
      subroutine gmuller (x,y,gx,gy)
      implicit none
      real*8 x,y,gx,gy
      real*8 t1,t2,t3,t4
c
c
c     compute the derivatives with respect to x
c
      t1 = -200.d0 * (-2.0d0*x+2.0d0) * exp(-(x-1.0d0)**2-10.0d0*y**2)
      t2 = 200.0d0 * x * exp(-x**2-10.0d0*(y-0.5d0)**2)
      t3 = -170.0d0 * (-13.0d0*x-23.0d0+11.0d0*y)
     &                     * exp(-6.5d0*(x+0.5d0)**2
     &                           +11.0d0*(x+0.5d0)*(y-1.5d0)
     &                           -6.5d0*(y-1.5d0)**2)
      t4 = 15.0d0 * (1.4d0*x+0.8d0+0.6d0*y)
     &                     * exp(0.7d0*(x+1.0d0)**2
     &                           +0.6d0*(x+1.0d0)*(y-1.0d0)
     &                           +0.7d0*(y-1.0d0)**2)
      gx = t1 + t2 + t3 + t4
c
c     compute the derivatives with respect to y
c
      t1 = 4000.0d0 * y * exp(-(x-1.0d0)**2-10.0d0*y**2)
      t2 = -100.0d0 * (-20.0d0*y+10.0d0)
     &                     * exp(-x**2-10.0d0*(y-0.5d0)**2)
      t3 = -170.0d0 * (11.0d0*x+25.0d0-13.0d0*y)
     &                     * exp(-6.5d0*(x+0.5d0)**2
     &                           +11.0d0*(x+0.5d0)*(y-1.5d0)
     &                           -6.5d0*(y-1.5d0)**2)
      t4 = 15.0d0 * (0.6d0*x-0.8d0+1.4d0*y)
     &                     * exp(0.7d0*(x+1.0d0)**2
     &                           +0.6d0*(x+1.0d0)*(y-1.0d0)
     &                           +0.7d0*(y-1.0d0)**2)
      gy = t1 + t2 + t3 + t4
c
c     write out the Muller-Brown function value and gradient
c
c     write (*,10)  gx,gy
c  10 format (' GMULLER  --  ',5x,2f16.4)
      return
      end

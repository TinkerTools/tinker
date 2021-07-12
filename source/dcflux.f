c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2020 by Chengwen Liu & Jay W. Ponder  ##
c     ##                 All Rights Reserved                  ##
c     ##########################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine dcflux  --  charge flux gradient chain rule  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "dcflux" takes as input the electrostatic potential at each
c     atomic site and calculates gradient chain rule terms due to
c     charge flux coupling with bond stretching and angle bending
c
c     literature reference:
c
c     C. Liu, J.-P. Piquemal and P. Ren, "Implementation of Geometry-
c     Dependent Charge Flux into the Polarizable AMOEBA+ Potential",
c     Journal of Physical Chemistry Letters, 11, 419-426 (2020)
c
c
      subroutine dcflux (pot,dcfx,dcfy,dcfz)
      use sizes
      use angbnd
      use atoms
      use bndstr
      use bound
      use cflux
      implicit none
      integer i,ia,ib,ic
      real*8 xa,ya,za
      real*8 xb,yb,zb
      real*8 xc,yc,zc
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab,rab2,rab3
      real*8 rcb,rcb2,rcb3
      real*8 dpot,dpota,dpotc
      real*8 ddqdx,ddqdy,ddqdz
      real*8 fx,fy,fz
      real*8 fxa1,fya1,fza1
      real*8 fxb1,fyb1,fzb1
      real*8 fxc1,fyc1,fzc1
      real*8 fxa2,fya2,fza2
      real*8 fxb2,fyb2,fzb2
      real*8 fxc2,fyc2,fzc2
      real*8 pb,pb1,pb2
      real*8 pa1,pa2
      real*8 eps,dot
      real*8 term,fterm
      real*8 termxa,termxc
      real*8 termya,termyc
      real*8 termza,termzc
      real*8 pot(*)
      real*8 dcfx(*)
      real*8 dcfy(*)
      real*8 dcfz(*)
c
c
c     zero out the charge flux correction forces
c
      do i = 1, n
         dcfx(i) = 0.0d0
         dcfy(i) = 0.0d0
         dcfz(i) = 0.0d0
      end do
c
c     set tolerance for minimum distance and angle values
c
      eps = 0.0001d0
c
c     calculate the charge flux forces due to bond stretches
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         pb = bflx(i)
         xa = x(ia)
         ya = y(ia)
         za = z(ia)
         xb = x(ib)
         yb = y(ib)
         zb = z(ib)
         xab = xa - xb
         yab = ya - yb
         zab = za - zb
         if (use_polymer)  call image (xab,yab,zab)
         rab2 = max(xab*xab+yab*yab+zab*zab,eps)
         pb = pb / sqrt(rab2)
         dpot = pot(ib) - pot(ia)
         ddqdx = pb * (xa-xb)
         ddqdy = pb * (ya-yb)
         ddqdz = pb * (za-zb)
         fx = dpot * ddqdx
         fy = dpot * ddqdy
         fz = dpot * ddqdz
         dcfx(ia) = dcfx(ia) + fx
         dcfy(ia) = dcfy(ia) + fy
         dcfz(ia) = dcfz(ia) + fz
         dcfx(ib) = dcfx(ib) - fx
         dcfy(ib) = dcfy(ib) - fy
         dcfz(ib) = dcfz(ib) - fz
      end do
c
c     calculate the charge flux forces due to angle bends
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         pa1 = aflx(1,i)
         pa2 = aflx(2,i)
         pb1 = abflx(1,i)
         pb2 = abflx(2,i)
         xa = x(ia)
         ya = y(ia)
         za = z(ia)
         xb = x(ib)
         yb = y(ib)
         zb = z(ib)
         xc = x(ic)
         yc = y(ic)
         zc = z(ic)
         xab = xa - xb
         yab = ya - yb
         zab = za - zb
         xcb = xc - xb
         ycb = yc - yb
         zcb = zc - zb
         if (use_polymer) then
            call image (xab,yab,zab)
            call image (xcb,ycb,zcb)
         end if
         rab2 = max(xab*xab+yab*yab+zab*zab,eps)
         rcb2 = max(xcb*xcb+ycb*ycb+zcb*zcb,eps)
         rab  = sqrt(rab2)
         rab3 = rab2 * rab
         rcb  = sqrt(rcb2)
         rcb3 = rcb2 * rcb
c
c     get terms corresponding to asymmetric bond stretches
c
         dpota = pot(ia) - pot(ib)
         dpotc = pot(ic) - pot(ib)
         pb1 = dpota * pb1
         pb2 = dpotc * pb2
         fxa1 = pb2 * xab/rab
         fya1 = pb2 * yab/rab
         fza1 = pb2 * zab/rab
         fxc1 = pb1 * xcb/rcb
         fyc1 = pb1 * ycb/rcb
         fzc1 = pb1 * zcb/rcb
         fxb1 = -fxa1 - fxc1
         fyb1 = -fya1 - fyc1
         fzb1 = -fza1 - fzc1
c
c     get terms corresponding to bond angle bending
c
         dot = xab*xcb + yab*ycb + zab*zcb
         term = -rab*rcb / max(sqrt(rab2*rcb2-dot*dot),eps)
         fterm = term * (dpota*pa1+dpotc*pa2)
         termxa = xcb/(rab*rcb) - xab*dot/(rab3*rcb)
         termya = ycb/(rab*rcb) - yab*dot/(rab3*rcb)
         termza = zcb/(rab*rcb) - zab*dot/(rab3*rcb)
         termxc = xab/(rab*rcb) - xcb*dot/(rab*rcb3)
         termyc = yab/(rab*rcb) - ycb*dot/(rab*rcb3)
         termzc = zab/(rab*rcb) - zcb*dot/(rab*rcb3)
         fxa2 = fterm * termxa
         fya2 = fterm * termya
         fza2 = fterm * termza
         fxc2 = fterm * termxc
         fyc2 = fterm * termyc
         fzc2 = fterm * termzc
         fxb2 = -fxa2 - fxc2
         fyb2 = -fya2 - fyc2
         fzb2 = -fza2 - fzc2
         dcfx(ia) = dcfx(ia) + fxa1 + fxa2
         dcfy(ia) = dcfy(ia) + fya1 + fya2
         dcfz(ia) = dcfz(ia) + fza1 + fza2
         dcfx(ib) = dcfx(ib) + fxb1 + fxb2
         dcfy(ib) = dcfy(ib) + fyb1 + fyb2
         dcfz(ib) = dcfz(ib) + fzb1 + fzb2
         dcfx(ic) = dcfx(ic) + fxc1 + fxc2
         dcfy(ic) = dcfy(ic) + fyc1 + fyc2
         dcfz(ic) = dcfz(ic) + fzc1 + fzc2
      end do
      return
      end

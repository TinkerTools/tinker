c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2020 by Chengwen Liu & Jay W. Ponder  ##
c     ##                 All Rights Reserved                  ##
c     ##########################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine alterchg  --  modification of partial charges  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "alterchg" calculates the change in atomic partial charge or
c     monopole values due to bond and angle charge flux coupling
c
c     literature reference:
c
c     C. Liu, J.-P. Piquemal and P. Ren, "Implementation of Geometry-
c     Dependent Charge Flux into the Polarizable AMOEBA+ Potential",
c     Journal of Physical Chemistry Letters, 11, 419-426 (2020)
c
c
      subroutine alterchg
      use atoms
      use charge
      use chgpen
      use inform
      use iounit
      use mplpot
      use mpole
      implicit none
      integer i,k
      real*8, allocatable :: pdelta(:)
      logical header
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pdelta(n))
c
c     zero out the change in charge value at each site
c
      do i = 1, n
         pdelta(i) = 0.0d0
      end do
c
c     find charge modifications due to charge flux
c
      call bndchg (pdelta)
      call angchg (pdelta)
c
c     alter atomic partial charge values for charge flux
c
      header = .true.
      do i = 1, nion
         k = iion(i)
         pchg(i) = pchg0(i) + pdelta(k)
         if (debug.ne.0 .and. pdelta(k).ne.0.0d0) then
            if (header) then
               header = .false.
               write (iout,10)
   10          format (/,' Charge Flux Modification of Partial',
     &                    ' Charges :',
     &                 //,4x,'Atom',14x,'Base Value',7x,'Actual',/)
            end if
            write (iout,20)  k,pchg0(i),pchg(i)
   20       format (i8,9x,2f14.5)
         end if
      end do
c
c     alter monopoles and charge penetration for charge flux
c
      header = .true.
      do i = 1, npole
         k = ipole(i)
         pole(1,i) = mono0(i) + pdelta(k)
         if (use_chgpen)  pval(i) = pval0(i) + pdelta(k)
         if (debug.ne.0 .and. pdelta(k).ne.0.0d0) then
            if (header) then
               header = .false.
               write (iout,30)
   30          format (/,' Charge Flux Modification of Atomic',
     &                    ' Monopoles :',
     &                 //,4x,'Atom',14x,'Base Value',7x,'Actual',/)
            end if
            write (iout,40)  k,mono0(i),pole(1,i)
   40       format (i8,9x,2f14.5)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pdelta)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine bndchg  --  charge flux bond stretch coupling  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "bndchg" computes modifications to atomic partial charges or
c     monopoles due to bond stretch using a charge flux formulation
c
c
      subroutine bndchg (pdelta)
      use sizes
      use atoms
      use bndstr
      use bound
      use cflux
      implicit none
      integer i,ia,ib
      real*8 xab,yab,zab
      real*8 rab,rab0
      real*8 pb,dq
      real*8 pdelta(*)
c
c
c     loop over all the bond distances in the system
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         pb = bflx(i)
c
c     compute the bond length value for the current bond
c
         xab = x(ia) - x(ib)
         yab = y(ia) - y(ib)
         zab = z(ia) - z(ib)
         if (use_polymer)  call image (xab,yab,zab)
         rab = sqrt(xab*xab + yab*yab + zab*zab)
c
c     find the charge flux increment for the current bond
c
         rab0 = bl(i)
         dq = pb * (rab-rab0)
         pdelta(ia) = pdelta(ia) - dq
         pdelta(ib) = pdelta(ib) + dq
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine angchg  --  charge flux angle bend coupling  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "angchg" computes modifications to atomic partial charges or
c     monopoles due to angle bending using a charge flux formulation
c
c
      subroutine angchg (pdelta)
      use sizes
      use angbnd
      use atmlst
      use atoms
      use bndstr
      use bound
      use cflux
      use math
      implicit none
      integer i,ia,ib,ic
      real*8 angle,eps
      real*8 rab,rcb
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 dot,cosine
      real*8 pa1,pa2
      real*8 pb1,pb2
      real*8 theta0
      real*8 rab0,rcb0
      real*8 dq1,dq2
      real*8 pdelta(*)
c
c
c     loop over all the bond angles in the system
c
      eps = 0.0001d0
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         pa1 = aflx(1,i)
         pa2 = aflx(2,i)
         pb1 = abflx(1,i)
         pb2 = abflx(2,i)
c
c     calculate the angle values and included bond lengths
c
         xia = x(ia)
         yia = y(ia)
         zia = z(ia)
         xib = x(ib)
         yib = y(ib)
         zib = z(ib)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xab = xia - xib
         yab = yia - yib
         zab = zia - zib
         xcb = xic - xib
         ycb = yic - yib
         zcb = zic - zib
         if (use_polymer) then
            call image (xab,yab,zab)
            call image (xcb,ycb,zcb)
         end if
         rab = sqrt(max(xab*xab+yab*yab+zab*zab,eps))
         rcb = sqrt(max(xcb*xcb+ycb*ycb+zcb*zcb,eps))
         dot = xab*xcb + yab*ycb + zab*zcb
         cosine = dot / (rab*rcb)
         cosine = min(1.0d0,max(-1.0d0,cosine))
         angle = radian * acos(cosine)
c
c     find the charge flux increment for the current angle
c
         theta0 = anat(i)
         rab0 = bl(balist(1,i))
         rcb0 = bl(balist(2,i))
         dq1 = pb1*(rcb-rcb0) + pa1*(angle-theta0)/radian
         dq2 = pb2*(rab-rab0) + pa2*(angle-theta0)/radian
         pdelta(ia) = pdelta(ia) + dq1
         pdelta(ib) = pdelta(ib) - dq1 - dq2
         pdelta(ic) = pdelta(ic) + dq2
      end do
      return
      end

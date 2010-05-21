c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Thomas Darden & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  routines below implement various coordinate and B-spline  ##
c     ##  manipulations for particle mesh Ewald summation applied   ##
c     ##  to polarizable atomic multipoles; modified from original  ##
c     ##  code due to Thomas Darden, NIEHS, Research Triangle, NC   ##
c     ##                                                            ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine bspline_fill  --  set multipole B-spline coeffs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "bspline_fill" finds B-spline coefficients and derivatives
c     for multipole sites along the fractional coordinate axes
c
c
      subroutine bspline_fill
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'mpole.i'
      include 'pme.i'
      integer i,ii,ifr
      real*8 w,fr
c
c
c     get the B-spline coefficients for each multipole site
c
      do i = 1, npole
         ii = ipole(i)
         w = x(ii)*recip(1,1) + y(ii)*recip(2,1) + z(ii)*recip(3,1)
         fr = dble(nfft1) * (w-anint(w)+0.5d0)
         ifr = int(fr)
         w = fr - dble(ifr)
         igrid(1,i) = ifr - bsorder
         call bspline_gen (w,thetai1(1,1,i))
         w = x(ii)*recip(1,2) + y(ii)*recip(2,2) + z(ii)*recip(3,2)
         fr = dble(nfft2) * (w-anint(w)+0.5d0)
         ifr = int(fr)
         w = fr - dble(ifr)
         igrid(2,i) = ifr - bsorder
         call bspline_gen (w,thetai2(1,1,i))
         w = x(ii)*recip(1,3) + y(ii)*recip(2,3) + z(ii)*recip(3,3)
         fr = dble(nfft3) * (w-anint(w)+0.5d0)
         ifr = int(fr)
         w = fr - dble(ifr)
         igrid(3,i) = ifr - bsorder
         call bspline_gen (w,thetai3(1,1,i))
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine bspline_gen  --  single-site B-spline coeffs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "bspline_gen" gets B-spline coefficients and derivatives for
c     a single multipole site along a particular direction
c
c
      subroutine bspline_gen (w,thetai)
      implicit none
      include 'sizes.i'
      include 'pme.i'
      integer i,j,k
      real*8 w,denom
      real*8 thetai(4,maxorder)
      real*8 array(maxorder,maxorder)
c
c
c     initialization to get to 2nd order recursion
c
      array(2,2) = w
      array(2,1) = 1.0d0 - w
c
c     perform one pass to get to 3rd order recursion
c
      array(3,3) = 0.5d0 * w * array(2,2)
      array(3,2) = 0.5d0 * ((1.0d0+w)*array(2,1)+(2.0d0-w)*array(2,2))
      array(3,1) = 0.5d0 * (1.0d0-w) * array(2,1)
c
c     compute standard B-spline recursion to desired order
c
      do i = 4, bsorder
         k = i - 1
         denom = 1.0d0 / dble(k)
         array(i,i) = denom * w * array(k,k)
         do j = 1, i-2
            array(i,i-j) = denom * ((w+dble(j))*array(k,i-j-1)
     &                             +(dble(i-j)-w)*array(k,i-j))
         end do
         array(i,1) = denom * (1.0d0-w) * array(k,1)
      end do
c
c     get coefficients for the B-spline first derivative
c
      k = bsorder - 1
      array(k,bsorder) = array(k,bsorder-1)
      do i = bsorder-1, 2, -1
         array(k,i) = array(k,i-1) - array(k,i)
      end do
      array(k,1) = -array(k,1)
c
c     get coefficients for the B-spline second derivative
c
      k = bsorder - 2
      array(k,bsorder-1) = array(k,bsorder-2)
      do i = bsorder-2, 2, -1
         array(k,i) = array(k,i-1) - array(k,i)
      end do
      array(k,1) = -array(k,1)
      array(k,bsorder) = array(k,bsorder-1)
      do i = bsorder-1, 2, -1
         array(k,i) = array(k,i-1) - array(k,i)
      end do
      array(k,1) = -array(k,1)
c
c     get coefficients for the B-spline third derivative
c
      k = bsorder - 3
      array(k,bsorder-2) = array(k,bsorder-3)
      do i = bsorder-3, 2, -1
         array(k,i) = array(k,i-1) - array(k,i)
      end do
      array(k,1) = -array(k,1)
      array(k,bsorder-1) = array(k,bsorder-2)
      do i = bsorder-2, 2, -1
         array(k,i) = array(k,i-1) - array(k,i)
      end do
      array(k,1) = -array(k,1)
      array(k,bsorder) = array(k,bsorder-1)
      do i = bsorder-1, 2, -1
         array(k,i) = array(k,i-1) - array(k,i)
      end do
      array(k,1) = -array(k,1)
c
c     copy coefficients from temporary to permanent storage
c
      do i = 1, bsorder
         do j = 1, 4
            thetai(j,i) = array(bsorder-j+1,i)
         end do
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine grid_mpole  --  put multipoles on PME grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "grid_mpole" places the fractional atomic multipoles onto
c     the particle mesh Ewald grid
c
c
      subroutine grid_mpole (fmp)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'pme.i'
      integer i,j,k,m
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,u0,t0
      real*8 v1,u1,t1
      real*8 v2,u2,t2
      real*8 term0,term1,term2
      real*8 fmp(10,maxatm)
c
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     put the permanent multipole moments onto the grid
c
      do m = 1, npole
         igrd0 = igrid(1,m)
         jgrd0 = igrid(2,m)
         kgrd0 = igrid(3,m)
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,m)
            v1 = thetai3(2,it3,m)
            v2 = thetai3(3,it3,m)
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,m)
               u1 = thetai2(2,it2,m)
               u2 = thetai2(3,it2,m)
               term0 = fmp(1,m)*u0*v0 + fmp(3,m)*u1*v0
     &                    + fmp(4,m)*u0*v1 + fmp(6,m)*u2*v0
     &                    + fmp(7,m)*u0*v2 + fmp(10,m)*u1*v1
               term1 = fmp(2,m)*u0*v0 + fmp(8,m)*u1*v0
     &                    + fmp(9,m)*u0*v1
               term2 = fmp(5,m) * u0 * v0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  t0 = thetai1(1,it1,m)
                  t1 = thetai1(2,it1,m)
                  t2 = thetai1(3,it1,m)
                  qgrid(1,i,j,k) = qgrid(1,i,j,k) + term0*t0
     &                                + term1*t1 + term2*t2
               end do
            end do
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine grid_uind  --  put induced dipoles on PME grid  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "grid_uind" places the fractional induced dipoles onto the
c     particle mesh Ewald grid
c
c
      subroutine grid_uind (fuind,fuinp)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'pme.i'
      integer i,j,k,m
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,u0,t0
      real*8 v1,u1,t1
      real*8 term01,term11
      real*8 term02,term12
      real*8 fuind(3,maxatm)
      real*8 fuinp(3,maxatm)
c
c
c     zero out the particle mesh Ewald charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               qgrid(1,i,j,k) = 0.0d0
               qgrid(2,i,j,k) = 0.0d0
            end do
         end do
      end do
c
c     put the induced dipole moments onto the grid
c
      do m = 1, npole
         igrd0 = igrid(1,m)
         jgrd0 = igrid(2,m)
         kgrd0 = igrid(3,m)
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,m)
            v1 = thetai3(2,it3,m)
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,m)
               u1 = thetai2(2,it2,m)
               term01 = fuind(2,m)*u1*v0 + fuind(3,m)*u0*v1
               term11 = fuind(1,m)*u0*v0
               term02 = fuinp(2,m)*u1*v0 + fuinp(3,m)*u0*v1
               term12 = fuinp(1,m)*u0*v0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  t0 = thetai1(1,it1,m)
                  t1 = thetai1(2,it1,m)
                  qgrid(1,i,j,k) = qgrid(1,i,j,k) + term01*t0
     &                                + term11*t1
                  qgrid(2,i,j,k) = qgrid(2,i,j,k) + term02*t0
     &                                + term12*t1
               end do
            end do
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_mpole  --  multipole potential from grid  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_mpole" extracts the permanent multipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_mpole (fphi)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'pme.i'
      integer i,j,k,m
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3,tq
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu21,tu12,tu30,tu03
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111
      real*8 fphi(20,maxatm)
c
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,igrid,bsorder,
!$OMP& nfft3,thetai3,nfft2,thetai2,nfft1,thetai1,qgrid,fphi)
!$OMP DO
c
c     extract the permanent multipole field at each site
c
      do m = 1, npole
         igrd0 = igrid(1,m)
         jgrd0 = igrid(2,m)
         kgrd0 = igrid(3,m)
         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
         tuv200 = 0.0d0
         tuv020 = 0.0d0
         tuv002 = 0.0d0
         tuv110 = 0.0d0
         tuv101 = 0.0d0
         tuv011 = 0.0d0
         tuv300 = 0.0d0
         tuv030 = 0.0d0
         tuv003 = 0.0d0
         tuv210 = 0.0d0
         tuv201 = 0.0d0
         tuv120 = 0.0d0
         tuv021 = 0.0d0
         tuv102 = 0.0d0
         tuv012 = 0.0d0
         tuv111 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,m)
            v1 = thetai3(2,it3,m)
            v2 = thetai3(3,it3,m)
            v3 = thetai3(4,it3,m)
            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
            tu20 = 0.0d0
            tu11 = 0.0d0
            tu02 = 0.0d0
            tu30 = 0.0d0
            tu21 = 0.0d0
            tu12 = 0.0d0
            tu03 = 0.0d0
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,m)
               u1 = thetai2(2,it2,m)
               u2 = thetai2(3,it2,m)
               u3 = thetai2(4,it2,m)
               t0 = 0.0d0
               t1 = 0.0d0
               t2 = 0.0d0
               t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq = qgrid(1,i,j,k)
                  t0 = t0 + tq*thetai1(1,it1,m)
                  t1 = t1 + tq*thetai1(2,it1,m)
                  t2 = t2 + tq*thetai1(3,it1,m)
                  t3 = t3 + tq*thetai1(4,it1,m)
               end do
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0 
               tu21 = tu21 + t2*u1 
               tu12 = tu12 + t1*u2 
               tu03 = tu03 + t0*u3
            end do
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1
         end do
         fphi(1,m) = tuv000
         fphi(2,m) = tuv100
         fphi(3,m) = tuv010
         fphi(4,m) = tuv001
         fphi(5,m) = tuv200
         fphi(6,m) = tuv020
         fphi(7,m) = tuv002
         fphi(8,m) = tuv110
         fphi(9,m) = tuv101
         fphi(10,m) = tuv011
         fphi(11,m) = tuv300
         fphi(12,m) = tuv030
         fphi(13,m) = tuv003
         fphi(14,m) = tuv210
         fphi(15,m) = tuv201
         fphi(16,m) = tuv120
         fphi(17,m) = tuv021
         fphi(18,m) = tuv102
         fphi(19,m) = tuv012
         fphi(20,m) = tuv111
      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine fphi_uind  --  induced potential from grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "fphi_uind" extracts the induced dipole potential from
c     the particle mesh Ewald grid
c
c
      subroutine fphi_uind (fdip_phi1,fdip_phi2,fdip_sum_phi)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'pme.i'
      integer i,j,k,m
      integer i0,j0,k0
      integer it1,it2,it3
      integer igrd0,jgrd0,kgrd0
      real*8 v0,v1,v2,v3
      real*8 u0,u1,u2,u3
      real*8 t0,t1,t2,t3
      real*8 t0_1,t0_2,t1_1,t1_2
      real*8 t2_1,t2_2,tq_1,tq_2
      real*8 tu00,tu10,tu01,tu20,tu11
      real*8 tu02,tu30,tu21,tu12,tu03
      real*8 tu00_1,tu01_1,tu10_1
      real*8 tu00_2,tu01_2,tu10_2
      real*8 tu20_1,tu11_1,tu02_1
      real*8 tu20_2,tu11_2,tu02_2
      real*8 tuv100_1,tuv010_1,tuv001_1
      real*8 tuv100_2,tuv010_2,tuv001_2
      real*8 tuv200_1,tuv020_1,tuv002_1
      real*8 tuv110_1,tuv101_1,tuv011_1
      real*8 tuv200_2,tuv020_2,tuv002_2
      real*8 tuv110_2,tuv101_2,tuv011_2
      real*8 tuv000,tuv100,tuv010,tuv001
      real*8 tuv200,tuv020,tuv002,tuv110
      real*8 tuv101,tuv011,tuv300,tuv030
      real*8 tuv003,tuv210,tuv201,tuv120
      real*8 tuv021,tuv102,tuv012,tuv111
      real*8 fdip_phi1(10,maxatm)
      real*8 fdip_phi2(10,maxatm)
      real*8 fdip_sum_phi(20,maxatm)
c
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,igrid,
!$OMP& bsorder,nfft3,thetai3,nfft2,thetai2,nfft1,
!$OMP& thetai1,qgrid,fdip_phi1,fdip_phi2,fdip_sum_phi)
!$OMP DO
c
c     extract the induced dipole field at each site
c
      do m = 1, npole
         igrd0 = igrid(1,m)
         jgrd0 = igrid(2,m)
         kgrd0 = igrid(3,m)
         tuv100_1 = 0.0d0
         tuv010_1 = 0.0d0
         tuv001_1 = 0.0d0
         tuv200_1 = 0.0d0
         tuv020_1 = 0.0d0
         tuv002_1 = 0.0d0
         tuv110_1 = 0.0d0
         tuv101_1 = 0.0d0
         tuv011_1 = 0.0d0
         tuv100_2 = 0.0d0
         tuv010_2 = 0.0d0
         tuv001_2 = 0.0d0
         tuv200_2 = 0.0d0
         tuv020_2 = 0.0d0
         tuv002_2 = 0.0d0
         tuv110_2 = 0.0d0
         tuv101_2 = 0.0d0
         tuv011_2 = 0.0d0
         tuv000 = 0.0d0
         tuv001 = 0.0d0
         tuv010 = 0.0d0
         tuv100 = 0.0d0
         tuv200 = 0.0d0
         tuv020 = 0.0d0
         tuv002 = 0.0d0
         tuv110 = 0.0d0
         tuv101 = 0.0d0
         tuv011 = 0.0d0
         tuv300 = 0.0d0
         tuv030 = 0.0d0
         tuv003 = 0.0d0
         tuv210 = 0.0d0
         tuv201 = 0.0d0
         tuv120 = 0.0d0
         tuv021 = 0.0d0
         tuv102 = 0.0d0
         tuv012 = 0.0d0
         tuv111 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,m)
            v1 = thetai3(2,it3,m)
            v2 = thetai3(3,it3,m)
            v3 = thetai3(4,it3,m)
            tu00_1 = 0.0d0
            tu01_1 = 0.0d0
            tu10_1 = 0.0d0
            tu20_1 = 0.0d0
            tu11_1 = 0.0d0
            tu02_1 = 0.0d0
            tu00_2 = 0.0d0
            tu01_2 = 0.0d0
            tu10_2 = 0.0d0
            tu20_2 = 0.0d0
            tu11_2 = 0.0d0
            tu02_2 = 0.0d0
            tu00 = 0.0d0
            tu10 = 0.0d0
            tu01 = 0.0d0
            tu20 = 0.0d0
            tu11 = 0.0d0
            tu02 = 0.0d0
            tu30 = 0.0d0
            tu21 = 0.0d0
            tu12 = 0.0d0
            tu03 = 0.0d0
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,m)
               u1 = thetai2(2,it2,m)
               u2 = thetai2(3,it2,m)
               u3 = thetai2(4,it2,m)
               t0_1 = 0.0d0
               t1_1 = 0.0d0
               t2_1 = 0.0d0
               t0_2 = 0.0d0
               t1_2 = 0.0d0
               t2_2 = 0.0d0
               t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2
                  tq_1 = qgrid(1,i,j,k)
                  tq_2 = qgrid(2,i,j,k)
                  t0_1 = t0_1 + tq_1*thetai1(1,it1,m)
                  t1_1 = t1_1 + tq_1*thetai1(2,it1,m)
                  t2_1 = t2_1 + tq_1*thetai1(3,it1,m)
                  t0_2 = t0_2 + tq_2*thetai1(1,it1,m)
                  t1_2 = t1_2 + tq_2*thetai1(2,it1,m)
                  t2_2 = t2_2 + tq_2*thetai1(3,it1,m)
                  t3 = t3 + (tq_1+tq_2)*thetai1(4,it1,m)
               end do
               tu00_1 = tu00_1 + t0_1*u0
               tu10_1 = tu10_1 + t1_1*u0
               tu01_1 = tu01_1 + t0_1*u1
               tu20_1 = tu20_1 + t2_1*u0
               tu11_1 = tu11_1 + t1_1*u1
               tu02_1 = tu02_1 + t0_1*u2
               tu00_2 = tu00_2 + t0_2*u0
               tu10_2 = tu10_2 + t1_2*u0
               tu01_2 = tu01_2 + t0_2*u1
               tu20_2 = tu20_2 + t2_2*u0
               tu11_2 = tu11_2 + t1_2*u1
               tu02_2 = tu02_2 + t0_2*u2
               t0 = t0_1 + t0_2
               t1 = t1_1 + t1_2
               t2 = t2_1 + t2_2
               tu00 = tu00 + t0*u0
               tu10 = tu10 + t1*u0
               tu01 = tu01 + t0*u1
               tu20 = tu20 + t2*u0
               tu11 = tu11 + t1*u1
               tu02 = tu02 + t0*u2
               tu30 = tu30 + t3*u0 
               tu21 = tu21 + t2*u1 
               tu12 = tu12 + t1*u2 
               tu03 = tu03 + t0*u3
            end do
            tuv100_1 = tuv100_1 + tu10_1*v0
            tuv010_1 = tuv010_1 + tu01_1*v0
            tuv001_1 = tuv001_1 + tu00_1*v1
            tuv200_1 = tuv200_1 + tu20_1*v0
            tuv020_1 = tuv020_1 + tu02_1*v0
            tuv002_1 = tuv002_1 + tu00_1*v2
            tuv110_1 = tuv110_1 + tu11_1*v0
            tuv101_1 = tuv101_1 + tu10_1*v1
            tuv011_1 = tuv011_1 + tu01_1*v1
            tuv100_2 = tuv100_2 + tu10_2*v0
            tuv010_2 = tuv010_2 + tu01_2*v0
            tuv001_2 = tuv001_2 + tu00_2*v1
            tuv200_2 = tuv200_2 + tu20_2*v0
            tuv020_2 = tuv020_2 + tu02_2*v0
            tuv002_2 = tuv002_2 + tu00_2*v2
            tuv110_2 = tuv110_2 + tu11_2*v0
            tuv101_2 = tuv101_2 + tu10_2*v1
            tuv011_2 = tuv011_2 + tu01_2*v1
            tuv000 = tuv000 + tu00*v0
            tuv100 = tuv100 + tu10*v0
            tuv010 = tuv010 + tu01*v0
            tuv001 = tuv001 + tu00*v1
            tuv200 = tuv200 + tu20*v0
            tuv020 = tuv020 + tu02*v0
            tuv002 = tuv002 + tu00*v2
            tuv110 = tuv110 + tu11*v0
            tuv101 = tuv101 + tu10*v1
            tuv011 = tuv011 + tu01*v1
            tuv300 = tuv300 + tu30*v0
            tuv030 = tuv030 + tu03*v0
            tuv003 = tuv003 + tu00*v3
            tuv210 = tuv210 + tu21*v0
            tuv201 = tuv201 + tu20*v1
            tuv120 = tuv120 + tu12*v0
            tuv021 = tuv021 + tu02*v1
            tuv102 = tuv102 + tu10*v2
            tuv012 = tuv012 + tu01*v2
            tuv111 = tuv111 + tu11*v1
         end do
         fdip_phi1(2,m) = tuv100_1
         fdip_phi1(3,m) = tuv010_1
         fdip_phi1(4,m) = tuv001_1
         fdip_phi1(5,m) = tuv200_1
         fdip_phi1(6,m) = tuv020_1
         fdip_phi1(7,m) = tuv002_1
         fdip_phi1(8,m) = tuv110_1
         fdip_phi1(9,m) = tuv101_1
         fdip_phi1(10,m) = tuv011_1
         fdip_phi2(2,m) = tuv100_2
         fdip_phi2(3,m) = tuv010_2
         fdip_phi2(4,m) = tuv001_2
         fdip_phi2(5,m) = tuv200_2
         fdip_phi2(6,m) = tuv020_2
         fdip_phi2(7,m) = tuv002_2
         fdip_phi2(8,m) = tuv110_2
         fdip_phi2(9,m) = tuv101_2
         fdip_phi2(10,m) = tuv011_2
         fdip_sum_phi(1,m) = tuv000
         fdip_sum_phi(2,m) = tuv100
         fdip_sum_phi(3,m) = tuv010
         fdip_sum_phi(4,m) = tuv001
         fdip_sum_phi(5,m) = tuv200
         fdip_sum_phi(6,m) = tuv020
         fdip_sum_phi(7,m) = tuv002
         fdip_sum_phi(8,m) = tuv110
         fdip_sum_phi(9,m) = tuv101
         fdip_sum_phi(10,m) = tuv011
         fdip_sum_phi(11,m) = tuv300
         fdip_sum_phi(12,m) = tuv030
         fdip_sum_phi(13,m) = tuv003
         fdip_sum_phi(14,m) = tuv210
         fdip_sum_phi(15,m) = tuv201
         fdip_sum_phi(16,m) = tuv120
         fdip_sum_phi(17,m) = tuv021
         fdip_sum_phi(18,m) = tuv102
         fdip_sum_phi(19,m) = tuv012
         fdip_sum_phi(20,m) = tuv111
      end do
c
c     end OpenMP directive for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine cmp_to_fmp  --  transformation of multipoles  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "cmp_to_fmp" transforms the atomic multipoles from Cartesian
c     to fractional coordinates
c
c
      subroutine cmp_to_fmp (cmp,fmp)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      integer i,j,k
      real*8 ctf(10,10)
      real*8 cmp(10,maxatm)
      real*8 fmp(10,maxatm)
c
c
c     find the matrix to convert Cartesian to fractional
c
      call cart_to_frac (ctf)
c
c     apply the transformation to get the fractional multipoles
c
      do i = 1, npole
         fmp(1,i) = ctf(1,1) * cmp(1,i)
         do j = 2, 4
            fmp(j,i) = 0.0d0
            do k = 2, 4
               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
            end do
         end do
         do j = 5, 10
            fmp(j,i) = 0.0d0
            do k = 5, 10
               fmp(j,i) = fmp(j,i) + ctf(j,k)*cmp(k,i)
            end do
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine cart_to_frac  --  Cartesian to fractional  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "cart_to_frac" computes a transformation matrix to convert
c     a multipole object in Cartesian coordinates to fractional
c
c     note the multipole components are stored in the condensed
c     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
c
c
      subroutine cart_to_frac (ctf)
      implicit none
      include 'sizes.i'
      include 'boxes.i'
      include 'pme.i'
      integer i,j,k,m
      integer i1,i2
      integer qi1(6)
      integer qi2(6)
      real*8 a(3,3)
      real*8 ctf(10,10)
      data qi1  / 1, 2, 3, 1, 1, 2 /
      data qi2  / 1, 2, 3, 2, 3, 3 /
c
c
c     set the reciprocal vector transformation matrix
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
c
c     get the Cartesian to fractional conversion matrix
c
      do i = 1, 10
         do j = 1, 10
            ctf(j,i) = 0.0d0
         end do
      end do
      ctf(1,1) = 1.0d0
      do i = 2, 4
         do j = 2, 4
            ctf(i,j) = a(i-1,j-1)
         end do
      end do
      do i1 = 1, 3
         k = qi1(i1)
         do i2 = 1, 6
            i = qi1(i2)
            j = qi2(i2)
            ctf(i1+4,i2+4) = a(k,i) * a(k,j)
         end do
      end do
      do i1 = 4, 6
         k = qi1(i1)
         m = qi2(i1)
         do i2 = 1, 6
            i = qi1(i2)
            j = qi2(i2)
            ctf(i1+4,i2+4) = a(k,i)*a(m,j) + a(k,j)*a(m,i)
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fphi_to_cphi  --  transformation of potential  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fphi_to_cphi" transforms the reciprocal space potential from
c     fractional to Cartesian coordinates
c
c
      subroutine fphi_to_cphi (fphi,cphi)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      integer i,j,k
      real*8 ftc(10,10)
      real*8 cphi(10,maxatm)
      real*8 fphi(20,maxatm)
c
c
c     find the matrix to convert fractional to Cartesian
c
      call frac_to_cart (ftc)
c
c     apply the transformation to get the Cartesian potential
c
      do i = 1, npole
         cphi(1,i) = ftc(1,1) * fphi(1,i)
         do j = 2, 4
            cphi(j,i) = 0.0d0
            do k = 2, 4
               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
            end do
         end do
         do j = 5, 10
            cphi(j,i) = 0.0d0
            do k = 5, 10
               cphi(j,i) = cphi(j,i) + ftc(j,k)*fphi(k,i)
            end do
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine frac_to_cart  --  fractional to Cartesian  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "frac_to_cart" computes a transformation matrix to convert
c     a multipole object in fraction coordinates to Cartesian
c
c     note the multipole components are stored in the condensed
c     order (m,dx,dy,dz,qxx,qyy,qzz,qxy,qxz,qyz)
c
c
      subroutine frac_to_cart (ftc)
      implicit none
      include 'sizes.i'
      include 'boxes.i'
      include 'pme.i'
      integer i,j,k,m
      integer i1,i2
      integer qi1(6)
      integer qi2(6)
      real*8 a(3,3)
      real*8 ftc(10,10)
      data qi1  / 1, 2, 3, 1, 1, 2 /
      data qi2  / 1, 2, 3, 2, 3, 3 /
c
c
c     set the reciprocal vector transformation matrix
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
      end do
c
c     get the fractional to Cartesian conversion matrix
c
      do i = 1, 10
         do j = 1, 10
            ftc(j,i) = 0.0d0
         end do
      end do
      ftc(1,1) = 1.0d0
      do i = 2, 4
         do j = 2, 4
            ftc(i,j) = a(i-1,j-1)
         end do
      end do
      do i1 = 1, 3
         k = qi1(i1)
         do i2 = 1, 3
            i = qi1(i2)
            ftc(i1+4,i2+4) = a(k,i) * a(k,i)
         end do
         do i2 = 4, 6
            i = qi1(i2)
            j = qi2(i2)
            ftc(i1+4,i2+4) = 2.0d0 * a(k,i) * a(k,j)
         end do
      end do
      do i1 = 4, 6
         k = qi1(i1)
         m = qi2(i1)
         do i2 = 1, 3
            i = qi1(i2)
            ftc(i1+4,i2+4) = a(k,i) * a(m,i)
         end do
         do i2 = 4, 6
            i = qi1(i2)
            j = qi2(i2)
            ftc(i1+4,i2+4) = a(k,i)*a(m,j) + a(m,i)*a(k,j)
         end do
      end do
      return
      end

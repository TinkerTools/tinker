c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge2  --  atomwise charge-charge Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge2" calculates second derivatives of the
c     charge-charge interaction energy for a single atom
c
c
      subroutine echarge2 (i)
      use limits
      use warp
      implicit none
      integer i
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_smooth) then
         call echarge2d (i)
      else if (use_ewald) then
         call num_echarge2c (i)
c         call echarge2c (i)
      else if (use_clist) then
         call echarge2b (i)
      else
         call echarge2a (i)
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine echarge2a  --  charge Hessian via pairwise loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echarge2a" calculates second derivatives of the charge-charge
c     interaction energy for a single atom using a pairwise loop
c
c
      subroutine echarge2a (i)
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
      use group
      use hessn
      use shunt
      implicit none
      integer i,j,k,kk
      integer in,ic,kn,kc
      integer jcell
      real*8 e,de,d2e
      real*8 fi,fik,fgrp
      real*8 r,r2,rb,rb2
      real*8 d2edx,d2edy,d2edz
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 shift,taper,trans
      real*8 dtaper,dtrans
      real*8 d2taper,d2trans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8 term(3,3)
      real*8, allocatable :: cscale(:)
      logical proceed
      character*6 mode
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = electric * pchg(k) / dielec
            in = jion(k)
            ic = kion(k)
            goto 10
         end if
      end do
      return
   10 continue
c
c     store the coordinates of the atom of interest
c
      xic = x(ic)
      yic = y(ic)
      zic = z(ic)
      xi = x(i) - xic
      yi = y(i) - yic
      zi = z(i) - zic
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do j = 1, nion
         cscale(iion(j)) = 1.0d0
      end do
      do j = 1, n12(in)
         cscale(i12(j,in)) = c2scale
      end do
      do j = 1, n13(in)
         cscale(i13(j,in)) = c3scale
      end do
      do j = 1, n14(in)
         cscale(i14(j,in)) = c4scale
      end do
      do j = 1, n15(in)
         cscale(i15(j,in)) = c5scale
      end do
c
c     set cutoff distances and switching function coefficients
c
      mode = 'CHARGE'
      call switch (mode)
c
c     calculate the charge interaction energy Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         kc = kion(kk)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         if (proceed)  proceed = (kn .ne. i)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            xc = xic - x(kc)
            yc = yic - y(kc)
            zc = zic - z(kc)
            if (use_bounds)  call image (xc,yc,zc)
            rc2 = xc*xc + yc*yc + zc*zc
            if (rc2 .le. off2) then
               xr = xc + xi - x(k) + x(kc)
               yr = yc + yi - y(k) + y(kc)
               zr = zc + zi - z(k) + z(kc)
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
               rb = r + ebuffer
               rb2 = rb * rb
               fik = fi * pchg(kk) * cscale(kn)
c
c     compute chain rule terms for Hessian matrix elements
c
               de = -fik / rb2
               d2e = -2.0d0 * de/rb
c
c     use shifted energy switching if near the cutoff distance
c
               if (rc2 .gt. cut2) then
                  e = fik / rb
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  rc = sqrt(rc2)
                  rc3 = rc2 * rc
                  rc4 = rc2 * rc2
                  rc5 = rc2 * rc3
                  rc6 = rc3 * rc3
                  rc7 = rc3 * rc4
                  taper = c5*rc5 + c4*rc4 + c3*rc3
     &                       + c2*rc2 + c1*rc + c0
                  dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                        + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                  d2taper = 20.0d0*c5*rc3 + 12.0d0*c4*rc2
     &                         + 6.0d0*c3*rc + 2.0d0*c2
                  trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                            + f3*rc3 + f2*rc2 + f1*rc + f0)
                  dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                            + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                            + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                  d2trans = fik * (42.0d0*f7*rc5 + 30.0d0*f6*rc4
     &                             + 20.0d0*f5*rc3 + 12.0d0*f4*rc2
     &                             + 6.0d0*f3*rc + 2.0d0*f2)
                  d2e = e*d2taper + 2.0d0*de*dtaper
     &                     + d2e*taper + d2trans
                  de = e*dtaper + de*taper + dtrans
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  de = de * fgrp
                  d2e = d2e * fgrp
               end if
c
c     form the individual Hessian element components
c
               de = de / r
               d2e = (d2e-de) / r2
               d2edx = d2e * xr
               d2edy = d2e * yr
               d2edz = d2e * zr
               term(1,1) = d2edx*xr + de
               term(1,2) = d2edx*yr
               term(1,3) = d2edx*zr
               term(2,1) = term(1,2)
               term(2,2) = d2edy*yr + de
               term(2,3) = d2edy*zr
               term(3,1) = term(1,3)
               term(3,2) = term(2,3)
               term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + term(1,j)
                  hessy(j,i) = hessy(j,i) + term(2,j)
                  hessz(j,i) = hessz(j,i) + term(3,j)
                  hessx(j,k) = hessx(j,k) - term(1,j)
                  hessy(j,k) = hessy(j,k) - term(2,j)
                  hessz(j,k) = hessz(j,k) - term(3,j)
               end do
            end if
         end if
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         kc = kion(kk)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            do jcell = 2, ncell
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               call imager (xc,yc,zc,jcell)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rb2 = rb * rb
                  fik = fi * pchg(kk)
                  if (use_polymer) then
                     if (rc2 .le. polycut2)  fik = fik * cscale(kn)
                  end if
c
c     compute chain rule terms for Hessian matrix elements
c
                  de = -fik / rb2
                  d2e = -2.0d0 * de/rb
c
c     use shifted energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     e = fik / rb
                     shift = fik / (0.5d0*(off+cut))
                     e = e - shift
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                           + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                     d2taper = 20.0d0*c5*rc3 + 12.0d0*c4*rc2
     &                            + 6.0d0*c3*rc + 2.0d0*c2
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                               + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                               + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                     d2trans = fik * (42.0d0*f7*rc5 + 30.0d0*f6*rc4
     &                                + 20.0d0*f5*rc3 + 12.0d0*f4*rc2
     &                                + 6.0d0*f3*rc + 2.0d0*f2)
                     d2e = e*d2taper + 2.0d0*de*dtaper
     &                        + d2e*taper + d2trans
                     de = e*dtaper + de*taper + dtrans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     de = de * fgrp
                     d2e = d2e * fgrp
                  end if
c
c     form the individual Hessian element components
c
                  de = de / r
                  d2e = (d2e-de) / r2
                  d2edx = d2e * xr
                  d2edy = d2e * yr
                  d2edz = d2e * zr
                  term(1,1) = d2edx*xr + de
                  term(1,2) = d2edx*yr
                  term(1,3) = d2edx*zr
                  term(2,1) = term(1,2)
                  term(2,2) = d2edy*yr + de
                  term(2,3) = d2edy*zr
                  term(3,1) = term(1,3)
                  term(3,2) = term(2,3)
                  term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                  do j = 1, 3
                     hessx(j,i) = hessx(j,i) + term(1,j)
                     hessy(j,i) = hessy(j,i) + term(2,j)
                     hessz(j,i) = hessz(j,i) + term(3,j)
                     hessx(j,k) = hessx(j,k) - term(1,j)
                     hessy(j,k) = hessy(j,k) - term(2,j)
                     hessz(j,k) = hessz(j,k) - term(3,j)
                  end do
               end if
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine echarge2b  --  charge Hessian via neighbor list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echarge2b" calculates second derivatives of the charge-charge
c     interaction energy for a single atom using a neighbor list
c
c
      subroutine echarge2b (i)
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
      use group
      use hessn
      use neigh
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer in,ic,kn,kc
      integer jcell
      real*8 e,de,d2e
      real*8 fi,fik,fgrp
      real*8 r,r2,rb,rb2
      real*8 d2edx,d2edy,d2edz
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 shift,taper,trans
      real*8 dtaper,dtrans
      real*8 d2taper,d2trans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8 term(3,3)
      real*8, allocatable :: cscale(:)
      logical proceed
      character*6 mode
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = electric * pchg(k) / dielec
            ii = k
            in = jion(k)
            ic = kion(k)
            goto 10
         end if
      end do
      return
   10 continue
c
c     store the coordinates of the atom of interest
c
      xic = x(ic)
      yic = y(ic)
      zic = z(ic)
      xi = x(i) - xic
      yi = y(i) - yic
      zi = z(i) - zic
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do j = 1, nion
         cscale(iion(j)) = 1.0d0
      end do
      do j = 1, n12(in)
         cscale(i12(j,in)) = c2scale
      end do
      do j = 1, n13(in)
         cscale(i13(j,in)) = c3scale
      end do
      do j = 1, n14(in)
         cscale(i14(j,in)) = c4scale
      end do
      do j = 1, n15(in)
         cscale(i15(j,in)) = c5scale
      end do
c
c     set cutoff distances and switching function coefficients
c
      mode = 'CHARGE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(i,ii,iion,jion,kion,x,y,z,
!$OMP& xi,yi,zi,xic,yic,zic,fi,pchg,nelst,elst,use_group,use_bounds,
!$OMP& off,off2,cut,cut2,c0,c1,c2,c3,c4,c5,f0,f1,f2,f3,f4,f5,f6,f7,
!$OMP& ebuffer,cscale)
!$OMP& shared (hessx,hessy,hessz)
!$OMP DO reduction(+:hessx,hessy,hessz) schedule(guided)
c
c     calculate the charge interaction energy Hessian elements
c
      do kkk = 1, nelst(ii)
         kk = elst(kkk,ii)
         k = iion(kk)
         kn = jion(kk)
         kc = kion(kk)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         if (proceed)  proceed = (kn .ne. i)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            xc = xic - x(kc)
            yc = yic - y(kc)
            zc = zic - z(kc)
            if (use_bounds)  call image (xc,yc,zc)
            rc2 = xc*xc + yc*yc + zc*zc
            if (rc2 .le. off2) then
               xr = xc + xi - x(k) + x(kc)
               yr = yc + yi - y(k) + y(kc)
               zr = zc + zi - z(k) + z(kc)
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
               rb = r + ebuffer
               rb2 = rb * rb
               fik = fi * pchg(kk) * cscale(kn)
c
c     compute chain rule terms for Hessian matrix elements
c
               de = -fik / rb2
               d2e = -2.0d0 * de/rb
c
c     use shifted energy switching if near the cutoff distance
c
               if (r2 .gt. cut2) then
                  e = fik / rb
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  rc = sqrt(rc2)
                  rc3 = rc2 * rc
                  rc4 = rc2 * rc2
                  rc5 = rc2 * rc3
                  rc6 = rc3 * rc3
                  rc7 = rc3 * rc4
                  taper = c5*rc5 + c4*rc4 + c3*rc3
     &                       + c2*rc2 + c1*rc + c0
                  dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                        + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                  d2taper = 20.0d0*c5*rc3 + 12.0d0*c4*rc2
     &                         + 6.0d0*c3*rc + 2.0d0*c2
                  trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                            + f3*rc3 + f2*rc2 + f1*rc + f0)
                  dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                            + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                            + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                  d2trans = fik * (42.0d0*f7*rc5 + 30.0d0*f6*rc4
     &                             + 20.0d0*f5*rc3 + 12.0d0*f4*rc2
     &                             + 6.0d0*f3*rc + 2.0d0*f2)
                  d2e = e*d2taper + 2.0d0*de*dtaper
     &                     + d2e*taper + d2trans
                  de = e*dtaper + de*taper + dtrans
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  de = de * fgrp
                  d2e = d2e * fgrp
               end if
c
c     form the individual Hessian element components
c
               de = de / r
               d2e = (d2e-de) / r2
               d2edx = d2e * xr
               d2edy = d2e * yr
               d2edz = d2e * zr
               term(1,1) = d2edx*xr + de
               term(1,2) = d2edx*yr
               term(1,3) = d2edx*zr
               term(2,1) = term(1,2)
               term(2,2) = d2edy*yr + de
               term(2,3) = d2edy*zr
               term(3,1) = term(1,3)
               term(3,2) = term(2,3)
               term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + term(1,j)
                  hessy(j,i) = hessy(j,i) + term(2,j)
                  hessz(j,i) = hessz(j,i) + term(3,j)
                  hessx(j,k) = hessx(j,k) - term(1,j)
                  hessy(j,k) = hessy(j,k) - term(2,j)
                  hessz(j,k) = hessz(j,k) - term(3,j)
               end do
            end if
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge2c  --  loop Ewald sum charge Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge2c" calculates second derivatives of the charge-charge
c     interaction energy for a single atom using a particle mesh
c     Ewald summation and pairwise loop
c
c     note this version does not incorporate the reciprocal space
c     contribution to the Hessian calculation
c
c
      subroutine echarge2c (i)
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
      use ewald
      use group
      use hessn
      use limits
      use math
      implicit none
      integer i,j,k,kk
      integer in,kn,jcell
      real*8 fi,fik,fgrp
      real*8 r,r2,rb,rb2
      real*8 de,d2e
      real*8 d2edx,d2edy,d2edz
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rew,erfc,cut2
      real*8 erfterm,expterm
      real*8 scale
      real*8 term(3,3)
      real*8, allocatable :: cscale(:)
      logical proceed
      external erfc
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = electric * pchg(k) / dielec
            in = jion(k)
            cut2 = ewaldcut * ewaldcut
            goto 10
         end if
      end do
      return
   10 continue
c
c     store the coordinates of the atom of interest
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
c
c     set Ewald coefficient for electrostatic interactions
c
      aewald = aeewald
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do j = 1, nion
         cscale(iion(j)) = 1.0d0
      end do
      do j = 1, n12(in)
         cscale(i12(j,in)) = c2scale
      end do
      do j = 1, n13(in)
         cscale(i13(j,in)) = c3scale
      end do
      do j = 1, n14(in)
         cscale(i14(j,in)) = c4scale
      end do
      do j = 1, n15(in)
         cscale(i15(j,in)) = c5scale
      end do
c
c     calculate the real space Ewald interaction Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         proceed = .true.
         if (proceed)  proceed = (kn .ne. i)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. cut2) then
               r = sqrt(r2)
               rb = r + ebuffer
               rb2 = rb * rb
               fik = fi * pchg(kk) * cscale(kn)
               rew = aewald * r
               erfterm = erfc (rew)
               expterm = exp(-rew**2)
               scale = cscale(kn)
               if (use_group)  scale = scale * fgrp
               scale = scale - 1.0d0
c
c     compute chain rule terms for Hessian matrix elements
c
               de = -fik * ((erfterm+scale)/rb2
     &                        + (2.0d0*aewald/sqrtpi)*expterm/r)
               d2e = -2.0d0*de/rb + 2.0d0*(fik/(rb*rb2))*scale
     &                  + (4.0d0*fik*aewald**3/sqrtpi)*expterm
     &                  + 2.0d0*(fik/(rb*rb2))*scale
c
c     form the individual Hessian element components
c
               de = de / r
               d2e = (d2e-de) / r2
               d2edx = d2e * xr
               d2edy = d2e * yr
               d2edz = d2e * zr
               term(1,1) = d2edx*xr + de
               term(1,2) = d2edx*yr
               term(1,3) = d2edx*zr
               term(2,1) = term(1,2)
               term(2,2) = d2edy*yr + de
               term(2,3) = d2edy*zr
               term(3,1) = term(1,3)
               term(3,2) = term(2,3)
               term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
               do j = 1, 3
                  hessx(j,i) = hessx(j,i) + term(1,j)
                  hessy(j,i) = hessy(j,i) + term(2,j)
                  hessz(j,i) = hessz(j,i) + term(3,j)
                  hessx(j,k) = hessx(j,k) - term(1,j)
                  hessy(j,k) = hessy(j,k) - term(2,j)
                  hessz(j,k) = hessz(j,k) - term(3,j)
               end do
            end if
         end if
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         proceed = .true.
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            do jcell = 2, ncell
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call imager (xr,yr,zr,jcell)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. cut2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rb2 = rb * rb
                  fik = fi * pchg(kk)
                  rew = aewald * r
                  erfterm = erfc (rew)
                  expterm = exp(-rew**2)
                  scale = 1.0d0
                  if (use_group)  scale = scale * fgrp
                  if (use_polymer) then
                     if (r2 .le. polycut2) then
                        scale = scale * cscale(kn)
                     end if
                  end if
                  scale = scale - 1.0d0
c
c     compute chain rule terms for Hessian matrix elements
c
                  de = -fik * ((erfterm+scale)/rb2
     &                    + (2.0d0*aewald/sqrtpi)*exp(-rew**2)/r)
                  d2e = -2.0d0*de/rb + 2.0d0*(fik/(rb*rb2))*scale
     &                     + (4.0d0*fik*aewald**3/sqrtpi)*expterm
     &                     + 2.0d0*(fik/(rb*rb2))*scale
c
c     form the individual Hessian element components
c
                  de = de / r
                  d2e = (d2e-de) / r2
                  d2edx = d2e * xr
                  d2edy = d2e * yr
                  d2edz = d2e * zr
                  term(1,1) = d2edx*xr + de
                  term(1,2) = d2edx*yr
                  term(1,3) = d2edx*zr
                  term(2,1) = term(1,2)
                  term(2,2) = d2edy*yr + de
                  term(2,3) = d2edy*zr
                  term(3,1) = term(1,3)
                  term(3,2) = term(2,3)
                  term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
                  do j = 1, 3
                     hessx(j,i) = hessx(j,i) + term(1,j)
                     hessy(j,i) = hessy(j,i) + term(2,j)
                     hessz(j,i) = hessz(j,i) + term(3,j)
                     hessx(j,k) = hessx(j,k) - term(1,j)
                     hessy(j,k) = hessy(j,k) - term(2,j)
                     hessz(j,k) = hessz(j,k) - term(3,j)
                  end do
               end if
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge2d  --  charge Hessian for smoothing  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge2d" calculates second derivatives of the charge-charge
c     interaction energy for a single atom for use with potential
c     smoothing methods
c
c
      subroutine echarge2d (i)
      use atoms
      use charge
      use chgpot
      use couple
      use group
      use hessn
      use math
      use warp
      implicit none
      integer i,j,k,kk
      integer in,kn
      real*8 fi,fik,fgrp
      real*8 r,r2,rb,rb2
      real*8 de,d2e
      real*8 d2edx,d2edy,d2edz
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 erf,erfterm
      real*8 expcut,expterm
      real*8 wterm,width
      real*8 width2,width3
      real*8 term(3,3)
      real*8, allocatable :: cscale(:)
      logical proceed
      external erf
c
c
c     first see if the atom of interest carries a charge
c
      do k = 1, nion
         if (iion(k) .eq. i) then
            fi = electric * pchg(k) / dielec
            in = jion(k)
            goto 10
         end if
      end do
      return
   10 continue
c
c     store the coordinates of the atom of interest
c
      xi = x(i)
      yi = y(i)
      zi = z(i)
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do j = 1, nion
         cscale(iion(j)) = 1.0d0
      end do
      do j = 1, n12(in)
         cscale(i12(j,in)) = c2scale
      end do
      do j = 1, n13(in)
         cscale(i13(j,in)) = c3scale
      end do
      do j = 1, n14(in)
         cscale(i14(j,in)) = c4scale
      end do
      do j = 1, n15(in)
         cscale(i15(j,in)) = c5scale
      end do
c
c     set the smallest exponential terms to be calculated
c
      expcut = -50.0d0
c
c     set the extent of smoothing to be performed
c
      width = deform * diffc
      if (use_dem) then
         if (width .gt. 0.0d0)  width = 0.5d0 / sqrt(width)
      else if (use_gda) then
         wterm = sqrt(3.0d0/(2.0d0*diffc))
      end if
      width2 = width * width
      width3 = width * width2
c
c     calculate the charge interaction energy Hessian elements
c
      do kk = 1, nion
         k = iion(kk)
         kn = jion(kk)
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
         if (proceed)  proceed = (kn .ne. i)
c
c     compute the energy contribution for this interaction
c
         if (proceed) then
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            r2 = xr*xr + yr*yr + zr*zr
            r = sqrt(r2)
            rb = r + ebuffer
            rb2 = rb * rb
            fik = fi * pchg(kk) * cscale(kn)
c
c     compute chain rule terms for Hessian matrix elements
c
            de = -fik / rb2
            d2e = -2.0d0 * de/rb
c
c     transform the potential function via smoothing
c
            if (use_dem) then
               if (width .gt. 0.0d0) then
                  erfterm = erf(width*rb)
                  expterm = -rb2 * width2
                  if (expterm .gt. expcut) then
                     expterm = 2.0d0*fik*width*exp(expterm)
     &                            / (sqrtpi*rb)
                  else
                     expterm = 0.0d0
                  end if
                  de = de*erfterm + expterm
                  d2e = -2.0d0 * (de/rb + expterm*rb*width2)
               end if
            else if (use_gda) then
               width = m2(i) + m2(k)
               if (width .gt. 0.0d0) then
                  width = wterm / sqrt(width)
                  width2 = width * width
                  erfterm = erf(width*rb)
                  expterm = -rb2 * width2
                  if (expterm .gt. expcut) then
                     expterm = 2.0d0*fik*width*exp(expterm)
     &                            / (sqrtpi*rb)
                  else
                     expterm = 0.0d0
                  end if
                  de = de*erfterm + expterm
                  d2e = -2.0d0 * (de/rb + expterm*r*width2)
               end if
            else if (use_tophat) then
               if (width .gt. rb) then
                  d2e = -fik / width3
                  de = d2e * rb
               end if
            else if (use_stophat) then
               wterm = rb + width
               de = -fik / (wterm*wterm)
               d2e = -2.0d0 * de / wterm
            end if
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               de = de * fgrp
               d2e = d2e * fgrp
            end if
c
c     form the individual Hessian element components
c
            de = de / r
            d2e = (d2e-de) / r2
            d2edx = d2e * xr
            d2edy = d2e * yr
            d2edz = d2e * zr
            term(1,1) = d2edx*xr + de
            term(1,2) = d2edx*yr
            term(1,3) = d2edx*zr
            term(2,1) = term(1,2)
            term(2,2) = d2edy*yr + de
            term(2,3) = d2edy*zr
            term(3,1) = term(1,3)
            term(3,2) = term(2,3)
            term(3,3) = d2edz*zr + de
c
c     increment diagonal and non-diagonal Hessian elements
c
            do j = 1, 3
               hessx(j,i) = hessx(j,i) + term(1,j)
               hessy(j,i) = hessy(j,i) + term(2,j)
               hessz(j,i) = hessz(j,i) + term(3,j)
               hessx(j,k) = hessx(j,k) - term(1,j)
               hessy(j,k) = hessy(j,k) - term(2,j)
               hessz(j,k) = hessz(j,k) - term(3,j)
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c     ################################################################
c     ##                                                            ##
c     ##  program num_echarge2c - numerical Ewald Hessian summation ##
c     ##                                                            ##
c     ################################################################
c
c
c     "num_echarge2c" computes the numerical Hessian matrices of the
c     ewald energy summation with respect to Cartesian coordinates using
c     a finite difference approach
c
c     This subroutine is modfied after the testhess.f file. It uses the
c     part to run finite differences, and outputs a matrix that contains
c     ewald hessian elements to be added in.
      subroutine num_echarge2c (i)
      use atoms
      use inform
      use hessn
      implicit none
      integer i,j,k
      integer ii,jj
      integer index
      real*8 old,eps
      real*8 delta
      real*8, allocatable :: g(:,:)
      real*8, allocatable :: g0(:,:)
      character*120 xyzfile
c
c     set difference threshhold via the energy precision
c
      delta = 1.0d-5
      if (digits .ge. 6)  delta = 1.0d-6
      if (digits .ge. 8)  delta = 1.0d-8
c
c     setting the stepsize for numerical Hessian calculation
c
      eps = 1.0d-5
c
c     perform dynamic allocation of some local arrays
c
      allocate (g(3,n))
      allocate (g0(3,n))
c
c     get the two-sided numerical Hessian matrix elements
c
      old = x(i)
      x(i) = x(i) - 0.5d0*eps
      call gradient_echarge (g)
      do k = 1, n
         do j = 1, 3
            g0(j,k) = g(j,k)
         end do
      end do
      x(i) = x(i) + eps
      call gradient_echarge (g)
      x(i) = old
      do k = 1, n
         do j = 1, 3
            hessx(j,k) = hessx(j,k) + (g(j,k) - g0(j,k)) / eps
         end do
      end do
      old = y(i)
      y(i) = y(i) - 0.5d0*eps
      call gradient_echarge (g)
      do k = 1, n
         do j = 1, 3
            g0(j,k) = g(j,k)
         end do
      end do
      y(i) = y(i) + eps
      call gradient_echarge (g)
      y(i) = old
      do k = 1, n
         do j = 1, 3
            hessy(j,k) = hessy(j,k) + (g(j,k) - g0(j,k)) / eps
         end do
      end do
      old = z(i)
      z(i) = z(i) - 0.5d0*eps
      call gradient_echarge (g)
      do k = 1, n
         do j = 1, 3
            g0(j,k) = g(j,k)
         end do
      end do
      z(i) = z(i) + eps
      call gradient_echarge (g)
      z(i) = old
      do k = 1, n
         do j = 1, 3
            hessz(j,k) = hessz(j,k) + (g(j,k) - g0(j,k)) / eps
         end do
      end do
      deallocate (g)
      deallocate (g0)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine gradient_echarge -- find gradient components     ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "gradient_echarge" calls subroutines to calculate the potential energy
c     and first derivatives with respect to Cartesian coordinates
c
c     Modeled after the gradient.f file to return the gradient of the
c     first derivative of the ewald summation. Helps to support the
c     num_echarge2c subroutine above.
c
      subroutine gradient_echarge (derivs)
      use atoms
      use deriv
      use energi
      use potent
      implicit none
      integer l,m
      real*8 energy,cutoff
      real*8 derivs(3,*)
      logical prior
c
c
c     perform dynamic allocation of some global arrays
c
      prior = .false.
      if (allocated(dec)) then
         prior = .true.
         if (size(dec) .lt. 3*n)  deallocate (dec)
      end if
      if (.not. allocated(dec))  allocate (dec(3,n))
c
c     zero out each of the first derivative components
c
      do l = 1, n
         do m = 1, 3
            dec(m,l) = 0.0d0
         end do
      end do
c
c     call the electrostatic energy and gradient routines
c
      if (use_charge)  call echarge1
c
c     sum up to get the total energy and first derivatives
c
      esum = ec
      energy = esum
      do l = 1, n
         do m = 1, 3
            derivs(m,l) = dec(m,l)
         end do
      end do
c
c     deallocate some of the global arrays
c
      if (.not. prior) then
         deallocate (dec)
      end if
      return
      end


c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine edisp2  --  atomwise damped dispersion Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "edisp2" calculates the damped dispersion second derivatives
c     for a single atom at a time
c
c     literature reference:
c
c     J. A. Rackers, C. Liu, P. Ren and J. W. Ponder, "A Physically
c     Grounded Damped Dispersion Model with Particle Mesh Ewald
c     Summation", Journal of Chemical Physics, 149, 084115 (2018)
c
c
      subroutine edisp2 (iatom)
      use atoms
      use bound
      use cell
      use couple
      use disp
      use dsppot
      use group
      use hessn
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer iatom,jcell
      integer nlist,list(5)
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,de,d2e
      real*8 ci,ck,fgrp
      real*8 r,r2,r3
      real*8 r4,r5,r6
      real*8 ai,ai2
      real*8 ai3,ai4
      real*8 ak,ak2
      real*8 ak3,ak4
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 damp3,damp5
      real*8 damp,ddamp
      real*8 ddampi,ddampk
      real*8 d2damp
      real*8 taper,dtaper
      real*8 d2taper
      real*8 d2edx,d2edy,d2edz
      real*8 term(3,3)
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (dspscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         dspscale(i) = 1.0d0
      end do
c
c     set cutoff and switching coefficients
c
      mode = 'DISP'
      call switch (mode)
c
c     check to see if atom of interest is a dispersion site
c
      nlist = 0
      do k = 1, ndisp
         if (idisp(k) .eq. iatom) then
            nlist = nlist + 1
            list(nlist) = iatom
            goto 10
         end if
      end do
      return
   10 continue
c
c     calculate the dispersion energy Hessian elements
c
      do ii = 1, nlist
         i = list(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         ci = csix(i)
         ai = adisp(i)
         usei = use(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            dspscale(i12(j,i)) = dsp2scale
         end do
         do j = 1, n13(i)
            dspscale(i13(j,i)) = dsp3scale
         end do
         do j = 1, n14(i)
            dspscale(i14(j,i)) = dsp4scale
         end do
         do j = 1, n15(i)
            dspscale(i15(j,i)) = dsp5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = 1, ndisp
            k = idisp(kk)
            ck = csix(k)
            ak = adisp(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (k .ne. i)
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  r6 = r2**3
                  e = -ci * ck / r6
                  de = -6.0d0 * e / r
                  d2e = -7.0d0 * de / r
c
c     find the damping factor for the dispersion interaction
c
                  ai2 = ai * ai
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (ai .ne. ak) then
                     ai3 = ai * ai2
                     ai4 = ai2 * ai2
                     ak2 = ak * ak
                     ak3 = ak * ak2
                     ak4 = ak2 * ak2
                     dk2 = dk * dk
                     dk3 = dk * dk2
                     ti = ak2 / (ak2-ai2)
                     tk = ai2 / (ai2-ak2)
                     ti2 = ti * ti
                     tk2 = tk * tk
                     damp3 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2)*expi
     &                          - tk2*(1.0d0+dk+0.5d0*dk2)*expk
     &                          - 2.0d0*ti2*tk*(1.0d0+di)*expi
     &                          - 2.0d0*tk2*ti*(1.0d0+dk)*expk
                     damp5 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2
     &                                       +di3/6.0d0)*expi
     &                          - tk2*(1.0d0+dk+0.5d0*dk2
     &                                    +dk3/6.0d0)*expk
     &                          - 2.0d0*ti2*tk*(1.0+di+di2/3.0d0)*expi
     &                          - 2.0d0*tk2*ti*(1.0+dk+dk2/3.0d0)*expk
                     ddampi = 0.25d0 * di2 * ti2 * ai * expi
     &                           * (r*ai+4.0d0*tk-1.0d0)
                     ddampk = 0.25d0 * dk2 * tk2 * ak * expk
     &                           * (r*ak+4.0d0*ti-1.0d0)
                     ddamp = ddampi + ddampk
                     d2damp = 2.0d0*ddamp/r - ai*ddampi - ak*ddampk
     &                           + 0.25d0*di2*ti2*ai2*expi
     &                           + 0.25d0*dk2*tk2*ak2*expk
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                     ddamp = (di5-3.0d0*di3-3.0d0*di2)
     &                           *ai*expi/96.0d0
                     d2damp = (5.0d0*di4-9.0d0*di2-6.0d0*di)
     &                            *ai2*expi/96.0d0 - ai*ddamp
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  d2e = d2e*damp**2 + 4.0d0*de*damp*ddamp
     &                     + 2.0d0*e*(ddamp**2+damp*d2damp)
                  de = de*damp**2 + 2.0d0*e*damp*ddamp
                  e = e * damp**2
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                            + 6.0d0*c3*r + 2.0d0*c2
                     d2e = e*d2taper + 2.0d0*de*dtaper + d2e*taper
                     de = e*dtaper + de*taper
                  end if
c
c     scale the interaction based on its group membership
c
                  de = de * dspscale(k)
                  d2e = d2e * dspscale(k)
                  if (use_group) then
                     de = de * fgrp
                     d2e = d2e * fgrp
                  end if
c
c     get chain rule terms for damped dispersion Hessian elements
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
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            dspscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            dspscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            dspscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            dspscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nlist
         i = list(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         ci = csix(i)
         ai = adisp(i)
         usei = use(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            dspscale(i12(j,i)) = dsp2scale
         end do
         do j = 1, n13(i)
            dspscale(i13(j,i)) = dsp3scale
         end do
         do j = 1, n14(i)
            dspscale(i14(j,i)) = dsp4scale
         end do
         do j = 1, n15(i)
            dspscale(i15(j,i)) = dsp5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = 1, ndisp
            k = idisp(kk)
            ck = csix(k)
            ak = adisp(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed) then
               do jcell = 2, ncell
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  call image (xr,yr,zr)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     r6 = r2**3
                     e = -ci * ck / r6
                     de = -6.0d0 * e / r
                     d2e = -7.0d0 * de / r
c
c     find the damping factor for the dispersion interaction
c
                     ai2 = ai * ai
                     di = ai * r
                     di2 = di * di
                     di3 = di * di2
                     dk = ak * r
                     expi = exp(-di)
                     expk = exp(-dk)
                     if (ai .ne. ak) then
                        ai3 = ai * ai2
                        ai4 = ai2 * ai2
                        ak2 = ak * ak
                        ak3 = ak * ak2
                        ak4 = ak2 * ak2
                        dk2 = dk * dk
                        dk3 = dk * dk2
                        ti = ak2 / (ak2-ai2)
                        tk = ai2 / (ai2-ak2)
                        ti2 = ti * ti
                        tk2 = tk * tk
                        damp3 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2)*expi
     &                             - tk2*(1.0d0+dk+0.5d0*dk2)*expk
     &                             - 2.0d0*ti2*tk*(1.0d0+di)*expi
     &                             - 2.0d0*tk2*ti*(1.0d0+dk)*expk
                        damp5 = 1.0d0 - ti2*(1.0d0+di+0.5d0*di2
     &                                          +di3/6.0d0)*expi
     &                             - tk2*(1.0d0+dk+0.5d0*dk2
     &                                       +dk3/6.0d0)*expk
     &                            - 2.0d0*ti2*tk*(1.0+di+di2/3.0d0)*expi
     &                            - 2.0d0*tk2*ti*(1.0+dk+dk2/3.0d0)*expk
                        ddampi = 0.25d0 * di2 * ti2 * ai * expi
     &                              * (r*ai+4.0d0*tk-1.0d0)
                        ddampk = 0.25d0 * dk2 * tk2 * ak * expk
     &                              * (r*ak+4.0d0*ti-1.0d0)
                        ddamp = ddampi + ddampk
                        d2damp = 2.0d0*ddamp/r - ai*ddampi - ak*ddampk
     &                              + 0.25d0*di2*ti2*ai2*expi
     &                              + 0.25d0*dk2*tk2*ak2*expk
                     else
                        di4 = di2 * di2
                        di5 = di2 * di3
                        damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                             +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                        damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                           +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                        ddamp = (di5-3.0d0*di3-3.0d0*di2)
     &                              *ai*expi/96.0d0
                        d2damp = (5.0d0*di4-9.0d0*di2-6.0d0*di)
     &                               *ai2*expi/96.0d0 - ai*ddamp
                     end if
                     damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                     d2e = d2e*damp**2 + 4.0d0*de*damp*ddamp
     &                        + 2.0d0*e*(ddamp**2+damp*d2damp)
                     de = de*damp**2 + 2.0d0*e*damp*ddamp
                     e = e * damp**2
c
c     use energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                               + 6.0d0*c3*r + 2.0d0*c2
                        d2e = e*d2taper + 2.0d0*de*dtaper + d2e*taper
                        de = e*dtaper + de*taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           de = de * dspscale(k)
                           d2e = d2e * dspscale(k)
                        end if
                     end if
                     if (use_group) then
                        de = de * fgrp
                        d2e = d2e * fgrp
                     end if
                     if (i .eq. k) then
                        de = 0.5d0 * de
                        d2e = 0.5d0 * d2e
                     end if
c
c     get chain rule terms for damped dispersion Hessian elements
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
                        hessx(j,ii) = hessx(j,ii) + term(1,j)
                        hessy(j,ii) = hessy(j,ii) + term(2,j)
                        hessz(j,ii) = hessz(j,ii) + term(3,j)
                        hessx(j,kk) = hessx(j,kk) - term(1,j)
                        hessy(j,kk) = hessy(j,kk) - term(2,j)
                        hessz(j,kk) = hessz(j,kk) - term(3,j)
                     end do
                  end if
               end do
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            dspscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            dspscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            dspscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            dspscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dspscale)
      return
      end

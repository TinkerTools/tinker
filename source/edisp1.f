c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edisp1  --  damped dispersion energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edisp1" calculates the damped dispersion energy and first
c     derivatives with respect to Cartesian coordinates
c
c     literature reference:
c
c     J. A. Rackers, C. Liu, P. Ren and J. W. Ponder, "A Physically
c     Grounded Damped Dispersion Model with Particle Mesh Ewald
c     Summation", Journal of Chemical Physics, 149, 084115 (2018)
c
c
      subroutine edisp1
      use dsppot
      use energi
      use ewald
      use limits
      use virial
      implicit none
      real*8 elrc,vlrc
      character*6 mode
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_dewald) then
         if (use_dlist) then
            call edisp1d
         else
            call edisp1c
         end if
      else
         if (use_dlist) then
            call edisp1b
         else
            call edisp1a
         end if
      end if
c
c     apply long range dispersion correction if desired
c
      if (use_dcorr .and. .not.use_dewald) then
         mode = 'DISP'
         call evcorr1 (mode,elrc,vlrc)
         edsp = edsp + elrc
         vir(1,1) = vir(1,1) + vlrc
         vir(2,2) = vir(2,2) + vlrc
         vir(3,3) = vir(3,3) + vlrc
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine edisp1a  --  double loop dispersion derivs  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "edisp1a" calculates the damped dispersion energy and
c     derivatives with respect to Cartesian coordinates using
c     a pairwise double loop
c
c
      subroutine edisp1a
      use atoms
      use bound
      use cell
      use couple
      use deriv
      use disp
      use dsppot
      use energi
      use group
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,de,fgrp
      real*8 ci,ck
      real*8 r,r2,r3
      real*8 r4,r5,r6
      real*8 ai,ai2,ai3
      real*8 ak,ak2,ak3
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 damp3,damp5
      real*8 damp,ddamp
      real*8 taper,dtaper
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      character*6 mode
c
c
c     zero out the damped dispersion energy and derivatives
c
      edsp = 0.0d0
      do i = 1, n
         dedsp(1,i) = 0.0d0
         dedsp(2,i) = 0.0d0
         dedsp(3,i) = 0.0d0
      end do
      if (ndisp .eq. 0) return
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
c     set conversion factor, cutoff and switching coefficients
c
      mode = 'DISP'
      call switch (mode)
c
c     find dispersion energy and derivatives via double loop
c
      do ii = 1, ndisp-1
         i = idisp(ii)
         ci = csix(ii)
         ai = adisp(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
c     decide whether to compute the current interaction
c
         do kk = ii+1, ndisp
            k = idisp(kk)
            ck = csix(kk)
            ak = adisp(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  r6 = r2**3
                  e = -ci * ck / r6
                  de = -6.0d0 * e / r
c
c     find the damping factor for the dispersion interaction
c
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (ai .ne. ak) then
                     ai2 = ai * ai
                     ai3 = ai * ai2
                     ak2 = ak * ak
                     ak3 = ak * ak2
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
                     ddamp = 0.25d0 * di2 * ti2 * ai * expi
     &                          * (r*ai+4.0d0*tk-1.0d0)
     &                       + 0.25d0 * dk2 * tk2 * ak * expk
     &                            * (r*ak+4.0d0*ti-1.0d0)
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                     ddamp = ai * expi * (di5-3.0d0*di3-3.0d0*di2)
     &                          / 96.0d0
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  de = de*damp**2 + 2.0d0*e*damp*ddamp
                  e = e * damp**2
                  e = e * dspscale(k)
                  de = de * dspscale(k)
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
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
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     increment the overall damped dispersion energy component
c
                  edsp = edsp + e
c
c     increment the damped dispersion derivative components
c
                  dedx = de * xr/r
                  dedy = de * yr/r
                  dedz = de * zr/r
                  dedsp(1,i) = dedsp(1,i) + dedx
                  dedsp(2,i) = dedsp(2,i) + dedy
                  dedsp(3,i) = dedsp(3,i) + dedz
                  dedsp(1,k) = dedsp(1,k) - dedx
                  dedsp(2,k) = dedsp(2,k) - dedy
                  dedsp(3,k) = dedsp(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
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
      do ii = 1, ndisp
         i = idisp(ii)
         ci = csix(ii)
         ai = adisp(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
c     decide whether to compute the current interaction
c
         do kk = ii, ndisp
            k = idisp(kk)
            ck = csix(kk)
            ak = adisp(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
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
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     r6 = r2**3
                     e = -ci * ck / r6
                     de = -6.0d0 * e / r
c
c     find the damping factor for the dispersion interaction
c
                     di = ai * r
                     di2 = di * di
                     di3 = di * di2
                     dk = ak * r
                     expi = exp(-di)
                     expk = exp(-dk)
                     if (ai .ne. ak) then
                        ai2 = ai * ai
                        ai3 = ai * ai2
                        ak2 = ak * ak
                        ak3 = ak * ak2
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
                        ddamp = 0.25d0 * di2 * ti2 * ai * expi
     &                             * (r*ai+4.0d0*tk-1.0d0)
     &                          + 0.25d0 * dk2 * tk2 * ak * expk
     &                               * (r*ak+4.0d0*ti-1.0d0)
                     else
                        di4 = di2 * di2
                        di5 = di2 * di3
                        damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                             +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                        damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                           +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                        ddamp = ai * expi * (di5-3.0d0*di3-3.0d0*di2)
     &                             / 96.0d0
                     end if
                     damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                     de = de*damp**2 + 2.0d0*e*damp*ddamp
                     e = e * damp**2
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           e = e * dspscale(k)
                           de = de * dspscale(k)
                        end if
                     end if
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                     end if
                     if (ii .eq. kk) then
                        e = 0.5d0 * e
                        de = 0.5d0 * de
                     end if
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
                        de = e*dtaper + de*taper
                        e = e * taper
                     end if
c
c     increment the overall damped dispersion energy component
c
                     edsp = edsp + e
c
c     increment the damped dispersion derivative components
c
                     dedx = de * xr/r
                     dedy = de * yr/r
                     dedz = de * zr/r
                     dedsp(1,i) = dedsp(1,i) + dedx
                     dedsp(2,i) = dedsp(2,i) + dedy
                     dedsp(3,i) = dedsp(3,i) + dedz
                     dedsp(1,k) = dedsp(1,k) - dedx
                     dedsp(2,k) = dedsp(2,k) - dedy
                     dedsp(3,k) = dedsp(3,k) - dedz
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
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
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine edisp1b  --  neighbor list dispersion derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "edisp1b" calculates the damped dispersion energy and
c     derivatives with respect to Cartesian coordinates using
c     a pairwise neighbor list
c
c
      subroutine edisp1b
      use atoms
      use bound
      use cell
      use couple
      use deriv
      use disp
      use dsppot
      use energi
      use group
      use neigh
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,de,fgrp
      real*8 ci,ck
      real*8 r,r2,r3
      real*8 r4,r5,r6
      real*8 ai,ai2,ai3
      real*8 ak,ak2,ak3
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 damp3,damp5
      real*8 damp,ddamp
      real*8 taper,dtaper
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      character*6 mode
c
c
c     zero out the damped dispersion energy and derivatives
c
      edsp = 0.0d0
      do i = 1, n
         dedsp(1,i) = 0.0d0
         dedsp(2,i) = 0.0d0
         dedsp(3,i) = 0.0d0
      end do
      if (ndisp .eq. 0) return
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
c     set conversion factor, cutoff and switching coefficients
c
      mode = 'DISP'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ndisp,idisp,csix,adisp,use,
!$OMP& x,y,z,n12,n13,n14,n15,i12,i13,i14,i15,nvlst,vlst,use_group,
!$OMP& dsp2scale,dsp3scale,dsp4scale,dsp5scale,off2,cut2)
!$OMP& firstprivate(dspscale) shared(edsp,dedsp,vir)
!$OMP DO reduction(+:edsp,dedsp,vir) schedule(guided)
c
c     find dispersion energy and derivatives via neighbor list
c
      do ii = 1, ndisp-1
         i = idisp(ii)
         ci = csix(ii)
         ai = adisp(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
c     decide whether to compute the current interaction
c
         do kkk = 1, nvlst(ii)
            kk = vlst(kkk,ii)
            k = idisp(kk)
            ck = csix(kk)
            ak = adisp(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  r6 = r2**3
                  e = -ci * ck / r6
                  de = -6.0d0 * e / r
c
c     find the damping factor for the dispersion interaction
c
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (ai .ne. ak) then
                     ai2 = ai * ai
                     ai3 = ai * ai2
                     ak2 = ak * ak
                     ak3 = ak * ak2
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
                     ddamp = 0.25d0 * di2 * ti2 * ai * expi
     &                          * (r*ai+4.0d0*tk-1.0d0)
     &                       + 0.25d0 * dk2 * tk2 * ak * expk
     &                            * (r*ak+4.0d0*ti-1.0d0)
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                     ddamp = ai * expi * (di5-3.0d0*di3-3.0d0*di2)
     &                          / 96.0d0
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  de = de*damp**2 + 2.0d0*e*damp*ddamp
                  e = e * damp**2
                  e = e * dspscale(k)
                  de = de * dspscale(k)
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
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
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     increment the overall damped dispersion energy component
c
                  edsp = edsp + e
c
c     increment the damped dispersion derivative components
c
                  dedx = de * xr/r
                  dedy = de * yr/r
                  dedz = de * zr/r
                  dedsp(1,i) = dedsp(1,i) + dedx
                  dedsp(2,i) = dedsp(2,i) + dedy
                  dedsp(3,i) = dedsp(3,i) + dedz
                  dedsp(1,k) = dedsp(1,k) - dedx
                  dedsp(2,k) = dedsp(2,k) - dedy
                  dedsp(3,k) = dedsp(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (dspscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edisp1c  --  Ewald dispersion derivs via loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edisp1c" calculates the damped dispersion energy and
c     derivatives with respect to Cartesian coordinates using
c     particle mesh Ewald summation and a double loop
c
c
      subroutine edisp1c
      use atoms
      use deriv
      use disp
      use energi
      use ewald
      use pme
      implicit none
      integer i
      real*8 term
c
c
c     zero out the damped dispersion energy and derivatives
c
      edsp = 0.0d0
      do i = 1, n
         dedsp(1,i) = 0.0d0
         dedsp(2,i) = 0.0d0
         dedsp(3,i) = 0.0d0
      end do
      if (ndisp .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = ndfft1
      nfft2 = ndfft2
      nfft3 = ndfft3
      bsorder = bsdorder
      aewald = adewald
c
c     compute the real space portion of the Ewald summation
c
      call edreal1c
c
c     compute the reciprocal space part of the Ewald summation
c
      call edrecip1
c
c     compute the self-energy portion of the Ewald summation
c
      do i = 1, ndisp
         term = aewald**6 / 12.0d0
         edsp = edsp + term*csix(i)*csix(i)
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edreal1c  --  Ewald real disp derivs via loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to damped dispersion
c     interactions via a double loop
c
c
      subroutine edreal1c
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use dsppot
      use deriv
      use energi
      use ewald
      use group
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,de,fgrp
      real*8 ci,ck
      real*8 r,r2,r6,r7
      real*8 ai,ai2
      real*8 ak,ak2
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 damp3,damp5
      real*8 damp,ddamp
      real*8 ralpha2,scale
      real*8 expterm,term
      real*8 expa,rterm
      real*8 dedx,dedy,dedz
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
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
c     set conversion factor, cutoff and switching coefficients
c
      mode = 'DEWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, ndisp-1
         i = idisp(ii)
         ci = csix(ii)
         ai = adisp(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
c     decide whether to compute the current interaction
c
         do kk = ii+1, ndisp
            k = idisp(kk)
            ck = csix(kk)
            ak = adisp(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r6 = r2**3
                  ralpha2 = r2 * aewald**2
                  term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                  expterm = exp(-ralpha2)
                  expa = expterm * term
c
c     find the damping factor for the dispersion interaction
c
                  r = sqrt(r2)
                  r7 = r6 * r
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (ai .ne. ak) then
                     ai2 = ai * ai
                     ak2 = ak * ak
                     dk2 = dk * dk
                     dk3 = dk * dk2
                     ti = ak2 / (ak2-ai2)
                     ti2 = ti * ti
                     tk = ai2 / (ai2-ak2)
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
                     ddamp = 0.25d0 * di2 * ti2 * ai * expi
     &                          * (r*ai+4.0d0*tk-1.0d0)
     &                       + 0.25d0 * dk2 * tk2 * ak * expk
     &                            * (r*ak+4.0d0*ti-1.0d0)
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                     ddamp = ai * expi * (di5-3.0d0*di3-3.0d0*di2)
     &                          / 96.0d0
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  scale = dspscale(k) * damp**2
                  if (use_group)  scale = scale * fgrp
                  scale = scale - 1.0d0
                  e = -ci * ck * (expa+scale) / r6
                  rterm = -(ralpha2**3) * expterm / r
                  de = -6.0d0*e/r2 - ci*ck*rterm/r7
     &                    - 2.0d0*ci*ck*dspscale(k)*damp*ddamp/r7
c
c     increment the overall damped dispersion energy component
c
                  edsp = edsp + e
c
c     increment the damped dispersion derivative components
c
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  dedsp(1,i) = dedsp(1,i) + dedx
                  dedsp(2,i) = dedsp(2,i) + dedy
                  dedsp(3,i) = dedsp(3,i) + dedz
                  dedsp(1,k) = dedsp(1,k) - dedx
                  dedsp(2,k) = dedsp(2,k) - dedy
                  dedsp(3,k) = dedsp(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
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
      do ii = 1, ndisp
         i = idisp(ii)
         ci = csix(ii)
         ai = adisp(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
c     decide whether to compute the current interaction
c
         do kk = ii, ndisp
            k = idisp(kk)
            ck = csix(kk)
            ak = adisp(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
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
                  if (r2 .le. off2) then
                     r6 = r2**3
                     ralpha2 = r2 * aewald**2
                     term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                     expterm = exp(-ralpha2)
                     expa = expterm * term
c
c     find the damping factor for the dispersion interaction
c
                     r = sqrt(r2)
                     r7 = r6 * r
                     di = ai * r
                     di2 = di * di
                     di3 = di * di2
                     dk = ak * r
                     expi = exp(-di)
                     expk = exp(-dk)
                     if (ai .ne. ak) then
                        ai2 = ai * ai
                        ak2 = ak * ak
                        dk2 = dk * dk
                        dk3 = dk * dk2
                        ti = ak2 / (ak2-ai2)
                        ti2 = ti * ti
                        tk = ai2 / (ai2-ak2)
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
                        ddamp = 0.25d0 * di2 * ti2 * ai * expi
     &                             * (r*ai+4.0d0*tk-1.0d0)
     &                          + 0.25d0 * dk2 * tk2 * ak * expk
     &                               * (r*ak+4.0d0*ti-1.0d0)
                     else
                        di4 = di2 * di2
                        di5 = di2 * di3
                        damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                             +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                        damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                           +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                        ddamp = ai * expi * (di5-3.0d0*di3-3.0d0*di2)
     &                             / 96.0d0
                     end if
                     damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                     scale = dspscale(k) * damp**2
                     if (use_group)  scale = scale * fgrp
                     scale = scale - 1.0d0
                     e = -ci * ck * (expa+scale) / r6
                     rterm = -(ralpha2**3) * expterm / r
                     de = -6.0d0*e/r2 - ci*ck*rterm/r7
     &                       - 2.0d0*ci*ck*dspscale(k)*damp*ddamp/r7
                     if (ii .eq. kk) then
                        e = 0.5d0 * e
                        de = 0.5d0 * de
                     end if
c
c     increment the overall damped dispersion energy component
c
                     edsp = edsp + e
c
c     increment the damped dispersion derivative components
c
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     dedsp(1,i) = dedsp(1,i) + dedx
                     dedsp(2,i) = dedsp(2,i) + dedy
                     dedsp(3,i) = dedsp(3,i) + dedz
                     dedsp(1,k) = dedsp(1,k) - dedx
                     dedsp(2,k) = dedsp(2,k) - dedy
                     dedsp(3,k) = dedsp(3,k) - dedz
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     vir(1,1) = vir(1,1) + vxx
                     vir(2,1) = vir(2,1) + vyx
                     vir(3,1) = vir(3,1) + vzx
                     vir(1,2) = vir(1,2) + vyx
                     vir(2,2) = vir(2,2) + vyy
                     vir(3,2) = vir(3,2) + vzy
                     vir(1,3) = vir(1,3) + vzx
                     vir(2,3) = vir(2,3) + vzy
                     vir(3,3) = vir(3,3) + vzz
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
      stop
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edisp1d  --  Ewald dispersion derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edisp1d" calculates the damped dispersion energy and
c     derivatives with respect to Cartesian coordinates using
c     particle mesh Ewald summation and a neighbor list
c
c
      subroutine edisp1d
      use atoms
      use deriv
      use disp
      use energi
      use ewald
      use pme
      implicit none
      integer i
      real*8 term
c
c
c     zero out the damped dispersion energy and derivatives
c
      edsp = 0.0d0
      do i = 1, n
         dedsp(1,i) = 0.0d0
         dedsp(2,i) = 0.0d0
         dedsp(3,i) = 0.0d0
      end do
      if (ndisp .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = ndfft1
      nfft2 = ndfft2
      nfft3 = ndfft3
      bsorder = bsdorder
      aewald = adewald
c
c     compute the real space portion of the Ewald summation
c
      call edreal1d
c
c     compute the reciprocal space part of the Ewald summation
c
      call edrecip1
c
c     compute the self-energy portion of the Ewald summation
c
      do i = 1, ndisp
         term = aewald**6 / 12.0d0
         edsp = edsp + term*csix(i)*csix(i)
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edreal1d  --  Ewald real disp derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edreal1d" evaluates the real space portion of the Ewald
c     summation energy and gradient due to damped dispersion
c     interactions via a neighbor list
c
c
      subroutine edreal1d
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use deriv
      use disp
      use dsppot
      use energi
      use ewald
      use group
      use neigh
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,de,fgrp
      real*8 ci,ck
      real*8 r,r2,r6,r7
      real*8 ai,ai2
      real*8 ak,ak2
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 damp3,damp5
      real*8 damp,ddamp
      real*8 ralpha2,scale
      real*8 expterm,term
      real*8 expa,rterm
      real*8 dedx,dedy,dedz
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
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
c     set conversion factor, cutoff and switching coefficients
c
      mode = 'DEWALD'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ndisp,idisp,csix,adisp,use,
!$OMP& x,y,z,n12,n13,n14,n15,i12,i13,i14,i15,nvlst,vlst,use_group,
!$OMP& dsp2scale,dsp3scale,dsp4scale,dsp5scale,off2,aewald)
!$OMP& firstprivate(dspscale) shared(edsp,dedsp,vir)
!$OMP DO reduction(+:edsp,dedsp,vir) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, ndisp
         i = idisp(ii)
         ci = csix(ii)
         ai = adisp(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
c     decide whether to compute the current interaction
c
         do kkk = 1, nvlst(ii)
            kk = vlst(kkk,ii)
            k = idisp(kk)
            ck = csix(kk)
            ak = adisp(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r6 = r2**3
                  ralpha2 = r2 * aewald**2
                  term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                  expterm = exp(-ralpha2)
                  expa = expterm * term
c
c     find the damping factor for the dispersion interaction
c
                  r = sqrt(r2)
                  r7 = r6 * r
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (ai .ne. ak) then
                     ai2 = ai * ai
                     ak2 = ak * ak
                     dk2 = dk * dk
                     dk3 = dk * dk2
                     ti = ak2 / (ak2-ai2)
                     ti2 = ti * ti
                     tk = ai2 / (ai2-ak2)
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
                     ddamp = 0.25d0 * di2 * ti2 * ai * expi
     &                          * (r*ai+4.0d0*tk-1.0d0)
     &                       + 0.25d0 * dk2 * tk2 * ak * expk
     &                            * (r*ak+4.0d0*ti-1.0d0)
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                     ddamp = ai * expi * (di5-3.0d0*di3-3.0d0*di2)
     &                          / 96.0d0
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  scale = dspscale(k) * damp**2
                  if (use_group)  scale = scale * fgrp
                  scale = scale - 1.0d0
                  e = -ci * ck * (expa+scale) / r6
                  rterm = -(ralpha2**3) * expterm / r
                  de = -6.0d0*e/r2 - ci*ck*rterm/r7
     &                    - 2.0d0*ci*ck*dspscale(k)*damp*ddamp/r7
c
c     increment the overall damped dispersion energy component
c
                  edsp = edsp + e
c
c     increment the damped dispersion derivative components
c
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  dedsp(1,i) = dedsp(1,i) + dedx
                  dedsp(2,i) = dedsp(2,i) + dedy
                  dedsp(3,i) = dedsp(3,i) + dedz
                  dedsp(1,k) = dedsp(1,k) - dedx
                  dedsp(2,k) = dedsp(2,k) - dedy
                  dedsp(3,k) = dedsp(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  vir(1,1) = vir(1,1) + vxx
                  vir(2,1) = vir(2,1) + vyx
                  vir(3,1) = vir(3,1) + vzx
                  vir(1,2) = vir(1,2) + vyx
                  vir(2,2) = vir(2,2) + vyy
                  vir(3,2) = vir(3,2) + vzy
                  vir(1,3) = vir(1,3) + vzx
                  vir(2,3) = vir(2,3) + vzy
                  vir(3,3) = vir(3,3) + vzz
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (dspscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine edrecip1  --  PME recip disp energy & derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "edrecip1" evaluates the reciprocal space portion of particle
c     mesh Ewald energy and gradient due to damped dispersion
c
c
      subroutine edrecip1
      use sizes
      use boxes
      use bound
      use disp
      use deriv
      use energi
      use ewald
      use math
      use pme
      use virial
      implicit none
      integer i,j,k,ii
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff,ntot
      integer i0,iatm,igrd0
      integer it1,it2,it3
      integer j0,jgrd0
      integer k0,kgrd0
      real*8 e,fi,denom
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 term,vterm
      real*8 expterm
      real*8 erfcterm
      real*8 hsq,struc2
      real*8 h,hhh,b,bfac
      real*8 term1,denom0
      real*8 fac1,fac2,fac3
      real*8 de1,de2,de3
      real*8 dn1,dn2,dn3
      real*8 dt1,dt2,dt3
      real*8 t1,t2,t3
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some global arrays
c
      ntot = nfft1 * nfft2 * nfft3
      if (allocated(qgrid)) then
         if (size(qgrid) .ne. 2*ntot) then
            call fftclose
            deallocate (qgrid)
         end if
      end if
      if (.not. allocated(qgrid)) then
         allocate (qgrid(2,nfft1,nfft2,nfft3))
         call fftsetup
      end if
c
c     setup spatial decomposition and B-spline coefficients
c
      call getchunk
      call moduli
      call bspline_fill
      call table_fill
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_disp
      call fftfront
c
c     use scalar sum to get the reciprocal space energy
c
      qgrid(1,1,1,1) = 0.0d0
      qgrid(2,1,1,1) = 0.0d0
      bfac = pi / aewald
      fac1 = 2.0d0*pi**(3.5d0)
      fac2 = aewald**3
      fac3 = -2.0d0*aewald*pi**2
      denom0 = (6.0d0*volbox)/(pi**1.5d0)
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      nff = nfft1 * nfft2
      ntot = nff * nfft3
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/ndfft1 + 1
         k1 = j - (k2-1)*ndfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - ndfft1
         if (k2 .gt. nf2)  m2 = m2 - ndfft2
         if (k3 .gt. nf3)  m3 = m3 - ndfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         h = sqrt(hsq)
         b = h*bfac
         hhh = h*hsq
         term = -b*b
         expterm = 0.0d0
         erfcterm = erfc(b)
         denom = denom0*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
         if (term .gt. -50.0d0) then
            expterm = exp(term)
            erfcterm = erfc(b)
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               erfcterm = erfcterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               if (mod(m1+m2+m3,2) .ne. 0)  erfcterm = 0.0d0
            end if
            term1 = fac1*erfcterm*hhh + expterm*(fac2 + fac3*hsq)
            struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
            e = -(term1 / denom) * struc2
            edsp = edsp + e
            vterm = 3.0d0*(fac1*erfcterm*h + fac3*expterm)*struc2/denom
            vir(1,1) = vir(1,1) + h1*h1*vterm - e
            vir(2,1) = vir(2,1) + h1*h2*vterm
            vir(3,1) = vir(3,1) + h1*h3*vterm
            vir(1,2) = vir(1,2) + h2*h1*vterm
            vir(2,2) = vir(2,2) + h2*h2*vterm - e
            vir(3,2) = vir(3,2) + h2*h3*vterm
            vir(1,3) = vir(1,3) + h3*h1*vterm
            vir(2,3) = vir(2,3) + h3*h2*vterm
            vir(3,3) = vir(3,3) + h3*h3*vterm - e
         end if
         qgrid(1,k1,k2,k3) = -(term1/denom) * qgrid(1,k1,k2,k3) 
         qgrid(2,k1,k2,k3) = -(term1/denom) * qgrid(2,k1,k2,k3)
      end do
c
c     perform the 3-D FFT backward transformation
c
      call fftback
c
c     get first derivatives of the reciprocal space energy 
c
      dn1 = dble(ndfft1)
      dn2 = dble(ndfft2)
      dn3 = dble(ndfft3)
      do ii = 1, ndisp
         iatm = idisp(ii)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
         fi = csix(ii)
         de1 = 0.0d0
         de2 = 0.0d0
         de3 = 0.0d0
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (ndfft3-sign(ndfft3,k0))/2
            t3 = thetai3(1,it3,iatm)
            dt3 = dn3 * thetai3(2,it3,iatm)
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (ndfft2-sign(ndfft2,j0))/2
               t2 = thetai2(1,it2,iatm)
               dt2 = dn2 * thetai2(2,it2,iatm)
               i0 = igrd0
               do it1 = 1, bsorder
                  i0 = i0 + 1
                  i = i0 + 1 + (ndfft1-sign(ndfft1,i0))/2
                  t1 = thetai1(1,it1,iatm)
                  dt1 = dn1 * thetai1(2,it1,iatm)
                  term = qgrid(1,i,j,k)
                  de1 = de1 + 2.0d0*term*dt1*t2*t3
                  de2 = de2 + 2.0d0*term*dt2*t1*t3
                  de3 = de3 + 2.0d0*term*dt3*t1*t2
               end do
            end do
         end do
         dedsp(1,iatm) = dedsp(1,iatm) + fi*(recip(1,1)*de1
     &                            +recip(1,2)*de2+recip(1,3)*de3)
         dedsp(2,iatm) = dedsp(2,iatm) + fi*(recip(2,1)*de1
     &                            +recip(2,2)*de2+recip(2,3)*de3)
         dedsp(3,iatm) = dedsp(3,iatm) + fi*(recip(3,1)*de1
     &                            +recip(3,2)*de2+recip(3,3)*de3)
      end do
c
c     account for the energy and virial correction terms
c
      term = csixpr * aewald**3 / denom0
      edsp = edsp - term
      vir(1,1) = vir(1,1) + term
      vir(2,2) = vir(2,2) + term
      vir(3,3) = vir(3,3) + term
      return
      end

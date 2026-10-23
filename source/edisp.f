c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edisp  --  damped dispersion potential energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edisp" calculates the damped dispersion potential energy
c
c     literature reference:
c
c     J. A. Rackers, C. Liu, P. Ren and J. W. Ponder, "A Physically
c     Grounded Damped Dispersion Model with Particle Mesh Ewald
c     Summation", Journal of Chemical Physics, 149, 084115 (2018)
c
c
      subroutine edisp
      use dsppot
      use energi
      use ewald
      use limits
      implicit none
      real*8 elrc
      character*6 mode
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_dewald) then
         if (use_dlist) then
            call edisp0d
         else
            call edisp0c
         end if
      else
         if (use_dlist) then
            call edisp0b
         else
            call edisp0a
         end if
      end if
c
c     apply long range dispersion correction if desired
c
      if (use_dcorr .and. .not.use_dewald) then
         mode = 'DISP'
         call evcorr (mode,elrc)
         edsp = edsp + elrc
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine edisp0a  --  damped dispersion via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "edisp0a" calculates the damped dispersion potential energy
c     using a pairwise double loop
c
c
      subroutine edisp0a
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use dsppot
      use energi
      use group
      use mutant
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,fgrp
      real*8 ci,ck
      real*8 r,r2,r3
      real*8 r4,r5,r6
      real*8 ai,ai2
      real*8 ak,ak2
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 damp3,damp5
      real*8 damp,taper
      real*8 dlsc,dlsq
      real*8 soft,dsoft,eps
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
      logical mutik
      character*6 mode
c
c
c     zero out the total damped dispersion energy
c
      edsp = 0.0d0
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
c     set soft core lambda scaling values for mutated interactions
c
      if (nmut .ne. 0)  call dspsc (dlsc,dlsq)
      rs2 = scrdsp * scrdsp
c
c     set damping tolerance, cutoff and switching coefficients
c
      eps = 0.0001d0
      mode = 'DISP'
      call switch (mode)
c
c     find the damped dispersion energy via double loop search
c
      do ii = 1, ndisp-1
         i = idisp(ii)
         ci = csix(i)
         ai = adisp(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         usei = use(i)
         muti = (mutg(i) .ne. 0)
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
            ck = csix(k)
            ak = adisp(k)
            mutk = (mutg(k) .ne. 0)
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
c
c     find the damping factor for the dispersion interaction
c
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (abs(ai-ak) .lt. eps) then
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                  else
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
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  e = e * dspscale(k) * damp**2
                  if (use_group)  e = e * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (muti.or.mutk) .and.
     &                    (vcouple.eq.1 .or. .not.(muti.and.mutk))
                  if (mutik) then
                     soft = dlsc
                     if (r2 .lt. rs2) then
                        call edspsc (r,dlsc,dlsq,soft,dsoft)
                     end if
                     e = e * soft
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
                     e = e * taper
                  end if
c
c     increment the overall damped dispersion energy component
c
                  edsp = edsp + e
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
         ci = csix(i)
         ai = adisp(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         usei = use(i)
         muti = (mutg(i) .ne. 0)
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
            ck = csix(k)
            ak = adisp(k)
            mutk = (mutg(k) .ne. 0)
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
c
c     find the damping factor for the dispersion interaction
c
                     di = ai * r
                     di2 = di * di
                     di3 = di * di2
                     dk = ak * r
                     expi = exp(-di)
                     expk = exp(-dk)
                     if (abs(ai-ak) .lt. eps) then
                        di4 = di2 * di2
                        di5 = di2 * di3
                        damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                             +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                        damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                           +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                     else
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
     &                           - 2.0d0*ti2*tk*(1.0+di+di2/3.0d0)*expi
     &                           - 2.0d0*tk2*ti*(1.0+dk+dk2/3.0d0)*expk
                     end if
                     damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                     e = e * damp**2
                     if (use_polymer) then
                        if (r2 .le. polycut2)  e = e * dspscale(k)
                     end if
                     if (use_group)  e = e * fgrp
                     if (i .eq. k)  e = 0.5d0 * e
c
c     set use of lambda scaling for decoupling or annihilation
c
                     mutik = (muti.or.mutk) .and.
     &                       (vcouple.eq.1 .or. .not.(muti.and.mutk))
                     if (mutik) then
                        soft = dlsc
                        if (r2 .lt. rs2) then
                           call edspsc (r,dlsc,dlsq,soft,dsoft)
                        end if
                        e = e * soft
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
                        e = e * taper
                     end if
c
c     increment the overall damped dispersion energy component
c
                     edsp = edsp + e
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine edisp0b  --  damp dispersion via neighbor list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "edisp0b" calculates the damped dispersion potential energy
c     using a pairwise neighbor list
c
c
      subroutine edisp0b
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use dsppot
      use energi
      use group
      use mutant
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,fgrp
      real*8 ci,ck
      real*8 r,r2,r3
      real*8 r4,r5,r6
      real*8 ai,ai2
      real*8 ak,ak2
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 damp3,damp5
      real*8 damp,taper
      real*8 dlsc,dlsq
      real*8 soft,dsoft,eps
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
      logical mutik
      character*6 mode
c
c
c     zero out the total damped dispersion energy
c
      edsp = 0.0d0
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
c     set soft core lambda scaling values for mutated interactions
c
      if (nmut .ne. 0)  call dspsc (dlsc,dlsq)
      rs2 = scrdsp * scrdsp
c
c     set damping tolerance, cutoff and switching coefficients
c
      eps = 0.0001d0
      mode = 'DISP'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ndisp,idisp,csix,adisp,use,
!$OMP& x,y,z,n12,n13,n14,n15,i12,i13,i14,i15,nvlst,vlst,use_group,
!$OMP& dsp2scale,dsp3scale,dsp4scale,dsp5scale,mutg,off2,cut2,
!$OMP& c0,c1,c2,c3,c4,c5,vcouple,dlsc,dlsq,rs2,eps)
!$OMP& firstprivate(dspscale) shared(edsp)
!$OMP DO reduction(+:edsp)
c
c     find the damped dispersion energy via neighbor list search
c
      do ii = 1, ndisp
         i = idisp(ii)
         ci = csix(i)
         ai = adisp(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         usei = use(i)
         muti = (mutg(i) .ne. 0)
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
         do kk = 1, nvlst(i)
            k = vlst(kk,i)
            ck = csix(k)
            ak = adisp(k)
            mutk = (mutg(k) .ne. 0)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
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
c
c     find the damping factor for the dispersion interaction
c
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (abs(ai-ak) .lt. eps) then
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                  else
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
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  e = e * dspscale(k) * damp**2
                  if (use_group)  e = e * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (muti.or.mutk) .and.
     &                    (vcouple.eq.1 .or. .not.(muti.and.mutk))
                  if (mutik) then
                     soft = dlsc
                     if (r2 .lt. rs2) then
                        call edspsc (r,dlsc,dlsq,soft,dsoft)
                     end if
                     e = e * soft
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
                     e = e * taper
                  end if
c
c     increment the overall dispersion energy component
c
                  edsp = edsp + e
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
c     ##  subroutine edisp0c  --  Ewald dispersion energy via loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edisp0c" calculates the dispersion interaction energy using
c     particle mesh Ewald summation and a double loop
c
c
      subroutine edisp0c
      use disp
      use energi
      use ewald
      use pme
      implicit none
c
c
c     zero out the total damped dispersion energy
c
      edsp = 0.0d0
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
      call edreal0c
c
c     compute the reciprocal space part of the Ewald summation
c
      call edrecip
c
c     compute the self-energy portion of the Ewald summation
c
      call dspself (.false.,.false.)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine edreal0c  --  real space dispersion via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "edreal0c" calculates the damped dispersion potential energy
c     using a particle mesh Ewald sum and pairwise double loop
c
c
      subroutine edreal0c
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use dsppot
      use energi
      use ewald
      use group
      use mutant
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,fgrp
      real*8 ci,ck
      real*8 r,r2,r6
      real*8 ai,ai2
      real*8 ak,ak2
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 ralpha2
      real*8 expa,term
      real*8 damp3,damp5
      real*8 damp,scale
      real*8 dlsc,dlsq,wsc
      real*8 soft,dsoft,eps
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
      logical mutik
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
c     set soft core lambda scaling values for mutated interactions
c
      if (nmut .ne. 0)  call dspsc (dlsc,dlsq)
      rs2 = scrdsp * scrdsp
c
c     set damping tolerance, cutoff and switching coefficients
c
      eps = 0.0001d0
      mode = 'DEWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, ndisp-1
         i = idisp(ii)
         ci = csix(i)
         ai = adisp(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         usei = use(i)
         muti = (mutg(i) .ne. 0)
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
            ck = csix(k)
            ak = adisp(k)
            mutk = (mutg(k) .ne. 0)
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
                  ralpha2 = r2 * aewald**2
                  term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                  expa = exp(-ralpha2) * term
c
c     find the damping factor for the dispersion interaction
c
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (abs(ai-ak) .lt. eps) then
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                  else
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
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  scale = dspscale(k) * damp**2
                  if (use_group)  scale = scale * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (muti.or.mutk) .and.
     &                    (vcouple.eq.1 .or. .not.(muti.and.mutk))
                  wsc = 1.0d0
                  if (mutik) then
                     soft = dlsc
                     if (r2 .lt. rs2) then
                        call edspsc (r,dlsc,dlsq,soft,dsoft)
                     end if
                     wsc = dlsc
                     scale = scale * soft
                  end if
c
c     increment the overall dispersion energy component;
c     the lambda weight also scales the reciprocal space share
c
                  scale = scale - wsc
                  e = e * (wsc*expa+scale)
                  edsp = edsp + e
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
         ci = csix(i)
         ai = adisp(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         usei = use(i)
         muti = (mutg(i) .ne. 0)
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
            ck = csix(k)
            ak = adisp(k)
            mutk = (mutg(k) .ne. 0)
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
                     ralpha2 = r2 * aewald**2
                     term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                     expa = exp(-ralpha2) * term
c
c     find the damping factor for the dispersion interaction
c
                     di = ai * r
                     di2 = di * di
                     di3 = di * di2
                     dk = ak * r
                     expi = exp(-di)
                     expk = exp(-dk)
                     if (abs(ai-ak) .lt. eps) then
                        di4 = di2 * di2
                        di5 = di2 * di3
                        damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                             +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                        damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                           +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                     else
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
                     end if
                     damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                     scale = damp**2
                     if (use_group)  scale = scale * fgrp
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           scale = scale * dspscale(k)
                        end if
                     end if
c
c     set use of lambda scaling for decoupling or annihilation
c
                     mutik = (muti.or.mutk) .and.
     &                       (vcouple.eq.1 .or. .not.(muti.and.mutk))
                     wsc = 1.0d0
                     if (mutik) then
                        soft = dlsc
                        if (r2 .lt. rs2) then
                           call edspsc (r,dlsc,dlsq,soft,dsoft)
                        end if
                        wsc = dlsc
                        scale = scale * soft
                     end if
c
c     increment the overall dispersion energy component;
c     the lambda weight also scales the reciprocal space share
c
                     scale = scale - wsc
                     e = e * (wsc*expa+scale)
                     if (i .eq. k)  e = 0.5d0 * e
                     edsp = edsp + e
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine edisp0d  --  Ewald dispersion energy via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "edisp0d" calculates the dispersion interaction energy using
c     particle mesh Ewald summation and a neighbor list
c
c
      subroutine edisp0d
      use disp
      use energi
      use ewald
      use pme
      implicit none
c
c
c     zero out the total damped dispersion energy
c
      edsp = 0.0d0
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
      call edreal0d
c
c     compute the reciprocal space part of the Ewald summation
c
      call edrecip
c
c     compute the self-energy portion of the Ewald summation
c
      call dspself (.false.,.false.)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine edreal0d  --  real space dispersion via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "edreal0d" calculates the real space portion of the damped
c     dispersion energy using a neighbor list
c
c
      subroutine edreal0d
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use dsppot
      use energi
      use ewald
      use group
      use mutant
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,fgrp
      real*8 ci,ck
      real*8 r,r2,r6
      real*8 ai,ai2
      real*8 ak,ak2
      real*8 di,di2,di3
      real*8 di4,di5
      real*8 dk,dk2,dk3
      real*8 ti,ti2
      real*8 tk,tk2
      real*8 expi,expk
      real*8 ralpha2
      real*8 expa,term
      real*8 damp3,damp5
      real*8 damp,scale
      real*8 dlsc,dlsq,wsc
      real*8 soft,dsoft,eps
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
      logical mutik
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
c     set soft core lambda scaling values for mutated interactions
c
      if (nmut .ne. 0)  call dspsc (dlsc,dlsq)
      rs2 = scrdsp * scrdsp
c
c     set damping tolerance, cutoff and switching coefficients
c
      eps = 0.0001d0
      mode = 'DEWALD'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ndisp,idisp,csix,adisp,use,
!$OMP& x,y,z,n12,n13,n14,n15,i12,i13,i14,i15,nvlst,vlst,use_group,
!$OMP& dsp2scale,dsp3scale,dsp4scale,dsp5scale,mutg,off2,aewald,
!$OMP& vcouple,dlsc,dlsq,rs2,eps)
!$OMP& firstprivate(dspscale) shared(edsp)
!$OMP DO reduction(+:edsp)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, ndisp
         i = idisp(ii)
         ci = csix(i)
         ai = adisp(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         usei = use(i)
         muti = (mutg(i) .ne. 0)
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
         do kk = 1, nvlst(i)
            k = vlst(kk,i)
            ck = csix(k)
            ak = adisp(k)
            mutk = (mutg(k) .ne. 0)
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
                  ralpha2 = r2 * aewald**2
                  term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                  expa = exp(-ralpha2) * term
c
c     find the damping factor for the dispersion interaction
c
                  di = ai * r
                  di2 = di * di
                  di3 = di * di2
                  dk = ak * r
                  expi = exp(-di)
                  expk = exp(-dk)
                  if (abs(ai-ak) .lt. eps) then
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
                  else
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
                  end if
                  damp = 1.5d0*damp5 - 0.5d0*damp3
c
c     apply damping and scaling factors for this interaction
c
                  scale = dspscale(k) * damp**2
                  if (use_group)  scale = scale * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (muti.or.mutk) .and.
     &                    (vcouple.eq.1 .or. .not.(muti.and.mutk))
                  wsc = 1.0d0
                  if (mutik) then
                     soft = dlsc
                     if (r2 .lt. rs2) then
                        call edspsc (r,dlsc,dlsq,soft,dsoft)
                     end if
                     wsc = dlsc
                     scale = scale * soft
                  end if
c
c     increment the overall dispersion energy component;
c     the lambda weight also scales the reciprocal space share
c
                  scale = scale - wsc
                  e = e * (wsc*expa+scale)
                  edsp = edsp + e
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine edrecip  --  PME recip space damped dispersion  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "edrecip" evaluates the reciprocal space portion of the particle
c     mesh Ewald energy due to damped dispersion
c
c     for a mutated system the environment and ligand sites fill the
c     real and imaginary parts of one grid, and their transforms are
c     separated by symmetry to weight each pair type by lambda
c
c
      subroutine edrecip
      use boxes
      use bound
      use disp
      use energi
      use ewald
      use math
      use pme
      implicit none
      integer i,j
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      real*8 e,denom
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 term,expterm
      real*8 eterm,denom0
      real*8 hsq,struc2
      real*8 h,hhh,b,bfac
      real*8 fac1,fac2,fac3
      real*8 erfcterm
      real*8 wel,wll,cpr
      real*8 ser,sei,slr,sli
      logical split
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
         if (size(qgrid) .ne. 2*ntot)  call fftclose
      end if
      if (.not. allocated(qgrid))  call fftsetup
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
      call dspmutwt (split,wel,wll,cpr)
      call grid_disp (split)
      call fftfront (qgrid)
c
c     use scalar sum to get the reciprocal space energy
c
      bfac = pi / aewald
      fac1 = 2.0d0 * pi**(3.5d0)
      fac2 = aewald**3
      fac3 = -2.0d0 * aewald * pi**2
      denom0 = (6.0d0*volbox) / (pi**1.5d0)
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      nff = nfft1 * nfft2
      ntot = nff * nfft3
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         h = sqrt(hsq)
         b = h * bfac
         hhh = h * hsq
         term = -b * b
         if (term .gt. -50.0d0) then
            denom = denom0 * bsmod1(k1) * bsmod2(k2) * bsmod3(k3)
            expterm = exp(term)
            erfcterm = erfc(b)
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               erfcterm = erfcterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (nonprism) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               if (mod(m1+m2+m3,2) .ne. 0)  erfcterm = 0.0d0
            end if
            eterm = (-fac1*erfcterm*hhh-expterm*(fac2+fac3*hsq))/denom
            if (split) then
               call dspsplit (qgrid,k1,k2,k3,ser,sei,slr,sli)
               struc2 = ser*ser + sei*sei + wll*(slr*slr+sli*sli)
     &                     + 2.0d0*wel*(ser*slr+sei*sli)
            else
               struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
            end if
            e = eterm * struc2
            edsp = edsp + e
         end if
      end do
c
c     account for the total energy correction term
c
      e = -cpr * aewald**3 / denom0
      edsp = edsp + e
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dspmutwt  --  mutated dispersion PME weights  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dspmutwt" finds the lambda weights of the ligand to environment
c     and the ligand to ligand dispersion pairs, whether the PME sum
c     has to carry the environment and ligand separately, and the
c     weighted pairwise sum of C6 values for the zero wavevector term
c
c     the lambda power of the soft core scales every pair with a ligand
c     atom under annihilation, but only the ligand to environment pairs
c     under decoupling, so the weights are not a product of site scales
c     and a single real density cannot reproduce them; the soft core
c     switch acts only inside its radius, so reciprocal space carries
c     the plain lambda power
c
c
      subroutine dspmutwt (split,wel,wll,cpr)
      use disp
      implicit none
      real*8 wel,wll,cpr
      real*8 dwel,dwll,dcpr
      real*8 d2wel,d2wll,d2cpr
      logical split
c
c
c     get the weights from the routine that also finds their lambda
c     derivatives, which are not needed here
c
      call dspmutwt4 (split,wel,wll,cpr,dwel,dwll,dcpr,
     &                d2wel,d2wll,d2cpr)
c
c     a fully coupled system needs no split
c
      if (wel.eq.1.0d0 .and. wll.eq.1.0d0) then
         split = .false.
         cpr = csixpr
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine dspmutwt4  --  dispersion PME weight derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "dspmutwt4" finds the lambda weights of the mutated dispersion
c     PME sum as in "dspmutwt", along with their first and second
c     derivatives with respect to lambda and those of the weighted
c     C6 pair sum, and always carries the environment and ligand
c     separately, since the weight derivatives are not zero at the
c     fully coupled end
c
c
      subroutine dspmutwt4 (split,wel,wll,cpr,dwel,dwll,dcpr,
     &                      d2wel,d2wll,d2cpr)
      use disp
      use dlmda
      use mutant
      implicit none
      integer i,ii
      real*8 v,v1,v2
      real*8 wel,wll,cpr
      real*8 dwel,dwll,dcpr
      real*8 d2wel,d2wll,d2cpr
      real*8 sume,suml
      logical split
c
c
c     an unmutated system needs no split
c
      split = .false.
      wel = 1.0d0
      wll = 1.0d0
      cpr = csixpr
      dwel = 0.0d0
      dwll = 0.0d0
      dcpr = 0.0d0
      d2wel = 0.0d0
      d2wll = 0.0d0
      d2cpr = 0.0d0
      if (nmut .eq. 0)  return
c
c     set the weights and their derivatives for annihilation or for
c     decoupling from the lambda power of the soft core
c
      call lmdapwr (vlambda,scexpd,use_d2lmda,v,v1,v2)
      split = .true.
      wel = v
      dwel = v1
      d2wel = v2
      if (vcouple .eq. 1) then
         wll = v
         dwll = v1
         d2wll = v2
      end if
c
c     weight the environment and ligand parts of the C6 pair sum
c
      sume = 0.0d0
      suml = 0.0d0
      do ii = 1, ndisp
         i = idisp(ii)
         if (mutg(i) .eq. 0) then
            sume = sume + csix(i)
         else
            suml = suml + csix(i)
         end if
      end do
      cpr = sume*sume + 2.0d0*wel*sume*suml + wll*suml*suml
      dcpr = 2.0d0*dwel*sume*suml + dwll*suml*suml
      d2cpr = 2.0d0*d2wel*sume*suml + d2wll*suml*suml
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine dspself  --  Ewald dispersion self energy  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "dspself" adds the self-energy portion of the Ewald summation
c     for damped dispersion, weighting each ligand site as the
c     reciprocal space sum does, along with its lambda derivatives
c     when "deriv" is set and its partitioning over atoms when "part"
c     is set
c
c
      subroutine dspself (deriv,part)
      use action
      use analyz
      use disp
      use dlmda
      use energi
      use ewald
      use mutant
      implicit none
      integer i,ii
      real*8 e,wel,wll,cpr
      real*8 dwel,dwll,dcpr
      real*8 d2wel,d2wll,d2cpr
      logical deriv,part
      logical split
      logical dodl,dod2
c
c
c     get the ligand weights, along with their lambda derivatives
c     while the van der Waals lambda map has slope or curvature
c
      dodl = .false.
      dod2 = .false.
      if (deriv)  call vdlflags (dodl,dod2)
      if (dodl) then
         call dspmutwt4 (split,wel,wll,cpr,dwel,dwll,dcpr,
     &                   d2wel,d2wll,d2cpr)
      else
         call dspmutwt (split,wel,wll,cpr)
      end if
c
c     add the weighted self energy of each dispersion site
c
      do ii = 1, ndisp
         i = idisp(ii)
         e = csix(i)**2 * aewald**6 / 12.0d0
         if (mutg(i) .ne. 0) then
            if (dodl)  devdl = devdl + dwll*e
            if (dod2)  d2evdl2 = d2evdl2 + d2wll*e
            e = e * wll
         end if
         edsp = edsp + e
         if (part) then
            nedsp = nedsp + 1
            aedsp(i) = aedsp(i) + e
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dspsplit  --  separate packed PME transforms  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dspsplit" recovers the transforms of the environment and the
c     ligand dispersion densities at one wavevector from a grid that
c     held them as its real and imaginary parts, using the conjugate
c     symmetry of the transform of a real density
c
c
      subroutine dspsplit (qg,k1,k2,k3,ser,sei,slr,sli)
      use pme
      implicit none
      integer k1,k2,k3
      integer j1,j2,j3
      real*8 zr,zi,yr,yi
      real*8 ser,sei,slr,sli
      real*8 qg(2,nfft1,nfft2,nfft3)
c
c
c     find the grid point of the opposite wavevector
c
      j1 = mod(nfft1+1-k1,nfft1) + 1
      j2 = mod(nfft2+1-k2,nfft2) + 1
      j3 = mod(nfft3+1-k3,nfft3) + 1
c
c     split the packed values into the two separate transforms
c
      zr = qg(1,k1,k2,k3)
      zi = qg(2,k1,k2,k3)
      yr = qg(1,j1,j2,j3)
      yi = qg(2,j1,j2,j3)
      ser = 0.5d0 * (zr+yr)
      sei = 0.5d0 * (zi-yi)
      slr = 0.5d0 * (zi+yi)
      sli = 0.5d0 * (yr-zr)
      return
      end

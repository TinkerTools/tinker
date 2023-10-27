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
      real*8 vterm
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
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
c     set lambda scaling values for mutated interactions
c
      if (nmut .ne. 0) then
         vterm = vlambda**4 / sqrt(1.0d0+vlambda**2-vlambda**3)
      end if
c
c     set conversion factor, cutoff and switching coefficients
c
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
         muti = mut(i)
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
            mutk = mut(k)
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
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
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
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        e = e * vterm
                     else if (.not.muti .or. .not.mutk) then
                        e = e * vterm
                     end if
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
         muti = mut(i)
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
            mutk = mut(k)
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
     &                           - 2.0d0*ti2*tk*(1.0+di+di2/3.0d0)*expi
     &                           - 2.0d0*tk2*ti*(1.0+dk+dk2/3.0d0)*expk
                     else
                        di4 = di2 * di2
                        di5 = di2 * di3
                        damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                             +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                        damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                           +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
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
                     if (muti .or. mutk) then
                        if (vcouple .eq. 1) then
                           e = e * vterm
                        else if (.not.muti .or. .not.mutk) then
                           e = e * vterm
                        end if
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
      real*8 vterm
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
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
c     set lambda scaling values for mutated interactions
c
      if (nmut .ne. 0) then
         vterm = vlambda**4 / sqrt(1.0d0+vlambda**2-vlambda**3)
      end if
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
!$OMP& dsp2scale,dsp3scale,dsp4scale,dsp5scale,mut,off2,cut2,
!$OMP& c0,c1,c2,c3,c4,c5,vcouple,vterm)
!$OMP& firstprivate(dspscale) shared(edsp)
!$OMP DO reduction(+:edsp) schedule(guided)
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
         muti = mut(i)
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
            mutk = mut(k)
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
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
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
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        e = e * vterm
                     else if (.not.muti .or. .not.mutk) then
                        e = e * vterm
                     end if
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
      integer i,ii
      real*8 e
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
      do ii = 1, ndisp
         i = idisp(i)
         e = csix(i)**2 * aewald**6 / 12.0d0
         edsp = edsp + e
      end do
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
      real*8 vterm
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
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
c     set lambda scaling values for mutated interactions
c
      if (nmut .ne. 0) then
         vterm = vlambda**4 / sqrt(1.0d0+vlambda**2-vlambda**3)
      end if
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
         ci = csix(i)
         ai = adisp(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         usei = use(i)
         muti = mut(i)
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
            mutk = mut(k)
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
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
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
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        scale = scale * vterm
                     else if (.not.muti .or. .not.mutk) then
                        scale = scale * vterm
                     end if
                  end if
c
c     increment the overall dispersion energy component
c
                  scale = scale - 1.0d0
                  e = e * (expa+scale)
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
         muti = mut(i)
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
            mutk = mut(k)
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
                     else
                        di4 = di2 * di2
                        di5 = di2 * di3
                        damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                             +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                        damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                           +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
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
                     if (muti .or. mutk) then
                        if (vcouple .eq. 1) then
                           scale = scale * vterm
                        else if (.not.muti .or. .not.mutk) then
                           scale = scale * vterm
                        end if
                     end if
c
c     increment the overall dispersion energy component
c
                     scale = scale - 1.0d0
                     e = e * (expa+scale)
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
      integer i,ii
      real*8 e
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
      do ii = 1, ndisp
         i = idisp(ii)
         e = csix(i)**2 * aewald**6 / 12.0d0
         edsp = edsp + e
      end do
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
c     "edreal0d" evaluated the real space portion of the damped
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
      real*8 vterm
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
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
c     set lambda scaling values for mutated interactions
c
      if (nmut .ne. 0) then
         vterm = vlambda**4 / sqrt(1.0d0+vlambda**2-vlambda**3)
      end if
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
!$OMP& dsp2scale,dsp3scale,dsp4scale,dsp5scale,mut,off2,aewald,
!$OMP& vcouple,vterm)
!$OMP& firstprivate(dspscale) shared(edsp)
!$OMP DO reduction(+:edsp) schedule(guided)
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
         muti = mut(i)
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
            mutk = mut(k)
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
                  else
                     di4 = di2 * di2
                     di5 = di2 * di3
                     damp3 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +7.0d0*di3/48.0d0+di4/48.0d0)*expi
                     damp5 = 1.0d0 - (1.0d0+di+0.5d0*di2
     &                          +di3/6.0d0+di4/24.0d0+di5/144.0d0)*expi
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
                  if (muti .or. mutk) then
                     if (vcouple .eq. 1) then
                        scale = scale * vterm
                     else if (.not.muti .or. .not.mutk) then
                        scale = scale * vterm
                     end if
                  end if
c
c     increment the overall dispersion energy component
c
                  scale = scale - 1.0d0
                  e = e * (expa+scale)
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
      call grid_disp
      call fftfront
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
            struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
            e = eterm * struc2
            edsp = edsp + e
         end if
      end do
c
c     account for the total energy correction term
c
      e = -csixpr * aewald**3 / denom0
      edsp = edsp + e
      return
      end

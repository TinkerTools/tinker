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
      integer mutgi,ig
      integer ii,kk
      integer jcell
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,fgrp
      real*8 ci,ck
      real*8 r,r2,r3
      real*8 r4,r5,r6
      real*8 ai,ak
      real*8 damp,ddamp,taper
      real*8 dlsc(2),dlsq(2)
      real*8 soft,dsoft
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
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
c     set cutoff and switching coefficients
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
         mutgi = mutg(i)
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
            ig = mutcls(mutg(k),mutgi)
c
c     the two ligands of a relative dual topology never interact
c
            proceed = (ig .ge. 0)
            if (use_group .and. proceed)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
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
                  call dampdisp (r,ai,ak,damp,ddamp)
c
c     apply damping and scaling factors for this interaction
c
                  e = e * dspscale(k) * damp**2
                  if (use_group)  e = e * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (ig .gt. 0)
                  if (mutik) then
                     soft = dlsc(ig)
                     if (r2 .lt. rs2) then
                        call edspsc (r,ig,dlsc,dlsq,soft,dsoft)
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
         mutgi = mutg(i)
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
            ig = mutcls(mutg(k),mutgi)
c
c     the two ligands of a relative dual topology never interact
c
            proceed = (ig .ge. 0)
            if (use_group .and. proceed)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
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
                     call dampdisp (r,ai,ak,damp,ddamp)
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
                     mutik = (ig .gt. 0)
                     if (mutik) then
                        soft = dlsc(ig)
                        if (r2 .lt. rs2) then
                           call edspsc (r,ig,dlsc,dlsq,soft,dsoft)
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
      integer mutgi,ig
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,fgrp
      real*8 ci,ck
      real*8 r,r2,r3
      real*8 r4,r5,r6
      real*8 ai,ak
      real*8 damp,ddamp,taper
      real*8 dlsc(2),dlsq(2)
      real*8 soft,dsoft
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
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
c     set cutoff and switching coefficients
c
      mode = 'DISP'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ndisp,idisp,csix,adisp,use,
!$OMP& x,y,z,n12,n13,n14,n15,i12,i13,i14,i15,ndlst,dlst,use_group,
!$OMP& dsp2scale,dsp3scale,dsp4scale,dsp5scale,mutg,off2,cut2,
!$OMP& c0,c1,c2,c3,c4,c5,mutcls,dlsc,dlsq,rs2)
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
         mutgi = mutg(i)
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
         do kk = 1, ndlst(i)
            k = dlst(kk,i)
            ck = csix(k)
            ak = adisp(k)
            ig = mutcls(mutg(k),mutgi)
c
c     the two ligands of a relative dual topology never interact
c
            proceed = (ig .ge. 0)
            if (use_group .and. proceed)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
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
                  call dampdisp (r,ai,ak,damp,ddamp)
c
c     apply damping and scaling factors for this interaction
c
                  e = e * dspscale(k) * damp**2
                  if (use_group)  e = e * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (ig .gt. 0)
                  if (mutik) then
                     soft = dlsc(ig)
                     if (r2 .lt. rs2) then
                        call edspsc (r,ig,dlsc,dlsq,soft,dsoft)
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
      integer mutgi,ig
      integer ii,kk
      integer jcell
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,fgrp
      real*8 ci,ck
      real*8 r,r2,r6
      real*8 ai,ak
      real*8 ralpha2
      real*8 expa,term
      real*8 damp,ddamp,scale
      real*8 dmp2,dhw,dhdw
      real*8 dlsc(2),dlsq(2),wsc
      real*8 soft,dsoft
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
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
c     set cutoff and switching coefficients
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
         mutgi = mutg(i)
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
            ig = mutcls(mutg(k),mutgi)
c
c     the two ligands of a relative dual topology never interact
c
            proceed = (ig .ge. 0)
            if (use_group .and. proceed)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
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
c     find the damping factor, which an undamped hal tail skips
c
                  damp = 1.0d0
                  if (dhdamp .or. .not.use_dhal) then
                     call dampdisp (r,ai,ak,damp,ddamp)
                  end if
                  dmp2 = damp * damp
c
c     an Ewald tail mapped from hal takes over across the vdw switch
c
                  if (use_dhal) then
                     call dhalwt (r,dhw,dhdw)
                     dmp2 = dhw * dmp2
                  end if
c
c     apply damping and scaling factors for this interaction
c
                  scale = dspscale(k) * dmp2
                  if (use_group)  scale = scale * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (ig .gt. 0)
                  wsc = 1.0d0
                  if (mutik) then
                     soft = dlsc(ig)
                     if (r2 .lt. rs2) then
                        call edspsc (r,ig,dlsc,dlsq,soft,dsoft)
                     end if
                     wsc = dlsc(ig)
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
         mutgi = mutg(i)
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
            ig = mutcls(mutg(k),mutgi)
c
c     the two ligands of a relative dual topology never interact
c
            proceed = (ig .ge. 0)
            if (use_group .and. proceed)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
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
c     find the damping factor, which an undamped hal tail skips
c
                     damp = 1.0d0
                     if (dhdamp .or. .not.use_dhal) then
                        call dampdisp (r,ai,ak,damp,ddamp)
                     end if
                     dmp2 = damp * damp
c
c     an Ewald tail mapped from hal takes over across the vdw switch
c
                     if (use_dhal) then
                        call dhalwt (r,dhw,dhdw)
                        dmp2 = dhw * dmp2
                     end if
c
c     apply damping and scaling factors for this interaction
c
                     scale = dmp2
                     if (use_group)  scale = scale * fgrp
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           scale = scale * dspscale(k)
                        end if
                     end if
c
c     set use of lambda scaling for decoupling or annihilation
c
                     mutik = (ig .gt. 0)
                     wsc = 1.0d0
                     if (mutik) then
                        soft = dlsc(ig)
                        if (r2 .lt. rs2) then
                           call edspsc (r,ig,dlsc,dlsq,soft,dsoft)
                        end if
                        wsc = dlsc(ig)
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
      integer mutgi,ig
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,fgrp
      real*8 ci,ck
      real*8 r,r2,r6
      real*8 ai,ak
      real*8 ralpha2
      real*8 expa,term
      real*8 damp,ddamp,scale
      real*8 dmp2,dhw,dhdw
      real*8 dlsc(2),dlsq(2),wsc
      real*8 soft,dsoft
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
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
c     set cutoff and switching coefficients
c
      mode = 'DEWALD'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ndisp,idisp,csix,adisp,use,
!$OMP& x,y,z,n12,n13,n14,n15,i12,i13,i14,i15,ndlst,dlst,use_group,
!$OMP& dsp2scale,dsp3scale,dsp4scale,dsp5scale,mutg,off2,aewald,
!$OMP& mutcls,dlsc,dlsq,rs2)
!$OMP& shared(use_dhal,dhdamp)
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
         mutgi = mutg(i)
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
         do kk = 1, ndlst(i)
            k = dlst(kk,i)
            ck = csix(k)
            ak = adisp(k)
            ig = mutcls(mutg(k),mutgi)
c
c     the two ligands of a relative dual topology never interact
c
            proceed = (ig .ge. 0)
            if (use_group .and. proceed)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
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
c     find the damping factor, which an undamped hal tail skips
c
                  damp = 1.0d0
                  if (dhdamp .or. .not.use_dhal) then
                     call dampdisp (r,ai,ak,damp,ddamp)
                  end if
                  dmp2 = damp * damp
c
c     an Ewald tail mapped from hal takes over across the vdw switch
c
                  if (use_dhal) then
                     call dhalwt (r,dhw,dhdw)
                     dmp2 = dhw * dmp2
                  end if
c
c     apply damping and scaling factors for this interaction
c
                  scale = dspscale(k) * dmp2
                  if (use_group)  scale = scale * fgrp
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (ig .gt. 0)
                  wsc = 1.0d0
                  if (mutik) then
                     soft = dlsc(ig)
                     if (r2 .lt. rs2) then
                        call edspsc (r,ig,dlsc,dlsq,soft,dsoft)
                     end if
                     wsc = dlsc(ig)
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
c     for a mutated system the environment and ligand groups fill
c     the density slots set by "dspslot", the first two as the real
c     and imaginary parts of one grid and a third on a second grid,
c     and their transforms are separated by symmetry to weight each
c     pair type by lambda
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
      integer j1,j2,j3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer nslot
      integer gslot(0:2)
      integer gpart(0:2)
      real*8 e,denom
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 term,expterm
      real*8 eterm,denom0
      real*8 hsq,struc2
      real*8 h,hhh,b,bfac
      real*8 fac1,fac2,fac3
      real*8 erfcterm
      real*8 cpr,dcpr,d2cpr
      real*8 dspstruc
      real*8 sr(3),si(3)
      real*8 qa(2,2),qb(2)
      real*8 ws(3,3)
      real*8 dws(3,3)
      real*8 d2ws(3,3)
      real*8, allocatable :: qgb(:,:,:,:)
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
c     assign the environment and ligand groups to density slots
c
      call dspslot (.false.,.false.,nslot,gslot,ws,dws,d2ws,
     &              cpr,dcpr,d2cpr)
c
c     no group left on the grid leaves no energy, as the pairwise
c     sum of C6 values is then zero as well
c
      if (nslot .eq. 0)  return
c
c     assign PME grids and perform 3-D FFT forward transforms, where
c     a third density slot takes a second grid
c
      do i = 0, 2
         gpart(i) = gslot(i)
         if (gslot(i) .eq. 3)  gpart(i) = 0
      end do
      call grid_disp (qgrid,gpart)
      call fftfront (qgrid)
      if (nslot .eq. 3) then
         allocate (qgb(2,nfft1,nfft2,nfft3))
         do i = 0, 2
            gpart(i) = 0
            if (gslot(i) .eq. 3)  gpart(i) = 1
         end do
         call grid_disp (qgb,gpart)
         call fftfront (qgb)
      end if
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
            call pmeopp (k1,k2,k3,j1,j2,j3)
            do j = 1, 2
               qa(j,1) = qgrid(j,k1,k2,k3)
               qa(j,2) = qgrid(j,j1,j2,j3)
               if (nslot .eq. 3)  qb(j) = qgb(j,k1,k2,k3)
            end do
            call dspsf (nslot,qa(1,1),qa(1,2),qb,sr,si)
            struc2 = dspstruc (nslot,ws,sr,si)
            e = eterm * struc2
            edsp = edsp + e
         end if
      end do
c
c     account for the total energy correction term
c
      e = -cpr * aewald**3 / denom0
      edsp = edsp + e
      if (nslot .eq. 3)  deallocate (qgb)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine dspslot  --  dispersion PME density slots  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "dspslot" assigns the dispersion sites of each group to the
c     real densities "slots" carried by the PME transforms, and finds
c     the pair weights between slots, their lambda derivatives when
c     "dodl" or "dod2" is set, and the weighted pairwise sum of C6
c     values for the zero wavevector term
c
c     the lambda power of the soft core scales every pair with a ligand
c     atom under annihilation, but only the ligand to environment pairs
c     under decoupling, so the weights are not a product of site scales
c     and a single real density cannot reproduce them; the soft core
c     switch acts only inside its radius, so reciprocal space carries
c     the plain lambda power
c
c     a group with no sites, or with no weight left to carry, is not
c     placed on the grid, and a group whose weights match those of an
c     earlier group shares its slot, so a fully coupled system and the
c     ligand legs of a relative schedule need a single transform; the
c     first two slots are the real and imaginary parts of the main PME
c     grid, and a third slot needs a second grid
c
c
      subroutine dspslot (dodl,dod2,nslot,gslot,ws,dws,d2ws,
     &                    cpr,dcpr,d2cpr)
      use disp
      use mutant
      implicit none
      integer i,ii
      integer g,h,f
      integer nslot
      integer gslot(0:2)
      integer grp(3)
      real*8 cpr,dcpr,d2cpr
      real*8 csum(0:2)
      real*8 wg(0:2,0:2)
      real*8 dwg(0:2,0:2)
      real*8 d2wg(0:2,0:2)
      real*8 ws(3,3)
      real*8 dws(3,3)
      real*8 d2ws(3,3)
      logical dodl,dod2
      logical same,live
      logical have(0:2)
c
c
c     get the pair weights of the environment and ligand groups
c
      call grpvwt (scexpd,wg,dwg,d2wg)
c
c     find the groups with sites and their sums of C6 values
c
      do g = 0, 2
         csum(g) = 0.0d0
         have(g) = .false.
      end do
      do ii = 1, ndisp
         i = idisp(ii)
         g = mutg(i)
         csum(g) = csum(g) + csix(i)
         have(g) = .true.
      end do
c
c     drop a group with no weight to any group on the grid
c
      do g = 0, 2
         if (have(g)) then
            live = .false.
            do h = 0, 2
               if (have(h)) then
                  if (wg(g,h) .ne. 0.0d0)  live = .true.
                  if (dodl .and. dwg(g,h).ne.0.0d0)  live = .true.
                  if (dod2 .and. d2wg(g,h).ne.0.0d0)  live = .true.
               end if
            end do
            have(g) = live
         end if
      end do
c
c     merge a group into an earlier one with the same weights
c
      nslot = 0
      do g = 0, 2
         gslot(g) = 0
         if (have(g)) then
            do h = 0, g-1
               if (gslot(h).ne.0 .and. gslot(g).eq.0) then
                  same = .true.
                  do f = 0, 2
                     if (have(f)) then
                        if (wg(g,f) .ne. wg(h,f))  same = .false.
                        if (dodl .and. dwg(g,f).ne.dwg(h,f))
     &                     same = .false.
                        if (dod2 .and. d2wg(g,f).ne.d2wg(h,f))
     &                     same = .false.
                     end if
                  end do
                  if (same)  gslot(g) = gslot(h)
               end if
            end do
            if (gslot(g) .eq. 0) then
               nslot = nslot + 1
               gslot(g) = nslot
               grp(nslot) = g
            end if
         end if
      end do
c
c     set the pair weights between slots
c
      do i = 1, 3
         do ii = 1, 3
            ws(ii,i) = 0.0d0
            dws(ii,i) = 0.0d0
            d2ws(ii,i) = 0.0d0
         end do
      end do
      do i = 1, nslot
         do ii = 1, nslot
            ws(ii,i) = wg(grp(ii),grp(i))
            dws(ii,i) = dwg(grp(ii),grp(i))
            d2ws(ii,i) = d2wg(grp(ii),grp(i))
         end do
      end do
c
c     weight the pairwise sum of C6 values over the groups
c
      cpr = 0.0d0
      dcpr = 0.0d0
      d2cpr = 0.0d0
      do g = 0, 2
         do h = 0, 2
            if (have(g) .and. have(h)) then
               cpr = cpr + wg(g,h)*csum(g)*csum(h)
               dcpr = dcpr + dwg(g,h)*csum(g)*csum(h)
               d2cpr = d2cpr + d2wg(g,h)*csum(g)*csum(h)
            end if
         end do
      end do
c
c     an unmutated system keeps the stored pairwise sum
c
      if (nmut .eq. 0)  cpr = csixpr
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine dspsf  --  dispersion slot transforms at k  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "dspsf" gets the transforms "sr" and "si" of each dispersion
c     density slot at one wavevector, from the value "za" there and
c     "ya" at the opposite wavevector of the main transform, which
c     holds the first two slots as its real and imaginary parts, and
c     from the value "zb" there of the second transform, which holds
c     the third slot
c
c
      subroutine dspsf (nslot,za,ya,zb,sr,si)
      implicit none
      integer nslot
      real*8 za(2),ya(2),zb(2)
      real*8 sr(3),si(3)
c
c
c     a single slot is the transform itself, while two slots are
c     packed and must be separated
c
      if (nslot .eq. 1) then
         sr(1) = za(1)
         si(1) = za(2)
      else if (nslot .ge. 2) then
         call pmesplit (za,ya,sr(1),si(1),sr(2),si(2))
      end if
      if (nslot .eq. 3) then
         sr(3) = zb(1)
         si(3) = zb(2)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function dspstruc  --  weighted dispersion structure sum  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dspstruc" contracts the slot transforms "sr" and "si" at one
c     wavevector with the slot pair weights "w"
c
c
      function dspstruc (nslot,w,sr,si)
      implicit none
      integer i,j,nslot
      real*8 dspstruc
      real*8 w(3,3)
      real*8 sr(3),si(3)
c
c
c     sum the weighted products of the slot transforms
c
      dspstruc = 0.0d0
      do i = 1, nslot
         do j = 1, nslot
            dspstruc = dspstruc + w(j,i)*(sr(i)*sr(j)+si(i)*si(j))
         end do
      end do
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
c     for damped dispersion, weighting each site by the pair weight
c     of its own group, along with its lambda derivatives when "deriv"
c     is set and its partitioning over atoms when "part" is set
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
      integer i,ii,g
      real*8 e
      real*8 wg(0:2,0:2)
      real*8 dwg(0:2,0:2)
      real*8 d2wg(0:2,0:2)
      logical deriv,part
      logical dodl,dod2
c
c
c     get the group weights, along with their lambda derivatives
c     while the van der Waals lambda map has slope or curvature
c
      dodl = .false.
      dod2 = .false.
      if (deriv)  call vdlflags (dodl,dod2)
      call grpvwt (scexpd,wg,dwg,d2wg)
c
c     add the weighted self energy of each dispersion site
c
      do ii = 1, ndisp
         i = idisp(ii)
         e = csix(i)**2 * aewald**6 / 12.0d0
         g = mutg(i)
         if (g .ne. 0) then
            if (dodl)  devdl = devdl + dwg(g,g)*e
            if (dod2)  d2evdl2 = d2evdl2 + d2wg(g,g)*e
            e = e * wg(g,g)
         end if
         edsp = edsp + e
         if (part) then
            nedsp = nedsp + 1
            aedsp(i) = aedsp(i) + e
         end if
      end do
      return
      end

c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine edisp3  --  damped dispersion energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "edisp3" calculates the dispersion energy; also partitions
c     the energy among the atoms
c
c     literature reference:
c
c     J. A. Rackers, C. Liu, P. Ren and J. W. Ponder, "A Physically
c     Grounded Damped Dispersion Model with Particle Mesh Ewald
c     Summation", Journal of Chemical Physics, 149, 084115 (2018)
c
c
      subroutine edisp3
      use analyz
      use atoms
      use dsppot
      use energi
      use ewald
      use inform
      use iounit
      use limits
      implicit none
      integer i
      real*8 elrc,aelrc
      character*6 mode
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_dewald) then
         if (use_dlist) then
            call edisp3d
         else
            call edisp3c
         end if
      else
         if (use_dlist) then
            call edisp3b
         else
            call edisp3a
         end if
      end if
c
c     apply long range dispersion correction if desired
c
      if (use_dcorr .and. .not.use_dewald) then
         mode = 'DISP'
         call evcorr (mode,elrc)
         edsp = edsp + elrc
         aelrc = elrc / dble(n)
         do i = 1, n
            aedsp(i) = aedsp(i) + aelrc
         end do
         if (verbose .and. elrc.ne.0.0d0) then
            if (digits .ge. 8) then
               write (iout,10)  elrc
   10          format (/,' Long-Range Dispersion :',9x,f16.8)
            else if (digits .ge. 6) then
               write (iout,20)  elrc
   20          format (/,' Long-Range Dispersion :',9x,f16.6)
            else
               write (iout,30)  elrc
   30          format (/,' Long-Range Dispersion :',9x,f16.4)
            end if
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine edisp3a  --  damp dispersion analysis via loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "edisp3a" calculates the dispersion potential energy and
c     also partitions the energy among the atoms using a pairwise
c     double loop
c
c
      subroutine edisp3a
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use dsppot
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
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
      real*8 ai,ak
      real*8 damp,ddamp,taper
      real*8 dlsc,dlsq
      real*8 soft,dsoft
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical header,huge
      logical muti,mutk
      logical mutik
      character*6 mode
c
c
c     zero out the dispersion energy and partitioning terms
c
      nedsp = 0
      edsp = 0.0d0
      do i = 1, n
         aedsp(i) = 0.0d0
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
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. ndisp.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dispersion Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
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
                  call dampdisp (r,ai,ak,damp,ddamp)
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
c     increment the overall dispersion energy components
c
                  if (e .ne. 0.0d0) then
                     edsp = edsp + e
                     nedsp = nedsp + 1
                     aedsp(i) = aedsp(i) + 0.5d0*e
                     aedsp(k) = aedsp(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 4.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Dispersion',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,e
   30                format (' Disper',4x,2(i7,'-',a3),
     &                          9x,f10.4,2x,f12.4)
                  end if
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
c     increment the overall dispersion energy components
c
                     if (e .ne. 0.0d0) then
                        edsp = edsp + e
                        nedsp = nedsp + 1
                        aedsp(i) = aedsp(i) + 0.5d0*e
                        aedsp(k) = aedsp(k) + 0.5d0*e
                     end if
c
c     increment the total intermolecular energy
c
                     if (molcule(i) .ne. molcule(k)) then
                        einter = einter + e
                     end if
c
c     print a message if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 4.0d0)
                     if ((debug.and.e.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Dispersion',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',8x,'Energy',/)
                        end if
                        write (iout,50)  i,name(i),k,name(k),r,e
   50                   format (' Disper',4x,2(i7,'-',a3),1x,
     &                             '(XTAL)',2x,f10.4,2x,f12.4)
                     end if
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
c     ##  subroutine edisp3b  --  damp dispersion analysis via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "edisp3b" calculates the damped dispersion potential energy
c     and also partitions the energy among the atoms using a pairwise
c     neighbor list
c
c
      subroutine edisp3b
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use cell
      use disp
      use dsppot
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
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
      real*8 ai,ak
      real*8 damp,ddamp,taper
      real*8 dlsc,dlsq
      real*8 soft,dsoft
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical header,huge
      logical muti,mutk
      logical mutik
      character*6 mode
c
c
c     zero out the dispersion energy and partitioning terms
c
      nedsp = 0
      edsp = 0.0d0
      do i = 1, n
         aedsp(i) = 0.0d0
      end do
      if (ndisp .eq. 0)  return
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
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. ndisp.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dispersion Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ndisp,idisp,csix,adisp,use,
!$OMP& x,y,z,n12,n13,n14,n15,i12,i13,i14,i15,ndlst,dlst,use_group,
!$OMP& dsp2scale,dsp3scale,dsp4scale,dsp5scale,mutg,off2,cut2,c0,c1,
!$OMP& c2,c3,c4,c5,vcouple,dlsc,dlsq,rs2,molcule,name,verbose,debug,
!$OMP& header,iout)
!$OMP& firstprivate(dspscale),shared(edsp,nedsp,aedsp,einter)
!$OMP DO reduction(+:edsp,nedsp,aedsp,einter)
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
         do kk = 1, ndlst(i)
            k = dlst(kk,i)
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
                  call dampdisp (r,ai,ak,damp,ddamp)
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
c     increment the overall dispersion energy components
c
                  if (e .ne. 0.0d0) then
                     edsp = edsp + e
                     nedsp = nedsp + 1
                     aedsp(i) = aedsp(i) + 0.5d0*e
                     aedsp(k) = aedsp(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 4.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Dispersion',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,e
   30                format (' Disper',4x,2(i7,'-',a3),
     &                          9x,f10.4,2x,f12.4)
                  end if
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine edisp3c  --  Ewald dispersion analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "edisp3c" calculates the dispersion interaction energy using
c     particle mesh Ewald summation and a double loop
c
c
      subroutine edisp3c
      use action
      use analyz
      use atoms
      use disp
      use energi
      use ewald
      use pme
      implicit none
      integer i
c
c
c     zero out the dispersion energy and partitioning terms
c
      nedsp = 0
      edsp = 0.0d0
      do i = 1, n
         aedsp(i) = 0.0d0
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
      call edreal3c
c
c     compute the reciprocal space part of the Ewald summation
c
      call edrecip
c
c     compute the self-energy portion of the Ewald summation,
c     along with its partitioning over atoms
c
      call dspself (.false.,.true.)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine edreal3c  --  real space dispersion via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "edreal3c" calculates the real space portion of the damped
c     dispersion energy and analysis using Ewald and a double loop
c
c
      subroutine edreal3c
      use action
      use analyz
      use atomid
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
      use inform
      use inter
      use iounit
      use molcul
      use mutant
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,efull
      real*8 fgrp,ci,ck
      real*8 r,r2,r6
      real*8 ai,ak
      real*8 ralpha2
      real*8 expa,term
      real*8 damp,ddamp,scale
      real*8 dmp2,dhw,dhdw
      real*8 dlsc,dlsq,wsc
      real*8 soft,dsoft
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
      logical mutik
      logical header,huge
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
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. ndisp.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dispersion Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
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
c     compute the full undamped energy for this interaction
c
                  efull = e * scale
                  if (efull .ne. 0.0d0) then
                     nedsp = nedsp + 1
                     aedsp(i) = aedsp(i) + 0.5d0*efull
                     aedsp(k) = aedsp(k) + 0.5d0*efull
                     if (molcule(i) .ne. molcule(k)) then
                        einter = einter + efull
                     end if
                  end if
c
c     increment the overall dispersion energy component;
c     the lambda weight also scales the reciprocal space share
c
                  scale = scale - wsc
                  e = e * (wsc*expa+scale)
                  edsp = edsp + e
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(efull) .gt. 4.0d0)
                  if ((debug.and.efull.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Dispersion',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,efull
   30                format (' Disper',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
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
c     compute the full undamped energy for this interaction
c
                     if (i .eq. k)  e = 0.5d0 * e
                     efull = e * scale
                     if (efull .ne. 0.0d0) then
                        nedsp = nedsp + 1
                        aedsp(i) = aedsp(i) + 0.5d0*efull
                        aedsp(k) = aedsp(k) + 0.5d0*efull
                        if (molcule(i) .ne. molcule(k)) then
                           einter = einter + efull
                        end if
                     end if
c
c     increment the overall dispersion energy component;
c     the lambda weight also scales the reciprocal space share
c
                     scale = scale - wsc
                     e = e * (wsc*expa+scale)
                     edsp = edsp + e
c
c     print a message if the energy of this interaction is large
c
                     huge = (abs(efull) .gt. 4.0d0)
                     if ((debug.and.efull.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Dispersion',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',8x,'Energy',/)
                        end if
                        write (iout,50)  i,name(i),k,name(k),r,efull
   50                   format (' Disper',4x,2(i7,'-',a3),1x,
     &                             '(XTAL)',2x,f10.4,2x,f12.4)
                     end if
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine edisp3d  --  Ewald dispersion analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "edisp3d" calculates the damped dispersion energy and analysis
c     using particle mesh Ewald summation and a neighbor list
c
c
      subroutine edisp3d
      use action
      use analyz
      use atoms
      use disp
      use energi
      use ewald
      use pme
      implicit none
      integer i
c
c
c     zero out the dispersion energy and partitioning terms
c
      nedsp = 0
      edsp = 0.0d0
      do i = 1, n
         aedsp(i) = 0.0d0
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
      call edreal3d
c
c     compute the reciprocal space part of the Ewald summation
c
      call edrecip
c
c     compute the self-energy portion of the Ewald summation,
c     along with its partitioning over atoms
c
      call dspself (.false.,.true.)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine edreal3d  --  real space disp analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "edreal3d" calculates the real space portion of the damped
c     dispersion energy and analysis using Ewald and a neighbor list
c
c
      subroutine edreal3d
      use action
      use analyz
      use atomid
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
      use inform
      use inter
      use iounit
      use molcul
      use mutant
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,efull
      real*8 fgrp,ci,ck
      real*8 r,r2,r6
      real*8 ai,ak
      real*8 ralpha2
      real*8 expa,term
      real*8 damp,ddamp,scale
      real*8 dmp2,dhw,dhdw
      real*8 dlsc,dlsq,wsc
      real*8 soft,dsoft
      real*8 rs2
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
      logical mutik
      logical header,huge
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
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. ndisp.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dispersion Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(ndisp,idisp,csix,adisp,use,
!$OMP& x,y,z,n12,n13,n14,n15,i12,i13,i14,i15,ndlst,dlst,use_group,
!$OMP& dsp2scale,dsp3scale,dsp4scale,dsp5scale,mutg,off2,aewald,
!$OMP& molcule,vcouple,dlsc,dlsq,rs2,name,verbose,debug,header,iout)
!$OMP& shared(use_dhal,dhdamp)
!$OMP& firstprivate(dspscale),shared(edsp,nedsp,aedsp,einter)
!$OMP DO reduction(+:edsp,nedsp,aedsp,einter)
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
         do kk = 1, ndlst(i)
            k = dlst(kk,i)
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
c     compute the full undamped energy for this interaction
c
                  efull = e * scale
                  if (efull .ne. 0.0d0) then
                     nedsp = nedsp + 1
                     aedsp(i) = aedsp(i) + 0.5d0*efull
                     aedsp(k) = aedsp(k) + 0.5d0*efull
                     if (molcule(i) .ne. molcule(k)) then
                        einter = einter + efull
                     end if
                  end if
c
c     increment the overall dispersion energy component;
c     the lambda weight also scales the reciprocal space share
c
                  scale = scale - wsc
                  e = e * (wsc*expa+scale)
                  edsp = edsp + e
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(efull) .gt. 4.0d0)
                  if ((debug.and.efull.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Dispersion',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,efull
   30                format (' Disper',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
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

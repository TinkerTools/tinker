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
c     derivatives with respect to Cartesian coordinates, and lambda
c     derivatives while the van der Waals lambda map has slope or
c     curvature
c
c     literature reference:
c
c     J. A. Rackers, C. Liu, P. Ren and J. W. Ponder, "A Physically
c     Grounded Damped Dispersion Model with Particle Mesh Ewald
c     Summation", Journal of Chemical Physics, 149, 084115 (2018)
c
c
      subroutine edisp1
      use dlmda
      use dsppot
      use energi
      use ewald
      use limits
      use virial
      implicit none
      integer i,j
      real*8 elrc,vlrc
      real*8 delrc,dvlrc
      real*8 d2elrc
      logical dodl,dod2
      character*6 mode
c
c
c     take lambda derivatives as "vdlflags" sets for the van der
c     Waals lambda
c
      call vdlflags (dodl,dod2)
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
c     apply long range dispersion correction if desired, along
c     with its lambda derivatives
c
      if (use_dcorr .and. .not.use_dewald) then
         mode = 'DISP'
         call evcorr1 (mode,elrc,vlrc,delrc,dvlrc,d2elrc)
         edsp = edsp + elrc
         if (dodl)  devdl = devdl + delrc
         if (dod2)  d2evdl2 = d2evdl2 + d2elrc
         do i = 1, 3
            edspvir(i,i) = edspvir(i,i) + vlrc
            if (dod2)  devvirdl(i,i) = devvirdl(i,i) + dvlrc
         end do
      end if
c
c     add the damped dispersion virial to main virial
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + edspvir(j,i)
         end do
      end do
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
      use dlmda
      use dsppot
      use energi
      use group
      use mutant
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
      real*8 ai,ak
      real*8 damp,ddamp
      real*8 dlsc,dlsq
      real*8 dlsc1,dlsq1
      real*8 dlsc2,dlsq2
      real*8 soft,dsoft,dsdl
      real*8 rs2
      real*8 d2sdl2,d2sdldr
      real*8 dle,dle2,dlde
      real*8 dldx,dldy,dldz
      real*8 taper,dtaper
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
      logical mutik,mutdl,mutd2
      logical dodl,dod2
      character*6 mode
c
c
c     zero out the dispersion energy and derivatives
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
c     set soft core lambda scaling values for mutated interactions,
c     along with their lambda derivative factors as "vdlflags" sets
c     for the van der Waals lambda
c
      call vdlflags (dodl,dod2)
      if (nmut .ne. 0) then
         if (dodl) then
            call dspsc4 (dlsc,dlsq,dlsc1,dlsq1,dlsc2,dlsq2)
         else
            call dspsc (dlsc,dlsq)
            dlsc1 = 0.0d0
            dlsq1 = 0.0d0
            dlsc2 = 0.0d0
            dlsq2 = 0.0d0
         end if
      end if
      rs2 = scrdsp * scrdsp
c
c     set cutoff and switching coefficients
c
      mode = 'DISP'
      call switch (mode)
c
c     find dispersion energy and derivatives via double loop
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
                  de = -6.0d0 * e / r
c
c     find the damping factor for the dispersion interaction
c
                  call dampdisp (r,ai,ak,damp,ddamp)
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
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (muti.or.mutk) .and.
     &                    (vcouple.eq.1 .or. .not.(muti.and.mutk))
                  mutdl = mutik .and. dodl
                  mutd2 = mutik .and. dod2
c
c     get energy and derivative via soft core lambda scaling, along
c     with their lambda derivatives from the unscaled pair terms
c
                  if (mutik) then
                     soft = dlsc
                     dsoft = 0.0d0
                     dsdl = dlsc1
                     d2sdl2 = dlsc2
                     d2sdldr = 0.0d0
                     if (r2 .lt. rs2) then
                        if (mutdl) then
                           call edspsc4 (r,dlsc,dlsq,dlsc1,dlsq1,
     &                                   dlsc2,dlsq2,soft,dsoft,
     &                                   dsdl,d2sdl2,d2sdldr)
                        else
                           call edspsc (r,dlsc,dlsq,soft,dsoft)
                        end if
                     end if
                  end if
                  if (mutdl) then
                     dle = dsdl * e
                     if (mutd2) then
                        dle2 = d2sdl2 * e
                        dlde = dsdl*de + d2sdldr*e
                     end if
                  end if
                  if (mutik) then
                     de = soft*de + dsoft*e
                     e = soft * e
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
                     if (mutd2) then
                        dlde = dle*dtaper + dlde*taper
                        dle2 = dle2 * taper
                     end if
                     if (mutdl)  dle = dle * taper
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
                  edspvir(1,1) = edspvir(1,1) + vxx
                  edspvir(2,1) = edspvir(2,1) + vyx
                  edspvir(3,1) = edspvir(3,1) + vzx
                  edspvir(1,2) = edspvir(1,2) + vyx
                  edspvir(2,2) = edspvir(2,2) + vyy
                  edspvir(3,2) = edspvir(3,2) + vzy
                  edspvir(1,3) = edspvir(1,3) + vzx
                  edspvir(2,3) = edspvir(2,3) + vzy
                  edspvir(3,3) = edspvir(3,3) + vzz
c
c     increment the lambda derivatives of the energy, gradient
c     and virial for a mutated pair
c
                  if (mutdl)  devdl = devdl + dle
                  if (mutd2) then
                     d2evdl2 = d2evdl2 + dle2
                     dldx = dlde * xr/r
                     dldy = dlde * yr/r
                     dldz = dlde * zr/r
                     call dlcpair (i,k,xr,yr,zr,dldx,dldy,dldz,dfvdl,
     &                              devvirdl)
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
                     de = -6.0d0 * e / r
c
c     find the damping factor for the dispersion interaction
c
                     call dampdisp (r,ai,ak,damp,ddamp)
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
                     if (i .eq. k) then
                        e = 0.5d0 * e
                        de = 0.5d0 * de
                     end if
c
c     set use of lambda scaling for decoupling or annihilation
c
                     mutik = (muti.or.mutk) .and.
     &                       (vcouple.eq.1 .or. .not.(muti.and.mutk))
                     mutdl = mutik .and. dodl
                     mutd2 = mutik .and. dod2
c
c     get energy and derivative via soft core lambda scaling, along
c     with their lambda derivatives from the unscaled pair terms
c
                     if (mutik) then
                        soft = dlsc
                        dsoft = 0.0d0
                        dsdl = dlsc1
                        d2sdl2 = dlsc2
                        d2sdldr = 0.0d0
                        if (r2 .lt. rs2) then
                           if (mutdl) then
                              call edspsc4 (r,dlsc,dlsq,dlsc1,dlsq1,
     &                                      dlsc2,dlsq2,soft,dsoft,
     &                                      dsdl,d2sdl2,d2sdldr)
                           else
                              call edspsc (r,dlsc,dlsq,soft,dsoft)
                           end if
                        end if
                     end if
                     if (mutdl) then
                        dle = dsdl * e
                        if (mutd2) then
                           dle2 = d2sdl2 * e
                           dlde = dsdl*de + d2sdldr*e
                        end if
                     end if
                     if (mutik) then
                        de = soft*de + dsoft*e
                        e = soft * e
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
                        if (mutd2) then
                           dlde = dle*dtaper + dlde*taper
                           dle2 = dle2 * taper
                        end if
                        if (mutdl)  dle = dle * taper
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
                     edspvir(1,1) = edspvir(1,1) + vxx
                     edspvir(2,1) = edspvir(2,1) + vyx
                     edspvir(3,1) = edspvir(3,1) + vzx
                     edspvir(1,2) = edspvir(1,2) + vyx
                     edspvir(2,2) = edspvir(2,2) + vyy
                     edspvir(3,2) = edspvir(3,2) + vzy
                     edspvir(1,3) = edspvir(1,3) + vzx
                     edspvir(2,3) = edspvir(2,3) + vzy
                     edspvir(3,3) = edspvir(3,3) + vzz
c
c     increment the lambda derivatives of the energy, gradient
c     and virial for a mutated pair
c
                     if (mutdl)  devdl = devdl + dle
                     if (mutd2) then
                        d2evdl2 = d2evdl2 + dle2
                        dldx = dlde * xr/r
                        dldy = dlde * yr/r
                        dldz = dlde * zr/r
                        call dlcpair (i,k,xr,yr,zr,dldx,dldy,dldz,dfvdl,
     &                                 devvirdl)
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
      use dlmda
      use dsppot
      use energi
      use group
      use mutant
      use neigh
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,de,fgrp
      real*8 ci,ck
      real*8 r,r2,r3
      real*8 r4,r5,r6
      real*8 ai,ak
      real*8 damp,ddamp
      real*8 dlsc,dlsq
      real*8 dlsc1,dlsq1
      real*8 dlsc2,dlsq2
      real*8 soft,dsoft,dsdl
      real*8 rs2
      real*8 d2sdl2,d2sdldr
      real*8 dle,dle2,dlde
      real*8 dldx,dldy,dldz
      real*8 taper,dtaper
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: dspscale(:)
      real*8, allocatable :: dfl(:,:)
      logical proceed,usei
      logical muti,mutk
      logical mutik,mutdl,mutd2
      logical dodl,dod2
      character*6 mode
c
c
c     zero out the dispersion energy and derivatives
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
c     set soft core lambda scaling values for mutated interactions,
c     along with their lambda derivative factors as "vdlflags" sets
c     for the van der Waals lambda
c
      call vdlflags (dodl,dod2)
c
c     the lambda derivative of the force is reduced over threads in
c     a local array that only has full size when it is needed
c
      if (dod2) then
         allocate (dfl(3,n))
      else
         allocate (dfl(3,1))
      end if
      dfl = 0.0d0
      if (nmut .ne. 0) then
         if (dodl) then
            call dspsc4 (dlsc,dlsq,dlsc1,dlsq1,dlsc2,dlsq2)
         else
            call dspsc (dlsc,dlsq)
            dlsc1 = 0.0d0
            dlsq1 = 0.0d0
            dlsc2 = 0.0d0
            dlsq2 = 0.0d0
         end if
      end if
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
!$OMP& c0,c1,c2,c3,c4,c5,vcouple,dlsc,dlsq,rs2,dlsc1,dlsq1,dlsc2,dlsq2,
!$OMP& dodl,dod2)
!$OMP& firstprivate(dspscale) shared(edsp,dedsp,edspvir,devdl,
!$OMP& d2evdl2,dfl,devvirdl)
!$OMP DO reduction(+:edsp,dedsp,edspvir,devdl,d2evdl2,dfl,devvirdl)
c
c     find dispersion energy and derivatives via neighbor list
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
                  de = -6.0d0 * e / r
c
c     find the damping factor for the dispersion interaction
c
                  call dampdisp (r,ai,ak,damp,ddamp)
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
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (muti.or.mutk) .and.
     &                    (vcouple.eq.1 .or. .not.(muti.and.mutk))
                  mutdl = mutik .and. dodl
                  mutd2 = mutik .and. dod2
c
c     get energy and derivative via soft core lambda scaling, along
c     with their lambda derivatives from the unscaled pair terms
c
                  if (mutik) then
                     soft = dlsc
                     dsoft = 0.0d0
                     dsdl = dlsc1
                     d2sdl2 = dlsc2
                     d2sdldr = 0.0d0
                     if (r2 .lt. rs2) then
                        if (mutdl) then
                           call edspsc4 (r,dlsc,dlsq,dlsc1,dlsq1,
     &                                   dlsc2,dlsq2,soft,dsoft,
     &                                   dsdl,d2sdl2,d2sdldr)
                        else
                           call edspsc (r,dlsc,dlsq,soft,dsoft)
                        end if
                     end if
                  end if
                  if (mutdl) then
                     dle = dsdl * e
                     if (mutd2) then
                        dle2 = d2sdl2 * e
                        dlde = dsdl*de + d2sdldr*e
                     end if
                  end if
                  if (mutik) then
                     de = soft*de + dsoft*e
                     e = soft * e
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
                     if (mutd2) then
                        dlde = dle*dtaper + dlde*taper
                        dle2 = dle2 * taper
                     end if
                     if (mutdl)  dle = dle * taper
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
                  edspvir(1,1) = edspvir(1,1) + vxx
                  edspvir(2,1) = edspvir(2,1) + vyx
                  edspvir(3,1) = edspvir(3,1) + vzx
                  edspvir(1,2) = edspvir(1,2) + vyx
                  edspvir(2,2) = edspvir(2,2) + vyy
                  edspvir(3,2) = edspvir(3,2) + vzy
                  edspvir(1,3) = edspvir(1,3) + vzx
                  edspvir(2,3) = edspvir(2,3) + vzy
                  edspvir(3,3) = edspvir(3,3) + vzz
c
c     increment the lambda derivatives of the energy, gradient
c     and virial for a mutated pair
c
                  if (mutdl)  devdl = devdl + dle
                  if (mutd2) then
                     d2evdl2 = d2evdl2 + dle2
                     dldx = dlde * xr/r
                     dldy = dlde * yr/r
                     dldz = dlde * zr/r
                     call dlcpair (i,k,xr,yr,zr,dldx,dldy,dldz,dfl,
     &                              devvirdl)
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
c     add the lambda derivative of the force to the global array
c
      if (dod2) then
         do i = 1, n
            dfvdl(1,i) = dfvdl(1,i) + dfl(1,i)
            dfvdl(2,i) = dfvdl(2,i) + dfl(2,i)
            dfvdl(3,i) = dfvdl(3,i) + dfl(3,i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dspscale)
      deallocate (dfl)
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
c
c
c     zero out the dispersion energy and derivatives
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
c     compute the self-energy portion of the Ewald summation,
c     along with its lambda derivatives
c
      call dspself (.true.,.false.)
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
      use dlmda
      use dsppot
      use deriv
      use energi
      use ewald
      use group
      use mutant
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
      real*8 ai,ak
      real*8 damp,ddamp
      real*8 dmp2,dmpd
      real*8 dhw,dhdw
      real*8 wsc
      real*8 dlsc,dlsq
      real*8 dlsc1,dlsq1
      real*8 dlsc2,dlsq2
      real*8 soft,dsoft,dsdl
      real*8 rs2
      real*8 d2sdl2,d2sdldr
      real*8 dle,dle2,dlde
      real*8 dldx,dldy,dldz
      real*8 sc,dscale
      real*8 ralpha2,scale
      real*8 expterm,term
      real*8 expa,rterm
      real*8 dedx,dedy,dedz
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
      real*8, allocatable :: dspscale(:)
      logical proceed,usei
      logical muti,mutk
      logical mutik,mutdl,mutd2
      logical dodl,dod2
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
c     set soft core lambda scaling values for mutated interactions,
c     along with their lambda derivative factors as "vdlflags" sets
c     for the van der Waals lambda
c
      call vdlflags (dodl,dod2)
      if (nmut .ne. 0) then
         if (dodl) then
            call dspsc4 (dlsc,dlsq,dlsc1,dlsq1,dlsc2,dlsq2)
         else
            call dspsc (dlsc,dlsq)
            dlsc1 = 0.0d0
            dlsq1 = 0.0d0
            dlsc2 = 0.0d0
            dlsq2 = 0.0d0
         end if
      end if
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
                  r6 = r2**3
                  e = -ci * ck / r6
                  ralpha2 = r2 * aewald**2
                  term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                  expterm = exp(-ralpha2)
                  expa = expterm * term
c
c     find the damping factor, which an undamped hal tail skips
c
                  r = sqrt(r2)
                  r7 = r6 * r
                  damp = 1.0d0
                  ddamp = 0.0d0
                  if (dhdamp .or. .not.use_dhal) then
                     call dampdisp (r,ai,ak,damp,ddamp)
                  end if
                  dmp2 = damp * damp
                  dmpd = damp * ddamp
c
c     an Ewald tail mapped from hal takes over across the vdw switch
c
                  if (use_dhal) then
                     call dhalwt (r,dhw,dhdw)
                     dmpd = dhw*dmpd + 0.5d0*dhdw*dmp2
                     dmp2 = dhw * dmp2
                  end if
c
c     apply damping and scaling factors for this interaction
c
                  sc = dspscale(k)
                  if (use_group)  sc = sc * fgrp
                  scale = sc * dmp2
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (muti.or.mutk) .and.
     &                    (vcouple.eq.1 .or. .not.(muti.and.mutk))
                  mutdl = mutik .and. dodl
                  mutd2 = mutik .and. dod2
c
c     get the soft core scaling and its lambda derivatives, where
c     the lambda power of the soft core weights the reciprocal space
c     share as it does in the reciprocal space sum
c
                  wsc = 1.0d0
                  soft = 1.0d0
                  dsoft = 0.0d0
                  if (mutik) then
                     wsc = dlsc
                     soft = dlsc
                     dsdl = dlsc1
                     d2sdl2 = dlsc2
                     d2sdldr = 0.0d0
                     if (r2 .lt. rs2) then
                        if (mutdl) then
                           call edspsc4 (r,dlsc,dlsq,dlsc1,dlsq1,
     &                                   dlsc2,dlsq2,soft,dsoft,
     &                                   dsdl,d2sdl2,d2sdldr)
                        else
                           call edspsc (r,dlsc,dlsq,soft,dsoft)
                        end if
                     end if
                  end if
                  rterm = -(ralpha2**3) * expterm / r
                  if (mutdl) then
                     dle = e * (dsdl*scale+dlsc1*(expa-1.0d0))
                     if (mutd2) then
                        dle2 = e * (d2sdl2*scale+dlsc2*(expa-1.0d0))
                        dlde = -6.0d0*dle/r2 - dlsc1*ci*ck*rterm/r7
     &                      - 2.0d0*ci*ck*sc*dmpd*dsdl/r7
     &                      - ci*ck*scale*d2sdldr/r7
                     end if
                  end if
                  dscale = scale * dsoft
c
c     increment the overall damped dispersion energy component;
c     the lambda weight also scales the reciprocal space share
c
                  scale = soft*scale - wsc
                  e = e * (wsc*expa+scale)
                  edsp = edsp + e
c
c     increment the damped dispersion derivative components
c
                  de = -6.0d0*e/r2 - wsc*ci*ck*rterm/r7
     &                    - 2.0d0*ci*ck*sc*dmpd*soft/r7
     &                    - ci*ck*dscale/r7
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
                  edspvir(1,1) = edspvir(1,1) + vxx
                  edspvir(2,1) = edspvir(2,1) + vyx
                  edspvir(3,1) = edspvir(3,1) + vzx
                  edspvir(1,2) = edspvir(1,2) + vyx
                  edspvir(2,2) = edspvir(2,2) + vyy
                  edspvir(3,2) = edspvir(3,2) + vzy
                  edspvir(1,3) = edspvir(1,3) + vzx
                  edspvir(2,3) = edspvir(2,3) + vzy
                  edspvir(3,3) = edspvir(3,3) + vzz
c
c     increment the lambda derivatives of the energy, gradient
c     and virial for a mutated pair
c
                  if (mutdl)  devdl = devdl + dle
                  if (mutd2) then
                     d2evdl2 = d2evdl2 + dle2
                     dldx = dlde * xr
                     dldy = dlde * yr
                     dldz = dlde * zr
                     call dlcpair (i,k,xr,yr,zr,dldx,dldy,dldz,dfvdl,
     &                              devvirdl)
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
                     r6 = r2**3
                     e = -ci * ck / r6
                     ralpha2 = r2 * aewald**2
                     term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                     expterm = exp(-ralpha2)
                     expa = expterm * term
c
c     find the damping factor, which an undamped hal tail skips
c
                     r = sqrt(r2)
                     r7 = r6 * r
                     damp = 1.0d0
                     ddamp = 0.0d0
                     if (dhdamp .or. .not.use_dhal) then
                        call dampdisp (r,ai,ak,damp,ddamp)
                     end if
                     dmp2 = damp * damp
                     dmpd = damp * ddamp
c
c     an Ewald tail mapped from hal takes over across the vdw switch
c
                     if (use_dhal) then
                        call dhalwt (r,dhw,dhdw)
                        dmpd = dhw*dmpd + 0.5d0*dhdw*dmp2
                        dmp2 = dhw * dmp2
                     end if
c
c     apply damping and scaling factors for this interaction
c
                     sc = 1.0d0
                     if (use_group)  sc = sc * fgrp
                     if (use_polymer) then
                        if (r2 .le. polycut2)  sc = sc * dspscale(k)
                     end if
                     scale = sc * dmp2
c
c     set use of lambda scaling for decoupling or annihilation
c
                     mutik = (muti.or.mutk) .and.
     &                       (vcouple.eq.1 .or. .not.(muti.and.mutk))
                     mutdl = mutik .and. dodl
                     mutd2 = mutik .and. dod2
c
c     get the soft core scaling and its lambda derivatives, where
c     the lambda power of the soft core weights the reciprocal space
c     share as it does in the reciprocal space sum
c
                     wsc = 1.0d0
                     soft = 1.0d0
                     dsoft = 0.0d0
                     if (mutik) then
                        wsc = dlsc
                        soft = dlsc
                        dsdl = dlsc1
                        d2sdl2 = dlsc2
                        d2sdldr = 0.0d0
                        if (r2 .lt. rs2) then
                           if (mutdl) then
                              call edspsc4 (r,dlsc,dlsq,dlsc1,dlsq1,
     &                                      dlsc2,dlsq2,soft,dsoft,
     &                                      dsdl,d2sdl2,d2sdldr)
                           else
                              call edspsc (r,dlsc,dlsq,soft,dsoft)
                           end if
                        end if
                     end if
                     rterm = -(ralpha2**3) * expterm / r
                     if (mutdl) then
                        dle = e * (dsdl*scale+dlsc1*(expa-1.0d0))
                        if (mutd2) then
                           dle2 = e * (d2sdl2*scale+dlsc2*(expa-1.0d0))
                           dlde = -6.0d0*dle/r2 - dlsc1*ci*ck*rterm/r7
     &                      - 2.0d0*ci*ck*sc*dmpd*dsdl/r7
     &                      - ci*ck*scale*d2sdldr/r7
                        end if
                     end if
                     dscale = scale * dsoft
c
c     increment the overall damped dispersion energy component;
c     the lambda weight also scales the reciprocal space share
c
                     scale = soft*scale - wsc
                     e = e * (wsc*expa+scale)
                     de = -6.0d0*e/r2 - wsc*ci*ck*rterm/r7
     &                      - 2.0d0*ci*ck*sc*dmpd*soft/r7
     &                      - ci*ck*dscale/r7
                     if (ii .eq. kk) then
                        e = 0.5d0 * e
                        de = 0.5d0 * de
                        if (mutdl)  dle = 0.5d0 * dle
                        if (mutd2) then
                           dle2 = 0.5d0 * dle2
                           dlde = 0.5d0 * dlde
                        end if
                     end if
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
                     edspvir(1,1) = edspvir(1,1) + vxx
                     edspvir(2,1) = edspvir(2,1) + vyx
                     edspvir(3,1) = edspvir(3,1) + vzx
                     edspvir(1,2) = edspvir(1,2) + vyx
                     edspvir(2,2) = edspvir(2,2) + vyy
                     edspvir(3,2) = edspvir(3,2) + vzy
                     edspvir(1,3) = edspvir(1,3) + vzx
                     edspvir(2,3) = edspvir(2,3) + vzy
                     edspvir(3,3) = edspvir(3,3) + vzz
c
c     increment the lambda derivatives of the energy, gradient
c     and virial for a mutated pair
c
                     if (mutdl)  devdl = devdl + dle
                     if (mutd2) then
                        d2evdl2 = d2evdl2 + dle2
                        dldx = dlde * xr
                        dldy = dlde * yr
                        dldz = dlde * zr
                        call dlcpair (i,k,xr,yr,zr,dldx,dldy,dldz,dfvdl,
     &                                 devvirdl)
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
c
c
c     zero out the dispersion energy and derivatives
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
c     compute the self-energy portion of the Ewald summation,
c     along with its lambda derivatives
c
      call dspself (.true.,.false.)
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
c     "edreal1d" calculates the real space portion of the Ewald
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
      use dlmda
      use dsppot
      use energi
      use ewald
      use group
      use mutant
      use neigh
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 e,de,fgrp
      real*8 ci,ck
      real*8 r,r2,r6,r7
      real*8 ai,ak
      real*8 damp,ddamp
      real*8 dmp2,dmpd
      real*8 dhw,dhdw
      real*8 wsc
      real*8 dlsc,dlsq
      real*8 dlsc1,dlsq1
      real*8 dlsc2,dlsq2
      real*8 soft,dsoft,dsdl
      real*8 rs2
      real*8 d2sdl2,d2sdldr
      real*8 dle,dle2,dlde
      real*8 dldx,dldy,dldz
      real*8 sc,dscale
      real*8 ralpha2,scale
      real*8 expterm,term
      real*8 expa,rterm
      real*8 dedx,dedy,dedz
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
      real*8, allocatable :: dspscale(:)
      real*8, allocatable :: dfl(:,:)
      logical proceed,usei
      logical muti,mutk
      logical mutik,mutdl,mutd2
      logical dodl,dod2
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
c     set soft core lambda scaling values for mutated interactions,
c     along with their lambda derivative factors as "vdlflags" sets
c     for the van der Waals lambda
c
      call vdlflags (dodl,dod2)
c
c     the lambda derivative of the force is reduced over threads in
c     a local array that only has full size when it is needed
c
      if (dod2) then
         allocate (dfl(3,n))
      else
         allocate (dfl(3,1))
      end if
      dfl = 0.0d0
      if (nmut .ne. 0) then
         if (dodl) then
            call dspsc4 (dlsc,dlsq,dlsc1,dlsq1,dlsc2,dlsq2)
         else
            call dspsc (dlsc,dlsq)
            dlsc1 = 0.0d0
            dlsq1 = 0.0d0
            dlsc2 = 0.0d0
            dlsq2 = 0.0d0
         end if
      end if
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
!$OMP& vcouple,dlsc,dlsq,rs2,dlsc1,dlsq1,dlsc2,dlsq2,dodl,dod2)
!$OMP& shared(use_dhal,dhdamp)
!$OMP& firstprivate(dspscale) shared(edsp,dedsp,edspvir,devdl,
!$OMP& d2evdl2,dfl,devvirdl)
!$OMP DO reduction(+:edsp,dedsp,edspvir,devdl,d2evdl2,dfl,devvirdl)
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
                  r6 = r2**3
                  e = -ci * ck / r6
                  ralpha2 = r2 * aewald**2
                  term = 1.0d0 + ralpha2 + 0.5d0*ralpha2**2
                  expterm = exp(-ralpha2)
                  expa = expterm * term
c
c     find the damping factor, which an undamped hal tail skips
c
                  r = sqrt(r2)
                  r7 = r6 * r
                  damp = 1.0d0
                  ddamp = 0.0d0
                  if (dhdamp .or. .not.use_dhal) then
                     call dampdisp (r,ai,ak,damp,ddamp)
                  end if
                  dmp2 = damp * damp
                  dmpd = damp * ddamp
c
c     an Ewald tail mapped from hal takes over across the vdw switch
c
                  if (use_dhal) then
                     call dhalwt (r,dhw,dhdw)
                     dmpd = dhw*dmpd + 0.5d0*dhdw*dmp2
                     dmp2 = dhw * dmp2
                  end if
c
c     apply damping and scaling factors for this interaction
c
                  sc = dspscale(k)
                  if (use_group)  sc = sc * fgrp
                  scale = sc * dmp2
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (muti.or.mutk) .and.
     &                    (vcouple.eq.1 .or. .not.(muti.and.mutk))
                  mutdl = mutik .and. dodl
                  mutd2 = mutik .and. dod2
c
c     get the soft core scaling and its lambda derivatives, where
c     the lambda power of the soft core weights the reciprocal space
c     share as it does in the reciprocal space sum
c
                  wsc = 1.0d0
                  soft = 1.0d0
                  dsoft = 0.0d0
                  if (mutik) then
                     wsc = dlsc
                     soft = dlsc
                     dsdl = dlsc1
                     d2sdl2 = dlsc2
                     d2sdldr = 0.0d0
                     if (r2 .lt. rs2) then
                        if (mutdl) then
                           call edspsc4 (r,dlsc,dlsq,dlsc1,dlsq1,
     &                                   dlsc2,dlsq2,soft,dsoft,
     &                                   dsdl,d2sdl2,d2sdldr)
                        else
                           call edspsc (r,dlsc,dlsq,soft,dsoft)
                        end if
                     end if
                  end if
                  rterm = -(ralpha2**3) * expterm / r
                  if (mutdl) then
                     dle = e * (dsdl*scale+dlsc1*(expa-1.0d0))
                     if (mutd2) then
                        dle2 = e * (d2sdl2*scale+dlsc2*(expa-1.0d0))
                        dlde = -6.0d0*dle/r2 - dlsc1*ci*ck*rterm/r7
     &                      - 2.0d0*ci*ck*sc*dmpd*dsdl/r7
     &                      - ci*ck*scale*d2sdldr/r7
                     end if
                  end if
                  dscale = scale * dsoft
c
c     increment the overall damped dispersion energy component;
c     the lambda weight also scales the reciprocal space share
c
                  scale = soft*scale - wsc
                  e = e * (wsc*expa+scale)
                  edsp = edsp + e
c
c     increment the damped dispersion derivative components
c
                  de = -6.0d0*e/r2 - wsc*ci*ck*rterm/r7
     &                    - 2.0d0*ci*ck*sc*dmpd*soft/r7
     &                    - ci*ck*dscale/r7
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
                  edspvir(1,1) = edspvir(1,1) + vxx
                  edspvir(2,1) = edspvir(2,1) + vyx
                  edspvir(3,1) = edspvir(3,1) + vzx
                  edspvir(1,2) = edspvir(1,2) + vyx
                  edspvir(2,2) = edspvir(2,2) + vyy
                  edspvir(3,2) = edspvir(3,2) + vzy
                  edspvir(1,3) = edspvir(1,3) + vzx
                  edspvir(2,3) = edspvir(2,3) + vzy
                  edspvir(3,3) = edspvir(3,3) + vzz
c
c     increment the lambda derivatives of the energy, gradient
c     and virial for a mutated pair
c
                  if (mutdl)  devdl = devdl + dle
                  if (mutd2) then
                     d2evdl2 = d2evdl2 + dle2
                     dldx = dlde * xr
                     dldy = dlde * yr
                     dldz = dlde * zr
                     call dlcpair (i,k,xr,yr,zr,dldx,dldy,dldz,dfl,
     &                              devvirdl)
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
c     add the lambda derivative of the force to the global array
c
      if (dod2) then
         do i = 1, n
            dfvdl(1,i) = dfvdl(1,i) + dfl(1,i)
            dfvdl(2,i) = dfvdl(2,i) + dfl(2,i)
            dfvdl(3,i) = dfvdl(3,i) + dfl(3,i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dspscale)
      deallocate (dfl)
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
      use boxes
      use bound
      use disp
      use deriv
      use dlmda
      use energi
      use ewald
      use math
      use mutant
      use pme
      use virial
      implicit none
      integer i,j,k,ii
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff,ntot
      integer i0,iatm,igrd0
      integer ipart
      integer it1,it2,it3
      integer j0,jgrd0
      integer k0,kgrd0
      real*8 e,fi,denom
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 term,denom0
      real*8 eterm,vterm
      real*8 expterm
      real*8 erfcterm
      real*8 hsq,struc2
      real*8 h,hhh,b,bfac
      real*8 fac1,fac2,fac3
      real*8 de1,de2,de3
      real*8 dn1,dn2,dn3
      real*8 dt1,dt2,dt3
      real*8 t1,t2,t3
      real*8 wel,wll,cpr
      real*8 dwel,dwll,dcpr
      real*8 d2wel,d2wll,d2cpr
      real*8 ser,sei,slr,sli
      real*8 strucd,strucd2
      real*8 edl,vtermd
      real*8 dl1,dl2,dl3
      real*8, allocatable :: qfwd(:,:,:,:)
      real*8, allocatable :: qdl(:,:,:,:)
      logical split
      logical dodl,dod2
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
c     get the mutated pair weights, along with their lambda
c     derivatives while the van der Waals lambda map has slope or
c     curvature, as in "edreal1c"
c
      call vdlflags (dodl,dod2)
      if (dodl) then
         call dspmutwt4 (split,wel,wll,cpr,dwel,dwll,dcpr,
     &                   d2wel,d2wll,d2cpr)
      else
         call dspmutwt (split,wel,wll,cpr)
      end if
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_disp (split)
      call fftfront (qgrid)
c
c     keep the packed transform, since the grid is overwritten
c     in place but each point needs its opposite wavevector
c
      if (split) then
         allocate (qfwd(2,nfft1,nfft2,nfft3))
         qfwd = qgrid
      end if
c
c     the lambda derivative of the potential uses a second grid
c
      if (dod2) then
         allocate (qdl(2,nfft1,nfft2,nfft3))
         qdl(1,1,1,1) = 0.0d0
         qdl(2,1,1,1) = 0.0d0
      end if
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
         eterm = 0.0d0
         denom = denom0*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
         if (split)  call dspsplit (qfwd,k1,k2,k3,ser,sei,slr,sli)
         if (term .gt. -50.0d0) then
            expterm = exp(term)
            erfcterm = erfc(b)
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               erfcterm = erfcterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (nonprism) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               if (mod(m1+m2+m3,2) .ne. 0)  erfcterm = 0.0d0
            end if
            eterm = (-fac1*erfcterm*hhh - expterm*(fac2+fac3*hsq))/denom
            if (split) then
               struc2 = ser*ser + sei*sei + wll*(slr*slr+sli*sli)
     &                     + 2.0d0*wel*(ser*slr+sei*sli)
            else
               struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
            end if
            vterm = 3.0d0*(fac1*erfcterm*h + fac3*expterm)*struc2/denom
            e = eterm * struc2
            edsp = edsp + e
            edspvir(1,1) = edspvir(1,1) + h1*h1*vterm - e
            edspvir(2,1) = edspvir(2,1) + h1*h2*vterm
            edspvir(3,1) = edspvir(3,1) + h1*h3*vterm
            edspvir(1,2) = edspvir(1,2) + h2*h1*vterm
            edspvir(2,2) = edspvir(2,2) + h2*h2*vterm - e
            edspvir(3,2) = edspvir(3,2) + h2*h3*vterm
            edspvir(1,3) = edspvir(1,3) + h3*h1*vterm
            edspvir(2,3) = edspvir(2,3) + h3*h2*vterm
            edspvir(3,3) = edspvir(3,3) + h3*h3*vterm - e
c
c     get the lambda derivatives of the energy and virial from the
c     derivatives of the pair weights
c
            if (dodl) then
               strucd = dwll*(slr*slr+sli*sli)
     &                     + 2.0d0*dwel*(ser*slr+sei*sli)
               devdl = devdl + eterm*strucd
            end if
            if (dod2) then
               strucd2 = d2wll*(slr*slr+sli*sli)
     &                      + 2.0d0*d2wel*(ser*slr+sei*sli)
               d2evdl2 = d2evdl2 + eterm*strucd2
               vtermd = 3.0d0 * (fac1*erfcterm*h+fac3*expterm)
     &                     * strucd / denom
               edl = eterm * strucd
               devvirdl(1,1) = devvirdl(1,1) + h1*h1*vtermd - edl
               devvirdl(2,1) = devvirdl(2,1) + h1*h2*vtermd
               devvirdl(3,1) = devvirdl(3,1) + h1*h3*vtermd
               devvirdl(1,2) = devvirdl(1,2) + h2*h1*vtermd
               devvirdl(2,2) = devvirdl(2,2) + h2*h2*vtermd - edl
               devvirdl(3,2) = devvirdl(3,2) + h2*h3*vtermd
               devvirdl(1,3) = devvirdl(1,3) + h3*h1*vtermd
               devvirdl(2,3) = devvirdl(2,3) + h3*h2*vtermd
               devvirdl(3,3) = devvirdl(3,3) + h3*h3*vtermd - edl
            end if
         end if
c
c     the environment potential takes the real part of the grid
c     and the ligand potential the imaginary part
c
         if (split) then
            qgrid(1,k1,k2,k3) = eterm * (ser-wel*sei+wel*slr-wll*sli)
            qgrid(2,k1,k2,k3) = eterm * (sei+wel*ser+wel*sli+wll*slr)
         else
            qgrid(1,k1,k2,k3) = eterm * qgrid(1,k1,k2,k3) 
            qgrid(2,k1,k2,k3) = eterm * qgrid(2,k1,k2,k3)
         end if
c
c     the lambda derivative of the potential packs the environment
c     and ligand parts with the derivatives of the pair weights
c
         if (dod2) then
            qdl(1,k1,k2,k3) = eterm * (dwel*slr-dwel*sei-dwll*sli)
            qdl(2,k1,k2,k3) = eterm * (dwel*ser+dwel*sli+dwll*slr)
         end if
      end do
      if (split)  deallocate (qfwd)
c
c     perform the 3-D FFT backward transformation
c
      call fftback (qgrid)
      if (dod2)  call fftback (qdl)
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
         fi = csix(iatm)
         ipart = 1
         if (split) then
            if (mutg(iatm) .ne. 0)  ipart = 2
         end if
         de1 = 0.0d0
         de2 = 0.0d0
         de3 = 0.0d0
         dl1 = 0.0d0
         dl2 = 0.0d0
         dl3 = 0.0d0
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
                  term = qgrid(ipart,i,j,k)
                  de1 = de1 + 2.0d0*term*dt1*t2*t3
                  de2 = de2 + 2.0d0*term*dt2*t1*t3
                  de3 = de3 + 2.0d0*term*dt3*t1*t2
                  if (dod2) then
                     term = qdl(ipart,i,j,k)
                     dl1 = dl1 + 2.0d0*term*dt1*t2*t3
                     dl2 = dl2 + 2.0d0*term*dt2*t1*t3
                     dl3 = dl3 + 2.0d0*term*dt3*t1*t2
                  end if
               end do
            end do
         end do
         dedsp(1,iatm) = dedsp(1,iatm) + fi*(recip(1,1)*de1
     &                            +recip(1,2)*de2+recip(1,3)*de3)
         dedsp(2,iatm) = dedsp(2,iatm) + fi*(recip(2,1)*de1
     &                            +recip(2,2)*de2+recip(2,3)*de3)
         dedsp(3,iatm) = dedsp(3,iatm) + fi*(recip(3,1)*de1
     &                            +recip(3,2)*de2+recip(3,3)*de3)
         if (dod2) then
            dfvdl(1,iatm) = dfvdl(1,iatm) + fi*(recip(1,1)*dl1
     &                         +recip(1,2)*dl2+recip(1,3)*dl3)
            dfvdl(2,iatm) = dfvdl(2,iatm) + fi*(recip(2,1)*dl1
     &                         +recip(2,2)*dl2+recip(2,3)*dl3)
            dfvdl(3,iatm) = dfvdl(3,iatm) + fi*(recip(3,1)*dl1
     &                         +recip(3,2)*dl2+recip(3,3)*dl3)
         end if
      end do
      if (dod2)  deallocate (qdl)
c
c     account for the energy and virial correction terms
c
      term = cpr * aewald**3 / denom0
      edsp = edsp - term
      edspvir(1,1) = edspvir(1,1) + term
      edspvir(2,2) = edspvir(2,2) + term
      edspvir(3,3) = edspvir(3,3) + term
      if (dodl) then
         term = dcpr * aewald**3 / denom0
         devdl = devdl - term
         if (dod2) then
            d2evdl2 = d2evdl2 - d2cpr*aewald**3/denom0
            devvirdl(1,1) = devvirdl(1,1) + term
            devvirdl(2,2) = devvirdl(2,2) + term
            devvirdl(3,3) = devvirdl(3,3) + term
         end if
      end if
      return
      end

c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emrecip1  --  mpole Ewald recip energy & derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to atomic multipole
c     interactions and dipole polarizability
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine emrecip1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'pme.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'virial.i'
      integer i,j,k,m
      integer ii,i1,i2
      integer j1,j2,j3
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer nframes,numfrlst
      integer frame_axis(3)
      integer qi1(6),qi2(6)
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      integer framelist(5,3*maxatm)
      real*8 e,r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 dfx,dfy,dfz
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 tfield(10)
      real*8 tdfield(10)
      real*8 tpfield(10)
      real*8 a(3,3),ftc(10,10)
      real*8 frc(3,maxatm)
      real*8 frci(3,maxatm)
      real*8 frct(3,maxatm)
      real*8 frcit(3,maxatm)
      real*8 fuind(3,maxatm)
      real*8 fuinp(3,maxatm)
      real*8 cmp(10,maxatm)
      real*8 fmp(10,maxatm)
      real*8 cphi(10,maxatm)
      real*8 fphi(20,maxatm)
      real*8 dipfield1(3,maxatm)
      real*8 dipfield2(3,maxatm)
      real*8 de_drotframe(3,maxatm)
      real*8 de_drotsite(3,maxatm)
      real*8 fdip_phi1(10,maxatm)
      real*8 fdip_phi2(10,maxatm)
      real*8 fdip_sum_phi(20,maxatm)
      real*8 frame(3,3,maxatm)
      real*8 frdefpts(3,2,maxatm)
      real*8 de_ddefpt(3,2,maxatm)
      real*8 qgrip(2,maxgrid,maxgrid,maxgrid)
c
c     quadrupole indices for fractional to Cartesian transform
c
      data qi1  / 1, 2, 3, 1, 1, 2 /
      data qi2  / 1, 2, 3, 2, 3, 3 /
c
c     derivative indices into the fphi and fdip_phi arrays
c
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vyx = 0.0d0
      vzx = 0.0d0
      vyy = 0.0d0
      vzy = 0.0d0
      vzz = 0.0d0
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
         frc(1,i) = 0.0d0
         frc(2,i) = 0.0d0
         frc(3,i) = 0.0d0
         frci(1,i) = 0.0d0
         frci(2,i) = 0.0d0
         frci(3,i) = 0.0d0
      end do
      do i = 1, n
         frct(1,i) = 0.0d0
         frct(2,i) = 0.0d0
         frct(3,i) = 0.0d0
         frcit(1,i) = 0.0d0
         frcit(2,i) = 0.0d0
         frcit(3,i) = 0.0d0
      end do
c
c     assemble fractional to Cartesian transformation matrix
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
      end do
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
c
c     handle the Z-then-X coordinate ordering
c
      frame_axis(1) = 3
      frame_axis(2) = 1
      frame_axis(3) = 2
c
c     compute the arrays of B-spline coefficients
c
      if (.not. use_polar)  call pme_bspline_fill
c
c     assign permanent and induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c
      if (use_polar) then
         do i = 1, npole
            do j = 2, 4
               cmp(j,i) = cmp(j,i) + uinp(j-1,i)
            end do
         end do
         call pme_cmp_to_fmp (cmp,fmp)
         call pme_grid_fmpole (fmp)
         call fftfront
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
               end do
            end do
         end do
         do i = 1, npole
            do j = 2, 4
               cmp(j,i) = cmp(j,i) + uind(j-1,i) - uinp(j-1,i)
            end do
         end do
         call pme_cmp_to_fmp (cmp,fmp)
         call pme_grid_fmpole (fmp)
         call fftfront
         do i = 1, npole
            do j = 2, 4
               cmp(j,i) = cmp(j,i) - uind(j-1,i)
            end do
         end do
      else
         call pme_cmp_to_fmp (cmp,fmp)
         call pme_grid_fmpole (fmp)
         call fftfront
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
               end do
            end do
         end do
      end if
c
c     make the scalar summation over reciprocal lattice
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
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
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm * hsq * bsmod1(k1) * bsmod2(k2)
     &                 * bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
            struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                  + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
            e = 0.5d0 * electric * expterm * struc2
            vterm = (2.0d0/hsq) * (1.0d0-term) * e
            vxx = vxx + h1*h1*vterm - e
            vyx = vyx + h2*h1*vterm
            vzx = vzx + h3*h1*vterm
            vyy = vyy + h2*h2*vterm - e
            vzy = vzy + h3*h2*vterm
            vzz = vzz + h3*h3*vterm - e
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     transform permanent multipoles without induced dipoles
c
      if (use_polar) then
         call pme_cmp_to_fmp (cmp,fmp)
         call pme_grid_fmpole (fmp)
         call fftfront
      end if
c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi * xboxi
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = 0.5d0 * expterm * struc2
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the PME grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call pme_get_fphi (fphi)
c
c     increment the permanent multipole energy and gradient
c
      e = 0.0d0
      do i = 1, npole
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 10
            e = e + fmp(k,i)*fphi(k,i)
            f1 = f1 + fmp(k,i)*fphi(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphi(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphi(deriv3(k),i)
         end do
         f1 = electric * dble(nfft1) * f1
         f2 = electric * dble(nfft2) * f2
         f3 = electric * dble(nfft3) * f3
         dfx = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         dfy = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         dfz = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         frc(1,i) = frc(1,i) + dfx
         frc(2,i) = frc(2,i) + dfy
         frc(3,i) = frc(3,i) + dfz
      end do
      e = 0.5d0 * electric * e
      em = em + e
      do i = 1, npole
         ii = ipole(i)
         dem(1,ii) = dem(1,ii) + frc(1,i)
         dem(2,ii) = dem(2,ii) + frc(2,i)
         dem(3,ii) = dem(3,ii) + frc(3,i)
      end do
c
c     distribute torques into the permanent multipole forces
c
      call pme_fphi_to_cphi (fphi,cphi)
      call build_frame (nframes,numfrlst,framelist,frdefpts,frame)
      call get_de_drot_mpole (cphi,cmp,de_drotsite)
      call accum_de_dframe_rot (nframes,frame_axis,de_drotsite,
     &                           frame,frdefpts,de_drotframe)
      call accum_de_ddefpts (nframes,frame_axis,de_drotframe,
     &                          frame,frdefpts,de_ddefpt)
      call trans_de_ddefpts_cent (numfrlst,framelist,frct,
     &                              nframes,de_ddefpt)
c
c     permanent torque contribution to the multipole derivatives
c
      do i = 1, n
         dem(1,i) = dem(1,i) - frct(1,i)
         dem(2,i) = dem(2,i) - frct(2,i)
         dem(3,i) = dem(3,i) - frct(3,i)
      end do
c
c     permanent torque contribution to the internal virial
c
      do i = 1, npole
         do j = 2, 4
            tfield(j) = 0.0d0
            do k = 2, 4
               tfield(j) = tfield(j) + ftc(j,k)*fphi(k,i)
            end do
            tfield(j) = electric * tfield(j)
         end do
         vxx = vxx - tfield(2)*cmp(2,i)
         vyx = vyx - 0.5d0*(tfield(2)*cmp(3,i)+tfield(3)*cmp(2,i))
         vzx = vzx - 0.5d0*(tfield(2)*cmp(4,i)+tfield(4)*cmp(2,i))
         vyy = vyy - tfield(3)*cmp(3,i)
         vzy = vzy - 0.5d0*(tfield(3)*cmp(4,i)+tfield(4)*cmp(3,i))
         vzz = vzz - tfield(4)*cmp(4,i)
         do j = 5, 10
            tfield(j) = 0.0d0
            do k = 5, 10
               tfield(j) = tfield(j) + ftc(j,k)*fphi(k,i)
            end do
            tfield(j) = electric * tfield(j)
         end do
         vxx = vxx - 2.0d0*tfield(5)*cmp(5,i) - tfield(8)*cmp(8,i)
     &                - tfield(9)*cmp(9,i)
         vyx = vyx - tfield(8)*(cmp(5,i)+cmp(6,i))
     &             - 0.5d0*(tfield(6)*cmp(8,i)+tfield(10)*cmp(9,i)
     &                     +tfield(5)*cmp(8,i)+tfield(9)*cmp(10,i))
         vzx = vzx - tfield(9)*(cmp(5,i)+cmp(7,i))
     &             - 0.5d0*(tfield(10)*cmp(8,i)+tfield(7)*cmp(9,i)
     &                     +tfield(5)*cmp(9,i)+tfield(8)*cmp(10,i))
         vyy = vyy - 2.0d0*tfield(6)*cmp(6,i) - tfield(8)*cmp(8,i)
     &                - tfield(10)*cmp(10,i)
         vzy = vzy - tfield(10)*(cmp(6,i)+cmp(7,i))
     &             - 0.5d0*(tfield(9)*cmp(8,i)+tfield(7)*cmp(10,i)
     &                     +tfield(8)*cmp(9,i)+tfield(6)*cmp(10,i))
         vzz = vzz - 2.0d0*tfield(7)*cmp(7,i) - tfield(9)*cmp(9,i)
     &                - tfield(10)*cmp(10,i)
      end do
c
c     convert Cartesian induced dipoles to fractional coordinates
c
      if (use_polar) then
         do i = 1, 3
            a(1,i) = dble(nfft1) * recip(i,1)
            a(2,i) = dble(nfft2) * recip(i,2)
            a(3,i) = dble(nfft3) * recip(i,3)
         end do
         do i = 1, npole
            do j = 1, 3
               fuind(j,i) = a(j,1)*uind(1,i) + a(j,2)*uind(2,i)
     &                          + a(j,3)*uind(3,i)
               fuinp(j,i) = a(j,1)*uinp(1,i) + a(j,2)*uinp(2,i)
     &                          + a(j,3)*uinp(3,i)
            end do
         end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
         call pme_grid_induce (fuind,fuinp)
         call fftfront
c
c     account for the zeroth grid point for a finite system
c
         if (.not. use_bounds) then
            expterm = 0.5d0 * pi * xboxi
            struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
            e = 0.5d0 * expterm * struc2
         end if
c
c     complete transformation of the PME grid
c
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  term = qfac(i,j,k)
                  qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
                  qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
               end do
            end do
         end do
c
c     perform 3-D FFT backward transform and get field
c
         call fftback
         call pme_get_induced (fdip_phi1,fdip_phi2,fdip_sum_phi)
c
c     convert the dipole fields from fractional to Cartesian
c
         do i = 1, 3
            a(i,1) = dble(nfft1) * recip(i,1)
            a(i,2) = dble(nfft2) * recip(i,2)
            a(i,3) = dble(nfft3) * recip(i,3)
         end do
         do i = 1, npole
            do j = 1, 3
               dipfield1(j,i) = a(j,1)*fdip_phi1(2,i)
     &                             + a(j,2)*fdip_phi1(3,i)
     &                             + a(j,3)*fdip_phi1(4,i)
               dipfield2(j,i) = a(j,1)*fdip_phi2(2,i)
     &                             + a(j,2)*fdip_phi2(3,i)
     &                             + a(j,3)*fdip_phi2(4,i)
            end do
         end do
c
c     increment the induced dipole energy and gradient
c
         e = 0.0d0
         do i = 1, npole
            f1 = 0.0d0
            f2 = 0.0d0
            f3 = 0.0d0
            do k = 1, 3
               j1 = deriv1(k+1)
               j2 = deriv2(k+1)
               j3 = deriv3(k+1)
               e = e + fuind(k,i)*fphi(k+1,i)
               f1 = f1 + (fuind(k,i)+fuinp(k,i))*fphi(j1,i)
     &                 + fuind(k,i)*fdip_phi2(j1,i)
     &                 + fuinp(k,i)*fdip_phi1(j1,i)
               f2 = f2 + (fuind(k,i)+fuinp(k,i))*fphi(j2,i)
     &                 + fuind(k,i)*fdip_phi2(j2,i)
     &                 + fuinp(k,i)*fdip_phi1(j2,i)
               f3 = f3 + (fuind(k,i)+fuinp(k,i))*fphi(j3,i)
     &                 + fuind(k,i)*fdip_phi2(j3,i)
     &                 + fuinp(k,i)*fdip_phi1(j3,i)
               if (poltyp .eq. 'DIRECT') then
                  f1 = f1 - fuind(k,i)*fdip_phi2(j1,i)
     &                    - fuinp(k,i)*fdip_phi1(j1,i)
                  f2 = f2 - fuind(k,i)*fdip_phi2(j2,i)
     &                    - fuinp(k,i)*fdip_phi1(j2,i)
                  f3 = f3 - fuind(k,i)*fdip_phi2(j3,i)
     &                    - fuinp(k,i)*fdip_phi1(j3,i)
               end if
            end do
            do k = 1, 10
               f1 = f1 + fmp(k,i)*fdip_sum_phi(deriv1(k),i)
               f2 = f2 + fmp(k,i)*fdip_sum_phi(deriv2(k),i)
               f3 = f3 + fmp(k,i)*fdip_sum_phi(deriv3(k),i)
            end do
            f1 = 0.5d0 * electric * dble(nfft1) * f1
            f2 = 0.5d0 * electric * dble(nfft2) * f2
            f3 = 0.5d0 * electric * dble(nfft3) * f3
            dfx = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
            dfy = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
            dfz = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
            frci(1,i) = frci(1,i) + dfx
            frci(2,i) = frci(2,i) + dfy
            frci(3,i) = frci(3,i) + dfz
         end do
         e = 0.5d0 * electric * e
         ep = ep + e
         do i = 1, npole
            ii = ipole(i)
            dep(1,ii) = dep(1,ii) + frci(1,i)
            dep(2,ii) = dep(2,ii) + frci(2,i)
            dep(3,ii) = dep(3,ii) + frci(3,i)
         end do
c
c     set the field to be the average induced dipole field
c
         do i = 1, npole
            do k = 1, 10
               fdip_sum_phi(k,i) = 0.5d0 * fdip_sum_phi(k,i)
            end do
         end do
c
c     distribute torques into the induced dipole forces
c
         call pme_fphi_to_cphi (fdip_sum_phi,cphi)
         call get_de_drot_mpole (cphi,cmp,de_drotsite)
         call accum_de_dframe_rot (nframes,frame_axis,de_drotsite,
     &                               frame,frdefpts,de_drotframe)
         call accum_de_ddefpts (nframes,frame_axis,de_drotframe,
     &                             frame,frdefpts,de_ddefpt)
         call trans_de_ddefpts_cent (numfrlst,framelist,frcit,
     &                                  nframes,de_ddefpt)
c
c     induced torque contribution to the polarization derivatives
c
         do i = 1, n
            dep(1,i) = dep(1,i) - frcit(1,i)
            dep(2,i) = dep(2,i) - frcit(2,i)
            dep(3,i) = dep(3,i) - frcit(3,i)
         end do
c
c     induced torque contribution to the internal virial
c
         do i = 1, npole
            do j = 2, 4
               tfield(j) = 0.0d0
               tdfield(j) = 0.0d0
               tpfield(j) = 0.0d0
               do k = 2, 4
                  tfield(j) = tfield(j) + ftc(j,k)*fphi(k,i)
                  tdfield(j) = tdfield(j) + ftc(j,k)*fdip_phi1(k,i)
                  tpfield(j) = tpfield(j) + ftc(j,k)*fdip_phi2(k,i)
               end do
               tfield(j) = electric * tfield(j)
               tdfield(j) = electric * tdfield(j)
               tpfield(j) = electric * tpfield(j)
            end do
            vxx = vxx - 0.5d0*((tdfield(2)+tpfield(2))*cmp(2,i)
     &                         +tfield(2)*(uind(1,i)+uinp(1,i))
     &                   +tdfield(2)*uinp(1,i)+tpfield(2)*uind(1,i))
            vyx = vyx - 0.25d0*((tdfield(2)+tpfield(2))*cmp(3,i)
     &                         +(tdfield(3)+tpfield(3))*cmp(2,i)
     &                          +tfield(2)*(uind(2,i)+uinp(2,i))
     &                          +tfield(3)*(uind(1,i)+uinp(1,i))
     &                   +tdfield(2)*uinp(2,i)+tpfield(2)*uind(2,i)
     &                   +tdfield(3)*uinp(1,i)+tpfield(3)*uind(1,i))
            vzx = vzx - 0.25d0*((tdfield(2)+tpfield(2))*cmp(4,i)
     &                         +(tdfield(4)+tpfield(4))*cmp(2,i)
     &                          +tfield(2)*(uind(3,i)+uinp(3,i))
     &                          +tfield(4)*(uind(1,i)+uinp(1,i))
     &                   +tdfield(2)*uinp(3,i)+tpfield(2)*uind(3,i)
     &                   +tdfield(4)*uinp(1,i)+tpfield(4)*uind(1,i))
            vyy = vyy - 0.5d0*((tdfield(3)+tpfield(3))*cmp(3,i)
     &                         +tfield(3)*(uind(2,i)+uinp(2,i))
     &                   +tdfield(3)*uinp(2,i)+tpfield(3)*uind(2,i))
            vzy = vzy - 0.25d0*((tdfield(3)+tpfield(3))*cmp(4,i)
     &                         +(tdfield(4)+tpfield(4))*cmp(3,i)
     &                          +tfield(3)*(uind(3,i)+uinp(3,i))
     &                          +tfield(4)*(uind(2,i)+uinp(2,i))
     &                   +tdfield(3)*uinp(3,i)+tpfield(3)*uind(3,i)
     &                   +tdfield(4)*uinp(2,i)+tpfield(4)*uind(2,i))
            vzz = vzz - 0.5d0*((tdfield(4)+tpfield(4))*cmp(4,i)
     &                         +tfield(4)*(uind(3,i)+uinp(3,i))
     &                   +tdfield(4)*uinp(3,i)+tpfield(4)*uind(3,i))
            do j = 5, 10
               tfield(j) = 0.0d0
               do k = 5, 10
                  tfield(j) = tfield(j) + ftc(j,k)*fdip_sum_phi(k,i)
               end do
               tfield(j) = electric * tfield(j)
            end do
            vxx = vxx - 2.0d0*tfield(5)*cmp(5,i) - tfield(8)*cmp(8,i)
     &                   - tfield(9)*cmp(9,i)
            vyx = vyx - tfield(8)*(cmp(5,i)+cmp(6,i))
     &                - 0.5d0*(tfield(6)*cmp(8,i)+tfield(10)*cmp(9,i)
     &                        +tfield(5)*cmp(8,i)+tfield(9)*cmp(10,i))
            vzx = vzx - tfield(9)*(cmp(5,i)+cmp(7,i))
     &                - 0.5d0*(tfield(10)*cmp(8,i)+tfield(7)*cmp(9,i)
     &                        +tfield(5)*cmp(9,i)+tfield(8)*cmp(10,i))
            vyy = vyy - 2.0d0*tfield(6)*cmp(6,i) - tfield(8)*cmp(8,i)
     &                   - tfield(10)*cmp(10,i)
            vzy = vzy - tfield(10)*(cmp(6,i)+cmp(7,i))
     &                - 0.5d0*(tfield(9)*cmp(8,i)+tfield(7)*cmp(10,i)
     &                        +tfield(8)*cmp(9,i)+tfield(6)*cmp(10,i))
            vzz = vzz - 2.0d0*tfield(7)*cmp(7,i) - tfield(9)*cmp(9,i)
     &                   - tfield(10)*cmp(10,i)
         end do
      end if
c
c     increment the internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vyx
      vir(3,1) = vir(3,1) + vzx
      vir(1,2) = vir(1,2) + vyx
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vzy
      vir(1,3) = vir(1,3) + vzx
      vir(2,3) = vir(2,3) + vzy
      vir(3,3) = vir(3,3) + vzz
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine build_frame  --  set local coordinate frames  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "build_frame" constructs the local coordinate frame axes at each
c     multipole site
c
c
      subroutine build_frame (nframes,numfrlst,framelist,frdefpt,frame)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,k,l,m,in
      integer k1,k2,k3
      integer numl,nframes
      integer numfr,numfrlst
      integer framelist(5,3*maxatm)
      real*8 wt,size,dot
      real*8 dx,dy,dz
      real*8 u(3),v(3),w(3)
      real*8 frdefpt(3,2,maxatm)
      real*8 frame(3,3,maxatm)
c
c
c     set local coordinate frame type at each multipole site
c
      numfr = 0
      numl = 0
      do i = 1, npole
         numfr = numfr + 1
         if (polaxe(i) .eq. 'Z-then-X') then
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = zaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 1
            framelist(5,numl) = 1
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = xaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 2
            framelist(5,numl) = 1
         else if (polaxe(i) .eq. 'Bisector') then
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = zaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 1
            framelist(5,numl) = 2
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = xaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 1
            framelist(5,numl) = 2
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = xaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 2
            framelist(5,numl) = 1
         else if (polaxe(i) .eq. 'Z-Bisect') then
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = zaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 1
            framelist(5,numl) = 1
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = xaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 2
            framelist(5,numl) = 2
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = yaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 2
            framelist(5,numl) = 2
         else if (polaxe(i) .eq. '3-Fold') then
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = zaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 1
            framelist(5,numl) = 3
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = xaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 1
            framelist(5,numl) = 3
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = yaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 1
            framelist(5,numl) = 3
            numl = numl + 1
            framelist(1,numl) = ipole(i)
            framelist(2,numl) = zaxis(i)
            framelist(3,numl) = numfr
            framelist(4,numl) = 2
            framelist(5,numl) = 1
         end if
      end do
      nframes = numfr
      numfrlst = numl
c
c     find frame defining points for each multipole site
c
      do k = 1, npole
         do j = 1, 2
            do i = 1, 3
               frdefpt(i,j,k) = 0.0d0
            end do
         end do
      end do
      do in = 1, numfrlst
         i = framelist(1,in)
         j = framelist(2,in)
         k = framelist(3,in)
         l = framelist(4,in)
         m = framelist(5,in)
         dx = x(j) - x(i)
         dy = y(j) - y(i)
         dz = z(j) - z(i)
         wt = m * sqrt(dx*dx+dy*dy+dz*dz)
         frdefpt(1,l,k) = frdefpt(1,l,k) + dx/wt
         frdefpt(2,l,k) = frdefpt(2,l,k) + dy/wt
         frdefpt(3,l,k) = frdefpt(3,l,k) + dz/wt
      end do
c
c     get the local axis vectors for each multipole site
c
      k1 = 3
      k2 = 1
      k3 = 2
      do in = 1, nframes
         do i = 1, 3
            u(i) = frdefpt(i,1,in)
         end do
         size = sqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))
         do i = 1, 3
            u(i) = u(i) / size
         end do
         do i = 1, 3
            v(i) = frdefpt(i,2,in)
         end do
         dot = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
         do i = 1, 3
            v(i) = v(i) - dot*u(i)
         end do
         size = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
         do i = 1, 3
            v(i) = v(i) / size
         end do
         w(1) = u(2)*v(3) - u(3)*v(2)
         w(2) = u(3)*v(1) - u(1)*v(3)
         w(3) = u(1)*v(2) - u(2)*v(1)
         do i = 1, 3
            frame(i,k1,in) = u(i)
            frame(i,k2,in) = v(i)
            frame(i,k3,in) = w(i)
         end do
      end do
      return
      end
c
c
c     ####################################
c     ##                                ##
c     ##  subroutine get_de_drot_mpole  ##
c     ##                                ##
c     ####################################
c
c
      subroutine get_de_drot_mpole (cphi,mpole_xyz,de_drotsite)
      implicit none
      include 'sizes.i'
      include 'chgpot.i'
      include 'mpole.i'
      integer i,j,k
      real*8 a(3,3)
      real*8 da(3,3)
      real*8 tmp_x(10)
      real*8 tmp_y(10)
      real*8 tmp_z(10)
      real*8 dmp_x(10,10)
      real*8 dmp_y(10,10)
      real*8 dmp_z(10,10)
      real*8 cphi(10,maxatm)
      real*8 mpole_xyz(10,maxatm)
      real*8 de_drotsite(3,maxatm)
c
c
c     to get de_drot we calculate the deriv of mpole with
c     respect to infinitesimal rotations about x, y and z
c
      do i = 1, 3
         do j = 1, 3
            a(i,j) = 0.0d0
         end do
         a(i,i) = 1.0d0
      end do
c
c     get x-axis rotation of dtheta
c
      do i = 1, 3
         do j = 1, 3
            da(i,j) = 0.0d0
         end do
      end do
      da(3,2) = 1.0d0
      da(2,3) = -1.0d0
      call xform_mpole_deriv (a,da,dmp_x)
c
c     get y-axis rotation of dtheta
c
      do i = 1, 3
         do j = 1, 3
            da(i,j) = 0.0d0
         end do
      end do
      da(3,1) = -1.0d0
      da(1,3) = 1.0d0
      call xform_mpole_deriv (a,da,dmp_y)
c
c     get z-axis rotation of dtheta
c
      do i = 1, 3
         do j = 1, 3
            da(i,j) = 0.0d0
         end do
      end do
      da(2,1) = 1.0d0
      da(1,2) = -1.0d0
      call xform_mpole_deriv (a,da,dmp_z)
c
c     set the values for each multipole site
c
      do i = 1, npole
         tmp_x(1) = dmp_x(1,1) * mpole_xyz(1,i)
         tmp_y(1) = dmp_y(1,1) * mpole_xyz(1,i)
         tmp_z(1) = dmp_z(1,1) * mpole_xyz(1,i)
         do j = 2, 4
            tmp_x(j) = 0.0d0
            tmp_y(j) = 0.0d0
            tmp_z(j) = 0.0d0
            do k = 2, 4
               tmp_x(j) = tmp_x(j) + dmp_x(j,k)*mpole_xyz(k,i)
               tmp_y(j) = tmp_y(j) + dmp_y(j,k)*mpole_xyz(k,i)
               tmp_z(j) = tmp_z(j) + dmp_z(j,k)*mpole_xyz(k,i)
            end do
         end do
         do j = 5, 10
            tmp_x(j) = 0.0d0
            tmp_y(j) = 0.0d0
            tmp_z(j) = 0.0d0
            do k = 5, 10
               tmp_x(j) = tmp_x(j) + dmp_x(j,k)*mpole_xyz(k,i)
               tmp_y(j) = tmp_y(j) + dmp_y(j,k)*mpole_xyz(k,i)
               tmp_z(j) = tmp_z(j) + dmp_z(j,k)*mpole_xyz(k,i)
            end do
         end do
         de_drotsite(1,i) = 0.0d0
         de_drotsite(2,i) = 0.0d0
         de_drotsite(3,i) = 0.0d0
         do j = 1, 10
            de_drotsite(1,i) = de_drotsite(1,i) + tmp_x(j)*cphi(j,i)
            de_drotsite(2,i) = de_drotsite(2,i) + tmp_y(j)*cphi(j,i)
            de_drotsite(3,i) = de_drotsite(3,i) + tmp_z(j)*cphi(j,i)
         end do
         de_drotsite(1,i) = electric * de_drotsite(1,i)
         de_drotsite(2,i) = electric * de_drotsite(2,i)
         de_drotsite(3,i) = electric * de_drotsite(3,i)
      end do
      return
      end
c
c
c     ####################################
c     ##                                ##
c     ##  subroutine xform_mpole_deriv  ##
c     ##                                ##
c     ####################################
c
c
      subroutine xform_mpole_deriv (a,da,dmp)
      implicit none
      integer i,j,k,m
      integer i1,i2
      integer qi1(6)
      integer qi2(6)
      real*8 a(3,3)
      real*8 da(3,3)
      real*8 dmp(10,10)
      data qi1  / 1, 2, 3, 1, 1, 2 /
      data qi2  / 1, 2, 3, 2, 3, 3 /
c
c
      do i = 1, 10
         do j = 1, 10
            dmp(j,i) = 0.0d0
         end do
      end do
      do i = 2, 4
         do j = 2, 4
            dmp(i,j) = da(i-1,j-1)
         end do
      end do
      do i1 = 1, 3
         k = qi1(i1)
         do i2 = 1, 6
            i = qi1(i2)
            j = qi2(i2)
            dmp(i1+4,i2+4) = da(k,i)*a(k,j) + a(k,i)*da(k,j)
         end do
      end do
      do i1 = 4, 6
         k = qi1(i1)
         m = qi2(i1)
         do i2 = 1, 6
            i = qi1(i2)
            j = qi2(i2)
            dmp(i1+4,i2+4) = da(k,i)*a(m,j) + da(k,j)*a(m,i)
     &                          + a(k,i)*da(m,j) + a(k,j)*da(m,i)
         end do
      end do
      return
      end
c
c
c     ######################################
c     ##                                  ##
c     ##  subroutine accum_de_dframe_rot  ##
c     ##                                  ##
c     ######################################
c
c
      subroutine accum_de_dframe_rot (nframes,frame_axis,de_drotsite,
     &                                 frame,frdefpt,de_drotframe)

      implicit none
      include 'sizes.i'
      include 'mpole.i'
      integer i,j,k
      integer k1,k2,k3
      integer nframes
      integer frame_axis(3)
      real*8 size
      real*8 p2unit(3)
      real*8 de_drotsite(3,npole)
      real*8 de_drotframe(3,nframes)
      real*8 frame(3,3,nframes)
      real*8 frdefpt(3,2,nframes)
c
c
      do i = 1, nframes
         do j = 1, 3
            de_drotframe(j,i) = 0.0d0
         end do
      end do
      k1 = frame_axis(1)
      k2 = frame_axis(2)
      k3 = frame_axis(3)
      k = 0
      do i = 1, npole
         k = k + 1
         size = sqrt(frdefpt(1,2,k)**2 + frdefpt(2,2,k)**2
     &                    + frdefpt(3,2,k)**2)
         do j = 1, 3
            p2unit(j) = frdefpt(j,2,k) / size
         end do
         de_drotframe(k1,k) = de_drotframe(k1,k)
     &                           + de_drotsite(1,i)*frame(1,k1,k)
     &                           + de_drotsite(2,i)*frame(2,k1,k)
     &                           + de_drotsite(3,i)*frame(3,k1,k)
         de_drotframe(k2,k) = de_drotframe(k2,k)
     &                           + de_drotsite(1,i)*p2unit(1)
     &                           + de_drotsite(2,i)*p2unit(2)
     &                           + de_drotsite(3,i)*p2unit(3)
         de_drotframe(k3,k) = de_drotframe(k3,k)
     &                           + de_drotsite(1,i)*frame(1,k3,k)
     &                           + de_drotsite(2,i)*frame(2,k3,k)
     &                           + de_drotsite(3,i)*frame(3,k3,k)
      end do
      return
      end
c
c
c     ###################################
c     ##                               ##
c     ##  subroutine accum_de_ddefpts  ##
c     ##                               ##
c     ###################################
c
c
      subroutine accum_de_ddefpts (nframes,frame_axis,de_drotframe,
     &                                 frame,defpts,de_ddefpt)
      implicit none
      integer i,j
      integer k1,k2,k3
      integer nframes
      integer frame_axis(3)
      real*8 dot12,dot21
      real*8 sizp1,sizp2
      real*8 sizp1perp2
      real*8 sizp2perp1
      real*8 dedv,dedw
      real*8 de_drotu
      real*8 de_drotv
      real*8 de_drotw
      real*8 u(3),v(3),w(3)
      real*8 p1(3),p2(3),p2unit(3)
      real*8 p2perp1(3),p1perp2(3)
      real*8 de_drotframe(3,nframes)
      real*8 defpts(3,2,nframes)
      real*8 de_ddefpt(3,2,nframes)
      real*8 frame(3,3,nframes)
c
c
      k1 = frame_axis(1)
      k2 = frame_axis(2)
      k3 = frame_axis(3)
      do i = 1, nframes
         do j = 1, 3
            p1(j) = defpts(j,1,i)
            p2(j) = defpts(j,2,i)
            u(j) = frame(j,k1,i)
            v(j) = frame(j,k2,i)
            w(j) = frame(j,k3,i)
         end do
         de_drotu = de_drotframe(k1,i)
         de_drotv = de_drotframe(k2,i)
         de_drotw = de_drotframe(k3,i)
         sizp1 = sqrt(p1(1)**2+p1(2)**2+p1(3)**2)
         sizp2 = sqrt(p2(1)**2+p2(2)**2+p2(3)**2)
         do j = 1, 3
            p2unit(j) = p2(j) / sizp2
         end do
         dot21 = u(1)*p2(1) + u(2)*p2(2) + u(3)*p2(3)
         dot12 = p1(1)*p2unit(1) + p1(2)*p2unit(2) + p1(3)*p2unit(3)
         do j = 1, 3
            p2perp1(j) = p2(j) - dot21*u(j)
            p1perp2(j) = p1(j) - dot12*p2unit(j)
         end do
         sizp2perp1 = sqrt(p2perp1(1)**2+p2perp1(2)**2+p2perp1(3)**2)
         sizp1perp2 = sqrt(p1perp2(1)**2+p1perp2(2)**2+p1perp2(3)**2)
         dedv = de_drotw / sizp1
         dedw = -de_drotv / sizp1perp2
         de_ddefpt(1,1,i) = dedv*v(1) + dedw*w(1)
         de_ddefpt(2,1,i) = dedv*v(2) + dedw*w(2)
         de_ddefpt(3,1,i) = dedv*v(3) + dedw*w(3)
         dedw = de_drotu / sizp2perp1
         de_ddefpt(1,2,i) = dedw*w(1)
         de_ddefpt(2,2,i) = dedw*w(2)
         de_ddefpt(3,2,i) = dedw*w(3)
      end do
      return
      end
c
c
c     ########################################
c     ##                                    ##
c     ##  subroutine trans_de_ddefpts_cent  ##
c     ##                                    ##
c     ########################################
c
c
      subroutine trans_de_ddefpts_cent (numfrlst,framelist,sitefrc,
     &                                      nframes,de_ddefpt)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,j,k,l,m
      integer in,nframes
      integer numfrlst
      integer framelist(5,3*maxatm)
      real*8 dx,dy,dz
      real*8 size,size2
      real*8 sizinv,sizinv3
      real*8 dedx,dedy,dedz
      real*8 dedux,deduy,deduz
      real*8 dux_dx,dux_dy,dux_dz
      real*8 duy_dy,duy_dz,duz_dz
      real*8 sitefrc(3,maxatm)
      real*8 de_ddefpt(3,2,nframes)
c
c
      do in = 1, numfrlst
         i = framelist(1,in)
         j = framelist(2,in)
         k = framelist(3,in)
         l = framelist(4,in)
         m = framelist(5,in)
         dx = x(j) - x(i)
         dy = y(j) - y(i)
         dz = z(j) - z(i)
         size2 = dx*dx + dy*dy + dz*dz
         size = sqrt(size2)
         sizinv = 1.0d0 / size
         sizinv3 = 1.0d0 / (size2*size)
         dux_dx = sizinv - dx*dx*sizinv3
         dux_dy = -dx * dy * sizinv3
         dux_dz = -dx * dz * sizinv3
         duy_dy = sizinv - dy*dy*sizinv3
         duy_dz = -dy * dz * sizinv3
         duz_dz = sizinv - dz*dz*sizinv3
         dedux = de_ddefpt(1,l,k) / dble(m)
         deduy = de_ddefpt(2,l,k) / dble(m)
         deduz = de_ddefpt(3,l,k) / dble(m)
         dedx = dedux*dux_dx + deduy*dux_dy + deduz*dux_dz
         dedy = dedux*dux_dy + deduy*duy_dy + deduz*duy_dz
         dedz = dedux*dux_dz + deduy*duy_dz + deduz*duz_dz
         sitefrc(1,i) = sitefrc(1,i) + dedx
         sitefrc(2,i) = sitefrc(2,i) + dedy
         sitefrc(3,i) = sitefrc(3,i) + dedz
         sitefrc(1,j) = sitefrc(1,j) - dedx
         sitefrc(2,j) = sitefrc(2,j) - dedy
         sitefrc(3,j) = sitefrc(3,j) - dedz
      end do
      return
      end

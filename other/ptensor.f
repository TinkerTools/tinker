c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder   ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine induce0a  --  induced dipoles via double loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "induce0a" computes the induced dipole moment at each
c     polarizable site using a pairwise double loop
c
c
      subroutine induce0a
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,m,skip(maxatm)
      integer ii,iz,ix,kk,kz,kx
      integer iter,maxiter
      real*8 xr,yr,zr,eps,taper
      real*8 r,r2,r3,r4,r5
      real*8 rpi(13),rpk(13)
      real*8 udir(3,maxatm)
      real*8 uold(3,maxatm)
      real*8 fieldi(3),fieldk(3)
      logical iuse,kuse
c
c
c     zero out induced dipoles and count the polarizable sites
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the switching function coefficients
c
      call switch ('MPOLE')
c
c     compute the direct induced dipole moment at each atom
c
      do ii = 1, npole-1
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * chg12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * chg13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * chg14use
         end do
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do kk = ii+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
            if (iuse .or. kuse) then
               if (skip(k) .ne. i) then
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  if (use_image)  call image (xr,yr,zr,0)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     do j = 1, maxpole
                        rpk(j) = rpole(j,kk)
                     end do
                     call udirect (ii,kk,xr,yr,zr,r,r2,rpi,
     &                                rpk,fieldi,fieldk)
                     if (skip(k) .eq. -i) then
                        do j = 1, 3
                           fieldi(j) = fieldi(j) / chgscale
                           fieldk(j) = fieldk(j) / chgscale
                        end do
                     end if
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        do j = 1, 3
                           fieldi(j) = fieldi(j) * taper
                           fieldk(j) = fieldk(j) * taper
                        end do
                     end if
                     do j = 1, 3
                        uind(j,ii) = uind(j,ii) + fieldi(j)
                        uind(j,kk) = uind(j,kk) + fieldk(j)
                     end do
                  end if
               end if
            end if
         end do
      end do
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do ii = 1, npole
            i = ipole(ii)
            iz = zaxis(ii)
            ix = xaxis(ii)
            iuse = (use(i) .or. use(iz) .or. use(ix))
            do j = 1, maxpole
               rpi(j) = rpole(j,ii)
            end do
            do kk = ii, npole
               k = ipole(kk)
               kz = zaxis(kk)
               kx = xaxis(kk)
               kuse = (use(k) .or. use(kz) .or. use(kx))
               if (iuse .or. kuse) then
                  do m = 1, ncell
                     xr = x(k) - x(i)
                     yr = y(k) - y(i)
                     zr = z(k) - z(i)
                     call image (xr,yr,zr,m)
                     r2 = xr*xr + yr*yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        do j = 1, maxpole
                           rpk(j) = rpole(j,kk)
                        end do
                        call udirect (ii,kk,xr,yr,zr,r,r2,rpi,
     &                                   rpk,fieldi,fieldk)
                        if (r2 .gt. cut2) then
                           r3 = r2 * r
                           r4 = r2 * r2
                           r5 = r2 * r3
                           taper = c5*r5 + c4*r4 + c3*r3
     &                                + c2*r2 + c1*r + c0
                           do j = 1, 3
                              fieldi(j) = fieldi(j) * taper
                              fieldk(j) = fieldk(j) * taper
                           end do
                        end if
                        do j = 1, 3
                           uind(j,ii) = uind(j,ii) + fieldi(j)
                           if (ii .ne. kk)
     &                        uind(j,kk) = uind(j,kk) + fieldk(j)
                        end do
                     end if
                  end do
               end if
            end do
         end do
      end if
c
c     convert total field at each atom to induced dipole moment
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = polarize(i) * uind(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         maxiter = 100
         iter = 0
         eps = 1.0d0
         do i = 1, npole
            do j = 1, 3
               udir(j,i) = uind(j,i)
            end do
         end do
c
c     compute mutual induced dipole moments by an iterative method
c
         dowhile (iter.lt.maxiter .and. eps.gt.poleps)
            do ii = 1, npole
               do j = 1, 3
                  uold(j,ii) = uind(j,ii)
                  uind(j,ii) = 0.0d0
               end do
            end do
            do ii = 1, npole-1
               i = ipole(ii)
               iz = zaxis(ii)
               ix = xaxis(ii)
               iuse = (use(i) .or. use(iz) .or. use(ix))
               do j = 1, n12(i)
                  skip(i12(j,i)) = i * chg12use
               end do
               do j = 1, n13(i)
                  skip(i13(j,i)) = i * chg13use
               end do
               do j = 1, n14(i)
                  skip(i14(j,i)) = i * chg14use
               end do
               do kk = ii+1, npole
                  k = ipole(kk)
                  kz = zaxis(kk)
                  kx = xaxis(kk)
                  kuse = (use(k) .or. use(kz) .or. use(kx))
                  if (iuse .or. kuse) then
                     if (skip(k) .ne. i) then
                        xr = x(k) - x(i)
                        yr = y(k) - y(i)
                        zr = z(k) - z(i)
                        if (use_image)  call image (xr,yr,zr,0)
                        r2 = xr*xr + yr*yr + zr*zr
                        if (r2 .le. off2) then
                           r = sqrt(r2)
                           call umutual (ii,kk,xr,yr,zr,r,r2,uold(1,ii),
     &                                      uold(1,kk),fieldi,fieldk)
                           if (skip(k) .eq. -i) then
                              do j = 1, 3
                                 fieldi(j) = fieldi(j) / chgscale
                                 fieldk(j) = fieldk(j) / chgscale
                              end do
                           end if
                           if (r2 .gt. cut2) then
                              r3 = r2 * r
                              r4 = r2 * r2
                              r5 = r2 * r3
                              taper = c5*r5 + c4*r4 + c3*r3
     &                                   + c2*r2 + c1*r + c0
                              do j = 1, 3
                                 fieldi(j) = fieldi(j) * taper
                                 fieldk(j) = fieldk(j) * taper
                              end do
                           end if
                           do j = 1, 3
                              uind(j,ii) = uind(j,ii) + fieldi(j)
                              uind(j,kk) = uind(j,kk) + fieldk(j)
                           end do
                        end if
                     end if
                  end if
               end do
            end do
c
c     periodic boundary for large cutoffs via replicates method
c
            if (use_replica) then
               do ii = 1, npole
                  i = ipole(ii)
                  iz = zaxis(ii)
                  ix = xaxis(ii)
                  iuse = (use(i) .or. use(iz) .or. use(ix))
                  do kk = ii, npole
                     k = ipole(kk)
                     kz = zaxis(kk)
                     kx = xaxis(kk)
                     kuse = (use(k) .or. use(kz) .or. use(kx))
                     if (iuse .or. kuse) then
                        do m = 1, ncell
                           xr = x(k) - x(i)
                           yr = y(k) - y(i)
                           zr = z(k) - z(i)
                           call image (xr,yr,zr,m)
                           r2 = xr*xr + yr*yr + zr*zr
                           if (r2 .le. off2) then
                              r = sqrt(r2)
                              call umutual (ii,kk,xr,yr,zr,r,r2,
     &                                      uold(1,ii),uold(1,kk),
     &                                         fieldi,fieldk)
                              if (r2 .gt. cut2) then
                                 r3 = r2 * r
                                 r4 = r2 * r2
                                 r5 = r2 * r3
                                 taper = c5*r5 + c4*r4 + c3*r3
     &                                      + c2*r2 + c1*r + c0
                                 do j = 1, 3
                                    fieldi(j) = fieldi(j) * taper
                                    fieldk(j) = fieldk(j) * taper
                                 end do
                              end if
                              do j = 1, 3
                                 uind(j,ii) = uind(j,ii) + fieldi(j)
                                 if (ii .ne. kk)
     &                              uind(j,kk) = uind(j,kk) + fieldk(j)
                              end do
                           end if
                        end do
                     end if
                  end do
               end do
            end if
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            eps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = udir(j,i) + polarize(i)*uind(j,i)
                  eps = eps + (uind(j,i)-uold(j,i))**2
               end do
            end do
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
         end do
c
c     terminate the calculation if dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,30)
   30       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine udirect  --  field for direct induced dipoles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "udirect" evaluates the electric field at a polarizable atom
c     due to permanent atomic multipoles at a second atom, and vice
c     versa, for use in computation of direct induced dipole moments
c
c
      subroutine udirect (ii,kk,xr,yr,zr,r,r2,rpi,rpk,fieldi,fieldk)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      integer i,ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,damp
      real*8 rr3,rr5,rr7
      real*8 xr5,yr5,xyz
      real*8 xr7,yr7,rr53
      real*8 t2,t3,t4,t5,t6,t7,t8,t9,t10,t11
      real*8 t12,t13,t14,t15,t16,t17,t18,t19
      real*8 i0,i1,i2,i3,i4,i5,i6,i7,i8,i9
      real*8 k0,k1,k2,k3,k4,k5,k6,k7,k8,k9
      real*8 i9i4,i9i7,k4k9,k7k9
      real*8 rpi(13),rpk(13)
      real*8 fieldi(3),fieldk(3)
      logical iquad,kquad
c
c
c     set extent of the multipole expansion at each site
c
      iquad = (polsiz(ii) .ge. 13)
      kquad = (polsiz(kk) .ge. 13)
c
c     compute the first order T2 matrix elements
c
      rr3 = -1.0d0 / (r2*r)
      t2 = xr * rr3
      t3 = yr * rr3
      t4 = zr * rr3
c
c     compute the second order T2 matrix elements
c
      rr5 = -3.0d0 * rr3 / r2
      xr5 = xr * rr5
      yr5 = yr * rr5
      t5 = rr3 + xr5*xr
      t6 = xr5 * yr
      t7 = xr5 * zr
      t8 = rr3 + yr5*yr
      t9 = yr5 * zr
      t10 = -t5 - t8
c
c     compute the third order T2 matrix elements
c
      if (iquad .or. kquad) then
         rr7 = -5.0d0 * rr5 / r2
         xyz = xr * yr * zr
         xr7 = xr * xr * rr7
         yr7 = yr * yr * rr7
         rr53 = 3.0d0 * rr5
         t11 = xr * (xr7+rr53)
         t12 = yr * (xr7+rr5)
         t13 = zr * (xr7+rr5)
         t14 = xr * (yr7+rr5)
         t15 = xyz * rr7
         t16 = -t11 - t14
         t17 = yr * (yr7+rr53)
         t18 = zr * (yr7+rr5)
         t19 = -t12 - t17
      end if
c
c     compute the field at site k due to multipoles at site i
c
      i0 = rpi(1)
      i1 = rpi(2)
      i2 = rpi(3)
      i3 = rpi(4)
      fieldk(1) = -i0*t2 + i1*t5 + i2*t6 + i3*t7
      fieldk(2) = -i0*t3 + i1*t6 + i2*t8 + i3*t9
      fieldk(3) = -i0*t4 + i1*t7 + i2*t9 + i3*t10
      if (iquad) then
         i4 = rpi(5)
         i5 = rpi(6)
         i6 = rpi(7)
         i7 = rpi(9)
         i8 = rpi(10)
         i9 = rpi(13)
         i9i4 = i9 - i4
         i9i7 = i9 - i7
         fieldk(1) = fieldk(1) + i9i4*t11 + i9i7*t14
     &                  - 2.0d0*(i5*t12 + i6*t13 + i8*t15)
         fieldk(2) = fieldk(2) + i9i4*t12 + i9i7*t17
     &                  - 2.0d0*(i5*t14 + i6*t15 + i8*t18)
         fieldk(3) = fieldk(3) + i9i4*t13 + i9i7*t18
     &                  - 2.0d0*(i5*t15 + i6*t16 + i8*t19)
      end if
c
c     compute the field at site i due to multipoles at site k
c
      k0 = rpk(1)
      k1 = rpk(2)
      k2 = rpk(3)
      k3 = rpk(4)
      fieldi(1) = k0*t2 + k1*t5 + k2*t6 + k3*t7
      fieldi(2) = k0*t3 + k1*t6 + k2*t8 + k3*t9
      fieldi(3) = k0*t4 + k1*t7 + k2*t9 + k3*t10
      if (kquad) then
         k4 = rpk(5)
         k5 = rpk(6)
         k6 = rpk(7)
         k7 = rpk(9)
         k8 = rpk(10)
         k9 = rpk(13)
         k4k9 = k4 - k9
         k7k9 = k7 - k9
         fieldi(1) = fieldi(1) + k4k9*t11 + k7k9*t14
     &                  + 2.0d0*(k5*t12 + k6*t13 + k8*t15)
         fieldi(2) = fieldi(2) + k4k9*t12 + k7k9*t17
     &                  + 2.0d0*(k5*t14 + k6*t15 + k8*t18)
         fieldi(3) = fieldi(3) + k4k9*t13 + k7k9*t18
     &                  + 2.0d0*(k5*t15 + k6*t16 + k8*t19)
      end if
c
c     apply a damping factor to reduce the field at short range
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         damp = -pgamma * (r/damp)**3
         if (damp .gt. -50.0d0) then
            damp = 1.0d0 - exp(damp)
            do i = 1, 3
               fieldi(i) = fieldi(i) * damp
               fieldk(i) = fieldk(i) * damp
            end do
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine umutual  --  field for mutual induced dipoles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "umutual" evaluates the electric field at a polarizable atom
c     due to the induced atomic dipoles at a second atom, and vice
c     versa, for use in computation of mutual induced dipole moments
c
c
      subroutine umutual (ii,kk,xr,yr,zr,r,r2,upi,upk,fieldi,fieldk)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      integer i,ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,damp
      real*8 rr3,rr5,xr5,yr5
      real*8 upi(3),upk(3)
      real*8 fieldi(3),fieldk(3)
      real*8 u1,u2,u3,v1,v2,v3
      real*8 t5,t6,t7,t8,t9,t10
c
c
c     set the current values of the induced dipoles
c
      u1 = upi(1)
      u2 = upi(2)
      u3 = upi(3)
      v1 = upk(1)
      v2 = upk(2)
      v3 = upk(3)
c
c     compute the second order T2 matrix elements
c
      rr3 = -1.0d0 / (r2*r)
      rr5 = -3.0d0 * rr3 / r2
      xr5 = xr * rr5
      yr5 = yr * rr5
      t5 = rr3 + xr5*xr
      t6 = xr5 * yr
      t7 = xr5 * zr
      t8 = rr3 + yr5*yr
      t9 = yr5 * zr
      t10 = -t5 - t8
c
c     compute the field at site i due to site k and vice versa
c
      fieldi(1) = v1*t5 + v2*t6 + v3*t7
      fieldi(2) = v1*t6 + v2*t8 + v3*t9
      fieldi(3) = v1*t7 + v2*t9 + v3*t10
      fieldk(1) = u1*t5 + u2*t6 + u3*t7
      fieldk(2) = u1*t6 + u2*t8 + u3*t9
      fieldk(3) = u1*t7 + u2*t9 + u3*t10
c
c     apply a damping factor to reduce the field at short range
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         damp = -pgamma * (r/damp)**3
         if (damp .gt. -50.0d0) then
            damp = 1.0d0 - exp(damp)
            do i = 1, 3
               fieldi(i) = fieldi(i) * damp
               fieldk(i) = fieldk(i) * damp
            end do
         end if
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole0a  --  double loop multipole energy  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole0a" calculates the electrostatic energy due to
c     atomic multipole interactions and dipole polarizability
c     using a pairwise double loop
c
c
      subroutine empole0a
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'mpole.i'
      include 'polar.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,m
      integer ii,iz,ix
      integer kk,kz,kx
      integer skip(maxatm)
      real*8 eik,ei,ek,a(3,3)
      real*8 shift,taper,trans
      real*8 f,fik,fgrp
      real*8 xr,yr,zr,r
      real*8 r2,r3,r4,r5,r6,r7
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      logical proceed,iuse,kuse
c
c
c     zero out the atomic multipole interaction energy
c
      em = 0.0d0
      ep = 0.0d0
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and cutoff values
c
      f = electric / dielec
      call switch ('MPOLE')
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each atom
c
      call induce
c
c     calculate the multipole interaction energy term
c
      do ii = 1, npole-1
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * chg12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * chg13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * chg14use
         end do
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. kuse)
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  do j = 1, maxpole
                     rpk(j) = rpole(j,kk)
                  end do
                  do j = 1, 3
                     indk(j) = uind(j,kk)
                  end do
                  r = sqrt(r2)
                  call empik (ii,kk,xr,yr,zr,r,rpi,rpk,
     &                          indi,indk,eik,ei,ek)
                  fik = f
                  if (skip(k) .eq. -i)  fik = fik / chgscale
                  eik = fik * eik
                  ei = fik * ei
                  ek = fik * ek
c
c     use shifted energy switching if near the cutoff distance
c
                  fik = fik * rpi(1) * rpk(1)
                  shift = fik / (0.5d0*(off+cut))
                  eik = eik - shift
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     eik = eik * taper + trans
                     ei = ei * taper
                     ek = ek * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     eik = eik * fgrp
                     ei = ei * fgrp
                     ek = ek * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  em = em + eik
                  ep = ep + ei + ek
               end if
            end if
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
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. kuse)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do m = 1, ncell
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  call image (xr,yr,zr,m)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     do j = 1, maxpole
                        rpk(j) = rpole(j,kk)
                     end do
                     do j = 1, 3
                        indk(j) = uind(j,kk)
                     end do
                     r = sqrt(r2)
                     call empik (ii,kk,xr,yr,zr,r,rpi,rpk,
     &                             indi,indk,eik,ei,ek)
                     fik = f
                     eik = fik * eik
                     ei = fik * ei
                     ek = fik * ek
c
c     use shifted energy switching if near the cutoff distance
c
                     fik = fik * rpi(1) * rpk(1)
                     shift = fik / (0.5d0*(off+cut))
                     eik = eik - shift
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        r6 = r3 * r3
                        r7 = r3 * r4
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                  + f3*r3 + f2*r2 + f1*r + f0)
                        eik = eik * taper + trans
                        ei = ei * taper
                        ek = ek * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        eik = eik * fgrp
                        ei = ei * fgrp
                        ek = ek * fgrp
                     end if
c
c     increment the overall multipole and polarization energies
c
                     if (i .eq. k) then
                        eik = 0.5d0 * eik
                        ei = 0.5d0 * ei
                        ek = 0.5d0 * ek
                     end if
                     em = em + eik
                     ep = ep + ei + ek
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empik  --  multipole & polarization pair energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empik" computes the permanent multipole and induced dipole
c     energies between a specified pair of atomic multipole sites
c
c
      subroutine empik (ii,kk,xr,yr,zr,r,rpi,rpk,indi,indk,eik,ei,ek)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      integer ii,kk
      real*8 eik,ei,ek
      real*8 xr,yr,zr,r,damp
      real*8 rr1,rr2,rr3,rr5,rr7,rr9
      real*8 xx,yy,xyz,xry,yrx,zrx,zry
      real*8 xr5,yr5,xr7,yr7
      real*8 xr9,xyr9,yr9,xr9x7,yr9y7
      real*8 rr53,xr73,yr73
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9,t11
      real*8 t12,t13,t14,t15,t17,t18,t21,t22
      real*8 t23,t24,t25,t27,t28,t31,t32
      real*8 i0,i1,i2,i3,i4,i5,i6,i7,i8,i9
      real*8 k0,k1,k2,k3,k4,k5,k6,k7,k8,k9
      real*8 u1,u2,u3,v1,v2,v3
      real*8 i3k3,i4i9,i7i9,k4k9,k7k9,i3k6
      real*8 i3k8,i6k3,i8k3,ik49,ik59,ik69
      real*8 ik79,ik89,ik68,i6k6,i8k8,i9k9
      real*8 i3v3,i6v3,i8v3,k3u3,k6u3,k8u3
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      logical iquad,kquad
c
c
c     zero out multipole and polarization energy components
c
      eik = 0.0d0
      ei = 0.0d0
      ek = 0.0d0
c
c     check for presence of quadrupole components at either site
c
      iquad = (polsiz(ii) .ge. 13)
      kquad = (polsiz(kk) .ge. 13)
c
c     set permanent and induced multipoles for first site
c
      i0 = rpi(1)
      i1 = rpi(2)
      i2 = rpi(3)
      i3 = rpi(4)
      i4 = rpi(5)
      i5 = rpi(6)
      i6 = rpi(7)
      i7 = rpi(9)
      i8 = rpi(10)
      i9 = rpi(13)
      u1 = indi(1)
      u2 = indi(2)
      u3 = indi(3)
c
c     set permanent and induced multipoles for second site
c
      k0 = rpk(1)
      k1 = rpk(2)
      k2 = rpk(3)
      k3 = rpk(4)
      k4 = rpk(5)
      k5 = rpk(6)
      k6 = rpk(7)
      k7 = rpk(9)
      k8 = rpk(10)
      k9 = rpk(13)
      v1 = indk(1)
      v2 = indk(2)
      v3 = indk(3)
c
c     compute the zeroth order T2 matrix element
c
      rr1 = 1.0d0 / r
      t1 = rr1
c
c     compute the first order T2 matrix elements
c
      rr2 = rr1 * rr1
      rr3 = rr1 * rr2
      t2 = -xr * rr3
      t3 = -yr * rr3
      t4 = -zr * rr3
c
c     compute the second order T2 matrix elements
c
      rr5 = 3.0d0 * rr3 * rr2
      xr5 = xr * rr5
      yr5 = yr * rr5
      t5 = -rr3 + xr5*xr
      t6 = xr5 * yr
      t7 = xr5 * zr
      t8 = -rr3 + yr5*yr
      t9 = yr5 * zr
c
c     compute the third order T2 matrix elements
c
      if (iquad .or. kquad) then
         rr7 = 5.0d0 * rr5 * rr2
         xx = xr * xr
         yy = yr * yr
         xyz = xr * yr * zr
         xr7 = xx * rr7
         yr7 = yy * rr7
         rr53 = 3.0d0 * rr5
         t11 = xr * (rr53-xr7)
         t12 = yr * (rr5-xr7)
         t13 = zr * (rr5-xr7)
         t14 = xr * (rr5-yr7)
         t15 = -xyz * rr7
         t17 = yr * (rr53-yr7)
         t18 = zr * (rr5-yr7)
c
c     compute the fourth order T2 matrix elements
c
         rr9 = 7.0d0 * rr7 * rr2
         if (xr .eq. 0.0d0) then
            yrx = 0.0d0
            zrx = 0.0d0
         else
            yrx = yr / xr
            zrx = zr / xr
         end if
         if (yr .eq. 0.0d0) then
            xry = 0.0d0
            zry = 0.0d0
         else
            xry = xr / yr
            zry = zr / yr
         end if
         xr9 = xx * xx * rr9
         xyr9 = xx * yy * rr9
         yr9 = yy * yy * rr9
         xr73 = 3.0d0 * xr7
         yr73 = 3.0d0 * yr7
         xr9x7 = xr9 - xr73
         yr9y7 = yr9 - yr73
         t21 = xr9x7 - xr73 + rr53
         t22 = yrx * xr9x7
         t23 = zrx * xr9x7
         t24 = xyr9 - xr7 - yr7 + rr5
         t25 = zry * (xyr9-yr7)
         t27 = xry * yr9y7
         t28 = zrx * (xyr9-xr7)
         t31 = yr9y7 - yr73 + rr53
         t32 = zry * yr9y7
      end if
c
c     get the M-M, M-D and D-D parts of the multipole energy
c
      i3k3 = i3 * k3
      eik = i0*k0*t1 + (i0*k1-i1*k0)*t2 + (i0*k2-i2*k0)*t3
     &         + (i0*k3-i3*k0)*t4 + (i3k3-i1*k1)*t5 - (i1*k2+i2*k1)*t6
     &         - (i1*k3+i3*k1)*t7 + (i3k3-i2*k2)*t8 - (i2*k3+i3*k2)*t9
c
c     get the M-indD and D-indD parts of the polarization energy
c
      i3v3 = i3 * v3
      k3u3 = k3 * u3
      ei = -k0*(u1*t2+u2*t3+u3*t4) + (k3u3-k1*u1)*t5 - (k1*u2+k2*u1)*t6
     &        - (k3*u1+k1*u3)*t7 + (k3u3-k2*u2)*t8 - (k2*u3+k3*u2)*t9
      ek = i0*(v1*t2+v2*t3+v3*t4) + (i3v3-i1*v1)*t5 - (i1*v2+i2*v1)*t6
     &        - (i3*v1+i1*v3)*t7 + (i3v3-i2*v2)*t8 - (i2*v3+i3*v2)*t9
c
c     get the M-Q, D-Q and Q-Dind energies, if necessary
c
      if (kquad) then
         k4k9 = k4 - k9
         k7k9 = k7 - k9
         i3k6 = 2.0d0 * i3 * k6
         i3k8 = 2.0d0 * i3 * k8
         eik = eik + i0*k4k9*t5 + 2.0d0*i0*k5*t6 + 2.0d0*i0*k6*t7
     &             + i0*k7k9*t8 + 2.0d0*i0*k8*t9 + (i3k6-i1*k4k9)*t11
     &             + (i3k8-i2*k4k9-2.0d0*i1*k5)*t12
     &             - (2.0d0*i1*k6+i3*k4k9)*t13
     &             + (i3k6-i1*k7k9-2.0d0*i2*k5)*t14
     &             - 2.0d0*(i1*k8+i2*k6+i3*k5)*t15
     &             + (i3k8-i2*k7k9)*t17
     &             - (2.0d0*i2*k8+i3*k7k9)*t18
         k6u3 = k6 * u3
         k8u3 = k8 * u3
         ei = ei + (-k4k9*u1+2.0d0*k6u3)*t11
     &           + (-k4k9*u2+2.0d0*(k8u3-k5*u1))*t12
     &           + (-k4k9*u3-2.0d0*k6*u1)*t13
     &           + (-k7k9*u1+2.0d0*(k6u3-k5*u2))*t14
     &           - 2.0d0*(k5*u3+k6*u2+k8*u1)*t15
     &           + (-k7k9*u2+2.0d0*k8u3)*t17
     &           + (-k7k9*u3-2.0d0*k8*u2)*t18
      end if
      if (iquad) then
         i4i9 = i4 - i9
         i7i9 = i7 - i9
         i6k3 = -2.0d0 * i6 * k3
         i8k3 = -2.0d0 * i8 * k3
         eik = eik + i4i9*k0*t5 + 2.0d0*i5*k0*t6 + 2.0d0*i6*k0*t7
     &             + i7i9*k0*t8 + 2.0d0*i8*k0*t9 + (i4i9*k1+i6k3)*t11
     &             + (i8k3+i4i9*k2+2.0d0*i5*k1)*t12
     &             + (2.0d0*i6*k1+i4i9*k3)*t13
     &             + (i6k3+i7i9*k1+2.0d0*i5*k2)*t14
     &             + 2.0d0*(i8*k1+i6*k2+i5*k3)*t15
     &             + (i8k3+i7i9*k2)*t17
     &             + (2.0d0*i8*k2+i7i9*k3)*t18
         i6v3 = i6 * v3
         i8v3 = i8 * v3
         ek = ek + (i4i9*v1-2.0d0*i6v3)*t11
     &           + (i4i9*v2+2.0d0*(i5*v1-i8v3))*t12
     &           + (i4i9*v3+2.0d0*i6*v1)*t13
     &           + (i7i9*v1+2.0d0*(i5*v2-i6v3))*t14
     &           + 2.0d0*(i5*v3+i6*v2+i8*v1)*t15
     &           + (i7i9*v2-2.0d0*i8v3)*t17
     &           + (i7i9*v3+2.0d0*i8*v2)*t18
      end if
c
c     get the Q-Q part of the multipole interaction energy
c
      if (iquad .and. kquad) then
         ik49 = i4*k9 + i9*k4
         ik59 = i5*k9 + i9*k5
         ik69 = i6*k9 + i9*k6
         ik79 = i7*k9 + i9*k7
         ik89 = i8*k9 + i9*k8
         ik68 = 2.0d0 * (i6*k8 + i8*k6)
         i6k6 = i6 * k6
         i8k8 = i8 * k8
         i9k9 = i9 * k9
         eik = eik + (i4*k4-ik49+i9k9-4.0d0*i6k6)*t21
     &             + 2.0d0*(i5*k4+i4*k5-ik68-ik59)*t22
     &             + 2.0d0*(i4*k6+i6*k4-ik69)*t23
     &             + (i4*k7+i7*k4-ik49-ik79+2.0d0*i9k9
     &                  +4.0d0*(i5*k5-i6k6-i8k8))*t24
     &             + 2.0d0*(i4*k8+i8*k4+2.0d0*(i5*k6+i6*k5)-ik89)*t25
     &             + 2.0d0*(i5*k7+i7*k5-ik68-ik59)*t27
     &             + 2.0d0*(i6*k7+i7*k6+2.0d0*(i5*k8+i8*k5)-ik69)*t28
     &             + (i7*k7-ik79+i9k9-4.0d0*i8k8)*t31
     &             + 2.0d0*(i7*k8+i8*k7-ik89)*t32
      end if
c
c     apply exponential damping to the polarization term
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         damp = -pgamma * (r/damp)**3
         if (damp .gt. -50.0d0) then
            damp = 1.0d0 - exp(damp)
            ei = ei * damp
            ek = ek * damp
         end if
      end if
c
c     final polarization energy is half of value computed above
c
      ei = 0.5d0 * ei
      ek = 0.5d0 * ek
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empik1  --  mpole & polarization pair gradient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empik1" computes the permanent multipole and induced dipole
c     energies and derivatives between a pair of multipole sites
c
c
      subroutine empik1 (ii,kk,xr,yr,zr,r,r2,rpi,rpk,indi,indk,
     &                    eik,ei,ek,dm,dmi,dmk,dp,dpi,dpk,utu)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      integer i,j,ii,kk
      real*8 eik,ei,ek
      real*8 damp,ddamp,de,term
      real*8 xr,yr,zr,r,r2,r4
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      real*8 dm(3),dp(3),utu
      real*8 dmi(3,3),dmk(3,3)
      real*8 dpi(3,3),dpk(3,3)
      real*8 rr1,rr2,rr3,rr5,rr7,rr9,rr11
      real*8 xx,yy,xyz,xry,yrx,zrx,zry
      real*8 xr5,yr5,xr7,yr7
      real*8 xr9,xyr9,yr9,xr9x7,yr9y7
      real*8 rr53,xr73,yr73
      real*8 x2,y2,z2,x4,y4,x2y2,x2r2,y2r2
      real*8 xrr11,yrr11,zrr11,xyzrr11
      real*8 r445,r4225,nnn1,nnn2
      real*8 n945x4,n945y4,n945x2y2
      real*8 i3k3,i6k6,i8k8,i9k9
      real*8 ik49,ik59,ik69,ik79,ik89,ik68
      real*8 k4k9,k7k9,i4i9,i7i9
      real*8 k3u3,i3k6,i3k8,k6u3,k8u3
      real*8 i3v3,i6k3,i8k3,i6v3,i8v3
      real*8 i0,i1,i2,i3,i4,i5,i6,i7,i8,i9
      real*8 k0,k1,k2,k3,k4,k5,k6,k7,k8,k9
      real*8 u1,u2,u3,v1,v2,v3
      real*8 di1,di2,di3,di4,di5,di6,di7,di8,di9
      real*8 dk1,dk2,dk3,dk4,dk5,dk6,dk7,dk8,dk9
      real*8 w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12
      real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13
      real*8 t14,t15,t16,t17,t18,t19,t21,t22,t23,t24,t25
      real*8 t26,t27,t28,t29,t31,t32,t33,t36,t37,t38,t39
      real*8 t40,t41,t42,t43,t44,t46,t47,t48,t51,t52,t53
      logical iquad,kquad
c
c
c     zero out the energy and first derivative components
c
      eik = 0.0d0
      ei = 0.0d0
      ek = 0.0d0
      utu = 0.0d0
      do i = 1, 3
         dm(i) = 0.0d0
         dp(i) = 0.0d0
         do j = 1, 3
            dmi(j,i) = 0.0d0
            dmk(j,i) = 0.0d0
            dpi(j,i) = 0.0d0
            dpk(j,i) = 0.0d0
         end do
      end do
c
c     check for presence of quadrupole components at either site
c
      iquad = (polsiz(ii) .ge. 13)
      kquad = (polsiz(kk) .ge. 13)
c
c     set permanent and induced multipoles for first site
c
      i0 = rpi(1)
      i1 = rpi(2)
      i2 = rpi(3)
      i3 = rpi(4)
      i4 = rpi(5)
      i5 = rpi(6)
      i6 = rpi(7)
      i7 = rpi(9)
      i8 = rpi(10)
      i9 = rpi(13)
      u1 = indi(1)
      u2 = indi(2)
      u3 = indi(3)
c
c     set permanent and induced multipoles for second site
c
      k0 = rpk(1)
      k1 = rpk(2)
      k2 = rpk(3)
      k3 = rpk(4)
      k4 = rpk(5)
      k5 = rpk(6)
      k6 = rpk(7)
      k7 = rpk(9)
      k8 = rpk(10)
      k9 = rpk(13)
      v1 = indk(1)
      v2 = indk(2)
      v3 = indk(3)
c
c     compute the zeroth order T2 matrix element
c
      rr1 = 1.0d0 / r
      t1 = rr1
c
c     compute the first order T2 matrix elements
c
      rr2 = rr1 * rr1
      rr3 = rr1 * rr2
      t2 = -xr * rr3
      t3 = -yr * rr3
      t4 = -zr * rr3
c
c     compute the second order T2 matrix elements
c
      rr5 = 3.0d0 * rr3 * rr2
      xr5 = xr * rr5
      yr5 = yr * rr5
      t5 = -rr3 + xr5*xr
      t6 = xr5 * yr
      t7 = xr5 * zr
      t8 = -rr3 + yr5*yr
      t9 = yr5 * zr
      t10 = -t5 - t8
c
c     compute the third order T2 matrix elements
c
      rr7 = 5.0d0 * rr5 * rr2
      xx = xr * xr
      yy = yr * yr
      xyz = xr * yr * zr
      xr7 = xx * rr7
      yr7 = yy * rr7
      rr53 = 3.0d0 * rr5
      t11 = xr * (rr53-xr7)
      t12 = yr * (rr5-xr7)
      t13 = zr * (rr5-xr7)
      t14 = xr * (rr5-yr7)
      t15 = -xyz * rr7
      t16 = -t11 - t14
      t17 = yr * (rr53-yr7)
      t18 = zr * (rr5-yr7)
      t19 = -t12 - t17
c
c     compute the fourth order T2 matrix elements
c
      if (iquad .or. kquad) then
         rr9 = 7.0d0 * rr7 * rr2
         if (xr .eq. 0.0d0) then
            yrx = 0.0d0
            zrx = 0.0d0
         else
            yrx = yr / xr
            zrx = zr / xr
         end if
         if (yr .eq. 0.0d0) then
            xry = 0.0d0
            zry = 0.0d0
         else
            xry = xr / yr
            zry = zr / yr
         end if
         xr9 = xx * xx * rr9
         xyr9 = xx * yy * rr9
         yr9 = yy * yy * rr9
         xr73 = 3.0d0 * xr7
         yr73 = 3.0d0 * yr7
         xr9x7 = xr9 - xr73
         yr9y7 = yr9 - yr73
         t21 = xr9x7 - xr73 + rr53
         t22 = yrx * xr9x7
         t23 = zrx * xr9x7
         t24 = xyr9 - xr7 - yr7 + rr5
         t25 = zry * (xyr9-yr7)
         t26 = -t21 - t24
         t27 = xry * yr9y7
         t28 = zrx * (xyr9-xr7)
         t29 = -t22 - t27
         t31 = yr9y7 - yr73 + rr53
         t32 = zry * yr9y7
         t33 = -t24 - t31
c
c     compute the fifth order T2 matrix elements
c
         r4 = r2 * r2
         rr5 = rr2 * rr3
         rr7 = rr2 * rr5
         rr9 = rr2 * rr7
         rr11 = rr2 * rr9
         x2 = xr * xr
         y2 = yr * yr
         z2 = zr * zr
         x4 = x2 * x2
         y4 = y2 * y2
         x2r2 = x2 * r2
         y2r2 = y2 * r2
         x2y2 = x2 * y2
         xrr11 = xr * rr11
         yrr11 = yr * rr11
         zrr11 = zr * rr11
         xyzrr11 = 315.0d0 * xyz * rr11
         r445 = 45.0d0 * r4
         n945x4 = -945.0d0 * x4
         n945y4 = -945.0d0 * y4
         n945x2y2 = -945.0d0 * x2y2
         nnn1 = n945x4 + 630.0d0*x2r2 - r445
         nnn2 = n945y4 + 630.0d0*y2r2 - r445
         r4225 = 225.0d0 * r4
         t36 = (n945x4+1050.0d0*x2r2-r4225) * xrr11
         t37 = nnn1 * yrr11
         t38 = nnn1 * zrr11
         t39 = (n945x2y2+105.0d0*x2r2+315.0d0*y2r2-r445) * xrr11
         t40 = (r2-3.0d0*x2) * xyzrr11
         t41 = -t36 - t39
         t42 = (n945x2y2+105.0d0*y2r2+315.0d0*x2r2-r445) * yrr11
         t43 = (n945x2y2+105.0d0*(x2r2+y2r2)-15.0d0*r4) * zrr11
         t44 = -t37 - t42
         t46 = nnn2 * xrr11
         t47 = (r2-3.0d0*y2) * xyzrr11
         t48 = -t39 - t46
         t51 = (n945y4+1050.0d0*y2r2-r4225) * yrr11
         t52 = nnn2 * zrr11
         t53 = -t42 - t51
      end if
c
c     get the M-M, M-D and D-D parts of the multipole energy
c
      i3k3 = i3 * k3
      w1 = i0 * k0
      w2 = i0*k1 - i1*k0
      w3 = i0*k2 - i2*k0
      w4 = i0*k3 - i3*k0
      w5 = i3k3 - i1*k1
      w6 = -i1*k2 - i2*k1
      w7 = -i1*k3 - i3*k1
      w8 = i3k3 - i2*k2
      w9 = -i2*k3 - i3*k2
      eik = w1*t1 + w2*t2 + w3*t3 + w4*t4 + w5*t5
     &         + w6*t6 + w7*t7 + w8*t8 + w9*t9
      dm(1) = w1*t2 + w2*t5 + w3*t6 + w4*t7 + w5*t11
     &           + w6*t12 + w7*t13 + w8*t14 + w9*t15
      dm(2) = w1*t3 + w2*t6 + w3*t8 + w4*t9 + w5*t12
     &           + w6*t14 + w7*t15 + w8*t17 + w9*t18
      dm(3) = w1*t4 + w2*t7 + w3*t9 + w4*t10 + w5*t13
     &           + w6*t15 + w7*t16 + w8*t18 + w9*t19
c
c     get the M-Q and D-Q parts of the multipole energy
c
      if (kquad) then
         k4k9 = k4 - k9
         k7k9 = k7 - k9
         i3k6 = 2.0d0 * i3 * k6
         i3k8 = 2.0d0 * i3 * k8
         w1 = i0 * k4k9
         w2 = 2.0d0 * i0 * k5
         w3 = 2.0d0 * i0 * k6
         w4 = i0 * k7k9
         w5 = 2.0d0 * i0 * k8
         w6 = i3k6 - i1*k4k9
         w7 = i3k8 - i2*k4k9 - 2.0d0*i1*k5
         w8 = -2.0d0*i1*k6 - i3*k4k9
         w9 = i3k6 - i1*k7k9 - 2.0d0*i2*k5
         w10 = -2.0d0 * (i1*k8 + i2*k6 + i3*k5)
         w11 = i3k8 - i2*k7k9
         w12 = -2.0d0*i2*k8 - i3*k7k9
         eik = eik + w1*t5 + w2*t6 + w3*t7 + w4*t8
     &            + w5*t9 + w6*t11 + w7*t12 + w8*t13
     &            + w9*t14 + w10*t15 + w11*t17 + w12*t18
         dm(1) = dm(1) + w1*t11 + w2*t12 + w3*t13 + w4*t14
     &              + w5*t15 + w6*t21 + w7*t22 + w8*t23
     &              + w9*t24 + w10*t25 + w11*t27 + w12*t28
         dm(2) = dm(2) + w1*t12 + w2*t14 + w3*t15 + w4*t17
     &              + w5*t18 + w6*t22 + w7*t24 + w8*t25
     &              + w9*t27 + w10*t28 + w11*t31 + w12*t32
         dm(3) = dm(3) + w1*t13 + w2*t15 + w3*t16 + w4*t18
     &              + w5*t19 + w6*t23 + w7*t25 + w8*t26
     &              + w9*t28 + w10*t29 + w11*t32 + w12*t33
      end if
c
c     get the M-Q and D-Q parts of the multipole energy
c
      if (iquad) then
         i4i9 = i4 - i9
         i7i9 = i7 - i9
         i6k3 = -2.0d0 * i6 * k3
         i8k3 = -2.0d0 * i8 * k3
         w1 = i4i9 * k0
         w2 = 2.0d0 * i5 * k0
         w3 = 2.0d0 * i6 * k0
         w4 = i7i9 * k0
         w5 = 2.0d0 * i8 * k0
         w6 = i4i9*k1 + i6k3
         w7 = i8k3 + i4i9*k2 + 2.0d0*i5*k1
         w8 = 2.0d0*i6*k1 + i4i9*k3
         w9 = i6k3 + i7i9*k1 + 2.0d0*i5*k2
         w10 = 2.0d0 * (i8*k1 + i6*k2 + i5*k3)
         w11 = i8k3 + i7i9*k2
         w12 = 2.0d0*i8*k2 + i7i9*k3
         eik = eik + w1*t5 + w2*t6 + w3*t7 + w4*t8
     &            + w5*t9 + w6*t11 + w7*t12 + w8*t13
     &            + w9*t14 + w10*t15 + w11*t17 + w12*t18
         dm(1) = dm(1) + w1*t11 + w2*t12 + w3*t13 + w4*t14
     &              + w5*t15 + w6*t21 + w7*t22 + w8*t23
     &              + w9*t24 + w10*t25 + w11*t27 + w12*t28
         dm(2) = dm(2) + w1*t12 + w2*t14 + w3*t15 + w4*t17
     &              + w5*t18 + w6*t22 + w7*t24 + w8*t25
     &              + w9*t27 + w10*t28 + w11*t31 + w12*t32
         dm(3) = dm(3) + w1*t13 + w2*t15 + w3*t16 + w4*t18
     &              + w5*t19 + w6*t23 + w7*t25 + w8*t26
     &              + w9*t28 + w10*t29 + w11*t32 + w12*t33
      end if
c
c     get the Q-Q part of the multipole interaction energy
c
      if (iquad .and. kquad) then
         ik49 = i4*k9 + i9*k4
         ik59 = i5*k9 + i9*k5
         ik69 = i6*k9 + i9*k6
         ik79 = i7*k9 + i9*k7
         ik89 = i8*k9 + i9*k8
         ik68 = 2.0d0 * (i6*k8 + i8*k6)
         i6k6 = i6 * k6
         i8k8 = i8 * k8
         i9k9 = i9 * k9
         w1 = i4*k4 - ik49 + i9k9 - 4.0d0*i6k6
         w2 = 2.0d0 * (i5*k4 + i4*k5 - ik68 - ik59)
         w3 = 2.0d0 * (i4*k6 + i6*k4 - ik69)
         w4 = i4*k7 + i7*k4 - ik49 - ik79 + 2.0d0*i9k9
     &           +4.0d0*(i5*k5-i6k6-i8k8)
         w5 = 2.0d0 * (i4*k8 + i8*k4 + 2.0d0*(i5*k6+i6*k5) - ik89)
         w6 = 2.0d0 * (i5*k7 + i7*k5 - ik68 - ik59)
         w7 = 2.0d0 * (i6*k7 + i7*k6 + 2.0d0*(i5*k8+i8*k5) - ik69)
         w8 = i7*k7 - ik79 + i9k9 - 4.0d0*i8k8
         w9 = 2.0d0 * (i7*k8 + i8*k7 - ik89)
         eik = eik + w1*t21 + w2*t22 + w3*t23 + w4*t24
     &            + w5*t25 + w6*t27 + w7*t28 + w8*t31 + w9*t32
         dm(1) = dm(1) + w1*t36 + w2*t37 + w3*t38 + w4*t39
     &              + w5*t40 + w6*t42 + w7*t43 + w8*t46 + w9*t47
         dm(2) = dm(2) + w1*t37 + w2*t39 + w3*t40 + w4*t42
     &              + w5*t43 + w6*t46 + w7*t47 + w8*t51 + w9*t52
         dm(3) = dm(3) + w1*t38 + w2*t40 + w3*t41 + w4*t43
     &              + w5*t44 + w6*t47 + w7*t48 + w8*t52 + w9*t53
      end if
c
c     get the (dM2/dx)*T*M1 terms for dipoles at both sites
c
      do i = 1, 3
         do j = 1, 3
            di1 = dpole(2,j,i,ii)
            di2 = dpole(3,j,i,ii)
            di3 = dpole(4,j,i,ii)
            dk1 = dpole(2,j,i,kk)
            dk2 = dpole(3,j,i,kk)
            dk3 = dpole(4,j,i,kk)
            w1 = k3 * di3
            dmi(j,i) = -k0*di1*t2 - k0*di2*t3 - k0*di3*t4
     &                    + (w1-k1*di1)*t5 - (k2*di1+k1*di2)*t6
     &                    - (k1*di3+k3*di1)*t7 + (w1-k2*di2)*t8
     &                    - (k2*di3+k3*di2)*t9
            w1 = dk3 * i3
            dmk(j,i) = dk1*i0*t2 + dk2*i0*t3 + dk3*i0*t4
     &                    + (w1-dk1*i1)*t5 - (dk2*i1+dk1*i2)*t6
     &                    - (dk1*i3+dk3*i1)*t7 + (w1-dk2*i2)*t8
     &                    - (dk2*i3+dk3*i2)*t9
            w1 = v3 * di3
            dpi(j,i) = (w1-v1*di1)*t5 - (v2*di1+v1*di2)*t6
     &                    - (v1*di3+v3*di1)*t7 + (w1-v2*di2)*t8
     &                    - (v3*di2+v2*di3)*t9
            w1 = u3 * dk3
            dpk(j,i) = (w1-u1*dk1)*t5 - (u2*dk1+u1*dk2)*t6
     &                    - (u1*dk3+u3*dk1)*t7 + (w1-u2*dk2)*t8
     &                    - (u3*dk2+u2*dk3)*t9
c
c     get the (dM2/dx)*T*M1 terms for quadrupole at first site
c
            if (iquad) then
               di4 = dpole(5,j,i,ii)
               di5 = dpole(6,j,i,ii)
               di6 = dpole(7,j,i,ii)
               di7 = dpole(9,j,i,ii)
               di8 = dpole(10,j,i,ii)
               di9 = dpole(13,j,i,ii)
               w1 = i9 - i4
               w2 = i9 - i7
               w3 = i6 * dk3
               w4 = i8 * dk3
               dmk(j,i) = dmk(j,i) - (w1*dk1+2.0d0*w3)*t11
     &                       - (w1*dk2+2.0d0*(w4-dk1*i5))*t12
     &                       - (w1*dk3-2.0d0*dk1*i6)*t13
     &                       - (w2*dk1+2.0d0*(w3-dk2*i5))*t14
     &                       + 2.0d0*(dk3*i5+dk2*i6+dk1*i8)*t15
     &                       - (w2*dk2+2.0d0*w4)*t17
     &                       - (w2*dk3-2.0d0*dk2*i8)*t18
               w1 = di9 - di4
               w2 = di9 - di7
               w3 = di6 * k3
               w4 = di8 * k3
               dmi(j,i) = dmi(j,i) - k0*(w1*t5+w2*t8)
     &                       + 2.0d0*k0*(di5*t6+di6*t7+di8*t9)
     &                       - (w1*k1+2.0d0*w3)*t11
     &                       - (w1*k2+2.0d0*(w4-k1*di5))*t12
     &                       - (w1*k3-2.0d0*k1*di6)*t13
     &                       - (w2*k1+2.0d0*(w3-k2*di5))*t14
     &                       + 2.0d0*(k3*di5+k2*di6+k1*di8)*t15
     &                       - (w2*k2+2.0d0*w4)*t17
     &                       - (w2*k3-2.0d0*k2*di8)*t18
c              w1 = di9 - di4
c              w2 = di9 - di7
               w3 = di6 * v3
               w4 = di8 * v3
               dpi(j,i) = dpi(j,i) - (w1*v1+2.0d0*w3)*t11
     &                       - (w1*v2+2.0d0*(w4-v1*di5))*t12
     &                       - (w1*v3-2.0d0*v1*di6)*t13
     &                       - (w2*v1+2.0d0*(w3-v2*di5))*t14
     &                       + 2.0d0*(v3*di5+v2*di6+v1*di8)*t15
     &                       - (w2*v2+2.0d0*w4)*t17
     &                       - (w2*v3-2.0d0*v2*di8)*t18
            end if
c
c     get the (dM2/dx)*T*M1 terms for quadrupole at second site
c
            if (kquad) then
               dk4 = dpole(5,j,i,kk)
               dk5 = dpole(6,j,i,kk)
               dk6 = dpole(7,j,i,kk)
               dk7 = dpole(9,j,i,kk)
               dk8 = dpole(10,j,i,kk)
               dk9 = dpole(13,j,i,kk)
               w1 = k9 - k4
               w2 = k9 - k7
               w3 = k6 * di3
               w4 = k8 * di3
               dmi(j,i) = dmi(j,i)
     &                       + (w1*di1+2.0d0*w3)*t11
     &                       + (w1*di2+2.0d0*(w4-k5*di1))*t12
     &                       + (w1*di3-2.0d0*k6*di1)*t13
     &                       + (w2*di1+2.0d0*(w3-k5*di2))*t14
     &                       - 2.0d0*(k5*di3+k6*di2+k8*di1)*t15
     &                       + (w2*di2+2.0d0*w4)*t17
     &                       + (w2*di3-2.0d0*k8*di2)*t18
               w1 = dk9 - dk4
               w2 = dk9 - dk7
               w3 = dk6 * i3
               w4 = dk8 * i3
               dmk(j,i) = dmk(j,i) - i0*(w1*t5+w2*t8)
     &                       + 2.0d0*i0*(dk5*t6+dk6*t7+dk8*t9)
     &                       + (w1*i1+2.0d0*w3)*t11
     &                       + (w1*i2+2.0d0*(w4-dk5*i1))*t12
     &                       + (w1*i3-2.0d0*dk6*i1)*t13
     &                       + (w2*i1+2.0d0*(w3-dk5*i2))*t14
     &                       - 2.0d0*(dk8*i1+dk5*i3+dk6*i2)*t15
     &                       + (w2*i2+2.0d0*w4)*t17
     &                       + (w2*i3-2.0d0*dk8*i2)*t18
c              w1 = dk9 - dk4
c              w2 = dk9 - dk7
               w3 = dk6 * u3
               w4 = dk8 * u3
               dpk(j,i) = dpk(j,i) + (w1*u1+2.0d0*w3)*t11
     &                       + (w1*u2+2.0d0*(w4-u1*dk5))*t12
     &                       + (w1*u3-2.0d0*u1*dk6)*t13
     &                       + (w2*u1+2.0d0*(w3-u2*dk5))*t14
     &                       - 2.0d0*(u1*dk8+u2*dk6+u3*dk5)*t15
     &                       + (w2*u2+2.0d0*w4)*t17
     &                       + (w2*u3-2.0d0*u2*dk8)*t18
            end if
c
c     get the (dM2/dx)*T*M1 terms for quadrupoles at both sites
c
            if (iquad .and. kquad) then
               w1 = di9 - di4
               w2 = di9 - di7
               w3 = di5*k9 + 2.0d0*(di6*k8+di8*k6)
               w4 = di6 * k6
               w5 = di6 * k9
               w6 = di8 * k8
               w7 = di8 * k9
               dmi(j,i) = dmi(j,i) + (w1*(k9-k4)-4.0d0*w4)*t21
     &                       + 2.0d0*(di5*k4-w1*k5-w3)*t22
     &                       + 2.0d0*(di6*k4-w1*k6-w5)*t23
     &                       + (4.0d0*(di5*k5-w4-w6)-w1*k7-w2*k4
     &                            +(2.0d0*di9-di4-di7)*k9)*t24
     &                       + 2.0d0*(di8*k4-w1*k8-w7
     &                            +2.0d0*(di6*k5+di5*k6))*t25
     &                       + 2.0d0*(di5*k7-w2*k5-w3)*t27
     &                       + 2.0d0*(di6*k7-w2*k6-w5
     &                            +2.0d0*(di8*k5+di5*k8))*t28
     &                       + (w2*(k9-k7)-4.0d0*w6)*t31
     &                       + 2.0d0*(di8*k7-w2*k8-w7)*t32
               w1 = dk9 - dk4
               w2 = dk9 - dk7
               w3 = dk5*i9 + 2.0d0*(dk6*i8+dk8*i6)
               w4 = dk6 * i6
               w5 = dk6 * i9
               w6 = dk8 * i8
               w7 = dk8 * i9
               dmk(j,i) = dmk(j,i) + (w1*(i9-i4)-4.0d0*w4)*t21
     &                       + 2.0d0*(dk5*i4-w1*i5-w3)*t22
     &                       + 2.0d0*(dk6*i4-w1*i6-w5)*t23
     &                       + (4.0d0*(dk5*i5-w4-w6)-w1*i7-w2*i4
     &                            +(2.0d0*dk9-dk4-dk7)*i9)*t24
     &                       + 2.0d0*(dk8*i4-w1*i8-w7
     &                            +2.0d0*(dk6*i5+dk5*i6))*t25
     &                       + 2.0d0*(dk5*i7-w2*i5-w3)*t27
     &                       + 2.0d0*(dk6*i7-w2*i6-w5
     &                            +2.0d0*(dk8*i5+dk5*i8))*t28
     &                       + (w2*(i9-i7)-4.0d0*w6)*t31
     &                       + 2.0d0*(dk8*i7-w2*i8-w7)*t32
            end if
         end do
      end do
c
c     get indD-M, indD-D and indD-Q polarization energy and gradients
c
      k3u3 = k3 * u3
      w1 = k0 * u1
      w2 = k0 * u2
      w3 = k0 * u3
      w4 = k3u3 - k1*u1
      w5 = k1*u2 + k2*u1
      w6 = k3*u1 + k1*u3
      w7 = k3u3 - k2*u2
      w8 = k2*u3 + k3*u2
      ei = -w1*t2 - w2*t3 - w3*t4 + w4*t5 - w5*t6
     &        - w6*t7 + w7*t8 - w8*t9
      dp(1) = -w1*t5 - w2*t6 - w3*t7 + w4*t11 - w5*t12
     &           - w6*t13 + w7*t14 - w8*t15
      dp(2) = -w1*t6 - w2*t8 - w3*t9 + w4*t12 - w5*t14
     &           - w6*t15 + w7*t17 - w8*t18
      dp(3) = -w1*t7 - w2*t9 - w3*t10 + w4*t13 - w5*t15
     &           - w6*t16 + w7*t18 - w8*t19
      if (kquad) then
         k6u3 = k6 * u3
         k8u3 = k8 * u3
         w1 = -k4k9*u1 + 2.0d0*k6u3
         w2 = -k4k9*u2 + 2.0d0*(k8u3-k5*u1)
         w3 = -k4k9*u3 - 2.0d0*k6*u1
         w4 = -k7k9*u1 + 2.0d0*(k6u3-k5*u2)
         w5 = 2.0d0 * (k5*u3 + k6*u2 + k8*u1)
         w6 = -k7k9*u2 + 2.0d0*k8u3
         w7 = -k7k9*u3 - 2.0d0*k8*u2
         ei = ei + w1*t11 + w2*t12 + w3*t13 + w4*t14
     &           - w5*t15 + w6*t17 + w7*t18
         dp(1) = dp(1) + w1*t21 + w2*t22 + w3*t23
     &              + w4*t24 - w5*t25 + w6*t27 + w7*t28
         dp(2) = dp(2) + w1*t22 + w2*t24 + w3*t25
     &              + w4*t27 - w5*t28 + w6*t31 + w7*t32
         dp(3) = dp(3) + w1*t23 + w2*t25 + w3*t26
     &              + w4*t28 - w5*t29 + w6*t32 + w7*t33
      end if
c
c     get M-indD, D-indD and Q-indD polarization energy and gradients
c
      i3v3 = i3 * v3
      w1 = i0 * v1
      w2 = i0 * v2
      w3 = i0 * v3
      w4 = i3v3 - i1*v1
      w5 = i1*v2 + i2*v1
      w6 = i3*v1 + i1*v3
      w7 = i3v3 - i2*v2
      w8 = i2*v3 + i3*v2
      ek = w1*t2 + w2*t3 + w3*t4 + w4*t5 - w5*t6
     &        - w6*t7 + w7*t8 - w8*t9
      dp(1) = dp(1) + w1*t5 + w2*t6 + w3*t7 + w4*t11 - w5*t12
     &           - w6*t13 + w7*t14 - w8*t15
      dp(2) = dp(2) + w1*t6 + w2*t8 + w3*t9 + w4*t12 - w5*t14
     &           - w6*t15 + w7*t17 - w8*t18
      dp(3) = dp(3) + w1*t7 + w2*t9 + w3*t10 + w4*t13 - w5*t15
     &           - w6*t16 + w7*t18 - w8*t19
      if (iquad) then
         i6v3 = i6 * v3
         i8v3 = i8 * v3
         w1 = i4i9*v1 - 2.0d0*i6v3
         w2 = i4i9*v2 + 2.0d0*(i5*v1-i8v3)
         w3 = i4i9*v3 + 2.0d0*i6*v1
         w4 = i7i9*v1 + 2.0d0*(i5*v2-i6v3)
         w5 = 2.0d0 * (i5*v3 + i6*v2 + i8*v1)
         w6 = i7i9*v2 - 2.0d0*i8v3
         w7 = i7i9*v3 + 2.0d0*i8*v2
         ek = ek + w1*t11 + w2*t12 + w3*t13 + w4*t14
     &           + w5*t15 + w6*t17 + w7*t18
         dp(1) = dp(1) + w1*t21 + w2*t22 + w3*t23
     &              + w4*t24 + w5*t25 + w6*t27 + w7*t28
         dp(2) = dp(2) + w1*t22 + w2*t24 + w3*t25
     &              + w4*t27 + w5*t28 + w6*t31 + w7*t32
         dp(3) = dp(3) + w1*t23 + w2*t25 + w3*t26
     &              + w4*t28 + w5*t29 + w6*t32 + w7*t33
      end if
c
c     compute the mutual polarization induced dipole gradient terms
c
      if (poltyp .eq. 'MUTUAL') then
         w1 = -u1*v1 + u3*v3
         w2 = -u2*v1 - u1*v2
         w3 = -u3*v1 - u1*v3
         w4 = -u2*v2 + u3*v3
         w5 = -u3*v2 - u2*v3
         utu = w1*t5 + w2*t6 + w3*t7 + w4*t8 + w5*t9
         dp(1) = dp(1) + w1*t11 + w2*t12 + w3*t13 + w4*t14 + w5*t15
         dp(2) = dp(2) + w1*t12 + w2*t14 + w3*t15 + w4*t17 + w5*t18
         dp(3) = dp(3) + w1*t13 + w2*t15 + w3*t16 + w4*t18 + w5*t19
      end if
c
c     apply a damping factor to polarization energy and derivatives
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         term = -pgamma * (r/damp)**3
         if (term .gt. -50.0d0) then
            term = exp(term)
            ddamp = (3.0d0*pgamma*r2/damp**3) * term
            damp = 1.0d0 - term
            de = (ei+ek+utu) * ddamp
            ei = ei * damp
            ek = ek * damp
            dp(1) = dp(1)*damp + de*(xr/r)
            dp(2) = dp(2)*damp + de*(yr/r)
            dp(3) = dp(3)*damp + de*(zr/r)
            do i = 1, 3
               do j = 1, 3
                  dpi(j,i) = dpi(j,i) * damp
                  dpk(j,i) = dpk(j,i) * damp
               end do
            end do
         end if
      end if
c
c     final polarization energy is half of value computed above
c
      ei = 0.5d0 * ei
      ek = 0.5d0 * ek
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole2  --  multipole & polarization Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole2" calculates second derivatives of the multipole
c     and dipole polarization energy for a single atom at a time
c
c
      subroutine empole2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'hessn.i'
      integer i,j,k
      real*8 eps,old
      real*8 d0(3,maxatm)
c
c
c     set the stepsize to be used for numerical derivatives
c
      eps = 1.0d-7
c
c     get multipole first derivatives for the base structure
c
      call empole2a (i)
      do k = 1, n
         do j = 1, 3
            d0(j,k) = dem(j,k) + dep(j,k)
         end do
      end do
c
c     find numerical x-components via perturbed structures
c
      old = x(i)
      x(i) = x(i) + eps
      call empole2a (i)
      x(i) = old
      do k = 1, n
         do j = 1, 3
            hessx(j,k) = hessx(j,k) + (dem(j,k)+dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical y-components via perturbed structures
c
      old = y(i)
      y(i) = y(i) + eps
      call empole2a (i)
      y(i) = old
      do k = 1, n
         do j = 1, 3
            hessy(j,k) = hessy(j,k) + (dem(j,k)+dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical z-components via perturbed structures
c
      old = z(i)
      z(i) = z(i) + eps
      call empole2a (i)
      z(i) = old
      do k = 1, n
         do j = 1, 3
            hessz(j,k) = hessz(j,k) + (dem(j,k)+dep(j,k)-d0(j,k))/eps
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empole2a  --  mpole & polar Hessian; numerical  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empole2a" computes multipole and dipole polarization first
c     derivatives for a single atom with respect to Cartesian
c     coordinates; used to get finite difference second derivatives
c
c     note that since polarization effects are many body, it is not
c     really correct to neglect interactions where "iatom" is not
c     directly involved as a multipole site or local coordinate axis;
c     however, other sites are neglected in this version via the
c     "ipart" and "kpart" checks to quickly get approximate values
c
c     also, in the current version, the induced dipoles are not
c     recomputed every time an atom is moved during computation of
c     the numerical Hessian resulting in a further approximation
c
c
      subroutine empole2a (iatom)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'deriv.i'
      include 'group.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,m
      integer jj,iatom
      integer ii,iz,ix
      integer kk,kz,kx
      real*8 eik,ei,ek,de
      real*8 f,fikm,fikp,fgrp
      real*8 shift,taper,dtaper
      real*8 trans,dtrans
      real*8 xi,yi,zi,xr,yr,zr
      real*8 r,r2,r3,r4,r5,r6,r7
      real*8 a(3,3),d(3,3,3,3)
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      real*8 dm(3),dp(3),utu
      real*8 dmi(3,3),dmk(3,3)
      real*8 dpi(3,3),dpk(3,3)
      real*8 mscale(maxatm)
      real*8 pscale(maxatm)
      logical proceed,iuse,kuse
      logical ipart,kpart
c
c
c     zero out the multipole and polarization first derivatives
c
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('MPOLE')
c
c     rotate multipole components and get rotation derivatives
c
      do i = 1, npole
         call rotmat (i,a)
         call rotsite (i,a)
         call drotmat (i,d)
         do j = 1, 3
            do k = 1, 3
               call drotpole (i,a,d,j,k)
            end do
         end do
      end do
c
c     compute the induced dipoles at each polarizable atom;
c     currently, we assume this was done in a calling routine
c
c     call induce
c
c     compute and partition multipole interaction energy
c
      do ii = 1, npole-1
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         iz = zaxis(ii)
         ix = xaxis(ii)
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         ipart = (i.eq.iatom .or. iz.eq.iatom .or. ix.eq.iatom)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = ii+1, npole
            mscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(i)
            mscale(i12(j,i)) = m2scale
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = m3scale
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = m4scale
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = m5scale
         end do
         do j = 1, np11(i)
            pscale(ip11(j,i)) = p1scale
         end do
         do j = 1, np12(i)
            pscale(ip12(j,i)) = p2scale
         end do
         do j = 1, np13(i)
            pscale(ip13(j,i)) = p3scale
         end do
         do j = 1, np14(i)
            pscale(ip14(j,i)) = p4scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kpart = (k.eq.iatom .or. kz.eq.iatom .or. kx.eq.iatom)
            kuse = (use(k) .or. use(kz) .or. use(kx))
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            proceed = .true.
            proceed = (ipart .or. kpart)
            if (proceed)  proceed = (iuse .or. kuse)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  do j = 1, maxpole
                     rpk(j) = rpole(j,kk)
                  end do
                  do j = 1, 3
                     indk(j) = uind(j,kk)
                  end do
                  r = sqrt(r2)
                  call empik1 (ii,kk,xr,yr,zr,r,r2,rpi,
     &                         rpk,indi,indk,eik,ei,ek,
     &                         dm,dmi,dmk,dp,dpi,dpk,utu)
                  fikm = f * mscale(k)
                  fikp = f * pscale(k)
                  eik = fikm * eik
                  ei = fikp * ei
                  ek = fikp * ek
                  utu = fikp * utu
                  do j = 1, 3
                     dm(j) = fikm * dm(j)
                     dp(j) = fikp * dp(j)
                     do jj = 1, 3
                        dmi(jj,j) = fikm * dmi(jj,j)
                        dmk(jj,j) = fikm * dmk(jj,j)
                        dpi(jj,j) = fikp * dpi(jj,j)
                        dpk(jj,j) = fikp * dpk(jj,j)
                     end do
                  end do
c
c     use shifted energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     fikm = fikm * rpi(1) * rpk(1)
                     shift = fikm / (0.5d0*(off+cut))
                     eik = eik - shift
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     trans = fikm * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                  + f3*r3 + f2*r2 + f1*r + f0)
                     dtrans = fikm * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                                  + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                              + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                     de = eik*dtaper + dtrans
                     dm(1) = dm(1)*taper + de*(xr/r)
                     dm(2) = dm(2)*taper + de*(yr/r)
                     dm(3) = dm(3)*taper + de*(zr/r)
                     de = (2.0d0*(ei+ek)+utu) * dtaper
                     dp(1) = dp(1)*taper + de*(xr/r)
                     dp(2) = dp(2)*taper + de*(yr/r)
                     dp(3) = dp(3)*taper + de*(zr/r)
                     do j = 1, 3
                        do jj = 1, 3
                           dmi(jj,j) = dmi(jj,j) * taper
                           dmk(jj,j) = dmk(jj,j) * taper
                           dpi(jj,j) = dpi(jj,j) * taper
                           dpk(jj,j) = dpk(jj,j) * taper
                        end do
                     end do
                  end if
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                  if (use_group) then
                     do j = 1, 3
                        dm(j) = dm(j) * fgrp
c                       dp(j) = dp(j) * fgrp
                        do jj = 1, 3
                           dmi(jj,j) = dmi(jj,j) * fgrp
                           dmk(jj,j) = dmk(jj,j) * fgrp
c                          dpi(jj,j) = dpi(jj,j) * fgrp
c                          dpk(jj,j) = dpk(jj,j) * fgrp
                        end do
                     end do
                  end if
c
c     increment the multipole first derivative expressions
c
                  dem(1,i) = dem(1,i) - dm(1) + dmi(1,1)
                  dem(2,i) = dem(2,i) - dm(2) + dmi(1,2)
                  dem(3,i) = dem(3,i) - dm(3) + dmi(1,3)
                  dem(1,iz) = dem(1,iz) + dmi(2,1)
                  dem(2,iz) = dem(2,iz) + dmi(2,2)
                  dem(3,iz) = dem(3,iz) + dmi(2,3)
                  dem(1,ix) = dem(1,ix) + dmi(3,1)
                  dem(2,ix) = dem(2,ix) + dmi(3,2)
                  dem(3,ix) = dem(3,ix) + dmi(3,3)
                  dem(1,k) = dem(1,k) + dm(1) + dmk(1,1)
                  dem(2,k) = dem(2,k) + dm(2) + dmk(1,2)
                  dem(3,k) = dem(3,k) + dm(3) + dmk(1,3)
                  dem(1,kz) = dem(1,kz) + dmk(2,1)
                  dem(2,kz) = dem(2,kz) + dmk(2,2)
                  dem(3,kz) = dem(3,kz) + dmk(2,3)
                  dem(1,kx) = dem(1,kx) + dmk(3,1)
                  dem(2,kx) = dem(2,kx) + dmk(3,2)
                  dem(3,kx) = dem(3,kx) + dmk(3,3)
c
c     increment the polarization first derivative expressions
c
                  dep(1,i) = dep(1,i) - dp(1) + dpi(1,1)
                  dep(2,i) = dep(2,i) - dp(2) + dpi(1,2)
                  dep(3,i) = dep(3,i) - dp(3) + dpi(1,3)
                  dep(1,iz) = dep(1,iz) + dpi(2,1)
                  dep(2,iz) = dep(2,iz) + dpi(2,2)
                  dep(3,iz) = dep(3,iz) + dpi(2,3)
                  dep(1,ix) = dep(1,ix) + dpi(3,1)
                  dep(2,ix) = dep(2,ix) + dpi(3,2)
                  dep(3,ix) = dep(3,ix) + dpi(3,3)
                  dep(1,k) = dep(1,k) + dp(1) + dpk(1,1)
                  dep(2,k) = dep(2,k) + dp(2) + dpk(1,2)
                  dep(3,k) = dep(3,k) + dp(3) + dpk(1,3)
                  dep(1,kz) = dep(1,kz) + dpk(2,1)
                  dep(2,kz) = dep(2,kz) + dpk(2,2)
                  dep(3,kz) = dep(3,kz) + dpk(2,3)
                  dep(1,kx) = dep(1,kx) + dpk(3,1)
                  dep(2,kx) = dep(2,kx) + dpk(3,2)
                  dep(3,kx) = dep(3,kx) + dpk(3,3)
               end if
            end if
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
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         iz = zaxis(ii)
         ix = xaxis(ii)
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         ipart = (i.eq.iatom .or. iz.eq.iatom .or. ix.eq.iatom)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = ii, npole
            mscale(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(i)
            mscale(i12(j,i)) = m2scale
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = m3scale
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = m4scale
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = m5scale
         end do
         do j = 1, np11(i)
            pscale(ip11(j,i)) = p1scale
         end do
         do j = 1, np12(i)
            pscale(ip12(j,i)) = p2scale
         end do
         do j = 1, np13(i)
            pscale(ip13(j,i)) = p3scale
         end do
         do j = 1, np14(i)
            pscale(ip14(j,i)) = p4scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kpart = (k.eq.iatom .or. kz.eq.iatom .or. kx.eq.iatom)
            kuse = (use(k) .or. use(kz) .or. use(kx))
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            proceed = .true.
            proceed = (ipart .or. kpart)
            if (proceed)  proceed = (iuse .or. kuse)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do m = 1, ncell
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  call image (xr,yr,zr,m)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     do j = 1, maxpole
                        rpk(j) = rpole(j,kk)
                     end do
                     do j = 1, 3
                        indk(j) = uind(j,kk)
                     end do
                     r = sqrt(r2)
                     call empik1 (ii,kk,xr,yr,zr,r,r2,rpi,
     &                            rpk,indi,indk,eik,ei,ek,
     &                            dm,dmi,dmk,dp,dpi,dpk,utu)
                     fikm = f
                     fikp = f
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           fikm = fikm * mscale(kk)
                           fikp = fikp * pscale(kk)
                        end if
                     end if
                     eik = fikm * eik
                     ei = fikp * ei
                     ek = fikp * ek
                     utu = fikp * utu
                     do j = 1, 3
                        dm(j) = fikm * dm(j)
                        dp(j) = fikp * dp(j)
                        do jj = 1, 3
                           dmi(jj,j) = fikm * dmi(jj,j)
                           dmk(jj,j) = fikm * dmk(jj,j)
                           dpi(jj,j) = fikp * dpi(jj,j)
                           dpk(jj,j) = fikp * dpk(jj,j)
                        end do
                     end do
c
c     use shifted energy switching if near the cutoff distance
c
                     if (r2 .gt. cut2) then
                        fikm = fikm * rpi(1) * rpk(1)
                        shift = fikm / (0.5d0*(off+cut))
                        eik = eik - shift
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        r6 = r3 * r3
                        r7 = r3 * r4
                        taper = c5*r5 + c4*r4 + c3*r3
     &                                + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                                 + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        trans = fikm * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                     + f3*r3 + f2*r2 + f1*r + f0)
                        dtrans = fikm * (7.0d0*f7*r6 + 6.0d0*f6*r5
     &                                     + 5.0d0*f5*r4 + 4.0d0*f4*r3
     &                                 + 3.0d0*f3*r2 + 2.0d0*f2*r + f1)
                        de = eik*dtaper + dtrans
                        dm(1) = dm(1)*taper + de*(xr/r)
                        dm(2) = dm(2)*taper + de*(yr/r)
                        dm(3) = dm(3)*taper + de*(zr/r)
                        de = (2.0d0*(ei+ek)+utu) * dtaper
                        dp(1) = dp(1)*taper + de*(xr/r)
                        dp(2) = dp(2)*taper + de*(yr/r)
                        dp(3) = dp(3)*taper + de*(zr/r)
                        do j = 1, 3
                           do jj = 1, 3
                              dmi(jj,j) = dmi(jj,j) * taper
                              dmk(jj,j) = dmk(jj,j) * taper
                              dpi(jj,j) = dpi(jj,j) * taper
                              dpk(jj,j) = dpk(jj,j) * taper
                           end do
                        end do
                     end if
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                     if (use_group) then
                        do j = 1, 3
                           dm(j) = dm(j) * fgrp
c                          dp(j) = dp(j) * fgrp
                           do jj = 1, 3
                              dmi(jj,j) = dmi(jj,j) * fgrp
                              dmk(jj,j) = dmk(jj,j) * fgrp
c                             dpi(jj,j) = dpi(jj,j) * fgrp
c                             dpk(jj,j) = dpk(jj,j) * fgrp
                           end do
                        end do
                     end if
c
c     increment the multipole first derivative expressions
c
                     dem(1,i) = dem(1,i) - dm(1) + dmi(1,1)
                     dem(2,i) = dem(2,i) - dm(2) + dmi(1,2)
                     dem(3,i) = dem(3,i) - dm(3) + dmi(1,3)
                     dem(1,iz) = dem(1,iz) + dmi(2,1)
                     dem(2,iz) = dem(2,iz) + dmi(2,2)
                     dem(3,iz) = dem(3,iz) + dmi(2,3)
                     dem(1,ix) = dem(1,ix) + dmi(3,1)
                     dem(2,ix) = dem(2,ix) + dmi(3,2)
                     dem(3,ix) = dem(3,ix) + dmi(3,3)
                     if (i .ne. k) then
                        dem(1,k) = dem(1,k) + dm(1) + dmk(1,1)
                        dem(2,k) = dem(2,k) + dm(2) + dmk(1,2)
                        dem(3,k) = dem(3,k) + dm(3) + dmk(1,3)
                        dem(1,kz) = dem(1,kz) + dmk(2,1)
                        dem(2,kz) = dem(2,kz) + dmk(2,2)
                        dem(3,kz) = dem(3,kz) + dmk(2,3)
                        dem(1,kx) = dem(1,kx) + dmk(3,1)
                        dem(2,kx) = dem(2,kx) + dmk(3,2)
                        dem(3,kx) = dem(3,kx) + dmk(3,3)
                     end if
c
c     increment the polarization first derivative expressions
c
                     dep(1,i) = dep(1,i) - dp(1) + dpi(1,1)
                     dep(2,i) = dep(2,i) - dp(2) + dpi(1,2)
                     dep(3,i) = dep(3,i) - dp(3) + dpi(1,3)
                     dep(1,iz) = dep(1,iz) + dpi(2,1)
                     dep(2,iz) = dep(2,iz) + dpi(2,2)
                     dep(3,iz) = dep(3,iz) + dpi(2,3)
                     dep(1,ix) = dep(1,ix) + dpi(3,1)
                     dep(2,ix) = dep(2,ix) + dpi(3,2)
                     dep(3,ix) = dep(3,ix) + dpi(3,3)
                     if (i .ne. k) then
                        dep(1,k) = dep(1,k) + dp(1) + dpk(1,1)
                        dep(2,k) = dep(2,k) + dp(2) + dpk(1,2)
                        dep(3,k) = dep(3,k) + dp(3) + dpk(1,3)
                        dep(1,kz) = dep(1,kz) + dpk(2,1)
                        dep(2,kz) = dep(2,kz) + dpk(2,2)
                        dep(3,kz) = dep(3,kz) + dpk(2,3)
                        dep(1,kx) = dep(1,kx) + dpk(3,1)
                        dep(2,kx) = dep(2,kx) + dpk(3,2)
                        dep(3,kx) = dep(3,kx) + dpk(3,3)
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3a  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole3a" calculates the electrostatic energy due to
c     atomic multipole interactions and dipole polarizability, and
c     partitions the energy among the atoms using a double loop
c
c
      subroutine empole3a
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moment.i'
      include 'mpole.i'
      include 'polar.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,m
      integer ii,iz,ix
      integer kk,kz,kx
      integer skip(maxatm)
      real*8 eik,ei,ek,a(3,3)
      real*8 shift,taper,trans
      real*8 f,fik,fgrp
      real*8 xr,yr,zr,r
      real*8 r2,r3,r4,r5,r6,r7
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      real*8 weight,xcenter,ycenter,zcenter
      real*8 totchg,xsum,ysum,zsum
      logical header,huge,proceed,iuse,kuse
c
c
c     zero out multipole and polarization energy and partitioning
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aep(i) = 0.0d0
      end do
      header = .true.
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set conversion factor and switching function coefficients
c
      f = electric / dielec
      call switch ('MPOLE')
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each atom
c
      call induce
c
c     find the center of mass of the set of active atoms
c
      weight = 0.0d0
      xcenter = 0.0d0
      ycenter = 0.0d0
      zcenter = 0.0d0
      do i = 1, n
         if (use(i)) then
            weight = weight + mass(i)
            xcenter = xcenter + x(i)*mass(i)
            ycenter = ycenter + y(i)*mass(i)
            zcenter = zcenter + z(i)*mass(i)
         end if
      end do
      if (weight .ne. 0.0d0) then
         xcenter = xcenter / weight
         ycenter = ycenter / weight
         zcenter = zcenter / weight
      end if
c
c     get net charge and dipole components over active atoms
c
      totchg = 0.0d0
      xsum = 0.0d0
      ysum = 0.0d0
      zsum = 0.0d0
      do i = 1, npole
         k = ipole(i)
         if (use(k)) then
            totchg = totchg + rpole(1,i)
            xsum = xsum + (x(k)-xcenter)*rpole(1,i)
            ysum = ysum + (y(k)-ycenter)*rpole(1,i)
            zsum = zsum + (z(k)-zcenter)*rpole(1,i)
            xsum = xsum + rpole(2,i) + uind(1,i)
            ysum = ysum + rpole(3,i) + uind(2,i)
            zsum = zsum + rpole(4,i) + uind(3,i)
         end if
      end do
      netchg = netchg + totchg
      xdipole = xdipole + debye*xsum
      ydipole = ydipole + debye*ysum
      zdipole = zdipole + debye*zsum
c
c     calculate the multipole interaction energy term
c
      do ii = 1, npole-1
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, n12(i)
            skip(i12(j,i)) = i * chg12use
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i * chg13use
         end do
         do j = 1, n14(i)
            skip(i14(j,i)) = i * chg14use
         end do
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. kuse)
            if (proceed)  proceed = (skip(k) .ne. i)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               if (use_image)  call image (xr,yr,zr,0)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  do j = 1, maxpole
                     rpk(j) = rpole(j,kk)
                  end do
                  do j = 1, 3
                     indk(j) = uind(j,kk)
                  end do
                  r = sqrt(r2)
                  call empik (ii,kk,xr,yr,zr,r,rpi,rpk,
     &                          indi,indk,eik,ei,ek)
                  fik = f
                  if (skip(k) .eq. -i)  fik = fik / chgscale
                  eik = fik * eik
                  ei = fik * ei
                  ek = fik * ek
c
c     use shifted energy switching if near the cutoff distance
c
                  fik = fik * rpi(1) * rpk(1)
                  shift = fik / (0.5d0*(off+cut))
                  eik = eik - shift
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     r6 = r3 * r3
                     r7 = r3 * r4
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                               + f3*r3 + f2*r2 + f1*r + f0)
                     eik = eik * taper + trans
                     ei = ei * taper
                     ek = ek * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     eik = eik * fgrp
                     ei = ei * fgrp
                     ek = ek * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  nem = nem + 1
                  em = em + eik
                  aem(i) = aem(i) + 0.5d0*eik
                  aem(k) = aem(k) + 0.5d0*eik
                  nep = nep + 2
                  ep = ep + ei + ek
                  aep(i) = aep(i) + 0.5d0*(ei+ek)
                  aep(k) = aep(k) + 0.5d0*(ei+ek)
c
c     print a warning if the energy of this interaction is large
c
                  huge = (max(abs(eik),abs(ei),abs(ek)) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Individual Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',13x,'Atom Names',
     &                             9x,'Distance',6x,'Energies',
     &                             ' (MPole, Pol1, Pol2)',/)
                     end if
                     write (iout,20)  i,name(i),k,name(k),r,eik,ei,ek
   20                format (' M-Pole',5x,i5,'-',a3,1x,i5,'-',a3,
     &                          5x,f8.4,3f12.4)
                  end if
               end if
            end if
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
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         do j = 1, maxpole
            rpi(j) = rpole(j,ii)
         end do
         do j = 1, 3
            indi(j) = uind(j,ii)
         end do
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            kuse = (use(k) .or. use(kz) .or. use(kx))
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,2,i,k,0,0,0)
            if (proceed)  proceed = (iuse .or. kuse)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do m = 1, ncell
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  call image (xr,yr,zr,m)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     do j = 1, maxpole
                        rpk(j) = rpole(j,kk)
                     end do
                     do j = 1, 3
                        indk(j) = uind(j,kk)
                     end do
                     r = sqrt(r2)
                     call empik (ii,kk,xr,yr,zr,r,rpi,rpk,
     &                             indi,indk,eik,ei,ek)
                     fik = f
                     eik = fik * eik
                     ei = fik * ei
                     ek = fik * ek
c
c     use shifted energy switching if near the cutoff distance
c
                     fik = fik * rpi(1) * rpk(1)
                     shift = fik / (0.5d0*(off+cut))
                     eik = eik - shift
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        r6 = r3 * r3
                        r7 = r3 * r4
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        trans = fik * (f7*r7 + f6*r6 + f5*r5 + f4*r4
     &                                  + f3*r3 + f2*r2 + f1*r + f0)
                        eik = eik * taper + trans
                        ei = ei * taper
                        ek = ek * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        eik = eik * fgrp
                        ei = ei * fgrp
                        ek = ek * fgrp
                     end if
c
c     increment the overall multipole and polarization energies
c
                     nem = nem + 1
                     nep = nep + 2
                     if (i .eq. k) then
                        em = em + 0.5d0*eik
                        aem(i) = aem(i) + 0.5d0*eik
                        ep = ep + ei
                        aep(i) = aep(i) + ei
                     else
                        em = em + eik
                        aem(i) = aem(i) + 0.5d0*eik
                        aem(k) = aem(k) + 0.5d0*eik
                        ep = ep + ei + ek
                        aep(i) = aep(i) + 0.5d0*(ei+ek)
                        aep(k) = aep(k) + 0.5d0*(ei+ek)
                     end if
c
c     print a warning if the energy of this interaction is large
c
                     huge = (max(abs(eik),abs(ei),abs(ek)) .gt. 100.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual Multipole and',
     &                                ' Polarization Interactions :',
     &                             //,' Type',13x,'Atom Names',
     &                                9x,'Distance',6x,'Energies',
     &                                ' (MPole, Pol1, Pol2)',/)
                        end if
                        write (iout,40)  i,name(i),k,name(k),r,eik,ei,ek
   40                   format (' M-Pole',5x,i5,'-',a3,1x,i5,'-',a3,
     &                             1x,'(X)',1x,f8.4,3f12.4)
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end

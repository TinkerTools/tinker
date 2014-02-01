c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine epolar  --  dipole polarization energy  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "epolar" calculates the induced dipole polarization energy
c
c
      subroutine epolar
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'energi.i'
      include 'mpole.i'
      include 'polar.i'
      include 'shunt.i'
      include 'usage.i'
      integer i,j,k,skip(maxatm)
      integer ii,kk,iz,ix
      real*8 e,epolik,a(3,3)
      real*8 xi,yi,zi,xr,yr,zr
      real*8 r,r2,r3,r4,r5,taper
      real*8 rpi(13),indk(3)
      logical iuse,kuse
c
c
c     zero out the induced dipole polarization energy
c
      ep = 0.0d0
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the switching function coefficients
c
      call switch ('CHGDPL')
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each atom
c
      call induce
c
c     calculate the induced dipole polarization energy term
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         xi = x(i)
         yi = y(i)
         zi = z(i)
         skip(i) = i
         do j = 1, n12(i)
            skip(i12(j,i)) = i
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i
         end do
         do j = 1, mdqsiz
            rpi(j) = rpole(j,ii)
         end do
         do kk = 1, npolar
            k = ipolar(kk)
            kuse = use(k)
            if (iuse .or. kuse) then
               if (skip(k) .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  if (use_image)  call image (xr,yr,zr,0)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     indk(1) = uind(1,k)
                     indk(2) = uind(2,k)
                     indk(3) = uind(3,k)
                     r = sqrt(r2)
                     e = epolik (r,xr,yr,zr,rpi,indk)
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        e = e * taper
                     end if
                     ep = ep + e
                  end if
               end if
            end if
         end do
      end do
      return
      end
c
c
c     #######################
c     ##                   ##
c     ##  function epolik  ##
c     ##                   ##
c     #######################
c
c
      function epolik (r,xr,yr,zr,rpi,indk)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'electr.i'
      include 'mpole.i'
      include 'units.i'
      integer i,j
      real*8 epolik,t2matrix
      real*8 r,r2,r3,r4
      real*8 xr,yr,zr,zrr
      real*8 xrr,xrr2,xrr3
      real*8 yrr,yrr2,yrr3
      real*8 rpi(13),indk(3)
      real*8 m2t2(13),t2(13,13)
      real*8 bx(0:5,0:5),by(11,0:5,0:1),bz(11,0:1)
c
c
c     zeroth order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 1) then
         bz(1,0) = c(1,0,0)
         by(1,0,0) = c(1,0,0) * bz(1,0)
         bx(0,0) = c(1,0,0)
c
c     first order coefficients to speed the T2 calculation
c
         r2 = r * r
         zrr = zr / r
         bz(3,0) = c(3,0,0)
         bz(1,1) = c(1,1,1) * zrr
         yrr = yr / r
         by(3,0,0) = c(3,0,0) * bz(3,0)
         by(1,0,1) = c(1,0,0) * bz(1,1)
         by(1,1,0) = c(1,1,1) * bz(3,0) * yrr
         xrr = xr / r
         bx(1,1) = c(1,1,1) * xrr
      end if
c
c     second order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 4) then
         r3 = r2 * r
         bz(5,0) = c(5,0,0)
         bz(3,1) = c(3,1,1) * zrr
         yrr2 = yrr * yrr
         by(5,0,0) = c(5,0,0) * bz(5,0)
         by(3,0,1) = c(3,0,0) * bz(3,1)
         by(3,1,0) = c(3,1,1) * bz(5,0) * yrr
         by(1,1,1) = c(1,1,1) * bz(3,1) * yrr
         by(1,2,0) = c(1,2,0)*bz(3,0) + c(1,2,2)*bz(5,0)*yrr2
         xrr2 = xrr * xrr
         bx(2,0) = c(1,2,0)
         bx(2,2) = c(1,2,2) * xrr2
      end if
c
c     third order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 13) then
         r4 = r2 * r2
         bz(7,0) = c(7,0,0)
         bz(5,1) = c(5,1,1) * zrr
         yrr3 = yrr2 * yrr
         by(7,0,0) = c(7,0,0)*bz(7,0)
         by(5,0,1) = c(5,0,0)*bz(5,1)
         by(5,1,0) = c(5,1,1)*bz(7,0)*yrr
         by(3,1,1) = c(3,1,1)*bz(5,1)*yrr
         by(3,2,0) = c(3,2,0)*bz(5,0) + c(3,2,2)*bz(7,0)*yrr2
         by(1,2,1) = c(1,2,0)*bz(3,1) + c(1,2,2)*bz(5,1)*yrr2
         by(1,3,0) = c(1,3,1)*bz(5,0)*yrr + c(1,3,3)*bz(7,0)*yrr3
         xrr3 = xrr2 * xrr
         bx(3,1) = c(1,3,1) * xrr
         bx(3,3) = c(1,3,3) * xrr3
      end if
c
c     first order T2 matrix elements
c
      if (mdqsiz .ge. 1) then
         t2(2,1) = t2matrix (0,0,1,r2,bx,by)
         t2(3,1) = t2matrix (0,1,0,r2,bx,by)
         t2(4,1) = t2matrix (1,0,0,r2,bx,by)
      end if
c
c     second order T2 matrix elements
c
      if (mdqsiz .ge. 4) then
         t2(2,2) = -t2matrix (0,0,2,r3,bx,by)
         t2(3,2) = -t2matrix (0,1,1,r3,bx,by)
         t2(4,2) = -t2matrix (1,0,1,r3,bx,by)
         t2(3,3) = -t2matrix (0,2,0,r3,bx,by)
         t2(4,3) = -t2matrix (1,1,0,r3,bx,by)
         t2(4,4) = -t2(2,2) - t2(3,3)
         t2(2,3) = t2(3,2)
         t2(2,4) = t2(4,2)
         t2(3,4) = t2(4,3)
      end if
c
c     third order T2 matrix elements
c
      if (mdqsiz .ge. 13) then
         t2(2,5) = t2matrix (0,0,3,r4,bx,by)
         t2(3,5) = t2matrix (0,1,2,r4,bx,by)
         t2(4,5) = t2matrix (1,0,2,r4,bx,by)
         t2(3,6) = t2matrix (0,2,1,r4,bx,by)
         t2(4,6) = t2matrix (1,1,1,r4,bx,by)
         t2(4,7) = -t2(2,5) - t2(3,6)
         t2(2,6) = t2(3,5)
         t2(2,7) = t2(4,5)
         t2(3,7) = t2(4,6)
         t2(2,8) = t2(2,6)
         t2(2,9) = t2(3,6)
         t2(2,10) = t2(4,6)
         t2(3,9) = t2matrix (0,3,0,r4,bx,by)
         t2(3,10) = t2matrix (1,2,0,r4,bx,by)
         t2(4,10) = -t2(2,8) - t2(3,9)
         t2(3,8) = t2(2,9)
         t2(4,8) = t2(2,10)
         t2(4,9) = t2(3,10)
         t2(2,11) = t2(2,7)
         t2(2,12) = t2(3,7)
         t2(2,13) = t2(4,7)
         t2(3,12) = t2(3,10)
         t2(3,13) = t2(4,10)
         t2(4,13) = -t2(2,11) - t2(3,12)
         t2(3,11) = t2(2,12)
         t2(4,11) = t2(2,13)
         t2(4,12) = t2(3,13)
      end if
c
c     compute energy of the induced dipole and multipole sites
c
      do i = 1, mdqsiz
         m2t2(i) = 0.0d0
         do j = 2, 4
            m2t2(i) = m2t2(i) + indk(j-1)*t2(j,i)
         end do
      end do
      epolik = 0.0d0
      do i = 1, mdqsiz
         epolik = epolik + m2t2(i)*rpi(i)
      end do
      epolik = 0.5d0 * (electric / dielec) * epolik
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epolar1  --  polarization energy & derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epolar1" calculates the induced dipole polarization energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine epolar1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'inter.i'
      include 'molcul.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'usage.i'
      include 'virial.i'
      integer i,j,k,skip(maxatm)
      integer ii,kk,iz,ix
      real*8 e,de,epolik1,dutu
      real*8 xi,yi,zi,xr,yr,zr
      real*8 xiz,yiz,ziz,xix,yix,zix
      real*8 taper,dtaper
      real*8 r,r2,r3,r4,r5
      real*8 a(3,3),d(3,3,3,3)
      real*8 rpi(13),indi(3),indk(3)
      real*8 d1(3),d2(3,3),d3(3),utu
      logical iuse,kuse
c
c
c     zero out the induced dipole energy and derivatives
c
      ep = 0.0d0
      do i = 1, n
         do j = 1, 3
            dep(j,i) = 0.0d0
         end do
      end do
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the switching function coefficients
c
      call switch ('CHGDPL')
c
c     rotate multipole components and get rotation derivatives
c
      do i = 1, npole
         call rotmat (i,a)
         call rotsite (i,a)
         do j = 1, 3
            do k = 1, 3
               call drotmat (i,d)
               call drotpole (i,a,d,j,k)
            end do
         end do
      end do
c
c     compute the induced dipoles at each atom
c
      call induce
c
c     calculate the induced dipole energy and derivatives
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         xi = x(i)
         yi = y(i)
         zi = z(i)
         skip(i) = i
         do j = 1, n12(i)
            skip(i12(j,i)) = i
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i
         end do
         do j = 1, mdqsiz
            rpi(j) = rpole(j,ii)
         end do
         indi(1) = uind(1,i)
         indi(2) = uind(2,i)
         indi(3) = uind(3,i)
         do kk = 1, npolar
            k = ipolar(kk)
            kuse = use(k)
            if (iuse .or. kuse) then
               if (skip(k) .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  if (use_image)  call image (xr,yr,zr,0)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     indk(1) = uind(1,k)
                     indk(2) = uind(2,k)
                     indk(3) = uind(3,k)
                     r = sqrt(r2)
                     e = epolik1 (ii,r,xr,yr,zr,rpi,indi,indk,
     &                                   d1,d2,d3,utu)
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                              + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                        de = 2.0d0 * e * dtaper
                        dutu = utu * dtaper
                        e = e * taper
                        d1(1) = d1(1)*taper + de*(xr/r)
                        d1(2) = d1(2)*taper + de*(yr/r)
                        d1(3) = d1(3)*taper + de*(zr/r)
                        d2(1,1) = d2(1,1) * taper
                        d2(1,2) = d2(1,2) * taper
                        d2(1,3) = d2(1,3) * taper
                        d2(2,1) = d2(2,1) * taper
                        d2(2,2) = d2(2,2) * taper
                        d2(2,3) = d2(2,3) * taper
                        d2(3,1) = d2(3,1) * taper
                        d2(3,2) = d2(3,2) * taper
                        d2(3,3) = d2(3,3) * taper
                        d3(1) = d3(1)*taper + dutu*(xr/r)
                        d3(2) = d3(2)*taper + dutu*(yr/r)
                        d3(3) = d3(3)*taper + dutu*(zr/r)
                     end if
                     ep = ep + e
                     dep(1,i) = dep(1,i) - d1(1) - d2(1,1)
                     dep(2,i) = dep(2,i) - d1(2) - d2(1,2)
                     dep(3,i) = dep(3,i) - d1(3) - d2(1,3)
                     if (poltyp .ne. 'DIRECT') then
                        dep(1,i) = dep(1,i) - d3(1)
                        dep(2,i) = dep(2,i) - d3(2)
                        dep(3,i) = dep(3,i) - d3(3)
                     end if
                     dep(1,iz) = dep(1,iz) - d2(2,1)
                     dep(2,iz) = dep(2,iz) - d2(2,2)
                     dep(3,iz) = dep(3,iz) - d2(2,3)
                     dep(1,ix) = dep(1,ix) - d2(3,1)
                     dep(2,ix) = dep(2,ix) - d2(3,2)
                     dep(3,ix) = dep(3,ix) - d2(3,3)
                     dep(1,k) = dep(1,k) + d1(1)
                     dep(2,k) = dep(2,k) + d1(2)
                     dep(3,k) = dep(3,k) + d1(3)
c
c     increment the total intermolecular energy
c
                     if (molcule(i) .ne. molcule(k)) then
                        einter = einter + e
                     end if
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        xiz = x(i) - x(iz)
                        yiz = y(i) - y(iz)
                        ziz = z(i) - z(iz)
                        xix = x(i) - x(ix)
                        yix = y(i) - y(ix)
                        zix = z(i) - z(ix)
                        virx = virx + xr*d1(1) + xiz*d2(2,1)
     &                            + xix*d2(3,1)
                        viry = viry + yr*d1(2) + yiz*d2(2,2)
     &                            + yix*d2(3,2)
                        virz = virz + zr*d1(3) + ziz*d2(2,3)
     &                            + zix*d2(3,3)
                        if (poltyp .ne. 'DIRECT') then
                           virx = virx + 0.5d0*xr*d3(1)
                           viry = viry + 0.5d0*yr*d3(2)
                           virz = virz + 0.5d0*zr*d3(3)
                        end if
                     end if
                  end if
               end if
            end if
         end do
      end do
      return
      end
c
c
c     ########################
c     ##                    ##
c     ##  function epolik1  ##
c     ##                    ##
c     ########################
c
c
      function epolik1 (ii,r,xr,yr,zr,rpi,indi,indk,d1,d2,d3,utu)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'electr.i'
      include 'mpole.i'
      include 'units.i'
      integer i,j,k,m,ii
      real*8 epolik1,t2matrix
      real*8 r,r2,r3,r4,r5
      real*8 xr,yr,zr,zrr,factor
      real*8 xrr,xrr2,xrr3,xrr4
      real*8 yrr,yrr2,yrr3,yrr4
      real*8 t3x(3,13),t3y(3,13),t3z(3,13),tt2(3,13)
      real*8 rpi(13),indi(3),indk(3)
      real*8 m2t2(13),t2(13,13)
      real*8 interx(3),intery(3),interz(3),m1t(13)
      real*8 bx(0:5,0:5),by(11,0:5,0:1),bz(11,0:1)
      real*8 d1(3),d2(3,3),d3(3),utu
c
c
c     zeroth order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 1) then
         bz(1,0) = c(1,0,0)
         by(1,0,0) = c(1,0,0) * bz(1,0)
         bx(0,0) = c(1,0,0)
c
c     first order coefficients to speed the T2 calculation
c
         r2 = r * r
         zrr = zr / r
         bz(3,0) = c(3,0,0)
         bz(1,1) = c(1,1,1) * zrr
         yrr = yr / r
         by(3,0,0) = c(3,0,0) * bz(3,0)
         by(1,0,1) = c(1,0,0) * bz(1,1)
         by(1,1,0) = c(1,1,1) * bz(3,0) * yrr
         xrr = xr / r
         bx(1,1) = c(1,1,1) * xrr
c
c     second order coefficients to speed the T2 calculation
c
         r3 = r2 * r
         bz(5,0) = c(5,0,0)
         bz(3,1) = c(3,1,1) * zrr
         yrr2 = yrr * yrr
         by(5,0,0) = c(5,0,0) * bz(5,0)
         by(3,0,1) = c(3,0,0) * bz(3,1)
         by(3,1,0) = c(3,1,1) * bz(5,0) * yrr
         by(1,1,1) = c(1,1,1) * bz(3,1) * yrr
         by(1,2,0) = c(1,2,0)*bz(3,0) + c(1,2,2)*bz(5,0)*yrr2
         xrr2 = xrr * xrr
         bx(2,0) = c(1,2,0)
         bx(2,2) = c(1,2,2) * xrr2
      end if
c
c     third order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 4) then
         r4 = r2 * r2
         bz(7,0) = c(7,0,0)
         bz(5,1) = c(5,1,1) * zrr
         yrr3 = yrr2 * yrr
         by(7,0,0) = c(7,0,0)*bz(7,0)
         by(5,0,1) = c(5,0,0)*bz(5,1)
         by(5,1,0) = c(5,1,1)*bz(7,0)*yrr
         by(3,1,1) = c(3,1,1)*bz(5,1)*yrr
         by(3,2,0) = c(3,2,0)*bz(5,0) + c(3,2,2)*bz(7,0)*yrr2
         by(1,2,1) = c(1,2,0)*bz(3,1) + c(1,2,2)*bz(5,1)*yrr2
         by(1,3,0) = c(1,3,1)*bz(5,0)*yrr + c(1,3,3)*bz(7,0)*yrr3
         xrr3 = xrr2 * xrr
         bx(3,1) = c(1,3,1) * xrr
         bx(3,3) = c(1,3,3) * xrr3
      end if
c
c     fourth order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 13) then
         r5 = r2 * r3
         bz(9,0) = c(9,0,0)
         bz(7,1) = c(7,1,1) * zrr
         yrr4 = yrr2 * yrr2
         by(9,0,0) = c(9,0,0)*bz(9,0)
         by(7,0,1) = c(7,0,0)*bz(7,1)
         by(7,1,0) = c(7,1,1)*bz(9,0)*yrr
         by(5,1,1) = c(5,1,1)*bz(7,1)*yrr
         by(5,2,0) = c(5,2,0)*bz(7,0) + c(5,2,2)*bz(9,0)*yrr2
         by(3,2,1) = c(3,2,0)*bz(5,1) + c(3,2,2)*bz(7,1)*yrr2
         by(3,3,0) = c(3,3,1)*bz(7,0)*yrr + c(3,3,3)*bz(9,0)*yrr3
         by(1,3,1) = c(1,3,1)*bz(5,1)*yrr + c(1,3,3)*bz(7,1)*yrr3
         by(1,4,0) = c(1,4,0)*bz(5,0) + c(1,4,2)*bz(7,0)*yrr2
     &                  + c(1,4,4)*bz(9,0)*yrr4
         xrr4 = xrr2 * xrr2
         bx(4,0) = c(1,4,0)
         bx(4,2) = c(1,4,2) * xrr2
         bx(4,4) = c(1,4,4) * xrr4
      end if
c
c     zeroth order T2 matrix elements
c
      if (mdqsiz .ge. 1) then
         t2(1,1) = t2matrix (0,0,0,r,bx,by)
c
c     first order T2 matrix elements
c
         t2(2,1) = t2matrix (0,0,1,r2,bx,by)
         t2(3,1) = t2matrix (0,1,0,r2,bx,by)
         t2(4,1) = t2matrix (1,0,0,r2,bx,by)
         t2(1,2) = -t2(2,1)
         t2(1,3) = -t2(3,1)
         t2(1,4) = -t2(4,1)
c
c     second order T2 matrix elements
c
         t2(2,2) = -t2matrix (0,0,2,r3,bx,by)
         t2(3,2) = -t2matrix (0,1,1,r3,bx,by)
         t2(4,2) = -t2matrix (1,0,1,r3,bx,by)
         t2(3,3) = -t2matrix (0,2,0,r3,bx,by)
         t2(4,3) = -t2matrix (1,1,0,r3,bx,by)
         t2(4,4) = -t2(2,2) - t2(3,3)
         t2(2,3) = t2(3,2)
         t2(2,4) = t2(4,2)
         t2(3,4) = t2(4,3)
         k = 5
         do i = 2, 4
            do j = 2, 4
               t2(1,k) = -t2(i,j)
               t2(k,1) = t2(1,k)
               k = k + 1
            end do
         end do
      end if
c
c     third order T2 matrix elements
c
      if (mdqsiz .ge. 4) then
         t2(2,5) = t2matrix (0,0,3,r4,bx,by)
         t2(3,5) = t2matrix (0,1,2,r4,bx,by)
         t2(4,5) = t2matrix (1,0,2,r4,bx,by)
         t2(3,6) = t2matrix (0,2,1,r4,bx,by)
         t2(4,6) = t2matrix (1,1,1,r4,bx,by)
         t2(4,7) = -t2(2,5) - t2(3,6)
         t2(2,6) = t2(3,5)
         t2(2,7) = t2(4,5)
         t2(3,7) = t2(4,6)
         t2(2,8) = t2(2,6)
         t2(2,9) = t2(3,6)
         t2(2,10) = t2(4,6)
         t2(3,9) = t2matrix (0,3,0,r4,bx,by)
         t2(3,10) = t2matrix (1,2,0,r4,bx,by)
         t2(4,10) = -t2(2,8) - t2(3,9)
         t2(3,8) = t2(2,9)
         t2(4,8) = t2(2,10)
         t2(4,9) = t2(3,10)
         t2(2,11) = t2(2,7)
         t2(2,12) = t2(3,7)
         t2(2,13) = t2(4,7)
         t2(3,12) = t2(3,10)
         t2(3,13) = t2(4,10)
         t2(4,13) = -t2(2,11) - t2(3,12)
         t2(3,11) = t2(2,12)
         t2(4,11) = t2(2,13)
         t2(4,12) = t2(3,13)
         do i = 5, 13
            do j = 2, 4
               t2(i,j) = -t2(j,i)
            end do
         end do
      end if
c
c     fourth order T2 matrix elements
c
      if (mdqsiz .ge.13) then
         t2(5,5) = t2matrix (0,0,4,r5,bx,by)
         t2(6,5) = t2matrix (0,1,3,r5,bx,by)
         t2(7,5) = t2matrix (1,0,3,r5,bx,by)
         t2(6,6) = t2matrix (0,2,2,r5,bx,by)
         t2(7,6) = t2matrix (1,1,2,r5,bx,by)
         t2(7,7) = -t2(5,5) - t2(6,6)
         t2(8,5) = t2(6,5)
         t2(9,5) = t2(6,6)
         t2(10,5) = t2(7,6)
         t2(9,6) = t2matrix (0,3,1,r5,bx,by)
         t2(10,6) = t2matrix (1,2,1,r5,bx,by)
         t2(10,7) = -t2(8,5) - t2(9,6)
         t2(8,6) = t2(9,5)
         t2(8,7) = t2(10,5)
         t2(9,7) = t2(10,6)
         t2(11,5) = t2(7,5)
         t2(12,5) = t2(7,6)
         t2(13,5) = t2(7,7)
         t2(12,6) = t2(10,6)
         t2(13,6) = t2(10,7)
         t2(13,7) = -t2(11,5) - t2(12,6)
         t2(11,6) = t2(12,5)
         t2(11,7) = t2(13,5)
         t2(12,7) = t2(13,6)
         t2(8,8) = t2(8,6)
         t2(9,8) = t2(9,6)
         t2(10,8) = t2(10,6)
         t2(9,9) = t2matrix (0,4,0,r5,bx,by)
         t2(10,9) = t2matrix (1,3,0,r5,bx,by)
         t2(10,10) = -t2(8,8) - t2(9,9)
         t2(11,8) = t2(11,6)
         t2(12,8) = t2(12,6)
         t2(13,8) = t2(13,6)
         t2(12,9) = t2(10,9)
         t2(13,9) = t2(10,10)
         t2(13,10) = -t2(11,8) - t2(12,9)
         t2(11,9) = t2(12,8)
         t2(11,10) = t2(13,8)
         t2(12,10) = t2(13,9)
         t2(11,11) = t2(11,7)
         t2(12,11) = t2(12,7)
         t2(13,11) = t2(13,7)
         t2(12,12) = t2(12,10)
         t2(13,12) = t2(13,10)
         t2(13,13) = -t2(11,11) - t2(12,12)
         do i = 5, 12
            do j = i+1, 13
               t2(i,j) = t2(j,i)
            end do
         end do
      end if
c
c     set some additional T matrix elements
c
      k = 1
      m = 1
      do i = 5, 13
         do j = 1, 13
            if (k .eq. 1) then
               t3x(m,j) = t2(i,j)
            else if (k .eq. 2) then
               t3y(m,j) = t2(i,j)
            else if (k .eq. 3) then
               t3z(m,j) = t2(i,j)
            end if
         end do
         k = k + 1
         if (k .eq. 4) then
            k = 1
            m = m + 1
         end if
      end do
      do i = 1, 13
         do j = 2, 4
            tt2(j-1,i) = t2(i,j)
         end do
      end do
      do i = 1, 3
         do j = 2, 4
            tt2(i,j) = -tt2(i,j)
         end do
      end do
c
c     compute the induced dipole polarization energy
c
      factor = electric / dielec
      do i = 1, mdqsiz
         m2t2(i) = 0.0d0
         do j = 2, 4
            m2t2(i) = m2t2(i) + indk(j-1)*t2(j,i)
         end do
      end do
      epolik1 = 0.0d0
      do i = 1, mdqsiz
         epolik1 = epolik1 + m2t2(i)*rpi(i)
      end do
      epolik1 = 0.5d0 * factor * epolik1
c
c     compute the d1 components for the induced dipole gradient
c
      do i = 1, 3
         interx(i) = 0.0d0
         intery(i) = 0.0d0
         interz(i) = 0.0d0
         do j = 1, mdqsiz
            interx(i) = interx(i) + rpi(j)*t3x(i,j)
            intery(i) = intery(i) + rpi(j)*t3y(i,j)
            interz(i) = interz(i) + rpi(j)*t3z(i,j)
         end do
      end do
      do i = 1, 3
         d1(i) = 0.0d0
      end do
      do i = 1, 3
         d1(1) = d1(1) + interx(i)*indk(i)
         d1(2) = d1(2) + intery(i)*indk(i)
         d1(3) = d1(3) + interz(i)*indk(i)
      end do
      do i = 1, 3
         d1(i) = factor * d1(i)
      end do
c
c     compute the d2 components for the induced dipole gradient
c
      do i = 1, mdqsiz
         m1t(i) = 0.0d0
         do j = 1, 3
            m1t(i) = m1t(i) + indk(j)*tt2(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            d2(i,j) = 0.0d0
            do k = 1, mdqsiz
               d2(i,j) = d2(i,j) + m1t(k)*dpole(k,i,j,ii)
            end do
            d2(i,j) = factor * d2(i,j)
         end do
      end do
c
c     compute the d3 components for the induced dipole gradient
c
      do i = 1, 3
         interx(i) = 0.0d0
         intery(i) = 0.0d0
         interz(i) = 0.0d0
         do j = 2, 4
            interx(i) = interx(i) + indk(j-1)*t3x(i,j)
            intery(i) = intery(i) + indk(j-1)*t3y(i,j)
            interz(i) = interz(i) + indk(j-1)*t3z(i,j)
         end do
      end do
      do i = 1, 3
         d3(i) = 0.0d0
      end do
      do i = 1, 3
         d3(1) = d3(1) + interx(i)*indi(i)
         d3(2) = d3(2) + intery(i)*indi(i)
         d3(3) = d3(3) + interz(i)*indi(i)
      end do
      do i = 1, 3
         d3(i) = factor * d3(i)
      end do
c
c     compute the utu component for the induced dipole gradient
c
      utu = 0.0d0
      do i = 1, 3
         utu = utu + m2t2(i+1)*indi(i)
      end do
      utu = factor * utu
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epolar2  --  atom-by-atom polarization Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epolar2" calculates second derivatives of the induced
c     dipole polarization energy for a single atom at a time
c
c
      subroutine epolar2 (i)
      implicit none
      integer i
c
c
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar3  --  polarization energy and analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine epolar3
      implicit none
      include 'sizes.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'electr.i'
      include 'energi.i'
      include 'inform.i'
      include 'iounit.i'
      include 'moment.i'
      include 'mpole.i'
      include 'polar.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k,skip(maxatm)
      integer ii,kk,iz,ix
      real*8 e,epolik,a(3,3)
      real*8 xi,yi,zi,xr,yr,zr
      real*8 r,r2,r3,r4,r5,taper
      real*8 rpi(13),indk(3)
      real*8 utotal,xsum,ysum,zsum
      logical iuse,kuse,huge,header
c
c
c     zero out the induced dipole energy and partitioning
c
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
      end do
      nplr = 0
      header = .true.
c
c     zero out the list of atoms to be skipped
c
      do i = 1, n
         skip(i) = 0
      end do
c
c     set the switching function coefficients
c
      call switch ('CHGDPL')
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each atom
c
      call induce
c
c     add induced dipole components to the permanent components
c
      xsum = 0.0d0
      ysum = 0.0d0
      zsum = 0.0d0
      do ii = 1, npolar
         i = ipolar(ii)
         xsum = xsum + uind(1,i)
         ysum = ysum + uind(2,i)
         zsum = zsum + uind(3,i)
      end do
      xdipole = xdipole + debye*xsum
      ydipole = ydipole + debye*ysum
      zdipole = zdipole + debye*zsum
c
c     compute and partition the induced dipole energy
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iuse = (use(i) .or. use(iz) .or. use(ix))
         xi = x(i)
         yi = y(i)
         zi = z(i)
         skip(i) = i
         do j = 1, n12(i)
            skip(i12(j,i)) = i
         end do
         do j = 1, n13(i)
            skip(i13(j,i)) = i
         end do
         do j = 1, mdqsiz
            rpi(j) = rpole(j,ii)
         end do
         do kk = 1, npolar
            k = ipolar(kk)
            kuse = use(k)
            if (iuse .or. kuse) then
               if (skip(k) .ne. i) then
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  if (use_image)  call image (xr,yr,zr,0)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     indk(1) = uind(1,k)
                     indk(2) = uind(2,k)
                     indk(3) = uind(3,k)
                     r = sqrt(r2)
                     e = epolik (r,xr,yr,zr,rpi,indk)
                     if (r2 .gt. cut2) then
                        r3 = r2 * r
                        r4 = r2 * r2
                        r5 = r2 * r3
                        taper = c5*r5 + c4*r4 + c3*r3
     &                             + c2*r2 + c1*r + c0
                        e = e * taper
                     end if
                     nplr = nplr + 1
                     ep = ep + e
                     aep(i) = aep(i) + 0.5d0*e
                     aep(k) = aep(k) + 0.5d0*e
c
c     print a warning if the energy of dipole interaction is large
c
                     huge = (abs(e) .gt. 1.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,10)
   10                      format (/,' Individual Induced Dipole',
     &                                ' Interactions :',
     &                             //,' Type',11x,'Atom Names',
     &                                13x,'IndDpl - Chg',3x,'Distance',
     &                                5x,'Energy',/)
                        end if
                        utotal = sqrt(uind(1,k)**2 + uind(2,k)**2
     &                                  + uind(3,k)**2) * debye
                        write (iout,20)  k,name(k),i,name(i),
     &                                   utotal,pole(1,ii),r,e
   20                   format (' Polarize ',i5,'-',a3,1x,i5,'-',
     &                             a3,8x,2f7.2,f10.4,f12.4)
                     end if
                  end if
               end if
            end if
         end do
      end do
      return
      end

c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole1  --  mpole/polar energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1" calculates the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c
c
      subroutine empole1
      use sizes
      use deriv
      use energi
      use limits
      use mpole
      use potent
      implicit none
      integer i,j,ii
c
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call empole1d
         else
            call empole1c
         end if
      else
         if (use_mlist) then
            call empole1b
         else
            call empole1a
         end if
      end if
c
c     zero out energy and derivative terms which are not in use
c
      if (.not. use_mpole) then
         em = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            do j = 1, 3
               dem(j,ii) = 0.0d0
            end do
         end do
      end if
      if (.not. use_polar) then
         ep = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            do j = 1, 3
               dep(j,ii) = 0.0d0
            end do
         end do
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole1a  --  double loop multipole derivatives  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole1a" calculates the multipole energy and derivatives with 
c     respect to Cartesian coordinates using a pairwise double loop
c
c
      subroutine empole1a
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use molcul
      use mplpot
      use mpole
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer ix,iy,iz
      integer kx,ky,kz
      integer iax,iay,iaz
      real*8 e,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qirx,qiry,qirz
      real*8 qkrx,qkry,qkrz
      real*8 qirxr,qiryr,qirzr
      real*8 qkrxr,qkryr,qkrzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 sc(9),ge(5),gf(7)
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
c
c     set arrays needed for connected atom scaling and torque
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the multipole interaction energy and gradient
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i+1, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 10
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f * mscale(kk) / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     construct several necessary additional variables
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
               qirx = qixx*xr + qixy*yr + qixz*zr
               qiry = qixy*xr + qiyy*yr + qiyz*zr
               qirz = qixz*xr + qiyz*yr + qizz*zr
               qkrx = qkxx*xr + qkxy*yr + qkxz*zr
               qkry = qkxy*xr + qkyy*yr + qkyz*zr
               qkrz = qkxz*xr + qkyz*yr + qkzz*zr
               qirxr = qirz*yr - qiry*zr
               qiryr = qirx*zr - qirz*xr
               qirzr = qiry*xr - qirx*yr
               qkrxr = qkrz*yr - qkry*zr
               qkryr = qkrx*zr - qkrz*xr
               qkrzr = qkry*xr - qkrx*yr
               qrrx = qkry*qirz - qkrz*qiry
               qrry = qkrz*qirx - qkrx*qirz
               qrrz = qkrx*qiry - qkry*qirx
               qikrx = qixx*qkrx + qixy*qkry + qixz*qkrz
               qikry = qixy*qkrx + qiyy*qkry + qiyz*qkrz
               qikrz = qixz*qkrx + qiyz*qkry + qizz*qkrz
               qkirx = qkxx*qirx + qkxy*qiry + qkxz*qirz
               qkiry = qkxy*qirx + qkyy*qiry + qkyz*qirz
               qkirz = qkxz*qirx + qkyz*qiry + qkzz*qirz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qkrz - diz*qkry + dky*qirz - dkz*qiry
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qkrx - dix*qkrz + dkz*qirx - dkx*qirz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qkry - diy*qkrx + dkx*qiry - dky*qirx
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate scalar products for multipole interactions
c
               sc(1) = dix*dkx + diy*dky + diz*dkz
               sc(2) = dix*xr + diy*yr + diz*zr
               sc(3) = dkx*xr + dky*yr + dkz*zr
               sc(4) = qirx*xr + qiry*yr + qirz*zr
               sc(5) = qkrx*xr + qkry*yr + qkrz*zr
               sc(6) = dkx*qirx + dky*qiry + dkz*qirz
               sc(7) = dix*qkrx + diy*qkry + diz*qkrz
               sc(8) = qirx*qkrx + qiry*qkry + qirz*qkrz
               sc(9) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                    + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     construct auxiliary variables for multipole energy
c
               ge(1) = ci*ck
               ge(2) = ck*sc(2) - ci*sc(3) + sc(1)
               ge(3) = ci*sc(5) + ck*sc(4) - sc(2)*sc(3)
     &                    + 2.0d0*(sc(6)-sc(7)+sc(9))
               ge(4) = sc(2)*sc(5) - sc(3)*sc(4) - 4.0d0*sc(8)
               ge(5) = sc(4)*sc(5)
c
c     compute the energy contribution for this interaction
c
               e = rr1*ge(1) + rr3*ge(2) + rr5*ge(3)
     &                + rr7*ge(4) + rr9*ge(5)
               if (ii .eq. kk)  e = 0.5d0 * e
               em = em + e
c
c     increment the total intramolecular energy
c
               if (molcule(ii) .ne. molcule(kk))
     &            einter = einter + e
c
c     construct auxiliary variables for force and torque
c
               gf(1) = rr3*ge(1) + rr5*ge(2) + rr7*ge(3)
     &                    + rr9*ge(4) + rr11*ge(5)
               gf(2) = -ck*rr3 + sc(3)*rr5 - sc(5)*rr7
               gf(3) = ci*rr3 + sc(2)*rr5 + sc(4)*rr7
               gf(4) = 2.0d0 * rr5
               gf(5) = 2.0d0 * (-ck*rr5+sc(3)*rr7-sc(5)*rr9)
               gf(6) = 2.0d0 * (-ci*rr5-sc(2)*rr7-sc(4)*rr9)
               gf(7) = 4.0d0 * rr7
c
c     compute the force components for both sites
c
               frcx = gf(1)*xr + gf(2)*dix + gf(3)*dkx
     &                   + gf(4)*(diqkx-dkqix) + gf(5)*qirx
     &                   + gf(6)*qkrx + gf(7)*(qikrx+qkirx)
               frcy = gf(1)*yr + gf(2)*diy + gf(3)*dky
     &                   + gf(4)*(diqky-dkqiy) + gf(5)*qiry
     &                   + gf(6)*qkry + gf(7)*(qikry+qkiry)
               frcz = gf(1)*zr + gf(2)*diz + gf(3)*dkz
     &                   + gf(4)*(diqkz-dkqiz) + gf(5)*qirz
     &                   + gf(6)*qkrz + gf(7)*(qikrz+qkirz)
c
c     compute the torque components for both sites
c
               ttmi(1) = -rr3*dikx + gf(2)*dirx
     &                      + gf(4)*(dqiqkx+dkqixr)
     &                      - gf(5)*qirxr - gf(7)*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + gf(2)*diry
     &                      + gf(4)*(dqiqky+dkqiyr)
     &                      - gf(5)*qiryr - gf(7)*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + gf(2)*dirz
     &                      + gf(4)*(dqiqkz+dkqizr)
     &                      - gf(5)*qirzr - gf(7)*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + gf(3)*dkrx
     &                      - gf(4)*(dqiqkx+diqkxr)
     &                      - gf(6)*qkrxr - gf(7)*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + gf(3)*dkry
     &                      - gf(4)*(dqiqky+diqkyr)
     &                      - gf(6)*qkryr - gf(7)*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + gf(3)*dkrz
     &                      - gf(4)*(dqiqkz+diqkzr)
     &                      - gf(6)*qkrzr - gf(7)*(qkirzr-qrrz)
c
c     force and torque components scaled by group membership
c
               if (use_group) then
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
                  do j = 1, 3
                     ttmi(j) = fgrp * ttmi(j)
                     ttmk(j) = fgrp * ttmk(j)
                  end do
               end if
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
   10       continue
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction with other unit cells
c
      do i = 1, npole
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 20
            do jcell = 1, ncell
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               mscale(kk) = 1.0d0
            end if
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f * mscale(kk) / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     construct several necessary additional variables
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
               qirx = qixx*xr + qixy*yr + qixz*zr
               qiry = qixy*xr + qiyy*yr + qiyz*zr
               qirz = qixz*xr + qiyz*yr + qizz*zr
               qkrx = qkxx*xr + qkxy*yr + qkxz*zr
               qkry = qkxy*xr + qkyy*yr + qkyz*zr
               qkrz = qkxz*xr + qkyz*yr + qkzz*zr
               qirxr = qirz*yr - qiry*zr
               qiryr = qirx*zr - qirz*xr
               qirzr = qiry*xr - qirx*yr
               qkrxr = qkrz*yr - qkry*zr
               qkryr = qkrx*zr - qkrz*xr
               qkrzr = qkry*xr - qkrx*yr
               qrrx = qkry*qirz - qkrz*qiry
               qrry = qkrz*qirx - qkrx*qirz
               qrrz = qkrx*qiry - qkry*qirx
               qikrx = qixx*qkrx + qixy*qkry + qixz*qkrz
               qikry = qixy*qkrx + qiyy*qkry + qiyz*qkrz
               qikrz = qixz*qkrx + qiyz*qkry + qizz*qkrz
               qkirx = qkxx*qirx + qkxy*qiry + qkxz*qirz
               qkiry = qkxy*qirx + qkyy*qiry + qkyz*qirz
               qkirz = qkxz*qirx + qkyz*qiry + qkzz*qirz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qkrz - diz*qkry + dky*qirz - dkz*qiry
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qkrx - dix*qkrz + dkz*qirx - dkx*qirz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qkry - diy*qkrx + dkx*qiry - dky*qirx
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate scalar products for multipole interactions
c
               sc(1) = dix*dkx + diy*dky + diz*dkz
               sc(2) = dix*xr + diy*yr + diz*zr
               sc(3) = dkx*xr + dky*yr + dkz*zr
               sc(4) = qirx*xr + qiry*yr + qirz*zr
               sc(5) = qkrx*xr + qkry*yr + qkrz*zr
               sc(6) = dkx*qirx + dky*qiry + dkz*qirz
               sc(7) = dix*qkrx + diy*qkry + diz*qkrz
               sc(8) = qirx*qkrx + qiry*qkry + qirz*qkrz
               sc(9) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                    + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     construct auxiliary variables for multipole energy
c
               ge(1) = ci*ck
               ge(2) = ck*sc(2) - ci*sc(3) + sc(1)
               ge(3) = ci*sc(5) + ck*sc(4) - sc(2)*sc(3)
     &                    + 2.0d0*(sc(6)-sc(7)+sc(9))
               ge(4) = sc(2)*sc(5) - sc(3)*sc(4) - 4.0d0*sc(8)
               ge(5) = sc(4)*sc(5)
c
c     compute the energy contribution for this interaction
c
               e = rr1*ge(1) + rr3*ge(2) + rr5*ge(3)
     &                + rr7*ge(4) + rr9*ge(5)
               if (ii .eq. kk)  e = 0.5d0 * e
               em = em + e
c
c     increment the total intramolecular energy
c
               einter = einter + e
c
c     construct auxiliary variables for force and torque
c
               gf(1) = rr3*ge(1) + rr5*ge(2) + rr7*ge(3)
     &                    + rr9*ge(4) + rr11*ge(5)
               gf(2) = -ck*rr3 + sc(3)*rr5 - sc(5)*rr7
               gf(3) = ci*rr3 + sc(2)*rr5 + sc(4)*rr7
               gf(4) = 2.0d0 * rr5
               gf(5) = 2.0d0 * (-ck*rr5+sc(3)*rr7-sc(5)*rr9)
               gf(6) = 2.0d0 * (-ci*rr5-sc(2)*rr7-sc(4)*rr9)
               gf(7) = 4.0d0 * rr7
c
c     compute the force components for both sites
c
               frcx = gf(1)*xr + gf(2)*dix + gf(3)*dkx
     &                   + gf(4)*(diqkx-dkqix) + gf(5)*qirx
     &                   + gf(6)*qkrx + gf(7)*(qikrx+qkirx)
               frcy = gf(1)*yr + gf(2)*diy + gf(3)*dky
     &                   + gf(4)*(diqky-dkqiy) + gf(5)*qiry
     &                   + gf(6)*qkry + gf(7)*(qikry+qkiry)
               frcz = gf(1)*zr + gf(2)*diz + gf(3)*dkz
     &                   + gf(4)*(diqkz-dkqiz) + gf(5)*qirz
     &                   + gf(6)*qkrz + gf(7)*(qikrz+qkirz)
c
c     compute the torque components for both sites
c
               ttmi(1) = -rr3*dikx + gf(2)*dirx
     &                      + gf(4)*(dqiqkx+dkqixr)
     &                      - gf(5)*qirxr - gf(7)*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + gf(2)*diry
     &                      + gf(4)*(dqiqky+dkqiyr)
     &                      - gf(5)*qiryr - gf(7)*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + gf(2)*dirz
     &                      + gf(4)*(dqiqkz+dkqizr)
     &                      - gf(5)*qirzr - gf(7)*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + gf(3)*dkrx
     &                      - gf(4)*(dqiqkx+diqkxr)
     &                      - gf(6)*qkrxr - gf(7)*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + gf(3)*dkry
     &                      - gf(4)*(dqiqky+diqkyr)
     &                      - gf(6)*qkryr - gf(7)*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + gf(3)*dkrz
     &                      - gf(4)*(dqiqkz+diqkzr)
     &                      - gf(6)*qkrzr - gf(7)*(qkirzr-qrrz)
c
c     force and torque scaled for self-interactions and groups
c
               if (ii .eq. kk) then
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  do j = 1, 3
                     ttmi(j) = 0.5d0 * ttmi(j)
                     ttmk(j) = 0.5d0 * ttmk(j)
                  end do
               end if
               if (use_group) then
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
                  do j = 1, 3
                     ttmi(j) = fgrp * ttmi(j)
                     ttmk(j) = fgrp * ttmk(j)
                  end do
               end if
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
            end do
   20       continue
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
      end if
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
         vxz = zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole1b  --  neighbor list multipole derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole1b" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using a neighbor list
c
c
      subroutine empole1b
      use sizes
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use inter
      use molcul
      use mplpot
      use mpole
      use neigh
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      integer iax,iay,iaz
      real*8 e,f,fgrp
      real*8 emo,eintero
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qirx,qiry,qirz
      real*8 qkrx,qkry,qkrz
      real*8 qirxr,qiryr,qirzr
      real*8 qkrxr,qkryr,qkrzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 sc(9),ge(5),gf(7)
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 viro(3,3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: demo(:,:)
      real*8, allocatable :: temo(:,:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (demo(3,n))
      allocate (temo(3,n))
c
c     set arrays needed to scale connected atom interactions
c
      if (npole .eq. 0)  return
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and scaling coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     initialize local variables for OpenMP calculation
c
      emo = 0.0d0
      eintero = einter
      do i = 1, n
         do j = 1, 3
            demo(j,i) = 0.0d0
            temo(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            viro(j,i) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,xaxis,yaxis,zaxis,rpole,use,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,nelst,
!$OMP& elst,use_group,use_intra,use_bounds,off2,f,molcule,emo,eintero,
!$OMP& demo,temo,viro)
!$OMP& firstprivate(mscale)
!$OMP DO reduction(+:emo,eintero,demo,temo,viro) schedule(guided)
c
c     compute the multipole interaction energy and gradient
c
      do i = 1, npole
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 10
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f * mscale(kk) / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     construct several necessary additional variables
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
               qirx = qixx*xr + qixy*yr + qixz*zr
               qiry = qixy*xr + qiyy*yr + qiyz*zr
               qirz = qixz*xr + qiyz*yr + qizz*zr
               qkrx = qkxx*xr + qkxy*yr + qkxz*zr
               qkry = qkxy*xr + qkyy*yr + qkyz*zr
               qkrz = qkxz*xr + qkyz*yr + qkzz*zr
               qirxr = qirz*yr - qiry*zr
               qiryr = qirx*zr - qirz*xr
               qirzr = qiry*xr - qirx*yr
               qkrxr = qkrz*yr - qkry*zr
               qkryr = qkrx*zr - qkrz*xr
               qkrzr = qkry*xr - qkrx*yr
               qrrx = qkry*qirz - qkrz*qiry
               qrry = qkrz*qirx - qkrx*qirz
               qrrz = qkrx*qiry - qkry*qirx
               qikrx = qixx*qkrx + qixy*qkry + qixz*qkrz
               qikry = qixy*qkrx + qiyy*qkry + qiyz*qkrz
               qikrz = qixz*qkrx + qiyz*qkry + qizz*qkrz
               qkirx = qkxx*qirx + qkxy*qiry + qkxz*qirz
               qkiry = qkxy*qirx + qkyy*qiry + qkyz*qirz
               qkirz = qkxz*qirx + qkyz*qiry + qkzz*qirz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qkrz - diz*qkry + dky*qirz - dkz*qiry
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qkrx - dix*qkrz + dkz*qirx - dkx*qirz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qkry - diy*qkrx + dkx*qiry - dky*qirx
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate scalar products for multipole interactions
c
               sc(1) = dix*dkx + diy*dky + diz*dkz
               sc(2) = dix*xr + diy*yr + diz*zr
               sc(3) = dkx*xr + dky*yr + dkz*zr
               sc(4) = qirx*xr + qiry*yr + qirz*zr
               sc(5) = qkrx*xr + qkry*yr + qkrz*zr
               sc(6) = dkx*qirx + dky*qiry + dkz*qirz
               sc(7) = dix*qkrx + diy*qkry + diz*qkrz
               sc(8) = qirx*qkrx + qiry*qkry + qirz*qkrz
               sc(9) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                    + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     construct auxiliary variables for multipole energy
c
               ge(1) = ci*ck
               ge(2) = ck*sc(2) - ci*sc(3) + sc(1)
               ge(3) = ci*sc(5) + ck*sc(4) - sc(2)*sc(3)
     &                    + 2.0d0*(sc(6)-sc(7)+sc(9))
               ge(4) = sc(2)*sc(5) - sc(3)*sc(4) - 4.0d0*sc(8)
               ge(5) = sc(4)*sc(5)
c
c     compute the energy contribution for this interaction
c
               e = rr1*ge(1) + rr3*ge(2) + rr5*ge(3)
     &                + rr7*ge(4) + rr9*ge(5)
               if (use_group)  e = e * fgrp
               emo = emo + e
c
c     increment the total intermolecular energy
c
               if (molcule(ii) .ne. molcule(kk))
     &            eintero = eintero + e
c
c     intermediate variables for the permanent multipoles
c
               gf(1) = rr3*ge(1) + rr5*ge(2) + rr7*ge(3)
     &                    + rr9*ge(4) + rr11*ge(5)
               gf(2) = -ck*rr3 + sc(3)*rr5 - sc(5)*rr7
               gf(3) = ci*rr3 + sc(2)*rr5 + sc(4)*rr7
               gf(4) = 2.0d0 * rr5
               gf(5) = 2.0d0 * (-ck*rr5+sc(3)*rr7-sc(5)*rr9)
               gf(6) = 2.0d0 * (-ci*rr5-sc(2)*rr7-sc(4)*rr9)
               gf(7) = 4.0d0 * rr7
c
c     get the permanent multipole force components
c
               frcx = gf(1)*xr + gf(2)*dix + gf(3)*dkx
     &                   + gf(4)*(diqkx-dkqix) + gf(5)*qirx
     &                   + gf(6)*qkrx + gf(7)*(qikrx+qkirx)
               frcy = gf(1)*yr + gf(2)*diy + gf(3)*dky
     &                   + gf(4)*(diqky-dkqiy) + gf(5)*qiry
     &                   + gf(6)*qkry + gf(7)*(qikry+qkiry)
               frcz = gf(1)*zr + gf(2)*diz + gf(3)*dkz
     &                   + gf(4)*(diqkz-dkqiz) + gf(5)*qirz
     &                   + gf(6)*qkrz + gf(7)*(qikrz+qkirz)
c
c     get the permanent multipole torque components
c
               ttmi(1) = -rr3*dikx + gf(2)*dirx - gf(5)*qirxr
     &                      + gf(4)*(dqiqkx+dkqixr)
     &                      - gf(7)*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + gf(2)*diry - gf(5)*qiryr
     &                      + gf(4)*(dqiqky+dkqiyr)
     &                      - gf(7)*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + gf(2)*dirz - gf(5)*qirzr
     &                      + gf(4)*(dqiqkz+dkqizr)
     &                      - gf(7)*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + gf(3)*dkrx - gf(6)*qkrxr
     &                      - gf(4)*(dqiqkx+diqkxr)
     &                      - gf(7)*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + gf(3)*dkry - gf(6)*qkryr
     &                      - gf(4)*(dqiqky+diqkyr)
     &                      - gf(7)*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + gf(3)*dkrz - gf(6)*qkrzr
     &                      - gf(4)*(dqiqkz+diqkzr)
     &                      - gf(7)*(qkirzr-qrrz)
c
c     increment force-based gradient and torque on first site
c
               demo(1,ii) = demo(1,ii) + frcx
               demo(2,ii) = demo(2,ii) + frcy
               demo(3,ii) = demo(3,ii) + frcz
               temo(1,i) = temo(1,i) + ttmi(1)
               temo(2,i) = temo(2,i) + ttmi(2)
               temo(3,i) = temo(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               demo(1,kk) = demo(1,kk) - frcx
               demo(2,kk) = demo(2,kk) - frcy
               demo(3,kk) = demo(3,kk) - frcz
               temo(1,k) = temo(1,k) + ttmk(1)
               temo(2,k) = temo(2,k) + ttmk(2)
               temo(3,k) = temo(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               viro(1,1) = viro(1,1) + vxx
               viro(2,1) = viro(2,1) + vxy
               viro(3,1) = viro(3,1) + vxz
               viro(1,2) = viro(1,2) + vxy
               viro(2,2) = viro(2,2) + vyy
               viro(3,2) = viro(3,2) + vyz
               viro(1,3) = viro(1,3) + vxz
               viro(2,3) = viro(2,3) + vyz
               viro(3,3) = viro(3,3) + vzz
            end if
   10       continue
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO reduction(+:demo,viro) schedule(guided)
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,temo(1,i),fix,fiy,fiz,demo)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
         vxz = zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         viro(1,1) = viro(1,1) + vxx
         viro(2,1) = viro(2,1) + vxy
         viro(3,1) = viro(3,1) + vxz
         viro(1,2) = viro(1,2) + vxy
         viro(2,2) = viro(2,2) + vyy
         viro(3,2) = viro(3,2) + vyz
         viro(1,3) = viro(1,3) + vxz
         viro(2,3) = viro(2,3) + vyz
         viro(3,3) = viro(3,3) + vzz
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add local to global variables for OpenMP calculation
c
      em = em + emo
      einter = eintero
      do i = 1, n
         do j = 1, 3
            dem(j,i) = dem(j,i) + demo(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + viro(j,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (demo)
      deallocate (temo)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole1c  --  Ewald multipole derivs via loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1c" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh
c     Ewald summation and a double loop
c
c
      subroutine empole1c
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use inter
      use math
      use mpole
      use virial
      implicit none
      integer i,j,ii
      real*8 e,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 trq(3),frcx(3)
      real*8 frcy(3),frcz(3)
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the real space part of the Ewald summation
c
      call emreal1c (eintra)
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip1
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         do i = 1, npole
            ii = ipole(i)
            dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
            dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
            dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         do i = 1, npole
            trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
            trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
            trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
            call torque (i,trq,frcx,frcy,frcz,dem)
         end do
c
c     boundary correction to virial due to overall cell dipole
c
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xq = 0.0d0
         yq = 0.0d0
         zq = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i)
            yd = yd + rpole(3,i)
            zd = zd + rpole(4,i)
            xq = xq + rpole(1,i)*x(ii)
            yq = yq + rpole(1,i)*y(ii)
            zq = zq + rpole(1,i)*z(ii)
         end do
         xv = xd * xq
         yv = yd * yq
         zv = zd * zq
         vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                      + xq*xq + yq*yq + zq*zq)
         vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
         vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
         vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
         vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
         vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
         vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
         vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
         vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
         vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + em - eintra
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1c  --  Ewald real space derivs via loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a double loop
c
c
      subroutine emreal1c (eintra)
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use math
      use molcul
      use mplpot
      use mpole
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer iax,iay,iaz
      real*8 e,f,bfac,erfc
      real*8 eintra,efull
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qirx,qiry,qirz
      real*8 qkrx,qkry,qkrz
      real*8 qirxr,qiryr,qirzr
      real*8 qkrxr,qkryr,qkrzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 sc(9),ge(5),gf(7)
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 bn(0:5)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      character*6 mode
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
c
c     set arrays needed for connected atom scaling and torque
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 5
                  bn(j) = f * bn(j)
               end do
c
c     construct several necessary additional variables
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
               qirx = qixx*xr + qixy*yr + qixz*zr
               qiry = qixy*xr + qiyy*yr + qiyz*zr
               qirz = qixz*xr + qiyz*yr + qizz*zr
               qkrx = qkxx*xr + qkxy*yr + qkxz*zr
               qkry = qkxy*xr + qkyy*yr + qkyz*zr
               qkrz = qkxz*xr + qkyz*yr + qkzz*zr
               qirxr = qirz*yr - qiry*zr
               qiryr = qirx*zr - qirz*xr
               qirzr = qiry*xr - qirx*yr
               qkrxr = qkrz*yr - qkry*zr
               qkryr = qkrx*zr - qkrz*xr
               qkrzr = qkry*xr - qkrx*yr
               qrrx = qkry*qirz - qkrz*qiry
               qrry = qkrz*qirx - qkrx*qirz
               qrrz = qkrx*qiry - qkry*qirx
               qikrx = qixx*qkrx + qixy*qkry + qixz*qkrz
               qikry = qixy*qkrx + qiyy*qkry + qiyz*qkrz
               qikrz = qixz*qkrx + qiyz*qkry + qizz*qkrz
               qkirx = qkxx*qirx + qkxy*qiry + qkxz*qirz
               qkiry = qkxy*qirx + qkyy*qiry + qkyz*qirz
               qkirz = qkxz*qirx + qkyz*qiry + qkzz*qirz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qkrz - diz*qkry + dky*qirz - dkz*qiry
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qkrx - dix*qkrz + dkz*qirx - dkx*qirz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qkry - diy*qkrx + dkx*qiry - dky*qirx
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate scalar products for multipole interactions
c
               sc(1) = dix*dkx + diy*dky + diz*dkz
               sc(2) = dix*xr + diy*yr + diz*zr
               sc(3) = dkx*xr + dky*yr + dkz*zr
               sc(4) = qirx*xr + qiry*yr + qirz*zr
               sc(5) = qkrx*xr + qkry*yr + qkrz*zr
               sc(6) = dkx*qirx + dky*qiry + dkz*qirz
               sc(7) = dix*qkrx + diy*qkry + diz*qkrz
               sc(8) = qirx*qkrx + qiry*qkry + qirz*qkrz
               sc(9) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                    + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     construct auxiliary variables for multipole energy
c
               ge(1) = ci*ck
               ge(2) = ck*sc(2) - ci*sc(3) + sc(1)
               ge(3) = ci*sc(5) + ck*sc(4) - sc(2)*sc(3)
     &                    + 2.0d0*(sc(6)-sc(7)+sc(9))
               ge(4) = sc(2)*sc(5) - sc(3)*sc(4) - 4.0d0*sc(8)
               ge(5) = sc(4)*sc(5)
c
c     compute the full energy without any Ewald scaling
c
               efull = rr1*ge(1) + rr3*ge(2) + rr5*ge(3)
     &                    + rr7*ge(4) + rr9*ge(5)
               efull = efull * mscale(kk)
c
c     modify distances to account for Ewald and exclusions
c
               scalekk = 1.0d0 - mscale(kk)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
               rr11 = bn(5) - scalekk*rr11
c
c     compute the energy contribution for this interaction
c
               e = rr1*ge(1) + rr3*ge(2) + rr5*ge(3)
     &                + rr7*ge(4) + rr9*ge(5)
               em = em + e
c
c     increment the total intramolecular energy
c
               if (molcule(ii) .eq. molcule(kk))
     &            eintra = eintra + efull
c
c     construct auxiliary variables for force and torque
c
               gf(1) = rr3*ge(1) + rr5*ge(2) + rr7*ge(3)
     &                    + rr9*ge(4) + rr11*ge(5)
               gf(2) = -ck*rr3 + sc(3)*rr5 - sc(5)*rr7
               gf(3) = ci*rr3 + sc(2)*rr5 + sc(4)*rr7
               gf(4) = 2.0d0 * rr5
               gf(5) = 2.0d0 * (-ck*rr5+sc(3)*rr7-sc(5)*rr9)
               gf(6) = 2.0d0 * (-ci*rr5-sc(2)*rr7-sc(4)*rr9)
               gf(7) = 4.0d0 * rr7
c
c     compute the force components for both sites
c
               frcx = gf(1)*xr + gf(2)*dix + gf(3)*dkx
     &                   + gf(4)*(diqkx-dkqix) + gf(5)*qirx
     &                   + gf(6)*qkrx + gf(7)*(qikrx+qkirx)
               frcy = gf(1)*yr + gf(2)*diy + gf(3)*dky
     &                   + gf(4)*(diqky-dkqiy) + gf(5)*qiry
     &                   + gf(6)*qkry + gf(7)*(qikry+qkiry)
               frcz = gf(1)*zr + gf(2)*diz + gf(3)*dkz
     &                   + gf(4)*(diqkz-dkqiz) + gf(5)*qirz
     &                   + gf(6)*qkrz + gf(7)*(qikrz+qkirz)
c
c     compute the torque components for both sites
c
               ttmi(1) = -rr3*dikx + gf(2)*dirx
     &                      + gf(4)*(dqiqkx+dkqixr)
     &                      - gf(5)*qirxr - gf(7)*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + gf(2)*diry
     &                      + gf(4)*(dqiqky+dkqiyr)
     &                      - gf(5)*qiryr - gf(7)*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + gf(2)*dirz
     &                      + gf(4)*(dqiqkz+dkqizr)
     &                      - gf(5)*qirzr - gf(7)*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + gf(3)*dkrx
     &                      - gf(4)*(dqiqkx+diqkxr)
     &                      - gf(6)*qkrxr - gf(7)*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + gf(3)*dkry
     &                      - gf(4)*(dqiqky+diqkyr)
     &                      - gf(6)*qkryr - gf(7)*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + gf(3)*dkrz
     &                      - gf(4)*(dqiqkz+diqkzr)
     &                      - gf(6)*qkrzr - gf(7)*(qkirzr-qrrz)
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction with other unit cells
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i, npole
            kk = ipole(k)
            do jcell = 1, ncell
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               mscale(kk) = 1.0d0
            end if
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 5
                  bn(j) = f * bn(j)
               end do
c
c     construct several necessary additional variables
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
               qirx = qixx*xr + qixy*yr + qixz*zr
               qiry = qixy*xr + qiyy*yr + qiyz*zr
               qirz = qixz*xr + qiyz*yr + qizz*zr
               qkrx = qkxx*xr + qkxy*yr + qkxz*zr
               qkry = qkxy*xr + qkyy*yr + qkyz*zr
               qkrz = qkxz*xr + qkyz*yr + qkzz*zr
               qirxr = qirz*yr - qiry*zr
               qiryr = qirx*zr - qirz*xr
               qirzr = qiry*xr - qirx*yr
               qkrxr = qkrz*yr - qkry*zr
               qkryr = qkrx*zr - qkrz*xr
               qkrzr = qkry*xr - qkrx*yr
               qrrx = qkry*qirz - qkrz*qiry
               qrry = qkrz*qirx - qkrx*qirz
               qrrz = qkrx*qiry - qkry*qirx
               qikrx = qixx*qkrx + qixy*qkry + qixz*qkrz
               qikry = qixy*qkrx + qiyy*qkry + qiyz*qkrz
               qikrz = qixz*qkrx + qiyz*qkry + qizz*qkrz
               qkirx = qkxx*qirx + qkxy*qiry + qkxz*qirz
               qkiry = qkxy*qirx + qkyy*qiry + qkyz*qirz
               qkirz = qkxz*qirx + qkyz*qiry + qkzz*qirz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qkrz - diz*qkry + dky*qirz - dkz*qiry
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qkrx - dix*qkrz + dkz*qirx - dkx*qirz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qkry - diy*qkrx + dkx*qiry - dky*qirx
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate scalar products for multipole interactions
c
               sc(1) = dix*dkx + diy*dky + diz*dkz
               sc(2) = dix*xr + diy*yr + diz*zr
               sc(3) = dkx*xr + dky*yr + dkz*zr
               sc(4) = qirx*xr + qiry*yr + qirz*zr
               sc(5) = qkrx*xr + qkry*yr + qkrz*zr
               sc(6) = dkx*qirx + dky*qiry + dkz*qirz
               sc(7) = dix*qkrx + diy*qkry + diz*qkrz
               sc(8) = qirx*qkrx + qiry*qkry + qirz*qkrz
               sc(9) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                    + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     construct auxiliary variables for multipole energy
c
               ge(1) = ci*ck
               ge(2) = ck*sc(2) - ci*sc(3) + sc(1)
               ge(3) = ci*sc(5) + ck*sc(4) - sc(2)*sc(3)
     &                    + 2.0d0*(sc(6)-sc(7)+sc(9))
               ge(4) = sc(2)*sc(5) - sc(3)*sc(4) - 4.0d0*sc(8)
               ge(5) = sc(4)*sc(5)
c
c     compute the full energy without any Ewald scaling
c
               efull = rr1*ge(1) + rr3*ge(2) + rr5*ge(3)
     &                    + rr7*ge(4) + rr9*ge(5)
               efull = efull * mscale(kk)
               if (ii .eq. kk)  efull = 0.5d0 * e
c
c     modify distances to account for Ewald and exclusions
c
               scalekk = 1.0d0 - mscale(kk)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
               rr11 = bn(5) - scalekk*rr11
c
c     compute the energy contribution for this interaction
c
               e = rr1*ge(1) + rr3*ge(2) + rr5*ge(3)
     &                + rr7*ge(4) + rr9*ge(5)
               if (ii .eq. kk)  e = 0.5d0 * e
               em = em + e
c
c     increment the total intramolecular energy
c
               if (molcule(ii) .eq. molcule(kk))
     &            eintra = eintra + efull
c
c     construct auxiliary variables for force and torque
c
               gf(1) = rr3*ge(1) + rr5*ge(2) + rr7*ge(3)
     &                    + rr9*ge(4) + rr11*ge(5)
               gf(2) = -ck*rr3 + sc(3)*rr5 - sc(5)*rr7
               gf(3) = ci*rr3 + sc(2)*rr5 + sc(4)*rr7
               gf(4) = 2.0d0 * rr5
               gf(5) = 2.0d0 * (-ck*rr5+sc(3)*rr7-sc(5)*rr9)
               gf(6) = 2.0d0 * (-ci*rr5-sc(2)*rr7-sc(4)*rr9)
               gf(7) = 4.0d0 * rr7
c
c     compute the force components for both sites
c
               frcx = gf(1)*xr + gf(2)*dix + gf(3)*dkx
     &                   + gf(4)*(diqkx-dkqix) + gf(5)*qirx
     &                   + gf(6)*qkrx + gf(7)*(qikrx+qkirx)
               frcy = gf(1)*yr + gf(2)*diy + gf(3)*dky
     &                   + gf(4)*(diqky-dkqiy) + gf(5)*qiry
     &                   + gf(6)*qkry + gf(7)*(qikry+qkiry)
               frcz = gf(1)*zr + gf(2)*diz + gf(3)*dkz
     &                   + gf(4)*(diqkz-dkqiz) + gf(5)*qirz
     &                   + gf(6)*qkrz + gf(7)*(qikrz+qkirz)
c
c     compute the torque components for both sites
c
               ttmi(1) = -rr3*dikx + gf(2)*dirx
     &                      + gf(4)*(dqiqkx+dkqixr)
     &                      - gf(5)*qirxr - gf(7)*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + gf(2)*diry
     &                      + gf(4)*(dqiqky+dkqiyr)
     &                      - gf(5)*qiryr - gf(7)*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + gf(2)*dirz
     &                      + gf(4)*(dqiqkz+dkqizr)
     &                      - gf(5)*qirzr - gf(7)*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + gf(3)*dkrx
     &                      - gf(4)*(dqiqkx+diqkxr)
     &                      - gf(6)*qkrxr - gf(7)*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + gf(3)*dkry
     &                      - gf(4)*(dqiqky+diqkyr)
     &                      - gf(6)*qkryr - gf(7)*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + gf(3)*dkrz
     &                      - gf(4)*(dqiqkz+diqkzr)
     &                      - gf(6)*qkrzr - gf(7)*(qkirzr-qrrz)
c
c     force and torque components scaled for self-interactions
c
               if (ii .eq. kk) then
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  do j = 1, 3
                     ttmi(j) = 0.5d0 * ttmi(j)
                     ttmk(j) = 0.5d0 * ttmk(j)
                  end do
               end if
c
c     increment force-based gradient and torque on first site
c
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
            end do
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
      end if
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
         vxz = zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vxy
         vir(3,1) = vir(3,1) + vxz
         vir(1,2) = vir(1,2) + vxy
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vyz
         vir(1,3) = vir(1,3) + vxz
         vir(2,3) = vir(2,3) + vyz
         vir(3,3) = vir(3,3) + vzz
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole1d  --  Ewald multipole derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1d" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using particle mesh Ewald
c     summation and a neighbor list
c
c
      subroutine empole1d
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use inter
      use math
      use mpole
      use virial
      implicit none
      integer i,j,ii
      real*8 e,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 trq(3),frcx(3)
      real*8 frcy(3),frcz(3)
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the real space part of the Ewald summation
c
      call emreal1d (eintra)
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip1
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         do i = 1, npole
            ii = ipole(i)
            dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
            dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
            dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         do i = 1, npole
            trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
            trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
            trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
            call torque (i,trq,frcx,frcy,frcz,dem)
         end do
c
c     boundary correction to virial due to overall cell dipole
c
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xq = 0.0d0
         yq = 0.0d0
         zq = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i)
            yd = yd + rpole(3,i)
            zd = zd + rpole(4,i)
            xq = xq + rpole(1,i)*x(ii)
            yq = yq + rpole(1,i)*y(ii)
            zq = zq + rpole(1,i)*z(ii)
         end do
         xv = xd * xq
         yv = yd * yq
         zv = zd * zq
         vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                      + xq*xq + yq*yq + zq*zq)
         vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
         vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
         vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
         vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
         vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
         vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
         vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
         vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
         vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + em - eintra
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1d  --  Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1d" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emreal1d (eintra)
      use sizes
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer iax,iay,iaz
      real*8 e,f,bfac,erfc
      real*8 eintra,eintrao
      real*8 emo,efull
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalekk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dikx,diky,dikz
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 qirx,qiry,qirz
      real*8 qkrx,qkry,qkrz
      real*8 qirxr,qiryr,qirzr
      real*8 qkrxr,qkryr,qkrzr
      real*8 qrrx,qrry,qrrz
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 qikrxr,qikryr,qikrzr
      real*8 qkirxr,qkiryr,qkirzr
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkxr,diqkyr,diqkzr
      real*8 dkqixr,dkqiyr,dkqizr
      real*8 dqiqkx,dqiqky,dqiqkz
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 sc(9),ge(5),gf(7)
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 bn(0:5)
      real*8 viro(3,3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: demo(:,:)
      real*8, allocatable :: temo(:,:)
      character*6 mode
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (demo(3,n))
      allocate (temo(3,n))
c
c     set array needed for scaling connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     initialize local variables for OpenMP calculation
c
      emo = 0.0d0
      eintrao = eintra
      do i = 1, n
         do j = 1, 3
            demo(j,i) = 0.0d0
         end do
         temo(1,i) = 0.0d0
         temo(2,i) = 0.0d0
         temo(3,i) = 0.0d0
      end do
      do i = 1, 3
         do j = 1, 3
            viro(j,i) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,rpole,n12,i12,n13,i13,n14,i14,n15,i15,
!$OMP& m2scale,m3scale,m4scale,m5scale,nelst,elst,use_bounds,f,off2,
!$OMP& aewald,molcule,xaxis,yaxis,zaxis,emo,eintrao,demo,temo,viro)
!$OMP& firstprivate(mscale)
!$OMP DO reduction(+:emo,eintrao,demo,temo,viro) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 5
                  bn(j) = f * bn(j)
               end do
c
c     construct several necessary additional variables
c
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
               qirx = qixx*xr + qixy*yr + qixz*zr
               qiry = qixy*xr + qiyy*yr + qiyz*zr
               qirz = qixz*xr + qiyz*yr + qizz*zr
               qkrx = qkxx*xr + qkxy*yr + qkxz*zr
               qkry = qkxy*xr + qkyy*yr + qkyz*zr
               qkrz = qkxz*xr + qkyz*yr + qkzz*zr
               qirxr = qirz*yr - qiry*zr
               qiryr = qirx*zr - qirz*xr
               qirzr = qiry*xr - qirx*yr
               qkrxr = qkrz*yr - qkry*zr
               qkryr = qkrx*zr - qkrz*xr
               qkrzr = qkry*xr - qkrx*yr
               qrrx = qkry*qirz - qkrz*qiry
               qrry = qkrz*qirx - qkrx*qirz
               qrrz = qkrx*qiry - qkry*qirx
               qikrx = qixx*qkrx + qixy*qkry + qixz*qkrz
               qikry = qixy*qkrx + qiyy*qkry + qiyz*qkrz
               qikrz = qixz*qkrx + qiyz*qkry + qizz*qkrz
               qkirx = qkxx*qirx + qkxy*qiry + qkxz*qirz
               qkiry = qkxy*qirx + qkyy*qiry + qkyz*qirz
               qkirz = qkxz*qirx + qkyz*qiry + qkzz*qirz
               qikrxr = qikrz*yr - qikry*zr
               qikryr = qikrx*zr - qikrz*xr
               qikrzr = qikry*xr - qikrx*yr
               qkirxr = qkirz*yr - qkiry*zr
               qkiryr = qkirx*zr - qkirz*xr
               qkirzr = qkiry*xr - qkirx*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqkxr = diqkz*yr - diqky*zr
               diqkyr = diqkx*zr - diqkz*xr
               diqkzr = diqky*xr - diqkx*yr
               dkqixr = dkqiz*yr - dkqiy*zr
               dkqiyr = dkqix*zr - dkqiz*xr
               dkqizr = dkqiy*xr - dkqix*yr
               dqiqkx = diy*qkrz - diz*qkry + dky*qirz - dkz*qiry
     &                     - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                             -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiqky = diz*qkrx - dix*qkrz + dkz*qirx - dkx*qirz
     &                     - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                             -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqiqkz = dix*qkry - diy*qkrx + dkx*qiry - dky*qirx
     &                     - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                             -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     calculate scalar products for multipole interactions
c
               sc(1) = dix*dkx + diy*dky + diz*dkz
               sc(2) = dix*xr + diy*yr + diz*zr
               sc(3) = dkx*xr + dky*yr + dkz*zr
               sc(4) = qirx*xr + qiry*yr + qirz*zr
               sc(5) = qkrx*xr + qkry*yr + qkrz*zr
               sc(6) = dkx*qirx + dky*qiry + dkz*qirz
               sc(7) = dix*qkrx + diy*qkry + diz*qkrz
               sc(8) = qirx*qkrx + qiry*qkry + qirz*qkrz
               sc(9) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                    + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     construct auxiliary variables for multipole energy
c
               ge(1) = ci*ck
               ge(2) = ck*sc(2) - ci*sc(3) + sc(1)
               ge(3) = ci*sc(5) + ck*sc(4) - sc(2)*sc(3)
     &                    + 2.0d0*(sc(6)-sc(7)+sc(9))
               ge(4) = sc(2)*sc(5) - sc(3)*sc(4) - 4.0d0*sc(8)
               ge(5) = sc(4)*sc(5)
c
c     compute the full energy without any Ewald scaling
c
               efull = rr1*ge(1) + rr3*ge(2) + rr5*ge(3)
     &                    + rr7*ge(4) + rr9*ge(5)
               efull = efull * mscale(kk)
c
c     modify distances to account for Ewald and exclusions
c
               scalekk = 1.0d0 - mscale(kk)
               rr1 = bn(0) - scalekk*rr1
               rr3 = bn(1) - scalekk*rr3
               rr5 = bn(2) - scalekk*rr5
               rr7 = bn(3) - scalekk*rr7
               rr9 = bn(4) - scalekk*rr9
               rr11 = bn(5) - scalekk*rr11
c
c     compute the energy contributions for this interaction
c
               e = rr1*ge(1) + rr3*ge(2) + rr5*ge(3)
     &                + rr7*ge(4) + rr9*ge(5)
               emo = emo + e
c
c     increment the total intramolecular energy
c
               if (molcule(ii) .eq. molcule(kk))
     &            eintrao = eintrao + efull
c
c     construct auxiliary variables for force and torque
c
               gf(1) = rr3*ge(1) + rr5*ge(2) + rr7*ge(3)
     &                    + rr9*ge(4) + rr11*ge(5)
               gf(2) = -ck*rr3 + sc(3)*rr5 - sc(5)*rr7
               gf(3) = ci*rr3 + sc(2)*rr5 + sc(4)*rr7
               gf(4) = 2.0d0 * rr5
               gf(5) = 2.0d0 * (-ck*rr5+sc(3)*rr7-sc(5)*rr9)
               gf(6) = 2.0d0 * (-ci*rr5-sc(2)*rr7-sc(4)*rr9)
               gf(7) = 4.0d0 * rr7
c
c     compute the force components for both sites
c
               frcx = gf(1)*xr + gf(2)*dix + gf(3)*dkx
     &                   + gf(4)*(diqkx-dkqix) + gf(5)*qirx
     &                   + gf(6)*qkrx + gf(7)*(qikrx+qkirx)
               frcy = gf(1)*yr + gf(2)*diy + gf(3)*dky
     &                   + gf(4)*(diqky-dkqiy) + gf(5)*qiry
     &                   + gf(6)*qkry + gf(7)*(qikry+qkiry)
               frcz = gf(1)*zr + gf(2)*diz + gf(3)*dkz
     &                   + gf(4)*(diqkz-dkqiz) + gf(5)*qirz
     &                   + gf(6)*qkrz + gf(7)*(qikrz+qkirz)
c
c     compute the torque components for both sites
c
               ttmi(1) = -rr3*dikx + gf(2)*dirx
     &                      + gf(4)*(dqiqkx+dkqixr)
     &                      - gf(5)*qirxr - gf(7)*(qikrxr+qrrx)
               ttmi(2) = -rr3*diky + gf(2)*diry
     &                      + gf(4)*(dqiqky+dkqiyr)
     &                      - gf(5)*qiryr - gf(7)*(qikryr+qrry)
               ttmi(3) = -rr3*dikz + gf(2)*dirz
     &                      + gf(4)*(dqiqkz+dkqizr)
     &                      - gf(5)*qirzr - gf(7)*(qikrzr+qrrz)
               ttmk(1) = rr3*dikx + gf(3)*dkrx
     &                      - gf(4)*(dqiqkx+diqkxr)
     &                      - gf(6)*qkrxr - gf(7)*(qkirxr-qrrx)
               ttmk(2) = rr3*diky + gf(3)*dkry
     &                      - gf(4)*(dqiqky+diqkyr)
     &                      - gf(6)*qkryr - gf(7)*(qkiryr-qrry)
               ttmk(3) = rr3*dikz + gf(3)*dkrz
     &                      - gf(4)*(dqiqkz+diqkzr)
     &                      - gf(6)*qkrzr - gf(7)*(qkirzr-qrrz)
c
c     increment force-based gradient and torque on first site
c
               demo(1,ii) = demo(1,ii) + frcx
               demo(2,ii) = demo(2,ii) + frcy
               demo(3,ii) = demo(3,ii) + frcz
               temo(1,i) = temo(1,i) + ttmi(1)
               temo(2,i) = temo(2,i) + ttmi(2)
               temo(3,i) = temo(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               demo(1,kk) = demo(1,kk) - frcx
               demo(2,kk) = demo(2,kk) - frcy
               demo(3,kk) = demo(3,kk) - frcz
               temo(1,k) = temo(1,k) + ttmk(1)
               temo(2,k) = temo(2,k) + ttmk(2)
               temo(3,k) = temo(3,k) + ttmk(3)
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -yr * frcx
               vxz = -zr * frcx
               vyy = -yr * frcy
               vyz = -zr * frcy
               vzz = -zr * frcz
               viro(1,1) = viro(1,1) + vxx
               viro(2,1) = viro(2,1) + vxy
               viro(3,1) = viro(3,1) + vxz
               viro(1,2) = viro(1,2) + vxy
               viro(2,2) = viro(2,2) + vyy
               viro(3,2) = viro(3,2) + vyz
               viro(1,3) = viro(1,3) + vxz
               viro(2,3) = viro(2,3) + vyz
               viro(3,3) = viro(3,3) + vzz
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO reduction(+:demo,viro) schedule(guided)
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,temo(1,i),fix,fiy,fiz,demo)
         ii = ipole(i)
         iaz = zaxis(i)
         iax = xaxis(i)
         iay = yaxis(i)
         if (iaz .eq. 0)  iaz = ii
         if (iax .eq. 0)  iax = ii
         if (iay .eq. 0)  iay = ii
         xiz = x(iaz) - x(ii)
         yiz = y(iaz) - y(ii)
         ziz = z(iaz) - z(ii)
         xix = x(iax) - x(ii)
         yix = y(iax) - y(ii)
         zix = z(iax) - z(ii)
         xiy = x(iay) - x(ii)
         yiy = y(iay) - y(ii)
         ziy = z(iay) - z(ii)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
         vxz = zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
         vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
         viro(1,1) = viro(1,1) + vxx
         viro(2,1) = viro(2,1) + vxy
         viro(3,1) = viro(3,1) + vxz
         viro(1,2) = viro(1,2) + vxy
         viro(2,2) = viro(2,2) + vyy
         viro(3,2) = viro(3,2) + vyz
         viro(1,3) = viro(1,3) + vxz
         viro(2,3) = viro(2,3) + vyz
         viro(3,3) = viro(3,3) + vzz
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add local to global variables for OpenMP calculation
c
      em = em + emo
      eintra = eintrao
      do i = 1, n
         do j = 1, 3
            dem(j,i) = dem(j,i) + demo(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + viro(j,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (demo)
      deallocate (temo)
      return
      end
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
c     mesh Ewald summation energy and gradient due to multipoles
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
      use sizes
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use mrecip
      use pme
      use virial
      implicit none
      integer i,j,k,ii
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
c
c     indices into the electrostatic field array
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
c     perform dynamic allocation of some global arrays
c
      if (allocated(cmp)) then
         if (size(cmp) .lt. 10*npole) then
            deallocate (cmp)
            deallocate (fmp)
            deallocate (cphi)
            deallocate (fphi)
         end if
      end if
      if (.not. allocated(cmp)) then
         allocate (cmp(10,npole))
         allocate (fmp(10,npole))
         allocate (cphi(10,npole))
         allocate (fphi(20,npole))
      end if
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vxy = 0.0d0
      vxz = 0.0d0
      vyy = 0.0d0
      vyz = 0.0d0
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
      end do
c
c     compute the arrays of B-spline coefficients
c
      call bspline_fill
      call table_fill
c
c     assign permanent multipoles to PME grid and perform
c     the 3-D FFT forward transformation
c
      call cmp_to_fmp (cmp,fmp)
      call grid_mpole (fmp)
      call fftfront
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
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
            struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
            eterm = 0.5d0 * electric * expterm * struc2
            vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
            vxx = vxx + h1*h1*vterm - eterm
            vxy = vxy + h1*h2*vterm
            vxz = vxz + h1*h3*vterm
            vyy = vyy + h2*h2*vterm - eterm
            vyz = vyz + h2*h3*vterm
            vzz = vzz + h3*h3*vterm - eterm
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     save the virial for use in polarization computation
c
      vmxx = vxx
      vmxy = vxy
      vmxz = vxz
      vmyy = vyy
      vmyz = vyz
      vmzz = vzz
c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
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
c     perform 3-D FFT backward transform and get potential
c
      call fftback
      call fphi_mpole (fphi)
      do i = 1, npole
         do j = 1, 20
            fphi(j,i) = electric * fphi(j,i)
         end do
      end do
      call fphi_to_cphi (fphi,cphi)
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
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         ii = ipole(i)
         dem(1,ii) = dem(1,ii) + h1
         dem(2,ii) = dem(2,ii) + h2
         dem(3,ii) = dem(3,ii) + h3
      end do
      e = 0.5d0 * e
      em = em + e
c
c     distribute torques into the permanent multipole gradient
c
      do i = 1, npole
         trq(1) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
     &               + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
     &               + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
     &               - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
         trq(2) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
     &               + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
     &               + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
     &               - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
         trq(3) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
     &               + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
     &               + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
     &               - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)
         call torque (i,trq,fix,fiy,fiz,dem)
      end do
c
c     permanent multipole contribution to the internal virial
c
      do i = 1, npole
         vxx = vxx - cmp(2,i)*cphi(2,i) - 2.0d0*cmp(5,i)*cphi(5,i)
     &            - cmp(8,i)*cphi(8,i) - cmp(9,i)*cphi(9,i)
         vxy = vxy - 0.5d0*(cmp(3,i)*cphi(2,i)+cmp(2,i)*cphi(3,i))
     &            - (cmp(5,i)+cmp(6,i))*cphi(8,i)
     &            - 0.5d0*cmp(8,i)*(cphi(5,i)+cphi(6,i))
     &            - 0.5d0*(cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
         vxz = vxz - 0.5d0*(cmp(4,i)*cphi(2,i)+cmp(2,i)*cphi(4,i))
     &            - (cmp(5,i)+cmp(7,i))*cphi(9,i)
     &            - 0.5d0*cmp(9,i)*(cphi(5,i)+cphi(7,i))
     &            - 0.5d0*(cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
         vyy = vyy - cmp(3,i)*cphi(3,i) - 2.0d0*cmp(6,i)*cphi(6,i)
     &            - cmp(8,i)*cphi(8,i) - cmp(10,i)*cphi(10,i)
         vyz = vyz - 0.5d0*(cmp(4,i)*cphi(3,i)+cmp(3,i)*cphi(4,i))
     &            - (cmp(6,i)+cmp(7,i))*cphi(10,i)
     &            - 0.5d0*cmp(10,i)*(cphi(6,i)+cphi(7,i))
     &            - 0.5d0*(cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
         vzz = vzz - cmp(4,i)*cphi(4,i) - 2.0d0*cmp(7,i)*cphi(7,i)
     &            - cmp(9,i)*cphi(9,i) - cmp(10,i)*cphi(10,i)
      end do
c
c     increment the internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vxy
      vir(3,1) = vir(3,1) + vxz
      vir(1,2) = vir(1,2) + vxy
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vyz
      vir(1,3) = vir(1,3) + vxz
      vir(2,3) = vir(2,3) + vyz
      vir(3,3) = vir(3,3) + vzz
      return
      end

c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole2  --  permanent atomic multipole Hessian  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole2" calculates second derivatives of the multipole energy
c     for a single atom at a time
c
c     it is incorrect to neglect interactions with atoms not directly
c     involved as the multipole site or its local axis definition; to
c     get better accuracy, "list" should include all atoms, an option
c     available by setting "biglist" to "true"
c
c     the "twosided" flag controls use of one-sided vs. two-sided
c     numerical derivatives; setting the flag to "true" gives a more
c     accurate Hessian at the expense of increased computation time
c
c
      subroutine empole2 (i)
      use sizes
      use atoms
      use deriv
      use hessn
      use limits
      use mpole
      use potent
      implicit none
      integer i,j,k
      integer nlist
      integer, allocatable :: list(:)
      real*8 eps,old
      real*8, allocatable :: d0(:,:)
      logical prior
      logical biglist
      logical twosided
c
c
c     set the default stepsize and accuracy control flags
c
      eps = 1.0d-5
      biglist = .false.
      twosided = .false.
      if (n .le. 300) then
         biglist = .true.
         twosided = .true.
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(npole))
      allocate (d0(3,n))
c
c     perform dynamic allocation of some global arrays
c
      prior = .false.
      if (allocated(dem)) then
         prior = .true.
         if (size(dem) .lt. 3*n)  deallocate (dem)
      end if
      if (.not. allocated(dem))  allocate (dem(3,n))
c
c     find the multipole definitions involving the current atom;
c     results in a faster but approximate Hessian calculation
c
      nlist = 0
      do k = 1, npole
         if (biglist .or. ipole(k).eq.i .or. zaxis(k).eq.i
     &          .or. xaxis(k).eq.i .or. yaxis(k).eq.i) then
            nlist = nlist + 1
            list(nlist) = k
         end if
      end do
c
c     get multipole first derivatives for the base structure
c
      if (.not. twosided) then
         if (use_mlist) then
            call empole2b (nlist,list)
         else
            call empole2a (nlist,list)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dem(j,k)
            end do
         end do
      end if
c
c     find numerical x-components via perturbed structures
c
      old = x(i)
      if (twosided) then
         x(i) = x(i) - 0.5d0*eps
         if (use_mlist) then
            call empole2b (nlist,list)
         else
            call empole2a (nlist,list)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dem(j,k)
            end do
         end do
      end if
      x(i) = x(i) + eps
      if (use_mlist) then
         call empole2b (nlist,list)
      else
         call empole2a (nlist,list)
      end if
      x(i) = old
      do k = 1, n
         do j = 1, 3
            hessx(j,k) = hessx(j,k) + (dem(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical y-components via perturbed structures
c
      old = y(i)
      if (twosided) then
         y(i) = y(i) - 0.5d0*eps
         if (use_mlist) then
            call empole2b (nlist,list)
         else
            call empole2a (nlist,list)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dem(j,k)
            end do
         end do
      end if
      y(i) = y(i) + eps
      if (use_mlist) then
         call empole2b (nlist,list)
      else
         call empole2a (nlist,list)
      end if
      y(i) = old
      do k = 1, n
         do j = 1, 3
            hessy(j,k) = hessy(j,k) + (dem(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical z-components via perturbed structures
c
      old = z(i)
      if (twosided) then
         z(i) = z(i) - 0.5d0*eps
         if (use_mlist) then
            call empole2b (nlist,list)
         else
            call empole2a (nlist,list)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dem(j,k)
            end do
         end do
      end if
      z(i) = z(i) + eps
      if (use_mlist) then
         call empole2b (nlist,list)
      else
         call empole2a (nlist,list)
      end if
      z(i) = old
      do k = 1, n
         do j = 1, 3
            hessz(j,k) = hessz(j,k) + (dem(j,k)-d0(j,k))/eps
         end do
      end do
c
c     perform deallocation of some global arrays
c
      if (.not. prior) then
         deallocate (dem)
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (d0)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole2a  --  multipole Hessian; numer, loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole2a" computes multipole first derivatives for a single
c     atom with respect to Cartesian coordinates; used to get finite
c     difference second derivatives
c
c
      subroutine empole2a (nlist,list)
      use sizes
      use atoms
      use bound
      use boxes
      use cell
      use chgpot
      use couple
      use deriv
      use group
      use limits
      use molcul
      use mplpot
      use mpole
      use potent
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,iii
      integer nlist,jcell
      integer ix,iy,iz
      integer kx,ky,kz
      integer list(*)
      real*8 f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
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
      real*8 sc(9),ge(5),gf(7)
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      logical proceed
      logical usei,usek
      character*6 mode
c
c
c     zero out the multipole first derivative components
c
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
         end do
      end do
      if (nlist .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
c
c     initialize connected atom scaling and torque arrays
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
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute components of the multipole interaction gradient
c
      do iii = 1, nlist
         i = list(iii)
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
               tem(1,ii) = tem(1,ii) + ttmi(1)
               tem(2,ii) = tem(2,ii) + ttmi(2)
               tem(3,ii) = tem(3,ii) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,kk) = tem(1,kk) + ttmk(1)
               tem(2,kk) = tem(2,kk) + ttmk(2)
               tem(3,kk) = tem(3,kk) + ttmk(3)
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
      do iii = 1, nlist
         i = list(iii)
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
               tem(1,ii) = tem(1,ii) + ttmi(1)
               tem(2,ii) = tem(2,ii) + ttmi(2)
               tem(3,ii) = tem(3,ii) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,kk) = tem(1,kk) + ttmk(1)
               tem(2,kk) = tem(2,kk) + ttmk(2)
               tem(3,kk) = tem(3,kk) + ttmk(3)
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
c     resolve site torques then increment multipole forces
c
      do i = 1, npole
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole2b  --  multipole Hessian; numer, list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole2b" calculates the multipole energy and derivatives
c     with respect to Cartesian coordinates using a neighbor list
c
c
      subroutine empole2b (nlist,list)
      use sizes
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use group
      use mplpot
      use mpole
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer iii,kkk
      integer nlist
      integer ix,iy,iz
      integer kx,ky,kz
      integer list(*)
      real*8 f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
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
      real*8 sc(9),ge(5),gf(7)
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      logical proceed,usei,usek
      character*6 mode
c
c
c     zero out the atomic multipole derivative components
c
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
c     initialize connected atom scaling and torque arrays
c
      if (npole .eq. 0)  return
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and scaling coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(nlist,list,npole,ipole,x,y,z,xaxis,yaxis,zaxis,rpole,use,
!$OMP& n12,i12,n13,i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,
!$OMP& nelst,elst,use_group,use_intra,use_bounds,off2,f,dem,tem)
!$OMP& firstprivate(mscale)
!$OMP DO reduction(+:dem,tem) schedule(guided)
c
c     compute the multipole interaction energy and gradient
c
      do iii = 1, nlist
         i = list(iii)
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
               dem(1,ii) = dem(1,ii) + frcx
               dem(2,ii) = dem(2,ii) + frcy
               dem(3,ii) = dem(3,ii) + frcz
               tem(1,ii) = tem(1,ii) + ttmi(1)
               tem(2,ii) = tem(2,ii) + ttmi(2)
               tem(3,ii) = tem(3,ii) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,kk) = dem(1,kk) - frcx
               dem(2,kk) = dem(2,kk) - frcy
               dem(3,kk) = dem(3,kk) - frcz
               tem(1,kk) = tem(1,kk) + ttmk(1)
               tem(2,kk) = tem(2,kk) + ttmk(2)
               tem(3,kk) = tem(3,kk) + ttmk(3)
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
!$OMP DO reduction(+:dem) schedule(guided)
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      return
      end

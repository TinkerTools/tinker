c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empole2  --  atomic multipole Hessian elements  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empole2" calculates second derivatives of the multipole energy
c     for a single atom at a time
c
c
      subroutine empole2 (i)
      use atoms
      use deriv
      use hessn
      use mpole
      use potent
      implicit none
      integer i,j,k
      integer nlist
      integer, allocatable :: list(:)
      real*8 eps,old
      real*8, allocatable :: d0(:,:)
      logical prior
      logical twosided
c
c
c     set the default stepsize and accuracy control flags
c
      eps = 1.0d-5
      twosided = .false.
      if (n .le. 300)  twosided = .true.
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
c     find the multipole definitions involving the current atom
c
      nlist = 0
      do k = 1, npole
         if (ipole(k).eq.i .or. zaxis(k).eq.i
     &          .or. xaxis(k).eq.i .or. abs(yaxis(k)).eq.i) then
            nlist = nlist + 1
            list(nlist) = k
         end if
      end do
c
c     get multipole first derivatives for the base structure
c
      if (.not. twosided) then
         call empole2a (nlist,list)
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
         call empole2a (nlist,list)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dem(j,k)
            end do
         end do
      end if
      x(i) = x(i) + eps
      call empole2a (nlist,list)
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
         call empole2a (nlist,list)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dem(j,k)
            end do
         end do
      end if
      y(i) = y(i) + eps
      call empole2a (nlist,list)
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
         call empole2a (nlist,list)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dem(j,k)
            end do
         end do
      end if
      z(i) = z(i) + eps
      call empole2a (nlist,list)
      z(i) = old
      do k = 1, n
         do j = 1, 3
            hessz(j,k) = hessz(j,k) + (dem(j,k)-d0(j,k))/eps
         end do
      end do
c
c     perform deallocation of some global arrays
c
      if (.not. prior)  deallocate (dem)
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (d0)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine empole2a  --  numerical multipole Hessian  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "empole2a" computes multipole first derivatives for a single
c     atom; used to get finite difference second derivatives
c
c
      subroutine empole2a (nlist,list)
      use atoms
      use bound
      use boxes
      use cell
      use chgpen
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
      real*8 de,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 rr1i,rr3i,rr5i,rr7i
      real*8 rr1k,rr3k,rr5k,rr7k
      real*8 rr1ik,rr3ik,rr5ik
      real*8 rr7ik,rr9ik,rr11ik
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr,dik,qik
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 diqk,dkqi,qiqk
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 dikx,diky,dikz
      real*8 qirx,qiry,qirz
      real*8 qkrx,qkry,qkrz
      real*8 qikx,qiky,qikz
      real*8 qixk,qiyk,qizk
      real*8 qkxi,qkyi,qkzi
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkrx,diqkry,diqkrz
      real*8 dkqirx,dkqiry,dkqirz
      real*8 dqikx,dqiky,dqikz
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 term1i,term2i,term3i
      real*8 term1k,term2k,term3k
      real*8 term1ik,term2ik
      real*8 term3ik,term4ik
      real*8 term5ik,term6ik
      real*8 frcx,frcy,frcz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(11)
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
         ii = list(iii)
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = abs(yaxis(ii))
         xi = x(i)
         yi = y(i)
         zi = z(i)
         ci = rpole(1,ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
         usei = (use(i) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set exclusion coefficients for connected atoms
c
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
c
c     evaluate all sites within the cutoff distance
c
         do kk = 1, npole
            if (kk .eq. ii)  goto 10
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            ky = abs(yaxis(kk))
            usek = (use(k) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 10
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,kk)
               dkx = rpole(2,kk)
               dky = rpole(3,kk)
               dkz = rpole(4,kk)
               qkxx = rpole(5,kk)
               qkxy = rpole(6,kk)
               qkxz = rpole(7,kk)
               qkyy = rpole(9,kk)
               qkyz = rpole(10,kk)
               qkzz = rpole(13,kk)
c
c     intermediates involving moments and separation distance
c
               dir = dix*xr + diy*yr + diz*zr
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qir = qix*xr + qiy*yr + qiz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
               qkr = qkx*xr + qky*yr + qkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qik = qix*qkx + qiy*qky + qiz*qkz
               diqk = dix*qkx + diy*qky + diz*qkz
               dkqi = dkx*qix + dky*qiy + dkz*qiz
               qiqk = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                   + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     additional intermediates involving moments and distance
c
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               qirx = qiz*yr - qiy*zr
               qiry = qix*zr - qiz*xr
               qirz = qiy*xr - qix*yr
               qkrx = qkz*yr - qky*zr
               qkry = qkx*zr - qkz*xr
               qkrz = qky*xr - qkx*yr
               qikx = qky*qiz - qkz*qiy
               qiky = qkz*qix - qkx*qiz
               qikz = qkx*qiy - qky*qix
               qixk = qixx*qkx + qixy*qky + qixz*qkz
               qiyk = qixy*qkx + qiyy*qky + qiyz*qkz
               qizk = qixz*qkx + qiyz*qky + qizz*qkz
               qkxi = qkxx*qix + qkxy*qiy + qkxz*qiz
               qkyi = qkxy*qix + qkyy*qiy + qkyz*qiz
               qkzi = qkxz*qix + qkyz*qiy + qkzz*qiz
               qikrx = qizk*yr - qiyk*zr
               qikry = qixk*zr - qizk*xr
               qikrz = qiyk*xr - qixk*yr
               qkirx = qkzi*yr - qkyi*zr
               qkiry = qkxi*zr - qkzi*xr
               qkirz = qkyi*xr - qkxi*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqkrx = diqkz*yr - diqky*zr
               diqkry = diqkx*zr - diqkz*xr
               diqkrz = diqky*xr - diqkx*yr
               dkqirx = dkqiz*yr - dkqiy*zr
               dkqiry = dkqix*zr - dkqiz*xr
               dkqirz = dkqiy*xr - dkqix*yr
               dqikx = diy*qkz - diz*qky + dky*qiz - dkz*qiy
     &                 - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                         -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiky = diz*qkx - dix*qkz + dkz*qix - dkx*qiz
     &                 - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                         -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqikz = dix*qky - diy*qkx + dkx*qiy - dky*qix
     &                 - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                         -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f * mscale(k) / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     find damped multipole intermediates for force and torque
c
               if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  term1 = corei*corek
                  term1i = corek*vali
                  term2i = corek*dir
                  term3i = corek*qir
                  term1k = corei*valk
                  term2k = -corei*dkr
                  term3k = corei*qkr
                  term1ik = vali*valk
                  term2ik = valk*dir - vali*dkr + dik
                  term3ik = vali*qkr + valk*qir - dir*dkr
     &                         + 2.0d0*(dkqi-diqk+qiqk)
                  term4ik = dir*qkr - dkr*qir - 4.0d0*qik
                  term5ik = qir*qkr
                  call damppole (r,11,alphai,alphak,
     &                            dmpi,dmpk,dmpik)
                  rr3i = dmpi(3)*rr3
                  rr5i = dmpi(5)*rr5
                  rr7i = dmpi(7)*rr7
                  rr3k = dmpk(3)*rr3
                  rr5k = dmpk(5)*rr5
                  rr7k = dmpk(7)*rr7
                  rr3ik = dmpik(3)*rr3
                  rr5ik = dmpik(5)*rr5
                  rr7ik = dmpik(7)*rr7
                  rr9ik = dmpik(9)*rr9
                  rr11ik = dmpik(11)*rr11
                  de = term1*rr3 + term4ik*rr9ik + term5ik*rr11ik
     &                    + term1i*rr3i + term1k*rr3k + term1ik*rr3ik
     &                    + term2i*rr5i + term2k*rr5k + term2ik*rr5ik
     &                    + term3i*rr7i + term3k*rr7k + term3ik*rr7ik
                  term1 = -corek*rr3i - valk*rr3ik
     &                       + dkr*rr5ik - qkr*rr7ik
                  term2 = corei*rr3k + vali*rr3ik
     &                       + dir*rr5ik + qir*rr7ik
                  term3 = 2.0d0 * rr5ik
                  term4 = -2.0d0 * (corek*rr5i+valk*rr5ik
     &                                -dkr*rr7ik+qkr*rr9ik)
                  term5 = -2.0d0 * (corei*rr5k+vali*rr5ik
     &                                +dir*rr7ik+qir*rr9ik)
                  term6 = 4.0d0 * rr7ik
                  rr3 = rr3ik
c
c     find standard multipole intermediates for force and torque
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  de = term1*rr3 + term2*rr5 + term3*rr7
     &                    + term4*rr9 + term5*rr11
                  term1 = -ck*rr3 + dkr*rr5 - qkr*rr7
                  term2 = ci*rr3 + dir*rr5 + qir*rr7
                  term3 = 2.0d0 * rr5
                  term4 = 2.0d0 * (-ck*rr5+dkr*rr7-qkr*rr9)
                  term5 = 2.0d0 * (-ci*rr5-dir*rr7-qir*rr9)
                  term6 = 4.0d0 * rr7
               end if
c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qix
     &                   + term5*qkx + term6*(qixk+qkxi)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qiy
     &                   + term5*qky + term6*(qiyk+qkyi)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qiz
     &                   + term5*qkz + term6*(qizk+qkzi)
c
c     compute the torque components for this interaction
c
               ttmi(1) = -rr3*dikx + term1*dirx
     &                      + term3*(dqikx+dkqirx)
     &                      - term4*qirx - term6*(qikrx+qikx)
               ttmi(2) = -rr3*diky + term1*diry
     &                      + term3*(dqiky+dkqiry)
     &                      - term4*qiry - term6*(qikry+qiky)
               ttmi(3) = -rr3*dikz + term1*dirz
     &                      + term3*(dqikz+dkqirz)
     &                      - term4*qirz - term6*(qikrz+qikz)
               ttmk(1) = rr3*dikx + term2*dkrx
     &                      - term3*(dqikx+diqkrx)
     &                      - term5*qkrx - term6*(qkirx-qikx)
               ttmk(2) = rr3*diky + term2*dkry
     &                      - term3*(dqiky+diqkry)
     &                      - term5*qkry - term6*(qkiry-qiky)
               ttmk(3) = rr3*dikz + term2*dkrz
     &                      - term3*(dqikz+diqkrz)
     &                      - term5*qkrz - term6*(qkirz-qikz)
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
               dem(1,i) = dem(1,i) + frcx
               dem(2,i) = dem(2,i) + frcy
               dem(3,i) = dem(3,i) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,k) = dem(1,k) - frcx
               dem(2,k) = dem(2,k) - frcy
               dem(3,k) = dem(3,k) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
            end if
   10       continue
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = 1.0d0
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
         ii = list(iii)
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = abs(yaxis(ii))
         xi = x(i)
         yi = y(i)
         zi = z(i)
         ci = rpole(1,ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
         usei = (use(i) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set exclusion coefficients for connected atoms
c
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
c
c     evaluate all sites within the cutoff distance
c
         do kk = 1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            ky = abs(yaxis(kk))
            usek = (use(k) .or. use(kz) .or. use(kx) .or. use(ky))
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 20
            do jcell = 2, ncell
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               mscale(k) = 1.0d0
            end if
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,kk)
               dkx = rpole(2,kk)
               dky = rpole(3,kk)
               dkz = rpole(4,kk)
               qkxx = rpole(5,kk)
               qkxy = rpole(6,kk)
               qkxz = rpole(7,kk)
               qkyy = rpole(9,kk)
               qkyz = rpole(10,kk)
               qkzz = rpole(13,kk)
c
c     intermediates involving moments and separation distance
c
               dir = dix*xr + diy*yr + diz*zr
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qir = qix*xr + qiy*yr + qiz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
               qkr = qkx*xr + qky*yr + qkz*zr
               dik = dix*dkx + diy*dky + diz*dkz
               qik = qix*qkx + qiy*qky + qiz*qkz
               diqk = dix*qkx + diy*qky + diz*qkz
               dkqi = dkx*qix + dky*qiy + dkz*qiz
               qiqk = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                   + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     additional intermediates involving moments and distance
c
               dirx = diy*zr - diz*yr
               diry = diz*xr - dix*zr
               dirz = dix*yr - diy*xr
               dkrx = dky*zr - dkz*yr
               dkry = dkz*xr - dkx*zr
               dkrz = dkx*yr - dky*xr
               dikx = diy*dkz - diz*dky
               diky = diz*dkx - dix*dkz
               dikz = dix*dky - diy*dkx
               qirx = qiz*yr - qiy*zr
               qiry = qix*zr - qiz*xr
               qirz = qiy*xr - qix*yr
               qkrx = qkz*yr - qky*zr
               qkry = qkx*zr - qkz*xr
               qkrz = qky*xr - qkx*yr
               qikx = qky*qiz - qkz*qiy
               qiky = qkz*qix - qkx*qiz
               qikz = qkx*qiy - qky*qix
               qixk = qixx*qkx + qixy*qky + qixz*qkz
               qiyk = qixy*qkx + qiyy*qky + qiyz*qkz
               qizk = qixz*qkx + qiyz*qky + qizz*qkz
               qkxi = qkxx*qix + qkxy*qiy + qkxz*qiz
               qkyi = qkxy*qix + qkyy*qiy + qkyz*qiz
               qkzi = qkxz*qix + qkyz*qiy + qkzz*qiz
               qikrx = qizk*yr - qiyk*zr
               qikry = qixk*zr - qizk*xr
               qikrz = qiyk*xr - qixk*yr
               qkirx = qkzi*yr - qkyi*zr
               qkiry = qkxi*zr - qkzi*xr
               qkirz = qkyi*xr - qkxi*yr
               diqkx = dix*qkxx + diy*qkxy + diz*qkxz
               diqky = dix*qkxy + diy*qkyy + diz*qkyz
               diqkz = dix*qkxz + diy*qkyz + diz*qkzz
               dkqix = dkx*qixx + dky*qixy + dkz*qixz
               dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
               dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
               diqkrx = diqkz*yr - diqky*zr
               diqkry = diqkx*zr - diqkz*xr
               diqkrz = diqky*xr - diqkx*yr
               dkqirx = dkqiz*yr - dkqiy*zr
               dkqiry = dkqix*zr - dkqiz*xr
               dkqirz = dkqiy*xr - dkqix*yr
               dqikx = diy*qkz - diz*qky + dky*qiz - dkz*qiy
     &                 - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                         -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
               dqiky = diz*qkx - dix*qkz + dkz*qix - dkx*qiz
     &                 - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                         -qixx*qkxz-qixy*qkyz-qixz*qkzz)
               dqikz = dix*qky - diy*qkx + dkx*qiy - dky*qix
     &                 - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                         -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f * mscale(k) / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
c
c     find damped multipole intermediates for force and torque
c
               if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  term1 = corei*corek
                  term1i = corek*vali
                  term2i = corek*dir
                  term3i = corek*qir
                  term1k = corei*valk
                  term2k = -corei*dkr
                  term3k = corei*qkr
                  term1ik = vali*valk
                  term2ik = valk*dir - vali*dkr + dik
                  term3ik = vali*qkr + valk*qir - dir*dkr
     &                         + 2.0d0*(dkqi-diqk+qiqk)
                  term4ik = dir*qkr - dkr*qir - 4.0d0*qik
                  term5ik = qir*qkr
                  call damppole (r,11,alphai,alphak,
     &                            dmpi,dmpk,dmpik)
                  rr3i = dmpi(3)*rr3
                  rr5i = dmpi(5)*rr5
                  rr7i = dmpi(7)*rr7
                  rr3k = dmpk(3)*rr3
                  rr5k = dmpk(5)*rr5
                  rr7k = dmpk(7)*rr7
                  rr3ik = dmpik(3)*rr3
                  rr5ik = dmpik(5)*rr5
                  rr7ik = dmpik(7)*rr7
                  rr9ik = dmpik(9)*rr9
                  rr11ik = dmpik(11)*rr11
                  de = term1*rr3 + term4ik*rr9ik + term5ik*rr11ik
     &                    + term1i*rr3i + term1k*rr3k + term1ik*rr3ik
     &                    + term2i*rr5i + term2k*rr5k + term2ik*rr5ik
     &                    + term3i*rr7i + term3k*rr7k + term3ik*rr7ik
                  term1 = -corek*rr3i - valk*rr3ik
     &                       + dkr*rr5ik - qkr*rr7ik
                  term2 = corei*rr3k + vali*rr3ik
     &                       + dir*rr5ik + qir*rr7ik
                  term3 = 2.0d0 * rr5ik
                  term4 = -2.0d0 * (corek*rr5i+valk*rr5ik
     &                                -dkr*rr7ik+qkr*rr9ik)
                  term5 = -2.0d0 * (corei*rr5k+vali*rr5ik
     &                                +dir*rr7ik+qir*rr9ik)
                  term6 = 4.0d0 * rr7ik
                  rr3 = rr3ik
c
c     find standard multipole intermediates for force and torque
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  de = term1*rr3 + term2*rr5 + term3*rr7
     &                    + term4*rr9 + term5*rr11
                  term1 = -ck*rr3 + dkr*rr5 - qkr*rr7
                  term2 = ci*rr3 + dir*rr5 + qir*rr7
                  term3 = 2.0d0 * rr5
                  term4 = 2.0d0 * (-ck*rr5+dkr*rr7-qkr*rr9)
                  term5 = 2.0d0 * (-ci*rr5-dir*rr7-qir*rr9)
                  term6 = 4.0d0 * rr7
               end if
c
c     compute the force components for this interaction
c
               frcx = de*xr + term1*dix + term2*dkx
     &                   + term3*(diqkx-dkqix) + term4*qix
     &                   + term5*qkx + term6*(qixk+qkxi)
               frcy = de*yr + term1*diy + term2*dky
     &                   + term3*(diqky-dkqiy) + term4*qiy
     &                   + term5*qky + term6*(qiyk+qkyi)
               frcz = de*zr + term1*diz + term2*dkz
     &                   + term3*(diqkz-dkqiz) + term4*qiz
     &                   + term5*qkz + term6*(qizk+qkzi)
c
c     compute the torque components for this interaction
c
               ttmi(1) = -rr3*dikx + term1*dirx
     &                      + term3*(dqikx+dkqirx)
     &                      - term4*qirx - term6*(qikrx+qikx)
               ttmi(2) = -rr3*diky + term1*diry
     &                      + term3*(dqiky+dkqiry)
     &                      - term4*qiry - term6*(qikry+qiky)
               ttmi(3) = -rr3*dikz + term1*dirz
     &                      + term3*(dqikz+dkqirz)
     &                      - term4*qirz - term6*(qikrz+qikz)
               ttmk(1) = rr3*dikx + term2*dkrx
     &                      - term3*(dqikx+diqkrx)
     &                      - term5*qkrx - term6*(qkirx-qikx)
               ttmk(2) = rr3*diky + term2*dkry
     &                      - term3*(dqiky+diqkry)
     &                      - term5*qkry - term6*(qkiry-qiky)
               ttmk(3) = rr3*dikz + term2*dkrz
     &                      - term3*(dqikz+diqkrz)
     &                      - term5*qkrz - term6*(qkirz-qikz)
c
c     force and torque scaled for self-interactions and groups
c
               if (i .eq. k) then
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
               dem(1,i) = dem(1,i) + frcx
               dem(2,i) = dem(2,i) + frcy
               dem(3,i) = dem(3,i) + frcz
               tem(1,i) = tem(1,i) + ttmi(1)
               tem(2,i) = tem(2,i) + ttmi(2)
               tem(3,i) = tem(3,i) + ttmi(3)
c
c     increment force-based gradient and torque on second site
c
               dem(1,k) = dem(1,k) - frcx
               dem(2,k) = dem(2,k) - frcy
               dem(3,k) = dem(3,k) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
            end if
            end do
   20       continue
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = 1.0d0
         end do
      end do
      end if
c
c     resolve site torques then increment multipole forces
c
      do ii = 1, npole
         i = ipole(ii)
         call torque (ii,tem(1,i),fix,fiy,fiz,dem)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end

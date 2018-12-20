c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine empole1  --  multipole energy & derivatives  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "empole1" calculates the atomic multipole energy and first
c     derivatives with respect to Cartesian coordinates
c
c
      subroutine empole1
      use limits
      implicit none
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
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use deriv
      use energi
      use group
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
      real*8 e,de,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
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
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(11)
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
c     compute the multipole interaction energy and gradient
c
      do ii = 1, npole-1
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = yaxis(ii)
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
         do kk = i+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            ky = yaxis(kk)
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
c     find damped multipole intermediates and energy value
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
                  rr1i = dmpi(1)*rr1
                  rr3i = dmpi(3)*rr3
                  rr5i = dmpi(5)*rr5
                  rr7i = dmpi(7)*rr7
                  rr1k = dmpk(1)*rr1
                  rr3k = dmpk(3)*rr3
                  rr5k = dmpk(5)*rr5
                  rr7k = dmpk(7)*rr7
                  rr1ik = dmpik(1)*rr1
                  rr3ik = dmpik(3)*rr3
                  rr5ik = dmpik(5)*rr5
                  rr7ik = dmpik(7)*rr7
                  rr9ik = dmpik(9)*rr9
                  rr11ik = dmpik(11)*rr11
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find damped multipole intermediates for force and torque
c
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
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
c
c     find standard multipole intermediates for force and torque
c
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
c     energy, force and torque scaled by group membership
c
               if (use_group) then
                  e = fgrp * e
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
                  do j = 1, 3
                     ttmi(j) = fgrp * ttmi(j)
                     ttmk(j) = fgrp * ttmk(j)
                  end do
               end if
c
c     increment the overall atomic multipole energy component
c
               em = em + e
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
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -0.5d0 * (yr*frcx+xr*frcy)
               vxz = -0.5d0 * (zr*frcx+xr*frcz)
               vyy = -yr * frcy
               vyz = -0.5d0 * (zr*frcy+yr*frcz)
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
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = yaxis(ii)
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
         do kk = i, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            ky = yaxis(kk)
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
c     find damped multipole intermediates and energy value
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
                  rr1i = dmpi(1)*rr1
                  rr3i = dmpi(3)*rr3
                  rr5i = dmpi(5)*rr5
                  rr7i = dmpi(7)*rr7
                  rr1k = dmpk(1)*rr1
                  rr3k = dmpk(3)*rr3
                  rr5k = dmpk(5)*rr5
                  rr7k = dmpk(7)*rr7
                  rr1ik = dmpik(1)*rr1
                  rr3ik = dmpik(3)*rr3
                  rr5ik = dmpik(5)*rr5
                  rr7ik = dmpik(7)*rr7
                  rr9ik = dmpik(9)*rr9
                  rr11ik = dmpik(11)*rr11
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find damped multipole intermediates for force and torque
c
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
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
c
c     find standard multipole intermediates for force and torque
c
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
c     energy, force and torque scaled by group membership
c
               if (i .eq. k) then
                  e = 0.5d0 * e
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  do j = 1, 3
                     ttmi(j) = 0.5d0 * ttmi(j)
                     ttmk(j) = 0.5d0 * ttmk(j)
                  end do
               end if
               if (use_group) then
                  e = fgrp * e
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
                  do j = 1, 3
                     ttmi(j) = fgrp * ttmi(j)
                     ttmk(j) = fgrp * ttmk(j)
                  end do
               end if
c
c     increment the overall atomic multipole energy component
c
               em = em + e
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
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -0.5d0 * (yr*frcx+xr*frcy)
               vxz = -0.5d0 * (zr*frcx+xr*frcz)
               vyy = -yr * frcy
               vyz = -0.5d0 * (zr*frcy+yr*frcz)
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
c     resolve site torques then increment forces and virial
c
      do ii = 1, npole
         i = ipole(ii)
         call torque (ii,tem(1,i),fix,fiy,fiz,dem)
         iaz = zaxis(ii)
         iax = xaxis(ii)
         iay = yaxis(ii)
         if (iaz .eq. 0)  iaz = i
         if (iax .eq. 0)  iax = i
         if (iay .eq. 0)  iay = i
         xiz = x(iaz) - x(i)
         yiz = y(iaz) - y(i)
         ziz = z(iaz) - z(i)
         xix = x(iax) - x(i)
         yix = y(iax) - y(i)
         zix = z(iax) - z(i)
         xiy = x(iay) - x(i)
         yiy = y(iay) - y(i)
         ziy = z(iay) - z(i)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
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
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use deriv
      use energi
      use group
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
      real*8 e,de,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
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
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(11)
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
c     initialize connected atom scaling and torque arrays
c
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
!$OMP& shared(npole,ipole,x,y,z,xaxis,yaxis,zaxis,rpole,pcore,
!$OMP& pval,palpha,use,n12,i12,n13,i13,n14,i14,n15,i15,m2scale,
!$OMP& m3scale,m4scale,m5scale,nelst,elst,use_chgpen,use_group,
!$OMP& use_intra,use_bounds,off2,f)
!$OMP& firstprivate(mscale) shared (em,dem,tem,vir)
!$OMP DO reduction(+:em,dem,tem,vir) schedule(guided)
c
c     compute the multipole interaction energy and gradient
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = yaxis(ii)
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
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            ky = yaxis(kk)
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
c     find damped multipole intermediates and energy value
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
                  rr1i = dmpi(1)*rr1
                  rr3i = dmpi(3)*rr3
                  rr5i = dmpi(5)*rr5
                  rr7i = dmpi(7)*rr7
                  rr1k = dmpk(1)*rr1
                  rr3k = dmpk(3)*rr3
                  rr5k = dmpk(5)*rr5
                  rr7k = dmpk(7)*rr7
                  rr1ik = dmpik(1)*rr1
                  rr3ik = dmpik(3)*rr3
                  rr5ik = dmpik(5)*rr5
                  rr7ik = dmpik(7)*rr7
                  rr9ik = dmpik(9)*rr9
                  rr11ik = dmpik(11)*rr11
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find damped multipole intermediates for force and torque
c
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
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
c
c     find standard multipole intermediates for force and torque
c
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
c     energy, force and torque scaled by group membership
c
               if (use_group) then
                  e = fgrp * e
                  frcx = fgrp * frcx
                  frcy = fgrp * frcy
                  frcz = fgrp * frcz
                  do j = 1, 3
                     ttmi(j) = fgrp * ttmi(j)
                     ttmk(j) = fgrp * ttmk(j)
                  end do
               end if
c
c     increment the overall atomic multipole energy component
c
               em = em + e
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
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -0.5d0 * (yr*frcx+xr*frcy)
               vxz = -0.5d0 * (zr*frcx+xr*frcz)
               vyy = -yr * frcy
               vyz = -0.5d0 * (zr*frcy+yr*frcz)
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO reduction(+:dem,vir) schedule(guided)
c
c     resolve site torques then increment forces and virial
c
      do ii = 1, npole
         i = ipole(ii)
         call torque (ii,tem(1,i),fix,fiy,fiz,dem)
         iaz = zaxis(ii)
         iax = xaxis(ii)
         iay = yaxis(ii)
         if (iaz .eq. 0)  iaz = i
         if (iax .eq. 0)  iax = i
         if (iay .eq. 0)  iay = i
         xiz = x(iaz) - x(i)
         yiz = y(iaz) - y(i)
         ziz = z(iaz) - z(i)
         xix = x(iax) - x(i)
         yix = y(iax) - y(i)
         zix = z(iax) - z(i)
         xiy = x(iay) - x(i)
         yiy = y(iay) - y(i)
         ziy = z(iay) - z(i)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
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
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use pme
      use virial
      implicit none
      integer i,j,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 tem(3),frcx(3)
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
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bseorder
      aewald = aeewald
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
      call emreal1c
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip1
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do ii = 1, npole
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
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
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
         do ii = 1, npole
            i = ipole(ii)
            xd = xd + rpole(2,ii) + rpole(1,ii)*x(i)
            yd = yd + rpole(3,ii) + rpole(1,ii)*y(i)
            zd = zd + rpole(4,ii) + rpole(1,ii)*z(i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         do ii = 1, npole
            i = ipole(ii)
            dem(1,i) = dem(1,i) + 2.0d0*term*rpole(1,ii)*xd
            dem(2,i) = dem(2,i) + 2.0d0*term*rpole(1,ii)*yd
            dem(3,i) = dem(3,i) + 2.0d0*term*rpole(1,ii)*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         do ii = 1, npole
            tem(1) = rpole(3,ii)*zdfield - rpole(4,ii)*ydfield
            tem(2) = rpole(4,ii)*xdfield - rpole(2,ii)*zdfield
            tem(3) = rpole(2,ii)*ydfield - rpole(3,ii)*xdfield
            call torque (ii,tem,frcx,frcy,frcz,dem)
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
         do ii = 1, npole
            i = ipole(ii)
            xd = xd + rpole(2,ii)
            yd = yd + rpole(3,ii)
            zd = zd + rpole(4,ii)
            xq = xq + rpole(1,ii)*x(i)
            yq = yq + rpole(1,ii)*y(i)
            zq = zq + rpole(1,ii)*z(i)
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
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1c  --  Ewald real mpole derivs via loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a double loop
c
c
      subroutine emreal1c
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use math
      use mplpot
      use mpole
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer iax,iay,iaz
      real*8 e,de,f
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
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
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(11)
      real*8 bn(0:5)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      character*6 mode
      external erfc
c
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
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npole-1
         i = ipole(ii)
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
         do kk = ii+1, npole
            k = ipole(kk)
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
c     find damped multipole intermediates and energy value
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
                  scalek = mscale(k)
                  rr1i = bn(0) - (1.0d0-scalek*dmpi(1))*rr1
                  rr3i = bn(1) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = bn(2) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = bn(3) - (1.0d0-scalek*dmpi(7))*rr7
                  rr1k = bn(0) - (1.0d0-scalek*dmpk(1))*rr1
                  rr3k = bn(1) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = bn(2) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = bn(3) - (1.0d0-scalek*dmpk(7))*rr7
                  rr1ik = bn(0) - (1.0d0-scalek*dmpik(1))*rr1
                  rr3ik = bn(1) - (1.0d0-scalek*dmpik(3))*rr3
                  rr5ik = bn(2) - (1.0d0-scalek*dmpik(5))*rr5
                  rr7ik = bn(3) - (1.0d0-scalek*dmpik(7))*rr7
                  rr9ik = bn(4) - (1.0d0-scalek*dmpik(9))*rr9
                  rr11ik = bn(5) - (1.0d0-scalek*dmpik(11))*rr11
                  rr1 = bn(0) - (1.0d0-scalek)*rr1
                  rr3 = bn(1) - (1.0d0-scalek)*rr3
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find damped multipole intermediates for force and torque
c
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
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  scalek = 1.0d0 - mscale(k)
                  rr1 = bn(0) - scalek*rr1
                  rr3 = bn(1) - scalek*rr3
                  rr5 = bn(2) - scalek*rr5
                  rr7 = bn(3) - scalek*rr7
                  rr9 = bn(4) - scalek*rr9
                  rr11 = bn(5) - scalek*rr11
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
c
c     find standard multipole intermediates for force and torque
c
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
c     increment the overall atomic multipole energy component
c
               em = em + e
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
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -0.5d0 * (yr*frcx+xr*frcy)
               vxz = -0.5d0 * (zr*frcx+xr*frcz)
               vyy = -yr * frcy
               vyz = -0.5d0 * (zr*frcy+yr*frcz)
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
      do ii = 1, npole
         i = ipole(ii)
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
         do kk = i, npole
            k = ipole(kk)
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
c     find damped multipole intermediates and energy value
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
                  scalek = mscale(k)
                  rr1i = bn(0) - (1.0d0-scalek*dmpi(1))*rr1
                  rr3i = bn(1) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = bn(2) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = bn(3) - (1.0d0-scalek*dmpi(7))*rr7
                  rr1k = bn(0) - (1.0d0-scalek*dmpk(1))*rr1
                  rr3k = bn(1) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = bn(2) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = bn(3) - (1.0d0-scalek*dmpk(7))*rr7
                  rr1ik = bn(0) - (1.0d0-scalek*dmpik(1))*rr1
                  rr3ik = bn(1) - (1.0d0-scalek*dmpik(3))*rr3
                  rr5ik = bn(2) - (1.0d0-scalek*dmpik(5))*rr5
                  rr7ik = bn(3) - (1.0d0-scalek*dmpik(7))*rr7
                  rr9ik = bn(4) - (1.0d0-scalek*dmpik(9))*rr9
                  rr11ik = bn(5) - (1.0d0-scalek*dmpik(11))*rr11
                  rr1 = bn(0) - (1.0d0-scalek)*rr1
                  rr3 = bn(1) - (1.0d0-scalek)*rr3
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find damped multipole intermediates for force and torque
c
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
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  scalek = 1.0d0 - mscale(k)
                  rr1 = bn(0) - scalek*rr1
                  rr3 = bn(1) - scalek*rr3
                  rr5 = bn(2) - scalek*rr5
                  rr7 = bn(3) - scalek*rr7
                  rr9 = bn(4) - scalek*rr9
                  rr11 = bn(5) - scalek*rr11
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
c
c     find standard multipole intermediates for force and torque
c
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
c     energy, force and torque scaled for self-interactions
c
               if (i .eq. k) then
                  e = 0.5d0 * e
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  do j = 1, 3
                     ttmi(j) = 0.5d0 * ttmi(j)
                     ttmk(j) = 0.5d0 * ttmk(j)
                  end do
               end if
c
c     increment the overall atomic multipole energy component
c
               em = em + e
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
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -0.5d0 * (yr*frcx+xr*frcy)
               vxz = -0.5d0 * (zr*frcx+xr*frcz)
               vyy = -yr * frcy
               vyz = -0.5d0 * (zr*frcy+yr*frcz)
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
c     resolve site torques then increment forces and virial
c
      do ii = 1, npole
         i = ipole(ii)
         call torque (ii,tem(1,i),fix,fiy,fiz,dem)
         iaz = zaxis(ii)
         iax = xaxis(ii)
         iay = yaxis(ii)
         if (iaz .eq. 0)  iaz = i
         if (iax .eq. 0)  iax = i
         if (iay .eq. 0)  iay = i
         xiz = x(iaz) - x(i)
         yiz = y(iaz) - y(i)
         ziz = z(iaz) - z(i)
         xix = x(iax) - x(i)
         yix = y(iax) - y(i)
         zix = z(iax) - z(i)
         xiy = x(iay) - x(i)
         yiy = y(iay) - y(i)
         ziy = z(iay) - z(i)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
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
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use pme
      use virial
      implicit none
      integer i,j,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 tem(3),frcx(3)
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
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bseorder
      aewald = aeewald
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
      call emreal1d
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip1
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do ii = 1, npole
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
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
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
         do ii = 1, npole
            i = ipole(ii)
            xd = xd + rpole(2,ii) + rpole(1,ii)*x(i)
            yd = yd + rpole(3,ii) + rpole(1,ii)*y(i)
            zd = zd + rpole(4,ii) + rpole(1,ii)*z(i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         do ii = 1, npole
            i = ipole(ii)
            dem(1,i) = dem(1,i) + 2.0d0*term*rpole(1,ii)*xd
            dem(2,i) = dem(2,i) + 2.0d0*term*rpole(1,ii)*yd
            dem(3,i) = dem(3,i) + 2.0d0*term*rpole(1,ii)*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         do ii = 1, npole
            tem(1) = rpole(3,ii)*zdfield - rpole(4,ii)*ydfield
            tem(2) = rpole(4,ii)*xdfield - rpole(2,ii)*zdfield
            tem(3) = rpole(2,ii)*ydfield - rpole(3,ii)*xdfield
            call torque (ii,tem,frcx,frcy,frcz,dem)
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
         do ii = 1, npole
            i = ipole(ii)
            xd = xd + rpole(2,ii)
            yd = yd + rpole(3,ii)
            zd = zd + rpole(4,ii)
            xq = xq + rpole(1,ii)*x(i)
            yq = yq + rpole(1,ii)*y(i)
            zq = zq + rpole(1,ii)*z(i)
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
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine emreal1d  --  Ewald real mpole derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "emreal1d" evaluates the real space portion of the Ewald
c     summation energy and gradient due to multipole interactions
c     via a neighbor list
c
c
      subroutine emreal1d
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use math
      use mplpot
      use mpole
      use neigh
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer iax,iay,iaz
      real*8 e,de,f
      real*8 bfac,erfc
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
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
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 ttmi(3),ttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(11)
      real*8 bn(0:5)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      character*6 mode
      external erfc
c
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
      mode = 'EWALD'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,rpole,pcore,pval,palpha,n12,i12,n13,
!$OMP& i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,nelst,
!$OMP& elst,use_chgpen,use_bounds,f,off2,aewald,xaxis,yaxis,zaxis)
!$OMP& firstprivate(mscale) shared (em,dem,tem,vir)
!$OMP DO reduction(+:em,dem,tem,vir) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npole
         i = ipole(ii)
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
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
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
c     find damped multipole intermediates and energy value
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
                  scalek = mscale(k)
                  rr1i = bn(0) - (1.0d0-scalek*dmpi(1))*rr1
                  rr3i = bn(1) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = bn(2) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = bn(3) - (1.0d0-scalek*dmpi(7))*rr7
                  rr1k = bn(0) - (1.0d0-scalek*dmpk(1))*rr1
                  rr3k = bn(1) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = bn(2) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = bn(3) - (1.0d0-scalek*dmpk(7))*rr7
                  rr1ik = bn(0) - (1.0d0-scalek*dmpik(1))*rr1
                  rr3ik = bn(1) - (1.0d0-scalek*dmpik(3))*rr3
                  rr5ik = bn(2) - (1.0d0-scalek*dmpik(5))*rr5
                  rr7ik = bn(3) - (1.0d0-scalek*dmpik(7))*rr7
                  rr9ik = bn(4) - (1.0d0-scalek*dmpik(9))*rr9
                  rr11ik = bn(5) - (1.0d0-scalek*dmpik(11))*rr11
                  rr1 = bn(0) - (1.0d0-scalek)*rr1
                  rr3 = bn(1) - (1.0d0-scalek)*rr3
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find damped multipole intermediates for force and torque
c
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
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  scalek = 1.0d0 - mscale(k)
                  rr1 = bn(0) - scalek*rr1
                  rr3 = bn(1) - scalek*rr3
                  rr5 = bn(2) - scalek*rr5
                  rr7 = bn(3) - scalek*rr7
                  rr9 = bn(4) - scalek*rr9
                  rr11 = bn(5) - scalek*rr11
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
c
c     find standard multipole intermediates for force and torque
c
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
c     increment the overall atomic multipole energy component
c
               em = em + e
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
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = -xr * frcx
               vxy = -0.5d0 * (yr*frcx+xr*frcy)
               vxz = -0.5d0 * (zr*frcx+xr*frcz)
               vyy = -yr * frcy
               vyz = -0.5d0 * (zr*frcy+yr*frcz)
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO reduction(+:dem,vir) schedule(guided)
c
c     resolve site torques then increment forces and virial
c
      do ii = 1, npole
         i = ipole(ii)
         call torque (ii,tem(1,i),fix,fiy,fiz,dem)
         iaz = zaxis(ii)
         iax = xaxis(ii)
         iay = yaxis(ii)
         if (iaz .eq. 0)  iaz = i
         if (iax .eq. 0)  iax = i
         if (iay .eq. 0)  iay = i
         xiz = x(iaz) - x(i)
         yiz = y(iaz) - y(i)
         ziz = z(iaz) - z(i)
         xix = x(iax) - x(i)
         yix = y(iax) - y(i)
         zix = z(iax) - z(i)
         xiy = x(iay) - x(i)
         yiy = y(iay) - y(i)
         ziy = z(iay) - z(i)
         vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
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
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine emrecip1  --  PME recip mpole energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "emrecip1" evaluates the reciprocal space portion of particle
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
      integer iax,iay,iaz
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm,f
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 tem(3),fix(3)
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
      f = electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(cmp) .and. size(cmp).lt.10*npole)
     &   deallocate (cmp)
      if (allocated(fmp) .and. size(fmp).lt.10*npole)
     &   deallocate (fmp)
      if (allocated(cphi) .and. size(cphi).lt.10*npole)
     &   deallocate (cphi)
      if (allocated(fphi) .and. size(fphi).lt.20*npole)
     &   deallocate (fphi)
      if (.not. allocated(cmp))  allocate (cmp(10,npole))
      if (.not. allocated(fmp))  allocate (fmp(10,npole))
      if (.not. allocated(cphi))  allocate (cphi(10,npole))
      if (.not. allocated(fphi))  allocate (fphi(20,npole))
c
c     perform dynamic allocation of some global arrays
c
      ntot = nfft1 * nfft2 * nfft3
      if (allocated(qgrid) .and. size(qgrid).ne.2*ntot)
     &   deallocate(qgrid)
      if (allocated(qfac) .and. size(qfac).ne.ntot)
     &   deallocate(qfac)
      if (.not. allocated(qgrid))
     &   allocate (qgrid(2,nfft1,nfft2,nfft3))
      if (.not. allocated(qfac))
     &   allocate (qfac(nfft1,nfft2,nfft3))
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
c     setup spatial decomposition, B-splines and PME arrays
c
      call getchunk
      call moduli
      call fftsetup
c
c     compute B-spline coefficients and spatial decomposition
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
      qfac(1,1,1) = 0.0d0
      pterm = (pi/aewald)**2
      volterm = pi * volbox
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
            eterm = 0.5d0 * f * expterm * struc2
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
c     account for zeroth grid point for nonperiodic system
c
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qfac(1,1,1) = expterm
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = f * expterm * struc2
         em = em + e
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
            fphi(j,i) = f * fphi(j,i)
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
c     increment the permanent multipole virial contributions
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
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         tem(1) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
     &               + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
     &               + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
     &               - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
         tem(2) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
     &               + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
     &               + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
     &               - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
         tem(3) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
     &               + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
     &               + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
     &               - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)
         call torque (i,tem,fix,fiy,fiz,dem)
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
         vxx = vxx + xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = vxy + 0.5d0*(yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                        + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = vxz + 0.5d0*(zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                        + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3)) 
         vyy = vyy + yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = vyz + 0.5d0*(zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                        + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = vzz + zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
      end do
c
c     increment the total internal virial tensor components
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

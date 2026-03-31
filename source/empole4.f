c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses Chung, Pengyu Ren, Jay Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine empole4  --  multipole lambda derivatives  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "empole4" calculates the atomic multipole energy,
c     first derivatives with respect to Cartesian coordinates,
c     and lambda derivatives
c
c
      subroutine empole4
      use energi
      use extfld
      use limits
      implicit none
      real*8 exf
      character*6 mode
c
c
c     choose the method to sum over multipole interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call empole4d
         else
            call empole4c
         end if
      else
         if (use_mlist) then
            call empole4b
         else
            call empole4a
         end if
      end if
c
c     get contribution from external electric field if used
c
      if (use_exfld) then
         mode = 'MPOLE'
         call exfield1 (mode,exf)
         em = em + exf
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empole4a  --  double loop mpole lambda derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole4a" calculates the multipole energy, derivatives with 
c     respect to Cartesian coordinates, and lambda derivatives
c     using a pairwise double loop
c
c
      subroutine empole4a
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use deriv
      use dlmda
      use energi
      use group
      use mplpot
      use mpole
      use mutant
      use potent
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer ix,iy,iz
      integer kx,ky,kz
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
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 poti,potk
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 dlde
      real*8 dlfrcx,dlfrcy,dlfrcz
      real*8 dlambda,dlambda2
      real*8 scalelmda
      real*8 ttmi(3),ttmk(3)
      real*8 dlttmi(3),dlttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(11)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      real*8, allocatable :: dltem(:,:)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
      logical proceed,usei,usek
      logical muti,mutk
      character*6 mode
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      demlmda = 0.0d0
      demlmda2 = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dldem(j,i) = 0.0d0
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
      call rotpole ('MPORG')
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
      allocate (dltem(3,n))
      allocate (pot(n))
      allocate (decfx(n))
      allocate (decfy(n))
      allocate (decfz(n))
c
c     initialize scaling, torque and potential arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
            dltem(j,i) = 0.0d0
         end do
         pot(i) = 0.0d0
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
         iz = zaxis(i)
         ix = xaxis(i)
         iy = abs(yaxis(i))
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(i)
            vali = pval(i)
            alphai = palpha(i)
         end if
         usei = (use(i) .or. use(iz) .or. use(ix) .or. use(iy))
         muti = mut(i)
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
            kz = zaxis(k)
            kx = xaxis(k)
            ky = abs(yaxis(k))
            usek = (use(k) .or. use(kz) .or. use(kx) .or. use(ky))
            mutk = mut(k)
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
                  corek = pcore(k)
                  valk = pval(k)
                  alphak = palpha(k)
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
                  term4 = -2.0d0 * (ck*rr5-dkr*rr7+qkr*rr9)
                  term5 = -2.0d0 * (ci*rr5+dir*rr7+qir*rr9)
                  term6 = 4.0d0 * rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_chgpen) then
                     term1i = corek*dmpi(1) + valk*dmpik(1) 
                     term1k = corei*dmpk(1) + vali*dmpik(1) 
                     term2i = -dkr * dmpik(3)
                     term2k = dir * dmpik(3)
                     term3i = qkr * dmpik(5)
                     term3k = qir * dmpik(5)
                     poti = term1i*rr1 + term2i*rr3 + term3i*rr5
                     potk = term1k*rr1 + term2k*rr3 + term3k*rr5
                  else
                     poti = ck*rr1 - dkr*rr3 + qkr*rr5
                     potk = ci*rr1 + dir*rr3 + qir*rr5
                  end if
                  pot(i) = pot(i) + poti
                  pot(k) = pot(k) + potk
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
               if (use_chgpen)  rr3 = rr3ik
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
c     compute lambda derivative
c
               scalelmda = 1.0d0
               if (muti .and. mutk) then
                  dlambda = 2.0d0 * elambda * e
                  dlambda2 = 2.0d0 * e
                  dlfrcx = 2.0d0 * elambda * frcx
                  dlfrcy = 2.0d0 * elambda * frcy
                  dlfrcz = 2.0d0 * elambda * frcz
                  dlttmi(1) = 2.0d0 * elambda * ttmi(1)
                  dlttmi(2) = 2.0d0 * elambda * ttmi(2)
                  dlttmi(3) = 2.0d0 * elambda * ttmi(3)
                  dlttmk(1) = 2.0d0 * elambda * ttmk(1)
                  dlttmk(2) = 2.0d0 * elambda * ttmk(2)
                  dlttmk(3) = 2.0d0 * elambda * ttmk(3)
                  scalelmda = elambda * elambda
               else if (muti .or. mutk) then
                  dlambda = e
                  dlambda2 = 0.0d0
                  dlfrcx = frcx
                  dlfrcy = frcy
                  dlfrcz = frcz
                  dlttmi(1) = ttmi(1)
                  dlttmi(2) = ttmi(2)
                  dlttmi(3) = ttmi(3)
                  dlttmk(1) = ttmk(1)
                  dlttmk(2) = ttmk(2)
                  dlttmk(3) = ttmk(3)
                  scalelmda = elambda
               end if
               if (muti .or. mutk) then
                  demlmda = demlmda + dlambda
                  demlmda2 = demlmda2 + dlambda2
               end if
c
c     modify the energy, force, and torque by lambda
c
               e = scalelmda * e
               frcx = scalelmda * frcx
               frcy = scalelmda * frcy
               frcz = scalelmda * frcz
               ttmi(1) = scalelmda * ttmi(1)
               ttmi(2) = scalelmda * ttmi(2)
               ttmi(3) = scalelmda * ttmi(3)
               ttmk(1) = scalelmda * ttmk(1)
               ttmk(2) = scalelmda * ttmk(2)
               ttmk(3) = scalelmda * ttmk(3)
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
               if (muti .or. mutk) then
                  dldem(1,i) = dldem(1,i) + dlfrcx
                  dldem(2,i) = dldem(2,i) + dlfrcy
                  dldem(3,i) = dldem(3,i) + dlfrcz
                  dltem(1,i) = dltem(1,i) + dlttmi(1)
                  dltem(2,i) = dltem(2,i) + dlttmi(2)
                  dltem(3,i) = dltem(3,i) + dlttmi(3)
               end if
c
c     increment force-based gradient and torque on second site
c
               dem(1,k) = dem(1,k) - frcx
               dem(2,k) = dem(2,k) - frcy
               dem(3,k) = dem(3,k) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
               if (muti .or. mutk) then
                  dldem(1,k) = dldem(1,k) - dlfrcx
                  dldem(2,k) = dldem(2,k) - dlfrcy
                  dldem(3,k) = dldem(3,k) - dlfrcz
                  dltem(1,k) = dltem(1,k) + dlttmk(1)
                  dltem(2,k) = dltem(2,k) + dlttmk(2)
                  dltem(3,k) = dltem(3,k) + dlttmk(3)
               end if
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
         iz = zaxis(i)
         ix = xaxis(i)
         iy = abs(yaxis(i))
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(i)
            vali = pval(i)
            alphai = palpha(i)
         end if
         usei = (use(i) .or. use(iz) .or. use(ix) .or. use(iy))
         muti = mut(i)
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
         do kk = ii, npole
            k = ipole(kk)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = abs(yaxis(k))
            usek = (use(k) .or. use(kz) .or. use(kx) .or. use(ky))
            mutk = mut(k)
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
                  corek = pcore(k)
                  valk = pval(k)
                  alphak = palpha(k)
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
                  term4 = -2.0d0 * (ck*rr5-dkr*rr7+qkr*rr9)
                  term5 = -2.0d0 * (ci*rr5+dir*rr7+qir*rr9)
                  term6 = 4.0d0 * rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_chgpen) then
                     term1i = corek*dmpi(1) + valk*dmpik(1) 
                     term1k = corei*dmpk(1) + vali*dmpik(1) 
                     term2i = -dkr * dmpik(3)
                     term2k = dir * dmpik(3)
                     term3i = qkr * dmpik(5)
                     term3k = qir * dmpik(5)
                     poti = term1i*rr1 + term2i*rr3 + term3i*rr5
                     potk = term1k*rr1 + term2k*rr3 + term3k*rr5
                  else
                     poti = ck*rr1 - dkr*rr3 + qkr*rr5
                     potk = ci*rr1 + dir*rr3 + qir*rr5
                  end if
                  pot(i) = pot(i) + poti
                  pot(k) = pot(k) + potk
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
               if (use_chgpen)  rr3 = rr3ik
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
c     compute lambda derivative
c
               scalelmda = 1.0d0
               if (muti .and. mutk) then
                  dlambda = 2.0d0 * elambda * e
                  dlambda2 = 2.0d0 * e
                  dlfrcx = 2.0d0 * elambda * frcx
                  dlfrcy = 2.0d0 * elambda * frcy
                  dlfrcz = 2.0d0 * elambda * frcz
                  dlttmi(1) = 2.0d0 * elambda * ttmi(1)
                  dlttmi(2) = 2.0d0 * elambda * ttmi(2)
                  dlttmi(3) = 2.0d0 * elambda * ttmi(3)
                  dlttmk(1) = 2.0d0 * elambda * ttmk(1)
                  dlttmk(2) = 2.0d0 * elambda * ttmk(2)
                  dlttmk(3) = 2.0d0 * elambda * ttmk(3)
                  scalelmda = elambda * elambda
               else if (muti .or. mutk) then
                  dlambda = e
                  dlambda2 = 0.0d0
                  dlfrcx = frcx
                  dlfrcy = frcy
                  dlfrcz = frcz
                  dlttmi(1) = ttmi(1)
                  dlttmi(2) = ttmi(2)
                  dlttmi(3) = ttmi(3)
                  dlttmk(1) = ttmk(1)
                  dlttmk(2) = ttmk(2)
                  dlttmk(3) = ttmk(3)
                  scalelmda = elambda
               end if
               if (muti .or. mutk) then
                  if (i .eq. k) then
                     dlambda = 0.5d0 * dlambda
                     dlambda2 = 0.5d0 * dlambda2
                  end if
                  demlmda = demlmda + dlambda
                  demlmda2 = demlmda2 + dlambda2
               end if
c
c     modify the energy, force, and torque by lambda
c
               e = scalelmda * e
               frcx = scalelmda * frcx
               frcy = scalelmda * frcy
               frcz = scalelmda * frcz
               ttmi(1) = scalelmda * ttmi(1)
               ttmi(2) = scalelmda * ttmi(2)
               ttmi(3) = scalelmda * ttmi(3)
               ttmk(1) = scalelmda * ttmk(1)
               ttmk(2) = scalelmda * ttmk(2)
               ttmk(3) = scalelmda * ttmk(3)
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
               if (muti .or. mutk) then
                  dldem(1,i) = dldem(1,i) + dlfrcx
                  dldem(2,i) = dldem(2,i) + dlfrcy
                  dldem(3,i) = dldem(3,i) + dlfrcz
                  dltem(1,i) = dltem(1,i) + dlttmi(1)
                  dltem(2,i) = dltem(2,i) + dlttmi(2)
                  dltem(3,i) = dltem(3,i) + dlttmi(3)
               end if
c
c     increment force-based gradient and torque on second site
c
               dem(1,k) = dem(1,k) - frcx
               dem(2,k) = dem(2,k) - frcy
               dem(3,k) = dem(3,k) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
               if (muti .or. mutk) then
                  dldem(1,k) = dldem(1,k) - dlfrcx
                  dldem(2,k) = dldem(2,k) - dlfrcy
                  dldem(3,k) = dldem(3,k) - dlfrcz
                  dltem(1,k) = dltem(1,k) + dlttmk(1)
                  dltem(2,k) = dltem(2,k) + dlttmk(2)
                  dltem(3,k) = dltem(3,k) + dlttmk(3)
               end if
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
         call torque (i,dltem(1,i),fix,fiy,fiz,dldem)
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = abs(yaxis(i))
         if (iz .eq. 0)  iz = i
         if (ix .eq. 0)  ix = i
         if (iy .eq. 0)  iy = i
         xiz = x(iz) - x(i)
         yiz = y(iz) - y(i)
         ziz = z(iz) - z(i)
         xix = x(ix) - x(i)
         yix = y(ix) - y(i)
         zix = z(ix) - z(i)
         xiy = x(iy) - x(i)
         yiy = y(iy) - y(i)
         ziy = z(iy) - z(i)
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
c     modify the gradient and virial for charge flux
c
      if (use_chgflx) then
         call dcflux (pot,decfx,decfy,decfz)
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            frcx = decfx(i)
            frcy = decfy(i)
            frcz = decfz(i)
            dem(1,i) = dem(1,i) + frcx
            dem(2,i) = dem(2,i) + frcy
            dem(3,i) = dem(3,i) + frcz
            vxx = xi * frcx
            vxy = yi * frcx
            vxz = zi * frcx
            vyy = yi * frcy
            vyz = zi * frcy
            vzz = zi * frcz
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
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      deallocate (dltem)
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole4b  --  neighbor list mpole lambda derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole4b" calculates the multipole energy, derivatives
c     with respect to Cartesian coordinates, and lambda derivatives
c     using a neighbor list
c
c
      subroutine empole4b
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use deriv
      use dlmda
      use energi
      use group
      use mplpot
      use mpole
      use mutant
      use neigh
      use potent
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
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
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 poti,potk
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 dlde
      real*8 dlfrcx,dlfrcy,dlfrcz
      real*8 dlambda,dlambda2
      real*8 scalelmda
      real*8 ttmi(3),ttmk(3)
      real*8 dlttmi(3),dlttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(11)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      real*8, allocatable :: dltem(:,:)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
      logical proceed,usei,usek
      logical muti,mutk
      character*6 mode
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      demlmda = 0.0d0
      demlmda2 = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dldem(j,i) = 0.0d0
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
      call rotpole ('MPORG')
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
      allocate (pot(n))
      allocate (dltem(3,n))
      allocate (decfx(n))
      allocate (decfy(n))
      allocate (decfz(n))
c
c     initialize scaling, torque and potential arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
            dltem(j,i) = 0.0d0
         end do
         pot(i) = 0.0d0
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
!$OMP& m3scale,m4scale,m5scale,nelst,elst,use_chgpen,use_chgflx,
!$OMP& use_group,use_intra,use_bounds,off2,f,mut,elambda)
!$OMP& firstprivate(mscale) shared (em,dem,dldem,tem,dltem,pot,vir,
!$OMP& demlmda,demlmda2)
!$OMP DO reduction(+:em,dem,dldem,tem,dltem,pot,vir,demlmda,demlmda2)
c
c     compute the multipole interaction energy and gradient
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = abs(yaxis(i))
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(i)
            vali = pval(i)
            alphai = palpha(i)
         end if
         usei = (use(i) .or. use(iz) .or. use(ix) .or. use(iy))
         muti = mut(i)
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
            kz = zaxis(k)
            kx = xaxis(k)
            ky = abs(yaxis(k))
            usek = (use(k) .or. use(kz) .or. use(kx) .or. use(ky))
            mutk = mut(k)
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
                  corek = pcore(k)
                  valk = pval(k)
                  alphak = palpha(k)
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
                  term4 = -2.0d0 * (ck*rr5-dkr*rr7+qkr*rr9)
                  term5 = -2.0d0 * (ci*rr5+dir*rr7+qir*rr9)
                  term6 = 4.0d0 * rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_chgpen) then
                     term1i = corek*dmpi(1) + valk*dmpik(1) 
                     term1k = corei*dmpk(1) + vali*dmpik(1) 
                     term2i = -dkr * dmpik(3)
                     term2k = dir * dmpik(3)
                     term3i = qkr * dmpik(5)
                     term3k = qir * dmpik(5)
                     poti = term1i*rr1 + term2i*rr3 + term3i*rr5
                     potk = term1k*rr1 + term2k*rr3 + term3k*rr5
                  else
                     poti = ck*rr1 - dkr*rr3 + qkr*rr5
                     potk = ci*rr1 + dir*rr3 + qir*rr5
                  end if
                  pot(i) = pot(i) + poti
                  pot(k) = pot(k) + potk
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
               if (use_chgpen)  rr3 = rr3ik
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
c     compute lambda derivative
c
               scalelmda = 1.0d0
               if (muti .and. mutk) then
                  dlambda = 2.0d0 * elambda * e
                  dlambda2 = 2.0d0 * e
                  dlfrcx = 2.0d0 * elambda * frcx
                  dlfrcy = 2.0d0 * elambda * frcy
                  dlfrcz = 2.0d0 * elambda * frcz
                  dlttmi(1) = 2.0d0 * elambda * ttmi(1)
                  dlttmi(2) = 2.0d0 * elambda * ttmi(2)
                  dlttmi(3) = 2.0d0 * elambda * ttmi(3)
                  dlttmk(1) = 2.0d0 * elambda * ttmk(1)
                  dlttmk(2) = 2.0d0 * elambda * ttmk(2)
                  dlttmk(3) = 2.0d0 * elambda * ttmk(3)
                  scalelmda = elambda * elambda
               else if (muti .or. mutk) then
                  dlambda = e
                  dlambda2 = 0.0d0
                  dlfrcx = frcx
                  dlfrcy = frcy
                  dlfrcz = frcz
                  dlttmi(1) = ttmi(1)
                  dlttmi(2) = ttmi(2)
                  dlttmi(3) = ttmi(3)
                  dlttmk(1) = ttmk(1)
                  dlttmk(2) = ttmk(2)
                  dlttmk(3) = ttmk(3)
                  scalelmda = elambda
               end if
               if (muti .or. mutk) then
                  demlmda = demlmda + dlambda
                  demlmda2 = demlmda2 + dlambda2
               end if
c
c     modify the energy, force, and torque by lambda
c
               e = scalelmda * e
               frcx = scalelmda * frcx
               frcy = scalelmda * frcy
               frcz = scalelmda * frcz
               ttmi(1) = scalelmda * ttmi(1)
               ttmi(2) = scalelmda * ttmi(2)
               ttmi(3) = scalelmda * ttmi(3)
               ttmk(1) = scalelmda * ttmk(1)
               ttmk(2) = scalelmda * ttmk(2)
               ttmk(3) = scalelmda * ttmk(3)
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
               if (muti .or. mutk) then
                  dldem(1,i) = dldem(1,i) + dlfrcx
                  dldem(2,i) = dldem(2,i) + dlfrcy
                  dldem(3,i) = dldem(3,i) + dlfrcz
                  dltem(1,i) = dltem(1,i) + dlttmi(1)
                  dltem(2,i) = dltem(2,i) + dlttmi(2)
                  dltem(3,i) = dltem(3,i) + dlttmi(3)
               end if
c
c     increment force-based gradient and torque on second site
c
               dem(1,k) = dem(1,k) - frcx
               dem(2,k) = dem(2,k) - frcy
               dem(3,k) = dem(3,k) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
               if (muti .or. mutk) then
                  dldem(1,k) = dldem(1,k) - dlfrcx
                  dldem(2,k) = dldem(2,k) - dlfrcy
                  dldem(3,k) = dldem(3,k) - dlfrcz
                  dltem(1,k) = dltem(1,k) + dlttmk(1)
                  dltem(2,k) = dltem(2,k) + dlttmk(2)
                  dltem(3,k) = dltem(3,k) + dlttmk(3)
               end if
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
!$OMP DO reduction(+:dem,dldem,vir)
c
c     resolve site torques then increment forces and virial
c
      do ii = 1, npole
         i = ipole(ii)
         call torque (i,dltem(1,i),fix,fiy,fiz,dldem)
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = abs(yaxis(i))
         if (iz .eq. 0)  iz = i
         if (ix .eq. 0)  ix = i
         if (iy .eq. 0)  iy = i
         xiz = x(iz) - x(i)
         yiz = y(iz) - y(i)
         ziz = z(iz) - z(i)
         xix = x(ix) - x(i)
         yix = y(ix) - y(i)
         zix = z(ix) - z(i)
         xiy = x(iy) - x(i)
         yiy = y(iy) - y(i)
         ziy = z(iy) - z(i)
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
c
c     modify the gradient and virial for charge flux
c
      if (use_chgflx) then
         call dcflux (pot,decfx,decfy,decfz)
!$OMP    DO reduction(+:dem,vir)
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            frcx = decfx(i)
            frcy = decfy(i)
            frcz = decfz(i)
            dem(1,i) = dem(1,i) + frcx
            dem(2,i) = dem(2,i) + frcy
            dem(3,i) = dem(3,i) + frcz
            vxx = xi * frcx
            vxy = yi * frcx
            vxz = zi * frcx
            vyy = yi * frcy
            vyz = zi * frcy
            vzz = zi * frcz
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
!$OMP    END DO
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      deallocate (dltem)
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole4c  --  Ewald lambda derivs via loop  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole4c" calculates the multipole energy, derivatives
c     with respect to Cartesian coordinates, and lambda derivatives
c     using particle mesh Ewald summation and a double loop
c
c
      subroutine empole4c
      use atoms
      use boxes
      use chgpot
      use deriv
      use dlmda
      use energi
      use ewald
      use math
      use mpole
      use mutant
      use pme
      use potent
      use virial
      implicit none
      integer i,j,ii
      real*8 e,f,sum,dlsum
      real*8 dlxd,dlyd,dlzd
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xi,yi,zi
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 dlxdfield,dlydfield
      real*8 dlzdfield
      real*8 fx,fy,fz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 tem(3),frcx(3)
      real*8 frcy(3),frcz(3)
      real*8 dltem(3)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
      logical muti
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      demlmda = 0.0d0
      demlmda2 = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dldem(j,i) = 0.0d0
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
      call rotpole ('MPORG')
c
c     compute the real space part of the Ewald summation
c
      call emreal4c
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip4
c
c     perform dynamic allocation of some local arrays
c
      allocate (pot(n))
      allocate (decfx(n))
      allocate (decfy(n))
      allocate (decfz(n))
c
c     initialize Ewald self-energy potential array
c
      do i = 1, n
         pot(i) = 0.0d0
      end do
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / rootpi
      do ii = 1, npole
         i = ipole(ii)
         muti = mut(i)
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
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         if (muti) then
            demlmda = demlmda + 2.0d0 * elambda * e
            demlmda2 = demlmda2 + 2.0d0 * e
            e = e * elambda * elambda
         end if
         em = em + e
         pot(i) = 2.0d0 * fterm * ci
      end do
c
c     modify gradient and virial for charge flux self-energy
c
      if (use_chgflx) then
         call dcflux (pot,decfx,decfy,decfz)
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            fx = decfx(i)
            fy = decfy(i)
            fz = decfz(i)
            dem(1,i) = dem(1,i) + fx
            dem(2,i) = dem(2,i) + fy
            dem(3,i) = dem(3,i) + fz
            vxx = xi * fx
            vxy = yi * fx
            vxz = zi * fx
            vyy = yi * fy
            vyz = zi * fy
            vzz = zi * fz
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
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
c
c     compute the uniform background charge correction term
c
      fterm = -0.5d0 * f * pi / (volbox*aewald**2)
      sum = 0.0d0
      dlsum = 0.0d0
      do ii = 1, npole
         i = ipole(ii)
         ci = rpole(1,i)
         muti = mut(i)
         if (muti) then
            dlsum = dlsum + ci
            ci = ci * elambda
         end if
         sum = sum + ci
      end do
      e = fterm * sum**2
      em = em + e
      demlmda = demlmda + fterm * 2.0d0 * sum * dlsum
      demlmda2 = demlmda2 + fterm * 2.0d0 * dlsum**2
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         dlxd = 0.0d0
         dlyd = 0.0d0
         dlzd = 0.0d0
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            muti = mut(i)
            if (muti) then
               dlxd = dlxd + dix + ci*xi
               dlyd = dlyd + diy + ci*yi
               dlzd = dlzd + diz + ci*zi
               ci = ci * elambda
               dix = dix * elambda
               diy = diy * elambda
               diz = diz * elambda
            end if
            xd = xd + dix + ci*xi
            yd = yd + diy + ci*yi
            zd = zd + diz + ci*zi
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         demlmda = demlmda + 2.0d0*term*(xd*dlxd+yd*dlyd+zd*dlzd)
         demlmda2 = demlmda2 + 2.0d0*term*(dlxd**2+dlyd**2+dlzd**2)
         do ii = 1, npole
            i = ipole(ii)
            ci = rpole(1,i)
            muti = mut(i)
            if (muti) then
               dldem(1,i) = dldem(1,i) + 2.0d0*term*ci*(xd+elambda*dlxd)
               dldem(2,i) = dldem(2,i) + 2.0d0*term*ci*(yd+elambda*dlyd)
               dldem(3,i) = dldem(3,i) + 2.0d0*term*ci*(zd+elambda*dlzd)
               ci = ci * elambda
            else
               dldem(1,i) = dldem(1,i) + 2.0d0*term*ci*dlxd
               dldem(2,i) = dldem(2,i) + 2.0d0*term*ci*dlyd
               dldem(3,i) = dldem(3,i) + 2.0d0*term*ci*dlzd
            end if
            dem(1,i) = dem(1,i) + 2.0d0*term*ci*xd
            dem(2,i) = dem(2,i) + 2.0d0*term*ci*yd
            dem(3,i) = dem(3,i) + 2.0d0*term*ci*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         dlxdfield = -2.0d0 * term * dlxd
         dlydfield = -2.0d0 * term * dlyd
         dlzdfield = -2.0d0 * term * dlzd
         do ii = 1, npole
            i = ipole(ii)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            muti = mut(i)
            if (muti) then
               dltem(1) = diy*(zdfield+elambda*dlzdfield)
     &                    - diz*(ydfield+elambda*dlydfield)
               dltem(2) = diz*(xdfield+elambda*dlxdfield)
     &                    - dix*(zdfield+elambda*dlzdfield)
               dltem(3) = dix*(ydfield+elambda*dlydfield)
     &                    - diy*(xdfield+elambda*dlxdfield)
               dix = dix * elambda
               diy = diy * elambda
               diz = diz * elambda
            else
               dltem(1) = diy*dlzdfield - diz*dlydfield
               dltem(2) = diz*dlxdfield - dix*dlzdfield
               dltem(3) = dix*dlydfield - diy*dlxdfield
            end if
            tem(1) = diy*zdfield - diz*ydfield
            tem(2) = diz*xdfield - dix*zdfield
            tem(3) = dix*ydfield - diy*xdfield
            call torque (i,dltem,frcx,frcy,frcz,dldem)
            call torque (i,tem,frcx,frcy,frcz,dem)
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
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            muti = mut(i)
            if (muti) then
               ci = ci * elambda
               dix = dix * elambda
               diy = diy * elambda
               diz = diz * elambda
            end if
            xd = xd + dix
            yd = yd + diy
            zd = zd + diz
            xq = xq + ci*x(i)
            yq = yq + ci*y(i)
            zq = zq + ci*z(i)
         end do
         xv = xd * xq
         yv = yd * yq
         zv = zd * zq
         vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                      + xq*xq + yq*yq + zq*zq)
         vxx = 2.0d0*term*(xq*xq+xv) + vterm
         vxy = 2.0d0*term*(xq*yq+xv)
         vxz = 2.0d0*term*(xq*zq+xv)
         vyy = 2.0d0*term*(yq*yq+yv) + vterm
         vyz = 2.0d0*term*(yq*zq+yv)
         vzz = 2.0d0*term*(zq*zq+zv) + vterm
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
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emreal4c  --  Ewald real lambda derivs via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emreal4c" evaluates the real space portion of the Ewald
c     summation energy, gradient, and lambda derivatives due to
c     multipole interactions via a double loop
c
c
      subroutine emreal4c
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use deriv
      use dlmda
      use energi
      use math
      use mplpot
      use mpole
      use mutant
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer ix,iy,iz
      real*8 e,de,f
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
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 poti,potk
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 dlde
      real*8 dlfrcx,dlfrcy,dlfrcz
      real*8 dlambda,dlambda2
      real*8 scalelmda
      real*8 ttmi(3),ttmk(3)
      real*8 dlttmi(3),dlttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(11),dmpe(11)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      real*8, allocatable :: dltem(:,:)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
      logical muti,mutk
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
      allocate (dltem(3,n))
      allocate (pot(n))
      allocate (decfx(n))
      allocate (decfy(n))
      allocate (decfz(n))
c
c     initialize scaling, torque and potential arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
            dltem(j,i) = 0.0d0
         end do
         pot(i) = 0.0d0
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
         if (use_chgpen) then
            corei = pcore(i)
            vali = pval(i)
            alphai = palpha(i)
         end if
         muti = mut(i)
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
            mutk = mut(k)
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
c     calculate real space Ewald error function damping
c
               call dampewald (11,r,r2,f,dmpe)
c
c     find damped multipole intermediates and energy value
c
               if (use_chgpen) then
                  corek = pcore(k)
                  valk = pval(k)
                  alphak = palpha(k)
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
                  rr1i = dmpe(1) - (1.0d0-scalek*dmpi(1))*rr1
                  rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
                  rr1k = dmpe(1) - (1.0d0-scalek*dmpk(1))*rr1
                  rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
                  rr1ik = dmpe(1) - (1.0d0-scalek*dmpik(1))*rr1
                  rr3ik = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
                  rr5ik = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
                  rr9ik = dmpe(9) - (1.0d0-scalek*dmpik(9))*rr9
                  rr11ik = dmpe(11) - (1.0d0-scalek*dmpik(11))*rr11
                  rr1 = dmpe(1) - (1.0d0-scalek)*rr1
                  rr3 = dmpe(3) - (1.0d0-scalek)*rr3
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
                  rr1 = dmpe(1) - scalek*rr1
                  rr3 = dmpe(3) - scalek*rr3
                  rr5 = dmpe(5) - scalek*rr5
                  rr7 = dmpe(7) - scalek*rr7
                  rr9 = dmpe(9) - scalek*rr9
                  rr11 = dmpe(11) - scalek*rr11
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
                  term4 = -2.0d0 * (ck*rr5-dkr*rr7+qkr*rr9)
                  term5 = -2.0d0 * (ci*rr5+dir*rr7+qir*rr9)
                  term6 = 4.0d0 * rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_chgpen) then
                     term1i = corek*rr1i + valk*rr1ik
                     term1k = corei*rr1k + vali*rr1ik
                     term2i = -dkr * rr3ik
                     term2k = dir * rr3ik
                     term3i = qkr * rr5ik
                     term3k = qir * rr5ik
                     poti = term1i + term2i + term3i
                     potk = term1k + term2k + term3k
                  else
                     poti = ck*rr1 - dkr*rr3 + qkr*rr5
                     potk = ci*rr1 + dir*rr3 + qir*rr5
                  end if
                  pot(i) = pot(i) + poti
                  pot(k) = pot(k) + potk
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
               if (use_chgpen)  rr3 = rr3ik
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
c     compute lambda derivative
c
               scalelmda = 1.0d0
               if (muti .and. mutk) then
                  dlambda = 2.0d0 * elambda * e
                  dlambda2 = 2.0d0 * e
                  dlfrcx = 2.0d0 * elambda * frcx
                  dlfrcy = 2.0d0 * elambda * frcy
                  dlfrcz = 2.0d0 * elambda * frcz
                  dlttmi(1) = 2.0d0 * elambda * ttmi(1)
                  dlttmi(2) = 2.0d0 * elambda * ttmi(2)
                  dlttmi(3) = 2.0d0 * elambda * ttmi(3)
                  dlttmk(1) = 2.0d0 * elambda * ttmk(1)
                  dlttmk(2) = 2.0d0 * elambda * ttmk(2)
                  dlttmk(3) = 2.0d0 * elambda * ttmk(3)
                  scalelmda = elambda * elambda
               else if (muti .or. mutk) then
                  dlambda = e
                  dlambda2 = 0.0d0
                  dlfrcx = frcx
                  dlfrcy = frcy
                  dlfrcz = frcz
                  dlttmi(1) = ttmi(1)
                  dlttmi(2) = ttmi(2)
                  dlttmi(3) = ttmi(3)
                  dlttmk(1) = ttmk(1)
                  dlttmk(2) = ttmk(2)
                  dlttmk(3) = ttmk(3)
                  scalelmda = elambda
               end if
               if (muti .or. mutk) then
                  demlmda = demlmda + dlambda
                  demlmda2 = demlmda2 + dlambda2
               end if
c
c     modify the energy, force, and torque by lambda
c
               e = scalelmda * e
               frcx = scalelmda * frcx
               frcy = scalelmda * frcy
               frcz = scalelmda * frcz
               ttmi(1) = scalelmda * ttmi(1)
               ttmi(2) = scalelmda * ttmi(2)
               ttmi(3) = scalelmda * ttmi(3)
               ttmk(1) = scalelmda * ttmk(1)
               ttmk(2) = scalelmda * ttmk(2)
               ttmk(3) = scalelmda * ttmk(3)
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
               if (muti .or. mutk) then
                  dldem(1,i) = dldem(1,i) + dlfrcx
                  dldem(2,i) = dldem(2,i) + dlfrcy
                  dldem(3,i) = dldem(3,i) + dlfrcz
                  dltem(1,i) = dltem(1,i) + dlttmi(1)
                  dltem(2,i) = dltem(2,i) + dlttmi(2)
                  dltem(3,i) = dltem(3,i) + dlttmi(3)
               end if
c
c     increment force-based gradient and torque on second site
c
               dem(1,k) = dem(1,k) - frcx
               dem(2,k) = dem(2,k) - frcy
               dem(3,k) = dem(3,k) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
               if (muti .or. mutk) then
                  dldem(1,k) = dldem(1,k) - dlfrcx
                  dldem(2,k) = dldem(2,k) - dlfrcy
                  dldem(3,k) = dldem(3,k) - dlfrcz
                  dltem(1,k) = dltem(1,k) + dlttmk(1)
                  dltem(2,k) = dltem(2,k) + dlttmk(2)
                  dltem(3,k) = dltem(3,k) + dlttmk(3)
               end if
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
         if (use_chgpen) then
            corei = pcore(i)
            vali = pval(i)
            alphai = palpha(i)
         end if
         muti = mut(i)
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
         do kk = ii, npole
            k = ipole(kk)
            do jcell = 2, ncell
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
            mutk = mut(k)
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               mscale(k) = 1.0d0
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
c     calculate real space Ewald error function damping
c
               call dampewald (11,r,r2,f,dmpe)
c
c     find damped multipole intermediates and energy value
c
               if (use_chgpen) then
                  corek = pcore(k)
                  valk = pval(k)
                  alphak = palpha(k)
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
                  rr1i = dmpe(1) - (1.0d0-scalek*dmpi(1))*rr1
                  rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
                  rr1k = dmpe(1) - (1.0d0-scalek*dmpk(1))*rr1
                  rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
                  rr1ik = dmpe(1) - (1.0d0-scalek*dmpik(1))*rr1
                  rr3ik = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
                  rr5ik = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
                  rr9ik = dmpe(9) - (1.0d0-scalek*dmpik(9))*rr9
                  rr11ik = dmpe(11) - (1.0d0-scalek*dmpik(11))*rr11
                  rr1 = dmpe(1) - (1.0d0-scalek)*rr1
                  rr3 = dmpe(3) - (1.0d0-scalek)*rr3
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
                  rr1 = dmpe(1) - scalek*rr1
                  rr3 = dmpe(3) - scalek*rr3
                  rr5 = dmpe(5) - scalek*rr5
                  rr7 = dmpe(7) - scalek*rr7
                  rr9 = dmpe(9) - scalek*rr9
                  rr11 = dmpe(11) - scalek*rr11
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
                  term4 = -2.0d0 * (ck*rr5-dkr*rr7+qkr*rr9)
                  term5 = -2.0d0 * (ci*rr5+dir*rr7+qir*rr9)
                  term6 = 4.0d0 * rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_chgpen) then
                     term1i = corek*rr1i + valk*rr1ik
                     term1k = corei*rr1k + vali*rr1ik
                     term2i = -dkr * rr3ik
                     term2k = dir * rr3ik
                     term3i = qkr * rr5ik
                     term3k = qir * rr5ik
                     poti = term1i + term2i + term3i
                     potk = term1k + term2k + term3k
                  else
                     poti = ck*rr1 - dkr*rr3 + qkr*rr5
                     potk = ci*rr1 + dir*rr3 + qir*rr5
                  end if
                  pot(i) = pot(i) + poti
                  pot(k) = pot(k) + potk
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
               if (use_chgpen)  rr3 = rr3ik
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
c     compute lambda derivative
c
               scalelmda = 1.0d0
               if (muti .and. mutk) then
                  dlambda = 2.0d0 * elambda * e
                  dlambda2 = 2.0d0 * e
                  dlfrcx = 2.0d0 * elambda * frcx
                  dlfrcy = 2.0d0 * elambda * frcy
                  dlfrcz = 2.0d0 * elambda * frcz
                  dlttmi(1) = 2.0d0 * elambda * ttmi(1)
                  dlttmi(2) = 2.0d0 * elambda * ttmi(2)
                  dlttmi(3) = 2.0d0 * elambda * ttmi(3)
                  dlttmk(1) = 2.0d0 * elambda * ttmk(1)
                  dlttmk(2) = 2.0d0 * elambda * ttmk(2)
                  dlttmk(3) = 2.0d0 * elambda * ttmk(3)
                  scalelmda = elambda * elambda
               else if (muti .or. mutk) then
                  dlambda = e
                  dlambda2 = 0.0d0
                  dlfrcx = frcx
                  dlfrcy = frcy
                  dlfrcz = frcz
                  dlttmi(1) = ttmi(1)
                  dlttmi(2) = ttmi(2)
                  dlttmi(3) = ttmi(3)
                  dlttmk(1) = ttmk(1)
                  dlttmk(2) = ttmk(2)
                  dlttmk(3) = ttmk(3)
                  scalelmda = elambda
               end if
               if (muti .or. mutk) then
                  demlmda = demlmda + dlambda
                  demlmda2 = demlmda2 + dlambda2
               end if
c
c     modify the energy, force, and torque by lambda
c
               e = scalelmda * e
               frcx = scalelmda * frcx
               frcy = scalelmda * frcy
               frcz = scalelmda * frcz
               ttmi(1) = scalelmda * ttmi(1)
               ttmi(2) = scalelmda * ttmi(2)
               ttmi(3) = scalelmda * ttmi(3)
               ttmk(1) = scalelmda * ttmk(1)
               ttmk(2) = scalelmda * ttmk(2)
               ttmk(3) = scalelmda * ttmk(3)
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
               if (muti .or. mutk) then
                  dldem(1,i) = dldem(1,i) + dlfrcx
                  dldem(2,i) = dldem(2,i) + dlfrcy
                  dldem(3,i) = dldem(3,i) + dlfrcz
                  dltem(1,i) = dltem(1,i) + dlttmi(1)
                  dltem(2,i) = dltem(2,i) + dlttmi(2)
                  dltem(3,i) = dltem(3,i) + dlttmi(3)
               end if
c
c     increment force-based gradient and torque on second site
c
               dem(1,k) = dem(1,k) - frcx
               dem(2,k) = dem(2,k) - frcy
               dem(3,k) = dem(3,k) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
               if (muti .or. mutk) then
                  dldem(1,k) = dldem(1,k) - dlfrcx
                  dldem(2,k) = dldem(2,k) - dlfrcy
                  dldem(3,k) = dldem(3,k) - dlfrcz
                  dltem(1,k) = dltem(1,k) + dlttmk(1)
                  dltem(2,k) = dltem(2,k) + dlttmk(2)
                  dltem(3,k) = dltem(3,k) + dlttmk(3)
               end if
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
         call torque (i,dltem(1,i),fix,fiy,fiz,dldem)
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = abs(yaxis(i))
         if (iz .eq. 0)  iz = i
         if (ix .eq. 0)  ix = i
         if (iy .eq. 0)  iy = i
         xiz = x(iz) - x(i)
         yiz = y(iz) - y(i)
         ziz = z(iz) - z(i)
         xix = x(ix) - x(i)
         yix = y(ix) - y(i)
         zix = z(ix) - z(i)
         xiy = x(iy) - x(i)
         yiy = y(iy) - y(i)
         ziy = z(iy) - z(i)
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
c     modify the gradient and virial for charge flux
c
      if (use_chgflx) then
         call dcflux (pot,decfx,decfy,decfz)
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            frcx = decfx(i)
            frcy = decfy(i)
            frcz = decfz(i)
            dem(1,i) = dem(1,i) + frcx
            dem(2,i) = dem(2,i) + frcy
            dem(3,i) = dem(3,i) + frcz
            vxx = xi * frcx
            vxy = yi * frcx
            vxz = zi * frcx
            vyy = yi * frcy
            vyz = zi * frcy
            vzz = zi * frcz
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
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      deallocate (dltem)
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole4d  --  Ewald lambda derivs via list  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole4d" calculates the multipole energy, derivatives with
c     respect to Cartesian coordinates, and lambda derivatives using
c     particle mesh Ewald summation and a neighbor list
c
c
      subroutine empole4d
      use atoms
      use boxes
      use chgpot
      use deriv
      use dlmda
      use energi
      use ewald
      use math
      use mpole
      use mutant
      use pme
      use potent
      use virial
      implicit none
      integer i,j,ii
      real*8 e,f,sum,dlsum
      real*8 dlxd,dlyd,dlzd
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xi,yi,zi
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,ydfield
      real*8 zdfield
      real*8 dlxdfield,dlydfield
      real*8 dlzdfield
      real*8 fx,fy,fz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 tem(3),frcx(3)
      real*8 frcy(3),frcz(3)
      real*8 dltem(3)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
      logical muti
c
c
c     zero out the atomic multipole energy and derivatives
c
      em = 0.0d0
      demlmda = 0.0d0
      demlmda2 = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dldem(j,i) = 0.0d0
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
      call rotpole ('MPORG')
c
c     compute the real space part of the Ewald summation
c
      call emreal4d
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip4
c
c     perform dynamic allocation of some local arrays
c
      allocate (pot(n))
      allocate (decfx(n))
      allocate (decfy(n))
      allocate (decfz(n))
c
c     initialize Ewald self-energy potential array
c
      do i = 1, n
         pot(i) = 0.0d0
      end do
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / rootpi
      do ii = 1, npole
         i = ipole(ii)
         muti = mut(i)
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
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         if (muti) then
            demlmda = demlmda + 2.0d0 * elambda * e
            demlmda2 = demlmda2 + 2.0d0 * e
            e = e * elambda * elambda
         end if
         em = em + e
         pot(i) = 2.0d0 * fterm * ci
      end do
c
c     modify gradient and virial for charge flux self-energy
c
      if (use_chgflx) then
         call dcflux (pot,decfx,decfy,decfz)
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            fx = decfx(i)
            fy = decfy(i)
            fz = decfz(i)
            dem(1,i) = dem(1,i) + fx
            dem(2,i) = dem(2,i) + fy
            dem(3,i) = dem(3,i) + fz
            vxx = xi * fx
            vxy = yi * fx
            vxz = zi * fx
            vyy = yi * fy
            vyz = zi * fy
            vzz = zi * fz
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
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
c
c     compute the uniform background charge correction term
c
      fterm = -0.5d0 * f * pi / (volbox*aewald**2)
      sum = 0.0d0
      dlsum = 0.0d0
      do ii = 1, npole
         i = ipole(ii)
         ci = rpole(1,i)
         muti = mut(i)
         if (muti) then
            dlsum = dlsum + ci
            ci = ci * elambda
         end if
         sum = sum + ci
      end do
      e = fterm * sum**2
      em = em + e
      demlmda = demlmda + fterm * 2.0d0 * sum * dlsum
      demlmda2 = demlmda2 + fterm * 2.0d0 * dlsum**2
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         dlxd = 0.0d0
         dlyd = 0.0d0
         dlzd = 0.0d0
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            muti = mut(i)
            if (muti) then
               dlxd = dlxd + dix + ci*xi
               dlyd = dlyd + diy + ci*yi
               dlzd = dlzd + diz + ci*zi
               ci = ci * elambda
               dix = dix * elambda
               diy = diy * elambda
               diz = diz * elambda
            end if
            xd = xd + dix + ci*xi
            yd = yd + diy + ci*yi
            zd = zd + diz + ci*zi
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         demlmda = demlmda + 2.0d0*term*(xd*dlxd+yd*dlyd+zd*dlzd)
         demlmda2 = demlmda2 + 2.0d0*term*(dlxd**2+dlyd**2+dlzd**2)
         do ii = 1, npole
            i = ipole(ii)
            ci = rpole(1,i)
            muti = mut(i)
            if (muti) then
               dldem(1,i) = dldem(1,i) + 2.0d0*term*ci*(xd+elambda*dlxd)
               dldem(2,i) = dldem(2,i) + 2.0d0*term*ci*(yd+elambda*dlyd)
               dldem(3,i) = dldem(3,i) + 2.0d0*term*ci*(zd+elambda*dlzd)
               ci = ci * elambda
            else
               dldem(1,i) = dldem(1,i) + 2.0d0*term*ci*dlxd
               dldem(2,i) = dldem(2,i) + 2.0d0*term*ci*dlyd
               dldem(3,i) = dldem(3,i) + 2.0d0*term*ci*dlzd
            end if
            dem(1,i) = dem(1,i) + 2.0d0*term*ci*xd
            dem(2,i) = dem(2,i) + 2.0d0*term*ci*yd
            dem(3,i) = dem(3,i) + 2.0d0*term*ci*zd
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         dlxdfield = -2.0d0 * term * dlxd
         dlydfield = -2.0d0 * term * dlyd
         dlzdfield = -2.0d0 * term * dlzd
         do ii = 1, npole
            i = ipole(ii)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            muti = mut(i)
            if (muti) then
               dltem(1) = diy*(zdfield+elambda*dlzdfield)
     &                    - diz*(ydfield+elambda*dlydfield)
               dltem(2) = diz*(xdfield+elambda*dlxdfield)
     &                    - dix*(zdfield+elambda*dlzdfield)
               dltem(3) = dix*(ydfield+elambda*dlydfield)
     &                    - diy*(xdfield+elambda*dlxdfield)
               dix = dix * elambda
               diy = diy * elambda
               diz = diz * elambda
            else
               dltem(1) = diy*dlzdfield - diz*dlydfield
               dltem(2) = diz*dlxdfield - dix*dlzdfield
               dltem(3) = dix*dlydfield - diy*dlxdfield
            end if
            tem(1) = diy*zdfield - diz*ydfield
            tem(2) = diz*xdfield - dix*zdfield
            tem(3) = dix*ydfield - diy*xdfield
            call torque (i,dltem,frcx,frcy,frcz,dldem)
            call torque (i,tem,frcx,frcy,frcz,dem)
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
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            muti = mut(i)
            if (muti) then
               ci = ci * elambda
               dix = dix * elambda
               diy = diy * elambda
               diz = diz * elambda
            end if
            xd = xd + dix
            yd = yd + diy
            zd = zd + diz
            xq = xq + ci*x(i)
            yq = yq + ci*y(i)
            zq = zq + ci*z(i)
         end do
         xv = xd * xq
         yv = yd * yq
         zv = zd * zq
         vterm = term * (xd*xd + yd*yd + zd*zd + 2.0d0*(xv+yv+zv)
     &                      + xq*xq + yq*yq + zq*zq)
         vxx = 2.0d0*term*(xq*xq+xv) + vterm
         vxy = 2.0d0*term*(xq*yq+xv)
         vxz = 2.0d0*term*(xq*zq+xv)
         vyy = 2.0d0*term*(yq*yq+yv) + vterm
         vyz = 2.0d0*term*(yq*zq+yv)
         vzz = 2.0d0*term*(zq*zq+zv) + vterm
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
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emreal4d  --  Ewald real lambda derivs via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emreal4d" evaluates the real space portion of the Ewald
c     summation energy, gradient, and lambda derivatives due to
c     multipole interactions via a neighbor list
c
c
      subroutine emreal4d
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use deriv
      use dlmda
      use energi
      use math
      use mplpot
      use mpole
      use mutant
      use neigh
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      real*8 e,de,f
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
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 poti,potk
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 dlde
      real*8 dlfrcx,dlfrcy,dlfrcz
      real*8 dlambda,dlambda2
      real*8 scalelmda
      real*8 ttmi(3),ttmk(3)
      real*8 dlttmi(3),dlttmk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(11),dmpe(11)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: tem(:,:)
      real*8, allocatable :: dltem(:,:)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
      logical muti,mutk
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (tem(3,n))
      allocate (dltem(3,n))
      allocate (pot(n))
      allocate (decfx(n))
      allocate (decfy(n))
      allocate (decfz(n))
c
c     initialize scaling, torque and potential arrays
c
      do i = 1, n
         mscale(i) = 1.0d0
         do j = 1, 3
            tem(j,i) = 0.0d0
            dltem(j,i) = 0.0d0
         end do
         pot(i) = 0.0d0
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
!$OMP& shared(npole,ipole,x,y,z,rpole,pcore,pval,palpha,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,
!$OMP& nelst,elst,use_chgpen,use_chgflx,use_bounds,f,off2,xaxis,
!$OMP& yaxis,zaxis,elambda,mut)
!$OMP& firstprivate(mscale) shared (em,dem,tem,pot,vir,demlmda,
!$OMP& demlmda2,dldem,dltem)
!$OMP DO reduction(+:em,dem,tem,pot,vir,demlmda,demlmda2,dldem,dltem)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(i)
            vali = pval(i)
            alphai = palpha(i)
         end if
         muti = mut(i)
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
            mutk = mut(k)
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
c     calculate real space Ewald error function damping
c
               call dampewald (11,r,r2,f,dmpe)
c
c     find damped multipole intermediates and energy value
c
               if (use_chgpen) then
                  corek = pcore(k)
                  valk = pval(k)
                  alphak = palpha(k)
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
                  rr1i = dmpe(1) - (1.0d0-scalek*dmpi(1))*rr1
                  rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = dmpe(7) - (1.0d0-scalek*dmpi(7))*rr7
                  rr1k = dmpe(1) - (1.0d0-scalek*dmpk(1))*rr1
                  rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = dmpe(7) - (1.0d0-scalek*dmpk(7))*rr7
                  rr1ik = dmpe(1) - (1.0d0-scalek*dmpik(1))*rr1
                  rr3ik = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
                  rr5ik = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
                  rr9ik = dmpe(9) - (1.0d0-scalek*dmpik(9))*rr9
                  rr11ik = dmpe(11) - (1.0d0-scalek*dmpik(11))*rr11
                  rr1 = dmpe(1) - (1.0d0-scalek)*rr1
                  rr3 = dmpe(3) - (1.0d0-scalek)*rr3
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
                  rr1 = dmpe(1) - scalek*rr1
                  rr3 = dmpe(3) - scalek*rr3
                  rr5 = dmpe(5) - scalek*rr5
                  rr7 = dmpe(7) - scalek*rr7
                  rr9 = dmpe(9) - scalek*rr9
                  rr11 = dmpe(11) - scalek*rr11
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
                  term4 = -2.0d0 * (ck*rr5-dkr*rr7+qkr*rr9)
                  term5 = -2.0d0 * (ci*rr5+dir*rr7+qir*rr9)
                  term6 = 4.0d0 * rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_chgpen) then
                     term1i = corek*rr1i + valk*rr1ik
                     term1k = corei*rr1k + vali*rr1ik
                     term2i = -dkr * rr3ik
                     term2k = dir * rr3ik
                     term3i = qkr * rr5ik
                     term3k = qir * rr5ik
                     poti = term1i + term2i + term3i
                     potk = term1k + term2k + term3k
                  else
                     poti = ck*rr1 - dkr*rr3 + qkr*rr5
                     potk = ci*rr1 + dir*rr3 + qir*rr5
                  end if
                  pot(i) = pot(i) + poti
                  pot(k) = pot(k) + potk
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
               if (use_chgpen)  rr3 = rr3ik
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
c     compute lambda derivative
c
               scalelmda = 1.0d0
               if (muti .and. mutk) then
                  dlambda = 2.0d0 * elambda * e
                  dlambda2 = 2.0d0 * e
                  dlfrcx = 2.0d0 * elambda * frcx
                  dlfrcy = 2.0d0 * elambda * frcy
                  dlfrcz = 2.0d0 * elambda * frcz
                  dlttmi(1) = 2.0d0 * elambda * ttmi(1)
                  dlttmi(2) = 2.0d0 * elambda * ttmi(2)
                  dlttmi(3) = 2.0d0 * elambda * ttmi(3)
                  dlttmk(1) = 2.0d0 * elambda * ttmk(1)
                  dlttmk(2) = 2.0d0 * elambda * ttmk(2)
                  dlttmk(3) = 2.0d0 * elambda * ttmk(3)
                  scalelmda = elambda * elambda
               else if (muti .or. mutk) then
                  dlambda = e
                  dlambda2 = 0.0d0
                  dlfrcx = frcx
                  dlfrcy = frcy
                  dlfrcz = frcz
                  dlttmi(1) = ttmi(1)
                  dlttmi(2) = ttmi(2)
                  dlttmi(3) = ttmi(3)
                  dlttmk(1) = ttmk(1)
                  dlttmk(2) = ttmk(2)
                  dlttmk(3) = ttmk(3)
                  scalelmda = elambda
               end if
               if (muti .or. mutk) then
                  demlmda = demlmda + dlambda
                  demlmda2 = demlmda2 + dlambda2
               end if
c
c     modify the energy, force, and torque by lambda
c
               e = scalelmda * e
               frcx = scalelmda * frcx
               frcy = scalelmda * frcy
               frcz = scalelmda * frcz
               ttmi(1) = scalelmda * ttmi(1)
               ttmi(2) = scalelmda * ttmi(2)
               ttmi(3) = scalelmda * ttmi(3)
               ttmk(1) = scalelmda * ttmk(1)
               ttmk(2) = scalelmda * ttmk(2)
               ttmk(3) = scalelmda * ttmk(3)
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
               if (muti .or. mutk) then
                  dldem(1,i) = dldem(1,i) + dlfrcx
                  dldem(2,i) = dldem(2,i) + dlfrcy
                  dldem(3,i) = dldem(3,i) + dlfrcz
                  dltem(1,i) = dltem(1,i) + dlttmi(1)
                  dltem(2,i) = dltem(2,i) + dlttmi(2)
                  dltem(3,i) = dltem(3,i) + dlttmi(3)
               end if
c
c     increment force-based gradient and torque on second site
c
               dem(1,k) = dem(1,k) - frcx
               dem(2,k) = dem(2,k) - frcy
               dem(3,k) = dem(3,k) - frcz
               tem(1,k) = tem(1,k) + ttmk(1)
               tem(2,k) = tem(2,k) + ttmk(2)
               tem(3,k) = tem(3,k) + ttmk(3)
               if (muti .or. mutk) then
                  dldem(1,k) = dldem(1,k) - dlfrcx
                  dldem(2,k) = dldem(2,k) - dlfrcy
                  dldem(3,k) = dldem(3,k) - dlfrcz
                  dltem(1,k) = dltem(1,k) + dlttmk(1)
                  dltem(2,k) = dltem(2,k) + dlttmk(2)
                  dltem(3,k) = dltem(3,k) + dlttmk(3)
               end if
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
!$OMP DO reduction(+:dem,vir,dldem)
c
c     resolve site torques then increment forces and virial
c
      do ii = 1, npole
         i = ipole(ii)
         call torque (i,dltem(1,i),fix,fiy,fiz,dldem)
         call torque (i,tem(1,i),fix,fiy,fiz,dem)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = abs(yaxis(i))
         if (iz .eq. 0)  iz = i
         if (ix .eq. 0)  ix = i
         if (iy .eq. 0)  iy = i
         xiz = x(iz) - x(i)
         yiz = y(iz) - y(i)
         ziz = z(iz) - z(i)
         xix = x(ix) - x(i)
         yix = y(ix) - y(i)
         zix = z(ix) - z(i)
         xiy = x(iy) - x(i)
         yiy = y(iy) - y(i)
         ziy = z(iy) - z(i)
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
c
c     modify the gradient and virial for charge flux
c
      if (use_chgflx) then
         call dcflux (pot,decfx,decfy,decfz)
!$OMP    DO reduction(+:dem,vir)
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            frcx = decfx(i)
            frcy = decfy(i)
            frcz = decfz(i)
            dem(1,i) = dem(1,i) + frcx
            dem(2,i) = dem(2,i) + frcy
            dem(3,i) = dem(3,i) + frcz
            vxx = xi * frcx
            vxy = yi * frcx
            vxz = zi * frcx
            vyy = yi * frcy
            vyz = zi * frcy
            vzz = zi * frcz
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
!$OMP    END DO
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (tem)
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine emrecip4  --  PME recip mpole lambda derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "emrecip4" evaluates the reciprocal space portion of particle
c     mesh Ewald summation energy, gradient, and lambda derivatives
c     due to multipoles
c
c     literature references:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     W. Smith and D. Fincham, "The Ewald Sum in Truncated Octahedral
c     and Rhombic Dodecahedral Boundary Conditions", Molecular
c     Simulation, 10, 67-71 (1993)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine emrecip4
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use dlmda
      use energi
      use ewald
      use math
      use mpole
      use mrecip
      use mutant
      use pme
      use potent
      use virial
      implicit none
      integer i,j,k,ii
      integer k1,k2,k3
      integer m1,m2,m3
      integer ix,iy,iz
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm,f
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 dlh1,dlh2,dlh3
      real*8 dlf1,dlf2,dlf3
      real*8 xi,yi,zi
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 frcx,frcy,frcz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 dlambda,dlambda2
      real*8 tem(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8 dltem(3)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
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
      if (allocated(cmp)) then
         if (size(cmp) .lt. 10*n)  deallocate (cmp)
      end if
      if (allocated(fmp)) then
         if (size(fmp) .lt. 10*n)  deallocate (fmp)
      end if
      if (allocated(lcmp)) then
         if (size(lcmp) .lt. 10*n)  deallocate (lcmp)
      end if
      if (allocated(lfmp)) then
         if (size(lfmp) .lt. 10*n)  deallocate (lfmp)
      end if
      if (allocated(cphi)) then
         if (size(cphi) .lt. 10*n)  deallocate (cphi)
      end if
      if (allocated(lcphi)) then
         if (size(lcphi) .lt. 10*n)  deallocate (lcphi)
      end if
      if (allocated(fphi)) then
         if (size(fphi) .lt. 20*n)  deallocate (fphi)
      end if
      if (allocated(lfphi)) then
         if (size(lfphi) .lt. 20*n)  deallocate (lfphi)
      end if
      if (.not. allocated(cmp))  allocate (cmp(10,n))
      if (.not. allocated(fmp))  allocate (fmp(10,n))
      if (.not. allocated(lcmp))  allocate (lcmp(10,n))
      if (.not. allocated(lfmp))  allocate (lfmp(10,n))
      if (.not. allocated(cphi))  allocate (cphi(10,n))
      if (.not. allocated(lcphi))  allocate (lcphi(10,n))
      if (.not. allocated(fphi))  allocate (fphi(20,n))
      if (.not. allocated(lfphi))  allocate (lfphi(20,n))
c
c     perform dynamic allocation of some global arrays
c
      ntot = nfft1 * nfft2 * nfft3
      if (allocated(qgrid)) then
         if (size(qgrid) .ne. 2*ntot)  call fftclose
      end if
      if (allocated(lqgrid)) then
         if (size(lqgrid) .ne. 2*ntot)  deallocate (lqgrid)
      end if
      if (.not. allocated(qgrid))  call fftsetup
      if (.not. allocated(lqgrid)) then
         allocate (lqgrid(2,nfft1,nfft2,nfft3))
      end if
c
c     setup spatial decomposition and B-spline coefficients
c
      call getchunk
      call moduli
      call bspline_fill
      call table_fill
c
c     copy multipole moments and coordinates to local storage
c
      do ii = 1, npole
         i = ipole(ii)
         if (mut(i)) then
            cmp(1,i)  = elambda * rpole(1,i)
            cmp(2,i)  = elambda * rpole(2,i)
            cmp(3,i)  = elambda * rpole(3,i)
            cmp(4,i)  = elambda * rpole(4,i)
            cmp(5,i)  = elambda * rpole(5,i)
            cmp(6,i)  = elambda * rpole(9,i)
            cmp(7,i)  = elambda * rpole(13,i)
            cmp(8,i)  = 2.0d0 * elambda * rpole(6,i)
            cmp(9,i)  = 2.0d0 * elambda * rpole(7,i)
            cmp(10,i) = 2.0d0 * elambda * rpole(10,i)
            lcmp(1,i)  = rpole(1,i)
            lcmp(2,i)  = rpole(2,i)
            lcmp(3,i)  = rpole(3,i)
            lcmp(4,i)  = rpole(4,i)
            lcmp(5,i)  = rpole(5,i)
            lcmp(6,i)  = rpole(9,i)
            lcmp(7,i)  = rpole(13,i)
            lcmp(8,i)  = 2.0d0 * rpole(6,i)
            lcmp(9,i)  = 2.0d0 * rpole(7,i)
            lcmp(10,i) = 2.0d0 * rpole(10,i)
         else
            cmp(1,i)  = rpole(1,i)
            cmp(2,i)  = rpole(2,i)
            cmp(3,i)  = rpole(3,i)
            cmp(4,i)  = rpole(4,i)
            cmp(5,i)  = rpole(5,i)
            cmp(6,i)  = rpole(9,i)
            cmp(7,i)  = rpole(13,i)
            cmp(8,i)  = 2.0d0 * rpole(6,i)
            cmp(9,i)  = 2.0d0 * rpole(7,i)
            cmp(10,i) = 2.0d0 * rpole(10,i)
            lcmp(1,i)  = 0.0d0
            lcmp(2,i)  = 0.0d0
            lcmp(3,i)  = 0.0d0
            lcmp(4,i)  = 0.0d0
            lcmp(5,i)  = 0.0d0
            lcmp(6,i)  = 0.0d0
            lcmp(7,i)  = 0.0d0
            lcmp(8,i)  = 0.0d0
            lcmp(9,i)  = 0.0d0
            lcmp(10,i) = 0.0d0
         end if
      end do
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp (cmp,fmp)
      call cmp_to_fmp (lcmp,lfmp)
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_mpole (fmp,qgrid)
      call fftfront (qgrid)
      call grid_mpole (lfmp,lqgrid)
      call fftfront (lqgrid)
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
c     make the scalar summation over reciprocal lattice
c
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
            else if (nonprism) then
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
         qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
         qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
         lqgrid(1,k1,k2,k3) = expterm * lqgrid(1,k1,k2,k3)
         lqgrid(2,k1,k2,k3) = expterm * lqgrid(2,k1,k2,k3)
      end do
c
c     save the partial virial for the polarization computation
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
      qgrid(1,1,1,1) = 0.0d0
      qgrid(2,1,1,1) = 0.0d0
      lqgrid(1,1,1,1) = 0.0d0
      lqgrid(2,1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qgrid(1,1,1,1) = expterm * qgrid(1,1,1,1)
         qgrid(2,1,1,1) = expterm * qgrid(2,1,1,1)
         lqgrid(1,1,1,1) = expterm * lqgrid(1,1,1,1)
         lqgrid(2,1,1,1) = expterm * lqgrid(2,1,1,1)
      end if
c
c     perform 3-D FFT backward transform and get potential
c
      call fftback (qgrid)
      call fphi_mpole (fphi,qgrid)
      call fftback (lqgrid)
      call fphi_mpole (lfphi,lqgrid)
      do ii = 1, npole
         i = ipole(ii)
         do j = 1, 20
            fphi(j,i) = f * fphi(j,i)
            lfphi(j,i) = f * lfphi(j,i)
         end do
      end do
      call fphi_to_cphi (fphi,cphi)
      call fphi_to_cphi (lfphi,lcphi)
c
c     increment the permanent multipole energy and gradient
c
      e = 0.0d0
      dlambda = 0.0d0
      dlambda2 = 0.0d0
      do ii = 1, npole
         i = ipole(ii)
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         dlf1 = 0.0d0
         dlf2 = 0.0d0
         dlf3 = 0.0d0
         do k = 1, 10
            e = e + fmp(k,i)*fphi(k,i)
            dlambda = dlambda + lfmp(k,i) * fphi(k,i)
            dlambda2 = dlambda2 + lfmp(k,i) * lfphi(k,i)
            f1 = f1 + fmp(k,i)*fphi(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphi(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphi(deriv3(k),i)
            dlf1 = dlf1 + lfmp(k,i)*fphi(deriv1(k),i)
     &             + fmp(k,i)*lfphi(deriv1(k),i)
            dlf2 = dlf2 + lfmp(k,i)*fphi(deriv2(k),i)
     &             + fmp(k,i)*lfphi(deriv2(k),i)
            dlf3 = dlf3 + lfmp(k,i)*fphi(deriv3(k),i)
     &             + fmp(k,i)*lfphi(deriv3(k),i)
         end do
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         dlf1 = dble(nfft1) * dlf1
         dlf2 = dble(nfft2) * dlf2
         dlf3 = dble(nfft3) * dlf3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         dlh1 = recip(1,1)*dlf1 + recip(1,2)*dlf2 + recip(1,3)*dlf3
         dlh2 = recip(2,1)*dlf1 + recip(2,2)*dlf2 + recip(2,3)*dlf3
         dlh3 = recip(3,1)*dlf1 + recip(3,2)*dlf2 + recip(3,3)*dlf3
         dem(1,i) = dem(1,i) + h1
         dem(2,i) = dem(2,i) + h2
         dem(3,i) = dem(3,i) + h3
         dldem(1,i) = dldem(1,i) + dlh1
         dldem(2,i) = dldem(2,i) + dlh2
         dldem(3,i) = dldem(3,i) + dlh3
      end do
      e = 0.5d0 * e
      em = em + e
      demlmda = demlmda + dlambda
      demlmda2 = demlmda2 + dlambda2
c
c     increment the permanent multipole virial contributions
c
      do ii = 1, npole
         i = ipole(ii)
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
      do ii = 1, npole
         i = ipole(ii)
         dltem(1) = lcmp(4,i)*cphi(3,i) + cmp(4,i)*lcphi(3,i)
     &            - lcmp(3,i)*cphi(4,i) - cmp(3,i)*lcphi(4,i)
     &            + 2.0d0*(lcmp(7,i)-lcmp(6,i))*cphi(10,i)
     &            + 2.0d0*(cmp(7,i)-cmp(6,i))*lcphi(10,i)
     &            + lcmp(9,i)*cphi(8,i) + cmp(9,i)*lcphi(8,i)
     &            + lcmp(10,i)*cphi(6,i) + cmp(10,i)*lcphi(6,i)
     &            - lcmp(8,i)*cphi(9,i) - cmp(8,i)*lcphi(9,i)
     &            - lcmp(10,i)*cphi(7,i) - cmp(10,i)*lcphi(7,i)

         dltem(2) = lcmp(2,i)*cphi(4,i) + cmp(2,i)*lcphi(4,i)
     &            - lcmp(4,i)*cphi(2,i) - cmp(4,i)*lcphi(2,i)
     &            + 2.0d0*(lcmp(5,i)-lcmp(7,i))*cphi(9,i)
     &            + 2.0d0*(cmp(5,i)-cmp(7,i))*lcphi(9,i)
     &            + lcmp(8,i)*cphi(10,i) + cmp(8,i)*lcphi(10,i)
     &            + lcmp(9,i)*cphi(7,i) + cmp(9,i)*lcphi(7,i)
     &            - lcmp(9,i)*cphi(5,i) - cmp(9,i)*lcphi(5,i)
     &            - lcmp(10,i)*cphi(8,i) - cmp(10,i)*lcphi(8,i)

         dltem(3) = lcmp(3,i)*cphi(2,i) + cmp(3,i)*lcphi(2,i)
     &            - lcmp(2,i)*cphi(3,i) - cmp(2,i)*lcphi(3,i)
     &            + 2.0d0*(lcmp(6,i)-lcmp(5,i))*cphi(8,i)
     &            + 2.0d0*(cmp(6,i)-cmp(5,i))*lcphi(8,i)
     &            + lcmp(8,i)*cphi(5,i) + cmp(8,i)*lcphi(5,i)
     &            + lcmp(10,i)*cphi(9,i) + cmp(10,i)*lcphi(9,i)
     &            - lcmp(8,i)*cphi(6,i) - cmp(8,i)*lcphi(6,i)
     &            - lcmp(9,i)*cphi(10,i) - cmp(9,i)*lcphi(10,i)
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
         call torque (i,dltem,fix,fiy,fiz,dldem)
         call torque (i,tem,fix,fiy,fiz,dem)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = abs(yaxis(i))
         if (iz .eq. 0)  iz = ii
         if (ix .eq. 0)  ix = ii
         if (iy .eq. 0)  iy = ii
         xiz = x(iz) - x(i)
         yiz = y(iz) - y(i)
         ziz = z(iz) - z(i)
         xix = x(ix) - x(i)
         yix = y(ix) - y(i)
         zix = z(ix) - z(i)
         xiy = x(iy) - x(i)
         yiy = y(iy) - y(i)
         ziy = z(iy) - z(i)
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
c     perform dynamic allocation of some local arrays
c
      if (use_chgflx) then
         allocate (pot(n))
         allocate (decfx(n))
         allocate (decfy(n))
         allocate (decfz(n))
c
c     modify the gradient and virial for charge flux
c
         do i = 1, n
            pot(i) = 0.0d0
         end do
         do ii = 1, npole
            i = ipole(ii)
            pot(i) = cphi(1,i)
         end do
         call dcflux (pot,decfx,decfy,decfz)
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            frcx = decfx(i)
            frcy = decfy(i)
            frcz = decfz(i)
            dem(1,i) = dem(1,i) + frcx
            dem(2,i) = dem(2,i) + frcy
            dem(3,i) = dem(3,i) + frcz
            vxx = vxx + xi*frcx
            vxy = vxy + yi*frcx
            vxz = vxz + zi*frcx
            vyy = vyy + yi*frcy
            vyz = vyz + zi*frcy
            vzz = vzz + zi*frcz
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (pot)
         deallocate (decfx)
         deallocate (decfy)
         deallocate (decfz)
      end if
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

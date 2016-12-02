c
c
c     ##################################################
c     ##  COPYRIGHT (C) 2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved             ##
c     ##################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine epolar1  --  polarization energy & derivs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar1" calculates the induced dipole polarization energy
c     and derivatives with respect to Cartesian coordinates
c
c
      subroutine epolar1
      use limits
      implicit none
c
c
c     choose the method for summing over polarization interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call epolar1d
         else
            call epolar1c
         end if
      else
         if (use_mlist) then
            call epolar1b
         else
            call epolar1a
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar1a  --  double loop polarization derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar1a" calculates the dipole polarization energy and
c     derivatives with respect to Cartesian coordinates using a
c     pairwise double loop
c
c
      subroutine epolar1a
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use deriv
      use energi
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer iax,iay,iaz
      real*8 e,f,fi
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 temp3,temp5,temp7
      real*8 sc3,sc5,sc7
      real*8 sr3,sr5,sr7
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dri,drk
      real*8 qri1,qri2,qri3
      real*8 qrk1,qrk2,qrk3
      real*8 qrri,qrrk
      real*8 txxi,tyyi,tzzi
      real*8 txyi,txzi,tyzi
      real*8 txxk,tyyk,tzzk
      real*8 txyk,txzk,tyzk
      real*8 turi,duri,puri
      real*8 turk,durk,purk
      real*8 txi3,tyi3,tzi3
      real*8 txi5,tyi5,tzi5
      real*8 txk3,tyk3,tzk3
      real*8 txk5,tyk5,tzk5
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 term7,term8
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      character*6 mode
c
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0d0
      do i = 1, n
         do j = 1, 3
            dep(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (ufld(3,n))
      allocate (dufld(6,n))
c
c     set arrays needed to scale interactions and store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
         do j = 1, 3
            ufld(j,i) = 0.0d0
         end do
         do j = 1, 6
            dufld(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     get polarization energy via induced dipoles times field
c
      do i = 1, npole
         if (polarity(i) .ne. 0.0d0) then
            fi = -f / polarity(i)
            e = 0.0d0
            do j = 1, 3
               e = fi * uind(j,i) * udirp(j,i)
               ep = ep + e
            end do
         end if
      end do
c
c     compute the dipole polarization gradient components
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
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
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
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
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     sc3 = 1.0d0 - expdamp
                     sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                     sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                    *expdamp
                     temp3 = -damp * expdamp * rr5
                     temp5 = -3.0d0 * damp / r2
                     temp7 = -(1.0d0+3.0d0*damp) / r2
                     rc3(1) = xr * temp3
                     rc3(2) = yr * temp3
                     rc3(3) = zr * temp3
                     rc5(1) = rc3(1) * temp5
                     rc5(2) = rc3(2) * temp5
                     rc5(3) = rc3(3) * temp5
                     rc7(1) = rc5(1) * temp7
                     rc7(2) = rc5(2) * temp7
                     rc7(3) = rc5(3) * temp7
                  end if
               end if
c
c     intermediates involving Thole damping and scale factors
c
               sr3 = rr3 * sc3
               sr5 = rr5 * sc5
               sr7 = rr7 * sc7
               dsr3 = sr3 * dscale(kk)
               dsr5 = sr5 * dscale(kk)
               dsr7 = sr7 * dscale(kk)
               psr3 = sr3 * pscale(kk)
               psr5 = sr5 * pscale(kk)
               psr7 = sr7 * pscale(kk)
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               qri1 = qixx*xr + qixy*yr + qixz*zr
               qri2 = qixy*xr + qiyy*yr + qiyz*zr
               qri3 = qixz*xr + qiyz*yr + qizz*zr
               qrk1 = qkxx*xr + qkxy*yr + qkxz*zr
               qrk2 = qkxy*xr + qkyy*yr + qkyz*zr
               qrk3 = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = (qixx*xr+2.0d0*qixy*yr)*xr
     &                   + (qiyy*yr+2.0d0*qiyz*zr)*yr
     &                   + (qizz*zr+2.0d0*qixz*xr)*zr
               qrrk = (qkxx*xr+2.0d0*qkxy*yr)*xr
     &                   + (qkyy*yr+2.0d0*qkyz*zr)*yr
     &                   + (qkzz*zr+2.0d0*qkxz*xr)*zr
               duri = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               durk = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
               puri = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               purk = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
c
c     get the dE/dR field gradient for direct polarization force
c
               term1 = sc3*(rr3-rr5*xr*xr) + rc3(1)*xr
               term2 = (sc3+sc5)*rr5*xr - rc3(1)
               term3 = sc5*(rr7*xr*xr-rr5) - rc5(1)*xr
               term4 = 2.0d0 * sc5 * rr5
               term5 = 2.0d0 * (sc5*rr7*xr-rc5(1)+1.5d0*sc7*rr7*xr)
               term6 = xr * (sc7*rr9*xr-rc7(1))
               txxi = ci*term1 + dix*term2 - dri*term3
     &                   - qixx*term4 + qri1*term5 - qrri*term6
     &                   + (qri2*yr+qri3*zr)*sc7*rr7
               txxk = ck*term1 - dkx*term2 + drk*term3
     &                   - qkxx*term4 + qrk1*term5 - qrrk*term6
     &                   + (qrk2*yr+qrk3*zr)*sc7*rr7
               term1 = sc3*(rr3-rr5*yr*yr) + rc3(2)*yr
               term2 = (sc3+sc5)*rr5*yr - rc3(2)
               term3 = sc5*(rr7*yr*yr-rr5) - rc5(2)*yr
               term4 = 2.0d0 * sc5 * rr5
               term5 = 2.0d0 * (sc5*rr7*yr-rc5(2)+1.5d0*sc7*rr7*yr)
               term6 = yr * (sc7*rr9*yr-rc7(2))
               tyyi = ci*term1 + diy*term2 - dri*term3
     &                   - qiyy*term4 + qri2*term5 - qrri*term6
     &                   + (qri1*xr+qri3*zr)*sc7*rr7
               tyyk = ck*term1 - dky*term2 + drk*term3
     &                   - qkyy*term4 + qrk2*term5 - qrrk*term6
     &                   + (qrk1*xr+qrk3*zr)*sc7*rr7
               term1 = sc3*(rr3-rr5*zr*zr) + rc3(3)*zr
               term2 = (sc3+sc5)*rr5*zr - rc3(3)
               term3 = sc5*(rr7*zr*zr-rr5) - rc5(3)*zr
               term4 = 2.0d0 * sc5 * rr5
               term5 = 2.0d0 * (sc5*rr7*zr-rc5(3)+1.5d0*sc7*rr7*zr)
               term6 = zr * (sc7*rr9*zr-rc7(3))
               tzzi = ci*term1 + diz*term2 - dri*term3
     &                   - qizz*term4 + qri3*term5 - qrri*term6
     &                   + (qri1*xr+qri2*yr)*sc7*rr7
               tzzk = ck*term1 - dkz*term2 + drk*term3
     &                   - qkzz*term4 + qrk3*term5 - qrrk*term6
     &                   + (qrk1*xr+qrk2*yr)*sc7*rr7
               term2 = sc3*rr5*xr - rc3(1)
               term1 = yr * term2
               term3 = sc5 * rr5 * yr
               term4 = yr * (sc5*rr7*xr-rc5(1))
               term5 = 2.0d0 * sc5 * rr5
               term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
               term7 = 2.0d0 * sc7 * rr7 * yr
               term8 = yr * (sc7*rr9*xr-rc7(1))
               txyi = -ci*term1 + diy*term2 + dix*term3
     &                   - dri*term4 - qixy*term5 + qri2*term6
     &                   + qri1*term7 - qrri*term8
               txyk = -ck*term1 - dky*term2 - dkx*term3
     &                   + drk*term4 - qkxy*term5 + qrk2*term6
     &                   + qrk1*term7 - qrrk*term8
               term2 = sc3*rr5*xr - rc3(1)
               term1 = zr * term2
               term3 = sc5 * rr5 * zr
               term4 = zr * (sc5*rr7*xr-rc5(1))
               term5 = 2.0d0 * sc5 * rr5
               term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
               term7 = 2.0d0 * sc7 * rr7 * zr
               term8 = zr * (sc7*rr9*xr-rc7(1))
               txzi = -ci*term1 + diz*term2 + dix*term3
     &                   - dri*term4 - qixz*term5 + qri3*term6
     &                   + qri1*term7 - qrri*term8
               txzk = -ck*term1 - dkz*term2 - dkx*term3
     &                   + drk*term4 - qkxz*term5 + qrk3*term6
     &                   + qrk1*term7 - qrrk*term8
               term2 = sc3*rr5*yr - rc3(2)
               term1 = zr * term2
               term3 = sc5 * rr5 * zr
               term4 = zr * (sc5*rr7*yr-rc5(2))
               term5 = 2.0d0 * sc5 * rr5
               term6 = 2.0d0 * (sc5*rr7*yr-rc5(2))
               term7 = 2.0d0 * sc7 * rr7 * zr
               term8 = zr * (sc7*rr9*yr-rc7(2))
               tyzi = -ci*term1 + diz*term2 + diy*term3
     &                   - dri*term4 - qiyz*term5 + qri3*term6
     &                   + qri2*term7 - qrri*term8
               tyzk = -ck*term1 - dkz*term2 - dky*term3
     &                + drk*term4 - qkyz*term5 + qrk3*term6
     &                + qrk2*term7 - qrrk*term8
c
c     get the dEd/dR field gradient for direct polarization force
c
               depx = txxi*uinp(1,k) - txxk*uinp(1,i)
     &                   + txyi*uinp(2,k) - txyk*uinp(2,i)
     &                   + txzi*uinp(3,k) - txzk*uinp(3,i)
               depy = txyi*uinp(1,k) - txyk*uinp(1,i)
     &                   + tyyi*uinp(2,k) - tyyk*uinp(2,i)
     &                   + tyzi*uinp(3,k) - tyzk*uinp(3,i)
               depz = txzi*uinp(1,k) - txzk*uinp(1,i)
     &                   + tyzi*uinp(2,k) - tyzk*uinp(2,i)
     &                   + tzzi*uinp(3,k) - tzzk*uinp(3,i)
               frcx = dscale(kk) * depx
               frcy = dscale(kk) * depy
               frcz = dscale(kk) * depz
c
c     get the dEp/dR field gradient for direct polarization force
c
               depx = txxi*uind(1,k) - txxk*uind(1,i)
     &                   + txyi*uind(2,k) - txyk*uind(2,i)
     &                   + txzi*uind(3,k) - txzk*uind(3,i)
               depy = txyi*uind(1,k) - txyk*uind(1,i)
     &                   + tyyi*uind(2,k) - tyyk*uind(2,i)
     &                   + tyzi*uind(3,k) - tyzk*uind(3,i)
               depz = txzi*uind(1,k) - txzk*uind(1,i)
     &                   + tyzi*uind(2,k) - tyzk*uind(2,i)
     &                   + tzzi*uind(3,k) - tzzk*uind(3,i)
               frcx = frcx + pscale(kk)*depx
               frcy = frcy + pscale(kk)*depy
               frcz = frcz + pscale(kk)*depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp .eq. 'MUTUAL') then
                  term1 = (sc3+sc5) * rr5
                  term2 = term1*xr - rc3(1)
                  term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  txxi = uind(1,i)*term2 + duri*term3
                  txxk = uind(1,k)*term2 + durk*term3
                  term2 = term1*yr - rc3(2)
                  term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tyyi = uind(2,i)*term2 + duri*term3
                  tyyk = uind(2,k)*term2 + durk*term3
                  term2 = term1*zr - rc3(3)
                  term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tzzi = uind(3,i)*term2 + duri*term3
                  tzzk = uind(3,k)*term2 + durk*term3
                  term1 = sc5 * rr5 * yr
                  term2 = sc3*rr5*xr - rc3(1)
                  term3 = yr * (sc5*rr7*xr-rc5(1))
                  txyi = uind(1,i)*term1 + uind(2,i)*term2 - duri*term3
                  txyk = uind(1,k)*term1 + uind(2,k)*term2 - durk*term3
                  term1 = sc5 * rr5 * zr
                  term3 = zr * (sc5*rr7*xr-rc5(1))
                  txzi = uind(1,i)*term1 + uind(3,i)*term2 - duri*term3
                  txzk = uind(1,k)*term1 + uind(3,k)*term2 - durk*term3
                  term2 = sc3*rr5*yr - rc3(2)
                  term3 = zr * (sc5*rr7*yr-rc5(2))
                  tyzi = uind(2,i)*term1 + uind(3,i)*term2 - duri*term3
                  tyzk = uind(2,k)*term1 + uind(3,k)*term2 - durk*term3
                  depx = txxi*uinp(1,k) + txxk*uinp(1,i)
     &                      + txyi*uinp(2,k) + txyk*uinp(2,i)
     &                      + txzi*uinp(3,k) + txzk*uinp(3,i)
                  depy = txyi*uinp(1,k) + txyk*uinp(1,i)
     &                      + tyyi*uinp(2,k) + tyyk*uinp(2,i)
     &                      + tyzi*uinp(3,k) + tyzk*uinp(3,i)
                  depz = txzi*uinp(1,k) + txzk*uinp(1,i)
     &                      + tyzi*uinp(2,k) + tyzk*uinp(2,i)
     &                      + tzzi*uinp(3,k) + tzzk*uinp(3,i)
                  frcx = frcx + uscale(kk)*depx
                  frcy = frcy + uscale(kk)*depy
                  frcz = frcz + uscale(kk)*depz
               end if
c
c     increment gradient and virial due to Cartesian forces
c
               dep(1,ii) = dep(1,ii) + frcx
               dep(2,ii) = dep(2,ii) + frcy
               dep(3,ii) = dep(3,ii) + frcz
               dep(1,kk) = dep(1,kk) - frcx
               dep(2,kk) = dep(2,kk) - frcy
               dep(3,kk) = dep(3,kk) - frcz
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
c
c     get the induced dipole field used for dipole torques
c
               txi3 = psr3*uind(1,k) + dsr3*uinp(1,k)
               tyi3 = psr3*uind(2,k) + dsr3*uinp(2,k)
               tzi3 = psr3*uind(3,k) + dsr3*uinp(3,k)
               txk3 = psr3*uind(1,i) + dsr3*uinp(1,i)
               tyk3 = psr3*uind(2,i) + dsr3*uinp(2,i)
               tzk3 = psr3*uind(3,i) + dsr3*uinp(3,i)
               turi = -psr5*durk - dsr5*purk
               turk = -psr5*duri - dsr5*puri
               ufld(1,i) = ufld(1,i) + txi3 + xr*turi
               ufld(2,i) = ufld(2,i) + tyi3 + yr*turi
               ufld(3,i) = ufld(3,i) + tzi3 + zr*turi
               ufld(1,k) = ufld(1,k) + txk3 + xr*turk
               ufld(2,k) = ufld(2,k) + tyk3 + yr*turk
               ufld(3,k) = ufld(3,k) + tzk3 + zr*turk
c
c     get induced dipole field gradient used for quadrupole torques
c
               txi5 = 2.0d0 * (psr5*uind(1,k)+dsr5*uinp(1,k))
               tyi5 = 2.0d0 * (psr5*uind(2,k)+dsr5*uinp(2,k))
               tzi5 = 2.0d0 * (psr5*uind(3,k)+dsr5*uinp(3,k))
               txk5 = 2.0d0 * (psr5*uind(1,i)+dsr5*uinp(1,i))
               tyk5 = 2.0d0 * (psr5*uind(2,i)+dsr5*uinp(2,i))
               tzk5 = 2.0d0 * (psr5*uind(3,i)+dsr5*uinp(3,i))
               turi = -psr7*durk - dsr7*purk
               turk = -psr7*duri - dsr7*puri
               dufld(1,i) = dufld(1,i) + xr*txi5 + xr*xr*turi
               dufld(2,i) = dufld(2,i) + xr*tyi5 + yr*txi5
     &                         + 2.0d0*xr*yr*turi
               dufld(3,i) = dufld(3,i) + yr*tyi5 + yr*yr*turi
               dufld(4,i) = dufld(4,i) + xr*tzi5 + zr*txi5
     &                         + 2.0d0*xr*zr*turi
               dufld(5,i) = dufld(5,i) + yr*tzi5 + zr*tyi5
     &                         + 2.0d0*yr*zr*turi
               dufld(6,i) = dufld(6,i) + zr*tzi5 + zr*zr*turi
               dufld(1,k) = dufld(1,k) - xr*txk5 - xr*xr*turk
               dufld(2,k) = dufld(2,k) - xr*tyk5 - yr*txk5
     &                         - 2.0d0*xr*yr*turk
               dufld(3,k) = dufld(3,k) - yr*tyk5 - yr*yr*turk
               dufld(4,k) = dufld(4,k) - xr*tzk5 - zr*txk5
     &                         - 2.0d0*xr*zr*turk
               dufld(5,k) = dufld(5,k) - yr*tzk5 - zr*tyk5
     &                         - 2.0d0*yr*zr*turk
               dufld(6,k) = dufld(6,k) - zr*tzk5 - zr*zr*turk
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
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
         pdi = pdamp(i)
         pti = thole(i)
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
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
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
            if (use_bounds)  call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               pscale(kk) = 1.0d0
               dscale(kk) = 1.0d0
               uscale(kk) = 1.0d0
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
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     sc3 = 1.0d0 - expdamp
                     sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                     sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                    *expdamp
                     temp3 = -damp * expdamp * rr5
                     temp5 = -3.0d0 * damp / r2
                     temp7 = -(1.0d0+3.0d0*damp) / r2
                     rc3(1) = xr * temp3
                     rc3(2) = yr * temp3
                     rc3(3) = zr * temp3
                     rc5(1) = rc3(1) * temp5
                     rc5(2) = rc3(2) * temp5
                     rc5(3) = rc3(3) * temp5
                     rc7(1) = rc5(1) * temp7
                     rc7(2) = rc5(2) * temp7
                     rc7(3) = rc5(3) * temp7
                  end if
               end if
c
c     intermediates involving Thole damping and scale factors
c
               sr3 = rr3 * sc3
               sr5 = rr5 * sc5
               sr7 = rr7 * sc7
               psr3 = sr3 * pscale(kk)
               psr5 = sr5 * pscale(kk)
               psr7 = sr7 * pscale(kk)
               dsr3 = sr3 * dscale(kk)
               dsr5 = sr5 * dscale(kk)
               dsr7 = sr7 * dscale(kk)
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               qri1 = qixx*xr + qixy*yr + qixz*zr
               qri2 = qixy*xr + qiyy*yr + qiyz*zr
               qri3 = qixz*xr + qiyz*yr + qizz*zr
               qrk1 = qkxx*xr + qkxy*yr + qkxz*zr
               qrk2 = qkxy*xr + qkyy*yr + qkyz*zr
               qrk3 = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = (qixx*xr+2.0d0*qixy*yr)*xr
     &                   + (qiyy*yr+2.0d0*qiyz*zr)*yr
     &                   + (qizz*zr+2.0d0*qixz*xr)*zr
               qrrk = (qkxx*xr+2.0d0*qkxy*yr)*xr
     &                   + (qkyy*yr+2.0d0*qkyz*zr)*yr
     &                   + (qkzz*zr+2.0d0*qkxz*xr)*zr
               duri = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               durk = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
               puri = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               purk = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
c
c     get the dE/dR field gradient for direct polarization force
c
               term1 = sc3*(rr3-rr5*xr*xr) + rc3(1)*xr
               term2 = (sc3+sc5)*rr5*xr - rc3(1)
               term3 = sc5*(rr7*xr*xr-rr5) - rc5(1)*xr
               term4 = 2.0d0 * sc5 * rr5
               term5 = 2.0d0 * (sc5*rr7*xr-rc5(1)+1.5d0*sc7*rr7*xr)
               term6 = xr * (sc7*rr9*xr-rc7(1))
               txxi = ci*term1 + dix*term2 - dri*term3
     &                   - qixx*term4 + qri1*term5 - qrri*term6
     &                   + (qri2*yr+qri3*zr)*sc7*rr7
               txxk = ck*term1 - dkx*term2 + drk*term3
     &                   - qkxx*term4 + qrk1*term5 - qrrk*term6
     &                   + (qrk2*yr+qrk3*zr)*sc7*rr7
               term1 = sc3*(rr3-rr5*yr*yr) + rc3(2)*yr
               term2 = (sc3+sc5)*rr5*yr - rc3(2)
               term3 = sc5*(rr7*yr*yr-rr5) - rc5(2)*yr
               term4 = 2.0d0 * sc5 * rr5
               term5 = 2.0d0 * (sc5*rr7*yr-rc5(2)+1.5d0*sc7*rr7*yr)
               term6 = yr * (sc7*rr9*yr-rc7(2))
               tyyi = ci*term1 + diy*term2 - dri*term3
     &                   - qiyy*term4 + qri2*term5 - qrri*term6
     &                   + (qri1*xr+qri3*zr)*sc7*rr7
               tyyk = ck*term1 - dky*term2 + drk*term3
     &                   - qkyy*term4 + qrk2*term5 - qrrk*term6
     &                   + (qrk1*xr+qrk3*zr)*sc7*rr7
               term1 = sc3*(rr3-rr5*zr*zr) + rc3(3)*zr
               term2 = (sc3+sc5)*rr5*zr - rc3(3)
               term3 = sc5*(rr7*zr*zr-rr5) - rc5(3)*zr
               term4 = 2.0d0 * sc5 * rr5
               term5 = 2.0d0 * (sc5*rr7*zr-rc5(3)+1.5d0*sc7*rr7*zr)
               term6 = zr * (sc7*rr9*zr-rc7(3))
               tzzi = ci*term1 + diz*term2 - dri*term3
     &                   - qizz*term4 + qri3*term5 - qrri*term6
     &                   + (qri1*xr+qri2*yr)*sc7*rr7
               tzzk = ck*term1 - dkz*term2 + drk*term3
     &                   - qkzz*term4 + qrk3*term5 - qrrk*term6
     &                   + (qrk1*xr+qrk2*yr)*sc7*rr7
               term2 = sc3*rr5*xr - rc3(1)
               term1 = yr * term2
               term3 = sc5 * rr5 * yr
               term4 = yr * (sc5*rr7*xr-rc5(1))
               term5 = 2.0d0 * sc5 * rr5
               term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
               term7 = 2.0d0 * sc7 * rr7 * yr
               term8 = yr * (sc7*rr9*xr-rc7(1))
               txyi = -ci*term1 + diy*term2 + dix*term3
     &                   - dri*term4 - qixy*term5 + qri2*term6
     &                   + qri1*term7 - qrri*term8
               txyk = -ck*term1 - dky*term2 - dkx*term3
     &                   + drk*term4 - qkxy*term5 + qrk2*term6
     &                   + qrk1*term7 - qrrk*term8
               term2 = sc3*rr5*xr - rc3(1)
               term1 = zr * term2
               term3 = sc5 * rr5 * zr
               term4 = zr * (sc5*rr7*xr-rc5(1))
               term5 = 2.0d0 * sc5 * rr5
               term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
               term7 = 2.0d0 * sc7 * rr7 * zr
               term8 = zr * (sc7*rr9*xr-rc7(1))
               txzi = -ci*term1 + diz*term2 + dix*term3
     &                   - dri*term4 - qixz*term5 + qri3*term6
     &                   + qri1*term7 - qrri*term8
               txzk = -ck*term1 - dkz*term2 - dkx*term3
     &                   + drk*term4 - qkxz*term5 + qrk3*term6
     &                   + qrk1*term7 - qrrk*term8
               term2 = sc3*rr5*yr - rc3(2)
               term1 = zr * term2
               term3 = sc5 * rr5 * zr
               term4 = zr * (sc5*rr7*yr-rc5(2))
               term5 = 2.0d0 * sc5 * rr5
               term6 = 2.0d0 * (sc5*rr7*yr-rc5(2))
               term7 = 2.0d0 * sc7 * rr7 * zr
               term8 = zr * (sc7*rr9*yr-rc7(2))
               tyzi = -ci*term1 + diz*term2 + diy*term3
     &                   - dri*term4 - qiyz*term5 + qri3*term6
     &                   + qri2*term7 - qrri*term8
               tyzk = -ck*term1 - dkz*term2 - dky*term3
     &                + drk*term4 - qkyz*term5 + qrk3*term6
     &                + qrk2*term7 - qrrk*term8
c
c     get the dEd/dR field gradient for direct polarization force
c
               depx = txxi*uinp(1,k) - txxk*uinp(1,i)
     &                   + txyi*uinp(2,k) - txyk*uinp(2,i)
     &                   + txzi*uinp(3,k) - txzk*uinp(3,i)
               depy = txyi*uinp(1,k) - txyk*uinp(1,i)
     &                   + tyyi*uinp(2,k) - tyyk*uinp(2,i)
     &                   + tyzi*uinp(3,k) - tyzk*uinp(3,i)
               depz = txzi*uinp(1,k) - txzk*uinp(1,i)
     &                   + tyzi*uinp(2,k) - tyzk*uinp(2,i)
     &                   + tzzi*uinp(3,k) - tzzk*uinp(3,i)
               frcx = dscale(kk) * depx
               frcy = dscale(kk) * depy
               frcz = dscale(kk) * depz
c
c     get the dEp/dR field gradient for direct polarization force
c
               depx = txxi*uind(1,k) - txxk*uind(1,i)
     &                   + txyi*uind(2,k) - txyk*uind(2,i)
     &                   + txzi*uind(3,k) - txzk*uind(3,i)
               depy = txyi*uind(1,k) - txyk*uind(1,i)
     &                   + tyyi*uind(2,k) - tyyk*uind(2,i)
     &                   + tyzi*uind(3,k) - tyzk*uind(3,i)
               depz = txzi*uind(1,k) - txzk*uind(1,i)
     &                   + tyzi*uind(2,k) - tyzk*uind(2,i)
     &                   + tzzi*uind(3,k) - tzzk*uind(3,i)
               frcx = frcx + pscale(kk)*depx
               frcy = frcy + pscale(kk)*depy
               frcz = frcz + pscale(kk)*depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp .eq. 'MUTUAL') then
                  term1 = (sc3+sc5) * rr5
                  term2 = term1*xr - rc3(1)
                  term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  txxi = uind(1,i)*term2 + duri*term3
                  txxk = uind(1,k)*term2 + durk*term3
                  term2 = term1*yr - rc3(2)
                  term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tyyi = uind(2,i)*term2 + duri*term3
                  tyyk = uind(2,k)*term2 + durk*term3
                  term2 = term1*zr - rc3(3)
                  term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tzzi = uind(3,i)*term2 + duri*term3
                  tzzk = uind(3,k)*term2 + durk*term3
                  term1 = sc5 * rr5 * yr
                  term2 = sc3*rr5*xr - rc3(1)
                  term3 = yr * (sc5*rr7*xr-rc5(1))
                  txyi = uind(1,i)*term1 + uind(2,i)*term2 - duri*term3
                  txyk = uind(1,k)*term1 + uind(2,k)*term2 - durk*term3
                  term1 = sc5 * rr5 * zr
                  term3 = zr * (sc5*rr7*xr-rc5(1))
                  txzi = uind(1,i)*term1 + uind(3,i)*term2 - duri*term3
                  txzk = uind(1,k)*term1 + uind(3,k)*term2 - durk*term3
                  term2 = sc3*rr5*yr - rc3(2)
                  term3 = zr * (sc5*rr7*yr-rc5(2))
                  tyzi = uind(2,i)*term1 + uind(3,i)*term2 - duri*term3
                  tyzk = uind(2,k)*term1 + uind(3,k)*term2 - durk*term3
                  depx = txxi*uinp(1,k) + txxk*uinp(1,i)
     &                      + txyi*uinp(2,k) + txyk*uinp(2,i)
     &                      + txzi*uinp(3,k) + txzk*uinp(3,i)
                  depy = txyi*uinp(1,k) + txyk*uinp(1,i)
     &                      + tyyi*uinp(2,k) + tyyk*uinp(2,i)
     &                      + tyzi*uinp(3,k) + tyzk*uinp(3,i)
                  depz = txzi*uinp(1,k) + txzk*uinp(1,i)
     &                      + tyzi*uinp(2,k) + tyzk*uinp(2,i)
     &                      + tzzi*uinp(3,k) + tzzk*uinp(3,i)
                  frcx = frcx + uscale(kk)*depx
                  frcy = frcy + uscale(kk)*depy
                  frcz = frcz + uscale(kk)*depz
               end if
c
c     force and torque components scaled for self-interactions
c
               if (ii .eq. kk) then
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  psr3 = 0.5d0 * psr3
                  psr5 = 0.5d0 * psr5
                  psr7 = 0.5d0 * psr7
                  dsr3 = 0.5d0 * dsr3
                  dsr5 = 0.5d0 * dsr5
                  dsr7 = 0.5d0 * dsr7
               end if
c
c     increment gradient and virial due to Cartesian forces
c
               dep(1,ii) = dep(1,ii) + frcx
               dep(2,ii) = dep(2,ii) + frcy
               dep(3,ii) = dep(3,ii) + frcz
               dep(1,kk) = dep(1,kk) - frcx
               dep(2,kk) = dep(2,kk) - frcy
               dep(3,kk) = dep(3,kk) - frcz
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
c
c     get the induced dipole field used for dipole torques
c
               txi3 = psr3*uind(1,k) + dsr3*uinp(1,k)
               tyi3 = psr3*uind(2,k) + dsr3*uinp(2,k)
               tzi3 = psr3*uind(3,k) + dsr3*uinp(3,k)
               txk3 = psr3*uind(1,i) + dsr3*uinp(1,i)
               tyk3 = psr3*uind(2,i) + dsr3*uinp(2,i)
               tzk3 = psr3*uind(3,i) + dsr3*uinp(3,i)
               turi = -psr5*durk - dsr5*purk
               turk = -psr5*duri - dsr5*puri
               ufld(1,i) = ufld(1,i) + txi3 + xr*turi
               ufld(2,i) = ufld(2,i) + tyi3 + yr*turi
               ufld(3,i) = ufld(3,i) + tzi3 + zr*turi
               ufld(1,k) = ufld(1,k) + txk3 + xr*turk
               ufld(2,k) = ufld(2,k) + tyk3 + yr*turk
               ufld(3,k) = ufld(3,k) + tzk3 + zr*turk
c
c     get induced dipole field gradient used for quadrupole torques
c
               txi5 = 2.0d0 * (psr5*uind(1,k)+dsr5*uinp(1,k))
               tyi5 = 2.0d0 * (psr5*uind(2,k)+dsr5*uinp(2,k))
               tzi5 = 2.0d0 * (psr5*uind(3,k)+dsr5*uinp(3,k))
               txk5 = 2.0d0 * (psr5*uind(1,i)+dsr5*uinp(1,i))
               tyk5 = 2.0d0 * (psr5*uind(2,i)+dsr5*uinp(2,i))
               tzk5 = 2.0d0 * (psr5*uind(3,i)+dsr5*uinp(3,i))
               turi = -psr7*durk - dsr7*purk
               turk = -psr7*duri - dsr7*puri
               dufld(1,i) = dufld(1,i) + xr*txi5 + xr*xr*turi
               dufld(2,i) = dufld(2,i) + xr*tyi5 + yr*txi5
     &                         + 2.0d0*xr*yr*turi
               dufld(3,i) = dufld(3,i) + yr*tyi5 + yr*yr*turi
               dufld(4,i) = dufld(4,i) + xr*tzi5 + zr*txi5
     &                         + 2.0d0*xr*zr*turi
               dufld(5,i) = dufld(5,i) + yr*tzi5 + zr*tyi5
     &                         + 2.0d0*yr*zr*turi
               dufld(6,i) = dufld(6,i) + zr*tzi5 + zr*zr*turi
               dufld(1,k) = dufld(1,k) - xr*txk5 - xr*xr*turk
               dufld(2,k) = dufld(2,k) - xr*tyk5 - yr*txk5
     &                         - 2.0d0*xr*yr*turk
               dufld(3,k) = dufld(3,k) - yr*tyk5 - yr*yr*turk
               dufld(4,k) = dufld(4,k) - xr*tzk5 - zr*txk5
     &                         - 2.0d0*xr*zr*turk
               dufld(5,k) = dufld(5,k) - yr*tzk5 - zr*tyk5
     &                         - 2.0d0*yr*zr*turk
               dufld(6,k) = dufld(6,k) - zr*tzk5 - zr*zr*turk
            end if
            end do
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
      end if
c
c     torque is induced field and gradient cross permanent moments
c
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         trq(1) = diz*ufld(2,i) - diy*ufld(3,i)
     &               + qixz*dufld(2,i) - qixy*dufld(4,i)
     &               + 2.0d0*qiyz*(dufld(3,i)-dufld(6,i))
     &               + (qizz-qiyy)*dufld(5,i)
         trq(2) = dix*ufld(3,i) - diz*ufld(1,i)
     &               - qiyz*dufld(2,i) + qixy*dufld(5,i)
     &               + 2.0d0*qixz*(dufld(6,i)-dufld(1,i))
     &               + (qixx-qizz)*dufld(4,i)
         trq(3) = diy*ufld(1,i) - dix*ufld(2,i)
     &               + qiyz*dufld(4,i) - qixz*dufld(5,i)
     &               + 2.0d0*qixy*(dufld(1,i)-dufld(3,i))
     &               + (qiyy-qixx)*dufld(2,i)
         call torque (i,trq,fix,fiy,fiz,dep)
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
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (ufld)
      deallocate (dufld)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar1b  --  neighbor list polarization derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar1b" calculates the dipole polarization energy and
c     derivatives with respect to Cartesian coordinates using a
c     neighbor list
c
c
      subroutine epolar1b
      use sizes
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use energi
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer iax,iay,iaz
      real*8 e,epo,f,fi
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 temp3,temp5,temp7
      real*8 sc3,sc5,sc7
      real*8 sr3,sr5,sr7
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dri,drk
      real*8 qri1,qri2,qri3
      real*8 qrk1,qrk2,qrk3
      real*8 qrri,qrrk
      real*8 txxi,tyyi,tzzi
      real*8 txyi,txzi,tyzi
      real*8 txxk,tyyk,tzzk
      real*8 txyk,txzk,tyzk
      real*8 turi,duri,puri
      real*8 turk,durk,purk
      real*8 txi3,tyi3,tzi3
      real*8 txi5,tyi5,tzi5
      real*8 txk3,tyk3,tzk3
      real*8 txk5,tyk5,tzk5
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 term7,term8
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8 viro(3,3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: depo(:,:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      character*6 mode
c
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0d0
      do i = 1, n
         do j = 1, 3
            dep(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (ufld(3,n))
      allocate (dufld(6,n))
      allocate (depo(3,n))
c
c     set arrays needed to scale interactions and store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
         do j = 1, 3
            ufld(j,i) = 0.0d0
         end do
         do j = 1, 6
            dufld(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     initialize local variables for OpenMP calculation
c
      epo = 0.0d0
      do i = 1, n
         do j = 1, 3
            depo(j,i) = 0.0d0
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
!$OMP PARALLEL default(private) shared(npole,polarity,f,uind,udirp,epo)
!$OMP& shared(ipole,pdamp,thole,x,y,z,xaxis,yaxis,zaxis,rpole,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,
!$OMP& p2scale,p3scale,p4scale,p41scale,p5scale,d1scale,d2scale,d3scale,
!$OMP& d4scale,u1scale,u2scale,u3scale,u4scale,nelst,elst,use_bounds,
!$OMP& off2,uinp,poltyp,depo,viro,ufld,dufld)
!$OMP& firstprivate(pscale,dscale,uscale)
!$OMP DO reduction(+:epo) schedule(guided)
c
c     get polarization energy via induced dipoles times field
c
      do i = 1, npole
         if (polarity(i) .ne. 0.0d0) then
            fi = -f / polarity(i)
            e = 0.0d0
            do j = 1, 3
               e = fi * uind(j,i) * udirp(j,i)
               epo = epo + e
            end do
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO reduction(+:depo,viro,ufld,dufld) schedule(guided)
c
c     compute the dipole polarization gradient components
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
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
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
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
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     sc3 = 1.0d0 - expdamp
                     sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                     sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                    *expdamp
                     temp3 = -damp * expdamp * rr5
                     temp5 = -3.0d0 * damp / r2
                     temp7 = -(1.0d0+3.0d0*damp) / r2
                     rc3(1) = xr * temp3
                     rc3(2) = yr * temp3
                     rc3(3) = zr * temp3
                     rc5(1) = rc3(1) * temp5
                     rc5(2) = rc3(2) * temp5
                     rc5(3) = rc3(3) * temp5
                     rc7(1) = rc5(1) * temp7
                     rc7(2) = rc5(2) * temp7
                     rc7(3) = rc5(3) * temp7
                  end if
               end if
c
c     intermediates involving Thole damping and scale factors
c
               sr3 = rr3 * sc3
               sr5 = rr5 * sc5
               sr7 = rr7 * sc7
               psr3 = sr3 * pscale(kk)
               psr5 = sr5 * pscale(kk)
               psr7 = sr7 * pscale(kk)
               dsr3 = sr3 * dscale(kk)
               dsr5 = sr5 * dscale(kk)
               dsr7 = sr7 * dscale(kk)
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               qri1 = qixx*xr + qixy*yr + qixz*zr
               qri2 = qixy*xr + qiyy*yr + qiyz*zr
               qri3 = qixz*xr + qiyz*yr + qizz*zr
               qrk1 = qkxx*xr + qkxy*yr + qkxz*zr
               qrk2 = qkxy*xr + qkyy*yr + qkyz*zr
               qrk3 = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = (qixx*xr+2.0d0*qixy*yr)*xr
     &                   + (qiyy*yr+2.0d0*qiyz*zr)*yr
     &                   + (qizz*zr+2.0d0*qixz*xr)*zr
               qrrk = (qkxx*xr+2.0d0*qkxy*yr)*xr
     &                   + (qkyy*yr+2.0d0*qkyz*zr)*yr
     &                   + (qkzz*zr+2.0d0*qkxz*xr)*zr
               duri = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               durk = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
               puri = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               purk = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
c
c     get the dE/dR field gradient for direct polarization force
c
               term1 = sc3*(rr3-rr5*xr*xr) + rc3(1)*xr
               term2 = (sc3+sc5)*rr5*xr - rc3(1)
               term3 = sc5*(rr7*xr*xr-rr5) - rc5(1)*xr
               term4 = 2.0d0 * sc5 * rr5
               term5 = 2.0d0 * (sc5*rr7*xr-rc5(1)+1.5d0*sc7*rr7*xr)
               term6 = xr * (sc7*rr9*xr-rc7(1))
               txxi = ci*term1 + dix*term2 - dri*term3
     &                   - qixx*term4 + qri1*term5 - qrri*term6
     &                   + (qri2*yr+qri3*zr)*sc7*rr7
               txxk = ck*term1 - dkx*term2 + drk*term3
     &                   - qkxx*term4 + qrk1*term5 - qrrk*term6
     &                   + (qrk2*yr+qrk3*zr)*sc7*rr7
               term1 = sc3*(rr3-rr5*yr*yr) + rc3(2)*yr
               term2 = (sc3+sc5)*rr5*yr - rc3(2)
               term3 = sc5*(rr7*yr*yr-rr5) - rc5(2)*yr
               term4 = 2.0d0 * sc5 * rr5
               term5 = 2.0d0 * (sc5*rr7*yr-rc5(2)+1.5d0*sc7*rr7*yr)
               term6 = yr * (sc7*rr9*yr-rc7(2))
               tyyi = ci*term1 + diy*term2 - dri*term3
     &                   - qiyy*term4 + qri2*term5 - qrri*term6
     &                   + (qri1*xr+qri3*zr)*sc7*rr7
               tyyk = ck*term1 - dky*term2 + drk*term3
     &                   - qkyy*term4 + qrk2*term5 - qrrk*term6
     &                   + (qrk1*xr+qrk3*zr)*sc7*rr7
               term1 = sc3*(rr3-rr5*zr*zr) + rc3(3)*zr
               term2 = (sc3+sc5)*rr5*zr - rc3(3)
               term3 = sc5*(rr7*zr*zr-rr5) - rc5(3)*zr
               term4 = 2.0d0 * sc5 * rr5
               term5 = 2.0d0 * (sc5*rr7*zr-rc5(3)+1.5d0*sc7*rr7*zr)
               term6 = zr * (sc7*rr9*zr-rc7(3))
               tzzi = ci*term1 + diz*term2 - dri*term3
     &                   - qizz*term4 + qri3*term5 - qrri*term6
     &                   + (qri1*xr+qri2*yr)*sc7*rr7
               tzzk = ck*term1 - dkz*term2 + drk*term3
     &                   - qkzz*term4 + qrk3*term5 - qrrk*term6
     &                   + (qrk1*xr+qrk2*yr)*sc7*rr7
               term2 = sc3*rr5*xr - rc3(1)
               term1 = yr * term2
               term3 = sc5 * rr5 * yr
               term4 = yr * (sc5*rr7*xr-rc5(1))
               term5 = 2.0d0 * sc5 * rr5
               term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
               term7 = 2.0d0 * sc7 * rr7 * yr
               term8 = yr * (sc7*rr9*xr-rc7(1))
               txyi = -ci*term1 + diy*term2 + dix*term3
     &                   - dri*term4 - qixy*term5 + qri2*term6
     &                   + qri1*term7 - qrri*term8
               txyk = -ck*term1 - dky*term2 - dkx*term3
     &                   + drk*term4 - qkxy*term5 + qrk2*term6
     &                   + qrk1*term7 - qrrk*term8
               term2 = sc3*rr5*xr - rc3(1)
               term1 = zr * term2
               term3 = sc5 * rr5 * zr
               term4 = zr * (sc5*rr7*xr-rc5(1))
               term5 = 2.0d0 * sc5 * rr5
               term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
               term7 = 2.0d0 * sc7 * rr7 * zr
               term8 = zr * (sc7*rr9*xr-rc7(1))
               txzi = -ci*term1 + diz*term2 + dix*term3
     &                   - dri*term4 - qixz*term5 + qri3*term6
     &                   + qri1*term7 - qrri*term8
               txzk = -ck*term1 - dkz*term2 - dkx*term3
     &                   + drk*term4 - qkxz*term5 + qrk3*term6
     &                   + qrk1*term7 - qrrk*term8
               term2 = sc3*rr5*yr - rc3(2)
               term1 = zr * term2
               term3 = sc5 * rr5 * zr
               term4 = zr * (sc5*rr7*yr-rc5(2))
               term5 = 2.0d0 * sc5 * rr5
               term6 = 2.0d0 * (sc5*rr7*yr-rc5(2))
               term7 = 2.0d0 * sc7 * rr7 * zr
               term8 = zr * (sc7*rr9*yr-rc7(2))
               tyzi = -ci*term1 + diz*term2 + diy*term3
     &                   - dri*term4 - qiyz*term5 + qri3*term6
     &                   + qri2*term7 - qrri*term8
               tyzk = -ck*term1 - dkz*term2 - dky*term3
     &                + drk*term4 - qkyz*term5 + qrk3*term6
     &                + qrk2*term7 - qrrk*term8
c
c     get the dEd/dR field gradient for direct polarization force
c
               depx = txxi*uinp(1,k) - txxk*uinp(1,i)
     &                   + txyi*uinp(2,k) - txyk*uinp(2,i)
     &                   + txzi*uinp(3,k) - txzk*uinp(3,i)
               depy = txyi*uinp(1,k) - txyk*uinp(1,i)
     &                   + tyyi*uinp(2,k) - tyyk*uinp(2,i)
     &                   + tyzi*uinp(3,k) - tyzk*uinp(3,i)
               depz = txzi*uinp(1,k) - txzk*uinp(1,i)
     &                   + tyzi*uinp(2,k) - tyzk*uinp(2,i)
     &                   + tzzi*uinp(3,k) - tzzk*uinp(3,i)
               frcx = dscale(kk) * depx
               frcy = dscale(kk) * depy
               frcz = dscale(kk) * depz
c
c     get the dEp/dR field gradient for direct polarization force
c
               depx = txxi*uind(1,k) - txxk*uind(1,i)
     &                   + txyi*uind(2,k) - txyk*uind(2,i)
     &                   + txzi*uind(3,k) - txzk*uind(3,i)
               depy = txyi*uind(1,k) - txyk*uind(1,i)
     &                   + tyyi*uind(2,k) - tyyk*uind(2,i)
     &                   + tyzi*uind(3,k) - tyzk*uind(3,i)
               depz = txzi*uind(1,k) - txzk*uind(1,i)
     &                   + tyzi*uind(2,k) - tyzk*uind(2,i)
     &                   + tzzi*uind(3,k) - tzzk*uind(3,i)
               frcx = frcx + pscale(kk)*depx
               frcy = frcy + pscale(kk)*depy
               frcz = frcz + pscale(kk)*depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp .eq. 'MUTUAL') then
                  term1 = (sc3+sc5) * rr5
                  term2 = term1*xr - rc3(1)
                  term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  txxi = uind(1,i)*term2 + duri*term3
                  txxk = uind(1,k)*term2 + durk*term3
                  term2 = term1*yr - rc3(2)
                  term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tyyi = uind(2,i)*term2 + duri*term3
                  tyyk = uind(2,k)*term2 + durk*term3
                  term2 = term1*zr - rc3(3)
                  term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tzzi = uind(3,i)*term2 + duri*term3
                  tzzk = uind(3,k)*term2 + durk*term3
                  term1 = sc5 * rr5 * yr
                  term2 = sc3*rr5*xr - rc3(1)
                  term3 = yr * (sc5*rr7*xr-rc5(1))
                  txyi = uind(1,i)*term1 + uind(2,i)*term2 - duri*term3
                  txyk = uind(1,k)*term1 + uind(2,k)*term2 - durk*term3
                  term1 = sc5 * rr5 * zr
                  term3 = zr * (sc5*rr7*xr-rc5(1))
                  txzi = uind(1,i)*term1 + uind(3,i)*term2 - duri*term3
                  txzk = uind(1,k)*term1 + uind(3,k)*term2 - durk*term3
                  term2 = sc3*rr5*yr - rc3(2)
                  term3 = zr * (sc5*rr7*yr-rc5(2))
                  tyzi = uind(2,i)*term1 + uind(3,i)*term2 - duri*term3
                  tyzk = uind(2,k)*term1 + uind(3,k)*term2 - durk*term3
                  depx = txxi*uinp(1,k) + txxk*uinp(1,i)
     &                      + txyi*uinp(2,k) + txyk*uinp(2,i)
     &                      + txzi*uinp(3,k) + txzk*uinp(3,i)
                  depy = txyi*uinp(1,k) + txyk*uinp(1,i)
     &                      + tyyi*uinp(2,k) + tyyk*uinp(2,i)
     &                      + tyzi*uinp(3,k) + tyzk*uinp(3,i)
                  depz = txzi*uinp(1,k) + txzk*uinp(1,i)
     &                      + tyzi*uinp(2,k) + tyzk*uinp(2,i)
     &                      + tzzi*uinp(3,k) + tzzk*uinp(3,i)
                  frcx = frcx + uscale(kk)*depx
                  frcy = frcy + uscale(kk)*depy
                  frcz = frcz + uscale(kk)*depz
               end if
c
c     increment gradient and virial due to Cartesian forces
c
               depo(1,ii) = depo(1,ii) + frcx
               depo(2,ii) = depo(2,ii) + frcy
               depo(3,ii) = depo(3,ii) + frcz
               depo(1,kk) = depo(1,kk) - frcx
               depo(2,kk) = depo(2,kk) - frcy
               depo(3,kk) = depo(3,kk) - frcz
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
c
c     get the induced dipole field used for dipole torques
c
               txi3 = psr3*uind(1,k) + dsr3*uinp(1,k)
               tyi3 = psr3*uind(2,k) + dsr3*uinp(2,k)
               tzi3 = psr3*uind(3,k) + dsr3*uinp(3,k)
               txk3 = psr3*uind(1,i) + dsr3*uinp(1,i)
               tyk3 = psr3*uind(2,i) + dsr3*uinp(2,i)
               tzk3 = psr3*uind(3,i) + dsr3*uinp(3,i)
               turi = -psr5*durk - dsr5*purk
               turk = -psr5*duri - dsr5*puri
               ufld(1,i) = ufld(1,i) + txi3 + xr*turi
               ufld(2,i) = ufld(2,i) + tyi3 + yr*turi
               ufld(3,i) = ufld(3,i) + tzi3 + zr*turi
               ufld(1,k) = ufld(1,k) + txk3 + xr*turk
               ufld(2,k) = ufld(2,k) + tyk3 + yr*turk
               ufld(3,k) = ufld(3,k) + tzk3 + zr*turk
c
c     get induced dipole field gradient used for quadrupole torques
c
               txi5 = 2.0d0 * (psr5*uind(1,k)+dsr5*uinp(1,k))
               tyi5 = 2.0d0 * (psr5*uind(2,k)+dsr5*uinp(2,k))
               tzi5 = 2.0d0 * (psr5*uind(3,k)+dsr5*uinp(3,k))
               txk5 = 2.0d0 * (psr5*uind(1,i)+dsr5*uinp(1,i))
               tyk5 = 2.0d0 * (psr5*uind(2,i)+dsr5*uinp(2,i))
               tzk5 = 2.0d0 * (psr5*uind(3,i)+dsr5*uinp(3,i))
               turi = -psr7*duri - dsr7*purk
               turk = -psr7*duri - dsr7*puri
               dufld(1,i) = dufld(1,i) + xr*txi5 + xr*xr*turi
               dufld(2,i) = dufld(2,i) + xr*tyi5 + yr*txi5
     &                         + 2.0d0*xr*yr*turi
               dufld(3,i) = dufld(3,i) + yr*tyi5 + yr*yr*turi
               dufld(4,i) = dufld(4,i) + xr*tzi5 + zr*txi5
     &                         + 2.0d0*xr*zr*turi
               dufld(5,i) = dufld(5,i) + yr*tzi5 + zr*tyi5
     &                         + 2.0d0*yr*zr*turi
               dufld(6,i) = dufld(6,i) + zr*tzi5 + zr*zr*turi
               dufld(1,k) = dufld(1,k) - xr*txk5 - xr*xr*turk
               dufld(2,k) = dufld(2,k) - xr*tyk5 - yr*txk5
     &                         - 2.0d0*xr*yr*turk
               dufld(3,k) = dufld(3,k) - yr*tyk5 - yr*yr*turk
               dufld(4,k) = dufld(4,k) - xr*tzk5 - zr*txk5
     &                         - 2.0d0*xr*zr*turk
               dufld(5,k) = dufld(5,k) - yr*tzk5 - zr*tyk5
     &                         - 2.0d0*yr*zr*turk
               dufld(6,k) = dufld(6,k) - zr*tzk5 - zr*zr*turk
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO reduction(+:depo,viro) schedule(guided)
c
c     torque is induced field and gradient cross permanent moments
c
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         trq(1) = diz*ufld(2,i) - diy*ufld(3,i)
     &               + qixz*dufld(2,i) - qixy*dufld(4,i)
     &               + 2.0d0*qiyz*(dufld(3,i)-dufld(6,i))
     &               + (qizz-qiyy)*dufld(5,i)
         trq(2) = dix*ufld(3,i) - diz*ufld(1,i)
     &               - qiyz*dufld(2,i) + qixy*dufld(5,i)
     &               + 2.0d0*qixz*(dufld(6,i)-dufld(1,i))
     &               + (qixx-qizz)*dufld(4,i)
         trq(3) = diy*ufld(1,i) - dix*ufld(2,i)
     &               + qiyz*dufld(4,i) - qixz*dufld(5,i)
     &               + 2.0d0*qixy*(dufld(1,i)-dufld(3,i))
     &               + (qiyy-qixx)*dufld(2,i)
         call torque (i,trq,fix,fiy,fiz,depo)
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
c     add local copies to global variables for OpenMP calculation
c
      ep = ep + epo
      do i = 1, n
         do j = 1, 3
            dep(j,i) = dep(j,i) + depo(j,i)
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
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (depo)
      deallocate (ufld)
      deallocate (dufld)
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar1c  --  Ewald polarization derivs via loop  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar1c" calculates the dipole polarization energy and
c     derivatives with respect to Cartesian coordinates using
c     particle mesh Ewald summation and a double loop
c
c
      subroutine epolar1c
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use polar
      use polpot
      use potent
      use virial
      implicit none
      integer i,j,ii
      real*8 e,f,term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xv,yv,zv,vterm
      real*8 xufield,yufield
      real*8 zufield
      real*8 fix(3),fiy(3),fiz(3)
      real*8 trq(3)
c
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0d0
      do i = 1, n
         do j = 1, 3
            dep(j,i) = 0.0d0
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
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     compute the real space part of the Ewald summation
c
      call epreal1c
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip1
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * term * uii / 3.0d0
c        ep = ep + e
      end do
c
c     compute the self-energy torque term due to induced dipole
c
      term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = 0.5d0 * (uind(1,i)+uinp(1,i))
         uiy = 0.5d0 * (uind(2,i)+uinp(2,i))
         uiz = 0.5d0 * (uind(3,i)+uinp(3,i))
         trq(1) = term * (diy*uiz-diz*uiy)
         trq(2) = term * (diz*uix-dix*uiz)
         trq(3) = term * (dix*uiy-diy*uix)
         call torque (i,trq,fix,fiy,fiz,dep)
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         xup = 0.0d0
         yup = 0.0d0
         zup = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
            xu = xu + uind(1,i)
            yu = yu + uind(2,i)
            zu = zu + uind(3,i)
            xup = xup + uinp(1,i)
            yup = yup + uinp(2,i)
            zup = zup + uinp(3,i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
c        ep = ep + term*(xd*xu+yd*yu+zd*zu)
         do i = 1, npole
            ii = ipole(i)
            dep(1,ii) = dep(1,ii) + term*rpole(1,i)*(xu+xup)
            dep(2,ii) = dep(2,ii) + term*rpole(1,i)*(yu+yup)
            dep(3,ii) = dep(3,ii) + term*rpole(1,i)*(zu+zup)
         end do
         xufield = -term * (xu+xup)
         yufield = -term * (yu+yup)
         zufield = -term * (zu+zup)
         do i = 1, npole
            trq(1) = rpole(3,i)*zufield - rpole(4,i)*yufield
            trq(2) = rpole(4,i)*xufield - rpole(2,i)*zufield
            trq(3) = rpole(2,i)*yufield - rpole(3,i)*xufield
            call torque (i,trq,fix,fiy,fiz,dep)
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
         xv = xq * (xu+xup)
         yv = yq * (yu+yup)
         zv = zq * (zu+zup)
         vterm = xv + yv + zv + xu*xup + yu*yup + zu*zup
     &              + xd*(xu+xup) + yd*(yu+yup) + zd*(zu+zup)
         vterm = term * vterm
         vir(1,1) = vir(1,1) + term*xv + vterm
         vir(2,1) = vir(2,1) + term*xv
         vir(3,1) = vir(3,1) + term*xv
         vir(1,2) = vir(1,2) + term*yv
         vir(2,2) = vir(2,2) + term*yv + vterm
         vir(3,2) = vir(3,2) + term*yv
         vir(1,3) = vir(1,3) + term*zv
         vir(2,3) = vir(2,3) + term*zv
         vir(3,3) = vir(3,3) + term*zv + vterm
         if (poltyp .eq. 'DIRECT') then
            vterm = term * (xu*xup+yu*yup+zu*zup)
            vir(1,1) = vir(1,1) + vterm
            vir(2,2) = vir(2,2) + vterm
            vir(3,3) = vir(3,3) + vterm
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal1c  --  Ewald real space derivs via loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal1c" evaluates the real space portion of the Ewald
c     summation energy and gradient due to dipole polarization
c     via a double loop
c
c
      subroutine epreal1c
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
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      integer iax,iay,iaz
      real*8 e,epo,f,fi
      real*8 erfc,bfac
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 temp3,temp5,temp7
      real*8 scalekk
      real*8 sc3,sc5,sc7
      real*8 psc3,psc5,psc7
      real*8 dsc3,dsc5,dsc7
      real*8 usc3,usc5
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 usr3,usr5
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dri,drk
      real*8 qri1,qri2,qri3
      real*8 qrk1,qrk2,qrk3
      real*8 qrri,qrrk
      real*8 txxi,tyyi,tzzi
      real*8 txyi,txzi,tyzi
      real*8 txxk,tyyk,tzzk
      real*8 txyk,txzk,tyzk
      real*8 turi,duri,puri
      real*8 turk,durk,purk
      real*8 txi3,tyi3,tzi3
      real*8 txi5,tyi5,tzi5
      real*8 txk3,tyk3,tzk3
      real*8 txk5,tyk5,tzk5
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term6,term7
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 prc3(3),prc5(3),prc7(3)
      real*8 drc3(3),drc5(3),drc7(3)
      real*8 urc3(3),urc5(3)
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8 bn(0:4)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (ufld(3,n))
      allocate (dufld(6,n))
c
c     set arrays needed to scale interactions and store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
         do j = 1, 3
            ufld(j,i) = 0.0d0
         end do
         do j = 1, 6
            dufld(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     get polarization energy via induced dipoles times field
c
      do i = 1, npole
         if (polarity(i) .ne. 0.0d0) then
            fi = -f / polarity(i)
            e = 0.0d0
            do j = 1, 3
               e = fi * uind(j,i) * udirp(j,i)
               ep = ep + e
            end do
         end if
      end do
c
c     compute the dipole polarization gradient components
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
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
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
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
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     sc3 = 1.0d0 - expdamp
                     sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                     sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                    *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     rc3(1) = xr * temp3
                     rc3(2) = yr * temp3
                     rc3(3) = zr * temp3
                     rc5(1) = rc3(1) * temp5
                     rc5(2) = rc3(2) * temp5
                     rc5(3) = rc3(3) * temp5
                     rc7(1) = rc5(1) * temp7
                     rc7(2) = rc5(2) * temp7
                     rc7(3) = rc5(3) * temp7
                  end if
               end if
c
c     intermediates involving Thole damping and scale factors
c
               psc3 = 1.0d0 - sc3*pscale(kk)
               psc5 = 1.0d0 - sc5*pscale(kk)
               psc7 = 1.0d0 - sc7*pscale(kk)
               dsc3 = 1.0d0 - sc3*dscale(kk)
               dsc5 = 1.0d0 - sc5*dscale(kk)
               dsc7 = 1.0d0 - sc7*dscale(kk)
               usc3 = 1.0d0 - sc3*uscale(kk)
               usc5 = 1.0d0 - sc5*uscale(kk)
               psr3 = bn(1) - psc3*rr3
               psr5 = bn(2) - psc5*rr5
               psr7 = bn(3) - psc7*rr7
               dsr3 = bn(1) - dsc3*rr3
               dsr5 = bn(2) - dsc5*rr5
               dsr7 = bn(3) - dsc7*rr7
               usr3 = bn(1) - usc3*rr3
               usr5 = bn(2) - usc5*rr5
               do j = 1, 3
                  prc3(j) = rc3(j) * pscale(kk)
                  prc5(j) = rc5(j) * pscale(kk)
                  prc7(j) = rc7(j) * pscale(kk)
                  drc3(j) = rc3(j) * dscale(kk)
                  drc5(j) = rc5(j) * dscale(kk)
                  drc7(j) = rc7(j) * dscale(kk)
                  urc3(j) = rc3(j) * uscale(kk)
                  urc5(j) = rc5(j) * uscale(kk)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               qri1 = qixx*xr + qixy*yr + qixz*zr
               qri2 = qixy*xr + qiyy*yr + qiyz*zr
               qri3 = qixz*xr + qiyz*yr + qizz*zr
               qrk1 = qkxx*xr + qkxy*yr + qkxz*zr
               qrk2 = qkxy*xr + qkyy*yr + qkyz*zr
               qrk3 = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = (qixx*xr+2.0d0*qixy*yr)*xr
     &                   + (qiyy*yr+2.0d0*qiyz*zr)*yr
     &                   + (qizz*zr+2.0d0*qixz*xr)*zr
               qrrk = (qkxx*xr+2.0d0*qkxy*yr)*xr
     &                   + (qkyy*yr+2.0d0*qkyz*zr)*yr
     &                   + (qkzz*zr+2.0d0*qkxz*xr)*zr
               puri = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               purk = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               duri = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               durk = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
c
c     get the dEd/dR field gradient for direct polarization force
c
               term1 = bn(2) - dsc3*rr5
               term2 = bn(3) - dsc5*rr7
               term3 = -dsr3 + term1*xr*xr - rr3*xr*drc3(1)
               term4 = rr3*drc3(1) - term1*xr - dsr5*xr
               term5 = term2*xr*xr - dsr5 - rr5*xr*drc5(1)
               term6 = (bn(4)-dsc7*rr9)*xr*xr - bn(3) - rr7*xr*drc7(1)
               term7 = rr5*drc5(1) - 2.0d0*bn(3)*xr
     &                    + (dsc5+1.5d0*dsc7)*rr7*xr
               txxi = ci*term3 + dix*term4 + dri*term5
     &                   + 2.0d0*dsr5*qixx + (qri2*yr+qri3*zr)*dsc7*rr7
     &                   + 2.0d0*qri1*term7 + qrri*term6
               txxk = ck*term3 - dkx*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkxx + (qrk2*yr+qrk3*zr)*dsc7*rr7
     &                   + 2.0d0*qrk1*term7 + qrrk*term6
               term3 = -dsr3 + term1*yr*yr - rr3*yr*drc3(2)
               term4 = rr3*drc3(2) - term1*yr - dsr5*yr
               term5 = term2*yr*yr - dsr5 - rr5*yr*drc5(2)
               term6 = (bn(4)-dsc7*rr9)*yr*yr - bn(3) - rr7*yr*drc7(2)
               term7 = rr5*drc5(2) - 2.0d0*bn(3)*yr
     &                    + (dsc5+1.5d0*dsc7)*rr7*yr
               tyyi = ci*term3 + diy*term4 + dri*term5
     &                   + 2.0d0*dsr5*qiyy + (qri1*xr+qri3*zr)*dsc7*rr7
     &                   + 2.0d0*qri2*term7 + qrri*term6
               tyyk = ck*term3 - dky*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkyy + (qrk1*xr+qrk3*zr)*dsc7*rr7
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = -dsr3 + term1*zr*zr - rr3*zr*drc3(3)
               term4 = rr3*drc3(3) - term1*zr - dsr5*zr
               term5 = term2*zr*zr - dsr5 - rr5*zr*drc5(3)
               term6 = (bn(4)-dsc7*rr9)*zr*zr - bn(3) - rr7*zr*drc7(3)
               term7 = rr5*drc5(3) - 2.0d0*bn(3)*zr
     &                    + (dsc5+1.5d0*dsc7)*rr7*zr
               tzzi = ci*term3 + diz*term4 + dri*term5
     &                   + 2.0d0*dsr5*qizz + (qri1*xr+qri2*yr)*dsc7*rr7
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tzzk = ck*term3 - dkz*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkzz + (qrk1*xr+qrk2*yr)*dsc7*rr7
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*xr*yr - rr3*yr*drc3(1)
               term4 = rr3*drc3(1) - term1*xr
               term5 = term2*xr*yr - rr5*yr*drc5(1)
               term6 = (bn(4)-dsc7*rr9)*xr*yr - rr7*yr*drc7(1)
               term7 = rr5*drc5(1) - term2*xr
               txyi = ci*term3 - dsr5*dix*yr + diy*term4 + dri*term5
     &                   + 2.0d0*dsr5*qixy - 2.0d0*dsr7*yr*qri1
     &                   + 2.0d0*qri2*term7 + qrri*term6
               txyk = ck*term3 + dsr5*dkx*yr - dky*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkxy - 2.0d0*dsr7*yr*qrk1
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = term1*xr*zr - rr3*zr*drc3(1)
               term5 = term2*xr*zr - rr5*zr*drc5(1)
               term6 = (bn(4)-dsc7*rr9)*xr*zr - rr7*zr*drc7(1)
               txzi = ci*term3 - dsr5*dix*zr + diz*term4 + dri*term5
     &                   + 2.0d0*dsr5*qixz - 2.0d0*dsr7*zr*qri1
     &                   + 2.0d0*qri3*term7 + qrri*term6
               txzk = ck*term3 + dsr5*dkx*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkxz - 2.0d0*dsr7*zr*qrk1
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*yr*zr - rr3*zr*drc3(2)
               term4 = rr3*drc3(2) - term1*yr
               term5 = term2*yr*zr - rr5*zr*drc5(2)
               term6 = (bn(4)-dsc7*rr9)*yr*zr - rr7*zr*drc7(2)
               term7 = rr5*drc5(2) - term2*yr
               tyzi = ci*term3 - dsr5*diy*zr + diz*term4 + dri*term5
     &                   + 2.0d0*dsr5*qiyz - 2.0d0*dsr7*zr*qri2
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tyzk = ck*term3 + dsr5*dky*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkyz - 2.0d0*dsr7*zr*qrk2
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               depx = txxi*uinp(1,k) - txxk*uinp(1,i)
     &                   + txyi*uinp(2,k) - txyk*uinp(2,i)
     &                   + txzi*uinp(3,k) - txzk*uinp(3,i)
               depy = txyi*uinp(1,k) - txyk*uinp(1,i)
     &                   + tyyi*uinp(2,k) - tyyk*uinp(2,i)
     &                   + tyzi*uinp(3,k) - tyzk*uinp(3,i)
               depz = txzi*uinp(1,k) - txzk*uinp(1,i)
     &                   + tyzi*uinp(2,k) - tyzk*uinp(2,i)
     &                   + tzzi*uinp(3,k) - tzzk*uinp(3,i)
               frcx = depx
               frcy = depy
               frcz = depz
c
c     get the dEp/dR field gradient for direct polarization force
c
               term1 = bn(2) - psc3*rr5
               term2 = bn(3) - psc5*rr7
               term3 = -psr3 + term1*xr*xr - rr3*xr*prc3(1)
               term4 = rr3*prc3(1) - term1*xr - psr5*xr
               term5 = term2*xr*xr - psr5 - rr5*xr*prc5(1)
               term6 = (bn(4)-psc7*rr9)*xr*xr - bn(3) - rr7*xr*prc7(1)
               term7 = rr5*prc5(1) - 2.0d0*bn(3)*xr
     &                    + (psc5+1.5d0*psc7)*rr7*xr
               txxi = ci*term3 + dix*term4 + dri*term5
     &                   + 2.0d0*psr5*qixx + (qri2*yr+qri3*zr)*psc7*rr7
     &                   + 2.0d0*qri1*term7 + qrri*term6
               txxk = ck*term3 - dkx*term4 - drk*term5
     &                   + 2.0d0*psr5*qkxx + (qrk2*yr+qrk3*zr)*psc7*rr7
     &                   + 2.0d0*qrk1*term7 + qrrk*term6
               term3 = -psr3 + term1*yr*yr - rr3*yr*prc3(2)
               term4 = rr3*prc3(2) - term1*yr - psr5*yr
               term5 = term2*yr*yr - psr5 - rr5*yr*prc5(2)
               term6 = (bn(4)-psc7*rr9)*yr*yr - bn(3) - rr7*yr*prc7(2)
               term7 = rr5*prc5(2) - 2.0d0*bn(3)*yr
     &                    + (psc5+1.5d0*psc7)*rr7*yr
               tyyi = ci*term3 + diy*term4 + dri*term5
     &                   + 2.0d0*psr5*qiyy + (qri1*xr+qri3*zr)*psc7*rr7
     &                   + 2.0d0*qri2*term7 + qrri*term6
               tyyk = ck*term3 - dky*term4 - drk*term5
     &                   + 2.0d0*psr5*qkyy + (qrk1*xr+qrk3*zr)*psc7*rr7
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = -psr3 + term1*zr*zr - rr3*zr*prc3(3)
               term4 = rr3*prc3(3) - term1*zr - psr5*zr
               term5 = term2*zr*zr - psr5 - rr5*zr*prc5(3)
               term6 = (bn(4)-psc7*rr9)*zr*zr - bn(3) - rr7*zr*prc7(3)
               term7 = rr5*prc5(3) - 2.0d0*bn(3)*zr
     &                    + (psc5+1.5d0*psc7)*rr7*zr
               tzzi = ci*term3 + diz*term4 + dri*term5
     &                   + 2.0d0*psr5*qizz + (qri1*xr+qri2*yr)*psc7*rr7
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tzzk = ck*term3 - dkz*term4 - drk*term5
     &                   + 2.0d0*psr5*qkzz + (qrk1*xr+qrk2*yr)*psc7*rr7
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*xr*yr - rr3*yr*prc3(1)
               term4 = rr3*prc3(1) - term1*xr
               term5 = term2*xr*yr - rr5*yr*prc5(1)
               term6 = (bn(4)-psc7*rr9)*xr*yr - rr7*yr*prc7(1)
               term7 = rr5*prc5(1) - term2*xr
               txyi = ci*term3 - psr5*dix*yr + diy*term4 + dri*term5
     &                   + 2.0d0*psr5*qixy - 2.0d0*psr7*yr*qri1
     &                   + 2.0d0*qri2*term7 + qrri*term6
               txyk = ck*term3 + psr5*dkx*yr - dky*term4 - drk*term5
     &                   + 2.0d0*psr5*qkxy - 2.0d0*psr7*yr*qrk1
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = term1*xr*zr - rr3*zr*prc3(1)
               term5 = term2*xr*zr - rr5*zr*prc5(1)
               term6 = (bn(4)-psc7*rr9)*xr*zr - rr7*zr*prc7(1)
               txzi = ci*term3 - psr5*dix*zr + diz*term4 + dri*term5
     &                   + 2.0d0*psr5*qixz - 2.0d0*psr7*zr*qri1
     &                   + 2.0d0*qri3*term7 + qrri*term6
               txzk = ck*term3 + psr5*dkx*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*psr5*qkxz - 2.0d0*psr7*zr*qrk1
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*yr*zr - rr3*zr*prc3(2)
               term4 = rr3*prc3(2) - term1*yr
               term5 = term2*yr*zr - rr5*zr*prc5(2)
               term6 = (bn(4)-psc7*rr9)*yr*zr - rr7*zr*prc7(2)
               term7 = rr5*prc5(2) - term2*yr
               tyzi = ci*term3 - psr5*diy*zr + diz*term4 + dri*term5
     &                   + 2.0d0*psr5*qiyz - 2.0d0*psr7*zr*qri2
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tyzk = ck*term3 + psr5*dky*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*psr5*qkyz - 2.0d0*psr7*zr*qrk2
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               depx = txxi*uind(1,k) - txxk*uind(1,i)
     &                   + txyi*uind(2,k) - txyk*uind(2,i)
     &                   + txzi*uind(3,k) - txzk*uind(3,i)
               depy = txyi*uind(1,k) - txyk*uind(1,i)
     &                   + tyyi*uind(2,k) - tyyk*uind(2,i)
     &                   + tyzi*uind(3,k) - tyzk*uind(3,i)
               depz = txzi*uind(1,k) - txzk*uind(1,i)
     &                   + tyzi*uind(2,k) - tyzk*uind(2,i)
     &                   + tzzi*uind(3,k) - tzzk*uind(3,i)
               frcx = frcx + depx
               frcy = frcy + depy
               frcz = frcz + depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp .eq. 'MUTUAL') then
                  term1 = bn(2) - usc3*rr5
                  term2 = bn(3) - usc5*rr7
                  term3 = usr5 + term1
                  term4 = rr3 * uscale(kk)
                  term5 = -xr*term3 + rc3(1)*term4
                  term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                  txxi = uind(1,i)*term5 + duri*term6
                  txxk = uind(1,k)*term5 + durk*term6
                  term5 = -yr*term3 + rc3(2)*term4
                  term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                  tyyi = uind(2,i)*term5 + duri*term6
                  tyyk = uind(2,k)*term5 + durk*term6
                  term5 = -zr*term3 + rc3(3)*term4
                  term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                  tzzi = uind(3,i)*term5 + duri*term6
                  tzzk = uind(3,k)*term5 + durk*term6
                  term4 = -usr5 * yr
                  term5 = -xr*term1 + rr3*urc3(1)
                  term6 = xr*yr*term2 - rr5*yr*urc5(1)
                  txyi = uind(1,i)*term4 + uind(2,i)*term5 + duri*term6
                  txyk = uind(1,k)*term4 + uind(2,k)*term5 + durk*term6
                  term4 = -usr5 * zr
                  term6 = xr*zr*term2 - rr5*zr*urc5(1)
                  txzi = uind(1,i)*term4 + uind(3,i)*term5 + duri*term6
                  txzk = uind(1,k)*term4 + uind(3,k)*term5 + durk*term6
                  term5 = -yr*term1 + rr3*urc3(2)
                  term6 = yr*zr*term2 - rr5*zr*urc5(2)
                  tyzi = uind(2,i)*term4 + uind(3,i)*term5 + duri*term6
                  tyzk = uind(2,k)*term4 + uind(3,k)*term5 + durk*term6
                  depx = txxi*uinp(1,k) + txxk*uinp(1,i)
     &                      + txyi*uinp(2,k) + txyk*uinp(2,i)
     &                      + txzi*uinp(3,k) + txzk*uinp(3,i)
                  depy = txyi*uinp(1,k) + txyk*uinp(1,i)
     &                      + tyyi*uinp(2,k) + tyyk*uinp(2,i)
     &                      + tyzi*uinp(3,k) + tyzk*uinp(3,i)
                  depz = txzi*uinp(1,k) + txzk*uinp(1,i)
     &                      + tyzi*uinp(2,k) + tyzk*uinp(2,i)
     &                      + tzzi*uinp(3,k) + tzzk*uinp(3,i)
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
               end if
c
c     increment gradient and virial due to Cartesian forces
c
               dep(1,ii) = dep(1,ii) - frcx
               dep(2,ii) = dep(2,ii) - frcy
               dep(3,ii) = dep(3,ii) - frcz
               dep(1,kk) = dep(1,kk) + frcx
               dep(2,kk) = dep(2,kk) + frcy
               dep(3,kk) = dep(3,kk) + frcz
               vxx = xr * frcx
               vxy = yr * frcx
               vxz = zr * frcx
               vyy = yr * frcy
               vyz = zr * frcy
               vzz = zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
c
c     get the induced dipole field used for dipole torques
c
               txi3 = psr3*uind(1,k) + dsr3*uinp(1,k)
               tyi3 = psr3*uind(2,k) + dsr3*uinp(2,k)
               tzi3 = psr3*uind(3,k) + dsr3*uinp(3,k)
               txk3 = psr3*uind(1,i) + dsr3*uinp(1,i)
               tyk3 = psr3*uind(2,i) + dsr3*uinp(2,i)
               tzk3 = psr3*uind(3,i) + dsr3*uinp(3,i)
               turi = -psr5*durk - dsr5*purk
               turk = -psr5*duri - dsr5*puri
               ufld(1,i) = ufld(1,i) + txi3 + xr*turi
               ufld(2,i) = ufld(2,i) + tyi3 + yr*turi
               ufld(3,i) = ufld(3,i) + tzi3 + zr*turi
               ufld(1,k) = ufld(1,k) + txk3 + xr*turk
               ufld(2,k) = ufld(2,k) + tyk3 + yr*turk
               ufld(3,k) = ufld(3,k) + tzk3 + zr*turk
c
c     get induced dipole field gradient used for quadrupole torques
c
               txi5 = 2.0d0 * (psr5*uind(1,k)+dsr5*uinp(1,k))
               tyi5 = 2.0d0 * (psr5*uind(2,k)+dsr5*uinp(2,k))
               tzi5 = 2.0d0 * (psr5*uind(3,k)+dsr5*uinp(3,k))
               txk5 = 2.0d0 * (psr5*uind(1,i)+dsr5*uinp(1,i))
               tyk5 = 2.0d0 * (psr5*uind(2,i)+dsr5*uinp(2,i))
               tzk5 = 2.0d0 * (psr5*uind(3,i)+dsr5*uinp(3,i))
               turi = -psr7*durk - dsr7*purk 
               turk = -psr7*duri - dsr7*puri 
               dufld(1,i) = dufld(1,i) + xr*txi5 + xr*xr*turi
               dufld(2,i) = dufld(2,i) + xr*tyi5 + yr*txi5
     &                         + 2.0d0*xr*yr*turi
               dufld(3,i) = dufld(3,i) + yr*tyi5 + yr*yr*turi
               dufld(4,i) = dufld(4,i) + xr*tzi5 + zr*txi5
     &                         + 2.0d0*xr*zr*turi
               dufld(5,i) = dufld(5,i) + yr*tzi5 + zr*tyi5
     &                         + 2.0d0*yr*zr*turi
               dufld(6,i) = dufld(6,i) + zr*tzi5 + zr*zr*turi
               dufld(1,k) = dufld(1,k) - xr*txk5 - xr*xr*turk
               dufld(2,k) = dufld(2,k) - xr*tyk5 - yr*txk5
     &                         - 2.0d0*xr*yr*turk
               dufld(3,k) = dufld(3,k) - yr*tyk5 - yr*yr*turk
               dufld(4,k) = dufld(4,k) - xr*tzk5 - zr*txk5
     &                         - 2.0d0*xr*zr*turk
               dufld(5,k) = dufld(5,k) - yr*tzk5 - zr*tyk5
     &                         - 2.0d0*yr*zr*turk
               dufld(6,k) = dufld(6,k) - zr*tzk5 - zr*zr*turk
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
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
         pdi = pdamp(i)
         pti = thole(i)
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
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
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
               pscale(kk) = 1.0d0
               dscale(kk) = 1.0d0
               uscale(kk) = 1.0d0
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
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     sc3 = 1.0d0 - expdamp
                     sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                     sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                    *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     rc3(1) = xr * temp3
                     rc3(2) = yr * temp3
                     rc3(3) = zr * temp3
                     rc5(1) = rc3(1) * temp5
                     rc5(2) = rc3(2) * temp5
                     rc5(3) = rc3(3) * temp5
                     rc7(1) = rc5(1) * temp7
                     rc7(2) = rc5(2) * temp7
                     rc7(3) = rc5(3) * temp7
                  end if
               end if
c
c     intermediates involving Thole damping and scale factors
c
               psc3 = 1.0d0 - sc3*pscale(kk)
               psc5 = 1.0d0 - sc5*pscale(kk)
               psc7 = 1.0d0 - sc7*pscale(kk)
               dsc3 = 1.0d0 - sc3*dscale(kk)
               dsc5 = 1.0d0 - sc5*dscale(kk)
               dsc7 = 1.0d0 - sc7*dscale(kk)
               usc3 = 1.0d0 - sc3*uscale(kk)
               usc5 = 1.0d0 - sc5*uscale(kk)
               psr3 = bn(1) - psc3*rr3
               psr5 = bn(2) - psc5*rr5
               psr7 = bn(3) - psc7*rr7
               dsr3 = bn(1) - dsc3*rr3
               dsr5 = bn(2) - dsc5*rr5
               dsr7 = bn(3) - dsc7*rr7
               usr3 = bn(1) - usc3*rr3
               usr5 = bn(2) - usc5*rr5
               do j = 1, 3
                  prc3(j) = rc3(j) * pscale(kk)
                  prc5(j) = rc5(j) * pscale(kk)
                  prc7(j) = rc7(j) * pscale(kk)
                  drc3(j) = rc3(j) * dscale(kk)
                  drc5(j) = rc5(j) * dscale(kk)
                  drc7(j) = rc7(j) * dscale(kk)
                  urc3(j) = rc3(j) * uscale(kk)
                  urc5(j) = rc5(j) * uscale(kk)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               qri1 = qixx*xr + qixy*yr + qixz*zr
               qri2 = qixy*xr + qiyy*yr + qiyz*zr
               qri3 = qixz*xr + qiyz*yr + qizz*zr
               qrk1 = qkxx*xr + qkxy*yr + qkxz*zr
               qrk2 = qkxy*xr + qkyy*yr + qkyz*zr
               qrk3 = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = (qixx*xr+2.0d0*qixy*yr)*xr
     &                   + (qiyy*yr+2.0d0*qiyz*zr)*yr
     &                   + (qizz*zr+2.0d0*qixz*xr)*zr
               qrrk = (qkxx*xr+2.0d0*qkxy*yr)*xr
     &                   + (qkyy*yr+2.0d0*qkyz*zr)*yr
     &                   + (qkzz*zr+2.0d0*qkxz*xr)*zr
               puri = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               purk = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               duri = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               durk = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
c
c     get the dEd/dR field gradient for direct polarization force
c
               term1 = bn(2) - dsc3*rr5
               term2 = bn(3) - dsc5*rr7
               term3 = -dsr3 + term1*xr*xr - rr3*xr*drc3(1)
               term4 = rr3*drc3(1) - term1*xr - dsr5*xr
               term5 = term2*xr*xr - dsr5 - rr5*xr*drc5(1)
               term6 = (bn(4)-dsc7*rr9)*xr*xr - bn(3) - rr7*xr*drc7(1)
               term7 = rr5*drc5(1) - 2.0d0*bn(3)*xr
     &                    + (dsc5+1.5d0*dsc7)*rr7*xr
               txxi = ci*term3 + dix*term4 + dri*term5
     &                   + 2.0d0*dsr5*qixx + (qri2*yr+qri3*zr)*dsc7*rr7
     &                   + 2.0d0*qri1*term7 + qrri*term6
               txxk = ck*term3 - dkx*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkxx + (qrk2*yr+qrk3*zr)*dsc7*rr7
     &                   + 2.0d0*qrk1*term7 + qrrk*term6
               term3 = -dsr3 + term1*yr*yr - rr3*yr*drc3(2)
               term4 = rr3*drc3(2) - term1*yr - dsr5*yr
               term5 = term2*yr*yr - dsr5 - rr5*yr*drc5(2)
               term6 = (bn(4)-dsc7*rr9)*yr*yr - bn(3) - rr7*yr*drc7(2)
               term7 = rr5*drc5(2) - 2.0d0*bn(3)*yr
     &                    + (dsc5+1.5d0*dsc7)*rr7*yr
               tyyi = ci*term3 + diy*term4 + dri*term5
     &                   + 2.0d0*dsr5*qiyy + (qri1*xr+qri3*zr)*dsc7*rr7
     &                   + 2.0d0*qri2*term7 + qrri*term6
               tyyk = ck*term3 - dky*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkyy + (qrk1*xr+qrk3*zr)*dsc7*rr7
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = -dsr3 + term1*zr*zr - rr3*zr*drc3(3)
               term4 = rr3*drc3(3) - term1*zr - dsr5*zr
               term5 = term2*zr*zr - dsr5 - rr5*zr*drc5(3)
               term6 = (bn(4)-dsc7*rr9)*zr*zr - bn(3) - rr7*zr*drc7(3)
               term7 = rr5*drc5(3) - 2.0d0*bn(3)*zr
     &                    + (dsc5+1.5d0*dsc7)*rr7*zr
               tzzi = ci*term3 + diz*term4 + dri*term5
     &                   + 2.0d0*dsr5*qizz + (qri1*xr+qri2*yr)*dsc7*rr7
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tzzk = ck*term3 - dkz*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkzz + (qrk1*xr+qrk2*yr)*dsc7*rr7
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*xr*yr - rr3*yr*drc3(1)
               term4 = rr3*drc3(1) - term1*xr
               term5 = term2*xr*yr - rr5*yr*drc5(1)
               term6 = (bn(4)-dsc7*rr9)*xr*yr - rr7*yr*drc7(1)
               term7 = rr5*drc5(1) - term2*xr
               txyi = ci*term3 - dsr5*dix*yr + diy*term4 + dri*term5
     &                   + 2.0d0*dsr5*qixy - 2.0d0*dsr7*yr*qri1
     &                   + 2.0d0*qri2*term7 + qrri*term6
               txyk = ck*term3 + dsr5*dkx*yr - dky*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkxy - 2.0d0*dsr7*yr*qrk1
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = term1*xr*zr - rr3*zr*drc3(1)
               term5 = term2*xr*zr - rr5*zr*drc5(1)
               term6 = (bn(4)-dsc7*rr9)*xr*zr - rr7*zr*drc7(1)
               txzi = ci*term3 - dsr5*dix*zr + diz*term4 + dri*term5
     &                   + 2.0d0*dsr5*qixz - 2.0d0*dsr7*zr*qri1
     &                   + 2.0d0*qri3*term7 + qrri*term6
               txzk = ck*term3 + dsr5*dkx*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkxz - 2.0d0*dsr7*zr*qrk1
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*yr*zr - rr3*zr*drc3(2)
               term4 = rr3*drc3(2) - term1*yr
               term5 = term2*yr*zr - rr5*zr*drc5(2)
               term6 = (bn(4)-dsc7*rr9)*yr*zr - rr7*zr*drc7(2)
               term7 = rr5*drc5(2) - term2*yr
               tyzi = ci*term3 - dsr5*diy*zr + diz*term4 + dri*term5
     &                   + 2.0d0*dsr5*qiyz - 2.0d0*dsr7*zr*qri2
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tyzk = ck*term3 + dsr5*dky*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkyz - 2.0d0*dsr7*zr*qrk2
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               depx = txxi*uinp(1,k) - txxk*uinp(1,i)
     &                   + txyi*uinp(2,k) - txyk*uinp(2,i)
     &                   + txzi*uinp(3,k) - txzk*uinp(3,i)
               depy = txyi*uinp(1,k) - txyk*uinp(1,i)
     &                   + tyyi*uinp(2,k) - tyyk*uinp(2,i)
     &                   + tyzi*uinp(3,k) - tyzk*uinp(3,i)
               depz = txzi*uinp(1,k) - txzk*uinp(1,i)
     &                   + tyzi*uinp(2,k) - tyzk*uinp(2,i)
     &                   + tzzi*uinp(3,k) - tzzk*uinp(3,i)
               frcx = depx
               frcy = depy
               frcz = depz
c
c     get the dEp/dR field gradient for direct polarization force
c
               term1 = bn(2) - psc3*rr5
               term2 = bn(3) - psc5*rr7
               term3 = -psr3 + term1*xr*xr - rr3*xr*prc3(1)
               term4 = rr3*prc3(1) - term1*xr - psr5*xr
               term5 = term2*xr*xr - psr5 - rr5*xr*prc5(1)
               term6 = (bn(4)-psc7*rr9)*xr*xr - bn(3) - rr7*xr*prc7(1)
               term7 = rr5*prc5(1) - 2.0d0*bn(3)*xr
     &                    + (psc5+1.5d0*psc7)*rr7*xr
               txxi = ci*term3 + dix*term4 + dri*term5
     &                   + 2.0d0*psr5*qixx + (qri2*yr+qri3*zr)*psc7*rr7
     &                   + 2.0d0*qri1*term7 + qrri*term6
               txxk = ck*term3 - dkx*term4 - drk*term5
     &                   + 2.0d0*psr5*qkxx + (qrk2*yr+qrk3*zr)*psc7*rr7
     &                   + 2.0d0*qrk1*term7 + qrrk*term6
               term3 = -psr3 + term1*yr*yr - rr3*yr*prc3(2)
               term4 = rr3*prc3(2) - term1*yr - psr5*yr
               term5 = term2*yr*yr - psr5 - rr5*yr*prc5(2)
               term6 = (bn(4)-psc7*rr9)*yr*yr - bn(3) - rr7*yr*prc7(2)
               term7 = rr5*prc5(2) - 2.0d0*bn(3)*yr
     &                    + (psc5+1.5d0*psc7)*rr7*yr
               tyyi = ci*term3 + diy*term4 + dri*term5
     &                   + 2.0d0*psr5*qiyy + (qri1*xr+qri3*zr)*psc7*rr7
     &                   + 2.0d0*qri2*term7 + qrri*term6
               tyyk = ck*term3 - dky*term4 - drk*term5
     &                   + 2.0d0*psr5*qkyy + (qrk1*xr+qrk3*zr)*psc7*rr7
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = -psr3 + term1*zr*zr - rr3*zr*prc3(3)
               term4 = rr3*prc3(3) - term1*zr - psr5*zr
               term5 = term2*zr*zr - psr5 - rr5*zr*prc5(3)
               term6 = (bn(4)-psc7*rr9)*zr*zr - bn(3) - rr7*zr*prc7(3)
               term7 = rr5*prc5(3) - 2.0d0*bn(3)*zr
     &                    + (psc5+1.5d0*psc7)*rr7*zr
               tzzi = ci*term3 + diz*term4 + dri*term5
     &                   + 2.0d0*psr5*qizz + (qri1*xr+qri2*yr)*psc7*rr7
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tzzk = ck*term3 - dkz*term4 - drk*term5
     &                   + 2.0d0*psr5*qkzz + (qrk1*xr+qrk2*yr)*psc7*rr7
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*xr*yr - rr3*yr*prc3(1)
               term4 = rr3*prc3(1) - term1*xr
               term5 = term2*xr*yr - rr5*yr*prc5(1)
               term6 = (bn(4)-psc7*rr9)*xr*yr - rr7*yr*prc7(1)
               term7 = rr5*prc5(1) - term2*xr
               txyi = ci*term3 - psr5*dix*yr + diy*term4 + dri*term5
     &                   + 2.0d0*psr5*qixy - 2.0d0*psr7*yr*qri1
     &                   + 2.0d0*qri2*term7 + qrri*term6
               txyk = ck*term3 + psr5*dkx*yr - dky*term4 - drk*term5
     &                   + 2.0d0*psr5*qkxy - 2.0d0*psr7*yr*qrk1
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = term1*xr*zr - rr3*zr*prc3(1)
               term5 = term2*xr*zr - rr5*zr*prc5(1)
               term6 = (bn(4)-psc7*rr9)*xr*zr - rr7*zr*prc7(1)
               txzi = ci*term3 - psr5*dix*zr + diz*term4 + dri*term5
     &                   + 2.0d0*psr5*qixz - 2.0d0*psr7*zr*qri1
     &                   + 2.0d0*qri3*term7 + qrri*term6
               txzk = ck*term3 + psr5*dkx*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*psr5*qkxz - 2.0d0*psr7*zr*qrk1
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*yr*zr - rr3*zr*prc3(2)
               term4 = rr3*prc3(2) - term1*yr
               term5 = term2*yr*zr - rr5*zr*prc5(2)
               term6 = (bn(4)-psc7*rr9)*yr*zr - rr7*zr*prc7(2)
               term7 = rr5*prc5(2) - term2*yr
               tyzi = ci*term3 - psr5*diy*zr + diz*term4 + dri*term5
     &                   + 2.0d0*psr5*qiyz - 2.0d0*psr7*zr*qri2
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tyzk = ck*term3 + psr5*dky*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*psr5*qkyz - 2.0d0*psr7*zr*qrk2
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               depx = txxi*uind(1,k) - txxk*uind(1,i)
     &                   + txyi*uind(2,k) - txyk*uind(2,i)
     &                   + txzi*uind(3,k) - txzk*uind(3,i)
               depy = txyi*uind(1,k) - txyk*uind(1,i)
     &                   + tyyi*uind(2,k) - tyyk*uind(2,i)
     &                   + tyzi*uind(3,k) - tyzk*uind(3,i)
               depz = txzi*uind(1,k) - txzk*uind(1,i)
     &                   + tyzi*uind(2,k) - tyzk*uind(2,i)
     &                   + tzzi*uind(3,k) - tzzk*uind(3,i)
               frcx = frcx + depx
               frcy = frcy + depy
               frcz = frcz + depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp .eq. 'MUTUAL') then
                  term1 = bn(2) - usc3*rr5
                  term2 = bn(3) - usc5*rr7
                  term3 = usr5 + term1
                  term4 = rr3 * uscale(kk)
                  term5 = -xr*term3 + rc3(1)*term4
                  term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                  txxi = uind(1,i)*term5 + duri*term6
                  txxk = uind(1,k)*term5 + durk*term6
                  term5 = -yr*term3 + rc3(2)*term4
                  term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                  tyyi = uind(2,i)*term5 + duri*term6
                  tyyk = uind(2,k)*term5 + durk*term6
                  term5 = -zr*term3 + rc3(3)*term4
                  term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                  tzzi = uind(3,i)*term5 + duri*term6
                  tzzk = uind(3,k)*term5 + durk*term6
                  term4 = -usr5 * yr
                  term5 = -xr*term1 + rr3*urc3(1)
                  term6 = xr*yr*term2 - rr5*yr*urc5(1)
                  txyi = uind(1,i)*term4 + uind(2,i)*term5 + duri*term6
                  txyk = uind(1,k)*term4 + uind(2,k)*term5 + durk*term6
                  term4 = -usr5 * zr
                  term6 = xr*zr*term2 - rr5*zr*urc5(1)
                  txzi = uind(1,i)*term4 + uind(3,i)*term5 + duri*term6
                  txzk = uind(1,k)*term4 + uind(3,k)*term5 + durk*term6
                  term5 = -yr*term1 + rr3*urc3(2)
                  term6 = yr*zr*term2 - rr5*zr*urc5(2)
                  tyzi = uind(2,i)*term4 + uind(3,i)*term5 + duri*term6
                  tyzk = uind(2,k)*term4 + uind(3,k)*term5 + durk*term6
                  depx = txxi*uinp(1,k) + txxk*uinp(1,i)
     &                      + txyi*uinp(2,k) + txyk*uinp(2,i)
     &                      + txzi*uinp(3,k) + txzk*uinp(3,i)
                  depy = txyi*uinp(1,k) + txyk*uinp(1,i)
     &                      + tyyi*uinp(2,k) + tyyk*uinp(2,i)
     &                      + tyzi*uinp(3,k) + tyzk*uinp(3,i)
                  depz = txzi*uinp(1,k) + txzk*uinp(1,i)
     &                      + tyzi*uinp(2,k) + tyzk*uinp(2,i)
     &                      + tzzi*uinp(3,k) + tzzk*uinp(3,i)
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
               end if
c
c     force and torque components scaled for self-interactions
c
               if (ii .eq. kk) then
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  do j = 1, 3
                     psr3 = 0.5d0 * psr3
                     psr5 = 0.5d0 * psr5
                     psr7 = 0.5d0 * psr7
                     dsr3 = 0.5d0 * dsr3
                     dsr5 = 0.5d0 * dsr5
                     dsr7 = 0.5d0 * dsr7
                  end do
               end if
c
c     increment gradient and virial due to Cartesian forces
c
               dep(1,ii) = dep(1,ii) - frcx
               dep(2,ii) = dep(2,ii) - frcy
               dep(3,ii) = dep(3,ii) - frcz
               dep(1,kk) = dep(1,kk) + frcx
               dep(2,kk) = dep(2,kk) + frcy
               dep(3,kk) = dep(3,kk) + frcz
               vxx = xr * frcx
               vxy = yr * frcx
               vxz = zr * frcx
               vyy = yr * frcy
               vyz = zr * frcy
               vzz = zr * frcz
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
c
c     get the induced dipole field used for dipole torques
c
               txi3 = psr3*uind(1,k) + dsr3*uinp(1,k)
               tyi3 = psr3*uind(2,k) + dsr3*uinp(2,k)
               tzi3 = psr3*uind(3,k) + dsr3*uinp(3,k)
               txk3 = psr3*uind(1,i) + dsr3*uinp(1,i)
               tyk3 = psr3*uind(2,i) + dsr3*uinp(2,i)
               tzk3 = psr3*uind(3,i) + dsr3*uinp(3,i)
               turi = -psr5*durk - dsr5*purk
               turk = -psr5*duri - dsr5*puri
               ufld(1,i) = ufld(1,i) + txi3 + xr*turi
               ufld(2,i) = ufld(2,i) + tyi3 + yr*turi
               ufld(3,i) = ufld(3,i) + tzi3 + zr*turi
               ufld(1,k) = ufld(1,k) + txk3 + xr*turk
               ufld(2,k) = ufld(2,k) + tyk3 + yr*turk
               ufld(3,k) = ufld(3,k) + tzk3 + zr*turk
c
c     get induced dipole field gradient used for quadrupole torques
c
               txi5 = 2.0d0 * (psr5*uind(1,k)+dsr5*uinp(1,k))
               tyi5 = 2.0d0 * (psr5*uind(2,k)+dsr5*uinp(2,k))
               tzi5 = 2.0d0 * (psr5*uind(3,k)+dsr5*uinp(3,k))
               txk5 = 2.0d0 * (psr5*uind(1,i)+dsr5*uinp(1,i))
               tyk5 = 2.0d0 * (psr5*uind(2,i)+dsr5*uinp(2,i))
               tzk5 = 2.0d0 * (psr5*uind(3,i)+dsr5*uinp(3,i))
               turi = -psr7*durk - dsr7*purk
               turk = -psr7*duri - dsr7*puri
               dufld(1,i) = dufld(1,i) + xr*txi5 + xr*xr*turi
               dufld(2,i) = dufld(2,i) + xr*tyi5 + yr*txi5
     &                         + 2.0d0*xr*yr*turi
               dufld(3,i) = dufld(3,i) + yr*tyi5 + yr*yr*turi
               dufld(4,i) = dufld(4,i) + xr*tzi5 + zr*txi5
     &                         + 2.0d0*xr*zr*turi
               dufld(5,i) = dufld(5,i) + yr*tzi5 + zr*tyi5
     &                         + 2.0d0*yr*zr*turi
               dufld(6,i) = dufld(6,i) + zr*tzi5 + zr*zr*turi
               dufld(1,k) = dufld(1,k) - xr*txk5 - xr*xr*turk
               dufld(2,k) = dufld(2,k) - xr*tyk5 - yr*txk5
     &                         - 2.0d0*xr*yr*turk
               dufld(3,k) = dufld(3,k) - yr*tyk5 - yr*yr*turk
               dufld(4,k) = dufld(4,k) - xr*tzk5 - zr*txk5
     &                         - 2.0d0*xr*zr*turk
               dufld(5,k) = dufld(5,k) - yr*tzk5 - zr*tyk5
     &                         - 2.0d0*yr*zr*turk
               dufld(6,k) = dufld(6,k) - zr*tzk5 - zr*zr*turk
            end if
            end do
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
      end if
c
c     torque is induced field and gradient cross permanent moments
c
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         trq(1) = diz*ufld(2,i) - diy*ufld(3,i)
     &               + qixz*dufld(2,i) - qixy*dufld(4,i)
     &               + 2.0d0*qiyz*(dufld(3,i)-dufld(6,i))
     &               + (qizz-qiyy)*dufld(5,i)
         trq(2) = dix*ufld(3,i) - diz*ufld(1,i)
     &               - qiyz*dufld(2,i) + qixy*dufld(5,i)
     &               + 2.0d0*qixz*(dufld(6,i)-dufld(1,i))
     &               + (qixx-qizz)*dufld(4,i)
         trq(3) = diy*ufld(1,i) - dix*ufld(2,i)
     &               + qiyz*dufld(4,i) - qixz*dufld(5,i)
     &               + 2.0d0*qixy*(dufld(1,i)-dufld(3,i))
     &               + (qiyy-qixx)*dufld(2,i)
         call torque (i,trq,fix,fiy,fiz,dep)
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
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (ufld)
      deallocate (dufld)
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar1d  --  Ewald polarization derivs via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar1d" calculates the dipole polarization energy and
c     derivatives with respect to Cartesian coordinates using
c     particle mesh Ewald summation and a neighbor list
c
c
      subroutine epolar1d
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use polar
      use polpot
      use potent
      use virial
      implicit none
      integer i,j,ii
      real*8 e,f,term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xv,yv,zv,vterm
      real*8 xufield,yufield
      real*8 zufield
      real*8 fix(3),fiy(3),fiz(3)
      real*8 trq(3)
c
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0d0
      do i = 1, n
         do j = 1, 3
            dep(j,i) = 0.0d0
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
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     compute the real space part of the Ewald summation
c
      call epreal1d
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip1
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * term * uii / 3.0d0
c        ep = ep + e
      end do
c
c     compute the self-energy torque term due to induced dipole
c
      term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = 0.5d0 * (uind(1,i)+uinp(1,i))
         uiy = 0.5d0 * (uind(2,i)+uinp(2,i))
         uiz = 0.5d0 * (uind(3,i)+uinp(3,i))
         trq(1) = term * (diy*uiz-diz*uiy)
         trq(2) = term * (diz*uix-dix*uiz)
         trq(3) = term * (dix*uiy-diy*uix)
         call torque (i,trq,fix,fiy,fiz,dep)
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         xup = 0.0d0
         yup = 0.0d0
         zup = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
            xu = xu + uind(1,i)
            yu = yu + uind(2,i)
            zu = zu + uind(3,i)
            xup = xup + uinp(1,i)
            yup = yup + uinp(2,i)
            zup = zup + uinp(3,i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
c        ep = ep + term*(xd*xu+yd*yu+zd*zu)
         do i = 1, npole
            ii = ipole(i)
            dep(1,ii) = dep(1,ii) + term*rpole(1,i)*(xu+xup)
            dep(2,ii) = dep(2,ii) + term*rpole(1,i)*(yu+yup)
            dep(3,ii) = dep(3,ii) + term*rpole(1,i)*(zu+zup)
         end do
         xufield = -term * (xu+xup)
         yufield = -term * (yu+yup)
         zufield = -term * (zu+zup)
         do i = 1, npole
            trq(1) = rpole(3,i)*zufield - rpole(4,i)*yufield
            trq(2) = rpole(4,i)*xufield - rpole(2,i)*zufield
            trq(3) = rpole(2,i)*yufield - rpole(3,i)*xufield
            call torque (i,trq,fix,fiy,fiz,dep)
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
         xv = xq * (xu+xup)
         yv = yq * (yu+yup)
         zv = zq * (zu+zup)
         vterm = xv + yv + zv + xu*xup + yu*yup + zu*zup
     &              + xd*(xu+xup) + yd*(yu+yup) + zd*(zu+zup)
         vterm = term * vterm
         vir(1,1) = vir(1,1) + term*xv + vterm
         vir(2,1) = vir(2,1) + term*xv
         vir(3,1) = vir(3,1) + term*xv
         vir(1,2) = vir(1,2) + term*yv
         vir(2,2) = vir(2,2) + term*yv + vterm
         vir(3,2) = vir(3,2) + term*yv
         vir(1,3) = vir(1,3) + term*zv
         vir(2,3) = vir(2,3) + term*zv
         vir(3,3) = vir(3,3) + term*zv + vterm
         if (poltyp .eq. 'DIRECT') then
            vterm = term * (xu*xup+yu*yup+zu*zup)
            vir(1,1) = vir(1,1) + vterm
            vir(2,2) = vir(2,2) + vterm
            vir(3,3) = vir(3,3) + vterm
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal1d  --  Ewald real space derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal1d" evaluates the real space portion of the Ewald
c     summation energy and gradient due to dipole polarization
c     via a neighbor list
c
c
      subroutine epreal1d
      use sizes
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer iax,iay,iaz
      real*8 e,epo,f,fi
      real*8 erfc,bfac
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 temp3,temp5,temp7
      real*8 scalekk
      real*8 sc3,sc5,sc7
      real*8 psc3,psc5,psc7
      real*8 dsc3,dsc5,dsc7
      real*8 usc3,usc5
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 usr3,usr5
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dri,drk
      real*8 qri1,qri2,qri3
      real*8 qrk1,qrk2,qrk3
      real*8 qrri,qrrk
      real*8 txxi,tyyi,tzzi
      real*8 txyi,txzi,tyzi
      real*8 txxk,tyyk,tzzk
      real*8 txyk,txzk,tyzk
      real*8 turi,duri,puri
      real*8 turk,durk,purk
      real*8 txi3,tyi3,tzi3
      real*8 txi5,tyi5,tzi5
      real*8 txk3,tyk3,tzk3
      real*8 txk5,tyk5,tzk5
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term6,term7
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 prc3(3),prc5(3),prc7(3)
      real*8 drc3(3),drc5(3),drc7(3)
      real*8 urc3(3),urc5(3)
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8 bn(0:4)
      real*8 viro(3,3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: depo(:,:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (depo(3,n))
      allocate (ufld(3,n))
      allocate (dufld(6,n))
c
c     set arrays needed to scale interactions and store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
         do j = 1, 3
            ufld(j,i) = 0.0d0
         end do
         do j = 1, 6
            dufld(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     initialize local variables for OpenMP calculation
c
      epo = 0.0d0
      do i = 1, n
         do j = 1, 3
            depo(j,i) = 0.0d0
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
!$OMP PARALLEL default(private) shared(npole,polarity,f,uind,udirp,epo)
!$OMP& shared(ipole,pdamp,thole,x,y,z,xaxis,yaxis,zaxis,rpole,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,
!$OMP& p2scale,p3scale,p4scale,p41scale,p5scale,d1scale,d2scale,d3scale,
!$OMP& d4scale,u1scale,u2scale,u3scale,u4scale,nelst,elst,use_bounds,
!$OMP& off2,aewald,uinp,poltyp,depo,viro,ufld,dufld)
!$OMP& firstprivate(pscale,dscale,uscale)
!$OMP DO reduction(+:epo) schedule(guided)
c
c     get polarization energy via induced dipoles times field
c
      do i = 1, npole
         if (polarity(i) .ne. 0.0d0) then
            fi = -f / polarity(i)
            e = 0.0d0
            do j = 1, 3
               e = fi * uind(j,i) * udirp(j,i)
               epo = epo + e
            end do
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO reduction(+:depo,viro,ufld,dufld) schedule(guided)
c
c     compute the dipole polarization gradient components
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
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
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
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
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 4
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 4
                  bn(j) = f * bn(j)
               end do
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     sc3 = 1.0d0 - expdamp
                     sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                     sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                    *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     rc3(1) = xr * temp3
                     rc3(2) = yr * temp3
                     rc3(3) = zr * temp3
                     rc5(1) = rc3(1) * temp5
                     rc5(2) = rc3(2) * temp5
                     rc5(3) = rc3(3) * temp5
                     rc7(1) = rc5(1) * temp7
                     rc7(2) = rc5(2) * temp7
                     rc7(3) = rc5(3) * temp7
                  end if
               end if
c
c     intermediates involving Thole damping and scale factors
c
               psc3 = 1.0d0 - sc3*pscale(kk)
               psc5 = 1.0d0 - sc5*pscale(kk)
               psc7 = 1.0d0 - sc7*pscale(kk)
               dsc3 = 1.0d0 - sc3*dscale(kk)
               dsc5 = 1.0d0 - sc5*dscale(kk)
               dsc7 = 1.0d0 - sc7*dscale(kk)
               usc3 = 1.0d0 - sc3*uscale(kk)
               usc5 = 1.0d0 - sc5*uscale(kk)
               psr3 = bn(1) - psc3*rr3
               psr5 = bn(2) - psc5*rr5
               psr7 = bn(3) - psc7*rr7
               dsr3 = bn(1) - dsc3*rr3
               dsr5 = bn(2) - dsc5*rr5
               dsr7 = bn(3) - dsc7*rr7
               usr3 = bn(1) - usc3*rr3
               usr5 = bn(2) - usc5*rr5
               do j = 1, 3
                  prc3(j) = rc3(j) * pscale(kk)
                  prc5(j) = rc5(j) * pscale(kk)
                  prc7(j) = rc7(j) * pscale(kk)
                  drc3(j) = rc3(j) * dscale(kk)
                  drc5(j) = rc5(j) * dscale(kk)
                  drc7(j) = rc7(j) * dscale(kk)
                  urc3(j) = rc3(j) * uscale(kk)
                  urc5(j) = rc5(j) * uscale(kk)
               end do
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               qri1 = qixx*xr + qixy*yr + qixz*zr
               qri2 = qixy*xr + qiyy*yr + qiyz*zr
               qri3 = qixz*xr + qiyz*yr + qizz*zr
               qrk1 = qkxx*xr + qkxy*yr + qkxz*zr
               qrk2 = qkxy*xr + qkyy*yr + qkyz*zr
               qrk3 = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = (qixx*xr+2.0d0*qixy*yr)*xr
     &                   + (qiyy*yr+2.0d0*qiyz*zr)*yr
     &                   + (qizz*zr+2.0d0*qixz*xr)*zr
               qrrk = (qkxx*xr+2.0d0*qkxy*yr)*xr
     &                   + (qkyy*yr+2.0d0*qkyz*zr)*yr
     &                   + (qkzz*zr+2.0d0*qkxz*xr)*zr
               puri = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               purk = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               duri = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               durk = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
c
c     get the dEd/dR field gradient for direct polarization force
c
               term1 = bn(2) - dsc3*rr5
               term2 = bn(3) - dsc5*rr7
               term3 = -dsr3 + term1*xr*xr - rr3*xr*drc3(1)
               term4 = rr3*drc3(1) - term1*xr - dsr5*xr
               term5 = term2*xr*xr - dsr5 - rr5*xr*drc5(1)
               term6 = (bn(4)-dsc7*rr9)*xr*xr - bn(3) - rr7*xr*drc7(1)
               term7 = rr5*drc5(1) - 2.0d0*bn(3)*xr
     &                    + (dsc5+1.5d0*dsc7)*rr7*xr
               txxi = ci*term3 + dix*term4 + dri*term5
     &                   + 2.0d0*dsr5*qixx + (qri2*yr+qri3*zr)*dsc7*rr7
     &                   + 2.0d0*qri1*term7 + qrri*term6
               txxk = ck*term3 - dkx*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkxx + (qrk2*yr+qrk3*zr)*dsc7*rr7
     &                   + 2.0d0*qrk1*term7 + qrrk*term6
               term3 = -dsr3 + term1*yr*yr - rr3*yr*drc3(2)
               term4 = rr3*drc3(2) - term1*yr - dsr5*yr
               term5 = term2*yr*yr - dsr5 - rr5*yr*drc5(2)
               term6 = (bn(4)-dsc7*rr9)*yr*yr - bn(3) - rr7*yr*drc7(2)
               term7 = rr5*drc5(2) - 2.0d0*bn(3)*yr
     &                    + (dsc5+1.5d0*dsc7)*rr7*yr
               tyyi = ci*term3 + diy*term4 + dri*term5
     &                   + 2.0d0*dsr5*qiyy + (qri1*xr+qri3*zr)*dsc7*rr7
     &                   + 2.0d0*qri2*term7 + qrri*term6
               tyyk = ck*term3 - dky*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkyy + (qrk1*xr+qrk3*zr)*dsc7*rr7
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = -dsr3 + term1*zr*zr - rr3*zr*drc3(3)
               term4 = rr3*drc3(3) - term1*zr - dsr5*zr
               term5 = term2*zr*zr - dsr5 - rr5*zr*drc5(3)
               term6 = (bn(4)-dsc7*rr9)*zr*zr - bn(3) - rr7*zr*drc7(3)
               term7 = rr5*drc5(3) - 2.0d0*bn(3)*zr
     &                    + (dsc5+1.5d0*dsc7)*rr7*zr
               tzzi = ci*term3 + diz*term4 + dri*term5
     &                   + 2.0d0*dsr5*qizz + (qri1*xr+qri2*yr)*dsc7*rr7
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tzzk = ck*term3 - dkz*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkzz + (qrk1*xr+qrk2*yr)*dsc7*rr7
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*xr*yr - rr3*yr*drc3(1)
               term4 = rr3*drc3(1) - term1*xr
               term5 = term2*xr*yr - rr5*yr*drc5(1)
               term6 = (bn(4)-dsc7*rr9)*xr*yr - rr7*yr*drc7(1)
               term7 = rr5*drc5(1) - term2*xr
               txyi = ci*term3 - dsr5*dix*yr + diy*term4 + dri*term5
     &                   + 2.0d0*dsr5*qixy - 2.0d0*dsr7*yr*qri1
     &                   + 2.0d0*qri2*term7 + qrri*term6
               txyk = ck*term3 + dsr5*dkx*yr - dky*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkxy - 2.0d0*dsr7*yr*qrk1
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = term1*xr*zr - rr3*zr*drc3(1)
               term5 = term2*xr*zr - rr5*zr*drc5(1)
               term6 = (bn(4)-dsc7*rr9)*xr*zr - rr7*zr*drc7(1)
               txzi = ci*term3 - dsr5*dix*zr + diz*term4 + dri*term5
     &                   + 2.0d0*dsr5*qixz - 2.0d0*dsr7*zr*qri1
     &                   + 2.0d0*qri3*term7 + qrri*term6
               txzk = ck*term3 + dsr5*dkx*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkxz - 2.0d0*dsr7*zr*qrk1
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*yr*zr - rr3*zr*drc3(2)
               term4 = rr3*drc3(2) - term1*yr
               term5 = term2*yr*zr - rr5*zr*drc5(2)
               term6 = (bn(4)-dsc7*rr9)*yr*zr - rr7*zr*drc7(2)
               term7 = rr5*drc5(2) - term2*yr
               tyzi = ci*term3 - dsr5*diy*zr + diz*term4 + dri*term5
     &                   + 2.0d0*dsr5*qiyz - 2.0d0*dsr7*zr*qri2
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tyzk = ck*term3 + dsr5*dky*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*dsr5*qkyz - 2.0d0*dsr7*zr*qrk2
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               depx = txxi*uinp(1,k) - txxk*uinp(1,i)
     &                   + txyi*uinp(2,k) - txyk*uinp(2,i)
     &                   + txzi*uinp(3,k) - txzk*uinp(3,i)
               depy = txyi*uinp(1,k) - txyk*uinp(1,i)
     &                   + tyyi*uinp(2,k) - tyyk*uinp(2,i)
     &                   + tyzi*uinp(3,k) - tyzk*uinp(3,i)
               depz = txzi*uinp(1,k) - txzk*uinp(1,i)
     &                   + tyzi*uinp(2,k) - tyzk*uinp(2,i)
     &                   + tzzi*uinp(3,k) - tzzk*uinp(3,i)
               frcx = depx
               frcy = depy
               frcz = depz
c
c     get the dEp/dR field gradient for direct polarization force
c
               term1 = bn(2) - psc3*rr5
               term2 = bn(3) - psc5*rr7
               term3 = -psr3 + term1*xr*xr - rr3*xr*prc3(1)
               term4 = rr3*prc3(1) - term1*xr - psr5*xr
               term5 = term2*xr*xr - psr5 - rr5*xr*prc5(1)
               term6 = (bn(4)-psc7*rr9)*xr*xr - bn(3) - rr7*xr*prc7(1)
               term7 = rr5*prc5(1) - 2.0d0*bn(3)*xr
     &                    + (psc5+1.5d0*psc7)*rr7*xr
               txxi = ci*term3 + dix*term4 + dri*term5
     &                   + 2.0d0*psr5*qixx + (qri2*yr+qri3*zr)*psc7*rr7
     &                   + 2.0d0*qri1*term7 + qrri*term6
               txxk = ck*term3 - dkx*term4 - drk*term5
     &                   + 2.0d0*psr5*qkxx + (qrk2*yr+qrk3*zr)*psc7*rr7
     &                   + 2.0d0*qrk1*term7 + qrrk*term6
               term3 = -psr3 + term1*yr*yr - rr3*yr*prc3(2)
               term4 = rr3*prc3(2) - term1*yr - psr5*yr
               term5 = term2*yr*yr - psr5 - rr5*yr*prc5(2)
               term6 = (bn(4)-psc7*rr9)*yr*yr - bn(3) - rr7*yr*prc7(2)
               term7 = rr5*prc5(2) - 2.0d0*bn(3)*yr
     &                    + (psc5+1.5d0*psc7)*rr7*yr
               tyyi = ci*term3 + diy*term4 + dri*term5
     &                   + 2.0d0*psr5*qiyy + (qri1*xr+qri3*zr)*psc7*rr7
     &                   + 2.0d0*qri2*term7 + qrri*term6
               tyyk = ck*term3 - dky*term4 - drk*term5
     &                   + 2.0d0*psr5*qkyy + (qrk1*xr+qrk3*zr)*psc7*rr7
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = -psr3 + term1*zr*zr - rr3*zr*prc3(3)
               term4 = rr3*prc3(3) - term1*zr - psr5*zr
               term5 = term2*zr*zr - psr5 - rr5*zr*prc5(3)
               term6 = (bn(4)-psc7*rr9)*zr*zr - bn(3) - rr7*zr*prc7(3)
               term7 = rr5*prc5(3) - 2.0d0*bn(3)*zr
     &                    + (psc5+1.5d0*psc7)*rr7*zr
               tzzi = ci*term3 + diz*term4 + dri*term5
     &                   + 2.0d0*psr5*qizz + (qri1*xr+qri2*yr)*psc7*rr7
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tzzk = ck*term3 - dkz*term4 - drk*term5
     &                   + 2.0d0*psr5*qkzz + (qrk1*xr+qrk2*yr)*psc7*rr7
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*xr*yr - rr3*yr*prc3(1)
               term4 = rr3*prc3(1) - term1*xr
               term5 = term2*xr*yr - rr5*yr*prc5(1)
               term6 = (bn(4)-psc7*rr9)*xr*yr - rr7*yr*prc7(1)
               term7 = rr5*prc5(1) - term2*xr
               txyi = ci*term3 - psr5*dix*yr + diy*term4 + dri*term5
     &                   + 2.0d0*psr5*qixy - 2.0d0*psr7*yr*qri1
     &                   + 2.0d0*qri2*term7 + qrri*term6
               txyk = ck*term3 + psr5*dkx*yr - dky*term4 - drk*term5
     &                   + 2.0d0*psr5*qkxy - 2.0d0*psr7*yr*qrk1
     &                   + 2.0d0*qrk2*term7 + qrrk*term6
               term3 = term1*xr*zr - rr3*zr*prc3(1)
               term5 = term2*xr*zr - rr5*zr*prc5(1)
               term6 = (bn(4)-psc7*rr9)*xr*zr - rr7*zr*prc7(1)
               txzi = ci*term3 - psr5*dix*zr + diz*term4 + dri*term5
     &                   + 2.0d0*psr5*qixz - 2.0d0*psr7*zr*qri1
     &                   + 2.0d0*qri3*term7 + qrri*term6
               txzk = ck*term3 + psr5*dkx*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*psr5*qkxz - 2.0d0*psr7*zr*qrk1
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               term3 = term1*yr*zr - rr3*zr*prc3(2)
               term4 = rr3*prc3(2) - term1*yr
               term5 = term2*yr*zr - rr5*zr*prc5(2)
               term6 = (bn(4)-psc7*rr9)*yr*zr - rr7*zr*prc7(2)
               term7 = rr5*prc5(2) - term2*yr
               tyzi = ci*term3 - psr5*diy*zr + diz*term4 + dri*term5
     &                   + 2.0d0*psr5*qiyz - 2.0d0*psr7*zr*qri2
     &                   + 2.0d0*qri3*term7 + qrri*term6
               tyzk = ck*term3 + psr5*dky*zr - dkz*term4 - drk*term5
     &                   + 2.0d0*psr5*qkyz - 2.0d0*psr7*zr*qrk2
     &                   + 2.0d0*qrk3*term7 + qrrk*term6
               depx = txxi*uind(1,k) - txxk*uind(1,i)
     &                   + txyi*uind(2,k) - txyk*uind(2,i)
     &                   + txzi*uind(3,k) - txzk*uind(3,i)
               depy = txyi*uind(1,k) - txyk*uind(1,i)
     &                   + tyyi*uind(2,k) - tyyk*uind(2,i)
     &                   + tyzi*uind(3,k) - tyzk*uind(3,i)
               depz = txzi*uind(1,k) - txzk*uind(1,i)
     &                   + tyzi*uind(2,k) - tyzk*uind(2,i)
     &                   + tzzi*uind(3,k) - tzzk*uind(3,i)
               frcx = frcx + depx
               frcy = frcy + depy
               frcz = frcz + depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp .eq. 'MUTUAL') then
                  term1 = bn(2) - usc3*rr5
                  term2 = bn(3) - usc5*rr7
                  term3 = usr5 + term1
                  term4 = rr3 * uscale(kk)
                  term5 = -xr*term3 + rc3(1)*term4
                  term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                  txxi = uind(1,i)*term5 + duri*term6
                  txxk = uind(1,k)*term5 + durk*term6
                  term5 = -yr*term3 + rc3(2)*term4
                  term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                  tyyi = uind(2,i)*term5 + duri*term6
                  tyyk = uind(2,k)*term5 + durk*term6
                  term5 = -zr*term3 + rc3(3)*term4
                  term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                  tzzi = uind(3,i)*term5 + duri*term6
                  tzzk = uind(3,k)*term5 + durk*term6
                  term4 = -usr5 * yr
                  term5 = -xr*term1 + rr3*urc3(1)
                  term6 = xr*yr*term2 - rr5*yr*urc5(1)
                  txyi = uind(1,i)*term4 + uind(2,i)*term5 + duri*term6
                  txyk = uind(1,k)*term4 + uind(2,k)*term5 + durk*term6
                  term4 = -usr5 * zr
                  term6 = xr*zr*term2 - rr5*zr*urc5(1)
                  txzi = uind(1,i)*term4 + uind(3,i)*term5 + duri*term6
                  txzk = uind(1,k)*term4 + uind(3,k)*term5 + durk*term6
                  term5 = -yr*term1 + rr3*urc3(2)
                  term6 = yr*zr*term2 - rr5*zr*urc5(2)
                  tyzi = uind(2,i)*term4 + uind(3,i)*term5 + duri*term6
                  tyzk = uind(2,k)*term4 + uind(3,k)*term5 + durk*term6
                  depx = txxi*uinp(1,k) + txxk*uinp(1,i)
     &                      + txyi*uinp(2,k) + txyk*uinp(2,i)
     &                      + txzi*uinp(3,k) + txzk*uinp(3,i)
                  depy = txyi*uinp(1,k) + txyk*uinp(1,i)
     &                      + tyyi*uinp(2,k) + tyyk*uinp(2,i)
     &                      + tyzi*uinp(3,k) + tyzk*uinp(3,i)
                  depz = txzi*uinp(1,k) + txzk*uinp(1,i)
     &                      + tyzi*uinp(2,k) + tyzk*uinp(2,i)
     &                      + tzzi*uinp(3,k) + tzzk*uinp(3,i)
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
               end if
c
c     increment gradient and virial due to Cartesian forces
c
               depo(1,ii) = depo(1,ii) - frcx
               depo(2,ii) = depo(2,ii) - frcy
               depo(3,ii) = depo(3,ii) - frcz
               depo(1,kk) = depo(1,kk) + frcx
               depo(2,kk) = depo(2,kk) + frcy
               depo(3,kk) = depo(3,kk) + frcz
               vxx = xr * frcx
               vxy = yr * frcx
               vxz = zr * frcx
               vyy = yr * frcy
               vyz = zr * frcy
               vzz = zr * frcz
               viro(1,1) = viro(1,1) + vxx
               viro(2,1) = viro(2,1) + vxy
               viro(3,1) = viro(3,1) + vxz
               viro(1,2) = viro(1,2) + vxy
               viro(2,2) = viro(2,2) + vyy
               viro(3,2) = viro(3,2) + vyz
               viro(1,3) = viro(1,3) + vxz
               viro(2,3) = viro(2,3) + vyz
               viro(3,3) = viro(3,3) + vzz
c
c     get the induced dipole field used for dipole torques
c
               txi3 = psr3*uind(1,k) + dsr3*uinp(1,k)
               tyi3 = psr3*uind(2,k) + dsr3*uinp(2,k)
               tzi3 = psr3*uind(3,k) + dsr3*uinp(3,k)
               txk3 = psr3*uind(1,i) + dsr3*uinp(1,i)
               tyk3 = psr3*uind(2,i) + dsr3*uinp(2,i)
               tzk3 = psr3*uind(3,i) + dsr3*uinp(3,i)
               turi = -psr5*durk - dsr5*purk
               turk = -psr5*duri - dsr5*puri
               ufld(1,i) = ufld(1,i) + txi3 + xr*turi
               ufld(2,i) = ufld(2,i) + tyi3 + yr*turi
               ufld(3,i) = ufld(3,i) + tzi3 + zr*turi
               ufld(1,k) = ufld(1,k) + txk3 + xr*turk
               ufld(2,k) = ufld(2,k) + tyk3 + yr*turk
               ufld(3,k) = ufld(3,k) + tzk3 + zr*turk
c
c     get induced dipole field gradient used for quadrupole torques
c
               txi5 = 2.0d0 * (psr5*uind(1,k)+dsr5*uinp(1,k))
               tyi5 = 2.0d0 * (psr5*uind(2,k)+dsr5*uinp(2,k))
               tzi5 = 2.0d0 * (psr5*uind(3,k)+dsr5*uinp(3,k))
               txk5 = 2.0d0 * (psr5*uind(1,i)+dsr5*uinp(1,i))
               tyk5 = 2.0d0 * (psr5*uind(2,i)+dsr5*uinp(2,i))
               tzk5 = 2.0d0 * (psr5*uind(3,i)+dsr5*uinp(3,i))
               turi = -psr7*durk - dsr7*purk
               turk = -psr7*duri - dsr7*puri
               dufld(1,i) = dufld(1,i) + xr*txi5 + xr*xr*turi
               dufld(2,i) = dufld(2,i) + xr*tyi5 + yr*txi5
     &                         + 2.0d0*xr*yr*turi
               dufld(3,i) = dufld(3,i) + yr*tyi5 + yr*yr*turi
               dufld(4,i) = dufld(4,i) + xr*tzi5 + zr*txi5
     &                         + 2.0d0*xr*zr*turi
               dufld(5,i) = dufld(5,i) + yr*tzi5 + zr*tyi5
     &                         + 2.0d0*yr*zr*turi
               dufld(6,i) = dufld(6,i) + zr*tzi5 + zr*zr*turi
               dufld(1,k) = dufld(1,k) - xr*txk5 - xr*xr*turk
               dufld(2,k) = dufld(2,k) - xr*tyk5 - yr*txk5
     &                         - 2.0d0*xr*yr*turk
               dufld(3,k) = dufld(3,k) - yr*tyk5 - yr*yr*turk
               dufld(4,k) = dufld(4,k) - xr*tzk5 - zr*txk5
     &                         - 2.0d0*xr*zr*turk
               dufld(5,k) = dufld(5,k) - yr*tzk5 - zr*tyk5
     &                         - 2.0d0*yr*zr*turk
               dufld(6,k) = dufld(6,k) - zr*tzk5 - zr*zr*turk
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO reduction(+:depo,viro) schedule(guided)
c
c     torque is induced field and gradient cross permanent moments
c
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         trq(1) = diz*ufld(2,i) - diy*ufld(3,i)
     &               + qixz*dufld(2,i) - qixy*dufld(4,i)
     &               + 2.0d0*qiyz*(dufld(3,i)-dufld(6,i))
     &               + (qizz-qiyy)*dufld(5,i)
         trq(2) = dix*ufld(3,i) - diz*ufld(1,i)
     &               - qiyz*dufld(2,i) + qixy*dufld(5,i)
     &               + 2.0d0*qixz*(dufld(6,i)-dufld(1,i))
     &               + (qixx-qizz)*dufld(4,i)
         trq(3) = diy*ufld(1,i) - dix*ufld(2,i)
     &               + qiyz*dufld(4,i) - qixz*dufld(5,i)
     &               + 2.0d0*qixy*(dufld(1,i)-dufld(3,i))
     &               + (qiyy-qixx)*dufld(2,i)
         call torque (i,trq,fix,fiy,fiz,depo)
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
c     add local copies to global variables for OpenMP calculation
c
      ep = ep + epo
      do i = 1, n
         do j = 1, 3
            dep(j,i) = dep(j,i) + depo(j,i)
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
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (depo)
      deallocate (ufld)
      deallocate (dufld)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine eprecip1  --  polar Ewald recip energy & derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "eprecip1" evaluates the reciprocal space portion of the
c     particle mesh Ewald summation energy and gradient due to
c     dipole polarization
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
      subroutine eprecip1
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
      use polar
      use polpot
      use potent
      use virial
      implicit none
      integer i,j,k,ii
      integer j1,j2,j3
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
      real*8 cphim(4),cphid(4)
      real*8 cphip(4)
      real*8 a(3,3),ftc(10,10)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: fphid(:,:)
      real*8, allocatable :: fphip(:,:)
      real*8, allocatable :: fphidp(:,:)
      real*8, allocatable :: qgrip(:,:,:,:)
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
c     get the fractional to Cartesian transformation matrix
c
      call frac_to_cart (ftc)
c
c     initialize variables required for the scalar summation
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
c
c     remove scalar sum virial from prior multipole 3-D FFT
c
      if (use_mpole) then
         vxx = -vmxx
         vxy = -vmxy
         vxz = -vmxz
         vyy = -vmyy
         vyz = -vmyz
         vzz = -vmzz
c
c     compute the arrays of B-spline coefficients
c
      else
         call bspline_fill
         call table_fill
c
c     assign only the permanent multipoles to the PME grid
c     and perform the 3-D FFT forward transformation
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
         call cmp_to_fmp (cmp,fmp)
         call grid_mpole (fmp)
         call fftfront
c
c     make the scalar summation over reciprocal lattice
c
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
               vxx = vxx - h1*h1*vterm + eterm
               vxy = vxy - h1*h2*vterm
               vxz = vxz - h1*h3*vterm
               vyy = vyy - h2*h2*vterm + eterm
               vyz = vyz - h2*h3*vterm
               vzz = vzz - h3*h3*vterm + eterm
            end if
         end do
c
c     account for the zeroth grid point for a finite system
c
         qfac(1,1,1) = 0.0d0
         if (.not. use_bounds) then
            expterm = 0.5d0 * pi / xbox
            struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
c           e = 0.5d0 * expterm * struc2
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
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole))
      allocate (fuinp(3,npole))
      allocate (fphid(10,npole))
      allocate (fphip(10,npole))
      allocate (fphidp(20,npole))
c
c     convert Cartesian induced dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
      do i = 1, npole
         do j = 1, 3
            fuind(j,i) = a(j,1)*uind(1,i) + a(j,2)*uind(2,i)
     &                      + a(j,3)*uind(3,i)
            fuinp(j,i) = a(j,1)*uinp(1,i) + a(j,2)*uinp(2,i)
     &                      + a(j,3)*uinp(3,i)
         end do
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_uind (fuind,fuinp)
      call fftfront
c
c     account for the zeroth grid point for a finite system
c
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
c        e = 0.5d0 * expterm * struc2
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
      call fphi_uind (fphid,fphip,fphidp)
      do i = 1, npole
         do j = 1, 10
            fphid(j,i) = electric * fphid(j,i)
            fphip(j,i) = electric * fphip(j,i)
         end do
         do j = 1, 20
            fphidp(j,i) = electric * fphidp(j,i)
         end do
      end do
c
c     increment the induced dipole energy and gradient
c
c     e = 0.0d0
      do i = 1, npole
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 3
            j1 = deriv1(k+1)
            j2 = deriv2(k+1)
            j3 = deriv3(k+1)
c           e = e + fuind(k,i)*fphi(k+1,i)
            f1 = f1 + (fuind(k,i)+fuinp(k,i))*fphi(j1,i)
     &              + fuind(k,i)*fphip(j1,i)
     &              + fuinp(k,i)*fphid(j1,i)
            f2 = f2 + (fuind(k,i)+fuinp(k,i))*fphi(j2,i)
     &              + fuind(k,i)*fphip(j2,i)
     &              + fuinp(k,i)*fphid(j2,i)
            f3 = f3 + (fuind(k,i)+fuinp(k,i))*fphi(j3,i)
     &              + fuind(k,i)*fphip(j3,i)
     &              + fuinp(k,i)*fphid(j3,i)
            if (poltyp .eq. 'DIRECT') then
               f1 = f1 - fuind(k,i)*fphip(j1,i)
     &                 - fuinp(k,i)*fphid(j1,i)
               f2 = f2 - fuind(k,i)*fphip(j2,i)
     &                 - fuinp(k,i)*fphid(j2,i)
               f3 = f3 - fuind(k,i)*fphip(j3,i)
     &                 - fuinp(k,i)*fphid(j3,i)
            end if
         end do
         do k = 1, 10
            f1 = f1 + fmp(k,i)*fphidp(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphidp(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphidp(deriv3(k),i)
         end do
         f1 = 0.5d0 * dble(nfft1) * f1
         f2 = 0.5d0 * dble(nfft2) * f2
         f3 = 0.5d0 * dble(nfft3) * f3
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         ii = ipole(i)
         dep(1,ii) = dep(1,ii) + h1
         dep(2,ii) = dep(2,ii) + h2
         dep(3,ii) = dep(3,ii) + h3
      end do
c     e = 0.5d0 * e
c     ep = ep + e
c
c     set the potential to be the induced dipole average
c
      do i = 1, npole
         do k = 1, 10
            fphidp(k,i) = 0.5d0 * fphidp(k,i)
         end do
      end do
      call fphi_to_cphi (fphidp,cphi)
c
c     distribute torques into the induced dipole gradient
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
         call torque (i,trq,fix,fiy,fiz,dep)
      end do
c
c     induced dipole contribution to the internal virial
c
      do i = 1, npole
         do j = 2, 4
            cphim(j) = 0.0d0
            cphid(j) = 0.0d0
            cphip(j) = 0.0d0
            do k = 2, 4
               cphim(j) = cphim(j) + ftc(j,k)*fphi(k,i)
               cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
               cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
            end do
         end do
         vxx = vxx - cphi(2,i)*cmp(2,i)
     &             - 0.5d0*(cphim(2)*(uind(1,i)+uinp(1,i))
     &                     +cphid(2)*uinp(1,i)+cphip(2)*uind(1,i))
         vxy = vxy - 0.5d0*(cphi(2,i)*cmp(3,i)+cphi(3,i)*cmp(2,i))
     &             - 0.25d0*(cphim(2)*(uind(2,i)+uinp(2,i))
     &                      +cphim(3)*(uind(1,i)+uinp(1,i))
     &                      +cphid(2)*uinp(2,i)+cphip(2)*uind(2,i)
     &                      +cphid(3)*uinp(1,i)+cphip(3)*uind(1,i))
         vxz = vxz - 0.5d0*(cphi(2,i)*cmp(4,i)+cphi(4,i)*cmp(2,i))
     &             - 0.25d0*(cphim(2)*(uind(3,i)+uinp(3,i))
     &                      +cphim(4)*(uind(1,i)+uinp(1,i))
     &                      +cphid(2)*uinp(3,i)+cphip(2)*uind(3,i)
     &                      +cphid(4)*uinp(1,i)+cphip(4)*uind(1,i))
         vyy = vyy - cphi(3,i)*cmp(3,i)
     &             - 0.5d0*(cphim(3)*(uind(2,i)+uinp(2,i))
     &                     +cphid(3)*uinp(2,i)+cphip(3)*uind(2,i))
         vyz = vyz - 0.5d0*(cphi(3,i)*cmp(4,i)+cphi(4,i)*cmp(3,i))
     &             - 0.25d0*(cphim(3)*(uind(3,i)+uinp(3,i))
     &                      +cphim(4)*(uind(2,i)+uinp(2,i))
     &                      +cphid(3)*uinp(3,i)+cphip(3)*uind(3,i)
     &                      +cphid(4)*uinp(2,i)+cphip(4)*uind(2,i))
         vzz = vzz - cphi(4,i)*cmp(4,i)
     &             - 0.5d0*(cphim(4)*(uind(3,i)+uinp(3,i))
     &                     +cphid(4)*uinp(3,i)+cphip(4)*uind(3,i))
         vxx = vxx - 2.0d0*cmp(5,i)*cphi(5,i) - cmp(8,i)*cphi(8,i)
     &             - cmp(9,i)*cphi(9,i)
         vxy = vxy - (cmp(5,i)+cmp(6,i))*cphi(8,i)
     &             - 0.5d0*(cmp(8,i)*(cphi(6,i)+cphi(5,i))
     &                  +cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
         vxz = vxz - (cmp(5,i)+cmp(7,i))*cphi(9,i)
     &             - 0.5d0*(cmp(9,i)*(cphi(5,i)+cphi(7,i))
     &                  +cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
         vyy = vyy - 2.0d0*cmp(6,i)*cphi(6,i) - cmp(8,i)*cphi(8,i)
     &             - cmp(10,i)*cphi(10,i)
         vyz = vyz - (cmp(6,i)+cmp(7,i))*cphi(10,i)
     &             - 0.5d0*(cmp(10,i)*(cphi(6,i)+cphi(7,i))
     &                  +cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
         vzz = vzz - 2.0d0*cmp(7,i)*cphi(7,i) - cmp(9,i)*cphi(9,i)
     &             - cmp(10,i)*cphi(10,i)
         if (poltyp .eq. 'DIRECT') then
            vxx = vxx + 0.5d0*(cphid(2)*uinp(1,i)+cphip(2)*uind(1,i))
            vxy = vxy + 0.25d0*(cphid(2)*uinp(2,i)+cphip(2)*uind(2,i)
     &                        +cphid(3)*uinp(1,i)+cphip(3)*uind(1,i))
            vxz = vxz + 0.25d0*(cphid(2)*uinp(3,i)+cphip(2)*uind(3,i)
     &                        +cphid(4)*uinp(1,i)+cphip(4)*uind(1,i))
            vyy = vyy + 0.5d0*(cphid(3)*uinp(2,i)+cphip(3)*uind(2,i))
            vyz = vyz + 0.25d0*(cphid(3)*uinp(3,i)+cphip(3)*uind(3,i)
     &                        +cphid(4)*uinp(2,i)+cphip(4)*uind(2,i))
            vzz = vzz + 0.5d0*(cphid(4)*uinp(3,i)+cphip(4)*uind(3,i))
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (qgrip(2,nfft1,nfft2,nfft3))
c
c     assign permanent and induced multipoles to the PME grid
c     and perform the 3-D FFT forward transformation
c
      do i = 1, npole
         do j = 2, 4
            cmp(j,i) = cmp(j,i) + uinp(j-1,i)
         end do
      end do
      call cmp_to_fmp (cmp,fmp)
      call grid_mpole (fmp)
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
      call cmp_to_fmp (cmp,fmp)
      call grid_mpole (fmp)
      call fftfront
c
c     make the scalar summation over reciprocal lattice
c
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
            struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                  + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
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
c     assign only the induced dipoles to the PME grid
c     and perform the 3-D FFT forward transformation
c
      if (poltyp .eq. 'DIRECT') then
         do i = 1, npole
            do j = 1, 10
               cmp(j,i) = 0.0d0
            end do
            do j = 2, 4
               cmp(j,i) = uinp(j-1,i)
            end do
         end do
         call cmp_to_fmp (cmp,fmp)
         call grid_mpole (fmp)
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
               cmp(j,i) = uind(j-1,i)
            end do
         end do
         call cmp_to_fmp (cmp,fmp)
         call grid_mpole (fmp)
         call fftfront
c
c     make the scalar summation over reciprocal lattice
c
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
               struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                     + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
               eterm = 0.5d0 * electric * expterm * struc2
               vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
               vxx = vxx - h1*h1*vterm + eterm
               vxy = vxy - h1*h2*vterm
               vxz = vxz - h1*h3*vterm
               vyy = vyy - h2*h2*vterm + eterm
               vyz = vyz - h2*h3*vterm
               vzz = vzz - h3*h3*vterm + eterm
            end if
         end do
      end if
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
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fphid)
      deallocate (fphip)
      deallocate (fphidp)
      deallocate (qgrip)
      return
      end

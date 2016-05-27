c
c
c     ##################################################
c     ##  COPYRIGHT (C) 2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved             ##
c     ##################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine epolar1  --  induced dipole energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
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
c           call epolar1d
         else
c           call epolar1c
         end if
      else
         if (use_mlist) then
c           call epolar1b
         else
c           call epolar1a
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
      use chgpot
      use couple
      use deriv
      use energi
      use limits
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 f,damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale5
      real*8 scale7
      real*8 temp3,temp5,temp7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 usc3,usc5
      real*8 dsr3,dsr5,dsr7
      real*8 psr3,psr5,psr7
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
      real*8 term0,term1,term2
      real*8 term3,term4,term5
      real*8 term6,term7,term8
      real*8 depx,depy,depz
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
      real*8 rc3(3),drc3(3)
      real*8 prc3(3),urc3(3)
      real*8 rc5(3),drc5(3)
      real*8 prc5(3),urc5(3)
      real*8 rc7(3),drc7(3)
      real*8 prc7(3)
      real*8 trqi(3),frcxi(3)
      real*8 frcyi(3),frczi(3)
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
c     set arrays needed to scale connected atom interactions
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
c     set the permanent multipole and induced dipole values
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
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
c
c     set interaction scaling coefficients for connected atoms
c
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
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
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
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
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
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -damp * expdamp * rr5
                     temp5 = -damp * rr5 / rr3
                     temp7 = (-0.2d0 - 0.6d0*damp) * rr7 / rr5
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
               dsc3 = -scale3 * dscale(kk)
               dsc5 = -scale5 * dscale(kk)
               dsc7 = -scale7 * dscale(kk)
               psc3 = -scale3 * pscale(kk)
               psc5 = -scale5 * pscale(kk)
               psc7 = -scale7 * pscale(kk)
               usc3 = -scale3 * uscale(kk)
               usc5 = -scale5 * uscale(kk)
               do j = 1, 3
                  drc3(j) = rc3(j) * dscale(kk)
                  drc5(j) = rc5(j) * dscale(kk)
                  drc7(j) = rc7(j) * dscale(kk)
                  prc3(j) = rc3(j) * pscale(kk)
                  prc5(j) = rc5(j) * pscale(kk)
                  prc7(j) = rc7(j) * pscale(kk)
                  urc3(j) = rc3(j) * uscale(kk)
                  urc5(j) = rc5(j) * uscale(kk)
               end do
               dsr3 = rr3 * dsc3
               dsr5 = rr5 * dsc5
               dsr7 = rr7 * dsc7
               psr3 = rr3 * psc3
               psr5 = rr5 * psc5
               psr7 = rr7 * psc7
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
c     get the dEd/dR field gradient for direct polarization force
c
               term1 = dsc3*(rr3-rr5*xr*xr) - drc3(1)*xr
               term2 = (dsc3+dsc5)*rr5*xr + drc3(1)
               term3 = dsc5*(rr7*xr*xr-rr5) + drc5(1)*xr
               term4 = 2.0d0 * dsc5 * rr5
               term5 = 2.0d0 * (dsc5*rr7*xr+drc5(1)+1.5d0*dsc7*rr7*xr)
               term6 = xr * (dsc7*rr9*xr+drc7(1))
               txxi = ci*term1 + dix*term2 - dri*term3
     &                - qixx*term4 + qri1*term5 - qrri*term6
     &                + (qri2*yr+qri3*zr)*dsc7*rr7
               txxk = ck*term1 - dkx*term2 + drk*term3
     &                - qkxx*term4 + qrk1*term5 - qrrk*term6
     &                + (qrk2*yr+qrk3*zr)*dsc7*rr7
               term1 = dsc3*(rr3-rr5*yr*yr) - drc3(2)*yr
               term2 = (dsc3+dsc5)*rr5*yr + drc3(2)
               term3 = dsc5*(rr7*yr*yr-rr5) + drc5(2)*yr
               term4 = 2.0d0 * dsc5 * rr5
               term5 = 2.0d0 * (dsc5*rr7*yr+drc5(2)+1.5d0*dsc7*rr7*yr)
               term6 = yr * (dsc7*rr9*yr+drc7(2))
               tyyi = ci*term1 + diy*term2 - dri*term3
     &                - qiyy*term4 + qri2*term5 - qrri*term6
     &                + (qri1*xr+qri3*zr)*dsc7*rr7
               tyyk = ck*term1 - dky*term2 + drk*term3
     &                - qkyy*term4 + qrk2*term5 - qrrk*term6
     &                + (qrk1*xr+qrk3*zr)*dsc7*rr7
               term1 = dsc3*(rr3-rr5*zr*zr) - drc3(3)*zr
               term2 = (dsc3+dsc5)*rr5*zr + drc3(3)
               term3 = dsc5*(rr7*zr*zr-rr5) + drc5(3)*zr
               term4 = 2.0d0 * dsc5 * rr5
               term5 = 2.0d0 * (dsc5*rr7*zr+drc5(3)+1.5d0*dsc7*rr7*zr)
               term6 = zr * (dsc7*rr9*zr+drc7(3))
               tzzi = ci*term1 + diz*term2 - dri*term3
     &                - qizz*term4 + qri3*term5 - qrri*term6
     &                + (qri1*xr+qri2*yr)*dsc7*rr7
               tzzk = ck*term1 - dkz*term2 + drk*term3
     &                - qkzz*term4 + qrk3*term5 - qrrk*term6
     &                + (qrk1*xr+qrk2*yr)*dsc7*rr7
               term2 = dsc3*rr5*xr + drc3(1)
               term1 = yr * term2
               term3 = dsc5 * rr5 * yr
               term4 = yr * (dsc5*rr7*xr+drc5(1))
               term5 = 2.0d0 * dsc5 * rr5
               term6 = 2.0d0 * (dsc5*rr7*xr+drc5(1))
               term7 = 2.0d0 * dsc7 * rr7 * yr
               term8 = yr * (dsc7*rr9*xr+drc7(1))
               txyi = -ci*term1 + diy*term2 + dix*term3
     &                - dri*term4 - qixy*term5 + qri2*term6
     &                + qri1*term7 - qrri*term8
               txyk = -ck*term1 - dky*term2 - dkx*term3
     &                + drk*term4 - qkxy*term5 + qrk2*term6
     &                + qrk1*term7 - qrrk*term8
               term2 = dsc3*rr5*xr + drc3(1)
               term1 = zr * term2
               term3 = dsc5 * rr5 * zr
               term4 = zr * (dsc5*rr7*xr+drc5(1))
               term5 = 2.0d0 * dsc5 * rr5
               term6 = 2.0d0 * (dsc5*rr7*xr+drc5(1))
               term7 = 2.0d0 * dsc7 * rr7 * zr
               term8 = zr * (dsc7*rr9*xr+drc7(1))
               txzi = -ci*term1 + diz*term2 + dix*term3
     &                - dri*term4 - qixz*term5 + qri3*term6
     &                + qri1*term7 - qrri*term8
               txzk = -ck*term1 - dkz*term2 - dkx*term3
     &                + drk*term4 - qkxz*term5 + qrk3*term6
     &                + qrk1*term7 - qrrk*term8
               term2 = dsc3*rr5*yr + drc3(2)
               term1 = zr * term2
               term3 = dsc5 * rr5 * zr
               term4 = zr * (dsc5*rr7*yr+drc5(2))
               term5 = 2.0d0 * dsc5 * rr5
               term6 = 2.0d0 * (dsc5*rr7*yr+drc5(2))
               term7 = 2.0d0 * dsc7 * rr7 * zr
               term8 = zr * (dsc7*rr9*yr+drc7(2))
               tyzi = -ci*term1 + diz*term2 + diy*term3
     &                - dri*term4 - qiyz*term5 + qri3*term6
     &                + qri2*term7 - qrri*term8
               tyzk = -ck*term1 - dkz*term2 - dky*term3
     &                + drk*term4 - qkyz*term5 + qrk3*term6
     &                + qrk2*term7 - qrrk*term8
               depx = txxi*uinp(1,k) - txxk*uinp(1,i)
     &                   + txyi*uinp(2,k) - txyk*uinp(2,i)
     &                   + txzi*uinp(3,k) - txzk*uinp(3,i)
               depy = txyi*uinp(1,k) - txyk*uinp(1,i)
     &                   + tyyi*uinp(2,k) - tyyk*uinp(2,i)
     &                   + tyzi*uinp(3,k) - tyzk*uinp(3,i)
               depz =  txzi*uinp(1,k) - txzk*uinp(1,i)
     &                   + tyzi*uinp(2,k) - tyzk*uinp(2,i)
     &                   + tzzi*uinp(3,k) - tzzk*uinp(3,i)
               dep(1,ii) = dep(1,ii) - f*depx
               dep(2,ii) = dep(2,ii) - f*depy
               dep(3,ii) = dep(3,ii) - f*depz
               dep(1,kk) = dep(1,kk) + f*depx
               dep(2,kk) = dep(2,kk) + f*depy
               dep(3,kk) = dep(3,kk) + f*depz
c
c     get the dEp/dR field gradient for direct polarization force
c
               term1 = psc3*(rr3-rr5*xr*xr) - prc3(1)*xr
               term2 = (psc3+psc5)*rr5*xr + prc3(1)
               term3 = psc5*(rr7*xr*xr-rr5) + prc5(1)*xr
               term4 = 2.0d0 * psc5 * rr5
               term5 = 2.0d0 * (psc5*rr7*xr+prc5(1)+1.5d0*psc7*rr7*xr)
               term6 = xr * (psc7*rr9*xr+prc7(1))
               txxi = ci*term1 + dix*term2 - dri*term3
     &                - qixx*term4 + qri1*term5 - qrri*term6
     &                + (qri2*yr+qri3*zr)*psc7*rr7
               txxk = ck*term1 - dkx*term2 + drk*term3
     &                - qkxx*term4 + qrk1*term5 - qrrk*term6
     &                + (qrk2*yr+qrk3*zr)*psc7*rr7
               term1 = psc3*(rr3-rr5*yr*yr) - prc3(2)*yr
               term2 = (psc3+psc5)*rr5*yr + prc3(2)
               term3 = psc5*(rr7*yr*yr-rr5) + prc5(2)*yr
               term4 = 2.0d0 * psc5 * rr5
               term5 = 2.0d0 * (psc5*rr7*yr+prc5(2)+1.5d0*psc7*rr7*yr)
               term6 = yr * (psc7*rr9*yr+prc7(2))
               tyyi = ci*term1 + diy*term2 - dri*term3
     &                - qiyy*term4 + qri2*term5 - qrri*term6
     &                + (qri1*xr+qri3*zr)*psc7*rr7
               tyyk = ck*term1 - dky*term2 + drk*term3
     &                - qkyy*term4 + qrk2*term5 - qrrk*term6
     &                + (qrk1*xr+qrk3*zr)*psc7*rr7
               term1 = psc3*(rr3-rr5*zr*zr) - prc3(3)*zr
               term2 = (psc3+psc5)*rr5*zr + prc3(3)
               term3 = psc5*(rr7*zr*zr-rr5) + prc5(3)*zr
               term4 = 2.0d0 * psc5 * rr5
               term5 = 2.0d0 * (psc5*rr7*zr+prc5(3)+1.5d0*psc7*rr7*zr)
               term6 = zr * (psc7*rr9*zr+prc7(3))
               tzzi = ci*term1 + diz*term2 - dri*term3
     &                - qizz*term4 + qri3*term5 - qrri*term6
     &                + (qri1*xr+qri2*yr)*psc7*rr7
               tzzk = ck*term1 - dkz*term2 + drk*term3
     &                - qkzz*term4 + qrk3*term5 - qrrk*term6
     &                + (qrk1*xr+qrk2*yr)*psc7*rr7
               term2 = psc3*rr5*xr + prc3(1)
               term1 = yr * term2
               term3 = psc5 * rr5 * yr
               term4 = yr * (psc5*rr7*xr+prc5(1))
               term5 = 2.0d0 * psc5 * rr5
               term6 = 2.0d0 * (psc5*rr7*xr+prc5(1))
               term7 = 2.0d0 * psc7 * rr7 * yr
               term8 = yr * (psc7*rr9*xr+prc7(1))
               txyi = -ci*term1 + diy*term2 + dix*term3
     &                - dri*term4 - qixy*term5 + qri2*term6
     &                + qri1*term7 - qrri*term8
               txyk = -ck*term1 - dky*term2 - dkx*term3
     &                + drk*term4 - qkxy*term5 + qrk2*term6
     &                + qrk1*term7 - qrrk*term8
               term2 = psc3*rr5*xr + prc3(1)
               term1 = zr * term2
               term3 = psc5 * rr5 * zr
               term4 = zr * (psc5*rr7*xr+prc5(1))
               term5 = 2.0d0 * psc5 * rr5
               term6 = 2.0d0 * (psc5*rr7*xr+prc5(1))
               term7 = 2.0d0 * psc7 * rr7 * zr
               term8 = zr * (psc7*rr9*xr+prc7(1))
               txzi = -ci*term1 + diz*term2 + dix*term3
     &               - dri*term4 - qixz*term5 + qri3*term6
     &               + qri1*term7 - qrri*term8
               txzk = -ck*term1 - dkz*term2 - dkx*term3
     &               + drk*term4 - qkxz*term5 + qrk3*term6
     &               + qrk1*term7 - qrrk*term8
               term2 = psc3*rr5*yr + prc3(2)
               term1 = zr * term2
               term3 = psc5 * rr5 * zr
               term4 = zr * (psc5*rr7*yr+prc5(2))
               term5 = 2.0d0 * psc5 * rr5
               term6 = 2.0d0 * (psc5*rr7*yr+prc5(2))
               term7 = 2.0d0 * psc7 * rr7 * zr
               term8 = zr * (psc7*rr9*yr+prc7(2))
               tyzi = -ci*term1 + diz*term2 + diy*term3
     &               - dri*term4 - qiyz*term5 + qri3*term6
     &               + qri2*term7 - qrri*term8
               tyzk = -ck*term1 - dkz*term2 - dky*term3
     &               + drk*term4 - qkyz*term5 + qrk3*term6
     &               + qrk2*term7 - qrrk*term8
               depx = txxi*uind(1,k) - txxk*uind(1,i)
     &                   + txyi*uind(2,k) - txyk*uind(2,i)
     &                   + txzi*uind(3,k) - txzk*uind(3,i)
               depy = txyi*uind(1,k) - txyk*uind(1,i)
     &                   + tyyi*uind(2,k) - tyyk*uind(2,i)
     &                   + tyzi*uind(3,k) - tyzk*uind(3,i)
               depz =  txzi*uind(1,k) - txzk*uind(1,i)
     &                   + tyzi*uind(2,k) - tyzk*uind(2,i)
     &                   + tzzi*uind(3,k) - tzzk*uind(3,i)
               dep(1,ii) = dep(1,ii) - f*depx
               dep(2,ii) = dep(2,ii) - f*depy
               dep(3,ii) = dep(3,ii) - f*depz
               dep(1,kk) = dep(1,kk) + f*depx
               dep(2,kk) = dep(2,kk) + f*depy
               dep(3,kk) = dep(3,kk) + f*depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp .eq. 'MUTUAL') then
                  term0 = (usc3+usc5) * rr5
                  term1 = term0*xr + urc3(1)
                  term2 = usc5*(rr5-rr7*xr*xr) - urc5(1)*xr
                  txxi = uind(1,i)*term1 + duri*term2
                  txxk = uind(1,k)*term1 + durk*term2
                  term1 = term0*yr + urc3(2)
                  term2 = usc5*(rr5-rr7*yr*yr) - urc5(2)*yr
                  tyyi = uind(2,i)*term1 + duri*term2
                  tyyk = uind(2,k)*term1 + durk*term2
                  term1 = term0*zr + urc3(3)
                  term2 = usc5*(rr5-rr7*zr*zr) - urc5(3)*zr
                  tzzi = uind(3,i)*term1 + duri*term2
                  tzzk = uind(3,k)*term1 + durk*term2
                  term1 = usc5 * rr5 * yr
                  term2 = usc3*rr5*xr + urc3(1)
                  term3 = yr * (usc5*rr7*xr+urc5(1))
                  txyi = uind(1,i)*term1 + uind(2,i)*term2 - duri*term3
                  txyk = uind(1,k)*term1 + uind(2,k)*term2 - durk*term3
                  term1 = usc5 * rr5 * zr
                  term3 = zr * (usc5*rr7*xr+urc5(1))
                  txzi = uind(1,i)*term1 + uind(3,i)*term2 - duri*term3
                  txzk = uind(1,k)*term1 + uind(3,k)*term2 - durk*term3
                  term2 = usc3*rr5*yr + urc3(2)
                  term3 = zr * (usc5*rr7*yr+urc5(2))
                  tyzi = uind(2,i)*term1 + uind(3,i)*term2 - duri*term3
                  tyzk = uind(2,k)*term1 + uind(3,k)*term2 - durk*term3
                  depx = txxi*uinp(1,k) + txxk*uinp(1,i)
     &                      + txyi*uinp(2,k) + txyk*uinp(2,i)
     &                      + txzi*uinp(3,k) + txzk*uinp(3,i)
                  depy = txyi*uinp(1,k) + txyk*uinp(1,i)
     &                      + tyyi*uinp(2,k) + tyyk*uinp(2,i)
     &                      + tyzi*uinp(3,k) + tyzk*uinp(3,i)
                  depz =  txzi*uinp(1,k) + txzk*uinp(1,i)
     &                      + tyzi*uinp(2,k) + tyzk*uinp(2,i)
     &                      + tzzi*uinp(3,k) + tzzk*uinp(3,i)
                  dep(1,ii) = dep(1,ii) - f*depx
                  dep(2,ii) = dep(2,ii) - f*depy
                  dep(3,ii) = dep(3,ii) - f*depz
                  dep(1,kk) = dep(1,kk) + f*depx
                  dep(2,kk) = dep(2,kk) + f*depy
                  dep(3,kk) = dep(3,kk) + f*depz
               end if
c
c     get the induced dipole field used for dipole torques
c
               txi3 = psr3*uind(1,k) + dsr3*uinp(1,k)
               tyi3 = psr3*uind(2,k) + dsr3*uinp(2,k)
               tzi3 = psr3*uind(3,k) + dsr3*uinp(3,k)
               txk3 = psr3*uind(1,i) + dsr3*uinp(1,i)
               tyk3 = psr3*uind(2,i) + dsr3*uinp(2,i)
               tzk3 = psr3*uind(3,i) + dsr3*uinp(3,i)
               turi = -dsr5*purk - psr5*durk
               turk = -dsr5*puri - psr5*duri
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
               turi = -dsr7*purk - psr7*durk
               turk = -dsr7*puri - psr7*duri
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
c     reset interaction scaling coefficients for connected atoms
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
         trqi(1) = diy*ufld(3,i) - diz*ufld(2,i)
     &                + qixy*dufld(4,i) - qixz*dufld(2,i)
     &                + 2.0d0*qiyz*(dufld(6,i)-dufld(3,i))
     &                + (qiyy-qizz)*dufld(5,i)
         trqi(2) = diz*ufld(1,i) - dix*ufld(3,i)
     &                - qixy*dufld(5,i) + qiyz*dufld(2,i)
     &                + 2.0d0*qixz*(dufld(1,i)-dufld(6,i))
     &                + (qizz-qixx)*dufld(4,i)
         trqi(3) = dix*ufld(2,i) - diy*ufld(1,i)
     &                + qixz*dufld(5,i) - qiyz*dufld(4,i)
     &                + 2.0d0*qixy*(dufld(3,i)-dufld(1,i))
     &                + (qixx-qiyy)*dufld(2,i)
         trqi(1) = f * trqi(1)
         trqi(2) = f * trqi(2)
         trqi(3) = f * trqi(3)
         call torque4 (i,trqi,frcxi,frcyi,frczi,dep)
      end do
c
c     dipole polarization contribution to the internal virial
c
c     (below is the permanent multipoles with the induced dipole
c     field/gradient; full polarization virial also requires both
c     the permanent field/gradient and induced field/gradient with
c     the induced dipoles; see the original emrecip1 code)
c
      vxx = 0.0d0
      vyx = 0.0d0
      vzx = 0.0d0
      vyy = 0.0d0
      vzy = 0.0d0
      vzz = 0.0d0
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
         vxx = vxx + dix*ufld(1,i) + 2.0d0*qixx*dufld(1,i)
     &            + 2.0d0*qixy*dufld(2,i) + 2.0d0*qixz*dufld(4,i)
         vyx = vyx + 0.5d0*dix*ufld(2,i) + 0.5d0*diy*ufld(1,i)
     &            + (qixx+qiyy)*dufld(2,i)
     &            + qixy*(dufld(1,i)+dufld(3,i))
     &            + qixz*dufld(5,i) + qiyz*dufld(4,i)
         vzx = vzx + 0.5d0*dix*ufld(3,i) + 0.5d0*diz*ufld(1,i)
     &            + (qixx+qizz)*dufld(4,i)
     &            + qixz*(dufld(1,i)+dufld(6,i))
     &            + qixy*dufld(5,i) + qiyz*dufld(2,i)
         vyy = vyy + diy*ufld(2,i) + 2.0d0*qixy*dufld(2,i)
     &            + 2.0d0*qiyy*dufld(3,i) + 2.0d0*qiyz*dufld(5,i)
         vzy = vzy + 0.5d0*diy*ufld(3,i) + 0.5d0*diz*ufld(2,i)
     &            + (qiyy+qizz)*dufld(5,i)
     &            + qiyz*(dufld(3,i)+dufld(6,i))
     &            + qixy*dufld(4,i) + qixz*dufld(2,i)
         vzz = vzz + diz*ufld(3,i) + 2.0d0*qixz*dufld(4,i)
     &            + 2.0d0*qiyz*dufld(5,i) + 2.0d0*qizz*dufld(6,i)
      end do
      vir(1,1) = vir(1,1) + f*vxx
      vir(2,1) = vir(2,1) + f*vyx
      vir(3,1) = vir(3,1) + f*vzx
      vir(1,2) = vir(1,2) + f*vyx
      vir(2,2) = vir(2,2) + f*vyy
      vir(3,2) = vir(3,2) + f*vzy
      vir(1,3) = vir(1,3) + f*vzx
      vir(2,3) = vir(2,3) + f*vzy
      vir(3,3) = vir(3,3) + f*vzz
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

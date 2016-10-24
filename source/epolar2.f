c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar2  --  induced dipole polarization Hessian  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar2" calculates second derivatives of the dipole polarization
c     energy for a single atom at a time
c
c     it is incorrect to neglect interactions with atoms not directly
c     involved as the multipole site; to get better accuracy, "list"
c     should include all atoms by setting "biglist" to "true"
c
c     the "twosided" flag controls use of one-sided vs. two-sided
c     numerical derivatives; setting the flag to "true" gives a more
c     accurate Hessian at the expense of increased computation time
c
c     also, the "reinduce" flag controls whether the induced dipoles
c     are recomputed every time an atom is moved during computation
c     of the numerical Hessian; setting the flag to "true" produces a
c     much slower calculation, but aids convergence of minimizations,
c     accuracy of vibrational frequencies, etc.
c
c
      subroutine epolar2 (i)
      use sizes
      use atoms
      use deriv
      use hessn
      use limits
      use mpole
      implicit none
      integer i,j,k
      integer nlist
      integer, allocatable :: list(:)
      real*8 eps,old
      real*8, allocatable :: d0(:,:)
      logical prior
      logical biglist
      logical twosided
      logical reinduce
c
c
c     set the default stepsize and accuracy control flags
c
      eps = 1.0d-5
      biglist = .false.
      twosided = .false.
      reinduce = .false.
      if (n .le. 300) then
         biglist = .true.
         twosided = .true.
         reinduce = .true.
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
      if (allocated(dep)) then
         prior = .true.
         if (size(dep) .lt. 3*n)  deallocate (dep)
      end if
      if (.not. allocated(dep))  allocate (dep(3,n))
c
c     find the multipole definition involving the current atom;
c     results in a faster but approximate Hessian calculation
c
      nlist = 0
      do k = 1, npole
         if (biglist .or. ipole(k).eq.i) then
            nlist = nlist + 1
            list(nlist) = k
         end if
      end do
c
c     get multipole first derivatives for the base structure
c
      if (.not. twosided) then
         if (use_mlist) then
            call epolar2b (nlist,list,reinduce)
         else
            call epolar2a (nlist,list,reinduce)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dep(j,k)
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
            call epolar2b (nlist,list,reinduce)
         else
            call epolar2a (nlist,list,reinduce)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dep(j,k)
            end do
         end do
      end if
      x(i) = x(i) + eps
      if (use_mlist) then
         call epolar2b (nlist,list,reinduce)
      else
         call epolar2a (nlist,list,reinduce)
      end if
      x(i) = old
      do k = 1, n
         do j = 1, 3
            hessx(j,k) = hessx(j,k) + (dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical y-components via perturbed structures
c
      old = y(i)
      if (twosided) then
         y(i) = y(i) - 0.5d0*eps
         if (use_mlist) then
            call epolar2b (nlist,list,reinduce)
         else
            call epolar2a (nlist,list,reinduce)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dep(j,k)
            end do
         end do
      end if
      y(i) = y(i) + eps
      if (use_mlist) then
         call epolar2b (nlist,list,reinduce)
      else
         call epolar2a (nlist,list,reinduce)
      end if
      y(i) = old
      do k = 1, n
         do j = 1, 3
            hessy(j,k) = hessy(j,k) + (dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical z-components via perturbed structures
c
      old = z(i)
      if (twosided) then
         z(i) = z(i) - 0.5d0*eps
         if (use_mlist) then
            call epolar2b (nlist,list,reinduce)
         else
            call epolar2a (nlist,list,reinduce)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dep(j,k)
            end do
         end do
      end if
      z(i) = z(i) + eps
      if (use_mlist) then
         call epolar2b (nlist,list,reinduce)
      else
         call epolar2a (nlist,list,reinduce)
      end if
      z(i) = old
      do k = 1, n
         do j = 1, 3
            hessz(j,k) = hessz(j,k) + (dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     perform deallocation of some global arrays
c
      if (.not. prior) then
         deallocate (dep)
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar2a  --  polarization Hessian; numer, loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar2a" computes polarization first derivatives for a single
c     atom with respect to Cartesian coordinates; used to get finite
c     difference second derivatives
c
c
      subroutine epolar2a (nlist,list,reinduce)
      use sizes
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use deriv
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,iii
      integer iax,iay,iaz
      integer nlist,jcell
      integer list(*)
      real*8 f,fi,damp
      real*8 expdamp
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
      real*8 term0,term1,term2
      real*8 term3,term4,term5
      real*8 term6,term7,term8
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      logical reinduce
      character*6 mode
c
c
c     zero out the polarization derivative components
c
      do i = 1, n
         do j = 1, 3
            dep(j,i) = 0.0d0
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
c     compute the induced dipoles at each polarizable atom
c
      if (reinduce)  call induce
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
c     compute the dipole polarization gradient components
c
      do iii = 1, nlist
         i = list(iii)
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
               puri = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               purk = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               duri = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               durk = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
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
                  term0 = (sc3+sc5) * rr5
                  term1 = term0*xr - rc3(1)
                  term2 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  txxi = uind(1,i)*term1 + duri*term2
                  txxk = uind(1,k)*term1 + durk*term2
                  term1 = term0*yr - rc3(2)
                  term2 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tyyi = uind(2,i)*term1 + duri*term2
                  tyyk = uind(2,k)*term1 + durk*term2
                  term1 = term0*zr - rc3(3)
                  term2 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tzzi = uind(3,i)*term1 + duri*term2
                  tzzk = uind(3,k)*term1 + durk*term2
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
      do iii = 1, nlist
         i = list(iii)
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
               dscale(kk) = 1.0d0
               pscale(kk) = 1.0d0
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
               puri = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               purk = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               duri = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               durk = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
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
                  term0 = (sc3+sc5) * rr5
                  term1 = term0*xr - rc3(1)
                  term2 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  txxi = uind(1,i)*term1 + duri*term2
                  txxk = uind(1,k)*term1 + durk*term2
                  term1 = term0*yr - rc3(2)
                  term2 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tyyi = uind(2,i)*term1 + duri*term2
                  tyyk = uind(2,k)*term1 + durk*term2
                  term1 = term0*zr - rc3(3)
                  term2 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tzzi = uind(3,i)*term1 + duri*term2
                  tzzk = uind(3,k)*term1 + durk*term2
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
c     increment gradient components due to Cartesian forces
c
               dep(1,ii) = dep(1,ii) + frcx
               dep(2,ii) = dep(2,ii) + frcy
               dep(3,ii) = dep(3,ii) + frcz
               dep(1,kk) = dep(1,kk) - frcx
               dep(2,kk) = dep(2,kk) - frcy
               dep(3,kk) = dep(3,kk) - frcz
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
c     ##  subroutine epolar2b  --  polarization Hessian; numer, list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar2b" computes polarization first derivatives for a single
c     atom with respect to Cartesian coordinates; used to get finite
c     difference second derivatives
c
c
      subroutine epolar2b (nlist,list,reinduce)
      use sizes
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use limits
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk
      integer iii,kkk
      integer nlist
      integer list(*)
      real*8 f,fi,damp
      real*8 expdamp
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
      real*8 term0,term1,term2
      real*8 term3,term4,term5
      real*8 term6,term7,term8
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 trq(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: depo(:,:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      logical reinduce
      character*6 mode
c
c
c     zero out the polarization derivative components
c
      do i = 1, n
         do j = 1, 3
            dep(j,i) = 0.0d0
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
c     compute the induced dipoles at each polarizable atom
c
      if (reinduce)  call induce
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
      do i = 1, n
         do j = 1, 3
            depo(j,i) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(nlist,list,npole,ipole,pdamp,thole,x,y,z,rpole,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,
!$OMP& p2scale,p3scale,p4scale,p41scale,p5scale,d1scale,d2scale,d3scale,
!$OMP& d4scale,u1scale,u2scale,u3scale,u4scale,nelst,elst,use_bounds,f,
!$OMP& off2,uind,uinp,poltyp,depo,ufld,dufld)
!$OMP& firstprivate(pscale,dscale,uscale)
!$OMP DO reduction(+:depo,ufld,dufld) schedule(guided)
c
c     compute the dipole polarization gradient components
c
      do iii = 1, nlist
         i = list(iii)
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
               puri = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               purk = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               duri = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               durk = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
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
                  term0 = (sc3+sc5) * rr5
                  term1 = term0*xr - rc3(1)
                  term2 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  txxi = uind(1,i)*term1 + duri*term2
                  txxk = uind(1,k)*term1 + durk*term2
                  term1 = term0*yr - rc3(2)
                  term2 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tyyi = uind(2,i)*term1 + duri*term2
                  tyyk = uind(2,k)*term1 + durk*term2
                  term1 = term0*zr - rc3(3)
                  term2 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tzzi = uind(3,i)*term1 + duri*term2
                  tzzk = uind(3,k)*term1 + durk*term2
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
c     increment gradient components due to Cartesian forces
c
               depo(1,ii) = depo(1,ii) + frcx
               depo(2,ii) = depo(2,ii) + frcy
               depo(3,ii) = depo(3,ii) + frcz
               depo(1,kk) = depo(1,kk) - frcx
               depo(2,kk) = depo(2,kk) - frcy
               depo(3,kk) = depo(3,kk) - frcz
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
!$OMP DO reduction(+:depo) schedule(guided)
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
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
      do i = 1, n
         do j = 1, 3
            dep(j,i) = dep(j,i) + depo(j,i)
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

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
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine epolar1
      use iounit
      use limits
      use mplpot
      use polpot
      implicit none
c
c
c     check for use of TCG polarization with charge penetration
c
      if (poltyp.eq.'TCG' .and. use_chgpen) then
         write (iout,10)
   10    format (/,' EPOLAR1  --  TCG Polarization not Available',
     &              ' with Charge Penetration')
         call fatal
      end if
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
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use deriv
      use energi
      use molcul
      use mplpot
      use mpole
      use polar
      use polgrp
      use polopt
      use polpot
      use poltcg
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k,m
      integer ii,kk,jcell
      integer ix,iy,iz
      real*8 f,pgamma
      real*8 pdi,pti,ddi
      real*8 damp,expdamp
      real*8 temp3,temp5,temp7
      real*8 sc3,sc5,sc7
      real*8 sr3,sr5,sr7
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 dsr3i,dsr5i,dsr7i
      real*8 dsr3k,dsr5k,dsr7k
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 uixp,uiyp,uizp
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 ukxp,ukyp,ukzp
      real*8 dir,uir,uirp
      real*8 dkr,ukr,ukrp
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 uirm,ukrm
      real*8 uirt,ukrt
      real*8 tuir,tukr
      real*8 tixx,tiyy,tizz
      real*8 tixy,tixz,tiyz
      real*8 tkxx,tkyy,tkzz
      real*8 tkxy,tkxz,tkyz
      real*8 tix3,tiy3,tiz3
      real*8 tix5,tiy5,tiz5
      real*8 tkx3,tky3,tkz3
      real*8 tkx5,tky5,tkz5
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 term7,term8
      real*8 term1core
      real*8 term1i,term2i,term3i
      real*8 term4i,term5i,term6i
      real*8 term7i,term8i
      real*8 term1k,term2k,term3k
      real*8 term4k,term5k,term6k
      real*8 term7k,term8k
      real*8 poti,potk
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 tep(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8 uax(3),uay(3),uaz(3)
      real*8 ubx(3),uby(3),ubz(3)
      real*8 uaxp(3),uayp(3),uazp(3)
      real*8 ubxp(3),ubyp(3),ubzp(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
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
c     compute the total induced dipole polarization energy
c
      call epolar1e
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (wscale(n))
      allocate (ufld(3,n))
      allocate (dufld(6,n))
      allocate (pot(n))
      allocate (decfx(n))
      allocate (decfy(n))
      allocate (decfz(n))
c
c     set exclusion coefficients and arrays to store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
         wscale(i) = 1.0d0
         do j = 1, 3
            ufld(j,i) = 0.0d0
         end do
         do j = 1, 6
            dufld(j,i) = 0.0d0
         end do
         pot(i) = 0.0d0
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
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         uixp = uinp(1,ii)
         uiyp = uinp(2,ii)
         uizp = uinp(3,ii)
         do j = 1, tcgnab
            uax(j) = uad(1,ii,j)
            uay(j) = uad(2,ii,j)
            uaz(j) = uad(3,ii,j)
            uaxp(j) = uap(1,ii,j)
            uayp(j) = uap(2,ii,j)
            uazp(j) = uap(3,ii,j)
            ubx(j) = ubd(1,ii,j)
            uby(j) = ubd(2,ii,j)
            ubz(j) = ubd(3,ii,j)
            ubxp(j) = ubp(1,ii,j)
            ubyp(j) = ubp(2,ii,j)
            ubzp(j) = ubp(3,ii,j)
         end do
         if (use_thole) then
            pdi = pdamp(ii)
            pti = thole(ii)
            ddi = dirdamp(ii)
         else if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               dscale(i12(j,i)) = pscale(i12(j,i))
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = u4scale
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
               uscale(ip14(j,i)) = u4scale
            end do
         end if
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
               ukx = uind(1,kk)
               uky = uind(2,kk)
               ukz = uind(3,kk)
               ukxp = uinp(1,kk)
               ukyp = uinp(2,kk)
               ukzp = uinp(3,kk)
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
               uir = uix*xr + uiy*yr + uiz*zr
               uirp = uixp*xr + uiyp*yr + uizp*zr
               ukr = ukx*xr + uky*yr + ukz*zr
               ukrp = ukxp*xr + ukyp*yr + ukzp*zr
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     set initial values for tha damping scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
c
c     apply Thole polarization damping to scale factors
c
               if (use_thole) then
                  damp = pdi * pdamp(kk)
                  if (use_dirdamp) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(ddi,dirdamp(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           sc3 = 1.0d0 - expdamp 
                           sc5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                      +0.15d0*damp**2)
                           temp3 = 0.5d0 * damp * expdamp 
                           temp5 = 1.5d0 * (1.0d0+damp)
                           temp7 = 5.0d0*(1.5d0*damp*expdamp
     &                                *(0.35d0+0.35d0*damp
     &                                   +0.15d0*damp**2))/(temp3*temp5)
                           temp3 = temp3 * rr5
                           temp5 = temp5 / r2
                           temp7 = temp7 / r2
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
                  else
                     pgamma = min(pti,thole(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(pti,thole(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           sc3 = 1.0d0 - expdamp
                           sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                      +0.6d0*damp**2)
                           temp3 = damp * expdamp * rr5
                           temp5 = 3.0d0 * damp / r2
                           temp7 = (-1.0d0+3.0d0*damp) / r2
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
                  end if
                  sr3 = rr3 * sc3
                  sr5 = rr5 * sc5
                  sr7 = rr7 * sc7
                  dsr3 = sr3 * dscale(k)
                  dsr5 = sr5 * dscale(k)
                  dsr7 = sr7 * dscale(k)
                  psr3 = sr3 * pscale(k)
                  psr5 = sr5 * pscale(k)
                  psr7 = sr7 * pscale(k)
c
c     apply charge penetration damping to scale factors
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call damppole (r,9,alphai,alphak,dmpi,dmpk,dmpik)
                  dsr3i = 2.0d0 * rr3 * dmpi(3) * dscale(k)
                  dsr5i = 2.0d0 * rr5 * dmpi(5) * dscale(k)
                  dsr7i = 2.0d0 * rr7 * dmpi(7) * dscale(k)
                  dsr3k = 2.0d0 * rr3 * dmpk(3) * dscale(k)
                  dsr5k = 2.0d0 * rr5 * dmpk(5) * dscale(k)
                  dsr7k = 2.0d0 * rr7 * dmpk(7) * dscale(k)
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_thole) then
                     poti = -ukr*psr3 - ukrp*dsr3
                     potk = uir*psr3 + uirp*dsr3
                  else if (use_chgpen) then
                     poti = -ukr * dsr3i
                     potk = uir * dsr3k
                  end if
                  pot(i) = pot(i) + poti 
                  pot(k) = pot(k) + potk 
               end if
c
c     get the induced dipole field used for dipole torques
c
               if (use_thole) then
                  tix3 = psr3*ukx + dsr3*ukxp
                  tiy3 = psr3*uky + dsr3*ukyp
                  tiz3 = psr3*ukz + dsr3*ukzp
                  tkx3 = psr3*uix + dsr3*uixp
                  tky3 = psr3*uiy + dsr3*uiyp
                  tkz3 = psr3*uiz + dsr3*uizp
                  tuir = -psr5*ukr - dsr5*ukrp
                  tukr = -psr5*uir - dsr5*uirp
               else if (use_chgpen) then
                  tix3 = dsr3i*ukx
                  tiy3 = dsr3i*uky
                  tiz3 = dsr3i*ukz
                  tkx3 = dsr3k*uix
                  tky3 = dsr3k*uiy
                  tkz3 = dsr3k*uiz
                  tuir = -dsr5i*ukr
                  tukr = -dsr5k*uir
               end if
               ufld(1,i) = ufld(1,i) + tix3 + xr*tuir
               ufld(2,i) = ufld(2,i) + tiy3 + yr*tuir
               ufld(3,i) = ufld(3,i) + tiz3 + zr*tuir
               ufld(1,k) = ufld(1,k) + tkx3 + xr*tukr
               ufld(2,k) = ufld(2,k) + tky3 + yr*tukr
               ufld(3,k) = ufld(3,k) + tkz3 + zr*tukr
c
c     get induced dipole field gradient used for quadrupole torques
c
               if (use_thole) then
                  tix5 = 2.0d0 * (psr5*ukx+dsr5*ukxp)
                  tiy5 = 2.0d0 * (psr5*uky+dsr5*ukyp)
                  tiz5 = 2.0d0 * (psr5*ukz+dsr5*ukzp)
                  tkx5 = 2.0d0 * (psr5*uix+dsr5*uixp)
                  tky5 = 2.0d0 * (psr5*uiy+dsr5*uiyp)
                  tkz5 = 2.0d0 * (psr5*uiz+dsr5*uizp)
                  tuir = -psr7*ukr - dsr7*ukrp
                  tukr = -psr7*uir - dsr7*uirp
               else if (use_chgpen) then
                  tix5 = 2.0d0 * (dsr5i*ukx)
                  tiy5 = 2.0d0 * (dsr5i*uky)
                  tiz5 = 2.0d0 * (dsr5i*ukz)
                  tkx5 = 2.0d0 * (dsr5k*uix)
                  tky5 = 2.0d0 * (dsr5k*uiy)
                  tkz5 = 2.0d0 * (dsr5k*uiz)
                  tuir = -dsr7i*ukr
                  tukr = -dsr7k*uir
               end if
               dufld(1,i) = dufld(1,i) + xr*tix5 + xr*xr*tuir
               dufld(2,i) = dufld(2,i) + xr*tiy5 + yr*tix5
     &                         + 2.0d0*xr*yr*tuir
               dufld(3,i) = dufld(3,i) + yr*tiy5 + yr*yr*tuir
               dufld(4,i) = dufld(4,i) + xr*tiz5 + zr*tix5
     &                         + 2.0d0*xr*zr*tuir
               dufld(5,i) = dufld(5,i) + yr*tiz5 + zr*tiy5
     &                         + 2.0d0*yr*zr*tuir
               dufld(6,i) = dufld(6,i) + zr*tiz5 + zr*zr*tuir
               dufld(1,k) = dufld(1,k) - xr*tkx5 - xr*xr*tukr
               dufld(2,k) = dufld(2,k) - xr*tky5 - yr*tkx5
     &                         - 2.0d0*xr*yr*tukr
               dufld(3,k) = dufld(3,k) - yr*tky5 - yr*yr*tukr
               dufld(4,k) = dufld(4,k) - xr*tkz5 - zr*tkx5
     &                         - 2.0d0*xr*zr*tukr
               dufld(5,k) = dufld(5,k) - yr*tkz5 - zr*tky5
     &                         - 2.0d0*yr*zr*tukr
               dufld(6,k) = dufld(6,k) - zr*tkz5 - zr*zr*tukr
c
c     get the field gradient for direct polarization force
c
               if (use_thole) then
                  term1 = sc3*(rr3-rr5*xr*xr) + rc3(1)*xr
                  term2 = (sc3+sc5)*rr5*xr - rc3(1)
                  term3 = sc5*(rr7*xr*xr-rr5) - rc5(1)*xr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*xr-rc5(1)+1.5d0*sc7*rr7*xr)
                  term6 = xr * (sc7*rr9*xr-rc7(1))
                  tixx = ci*term1 + dix*term2 - dir*term3
     &                      - qixx*term4 + qix*term5 - qir*term6
     &                      + (qiy*yr+qiz*zr)*sc7*rr7
                  tkxx = ck*term1 - dkx*term2 + dkr*term3
     &                      - qkxx*term4 + qkx*term5 - qkr*term6
     &                      + (qky*yr+qkz*zr)*sc7*rr7
                  term1 = sc3*(rr3-rr5*yr*yr) + rc3(2)*yr
                  term2 = (sc3+sc5)*rr5*yr - rc3(2)
                  term3 = sc5*(rr7*yr*yr-rr5) - rc5(2)*yr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*yr-rc5(2)+1.5d0*sc7*rr7*yr)
                  term6 = yr * (sc7*rr9*yr-rc7(2))
                  tiyy = ci*term1 + diy*term2 - dir*term3
     &                      - qiyy*term4 + qiy*term5 - qir*term6
     &                      + (qix*xr+qiz*zr)*sc7*rr7
                  tkyy = ck*term1 - dky*term2 + dkr*term3
     &                      - qkyy*term4 + qky*term5 - qkr*term6
     &                      + (qkx*xr+qkz*zr)*sc7*rr7
                  term1 = sc3*(rr3-rr5*zr*zr) + rc3(3)*zr
                  term2 = (sc3+sc5)*rr5*zr - rc3(3)
                  term3 = sc5*(rr7*zr*zr-rr5) - rc5(3)*zr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*zr-rc5(3)+1.5d0*sc7*rr7*zr)
                  term6 = zr * (sc7*rr9*zr-rc7(3))
                  tizz = ci*term1 + diz*term2 - dir*term3
     &                      - qizz*term4 + qiz*term5 - qir*term6
     &                      + (qix*xr+qiy*yr)*sc7*rr7
                  tkzz = ck*term1 - dkz*term2 + dkr*term3
     &                      - qkzz*term4 + qkz*term5 - qkr*term6
     &                      + (qkx*xr+qky*yr)*sc7*rr7
                  term2 = sc3*rr5*xr - rc3(1)
                  term1 = yr * term2
                  term3 = sc5 * rr5 * yr
                  term4 = yr * (sc5*rr7*xr-rc5(1))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
                  term7 = 2.0d0 * sc7 * rr7 * yr
                  term8 = yr * (sc7*rr9*xr-rc7(1))
                  tixy = -ci*term1 + diy*term2 + dix*term3
     &                      - dir*term4 - qixy*term5 + qiy*term6
     &                      + qix*term7 - qir*term8
                  tkxy = -ck*term1 - dky*term2 - dkx*term3
     &                      + dkr*term4 - qkxy*term5 + qky*term6
     &                      + qkx*term7 - qkr*term8
                  term2 = sc3*rr5*xr - rc3(1)
                  term1 = zr * term2
                  term3 = sc5 * rr5 * zr
                  term4 = zr * (sc5*rr7*xr-rc5(1))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
                  term7 = 2.0d0 * sc7 * rr7 * zr
                  term8 = zr * (sc7*rr9*xr-rc7(1))
                  tixz = -ci*term1 + diz*term2 + dix*term3
     &                      - dir*term4 - qixz*term5 + qiz*term6
     &                      + qix*term7 - qir*term8
                  tkxz = -ck*term1 - dkz*term2 - dkx*term3
     &                      + dkr*term4 - qkxz*term5 + qkz*term6
     &                      + qkx*term7 - qkr*term8
                  term2 = sc3*rr5*yr - rc3(2)
                  term1 = zr * term2
                  term3 = sc5 * rr5 * zr
                  term4 = zr * (sc5*rr7*yr-rc5(2))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*yr-rc5(2))
                  term7 = 2.0d0 * sc7 * rr7 * zr
                  term8 = zr * (sc7*rr9*yr-rc7(2))
                  tiyz = -ci*term1 + diz*term2 + diy*term3
     &                      - dir*term4 - qiyz*term5 + qiz*term6
     &                      + qiy*term7 - qir*term8
                  tkyz = -ck*term1 - dkz*term2 - dky*term3
     &                      + dkr*term4 - qkyz*term5 + qkz*term6
     &                      + qky*term7 - qkr*term8
c
c     get the field gradient for direct polarization force
c
               else if (use_chgpen) then
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*xr*xr
                  term1core = rr3 - rr5*xr*xr
                  term2i = 2.0d0*rr5*dmpi(5)*xr 
                  term3i = rr7*dmpi(7)*xr*xr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*xr
                  term6i = rr9*dmpi(9)*xr*xr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*xr*xr
                  term2k = 2.0d0*rr5*dmpk(5)*xr
                  term3k = rr7*dmpk(7)*xr*xr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*xr
                  term6k = rr9*dmpk(9)*xr*xr
                  tixx = vali*term1i + corei*term1core  
     &                      + dix*term2i - dir*term3i
     &                      - qixx*term4i + qix*term5i - qir*term6i
     &                      + (qiy*yr+qiz*zr)*rr7*dmpi(7)
                  tkxx = valk*term1k + corek*term1core
     &                      - dkx*term2k + dkr*term3k
     &                      - qkxx*term4k + qkx*term5k - qkr*term6k
     &                      + (qky*yr+qkz*zr)*rr7*dmpk(7)
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*yr*yr
                  term1core = rr3 - rr5*yr*yr
                  term2i = 2.0d0*rr5*dmpi(5)*yr
                  term3i = rr7*dmpi(7)*yr*yr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*yr
                  term6i = rr9*dmpi(9)*yr*yr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*yr*yr
                  term2k = 2.0d0*rr5*dmpk(5)*yr
                  term3k = rr7*dmpk(7)*yr*yr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*yr
                  term6k = rr9*dmpk(9)*yr*yr
                  tiyy = vali*term1i + corei*term1core
     &                      + diy*term2i - dir*term3i
     &                      - qiyy*term4i + qiy*term5i - qir*term6i
     &                      + (qix*xr+qiz*zr)*rr7*dmpi(7)
                  tkyy = valk*term1k + corek*term1core
     &                      - dky*term2k + dkr*term3k
     &                      - qkyy*term4k + qky*term5k - qkr*term6k
     &                      + (qkx*xr+qkz*zr)*rr7*dmpk(7)
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*zr*zr
                  term1core = rr3 - rr5*zr*zr
                  term2i = 2.0d0*rr5*dmpi(5)*zr
                  term3i = rr7*dmpi(7)*zr*zr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*zr
                  term6i = rr9*dmpi(9)*zr*zr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*zr*zr
                  term2k = 2.0d0*rr5*dmpk(5)*zr
                  term3k = rr7*dmpk(7)*zr*zr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*zr
                  term6k = rr9*dmpk(9)*zr*zr
                  tizz = vali*term1i + corei*term1core
     &                      + diz*term2i - dir*term3i
     &                      - qizz*term4i + qiz*term5i - qir*term6i
     &                      + (qix*xr+qiy*yr)*rr7*dmpi(7)
                  tkzz = valk*term1k + corek*term1core
     &                      - dkz*term2k + dkr*term3k
     &                      - qkzz*term4k + qkz*term5k - qkr*term6k
     &                      + (qkx*xr+qky*yr)*rr7*dmpk(7)
                  term2i = rr5*dmpi(5)*xr 
                  term1i = yr * term2i
                  term1core = rr5*xr*yr
                  term3i = rr5*dmpi(5)*yr
                  term4i = yr * (rr7*dmpi(7)*xr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*xr
                  term7i = 2.0d0*rr7*dmpi(7)*yr
                  term8i = yr*rr9*dmpi(9)*xr
                  term2k = rr5*dmpk(5)*xr
                  term1k = yr * term2k
                  term3k = rr5*dmpk(5)*yr
                  term4k = yr * (rr7*dmpk(7)*xr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*xr
                  term7k = 2.0d0*rr7*dmpk(7)*yr
                  term8k = yr*rr9*dmpk(9)*xr
                  tixy = -vali*term1i - corei*term1core 
     &                      + diy*term2i + dix*term3i
     &                      - dir*term4i - qixy*term5i + qiy*term6i
     &                      + qix*term7i - qir*term8i
                  tkxy = -valk*term1k - corek*term1core 
     &                      - dky*term2k - dkx*term3k
     &                      + dkr*term4k - qkxy*term5k + qky*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5*dmpi(5)*xr
                  term1i = zr * term2i
                  term1core = rr5*xr*zr
                  term3i = rr5*dmpi(5)*zr
                  term4i = zr * (rr7*dmpi(7)*xr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*xr
                  term7i = 2.0d0*rr7*dmpi(7)*zr
                  term8i = zr*rr9*dmpi(9)*xr
                  term2k = rr5*dmpk(5)*xr
                  term1k = zr * term2k
                  term3k = rr5*dmpk(5)*zr
                  term4k = zr * (rr7*dmpk(7)*xr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*xr
                  term7k = 2.0d0*rr7*dmpk(7)*zr
                  term8k = zr*rr9*dmpk(9)*xr
                  tixz = -vali*term1i - corei*term1core
     &                      + diz*term2i + dix*term3i
     &                      - dir*term4i - qixz*term5i + qiz*term6i
     &                      + qix*term7i - qir*term8i
                  tkxz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dkx*term3k
     &                      + dkr*term4k - qkxz*term5k + qkz*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5*dmpi(5)*yr
                  term1i = zr * term2i
                  term1core = rr5*yr*zr
                  term3i = rr5*dmpi(5)*zr
                  term4i = zr * (rr7*dmpi(7)*yr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*yr
                  term7i = 2.0d0*rr7*dmpi(7)*zr
                  term8i = zr*rr9*dmpi(9)*yr
                  term2k = rr5*dmpk(5)*yr
                  term1k = zr * term2k
                  term3k = rr5*dmpk(5)*zr
                  term4k = zr * (rr7*dmpk(7)*yr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*yr
                  term7k = 2.0d0*rr7*dmpk(7)*zr
                  term8k = zr*rr9*dmpk(9)*yr
                  tiyz = -vali*term1i - corei*term1core
     &                      + diz*term2i + diy*term3i
     &                      - dir*term4i - qiyz*term5i + qiz*term6i
     &                      + qiy*term7i - qir*term8i
                  tkyz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dky*term3k
     &                      + dkr*term4k - qkyz*term5k + qkz*term6k
     &                      + qky*term7k - qkr*term8k
               end if
c
c     get the dEd/dR terms for Thole direct polarization force
c
               if (use_thole) then
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      - tkxx*uixp - tkxy*uiyp - tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      - tkxy*uixp - tkyy*uiyp - tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      - tkxz*uixp - tkyz*uiyp - tkzz*uizp
                  frcx = dscale(k) * depx
                  frcy = dscale(k) * depy
                  frcz = dscale(k) * depz
c
c     get the dEp/dR terms for Thole direct polarization force
c
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = frcx + pscale(k)*depx
                  frcy = frcy + pscale(k)*depy
                  frcz = frcz + pscale(k)*depz
c
c     get the dEp/dR terms for chgpen direct polarization force
c
               else if (use_chgpen) then
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = 2.0d0*dscale(k)*depx
                  frcy = 2.0d0*dscale(k)*depy
                  frcz = 2.0d0*dscale(k)*depz
               end if
c
c     reset Thole values if alternate direct damping was used
c
               if (use_dirdamp) then
                  sc3 = 1.0d0
                  sc5 = 1.0d0
                  do j = 1, 3
                     rc3(j) = 0.0d0
                     rc5(j) = 0.0d0
                  end do
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        sc3 = 1.0d0 - expdamp
                        sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                        temp3 = damp * expdamp * rr5
                        temp5 = 3.0d0 * damp / r2
                        rc3(1) = xr * temp3
                        rc3(2) = yr * temp3
                        rc3(3) = zr * temp3
                        rc5(1) = rc3(1) * temp5
                        rc5(2) = rc3(2) * temp5
                        rc5(3) = rc3(3) * temp5
                     end if
                  end if
               end if
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp.eq.'MUTUAL' .and. use_thole) then
                  term1 = (sc3+sc5) * rr5
                  term2 = term1*xr - rc3(1)
                  term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr - rc3(2)
                  term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr - rc3(3)
                  term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = sc5 * rr5 * yr
                  term2 = sc3*rr5*xr - rc3(1)
                  term3 = yr * (sc5*rr7*xr-rc5(1))
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = sc5 * rr5 * zr
                  term3 = zr * (sc5*rr7*xr-rc5(1))
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = sc3*rr5*yr - rc3(2)
                  term3 = zr * (sc5*rr7*yr-rc5(2))
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + uscale(k)*depx
                  frcy = frcy + uscale(k)*depy
                  frcz = frcz + uscale(k)*depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               else if (poltyp.eq.'MUTUAL' .and. use_chgpen) then
                  term1 = 2.0d0 * dmpik(5) * rr5
                  term2 = term1*xr
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*xr*xr 
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr 
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*yr*yr 
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr 
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*zr*zr 
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = rr5*dmpik(5)*yr
                  term2 = rr5*dmpik(5)*xr 
                  term3 = yr * (rr7*dmpik(7)*xr)
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = rr5 *dmpik(5) * zr
                  term3 = zr * (rr7*dmpik(7)*xr)
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = rr5*dmpik(5)*yr 
                  term3 = zr * (rr7*dmpik(7)*yr)
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + wscale(k)*depx
                  frcy = frcy + wscale(k)*depy
                  frcz = frcz + wscale(k)*depz
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_thole) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,ii)*xr + uopt(j,2,ii)*yr
     &                          + uopt(j,3,ii)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,kk)*xr + uopt(m,2,kk)*yr
     &                             + uopt(m,3,kk)*zr
                        term1 = (sc3+sc5) * rr5
                        term2 = term1*xr - rc3(1)
                        term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                        tixx = uopt(j,1,ii)*term2 + uirm*term3
                        tkxx = uopt(m,1,kk)*term2 + ukrm*term3
                        term2 = term1*yr - rc3(2)
                        term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                        tiyy = uopt(j,2,ii)*term2 + uirm*term3
                        tkyy = uopt(m,2,kk)*term2 + ukrm*term3
                        term2 = term1*zr - rc3(3)
                        term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                        tizz = uopt(j,3,ii)*term2 + uirm*term3
                        tkzz = uopt(m,3,kk)*term2 + ukrm*term3
                        term1 = sc5 * rr5 * yr
                        term2 = sc3*rr5*xr - rc3(1)
                        term3 = yr * (sc5*rr7*xr-rc5(1))
                        tixy = uopt(j,1,ii)*term1 + uopt(j,2,ii)*term2
     &                            - uirm*term3
                        tkxy = uopt(m,1,kk)*term1 + uopt(m,2,kk)*term2
     &                            - ukrm*term3
                        term1 = sc5 * rr5 * zr
                        term3 = zr * (sc5*rr7*xr-rc5(1))
                        tixz = uopt(j,1,ii)*term1 + uopt(j,3,ii)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,kk)*term1 + uopt(m,3,kk)*term2
     &                            - ukrm*term3
                        term2 = sc3*rr5*yr - rc3(2)
                        term3 = zr * (sc5*rr7*yr-rc5(2))
                        tiyz = uopt(j,2,ii)*term1 + uopt(j,3,ii)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,kk)*term1 + uopt(m,3,kk)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,kk) + tkxx*uoptp(j,1,ii)
     &                       + tixy*uoptp(m,2,kk) + tkxy*uoptp(j,2,ii)
     &                       + tixz*uoptp(m,3,kk) + tkxz*uoptp(j,3,ii)
                        depy = tixy*uoptp(m,1,kk) + tkxy*uoptp(j,1,ii)
     &                       + tiyy*uoptp(m,2,kk) + tkyy*uoptp(j,2,ii)
     &                       + tiyz*uoptp(m,3,kk) + tkyz*uoptp(j,3,ii)
                        depz = tixz*uoptp(m,1,kk) + tkxz*uoptp(j,1,ii)
     &                       + tiyz*uoptp(m,2,kk) + tkyz*uoptp(j,2,ii)
     &                       + tizz*uoptp(m,3,kk) + tkzz*uoptp(j,3,ii)
                        frcx = frcx + copm(j+m+1)*uscale(k)*depx
                        frcy = frcy + copm(j+m+1)*uscale(k)*depy
                        frcz = frcz + copm(j+m+1)*uscale(k)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_chgpen) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,i)*xr + uopt(j,2,i)*yr
     &                          + uopt(j,3,i)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,k)*xr + uopt(m,2,k)*yr
     &                             + uopt(m,3,k)*zr
                        term1 = 2.0d0 * dmpik(5) * rr5
                        term2 = term1*xr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*xr*xr
                        tixx = uopt(j,1,i)*term2 + uirm*term3
                        tkxx = uopt(m,1,k)*term2 + ukrm*term3
                        term2 = term1*yr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*yr*yr
                        tiyy = uopt(j,2,i)*term2 + uirm*term3
                        tkyy = uopt(m,2,k)*term2 + ukrm*term3
                        term2 = term1*zr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*zr*zr
                        tizz = uopt(j,3,i)*term2 + uirm*term3
                        tkzz = uopt(m,3,k)*term2 + ukrm*term3
                        term1 = rr5*dmpik(5)*yr
                        term2 = rr5*dmpik(5)*xr
                        term3 = yr * (rr7*dmpik(7)*xr)
                        tixy = uopt(j,1,i)*term1 + uopt(j,2,i)*term2
     &                            - uirm*term3
                        tkxy = uopt(m,1,k)*term1 + uopt(m,2,k)*term2
     &                            - ukrm*term3
                        term1 = rr5 *dmpik(5) * zr
                        term3 = zr * (rr7*dmpik(7)*xr)
                        tixz = uopt(j,1,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        term2 = rr5*dmpik(5)*yr
                        term3 = zr * (rr7*dmpik(7)*yr)
                        tiyz = uopt(j,2,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,k) + tkxx*uoptp(j,1,i)
     &                       + tixy*uoptp(m,2,k) + tkxy*uoptp(j,2,i)
     &                       + tixz*uoptp(m,3,k) + tkxz*uoptp(j,3,i)
                        depy = tixy*uoptp(m,1,k) + tkxy*uoptp(j,1,i)
     &                       + tiyy*uoptp(m,2,k) + tkyy*uoptp(j,2,i)
     &                       + tiyz*uoptp(m,3,k) + tkyz*uoptp(j,3,i)
                        depz = tixz*uoptp(m,1,k) + tkxz*uoptp(j,1,i)
     &                       + tiyz*uoptp(m,2,k) + tkyz*uoptp(j,2,i)
     &                       + tizz*uoptp(m,3,k) + tkzz*uoptp(j,3,i)
                        frcx = frcx + copm(j+m+1)*wscale(k)*depx
                        frcy = frcy + copm(j+m+1)*wscale(k)*depy
                        frcz = frcz + copm(j+m+1)*wscale(k)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for TCG polarization force
c
               else if (poltyp.eq.'TCG' .and. use_thole) then
                  do j = 1, tcgnab
                     ukx = ubd(1,kk,j)
                     uky = ubd(2,kk,j)
                     ukz = ubd(3,kk,j)
                     ukxp = ubp(1,kk,j)
                     ukyp = ubp(2,kk,j)
                     ukzp = ubp(3,kk,j)
                     uirt = uax(j)*xr + uay(j)*yr + uaz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = (sc3+sc5) * rr5
                     term2 = term1*xr - rc3(1)
                     term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                     tixx = uax(j)*term2 + uirt*term3
                     tkxx = ukx*term2 + ukrt*term3
                     term2 = term1*yr - rc3(2)
                     term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                     tiyy = uay(j)*term2 + uirt*term3
                     tkyy = uky*term2 + ukrt*term3
                     term2 = term1*zr - rc3(3)
                     term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                     tizz = uaz(j)*term2 + uirt*term3
                     tkzz = ukz*term2 + ukrt*term3
                     term1 = sc5 * rr5 * yr
                     term2 = sc3*rr5*xr - rc3(1)
                     term3 = yr * (sc5*rr7*xr-rc5(1))
                     tixy = uax(j)*term1 + uay(j)*term2 - uirt*term3
                     tkxy = ukx*term1 + uky*term2 - ukrt*term3
                     term1 = sc5 * rr5 * zr
                     term3 = zr * (sc5*rr7*xr-rc5(1))
                     tixz = uax(j)*term1 + uaz(j)*term2 - uirt*term3
                     tkxz = ukx*term1 + ukz*term2 - ukrt*term3
                     term2 = sc3*rr5*yr - rc3(2)
                     term3 = zr * (sc5*rr7*yr-rc5(2))
                     tiyz = uay(j)*term1 + uaz(j)*term2 - uirt*term3
                     tkyz = uky*term1 + ukz*term2 - ukrt*term3
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*uaxp(j) + tkxy*uayp(j)
     &                         + tkxz*uazp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*uaxp(j) + tkyy*uayp(j)
     &                         + tkyz*uazp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*uaxp(j) + tkyz*uayp(j)
     &                         + tkzz*uazp(j)
                     frcx = frcx + uscale(k)*depx
                     frcy = frcy + uscale(k)*depy
                     frcz = frcz + uscale(k)*depz
                     ukx = uad(1,kk,j)
                     uky = uad(2,kk,j)
                     ukz = uad(3,kk,j)
                     ukxp = uap(1,kk,j)
                     ukyp = uap(2,kk,j)
                     ukzp = uap(3,kk,j)
                     uirt = ubx(j)*xr + uby(j)*yr + ubz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = (sc3+sc5) * rr5
                     term2 = term1*xr - rc3(1)
                     term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                     tixx = ubx(j)*term2 + uirt*term3
                     tkxx = ukx*term2 + ukrt*term3
                     term2 = term1*yr - rc3(2)
                     term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                     tiyy = uby(j)*term2 + uirt*term3
                     tkyy = uky*term2 + ukrt*term3
                     term2 = term1*zr - rc3(3)
                     term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                     tizz = ubz(j)*term2 + uirt*term3
                     tkzz = ukz*term2 + ukrt*term3
                     term1 = sc5 * rr5 * yr
                     term2 = sc3*rr5*xr - rc3(1)
                     term3 = yr * (sc5*rr7*xr-rc5(1))
                     tixy = ubx(j)*term1 + uby(j)*term2 - uirt*term3
                     tkxy = ukx*term1 + uky*term2 - ukrt*term3
                     term1 = sc5 * rr5 * zr
                     term3 = zr * (sc5*rr7*xr-rc5(1))
                     tixz = ubx(j)*term1 + ubz(j)*term2 - uirt*term3
                     tkxz = ukx*term1 + ukz*term2 - ukrt*term3
                     term2 = sc3*rr5*yr - rc3(2)
                     term3 = zr * (sc5*rr7*yr-rc5(2))
                     tiyz = uby(j)*term1 + ubz(j)*term2 - uirt*term3
                     tkyz = uky*term1 + ukz*term2 - ukrt*term3
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*ubxp(j) + tkxy*ubyp(j)
     &                         + tkxz*ubzp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*ubxp(j) + tkyy*ubyp(j)
     &                         + tkyz*ubzp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*ubxp(j) + tkyz*ubyp(j)
     &                         + tkzz*ubzp(j)
                     frcx = frcx + uscale(k)*depx
                     frcy = frcy + uscale(k)*depy
                     frcz = frcz + uscale(k)*depz
                  end do
               end if
c
c     increment force-based gradient on the interaction sites
c
               dep(1,i) = dep(1,i) + frcx
               dep(2,i) = dep(2,i) + frcy
               dep(3,i) = dep(3,i) + frcz
               dep(1,k) = dep(1,k) - frcx
               dep(2,k) = dep(2,k) - frcy
               dep(3,k) = dep(3,k) - frcz
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
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = 1.0d0
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
               uscale(ip14(j,i)) = 1.0d0
            end do
         end if
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
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         uixp = uinp(1,ii)
         uiyp = uinp(2,ii)
         uizp = uinp(3,ii)
         do j = 1, tcgnab
            uax(j) = uad(1,ii,j)
            uay(j) = uad(2,ii,j)
            uaz(j) = uad(3,ii,j)
            uaxp(j) = uap(1,ii,j)
            uayp(j) = uap(2,ii,j)
            uazp(j) = uap(3,ii,j)
            ubx(j) = ubd(1,ii,j)
            uby(j) = ubd(2,ii,j)
            ubz(j) = ubd(3,ii,j)
            ubxp(j) = ubp(1,ii,j)
            ubyp(j) = ubp(2,ii,j)
            ubzp(j) = ubp(3,ii,j)
         end do
         if (use_thole) then
            pdi = pdamp(ii)
            pti = thole(ii)
            ddi = dirdamp(ii)
         else if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               dscale(i12(j,i)) = pscale(i12(j,i))
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = u4scale
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
               uscale(ip14(j,i)) = u4scale
            end do
         end if
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii, npole
            k = ipole(kk)
            do jcell = 2, ncell
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
            if (use_bounds)  call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               pscale(k) = 1.0d0
               dscale(k) = 1.0d0
               uscale(k) = 1.0d0
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
               ukx = uind(1,kk)
               uky = uind(2,kk)
               ukz = uind(3,kk)
               ukxp = uinp(1,kk)
               ukyp = uinp(2,kk)
               ukzp = uinp(3,kk)
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
               uir = uix*xr + uiy*yr + uiz*zr
               uirp = uixp*xr + uiyp*yr + uizp*zr
               ukr = ukx*xr + uky*yr + ukz*zr
               ukrp = ukxp*xr + ukyp*yr + ukzp*zr
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
c
c     apply Thole polarization damping to scale factors
c
               if (use_thole) then
                  damp = pdi * pdamp(kk)
                  if (use_dirdamp) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(ddi,dirdamp(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           sc3 = 1.0d0 - expdamp 
                           sc5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                      +0.15d0*damp**2)
                           temp3 = 0.5d0 * damp * expdamp 
                           temp5 = 1.5d0 * (1.0d0+damp)
                           temp7 = 5.0d0*(1.5d0*damp*expdamp
     &                                *(0.35d0+0.35d0*damp
     &                                   +0.15d0*damp**2))/(temp3*temp5)
                           temp3 = temp3 * rr5
                           temp5 = temp5 / r2
                           temp7 = temp7 / r2
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
                  else
                     pgamma = min(pti,thole(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(pti,thole(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           sc3 = 1.0d0 - expdamp
                           sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                      +0.6d0*damp**2)
                           temp3 = damp * expdamp * rr5
                           temp5 = 3.0d0 * damp / r2
                           temp7 = (-1.0d0+3.0d0*damp) / r2
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
                  end if
                  sr3 = rr3 * sc3
                  sr5 = rr5 * sc5
                  sr7 = rr7 * sc7
                  dsr3 = sr3 * dscale(k)
                  dsr5 = sr5 * dscale(k)
                  dsr7 = sr7 * dscale(k)
                  psr3 = sr3 * pscale(k)
                  psr5 = sr5 * pscale(k)
                  psr7 = sr7 * pscale(k)
c
c     apply charge penetration damping to scale factors
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call damppole (r,9,alphai,alphak,dmpi,dmpk,dmpik)
                  dsr3i = 2.0d0 * rr3 * dmpi(3) * dscale(k)
                  dsr5i = 2.0d0 * rr5 * dmpi(5) * dscale(k)
                  dsr7i = 2.0d0 * rr7 * dmpi(7) * dscale(k)
                  dsr3k = 2.0d0 * rr3 * dmpk(3) * dscale(k)
                  dsr5k = 2.0d0 * rr5 * dmpk(5) * dscale(k)
                  dsr7k = 2.0d0 * rr7 * dmpk(7) * dscale(k)
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_thole) then
                     poti = -ukr*psr3 - ukrp*dsr3
                     potk = uir*psr3 + uirp*dsr3
                  else if (use_chgpen) then
                     poti = -ukr * dsr3i
                     potk = uir * dsr3k
                  end if
                  pot(i) = pot(i) + poti 
                  pot(k) = pot(k) + potk 
               end if
c
c     get the induced dipole field used for dipole torques
c
               if (use_thole) then
                  tix3 = psr3*ukx + dsr3*ukxp
                  tiy3 = psr3*uky + dsr3*ukyp
                  tiz3 = psr3*ukz + dsr3*ukzp
                  tkx3 = psr3*uix + dsr3*uixp
                  tky3 = psr3*uiy + dsr3*uiyp
                  tkz3 = psr3*uiz + dsr3*uizp
                  tuir = -psr5*ukr - dsr5*ukrp
                  tukr = -psr5*uir - dsr5*uirp
               else if (use_chgpen) then
                  tix3 = dsr3i*ukx
                  tiy3 = dsr3i*uky
                  tiz3 = dsr3i*ukz
                  tkx3 = dsr3k*uix
                  tky3 = dsr3k*uiy
                  tkz3 = dsr3k*uiz
                  tuir = -dsr5i*ukr
                  tukr = -dsr5k*uir
               end if
               ufld(1,i) = ufld(1,i) + tix3 + xr*tuir
               ufld(2,i) = ufld(2,i) + tiy3 + yr*tuir
               ufld(3,i) = ufld(3,i) + tiz3 + zr*tuir
               ufld(1,k) = ufld(1,k) + tkx3 + xr*tukr
               ufld(2,k) = ufld(2,k) + tky3 + yr*tukr
               ufld(3,k) = ufld(3,k) + tkz3 + zr*tukr
c
c     get induced dipole field gradient used for quadrupole torques
c
               if (use_thole) then
                  tix5 = 2.0d0 * (psr5*ukx+dsr5*ukxp)
                  tiy5 = 2.0d0 * (psr5*uky+dsr5*ukyp)
                  tiz5 = 2.0d0 * (psr5*ukz+dsr5*ukzp)
                  tkx5 = 2.0d0 * (psr5*uix+dsr5*uixp)
                  tky5 = 2.0d0 * (psr5*uiy+dsr5*uiyp)
                  tkz5 = 2.0d0 * (psr5*uiz+dsr5*uizp)
                  tuir = -psr7*ukr - dsr7*ukrp
                  tukr = -psr7*uir - dsr7*uirp
               else if (use_chgpen) then
                  tix5 = 2.0d0 * (dsr5i*ukx)
                  tiy5 = 2.0d0 * (dsr5i*uky)
                  tiz5 = 2.0d0 * (dsr5i*ukz)
                  tkx5 = 2.0d0 * (dsr5k*uix)
                  tky5 = 2.0d0 * (dsr5k*uiy)
                  tkz5 = 2.0d0 * (dsr5k*uiz)
                  tuir = -dsr7i*ukr
                  tukr = -dsr7k*uir
               end if
               dufld(1,i) = dufld(1,i) + xr*tix5 + xr*xr*tuir
               dufld(2,i) = dufld(2,i) + xr*tiy5 + yr*tix5
     &                         + 2.0d0*xr*yr*tuir
               dufld(3,i) = dufld(3,i) + yr*tiy5 + yr*yr*tuir
               dufld(4,i) = dufld(4,i) + xr*tiz5 + zr*tix5
     &                         + 2.0d0*xr*zr*tuir
               dufld(5,i) = dufld(5,i) + yr*tiz5 + zr*tiy5
     &                         + 2.0d0*yr*zr*tuir
               dufld(6,i) = dufld(6,i) + zr*tiz5 + zr*zr*tuir
               dufld(1,k) = dufld(1,k) - xr*tkx5 - xr*xr*tukr
               dufld(2,k) = dufld(2,k) - xr*tky5 - yr*tkx5
     &                         - 2.0d0*xr*yr*tukr
               dufld(3,k) = dufld(3,k) - yr*tky5 - yr*yr*tukr
               dufld(4,k) = dufld(4,k) - xr*tkz5 - zr*tkx5
     &                         - 2.0d0*xr*zr*tukr
               dufld(5,k) = dufld(5,k) - yr*tkz5 - zr*tky5
     &                         - 2.0d0*yr*zr*tukr
               dufld(6,k) = dufld(6,k) - zr*tkz5 - zr*zr*tukr
c
c     get the field gradient for direct polarization force
c
               if (use_thole) then
                  term1 = sc3*(rr3-rr5*xr*xr) + rc3(1)*xr
                  term2 = (sc3+sc5)*rr5*xr - rc3(1)
                  term3 = sc5*(rr7*xr*xr-rr5) - rc5(1)*xr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*xr-rc5(1)+1.5d0*sc7*rr7*xr)
                  term6 = xr * (sc7*rr9*xr-rc7(1))
                  tixx = ci*term1 + dix*term2 - dir*term3
     &                      - qixx*term4 + qix*term5 - qir*term6
     &                      + (qiy*yr+qiz*zr)*sc7*rr7
                  tkxx = ck*term1 - dkx*term2 + dkr*term3
     &                      - qkxx*term4 + qkx*term5 - qkr*term6
     &                      + (qky*yr+qkz*zr)*sc7*rr7
                  term1 = sc3*(rr3-rr5*yr*yr) + rc3(2)*yr
                  term2 = (sc3+sc5)*rr5*yr - rc3(2)
                  term3 = sc5*(rr7*yr*yr-rr5) - rc5(2)*yr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*yr-rc5(2)+1.5d0*sc7*rr7*yr)
                  term6 = yr * (sc7*rr9*yr-rc7(2))
                  tiyy = ci*term1 + diy*term2 - dir*term3
     &                      - qiyy*term4 + qiy*term5 - qir*term6
     &                      + (qix*xr+qiz*zr)*sc7*rr7
                  tkyy = ck*term1 - dky*term2 + dkr*term3
     &                      - qkyy*term4 + qky*term5 - qkr*term6
     &                      + (qkx*xr+qkz*zr)*sc7*rr7
                  term1 = sc3*(rr3-rr5*zr*zr) + rc3(3)*zr
                  term2 = (sc3+sc5)*rr5*zr - rc3(3)
                  term3 = sc5*(rr7*zr*zr-rr5) - rc5(3)*zr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*zr-rc5(3)+1.5d0*sc7*rr7*zr)
                  term6 = zr * (sc7*rr9*zr-rc7(3))
                  tizz = ci*term1 + diz*term2 - dir*term3
     &                      - qizz*term4 + qiz*term5 - qir*term6
     &                      + (qix*xr+qiy*yr)*sc7*rr7
                  tkzz = ck*term1 - dkz*term2 + dkr*term3
     &                      - qkzz*term4 + qkz*term5 - qkr*term6
     &                      + (qkx*xr+qky*yr)*sc7*rr7
                  term2 = sc3*rr5*xr - rc3(1)
                  term1 = yr * term2
                  term3 = sc5 * rr5 * yr
                  term4 = yr * (sc5*rr7*xr-rc5(1))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
                  term7 = 2.0d0 * sc7 * rr7 * yr
                  term8 = yr * (sc7*rr9*xr-rc7(1))
                  tixy = -ci*term1 + diy*term2 + dix*term3
     &                      - dir*term4 - qixy*term5 + qiy*term6
     &                      + qix*term7 - qir*term8
                  tkxy = -ck*term1 - dky*term2 - dkx*term3
     &                      + dkr*term4 - qkxy*term5 + qky*term6
     &                      + qkx*term7 - qkr*term8
                  term2 = sc3*rr5*xr - rc3(1)
                  term1 = zr * term2
                  term3 = sc5 * rr5 * zr
                  term4 = zr * (sc5*rr7*xr-rc5(1))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
                  term7 = 2.0d0 * sc7 * rr7 * zr
                  term8 = zr * (sc7*rr9*xr-rc7(1))
                  tixz = -ci*term1 + diz*term2 + dix*term3
     &                      - dir*term4 - qixz*term5 + qiz*term6
     &                      + qix*term7 - qir*term8
                  tkxz = -ck*term1 - dkz*term2 - dkx*term3
     &                      + dkr*term4 - qkxz*term5 + qkz*term6
     &                      + qkx*term7 - qkr*term8
                  term2 = sc3*rr5*yr - rc3(2)
                  term1 = zr * term2
                  term3 = sc5 * rr5 * zr
                  term4 = zr * (sc5*rr7*yr-rc5(2))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*yr-rc5(2))
                  term7 = 2.0d0 * sc7 * rr7 * zr
                  term8 = zr * (sc7*rr9*yr-rc7(2))
                  tiyz = -ci*term1 + diz*term2 + diy*term3
     &                      - dir*term4 - qiyz*term5 + qiz*term6
     &                      + qiy*term7 - qir*term8
                  tkyz = -ck*term1 - dkz*term2 - dky*term3
     &                      + dkr*term4 - qkyz*term5 + qkz*term6
     &                      + qky*term7 - qkr*term8
c
c     get the field gradient for direct polarization force
c
               else if (use_chgpen) then
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*xr*xr
                  term1core = rr3 - rr5*xr*xr
                  term2i = 2.0d0*rr5*dmpi(5)*xr 
                  term3i = rr7*dmpi(7)*xr*xr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*xr
                  term6i = rr9*dmpi(9)*xr*xr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*xr*xr
                  term2k = 2.0d0*rr5*dmpk(5)*xr
                  term3k = rr7*dmpk(7)*xr*xr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*xr
                  term6k = rr9*dmpk(9)*xr*xr
                  tixx = vali*term1i + corei*term1core  
     &                      + dix*term2i - dir*term3i
     &                      - qixx*term4i + qix*term5i - qir*term6i
     &                      + (qiy*yr+qiz*zr)*rr7*dmpi(7)
                  tkxx = valk*term1k + corek*term1core
     &                      - dkx*term2k + dkr*term3k
     &                      - qkxx*term4k + qkx*term5k - qkr*term6k
     &                      + (qky*yr+qkz*zr)*rr7*dmpk(7)
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*yr*yr
                  term1core = rr3 - rr5*yr*yr
                  term2i = 2.0d0*rr5*dmpi(5)*yr
                  term3i = rr7*dmpi(7)*yr*yr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*yr
                  term6i = rr9*dmpi(9)*yr*yr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*yr*yr
                  term2k = 2.0d0*rr5*dmpk(5)*yr
                  term3k = rr7*dmpk(7)*yr*yr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*yr
                  term6k = rr9*dmpk(9)*yr*yr
                  tiyy = vali*term1i + corei*term1core
     &                      + diy*term2i - dir*term3i
     &                      - qiyy*term4i + qiy*term5i - qir*term6i
     &                      + (qix*xr+qiz*zr)*rr7*dmpi(7)
                  tkyy = valk*term1k + corek*term1core
     &                      - dky*term2k + dkr*term3k
     &                      - qkyy*term4k + qky*term5k - qkr*term6k
     &                      + (qkx*xr+qkz*zr)*rr7*dmpk(7)
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*zr*zr
                  term1core = rr3 - rr5*zr*zr
                  term2i = 2.0d0*rr5*dmpi(5)*zr
                  term3i = rr7*dmpi(7)*zr*zr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*zr
                  term6i = rr9*dmpi(9)*zr*zr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*zr*zr
                  term2k = 2.0d0*rr5*dmpk(5)*zr
                  term3k = rr7*dmpk(7)*zr*zr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*zr
                  term6k = rr9*dmpk(9)*zr*zr
                  tizz = vali*term1i + corei*term1core
     &                      + diz*term2i - dir*term3i
     &                      - qizz*term4i + qiz*term5i - qir*term6i
     &                      + (qix*xr+qiy*yr)*rr7*dmpi(7)
                  tkzz = valk*term1k + corek*term1core
     &                      - dkz*term2k + dkr*term3k
     &                      - qkzz*term4k + qkz*term5k - qkr*term6k
     &                      + (qkx*xr+qky*yr)*rr7*dmpk(7)
                  term2i = rr5*dmpi(5)*xr 
                  term1i = yr * term2i
                  term1core = rr5*xr*yr
                  term3i = rr5*dmpi(5)*yr
                  term4i = yr * (rr7*dmpi(7)*xr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*xr
                  term7i = 2.0d0*rr7*dmpi(7)*yr
                  term8i = yr*rr9*dmpi(9)*xr
                  term2k = rr5*dmpk(5)*xr
                  term1k = yr * term2k
                  term3k = rr5*dmpk(5)*yr
                  term4k = yr * (rr7*dmpk(7)*xr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*xr
                  term7k = 2.0d0*rr7*dmpk(7)*yr
                  term8k = yr*rr9*dmpk(9)*xr
                  tixy = -vali*term1i - corei*term1core 
     &                      + diy*term2i + dix*term3i
     &                      - dir*term4i - qixy*term5i + qiy*term6i
     &                      + qix*term7i - qir*term8i
                  tkxy = -valk*term1k - corek*term1core 
     &                      - dky*term2k - dkx*term3k
     &                      + dkr*term4k - qkxy*term5k + qky*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5*dmpi(5)*xr
                  term1i = zr * term2i
                  term1core = rr5*xr*zr
                  term3i = rr5*dmpi(5)*zr
                  term4i = zr * (rr7*dmpi(7)*xr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*xr
                  term7i = 2.0d0*rr7*dmpi(7)*zr
                  term8i = zr*rr9*dmpi(9)*xr
                  term2k = rr5*dmpk(5)*xr
                  term1k = zr * term2k
                  term3k = rr5*dmpk(5)*zr
                  term4k = zr * (rr7*dmpk(7)*xr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*xr
                  term7k = 2.0d0*rr7*dmpk(7)*zr
                  term8k = zr*rr9*dmpk(9)*xr
                  tixz = -vali*term1i - corei*term1core
     &                      + diz*term2i + dix*term3i
     &                      - dir*term4i - qixz*term5i + qiz*term6i
     &                      + qix*term7i - qir*term8i
                  tkxz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dkx*term3k
     &                      + dkr*term4k - qkxz*term5k + qkz*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5*dmpi(5)*yr
                  term1i = zr * term2i
                  term1core = rr5*yr*zr
                  term3i = rr5*dmpi(5)*zr
                  term4i = zr * (rr7*dmpi(7)*yr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*yr
                  term7i = 2.0d0*rr7*dmpi(7)*zr
                  term8i = zr*rr9*dmpi(9)*yr
                  term2k = rr5*dmpk(5)*yr
                  term1k = zr * term2k
                  term3k = rr5*dmpk(5)*zr
                  term4k = zr * (rr7*dmpk(7)*yr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*yr
                  term7k = 2.0d0*rr7*dmpk(7)*zr
                  term8k = zr*rr9*dmpk(9)*yr
                  tiyz = -vali*term1i - corei*term1core
     &                      + diz*term2i + diy*term3i
     &                      - dir*term4i - qiyz*term5i + qiz*term6i
     &                      + qiy*term7i - qir*term8i
                  tkyz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dky*term3k
     &                      + dkr*term4k - qkyz*term5k + qkz*term6k
     &                      + qky*term7k - qkr*term8k
               end if
c
c     get the dEd/dR terms for Thole direct polarization force
c
               if (use_thole) then
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      - tkxx*uixp - tkxy*uiyp - tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      - tkxy*uixp - tkyy*uiyp - tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      - tkxz*uixp - tkyz*uiyp - tkzz*uizp
                  frcx = dscale(k) * depx
                  frcy = dscale(k) * depy
                  frcz = dscale(k) * depz
c
c     get the dEp/dR terms for Thole direct polarization force
c
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = frcx + pscale(k)*depx
                  frcy = frcy + pscale(k)*depy
                  frcz = frcz + pscale(k)*depz
c
c     get the dEp/dR terms for chgpen direct polarization force
c
               else if (use_chgpen) then
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = 2.0d0*dscale(k)*depx
                  frcy = 2.0d0*dscale(k)*depy
                  frcz = 2.0d0*dscale(k)*depz
               end if
c
c     reset Thole values if alternate direct damping was used
c
               if (use_dirdamp) then
                  sc3 = 1.0d0
                  sc5 = 1.0d0
                  do j = 1, 3
                     rc3(j) = 0.0d0
                     rc5(j) = 0.0d0
                  end do
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        sc3 = 1.0d0 - expdamp
                        sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                        temp3 = damp * expdamp * rr5
                        temp5 = 3.0d0 * damp / r2
                        rc3(1) = xr * temp3
                        rc3(2) = yr * temp3
                        rc3(3) = zr * temp3
                        rc5(1) = rc3(1) * temp5
                        rc5(2) = rc3(2) * temp5
                        rc5(3) = rc3(3) * temp5
                     end if
                  end if
               end if
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp.eq.'MUTUAL' .and. use_thole) then
                  term1 = (sc3+sc5) * rr5
                  term2 = term1*xr - rc3(1)
                  term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr - rc3(2)
                  term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr - rc3(3)
                  term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = sc5 * rr5 * yr
                  term2 = sc3*rr5*xr - rc3(1)
                  term3 = yr * (sc5*rr7*xr-rc5(1))
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = sc5 * rr5 * zr
                  term3 = zr * (sc5*rr7*xr-rc5(1))
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = sc3*rr5*yr - rc3(2)
                  term3 = zr * (sc5*rr7*yr-rc5(2))
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + uscale(k)*depx
                  frcy = frcy + uscale(k)*depy
                  frcz = frcz + uscale(k)*depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               else if (poltyp.eq.'MUTUAL' .and. use_chgpen) then
                  term1 = 2.0d0 * dmpik(5) * rr5
                  term2 = term1*xr
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*xr*xr 
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr 
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*yr*yr 
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr 
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*zr*zr 
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = rr5*dmpik(5)*yr
                  term2 = rr5*dmpik(5)*xr 
                  term3 = yr * (rr7*dmpik(7)*xr)
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = rr5 *dmpik(5) * zr
                  term3 = zr * (rr7*dmpik(7)*xr)
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = rr5*dmpik(5)*yr 
                  term3 = zr * (rr7*dmpik(7)*yr)
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + wscale(k)*depx
                  frcy = frcy + wscale(k)*depy
                  frcz = frcz + wscale(k)*depz
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_thole) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,ii)*xr + uopt(j,2,ii)*yr
     &                          + uopt(j,3,ii)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,kk)*xr + uopt(m,2,kk)*yr
     &                             + uopt(m,3,kk)*zr
                        term1 = (sc3+sc5) * rr5
                        term2 = term1*xr - rc3(1)
                        term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                        tixx = uopt(j,1,ii)*term2 + uirm*term3
                        tkxx = uopt(m,1,kk)*term2 + ukrm*term3
                        term2 = term1*yr - rc3(2)
                        term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                        tiyy = uopt(j,2,ii)*term2 + uirm*term3
                        tkyy = uopt(m,2,kk)*term2 + ukrm*term3
                        term2 = term1*zr - rc3(3)
                        term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                        tizz = uopt(j,3,ii)*term2 + uirm*term3
                        tkzz = uopt(m,3,kk)*term2 + ukrm*term3
                        term1 = sc5 * rr5 * yr
                        term2 = sc3*rr5*xr - rc3(1)
                        term3 = yr * (sc5*rr7*xr-rc5(1))
                        tixy = uopt(j,1,ii)*term1 + uopt(j,2,ii)*term2
     &                            - uirm*term3
                        tkxy = uopt(m,1,kk)*term1 + uopt(m,2,kk)*term2
     &                            - ukrm*term3
                        term1 = sc5 * rr5 * zr
                        term3 = zr * (sc5*rr7*xr-rc5(1))
                        tixz = uopt(j,1,ii)*term1 + uopt(j,3,ii)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,kk)*term1 + uopt(m,3,kk)*term2
     &                            - ukrm*term3
                        term2 = sc3*rr5*yr - rc3(2)
                        term3 = zr * (sc5*rr7*yr-rc5(2))
                        tiyz = uopt(j,2,ii)*term1 + uopt(j,3,ii)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,kk)*term1 + uopt(m,3,kk)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,kk) + tkxx*uoptp(j,1,ii)
     &                       + tixy*uoptp(m,2,kk) + tkxy*uoptp(j,2,ii)
     &                       + tixz*uoptp(m,3,kk) + tkxz*uoptp(j,3,ii)
                        depy = tixy*uoptp(m,1,kk) + tkxy*uoptp(j,1,ii)
     &                       + tiyy*uoptp(m,2,kk) + tkyy*uoptp(j,2,ii)
     &                       + tiyz*uoptp(m,3,kk) + tkyz*uoptp(j,3,ii)
                        depz = tixz*uoptp(m,1,kk) + tkxz*uoptp(j,1,ii)
     &                       + tiyz*uoptp(m,2,kk) + tkyz*uoptp(j,2,ii)
     &                       + tizz*uoptp(m,3,kk) + tkzz*uoptp(j,3,ii)
                        frcx = frcx + copm(j+m+1)*uscale(k)*depx
                        frcy = frcy + copm(j+m+1)*uscale(k)*depy
                        frcz = frcz + copm(j+m+1)*uscale(k)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_chgpen) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,i)*xr + uopt(j,2,i)*yr
     &                          + uopt(j,3,i)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,k)*xr + uopt(m,2,k)*yr
     &                             + uopt(m,3,k)*zr
                        term1 = 2.0d0 * dmpik(5) * rr5
                        term2 = term1*xr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*xr*xr
                        tixx = uopt(j,1,i)*term2 + uirm*term3
                        tkxx = uopt(m,1,k)*term2 + ukrm*term3
                        term2 = term1*yr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*yr*yr
                        tiyy = uopt(j,2,i)*term2 + uirm*term3
                        tkyy = uopt(m,2,k)*term2 + ukrm*term3
                        term2 = term1*zr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*zr*zr
                        tizz = uopt(j,3,i)*term2 + uirm*term3
                        tkzz = uopt(m,3,k)*term2 + ukrm*term3
                        term1 = rr5*dmpik(5)*yr
                        term2 = rr5*dmpik(5)*xr
                        term3 = yr * (rr7*dmpik(7)*xr)
                        tixy = uopt(j,1,i)*term1 + uopt(j,2,i)*term2
     &                            - uirm*term3
                        tkxy = uopt(m,1,k)*term1 + uopt(m,2,k)*term2
     &                            - ukrm*term3
                        term1 = rr5 *dmpik(5) * zr
                        term3 = zr * (rr7*dmpik(7)*xr)
                        tixz = uopt(j,1,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        term2 = rr5*dmpik(5)*yr
                        term3 = zr * (rr7*dmpik(7)*yr)
                        tiyz = uopt(j,2,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,k) + tkxx*uoptp(j,1,i)
     &                       + tixy*uoptp(m,2,k) + tkxy*uoptp(j,2,i)
     &                       + tixz*uoptp(m,3,k) + tkxz*uoptp(j,3,i)
                        depy = tixy*uoptp(m,1,k) + tkxy*uoptp(j,1,i)
     &                       + tiyy*uoptp(m,2,k) + tkyy*uoptp(j,2,i)
     &                       + tiyz*uoptp(m,3,k) + tkyz*uoptp(j,3,i)
                        depz = tixz*uoptp(m,1,k) + tkxz*uoptp(j,1,i)
     &                       + tiyz*uoptp(m,2,k) + tkyz*uoptp(j,2,i)
     &                       + tizz*uoptp(m,3,k) + tkzz*uoptp(j,3,i)
                        frcx = frcx + copm(j+m+1)*wscale(k)*depx
                        frcy = frcy + copm(j+m+1)*wscale(k)*depy
                        frcz = frcz + copm(j+m+1)*wscale(k)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for TCG polarization force
c
               else if (poltyp.eq.'TCG' .and. use_thole) then
                  do j = 1, tcgnab
                     ukx = ubd(1,kk,j)
                     uky = ubd(2,kk,j)
                     ukz = ubd(3,kk,j)
                     ukxp = ubp(1,kk,j)
                     ukyp = ubp(2,kk,j)
                     ukzp = ubp(3,kk,j)
                     uirt = uax(j)*xr + uay(j)*yr + uaz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = (sc3+sc5) * rr5
                     term2 = term1*xr - rc3(1)
                     term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                     tixx = uax(j)*term2 + uirt*term3
                     tkxx = ukx*term2 + ukrt*term3
                     term2 = term1*yr - rc3(2)
                     term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                     tiyy = uay(j)*term2 + uirt*term3
                     tkyy = uky*term2 + ukrt*term3
                     term2 = term1*zr - rc3(3)
                     term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                     tizz = uaz(j)*term2 + uirt*term3
                     tkzz = ukz*term2 + ukrt*term3
                     term1 = sc5 * rr5 * yr
                     term2 = sc3*rr5*xr - rc3(1)
                     term3 = yr * (sc5*rr7*xr-rc5(1))
                     tixy = uax(j)*term1 + uay(j)*term2 - uirt*term3
                     tkxy = ukx*term1 + uky*term2 - ukrt*term3
                     term1 = sc5 * rr5 * zr
                     term3 = zr * (sc5*rr7*xr-rc5(1))
                     tixz = uax(j)*term1 + uaz(j)*term2 - uirt*term3
                     tkxz = ukx*term1 + ukz*term2 - ukrt*term3
                     term2 = sc3*rr5*yr - rc3(2)
                     term3 = zr * (sc5*rr7*yr-rc5(2))
                     tiyz = uay(j)*term1 + uaz(j)*term2 - uirt*term3
                     tkyz = uky*term1 + ukz*term2 - ukrt*term3
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*uaxp(j) + tkxy*uayp(j)
     &                         + tkxz*uazp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*uaxp(j) + tkyy*uayp(j)
     &                         + tkyz*uazp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*uaxp(j) + tkyz*uayp(j)
     &                         + tkzz*uazp(j)
                     frcx = frcx + uscale(k)*depx
                     frcy = frcy + uscale(k)*depy
                     frcz = frcz + uscale(k)*depz
                     ukx = uad(1,kk,j)
                     uky = uad(2,kk,j)
                     ukz = uad(3,kk,j)
                     ukxp = uap(1,kk,j)
                     ukyp = uap(2,kk,j)
                     ukzp = uap(3,kk,j)
                     uirt = ubx(j)*xr + uby(j)*yr + ubz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = (sc3+sc5) * rr5
                     term2 = term1*xr - rc3(1)
                     term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                     tixx = ubx(j)*term2 + uirt*term3
                     tkxx = ukx*term2 + ukrt*term3
                     term2 = term1*yr - rc3(2)
                     term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                     tiyy = uby(j)*term2 + uirt*term3
                     tkyy = uky*term2 + ukrt*term3
                     term2 = term1*zr - rc3(3)
                     term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                     tizz = ubz(j)*term2 + uirt*term3
                     tkzz = ukz*term2 + ukrt*term3
                     term1 = sc5 * rr5 * yr
                     term2 = sc3*rr5*xr - rc3(1)
                     term3 = yr * (sc5*rr7*xr-rc5(1))
                     tixy = ubx(j)*term1 + uby(j)*term2 - uirt*term3
                     tkxy = ukx*term1 + uky*term2 - ukrt*term3
                     term1 = sc5 * rr5 * zr
                     term3 = zr * (sc5*rr7*xr-rc5(1))
                     tixz = ubx(j)*term1 + ubz(j)*term2 - uirt*term3
                     tkxz = ukx*term1 + ukz*term2 - ukrt*term3
                     term2 = sc3*rr5*yr - rc3(2)
                     term3 = zr * (sc5*rr7*yr-rc5(2))
                     tiyz = uby(j)*term1 + ubz(j)*term2 - uirt*term3
                     tkyz = uky*term1 + ukz*term2 - ukrt*term3
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*ubxp(j) + tkxy*ubyp(j)
     &                         + tkxz*ubzp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*ubxp(j) + tkyy*ubyp(j)
     &                         + tkyz*ubzp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*ubxp(j) + tkyz*ubyp(j)
     &                         + tkzz*ubzp(j)
                     frcx = frcx + uscale(k)*depx
                     frcy = frcy + uscale(k)*depy
                     frcz = frcz + uscale(k)*depz
                  end do
               end if
c
c     force and torque components scaled for self-interactions
c
               if (i .eq. k) then
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
c     increment force-based gradient on the interaction sites
c
               dep(1,i) = dep(1,i) + frcx
               dep(2,i) = dep(2,i) + frcy
               dep(3,i) = dep(3,i) + frcz
               dep(1,k) = dep(1,k) - frcx
               dep(2,k) = dep(2,k) - frcy
               dep(3,k) = dep(3,k) - frcz
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
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = 1.0d0
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
               uscale(ip14(j,i)) = 1.0d0
            end do
         end if
      end do
      end if
c
c     torque is induced field and gradient cross permanent moments
c
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         tep(1) = diz*ufld(2,i) - diy*ufld(3,i)
     &               + qixz*dufld(2,i) - qixy*dufld(4,i)
     &               + 2.0d0*qiyz*(dufld(3,i)-dufld(6,i))
     &               + (qizz-qiyy)*dufld(5,i)
         tep(2) = dix*ufld(3,i) - diz*ufld(1,i)
     &               - qiyz*dufld(2,i) + qixy*dufld(5,i)
     &               + 2.0d0*qixz*(dufld(6,i)-dufld(1,i))
     &               + (qixx-qizz)*dufld(4,i)
         tep(3) = diy*ufld(1,i) - dix*ufld(2,i)
     &               + qiyz*dufld(4,i) - qixz*dufld(5,i)
     &               + 2.0d0*qixy*(dufld(1,i)-dufld(3,i))
     &               + (qiyy-qixx)*dufld(2,i)
         call torque (ii,tep,fix,fiy,fiz,dep)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = abs(yaxis(ii))
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
            dep(1,i) = dep(1,i) + frcx
            dep(2,i) = dep(2,i) + frcy
            dep(3,i) = dep(3,i) + frcz
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
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (wscale)
      deallocate (ufld)
      deallocate (dufld)
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
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use deriv
      use energi
      use molcul
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polopt
      use polpot
      use poltcg
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      integer ix,iy,iz
      real*8 f,pgamma
      real*8 pdi,pti,ddi
      real*8 damp,expdamp
      real*8 temp3,temp5,temp7
      real*8 sc3,sc5,sc7
      real*8 sr3,sr5,sr7
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 dsr3i,dsr5i,dsr7i
      real*8 dsr3k,dsr5k,dsr7k
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 uixp,uiyp,uizp
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 ukxp,ukyp,ukzp
      real*8 dir,uir,uirp
      real*8 dkr,ukr,ukrp
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 uirm,ukrm
      real*8 uirt,ukrt
      real*8 tuir,tukr
      real*8 tixx,tiyy,tizz
      real*8 tixy,tixz,tiyz
      real*8 tkxx,tkyy,tkzz
      real*8 tkxy,tkxz,tkyz
      real*8 tix3,tiy3,tiz3
      real*8 tix5,tiy5,tiz5
      real*8 tkx3,tky3,tkz3
      real*8 tkx5,tky5,tkz5
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 term7,term8
      real*8 term1core
      real*8 term1i,term2i,term3i
      real*8 term4i,term5i,term6i
      real*8 term7i,term8i
      real*8 term1k,term2k,term3k
      real*8 term4k,term5k,term6k
      real*8 term7k,term8k
      real*8 poti,potk
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 tep(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8 uax(3),uay(3),uaz(3)
      real*8 ubx(3),uby(3),ubz(3)
      real*8 uaxp(3),uayp(3),uazp(3)
      real*8 ubxp(3),ubyp(3),ubzp(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
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
c     compute the total induced dipole polarization energy
c
      call epolar1e
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (wscale(n))
      allocate (ufld(3,n))
      allocate (dufld(6,n))
      allocate (pot(n))
      allocate (decfx(n))
      allocate (decfy(n))
      allocate (decfz(n))
c
c     set exclusion coefficients and arrays to store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
         wscale(i) = 1.0d0
         do j = 1, 3
            ufld(j,i) = 0.0d0
         end do
         do j = 1, 6
            dufld(j,i) = 0.0d0
         end do
         pot(i) = 0.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,x,y,z,rpole,uind,
!$OMP& uinp,pdamp,thole,dirdamp,pcore,pval,palpha,n12,i12,n13,i13,n14,
!$OMP& i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,p2scale,
!$OMP& p3scale,p4scale,p5scale,p2iscale,p3iscale,p4iscale,p5iscale,
!$OMP& d1scale,d2scale,d3scale,d4scale,u1scale,u2scale,u3scale,u4scale,
!$OMP& w2scale,w3scale,w4scale,w5scale,nelst,elst,dpequal,use_thole,
!$OMP& use_dirdamp,use_chgpen,use_chgflx,use_bounds,off2,f,molcule,
!$OMP& optorder,copm,uopt,uoptp,poltyp,tcgnab,uad,uap,ubd,ubp,
!$OMP& xaxis,yaxis,zaxis)
!$OMP& shared (dep,ufld,dufld,pot,vir)
!$OMP& firstprivate(pscale,dscale,uscale,wscale)
!$OMP DO reduction(+:dep,ufld,dufld,pot,vir) schedule(guided)
c
c     compute the dipole polarization gradient components
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
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         uixp = uinp(1,ii)
         uiyp = uinp(2,ii)
         uizp = uinp(3,ii)
         do j = 1, tcgnab
            uax(j) = uad(1,ii,j)
            uay(j) = uad(2,ii,j)
            uaz(j) = uad(3,ii,j)
            uaxp(j) = uap(1,ii,j)
            uayp(j) = uap(2,ii,j)
            uazp(j) = uap(3,ii,j)
            ubx(j) = ubd(1,ii,j)
            uby(j) = ubd(2,ii,j)
            ubz(j) = ubd(3,ii,j)
            ubxp(j) = ubp(1,ii,j)
            ubyp(j) = ubp(2,ii,j)
            ubzp(j) = ubp(3,ii,j)
         end do
         if (use_thole) then
            pdi = pdamp(ii)
            pti = thole(ii)
            ddi = dirdamp(ii)
         else if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               dscale(i12(j,i)) = pscale(i12(j,i))
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = u4scale
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
               uscale(ip14(j,i)) = u4scale
            end do
         end if
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
               ukx = uind(1,kk)
               uky = uind(2,kk)
               ukz = uind(3,kk)
               ukxp = uinp(1,kk)
               ukyp = uinp(2,kk)
               ukzp = uinp(3,kk)
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
               uir = uix*xr + uiy*yr + uiz*zr
               uirp = uixp*xr + uiyp*yr + uizp*zr
               ukr = ukx*xr + uky*yr + ukz*zr
               ukrp = ukxp*xr + ukyp*yr + ukzp*zr
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     set initial values for tha damping scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
c
c     apply Thole polarization damping to scale factors
c
               if (use_thole) then
                  damp = pdi * pdamp(kk)
                  if (use_dirdamp) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(ddi,dirdamp(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           sc3 = 1.0d0 - expdamp 
                           sc5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                      +0.15d0*damp**2)
                           temp3 = 0.5d0 * damp * expdamp 
                           temp5 = 1.5d0 * (1.0d0+damp)
                           temp7 = 5.0d0*(1.5d0*damp*expdamp
     &                                *(0.35d0+0.35d0*damp
     &                                   +0.15d0*damp**2))/(temp3*temp5)
                           temp3 = temp3 * rr5
                           temp5 = temp5 / r2
                           temp7 = temp7 / r2
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
                  else
                     pgamma = min(pti,thole(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(pti,thole(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           sc3 = 1.0d0 - expdamp
                           sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                      +0.6d0*damp**2)
                           temp3 = damp * expdamp * rr5
                           temp5 = 3.0d0 * damp / r2
                           temp7 = (-1.0d0+3.0d0*damp) / r2
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
                  end if
                  sr3 = rr3 * sc3
                  sr5 = rr5 * sc5
                  sr7 = rr7 * sc7
                  dsr3 = sr3 * dscale(k)
                  dsr5 = sr5 * dscale(k)
                  dsr7 = sr7 * dscale(k)
                  psr3 = sr3 * pscale(k)
                  psr5 = sr5 * pscale(k)
                  psr7 = sr7 * pscale(k)
c
c     apply charge penetration damping to scale factors
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call damppole (r,9,alphai,alphak,dmpi,dmpk,dmpik)
                  dsr3i = 2.0d0 * rr3 * dmpi(3) * dscale(k)
                  dsr5i = 2.0d0 * rr5 * dmpi(5) * dscale(k)
                  dsr7i = 2.0d0 * rr7 * dmpi(7) * dscale(k)
                  dsr3k = 2.0d0 * rr3 * dmpk(3) * dscale(k)
                  dsr5k = 2.0d0 * rr5 * dmpk(5) * dscale(k)
                  dsr7k = 2.0d0 * rr7 * dmpk(7) * dscale(k)
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_thole) then
                     poti = -ukr*psr3 - ukrp*dsr3
                     potk = uir*psr3 + uirp*dsr3
                  else if (use_chgpen) then
                     poti = -ukr * dsr3i
                     potk = uir * dsr3k
                  end if
                  pot(i) = pot(i) + poti 
                  pot(k) = pot(k) + potk 
               end if
c
c     get the induced dipole field used for dipole torques
c
               if (use_thole) then
                  tix3 = psr3*ukx + dsr3*ukxp
                  tiy3 = psr3*uky + dsr3*ukyp
                  tiz3 = psr3*ukz + dsr3*ukzp
                  tkx3 = psr3*uix + dsr3*uixp
                  tky3 = psr3*uiy + dsr3*uiyp
                  tkz3 = psr3*uiz + dsr3*uizp
                  tuir = -psr5*ukr - dsr5*ukrp
                  tukr = -psr5*uir - dsr5*uirp
               else if (use_chgpen) then
                  tix3 = dsr3i*ukx
                  tiy3 = dsr3i*uky
                  tiz3 = dsr3i*ukz
                  tkx3 = dsr3k*uix
                  tky3 = dsr3k*uiy
                  tkz3 = dsr3k*uiz
                  tuir = -dsr5i*ukr
                  tukr = -dsr5k*uir
               end if
               ufld(1,i) = ufld(1,i) + tix3 + xr*tuir
               ufld(2,i) = ufld(2,i) + tiy3 + yr*tuir
               ufld(3,i) = ufld(3,i) + tiz3 + zr*tuir
               ufld(1,k) = ufld(1,k) + tkx3 + xr*tukr
               ufld(2,k) = ufld(2,k) + tky3 + yr*tukr
               ufld(3,k) = ufld(3,k) + tkz3 + zr*tukr
c
c     get induced dipole field gradient used for quadrupole torques
c
               if (use_thole) then
                  tix5 = 2.0d0 * (psr5*ukx+dsr5*ukxp)
                  tiy5 = 2.0d0 * (psr5*uky+dsr5*ukyp)
                  tiz5 = 2.0d0 * (psr5*ukz+dsr5*ukzp)
                  tkx5 = 2.0d0 * (psr5*uix+dsr5*uixp)
                  tky5 = 2.0d0 * (psr5*uiy+dsr5*uiyp)
                  tkz5 = 2.0d0 * (psr5*uiz+dsr5*uizp)
                  tuir = -psr7*ukr - dsr7*ukrp
                  tukr = -psr7*uir - dsr7*uirp
               else if (use_chgpen) then
                  tix5 = 2.0d0 * (dsr5i*ukx)
                  tiy5 = 2.0d0 * (dsr5i*uky)
                  tiz5 = 2.0d0 * (dsr5i*ukz)
                  tkx5 = 2.0d0 * (dsr5k*uix)
                  tky5 = 2.0d0 * (dsr5k*uiy)
                  tkz5 = 2.0d0 * (dsr5k*uiz)
                  tuir = -dsr7i*ukr
                  tukr = -dsr7k*uir
               end if
               dufld(1,i) = dufld(1,i) + xr*tix5 + xr*xr*tuir
               dufld(2,i) = dufld(2,i) + xr*tiy5 + yr*tix5
     &                         + 2.0d0*xr*yr*tuir
               dufld(3,i) = dufld(3,i) + yr*tiy5 + yr*yr*tuir
               dufld(4,i) = dufld(4,i) + xr*tiz5 + zr*tix5
     &                         + 2.0d0*xr*zr*tuir
               dufld(5,i) = dufld(5,i) + yr*tiz5 + zr*tiy5
     &                         + 2.0d0*yr*zr*tuir
               dufld(6,i) = dufld(6,i) + zr*tiz5 + zr*zr*tuir
               dufld(1,k) = dufld(1,k) - xr*tkx5 - xr*xr*tukr
               dufld(2,k) = dufld(2,k) - xr*tky5 - yr*tkx5
     &                         - 2.0d0*xr*yr*tukr
               dufld(3,k) = dufld(3,k) - yr*tky5 - yr*yr*tukr
               dufld(4,k) = dufld(4,k) - xr*tkz5 - zr*tkx5
     &                         - 2.0d0*xr*zr*tukr
               dufld(5,k) = dufld(5,k) - yr*tkz5 - zr*tky5
     &                         - 2.0d0*yr*zr*tukr
               dufld(6,k) = dufld(6,k) - zr*tkz5 - zr*zr*tukr
c
c     get the field gradient for direct polarization force
c
               if (use_thole) then
                  term1 = sc3*(rr3-rr5*xr*xr) + rc3(1)*xr
                  term2 = (sc3+sc5)*rr5*xr - rc3(1)
                  term3 = sc5*(rr7*xr*xr-rr5) - rc5(1)*xr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*xr-rc5(1)+1.5d0*sc7*rr7*xr)
                  term6 = xr * (sc7*rr9*xr-rc7(1))
                  tixx = ci*term1 + dix*term2 - dir*term3
     &                      - qixx*term4 + qix*term5 - qir*term6
     &                      + (qiy*yr+qiz*zr)*sc7*rr7
                  tkxx = ck*term1 - dkx*term2 + dkr*term3
     &                      - qkxx*term4 + qkx*term5 - qkr*term6
     &                      + (qky*yr+qkz*zr)*sc7*rr7
                  term1 = sc3*(rr3-rr5*yr*yr) + rc3(2)*yr
                  term2 = (sc3+sc5)*rr5*yr - rc3(2)
                  term3 = sc5*(rr7*yr*yr-rr5) - rc5(2)*yr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*yr-rc5(2)+1.5d0*sc7*rr7*yr)
                  term6 = yr * (sc7*rr9*yr-rc7(2))
                  tiyy = ci*term1 + diy*term2 - dir*term3
     &                      - qiyy*term4 + qiy*term5 - qir*term6
     &                      + (qix*xr+qiz*zr)*sc7*rr7
                  tkyy = ck*term1 - dky*term2 + dkr*term3
     &                      - qkyy*term4 + qky*term5 - qkr*term6
     &                      + (qkx*xr+qkz*zr)*sc7*rr7
                  term1 = sc3*(rr3-rr5*zr*zr) + rc3(3)*zr
                  term2 = (sc3+sc5)*rr5*zr - rc3(3)
                  term3 = sc5*(rr7*zr*zr-rr5) - rc5(3)*zr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*zr-rc5(3)+1.5d0*sc7*rr7*zr)
                  term6 = zr * (sc7*rr9*zr-rc7(3))
                  tizz = ci*term1 + diz*term2 - dir*term3
     &                      - qizz*term4 + qiz*term5 - qir*term6
     &                      + (qix*xr+qiy*yr)*sc7*rr7
                  tkzz = ck*term1 - dkz*term2 + dkr*term3
     &                      - qkzz*term4 + qkz*term5 - qkr*term6
     &                      + (qkx*xr+qky*yr)*sc7*rr7
                  term2 = sc3*rr5*xr - rc3(1)
                  term1 = yr * term2
                  term3 = sc5 * rr5 * yr
                  term4 = yr * (sc5*rr7*xr-rc5(1))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
                  term7 = 2.0d0 * sc7 * rr7 * yr
                  term8 = yr * (sc7*rr9*xr-rc7(1))
                  tixy = -ci*term1 + diy*term2 + dix*term3
     &                      - dir*term4 - qixy*term5 + qiy*term6
     &                      + qix*term7 - qir*term8
                  tkxy = -ck*term1 - dky*term2 - dkx*term3
     &                      + dkr*term4 - qkxy*term5 + qky*term6
     &                      + qkx*term7 - qkr*term8
                  term2 = sc3*rr5*xr - rc3(1)
                  term1 = zr * term2
                  term3 = sc5 * rr5 * zr
                  term4 = zr * (sc5*rr7*xr-rc5(1))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
                  term7 = 2.0d0 * sc7 * rr7 * zr
                  term8 = zr * (sc7*rr9*xr-rc7(1))
                  tixz = -ci*term1 + diz*term2 + dix*term3
     &                      - dir*term4 - qixz*term5 + qiz*term6
     &                      + qix*term7 - qir*term8
                  tkxz = -ck*term1 - dkz*term2 - dkx*term3
     &                      + dkr*term4 - qkxz*term5 + qkz*term6
     &                      + qkx*term7 - qkr*term8
                  term2 = sc3*rr5*yr - rc3(2)
                  term1 = zr * term2
                  term3 = sc5 * rr5 * zr
                  term4 = zr * (sc5*rr7*yr-rc5(2))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*yr-rc5(2))
                  term7 = 2.0d0 * sc7 * rr7 * zr
                  term8 = zr * (sc7*rr9*yr-rc7(2))
                  tiyz = -ci*term1 + diz*term2 + diy*term3
     &                      - dir*term4 - qiyz*term5 + qiz*term6
     &                      + qiy*term7 - qir*term8
                  tkyz = -ck*term1 - dkz*term2 - dky*term3
     &                      + dkr*term4 - qkyz*term5 + qkz*term6
     &                      + qky*term7 - qkr*term8
c
c     get the field gradient for direct polarization force
c
               else if (use_chgpen) then
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*xr*xr
                  term1core = rr3 - rr5*xr*xr
                  term2i = 2.0d0*rr5*dmpi(5)*xr 
                  term3i = rr7*dmpi(7)*xr*xr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*xr
                  term6i = rr9*dmpi(9)*xr*xr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*xr*xr
                  term2k = 2.0d0*rr5*dmpk(5)*xr
                  term3k = rr7*dmpk(7)*xr*xr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*xr
                  term6k = rr9*dmpk(9)*xr*xr
                  tixx = vali*term1i + corei*term1core  
     &                      + dix*term2i - dir*term3i
     &                      - qixx*term4i + qix*term5i - qir*term6i
     &                      + (qiy*yr+qiz*zr)*rr7*dmpi(7)
                  tkxx = valk*term1k + corek*term1core
     &                      - dkx*term2k + dkr*term3k
     &                      - qkxx*term4k + qkx*term5k - qkr*term6k
     &                      + (qky*yr+qkz*zr)*rr7*dmpk(7)
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*yr*yr
                  term1core = rr3 - rr5*yr*yr
                  term2i = 2.0d0*rr5*dmpi(5)*yr
                  term3i = rr7*dmpi(7)*yr*yr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*yr
                  term6i = rr9*dmpi(9)*yr*yr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*yr*yr
                  term2k = 2.0d0*rr5*dmpk(5)*yr
                  term3k = rr7*dmpk(7)*yr*yr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*yr
                  term6k = rr9*dmpk(9)*yr*yr
                  tiyy = vali*term1i + corei*term1core
     &                      + diy*term2i - dir*term3i
     &                      - qiyy*term4i + qiy*term5i - qir*term6i
     &                      + (qix*xr+qiz*zr)*rr7*dmpi(7)
                  tkyy = valk*term1k + corek*term1core
     &                      - dky*term2k + dkr*term3k
     &                      - qkyy*term4k + qky*term5k - qkr*term6k
     &                      + (qkx*xr+qkz*zr)*rr7*dmpk(7)
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*zr*zr
                  term1core = rr3 - rr5*zr*zr
                  term2i = 2.0d0*rr5*dmpi(5)*zr
                  term3i = rr7*dmpi(7)*zr*zr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*zr
                  term6i = rr9*dmpi(9)*zr*zr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*zr*zr
                  term2k = 2.0d0*rr5*dmpk(5)*zr
                  term3k = rr7*dmpk(7)*zr*zr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*zr
                  term6k = rr9*dmpk(9)*zr*zr
                  tizz = vali*term1i + corei*term1core
     &                      + diz*term2i - dir*term3i
     &                      - qizz*term4i + qiz*term5i - qir*term6i
     &                      + (qix*xr+qiy*yr)*rr7*dmpi(7)
                  tkzz = valk*term1k + corek*term1core
     &                      - dkz*term2k + dkr*term3k
     &                      - qkzz*term4k + qkz*term5k - qkr*term6k
     &                      + (qkx*xr+qky*yr)*rr7*dmpk(7)
                  term2i = rr5*dmpi(5)*xr 
                  term1i = yr * term2i
                  term1core = rr5*xr*yr
                  term3i = rr5*dmpi(5)*yr
                  term4i = yr * (rr7*dmpi(7)*xr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*xr
                  term7i = 2.0d0*rr7*dmpi(7)*yr
                  term8i = yr*rr9*dmpi(9)*xr
                  term2k = rr5*dmpk(5)*xr
                  term1k = yr * term2k
                  term3k = rr5*dmpk(5)*yr
                  term4k = yr * (rr7*dmpk(7)*xr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*xr
                  term7k = 2.0d0*rr7*dmpk(7)*yr
                  term8k = yr*rr9*dmpk(9)*xr
                  tixy = -vali*term1i - corei*term1core 
     &                      + diy*term2i + dix*term3i
     &                      - dir*term4i - qixy*term5i + qiy*term6i
     &                      + qix*term7i - qir*term8i
                  tkxy = -valk*term1k - corek*term1core 
     &                      - dky*term2k - dkx*term3k
     &                      + dkr*term4k - qkxy*term5k + qky*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5*dmpi(5)*xr
                  term1i = zr * term2i
                  term1core = rr5*xr*zr
                  term3i = rr5*dmpi(5)*zr
                  term4i = zr * (rr7*dmpi(7)*xr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*xr
                  term7i = 2.0d0*rr7*dmpi(7)*zr
                  term8i = zr*rr9*dmpi(9)*xr
                  term2k = rr5*dmpk(5)*xr
                  term1k = zr * term2k
                  term3k = rr5*dmpk(5)*zr
                  term4k = zr * (rr7*dmpk(7)*xr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*xr
                  term7k = 2.0d0*rr7*dmpk(7)*zr
                  term8k = zr*rr9*dmpk(9)*xr
                  tixz = -vali*term1i - corei*term1core
     &                      + diz*term2i + dix*term3i
     &                      - dir*term4i - qixz*term5i + qiz*term6i
     &                      + qix*term7i - qir*term8i
                  tkxz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dkx*term3k
     &                      + dkr*term4k - qkxz*term5k + qkz*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5*dmpi(5)*yr
                  term1i = zr * term2i
                  term1core = rr5*yr*zr
                  term3i = rr5*dmpi(5)*zr
                  term4i = zr * (rr7*dmpi(7)*yr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*yr
                  term7i = 2.0d0*rr7*dmpi(7)*zr
                  term8i = zr*rr9*dmpi(9)*yr
                  term2k = rr5*dmpk(5)*yr
                  term1k = zr * term2k
                  term3k = rr5*dmpk(5)*zr
                  term4k = zr * (rr7*dmpk(7)*yr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*yr
                  term7k = 2.0d0*rr7*dmpk(7)*zr
                  term8k = zr*rr9*dmpk(9)*yr
                  tiyz = -vali*term1i - corei*term1core
     &                      + diz*term2i + diy*term3i
     &                      - dir*term4i - qiyz*term5i + qiz*term6i
     &                      + qiy*term7i - qir*term8i
                  tkyz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dky*term3k
     &                      + dkr*term4k - qkyz*term5k + qkz*term6k
     &                      + qky*term7k - qkr*term8k
               end if
c
c     get the dEd/dR terms for Thole direct polarization force
c
               if (use_thole) then
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      - tkxx*uixp - tkxy*uiyp - tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      - tkxy*uixp - tkyy*uiyp - tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      - tkxz*uixp - tkyz*uiyp - tkzz*uizp
                  frcx = dscale(k) * depx
                  frcy = dscale(k) * depy
                  frcz = dscale(k) * depz
c
c     get the dEp/dR terms for Thole direct polarization force
c
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = frcx + pscale(k)*depx
                  frcy = frcy + pscale(k)*depy
                  frcz = frcz + pscale(k)*depz
c
c     get the dEp/dR terms for chgpen direct polarization force
c
               else if (use_chgpen) then
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = 2.0d0*dscale(k)*depx
                  frcy = 2.0d0*dscale(k)*depy
                  frcz = 2.0d0*dscale(k)*depz
               end if
c
c     reset Thole values if alternate direct damping was used
c
               if (use_dirdamp) then
                  sc3 = 1.0d0
                  sc5 = 1.0d0
                  do j = 1, 3
                     rc3(j) = 0.0d0
                     rc5(j) = 0.0d0
                  end do
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        sc3 = 1.0d0 - expdamp
                        sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                        temp3 = damp * expdamp * rr5
                        temp5 = 3.0d0 * damp / r2
                        rc3(1) = xr * temp3
                        rc3(2) = yr * temp3
                        rc3(3) = zr * temp3
                        rc5(1) = rc3(1) * temp5
                        rc5(2) = rc3(2) * temp5
                        rc5(3) = rc3(3) * temp5
                     end if
                  end if
               end if
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp.eq.'MUTUAL' .and. use_thole) then
                  term1 = (sc3+sc5) * rr5
                  term2 = term1*xr - rc3(1)
                  term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr - rc3(2)
                  term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr - rc3(3)
                  term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = sc5 * rr5 * yr
                  term2 = sc3*rr5*xr - rc3(1)
                  term3 = yr * (sc5*rr7*xr-rc5(1))
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = sc5 * rr5 * zr
                  term3 = zr * (sc5*rr7*xr-rc5(1))
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = sc3*rr5*yr - rc3(2)
                  term3 = zr * (sc5*rr7*yr-rc5(2))
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + uscale(k)*depx
                  frcy = frcy + uscale(k)*depy
                  frcz = frcz + uscale(k)*depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               else if (poltyp.eq.'MUTUAL' .and. use_chgpen) then
                  term1 = 2.0d0 * dmpik(5) * rr5
                  term2 = term1*xr
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*xr*xr 
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr 
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*yr*yr 
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr 
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*zr*zr 
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = rr5*dmpik(5)*yr
                  term2 = rr5*dmpik(5)*xr 
                  term3 = yr * (rr7*dmpik(7)*xr)
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = rr5 *dmpik(5) * zr
                  term3 = zr * (rr7*dmpik(7)*xr)
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = rr5*dmpik(5)*yr 
                  term3 = zr * (rr7*dmpik(7)*yr)
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + wscale(k)*depx
                  frcy = frcy + wscale(k)*depy
                  frcz = frcz + wscale(k)*depz
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_thole) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,ii)*xr + uopt(j,2,ii)*yr
     &                          + uopt(j,3,ii)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,kk)*xr + uopt(m,2,kk)*yr
     &                             + uopt(m,3,kk)*zr
                        term1 = (sc3+sc5) * rr5
                        term2 = term1*xr - rc3(1)
                        term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                        tixx = uopt(j,1,ii)*term2 + uirm*term3
                        tkxx = uopt(m,1,kk)*term2 + ukrm*term3
                        term2 = term1*yr - rc3(2)
                        term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                        tiyy = uopt(j,2,ii)*term2 + uirm*term3
                        tkyy = uopt(m,2,kk)*term2 + ukrm*term3
                        term2 = term1*zr - rc3(3)
                        term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                        tizz = uopt(j,3,ii)*term2 + uirm*term3
                        tkzz = uopt(m,3,kk)*term2 + ukrm*term3
                        term1 = sc5 * rr5 * yr
                        term2 = sc3*rr5*xr - rc3(1)
                        term3 = yr * (sc5*rr7*xr-rc5(1))
                        tixy = uopt(j,1,ii)*term1 + uopt(j,2,ii)*term2
     &                            - uirm*term3
                        tkxy = uopt(m,1,kk)*term1 + uopt(m,2,kk)*term2
     &                            - ukrm*term3
                        term1 = sc5 * rr5 * zr
                        term3 = zr * (sc5*rr7*xr-rc5(1))
                        tixz = uopt(j,1,ii)*term1 + uopt(j,3,ii)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,kk)*term1 + uopt(m,3,kk)*term2
     &                            - ukrm*term3
                        term2 = sc3*rr5*yr - rc3(2)
                        term3 = zr * (sc5*rr7*yr-rc5(2))
                        tiyz = uopt(j,2,ii)*term1 + uopt(j,3,ii)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,kk)*term1 + uopt(m,3,kk)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,kk) + tkxx*uoptp(j,1,ii)
     &                       + tixy*uoptp(m,2,kk) + tkxy*uoptp(j,2,ii)
     &                       + tixz*uoptp(m,3,kk) + tkxz*uoptp(j,3,ii)
                        depy = tixy*uoptp(m,1,kk) + tkxy*uoptp(j,1,ii)
     &                       + tiyy*uoptp(m,2,kk) + tkyy*uoptp(j,2,ii)
     &                       + tiyz*uoptp(m,3,kk) + tkyz*uoptp(j,3,ii)
                        depz = tixz*uoptp(m,1,kk) + tkxz*uoptp(j,1,ii)
     &                       + tiyz*uoptp(m,2,kk) + tkyz*uoptp(j,2,ii)
     &                       + tizz*uoptp(m,3,kk) + tkzz*uoptp(j,3,ii)
                        frcx = frcx + copm(j+m+1)*uscale(k)*depx
                        frcy = frcy + copm(j+m+1)*uscale(k)*depy
                        frcz = frcz + copm(j+m+1)*uscale(k)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_chgpen) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,i)*xr + uopt(j,2,i)*yr
     &                          + uopt(j,3,i)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,k)*xr + uopt(m,2,k)*yr
     &                             + uopt(m,3,k)*zr
                        term1 = 2.0d0 * dmpik(5) * rr5
                        term2 = term1*xr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*xr*xr
                        tixx = uopt(j,1,i)*term2 + uirm*term3
                        tkxx = uopt(m,1,k)*term2 + ukrm*term3
                        term2 = term1*yr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*yr*yr
                        tiyy = uopt(j,2,i)*term2 + uirm*term3
                        tkyy = uopt(m,2,k)*term2 + ukrm*term3
                        term2 = term1*zr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*zr*zr
                        tizz = uopt(j,3,i)*term2 + uirm*term3
                        tkzz = uopt(m,3,k)*term2 + ukrm*term3
                        term1 = rr5*dmpik(5)*yr
                        term2 = rr5*dmpik(5)*xr
                        term3 = yr * (rr7*dmpik(7)*xr)
                        tixy = uopt(j,1,i)*term1 + uopt(j,2,i)*term2
     &                            - uirm*term3
                        tkxy = uopt(m,1,k)*term1 + uopt(m,2,k)*term2
     &                            - ukrm*term3
                        term1 = rr5 *dmpik(5) * zr
                        term3 = zr * (rr7*dmpik(7)*xr)
                        tixz = uopt(j,1,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        term2 = rr5*dmpik(5)*yr
                        term3 = zr * (rr7*dmpik(7)*yr)
                        tiyz = uopt(j,2,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,k) + tkxx*uoptp(j,1,i)
     &                       + tixy*uoptp(m,2,k) + tkxy*uoptp(j,2,i)
     &                       + tixz*uoptp(m,3,k) + tkxz*uoptp(j,3,i)
                        depy = tixy*uoptp(m,1,k) + tkxy*uoptp(j,1,i)
     &                       + tiyy*uoptp(m,2,k) + tkyy*uoptp(j,2,i)
     &                       + tiyz*uoptp(m,3,k) + tkyz*uoptp(j,3,i)
                        depz = tixz*uoptp(m,1,k) + tkxz*uoptp(j,1,i)
     &                       + tiyz*uoptp(m,2,k) + tkyz*uoptp(j,2,i)
     &                       + tizz*uoptp(m,3,k) + tkzz*uoptp(j,3,i)
                        frcx = frcx + copm(j+m+1)*wscale(k)*depx
                        frcy = frcy + copm(j+m+1)*wscale(k)*depy
                        frcz = frcz + copm(j+m+1)*wscale(k)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for TCG polarization force
c
               else if (poltyp.eq.'TCG' .and. use_thole) then
                  do j = 1, tcgnab
                     ukx = ubd(1,kk,j)
                     uky = ubd(2,kk,j)
                     ukz = ubd(3,kk,j)
                     ukxp = ubp(1,kk,j)
                     ukyp = ubp(2,kk,j)
                     ukzp = ubp(3,kk,j)
                     uirt = uax(j)*xr + uay(j)*yr + uaz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = (sc3+sc5) * rr5
                     term2 = term1*xr - rc3(1)
                     term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                     tixx = uax(j)*term2 + uirt*term3
                     tkxx = ukx*term2 + ukrt*term3
                     term2 = term1*yr - rc3(2)
                     term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                     tiyy = uay(j)*term2 + uirt*term3
                     tkyy = uky*term2 + ukrt*term3
                     term2 = term1*zr - rc3(3)
                     term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                     tizz = uaz(j)*term2 + uirt*term3
                     tkzz = ukz*term2 + ukrt*term3
                     term1 = sc5 * rr5 * yr
                     term2 = sc3*rr5*xr - rc3(1)
                     term3 = yr * (sc5*rr7*xr-rc5(1))
                     tixy = uax(j)*term1 + uay(j)*term2 - uirt*term3
                     tkxy = ukx*term1 + uky*term2 - ukrt*term3
                     term1 = sc5 * rr5 * zr
                     term3 = zr * (sc5*rr7*xr-rc5(1))
                     tixz = uax(j)*term1 + uaz(j)*term2 - uirt*term3
                     tkxz = ukx*term1 + ukz*term2 - ukrt*term3
                     term2 = sc3*rr5*yr - rc3(2)
                     term3 = zr * (sc5*rr7*yr-rc5(2))
                     tiyz = uay(j)*term1 + uaz(j)*term2 - uirt*term3
                     tkyz = uky*term1 + ukz*term2 - ukrt*term3
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*uaxp(j) + tkxy*uayp(j)
     &                         + tkxz*uazp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*uaxp(j) + tkyy*uayp(j)
     &                         + tkyz*uazp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*uaxp(j) + tkyz*uayp(j)
     &                         + tkzz*uazp(j)
                     frcx = frcx + uscale(k)*depx
                     frcy = frcy + uscale(k)*depy
                     frcz = frcz + uscale(k)*depz
                     ukx = uad(1,kk,j)
                     uky = uad(2,kk,j)
                     ukz = uad(3,kk,j)
                     ukxp = uap(1,kk,j)
                     ukyp = uap(2,kk,j)
                     ukzp = uap(3,kk,j)
                     uirt = ubx(j)*xr + uby(j)*yr + ubz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = (sc3+sc5) * rr5
                     term2 = term1*xr - rc3(1)
                     term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                     tixx = ubx(j)*term2 + uirt*term3
                     tkxx = ukx*term2 + ukrt*term3
                     term2 = term1*yr - rc3(2)
                     term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                     tiyy = uby(j)*term2 + uirt*term3
                     tkyy = uky*term2 + ukrt*term3
                     term2 = term1*zr - rc3(3)
                     term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                     tizz = ubz(j)*term2 + uirt*term3
                     tkzz = ukz*term2 + ukrt*term3
                     term1 = sc5 * rr5 * yr
                     term2 = sc3*rr5*xr - rc3(1)
                     term3 = yr * (sc5*rr7*xr-rc5(1))
                     tixy = ubx(j)*term1 + uby(j)*term2 - uirt*term3
                     tkxy = ukx*term1 + uky*term2 - ukrt*term3
                     term1 = sc5 * rr5 * zr
                     term3 = zr * (sc5*rr7*xr-rc5(1))
                     tixz = ubx(j)*term1 + ubz(j)*term2 - uirt*term3
                     tkxz = ukx*term1 + ukz*term2 - ukrt*term3
                     term2 = sc3*rr5*yr - rc3(2)
                     term3 = zr * (sc5*rr7*yr-rc5(2))
                     tiyz = uby(j)*term1 + ubz(j)*term2 - uirt*term3
                     tkyz = uky*term1 + ukz*term2 - ukrt*term3
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*ubxp(j) + tkxy*ubyp(j)
     &                         + tkxz*ubzp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*ubxp(j) + tkyy*ubyp(j)
     &                         + tkyz*ubzp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*ubxp(j) + tkyz*ubyp(j)
     &                         + tkzz*ubzp(j)
                     frcx = frcx + uscale(k)*depx
                     frcy = frcy + uscale(k)*depy
                     frcz = frcz + uscale(k)*depz
                  end do
               end if
c
c     increment force-based gradient on the interaction sites
c
               dep(1,i) = dep(1,i) + frcx
               dep(2,i) = dep(2,i) + frcy
               dep(3,i) = dep(3,i) + frcz
               dep(1,k) = dep(1,k) - frcx
               dep(2,k) = dep(2,k) - frcy
               dep(3,k) = dep(3,k) - frcz
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
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = 1.0d0
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
               uscale(ip14(j,i)) = 1.0d0
            end do
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO reduction(+:dep,vir) schedule(guided)
c
c     torque is induced field and gradient cross permanent moments
c
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         tep(1) = diz*ufld(2,i) - diy*ufld(3,i)
     &               + qixz*dufld(2,i) - qixy*dufld(4,i)
     &               + 2.0d0*qiyz*(dufld(3,i)-dufld(6,i))
     &               + (qizz-qiyy)*dufld(5,i)
         tep(2) = dix*ufld(3,i) - diz*ufld(1,i)
     &               - qiyz*dufld(2,i) + qixy*dufld(5,i)
     &               + 2.0d0*qixz*(dufld(6,i)-dufld(1,i))
     &               + (qixx-qizz)*dufld(4,i)
         tep(3) = diy*ufld(1,i) - dix*ufld(2,i)
     &               + qiyz*dufld(4,i) - qixz*dufld(5,i)
     &               + 2.0d0*qixy*(dufld(1,i)-dufld(3,i))
     &               + (qiyy-qixx)*dufld(2,i)
         call torque (ii,tep,fix,fiy,fiz,dep)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = abs(yaxis(ii))
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
!$OMP    DO reduction(+:dep,vir) schedule(guided)
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            frcx = decfx(i)
            frcy = decfy(i)
            frcz = decfz(i)
            dep(1,i) = dep(1,i) + frcx
            dep(2,i) = dep(2,i) + frcy
            dep(3,i) = dep(3,i) + frcz
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
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (wscale)
      deallocate (ufld)
      deallocate (dufld)
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
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
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use poltcg
      use potent
      use virial
      implicit none
      integer i,j,ii
      integer ix,iy,iz
      real*8 f,term
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xv,yv,zv,vterm
      real*8 xufield,yufield
      real*8 zufield
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 fix(3),fiy(3),fiz(3)
      real*8 tep(3)
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
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bsporder
      aewald = apewald
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
c     compute the total induced dipole polarization energy
c
      call epolar1e
c
c     compute the real space part of the Ewald summation
c
      call epreal1c
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip1
c
c     compute the Ewald self-energy torque and virial terms
c
      term = (4.0d0/3.0d0) * f * aewald**3 / rootpi
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         uix = 0.5d0 * (uind(1,ii)+uinp(1,ii))
         uiy = 0.5d0 * (uind(2,ii)+uinp(2,ii))
         uiz = 0.5d0 * (uind(3,ii)+uinp(3,ii))
         tep(1) = term * (diy*uiz-diz*uiy)
         tep(2) = term * (diz*uix-dix*uiz)
         tep(3) = term * (dix*uiy-diy*uix)
         call torque (ii,tep,fix,fiy,fiz,dep)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = abs(yaxis(ii))
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
     &                     + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                     + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3))
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                     + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
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
         do ii = 1, npole
            i = ipole(ii)
            xd = xd + rpole(2,ii) + rpole(1,ii)*x(i)
            yd = yd + rpole(3,ii) + rpole(1,ii)*y(i)
            zd = zd + rpole(4,ii) + rpole(1,ii)*z(i)
            xu = xu + uind(1,ii)
            yu = yu + uind(2,ii)
            zu = zu + uind(3,ii)
            xup = xup + uinp(1,ii)
            yup = yup + uinp(2,ii)
            zup = zup + uinp(3,ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
         do ii = 1, npole
            i = ipole(ii)
            dep(1,i) = dep(1,i) + term*rpole(1,ii)*(xu+xup)
            dep(2,i) = dep(2,i) + term*rpole(1,ii)*(yu+yup)
            dep(3,i) = dep(3,i) + term*rpole(1,ii)*(zu+zup)
         end do
         xufield = -term * (xu+xup)
         yufield = -term * (yu+yup)
         zufield = -term * (zu+zup)
         do ii = 1, npole
            tep(1) = rpole(3,ii)*zufield - rpole(4,ii)*yufield
            tep(2) = rpole(4,ii)*xufield - rpole(2,ii)*zufield
            tep(3) = rpole(2,ii)*yufield - rpole(3,ii)*xufield
            call torque (ii,tep,fix,fiy,fiz,dep)
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
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use deriv
      use ewald
      use math
      use mplpot
      use molcul
      use mpole
      use polar
      use polgrp
      use polopt
      use polpot
      use poltcg
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k,m
      integer ii,kk,jcell
      integer ix,iy,iz
      real*8 f,pgamma
      real*8 pdi,pti,ddi
      real*8 damp,expdamp
      real*8 temp3,temp5,temp7
      real*8 sc3,sc5,sc7
      real*8 psc3,psc5,psc7
      real*8 dsc3,dsc5,dsc7
      real*8 usc3,usc5
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 usr3,usr5
      real*8 rr3core,rr5core
      real*8 rr3i,rr5i
      real*8 rr7i,rr9i
      real*8 rr3k,rr5k
      real*8 rr7k,rr9k
      real*8 rr5ik,rr7ik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 uixp,uiyp,uizp
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 ukxp,ukyp,ukzp
      real*8 dir,uir,uirp
      real*8 dkr,ukr,ukrp
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 uirm,ukrm
      real*8 uirt,ukrt
      real*8 tuir,tukr
      real*8 tixx,tiyy,tizz
      real*8 tixy,tixz,tiyz
      real*8 tkxx,tkyy,tkzz
      real*8 tkxy,tkxz,tkyz
      real*8 tix3,tiy3,tiz3
      real*8 tix5,tiy5,tiz5
      real*8 tkx3,tky3,tkz3
      real*8 tkx5,tky5,tkz5
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term6,term7
      real*8 term1core
      real*8 term1i,term2i,term3i
      real*8 term4i,term5i,term6i
      real*8 term7i,term8i
      real*8 term1k,term2k,term3k
      real*8 term4k,term5k,term6k
      real*8 term7k,term8k
      real*8 poti,potk
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
      real*8 urc3(3),urc5(3),tep(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 uax(3),uay(3),uaz(3)
      real*8 ubx(3),uby(3),ubz(3)
      real*8 uaxp(3),uayp(3),uazp(3)
      real*8 ubxp(3),ubyp(3),ubzp(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9),dmpe(9)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (wscale(n))
      allocate (ufld(3,n))
      allocate (dufld(6,n))
      allocate (pot(n))
      allocate (decfx(n))
      allocate (decfy(n))
      allocate (decfz(n))
c
c     set exclusion coefficients and arrays to store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
         wscale(i) = 1.0d0
         do j = 1, 3
            ufld(j,i) = 0.0d0
         end do
         do j = 1, 6
            dufld(j,i) = 0.0d0
         end do
         pot(i) = 0.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the dipole polarization gradient components
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
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         uixp = uinp(1,ii)
         uiyp = uinp(2,ii)
         uizp = uinp(3,ii)
         do j = 1, tcgnab
            uax(j) = uad(1,ii,j)
            uay(j) = uad(2,ii,j)
            uaz(j) = uad(3,ii,j)
            uaxp(j) = uap(1,ii,j)
            uayp(j) = uap(2,ii,j)
            uazp(j) = uap(3,ii,j)
            ubx(j) = ubd(1,ii,j)
            uby(j) = ubd(2,ii,j)
            ubz(j) = ubd(3,ii,j)
            ubxp(j) = ubp(1,ii,j)
            ubyp(j) = ubp(2,ii,j)
            ubzp(j) = ubp(3,ii,j)
         end do
         if (use_thole) then
            pdi = pdamp(ii)
            pti = thole(ii)
            ddi = dirdamp(ii)
         else if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               dscale(i12(j,i)) = pscale(i12(j,i))
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = u4scale
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
               uscale(ip14(j,i)) = u4scale
            end do
         end if
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
               ukx = uind(1,kk)
               uky = uind(2,kk)
               ukz = uind(3,kk)
               ukxp = uinp(1,kk)
               ukyp = uinp(2,kk)
               ukzp = uinp(3,kk)
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
               uir = uix*xr + uiy*yr + uiz*zr
               uirp = uixp*xr + uiyp*yr + uizp*zr
               ukr = ukx*xr + uky*yr + ukz*zr
               ukrp = ukxp*xr + ukyp*yr + ukzp*zr
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate real space Ewald error function damping
c
               call dampewald (9,r,r2,f,dmpe)
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
c
c     apply Thole polarization damping to scale factors
c
               if (use_thole) then
                  damp = pdi * pdamp(kk)
                  if (use_dirdamp) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(ddi,dirdamp(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           sc3 = 1.0d0 - expdamp 
                           sc5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                      +0.15d0*damp**2)
                           temp3 = 1.5d0 * damp * expdamp / r2
                           temp5 = 0.5d0 * (1.0d0+damp)
                           temp7 = 0.7d0 + 0.15d0*damp**2/temp5
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
                  else
                     pgamma = min(pti,thole(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(pti,thole(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           sc3 = 1.0d0 - expdamp
                           sc5 = 1.0d0 - (1.0d0+damp)*expdamp
                           sc7 = 1.0d0 - (1.0d0+damp+0.6d0*damp**2)
     &                                          *expdamp
                           temp3 = 3.0d0 * damp * expdamp / r2
                           temp5 = damp
                           temp7 = -0.2d0 + 0.6d0*damp
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
                  end if
                  psc3 = 1.0d0 - sc3*pscale(k)
                  psc5 = 1.0d0 - sc5*pscale(k)
                  psc7 = 1.0d0 - sc7*pscale(k)
                  dsc3 = 1.0d0 - sc3*dscale(k)
                  dsc5 = 1.0d0 - sc5*dscale(k)
                  dsc7 = 1.0d0 - sc7*dscale(k)
                  usc3 = 1.0d0 - sc3*uscale(k)
                  usc5 = 1.0d0 - sc5*uscale(k)
                  psr3 = dmpe(3) - psc3*rr3
                  psr5 = dmpe(5) - psc5*rr5
                  psr7 = dmpe(7) - psc7*rr7
                  dsr3 = dmpe(3) - dsc3*rr3
                  dsr5 = dmpe(5) - dsc5*rr5
                  dsr7 = dmpe(7) - dsc7*rr7
                  usr3 = dmpe(3) - usc3*rr3
                  usr5 = dmpe(5) - usc5*rr5
                  do j = 1, 3
                     prc3(j) = rc3(j) * pscale(k)
                     prc5(j) = rc5(j) * pscale(k)
                     prc7(j) = rc7(j) * pscale(k)
                     drc3(j) = rc3(j) * dscale(k)
                     drc5(j) = rc5(j) * dscale(k)
                     drc7(j) = rc7(j) * dscale(k)
                     urc3(j) = rc3(j) * uscale(k)
                     urc5(j) = rc5(j) * uscale(k)
                  end do
c
c     apply charge penetration damping to scale factors
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call damppole (r,9,alphai,alphak,dmpi,dmpk,dmpik)
                  rr3core = dmpe(3) - (1.0d0-dscale(k))*rr3
                  rr5core = dmpe(5) - (1.0d0-dscale(k))*rr5
                  rr3i = dmpe(3) - (1.0d0-dscale(k)*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-dscale(k)*dmpi(5))*rr5
                  rr7i = dmpe(7) - (1.0d0-dscale(k)*dmpi(7))*rr7
                  rr9i = dmpe(9) - (1.0d0-dscale(k)*dmpi(9))*rr9
                  rr3k = dmpe(3) - (1.0d0-dscale(k)*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-dscale(k)*dmpk(5))*rr5
                  rr7k = dmpe(7) - (1.0d0-dscale(k)*dmpk(7))*rr7
                  rr9k = dmpe(9) - (1.0d0-dscale(k)*dmpk(9))*rr9
                  rr5ik = dmpe(5) - (1.0d0-wscale(k)*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-wscale(k)*dmpik(7))*rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_thole) then
                     poti = -ukr*psr3 - ukrp*dsr3
                     potk = uir*psr3 + uirp*dsr3
                  else if (use_chgpen) then
                     poti = -2.0d0 * ukr * rr3i
                     potk = 2.0d0 * uir * rr3k
                  end if
                  pot(i) = pot(i) + poti 
                  pot(k) = pot(k) + potk 
               end if
c
c     get the induced dipole field used for dipole torques
c
               if (use_thole) then
                  tix3 = psr3*ukx + dsr3*ukxp
                  tiy3 = psr3*uky + dsr3*ukyp
                  tiz3 = psr3*ukz + dsr3*ukzp
                  tkx3 = psr3*uix + dsr3*uixp
                  tky3 = psr3*uiy + dsr3*uiyp
                  tkz3 = psr3*uiz + dsr3*uizp
                  tuir = -psr5*ukr - dsr5*ukrp
                  tukr = -psr5*uir - dsr5*uirp
               else if (use_chgpen) then
                  tix3 = 2.0d0*rr3i*ukx
                  tiy3 = 2.0d0*rr3i*uky
                  tiz3 = 2.0d0*rr3i*ukz
                  tkx3 = 2.0d0*rr3k*uix
                  tky3 = 2.0d0*rr3k*uiy
                  tkz3 = 2.0d0*rr3k*uiz
                  tuir = -2.0d0*rr5i*ukr
                  tukr = -2.0d0*rr5k*uir
               end if
               ufld(1,i) = ufld(1,i) + tix3 + xr*tuir
               ufld(2,i) = ufld(2,i) + tiy3 + yr*tuir
               ufld(3,i) = ufld(3,i) + tiz3 + zr*tuir
               ufld(1,k) = ufld(1,k) + tkx3 + xr*tukr
               ufld(2,k) = ufld(2,k) + tky3 + yr*tukr
               ufld(3,k) = ufld(3,k) + tkz3 + zr*tukr
c
c     get induced dipole field gradient used for quadrupole torques
c
               if (use_thole) then
                  tix5 = 2.0d0 * (psr5*ukx+dsr5*ukxp)
                  tiy5 = 2.0d0 * (psr5*uky+dsr5*ukyp)
                  tiz5 = 2.0d0 * (psr5*ukz+dsr5*ukzp)
                  tkx5 = 2.0d0 * (psr5*uix+dsr5*uixp)
                  tky5 = 2.0d0 * (psr5*uiy+dsr5*uiyp)
                  tkz5 = 2.0d0 * (psr5*uiz+dsr5*uizp)
                  tuir = -psr7*ukr - dsr7*ukrp
                  tukr = -psr7*uir - dsr7*uirp
               else if (use_chgpen) then
                  tix5 = 4.0d0 * (rr5i*ukx)
                  tiy5 = 4.0d0 * (rr5i*uky)
                  tiz5 = 4.0d0 * (rr5i*ukz)
                  tkx5 = 4.0d0 * (rr5k*uix)
                  tky5 = 4.0d0 * (rr5k*uiy)
                  tkz5 = 4.0d0 * (rr5k*uiz)
                  tuir = -2.0d0*rr7i*ukr 
                  tukr = -2.0d0*rr7k*uir 
               end if
               dufld(1,i) = dufld(1,i) + xr*tix5 + xr*xr*tuir
               dufld(2,i) = dufld(2,i) + xr*tiy5 + yr*tix5
     &                         + 2.0d0*xr*yr*tuir
               dufld(3,i) = dufld(3,i) + yr*tiy5 + yr*yr*tuir
               dufld(4,i) = dufld(4,i) + xr*tiz5 + zr*tix5
     &                         + 2.0d0*xr*zr*tuir
               dufld(5,i) = dufld(5,i) + yr*tiz5 + zr*tiy5
     &                         + 2.0d0*yr*zr*tuir
               dufld(6,i) = dufld(6,i) + zr*tiz5 + zr*zr*tuir
               dufld(1,k) = dufld(1,k) - xr*tkx5 - xr*xr*tukr
               dufld(2,k) = dufld(2,k) - xr*tky5 - yr*tkx5
     &                         - 2.0d0*xr*yr*tukr
               dufld(3,k) = dufld(3,k) - yr*tky5 - yr*yr*tukr
               dufld(4,k) = dufld(4,k) - xr*tkz5 - zr*tkx5
     &                         - 2.0d0*xr*zr*tukr
               dufld(5,k) = dufld(5,k) - yr*tkz5 - zr*tky5
     &                         - 2.0d0*yr*zr*tukr
               dufld(6,k) = dufld(6,k) - zr*tkz5 - zr*zr*tukr
c
c     get the dEd/dR terms used for direct polarization force
c
               if (use_thole) then
                  term1 = dmpe(5) - dsc3*rr5
                  term2 = dmpe(7) - dsc5*rr7
                  term3 = -dsr3 + term1*xr*xr - rr3*xr*drc3(1)
                  term4 = rr3*drc3(1) - term1*xr - dsr5*xr
                  term5 = term2*xr*xr - dsr5 - rr5*xr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*xr - dmpe(7)
     &                       - rr7*xr*drc7(1)
                  term7 = rr5*drc5(1) - 2.0d0*dmpe(7)*xr
     &                       + (dsc5+1.5d0*dsc7)*rr7*xr
                  tixx = ci*term3 + dix*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixx + (qiy*yr+qiz*zr)*dsc7*rr7
     &                      + 2.0d0*qix*term7 + qir*term6
                  tkxx = ck*term3 - dkx*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxx + (qky*yr+qkz*zr)*dsc7*rr7
     &                      + 2.0d0*qkx*term7 + qkr*term6
                  term3 = -dsr3 + term1*yr*yr - rr3*yr*drc3(2)
                  term4 = rr3*drc3(2) - term1*yr - dsr5*yr
                  term5 = term2*yr*yr - dsr5 - rr5*yr*drc5(2)
                  term6 = (dmpe(9)-dsc7*rr9)*yr*yr - dmpe(7)
     &                       - rr7*yr*drc7(2)
                  term7 = rr5*drc5(2) - 2.0d0*dmpe(7)*yr
     &                       + (dsc5+1.5d0*dsc7)*rr7*yr
                  tiyy = ci*term3 + diy*term4 + dir*term5
     &                      + 2.0d0*dsr5*qiyy + (qix*xr+qiz*zr)*dsc7*rr7
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkyy = ck*term3 - dky*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkyy + (qkx*xr+qkz*zr)*dsc7*rr7
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = -dsr3 + term1*zr*zr - rr3*zr*drc3(3)
                  term4 = rr3*drc3(3) - term1*zr - dsr5*zr
                  term5 = term2*zr*zr - dsr5 - rr5*zr*drc5(3)
                  term6 = (dmpe(9)-dsc7*rr9)*zr*zr - dmpe(7)
     &                       - rr7*zr*drc7(3)
                  term7 = rr5*drc5(3) - 2.0d0*dmpe(7)*zr
     &                       + (dsc5+1.5d0*dsc7)*rr7*zr
                  tizz = ci*term3 + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qizz + (qix*xr+qiy*yr)*dsc7*rr7
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkzz = ck*term3 - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkzz + (qkx*xr+qky*yr)*dsc7*rr7
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*xr*yr - rr3*yr*drc3(1)
                  term4 = rr3*drc3(1) - term1*xr
                  term5 = term2*xr*yr - rr5*yr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*yr - rr7*yr*drc7(1)
                  term7 = rr5*drc5(1) - term2*xr
                  tixy = ci*term3 - dsr5*dix*yr + diy*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixy - 2.0d0*dsr7*yr*qix
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkxy = ck*term3 + dsr5*dkx*yr - dky*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxy - 2.0d0*dsr7*yr*qkx
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = term1*xr*zr - rr3*zr*drc3(1)
                  term5 = term2*xr*zr - rr5*zr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*zr - rr7*zr*drc7(1)
                  tixz = ci*term3 - dsr5*dix*zr + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixz - 2.0d0*dsr7*zr*qix
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkxz = ck*term3 + dsr5*dkx*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxz - 2.0d0*dsr7*zr*qkx
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*yr*zr - rr3*zr*drc3(2)
                  term4 = rr3*drc3(2) - term1*yr
                  term5 = term2*yr*zr - rr5*zr*drc5(2)
                  term6 = (dmpe(9)-dsc7*rr9)*yr*zr - rr7*zr*drc7(2)
                  term7 = rr5*drc5(2) - term2*yr
                  tiyz = ci*term3 - dsr5*diy*zr + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qiyz - 2.0d0*dsr7*zr*qiy
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkyz = ck*term3 + dsr5*dky*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkyz - 2.0d0*dsr7*zr*qky
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      - tkxx*uixp - tkxy*uiyp - tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      - tkxy*uixp - tkyy*uiyp - tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      - tkxz*uixp - tkyz*uiyp - tkzz*uizp
                  frcx = depx
                  frcy = depy
                  frcz = depz
c
c     get the dEp/dR terms used for direct polarization force
c
                  term1 = dmpe(5) - psc3*rr5
                  term2 = dmpe(7) - psc5*rr7
                  term3 = -psr3 + term1*xr*xr - rr3*xr*prc3(1)
                  term4 = rr3*prc3(1) - term1*xr - psr5*xr
                  term5 = term2*xr*xr - psr5 - rr5*xr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*xr - dmpe(7)
     &                       - rr7*xr*prc7(1)
                  term7 = rr5*prc5(1) - 2.0d0*dmpe(7)*xr
     &                       + (psc5+1.5d0*psc7)*rr7*xr
                  tixx = ci*term3 + dix*term4 + dir*term5
     &                      + 2.0d0*psr5*qixx + (qiy*yr+qiz*zr)*psc7*rr7
     &                      + 2.0d0*qix*term7 + qir*term6
                  tkxx = ck*term3 - dkx*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxx + (qky*yr+qkz*zr)*psc7*rr7
     &                      + 2.0d0*qkx*term7 + qkr*term6
                  term3 = -psr3 + term1*yr*yr - rr3*yr*prc3(2)
                  term4 = rr3*prc3(2) - term1*yr - psr5*yr
                  term5 = term2*yr*yr - psr5 - rr5*yr*prc5(2)
                  term6 = (dmpe(9)-psc7*rr9)*yr*yr - dmpe(7)
     &                       - rr7*yr*prc7(2)
                  term7 = rr5*prc5(2) - 2.0d0*dmpe(7)*yr
     &                       + (psc5+1.5d0*psc7)*rr7*yr
                  tiyy = ci*term3 + diy*term4 + dir*term5
     &                      + 2.0d0*psr5*qiyy + (qix*xr+qiz*zr)*psc7*rr7
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkyy = ck*term3 - dky*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkyy + (qkx*xr+qkz*zr)*psc7*rr7
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = -psr3 + term1*zr*zr - rr3*zr*prc3(3)
                  term4 = rr3*prc3(3) - term1*zr - psr5*zr
                  term5 = term2*zr*zr - psr5 - rr5*zr*prc5(3)
                  term6 = (dmpe(9)-psc7*rr9)*zr*zr - dmpe(7)
     &                       - rr7*zr*prc7(3)
                  term7 = rr5*prc5(3) - 2.0d0*dmpe(7)*zr
     &                       + (psc5+1.5d0*psc7)*rr7*zr
                  tizz = ci*term3 + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qizz + (qix*xr+qiy*yr)*psc7*rr7
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkzz = ck*term3 - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkzz + (qkx*xr+qky*yr)*psc7*rr7
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*xr*yr - rr3*yr*prc3(1)
                  term4 = rr3*prc3(1) - term1*xr
                  term5 = term2*xr*yr - rr5*yr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*yr - rr7*yr*prc7(1)
                  term7 = rr5*prc5(1) - term2*xr
                  tixy = ci*term3 - psr5*dix*yr + diy*term4 + dir*term5
     &                      + 2.0d0*psr5*qixy - 2.0d0*psr7*yr*qix
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkxy = ck*term3 + psr5*dkx*yr - dky*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxy - 2.0d0*psr7*yr*qkx
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = term1*xr*zr - rr3*zr*prc3(1)
                  term5 = term2*xr*zr - rr5*zr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*zr - rr7*zr*prc7(1)
                  tixz = ci*term3 - psr5*dix*zr + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qixz - 2.0d0*psr7*zr*qix
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkxz = ck*term3 + psr5*dkx*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxz - 2.0d0*psr7*zr*qkx
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*yr*zr - rr3*zr*prc3(2)
                  term4 = rr3*prc3(2) - term1*yr
                  term5 = term2*yr*zr - rr5*zr*prc5(2)
                  term6 = (dmpe(9)-psc7*rr9)*yr*zr - rr7*zr*prc7(2)
                  term7 = rr5*prc5(2) - term2*yr
                  tiyz = ci*term3 - psr5*diy*zr + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qiyz - 2.0d0*psr7*zr*qiy
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkyz = ck*term3 + psr5*dky*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkyz - 2.0d0*psr7*zr*qky
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
c
c     get the field gradient for direct polarization force
c
               else if (use_chgpen) then
                  term1i = rr3i - rr5i*xr*xr
                  term1core = rr3core - rr5core*xr*xr
                  term2i = 2.0d0*rr5i*xr 
                  term3i = rr7i*xr*xr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*xr
                  term6i = rr9i*xr*xr
                  term1k = rr3k - rr5k*xr*xr
                  term2k = 2.0d0*rr5k*xr
                  term3k = rr7k*xr*xr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*xr
                  term6k = rr9k*xr*xr
                  tixx = vali*term1i + corei*term1core  
     &                      + dix*term2i - dir*term3i
     &                      - qixx*term4i + qix*term5i - qir*term6i
     &                      + (qiy*yr+qiz*zr)*rr7i
                  tkxx = valk*term1k + corek*term1core
     &                      - dkx*term2k + dkr*term3k
     &                      - qkxx*term4k + qkx*term5k - qkr*term6k
     &                      + (qky*yr+qkz*zr)*rr7k
                  term1i = rr3i - rr5i*yr*yr
                  term1core = rr3core - rr5core*yr*yr
                  term2i = 2.0d0*rr5i*yr
                  term3i = rr7i*yr*yr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*yr
                  term6i = rr9i*yr*yr
                  term1k = rr3k - rr5k*yr*yr
                  term2k = 2.0d0*rr5k*yr
                  term3k = rr7k*yr*yr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*yr
                  term6k = rr9k*yr*yr
                  tiyy = vali*term1i + corei*term1core
     &                      + diy*term2i - dir*term3i
     &                      - qiyy*term4i + qiy*term5i - qir*term6i
     &                      + (qix*xr+qiz*zr)*rr7i
                  tkyy = valk*term1k + corek*term1core
     &                      - dky*term2k + dkr*term3k
     &                      - qkyy*term4k + qky*term5k - qkr*term6k
     &                      + (qkx*xr+qkz*zr)*rr7k
                  term1i = rr3i - rr5i*zr*zr
                  term1core = rr3core - rr5core*zr*zr
                  term2i = 2.0d0*rr5i*zr
                  term3i = rr7i*zr*zr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*zr
                  term6i = rr9i*zr*zr
                  term1k = rr3k - rr5k*zr*zr
                  term2k = 2.0d0*rr5k*zr
                  term3k = rr7k*zr*zr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*zr
                  term6k = rr9k*zr*zr
                  tizz = vali*term1i + corei*term1core
     &                      + diz*term2i - dir*term3i
     &                      - qizz*term4i + qiz*term5i - qir*term6i
     &                      + (qix*xr+qiy*yr)*rr7i
                  tkzz = valk*term1k + corek*term1core
     &                      - dkz*term2k + dkr*term3k
     &                      - qkzz*term4k + qkz*term5k - qkr*term6k
     &                      + (qkx*xr+qky*yr)*rr7k
                  term2i = rr5i*xr 
                  term1i = yr * term2i
                  term1core = rr5core*xr*yr
                  term3i = rr5i*yr
                  term4i = yr * (rr7i*xr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*xr
                  term7i = 2.0d0*rr7i*yr
                  term8i = yr*rr9i*xr
                  term2k = rr5k*xr
                  term1k = yr * term2k
                  term3k = rr5k*yr
                  term4k = yr * (rr7k*xr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*xr
                  term7k = 2.0d0*rr7k*yr
                  term8k = yr*rr9k*xr
                  tixy = -vali*term1i - corei*term1core 
     &                      + diy*term2i + dix*term3i
     &                      - dir*term4i - qixy*term5i + qiy*term6i
     &                      + qix*term7i - qir*term8i
                  tkxy = -valk*term1k - corek*term1core 
     &                      - dky*term2k - dkx*term3k
     &                      + dkr*term4k - qkxy*term5k + qky*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5i*xr
                  term1i = zr * term2i
                  term1core = rr5core*xr*zr
                  term3i = rr5i*zr
                  term4i = zr * (rr7i*xr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*xr
                  term7i = 2.0d0*rr7i*zr
                  term8i = zr*rr9i*xr
                  term2k = rr5k*xr
                  term1k = zr * term2k
                  term3k = rr5k*zr
                  term4k = zr * (rr7k*xr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*xr
                  term7k = 2.0d0*rr7k*zr
                  term8k = zr*rr9k*xr
                  tixz = -vali*term1i - corei*term1core
     &                      + diz*term2i + dix*term3i
     &                      - dir*term4i - qixz*term5i + qiz*term6i
     &                      + qix*term7i - qir*term8i
                  tkxz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dkx*term3k
     &                      + dkr*term4k - qkxz*term5k + qkz*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5i*yr
                  term1i = zr * term2i
                  term1core = rr5core*yr*zr
                  term3i = rr5i*zr
                  term4i = zr * (rr7i*yr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*yr
                  term7i = 2.0d0*rr7i*zr
                  term8i = zr*rr9i*yr
                  term2k = rr5k*yr
                  term1k = zr * term2k
                  term3k = rr5k*zr
                  term4k = zr * (rr7k*yr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*yr
                  term7k = 2.0d0*rr7k*zr
                  term8k = zr*rr9k*yr
                  tiyz = -vali*term1i - corei*term1core
     &                      + diz*term2i + diy*term3i
     &                      - dir*term4i - qiyz*term5i + qiz*term6i
     &                      + qiy*term7i - qir*term8i
                  tkyz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dky*term3k
     &                      + dkr*term4k - qkyz*term5k + qkz*term6k
     &                      + qky*term7k - qkr*term8k
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = -2.0d0 * depx
                  frcy = -2.0d0 * depy
                  frcz = -2.0d0 * depz
               end if
c
c     reset Thole values if alternate direct damping was used
c
               if (use_dirdamp) then
                  sc3 = 1.0d0
                  sc5 = 1.0d0
                  do j = 1, 3
                     rc3(j) = 0.0d0
                     rc5(j) = 0.0d0
                  end do
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        sc3 = 1.0d0 - expdamp
                        sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                        temp3 = 3.0d0 * damp * expdamp / r2
                        temp5 = damp
                        rc3(1) = xr * temp3
                        rc3(2) = yr * temp3
                        rc3(3) = zr * temp3
                        rc5(1) = rc3(1) * temp5
                        rc5(2) = rc3(2) * temp5
                        rc5(3) = rc3(3) * temp5
                     end if
                  end if
                  usc3 = 1.0d0 - sc3*uscale(k)
                  usc5 = 1.0d0 - sc5*uscale(k)
                  usr3 = dmpe(3) - usc3*rr3
                  usr5 = dmpe(5) - usc5*rr5
                  do j = 1, 3
                     urc3(j) = rc3(j) * uscale(k)
                     urc5(j) = rc5(j) * uscale(k)
                  end do
               end if
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp.eq.'MUTUAL' .and. use_thole) then
                  term1 = dmpe(5) - usc3*rr5
                  term2 = dmpe(7) - usc5*rr7
                  term3 = usr5 + term1
                  term4 = rr3 * uscale(k)
                  term5 = -xr*term3 + rc3(1)*term4
                  term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                  tixx = uix*term5 + uir*term6
                  tkxx = ukx*term5 + ukr*term6
                  term5 = -yr*term3 + rc3(2)*term4
                  term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                  tiyy = uiy*term5 + uir*term6
                  tkyy = uky*term5 + ukr*term6
                  term5 = -zr*term3 + rc3(3)*term4
                  term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                  tizz = uiz*term5 + uir*term6
                  tkzz = ukz*term5 + ukr*term6
                  term4 = -usr5 * yr
                  term5 = -xr*term1 + rr3*urc3(1)
                  term6 = xr*yr*term2 - rr5*yr*urc5(1)
                  tixy = uix*term4 + uiy*term5 + uir*term6
                  tkxy = ukx*term4 + uky*term5 + ukr*term6
                  term4 = -usr5 * zr
                  term6 = xr*zr*term2 - rr5*zr*urc5(1)
                  tixz = uix*term4 + uiz*term5 + uir*term6
                  tkxz = ukx*term4 + ukz*term5 + ukr*term6
                  term5 = -yr*term1 + rr3*urc3(2)
                  term6 = yr*zr*term2 - rr5*zr*urc5(2)
                  tiyz = uiy*term4 + uiz*term5 + uir*term6
                  tkyz = uky*term4 + ukz*term5 + ukr*term6
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               else if (poltyp.eq.'MUTUAL' .and. use_chgpen) then
                  term1 = 2.0d0 * rr5ik
                  term2 = term1*xr
                  term3 = rr5ik - rr7ik*xr*xr
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr
                  term3 = rr5ik - rr7ik*yr*yr
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr
                  term3 = rr5ik - rr7ik*zr*zr
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = rr5ik*yr
                  term2 = rr5ik*xr
                  term3 = yr * (rr7ik*xr)
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = rr5ik * zr
                  term3 = zr * (rr7ik*xr)
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = rr5ik*yr
                  term3 = zr * (rr7ik*yr)
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx - depx
                  frcy = frcy - depy
                  frcz = frcz - depz
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_thole) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,ii)*xr + uopt(j,2,ii)*yr
     &                          + uopt(j,3,ii)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,kk)*xr + uopt(m,2,kk)*yr
     &                             + uopt(m,3,kk)*zr
                        term1 = dmpe(5) - usc3*rr5
                        term2 = dmpe(7) - usc5*rr7
                        term3 = usr5 + term1
                        term4 = rr3 * uscale(k)
                        term5 = -xr*term3 + rc3(1)*term4
                        term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                        tixx = uopt(j,1,ii)*term5 + uirm*term6
                        tkxx = uopt(m,1,kk)*term5 + ukrm*term6
                        term5 = -yr*term3 + rc3(2)*term4
                        term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                        tiyy = uopt(j,2,ii)*term5 + uirm*term6
                        tkyy = uopt(m,2,kk)*term5 + ukrm*term6
                        term5 = -zr*term3 + rc3(3)*term4
                        term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                        tizz = uopt(j,3,ii)*term5 + uirm*term6
                        tkzz = uopt(m,3,kk)*term5 + ukrm*term6
                        term4 = -usr5 * yr
                        term5 = -xr*term1 + rr3*urc3(1)
                        term6 = xr*yr*term2 - rr5*yr*urc5(1)
                        tixy = uopt(j,1,ii)*term4 + uopt(j,2,ii)*term5
     &                            + uirm*term6
                        tkxy = uopt(m,1,kk)*term4 + uopt(m,2,kk)*term5
     &                            + ukrm*term6
                        term4 = -usr5 * zr
                        term6 = xr*zr*term2 - rr5*zr*urc5(1)
                        tixz = uopt(j,1,ii)*term4 + uopt(j,3,ii)*term5
     &                            + uirm*term6
                        tkxz = uopt(m,1,kk)*term4 + uopt(m,3,kk)*term5
     &                            + ukrm*term6
                        term5 = -yr*term1 + rr3*urc3(2)
                        term6 = yr*zr*term2 - rr5*zr*urc5(2)
                        tiyz = uopt(j,2,ii)*term4 + uopt(j,3,ii)*term5
     &                            + uirm*term6
                        tkyz = uopt(m,2,kk)*term4 + uopt(m,3,kk)*term5
     &                            + ukrm*term6
                        depx = tixx*uoptp(m,1,kk) + tkxx*uoptp(j,1,ii)
     &                       + tixy*uoptp(m,2,kk) + tkxy*uoptp(j,2,ii)
     &                       + tixz*uoptp(m,3,kk) + tkxz*uoptp(j,3,ii)
                        depy = tixy*uoptp(m,1,kk) + tkxy*uoptp(j,1,ii)
     &                       + tiyy*uoptp(m,2,kk) + tkyy*uoptp(j,2,ii)
     &                       + tiyz*uoptp(m,3,kk) + tkyz*uoptp(j,3,ii)
                        depz = tixz*uoptp(m,1,kk) + tkxz*uoptp(j,1,ii)
     &                       + tiyz*uoptp(m,2,kk) + tkyz*uoptp(j,2,ii)
     &                       + tizz*uoptp(m,3,kk) + tkzz*uoptp(j,3,ii)
                        frcx = frcx + copm(j+m+1)*depx
                        frcy = frcy + copm(j+m+1)*depy
                        frcz = frcz + copm(j+m+1)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_chgpen) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,i)*xr + uopt(j,2,i)*yr
     &                          + uopt(j,3,i)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,k)*xr + uopt(m,2,k)*yr
     &                             + uopt(m,3,k)*zr
                        term1 = 2.0d0 * rr5ik
                        term2 = term1*xr
                        term3 = rr5ik - rr7ik*xr*xr
                        tixx = uopt(j,1,i)*term2 + uirm*term3
                        tkxx = uopt(m,1,k)*term2 + ukrm*term3
                        term2 = term1*yr
                        term3 = rr5ik - rr7ik*yr*yr
                        tiyy = uopt(j,2,i)*term2 + uirm*term3
                        tkyy = uopt(m,2,k)*term2 + ukrm*term3
                        term2 = term1*zr
                        term3 = rr5ik - rr7ik*zr*zr
                        tizz = uopt(j,3,i)*term2 + uirm*term3
                        tkzz = uopt(m,3,k)*term2 + ukrm*term3
                        term1 = rr5ik*yr
                        term2 = rr5ik*xr
                        term3 = yr * (rr7ik*xr)
                        tixy = uopt(j,1,i)*term1 + uopt(j,2,i)*term2 
     &                       - uirm*term3
                        tkxy = uopt(m,1,k)*term1 + uopt(m,2,k)*term2 
     &                       - ukrm*term3
                        term1 = rr5ik * zr
                        term3 = zr * (rr7ik*xr)
                        tixz = uopt(j,1,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        term2 = rr5ik*yr
                        term3 = zr * (rr7ik*yr)
                        tiyz = uopt(j,2,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,k) + tkxx*uoptp(j,1,i)
     &                       + tixy*uoptp(m,2,k) + tkxy*uoptp(j,2,i)
     &                       + tixz*uoptp(m,3,k) + tkxz*uoptp(j,3,i)
                        depy = tixy*uoptp(m,1,k) + tkxy*uoptp(j,1,i)
     &                       + tiyy*uoptp(m,2,k) + tkyy*uoptp(j,2,i)
     &                       + tiyz*uoptp(m,3,k) + tkyz*uoptp(j,3,i)
                        depz = tixz*uoptp(m,1,k) + tkxz*uoptp(j,1,i)
     &                       + tiyz*uoptp(m,2,k) + tkyz*uoptp(j,2,i)
     &                       + tizz*uoptp(m,3,k) + tkzz*uoptp(j,3,i)
                        frcx = frcx - copm(j+m+1)*depx
                        frcy = frcy - copm(j+m+1)*depy
                        frcz = frcz - copm(j+m+1)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for TCG polarization force
c
               else if (poltyp.eq.'TCG' .and. use_thole) then
                  do j = 1, tcgnab
                     ukx = ubd(1,kk,j)
                     uky = ubd(2,kk,j)
                     ukz = ubd(3,kk,j)
                     ukxp = ubp(1,kk,j)
                     ukyp = ubp(2,kk,j)
                     ukzp = ubp(3,kk,j)
                     uirt = uax(j)*xr + uay(j)*yr + uaz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = dmpe(5) - usc3*rr5
                     term2 = dmpe(7) - usc5*rr7
                     term3 = usr5 + term1
                     term4 = rr3 * uscale(k)
                     term5 = -xr*term3 + rc3(1)*term4
                     term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                     tixx = uax(j)*term5 + uirt*term6
                     tkxx = ukx*term5 + ukrt*term6
                     term5 = -yr*term3 + rc3(2)*term4
                     term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                     tiyy = uay(j)*term5 + uirt*term6
                     tkyy = uky*term5 + ukrt*term6
                     term5 = -zr*term3 + rc3(3)*term4
                     term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                     tizz = uaz(j)*term5 + uirt*term6
                     tkzz = ukz*term5 + ukrt*term6
                     term4 = -usr5 * yr
                     term5 = -xr*term1 + rr3*urc3(1)
                     term6 = xr*yr*term2 - rr5*yr*urc5(1)
                     tixy = uax(j)*term4 + uay(j)*term5 + uirt*term6
                     tkxy = ukx*term4 + uky*term5 + ukrt*term6
                     term4 = -usr5 * zr
                     term6 = xr*zr*term2 - rr5*zr*urc5(1)
                     tixz = uax(j)*term4 + uaz(j)*term5 + uirt*term6
                     tkxz = ukx*term4 + ukz*term5 + ukrt*term6
                     term5 = -yr*term1 + rr3*urc3(2)
                     term6 = yr*zr*term2 - rr5*zr*urc5(2)
                     tiyz = uay(j)*term4 + uaz(j)*term5 + uirt*term6
                     tkyz = uky*term4 + ukz*term5 + ukrt*term6
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*uaxp(j) + tkxy*uayp(j)
     &                         + tkxz*uazp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*uaxp(j) + tkyy*uayp(j)
     &                         + tkyz*uazp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*uaxp(j) + tkyz*uayp(j)
     &                         + tkzz*uazp(j)
                     frcx = frcx + depx
                     frcy = frcy + depy
                     frcz = frcz + depz
                     ukx = uad(1,kk,j)
                     uky = uad(2,kk,j)
                     ukz = uad(3,kk,j)
                     ukxp = uap(1,kk,j)
                     ukyp = uap(2,kk,j)
                     ukzp = uap(3,kk,j)
                     uirt = ubx(j)*xr + uby(j)*yr + ubz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = dmpe(5) - usc3*rr5
                     term2 = dmpe(7) - usc5*rr7
                     term3 = usr5 + term1
                     term4 = rr3 * uscale(k)
                     term5 = -xr*term3 + rc3(1)*term4
                     term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                     tixx = ubx(j)*term5 + uirt*term6
                     tkxx = ukx*term5 + ukrt*term6
                     term5 = -yr*term3 + rc3(2)*term4
                     term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                     tiyy = uby(j)*term5 + uirt*term6
                     tkyy = uky*term5 + ukrt*term6
                     term5 = -zr*term3 + rc3(3)*term4
                     term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                     tizz = ubz(j)*term5 + uirt*term6
                     tkzz = ukz*term5 + ukrt*term6
                     term4 = -usr5 * yr
                     term5 = -xr*term1 + rr3*urc3(1)
                     term6 = xr*yr*term2 - rr5*yr*urc5(1)
                     tixy = ubx(j)*term4 + uby(j)*term5 + uirt*term6
                     tkxy = ukx*term4 + uky*term5 + ukrt*term6
                     term4 = -usr5 * zr
                     term6 = xr*zr*term2 - rr5*zr*urc5(1)
                     tixz = ubx(j)*term4 + ubz(j)*term5 + uirt*term6
                     tkxz = ukx*term4 + ukz*term5 + ukrt*term6
                     term5 = -yr*term1 + rr3*urc3(2)
                     term6 = yr*zr*term2 - rr5*zr*urc5(2)
                     tiyz = uby(j)*term4 + ubz(j)*term5 + uirt*term6
                     tkyz = uky*term4 + ukz*term5 + ukrt*term6
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*ubxp(j) + tkxy*ubyp(j)
     &                         + tkxz*ubzp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*ubxp(j) + tkyy*ubyp(j)
     &                         + tkyz*ubzp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*ubxp(j) + tkyz*ubyp(j)
     &                         + tkzz*ubzp(j)
                     frcx = frcx + depx
                     frcy = frcy + depy
                     frcz = frcz + depz
                  end do
               end if
c
c     increment force-based gradient on the interaction sites
c
               dep(1,i) = dep(1,i) - frcx
               dep(2,i) = dep(2,i) - frcy
               dep(3,i) = dep(3,i) - frcz
               dep(1,k) = dep(1,k) + frcx
               dep(2,k) = dep(2,k) + frcy
               dep(3,k) = dep(3,k) + frcz
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = xr * frcx
               vxy = 0.5d0 * (yr*frcx+xr*frcy)
               vxz = 0.5d0 * (zr*frcx+xr*frcz)
               vyy = yr * frcy
               vyz = 0.5d0 * (zr*frcy+yr*frcz)
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
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = 1.0d0
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
               uscale(ip14(j,i)) = 1.0d0
            end do
         end if
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
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         uixp = uinp(1,ii)
         uiyp = uinp(2,ii)
         uizp = uinp(3,ii)
         do j = 1, tcgnab
            uax(j) = uad(1,ii,j)
            uay(j) = uad(2,ii,j)
            uaz(j) = uad(3,ii,j)
            uaxp(j) = uap(1,ii,j)
            uayp(j) = uap(2,ii,j)
            uazp(j) = uap(3,ii,j)
            ubx(j) = ubd(1,ii,j)
            uby(j) = ubd(2,ii,j)
            ubz(j) = ubd(3,ii,j)
            ubxp(j) = ubp(1,ii,j)
            ubyp(j) = ubp(2,ii,j)
            ubzp(j) = ubp(3,ii,j)
         end do
         if (use_thole) then
            pdi = pdamp(ii)
            pti = thole(ii)
            ddi = dirdamp(ii)
         else if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               dscale(i12(j,i)) = pscale(i12(j,i))
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = u4scale
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
               uscale(ip14(j,i)) = u4scale
            end do
         end if
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii, npole
            k = ipole(kk)
            do jcell = 2, ncell
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               pscale(k) = 1.0d0
               dscale(k) = 1.0d0
               uscale(k) = 1.0d0
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
               ukx = uind(1,kk)
               uky = uind(2,kk)
               ukz = uind(3,kk)
               ukxp = uinp(1,kk)
               ukyp = uinp(2,kk)
               ukzp = uinp(3,kk)
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
               uir = uix*xr + uiy*yr + uiz*zr
               uirp = uixp*xr + uiyp*yr + uizp*zr
               ukr = ukx*xr + uky*yr + ukz*zr
               ukrp = ukxp*xr + ukyp*yr + ukzp*zr
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate real space Ewald error function damping
c
               call dampewald (9,r,r2,f,dmpe)
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
c
c     apply Thole polarization damping to scale factors
c
               if (use_thole) then
                  damp = pdi * pdamp(kk)
                  if (use_dirdamp) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(ddi,dirdamp(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           sc3 = 1.0d0 - expdamp 
                           sc5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                      +0.15d0*damp**2)
                           temp3 = 1.5d0 * damp * expdamp / r2
                           temp5 = 0.5d0 * (1.0d0+damp)
                           temp7 = 0.7d0 + 0.15d0*damp**2/temp5
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
                  else
                     pgamma = min(pti,thole(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(pti,thole(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           sc3 = 1.0d0 - expdamp
                           sc5 = 1.0d0 - (1.0d0+damp)*expdamp
                           sc7 = 1.0d0 - (1.0d0+damp+0.6d0*damp**2)
     &                                          *expdamp
                           temp3 = 3.0d0 * damp * expdamp / r2
                           temp5 = damp
                           temp7 = -0.2d0 + 0.6d0*damp
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
                  end if
                  psc3 = 1.0d0 - sc3*pscale(k)
                  psc5 = 1.0d0 - sc5*pscale(k)
                  psc7 = 1.0d0 - sc7*pscale(k)
                  dsc3 = 1.0d0 - sc3*dscale(k)
                  dsc5 = 1.0d0 - sc5*dscale(k)
                  dsc7 = 1.0d0 - sc7*dscale(k)
                  usc3 = 1.0d0 - sc3*uscale(k)
                  usc5 = 1.0d0 - sc5*uscale(k)
                  psr3 = dmpe(3) - psc3*rr3
                  psr5 = dmpe(5) - psc5*rr5
                  psr7 = dmpe(7) - psc7*rr7
                  dsr3 = dmpe(3) - dsc3*rr3
                  dsr5 = dmpe(5) - dsc5*rr5
                  dsr7 = dmpe(7) - dsc7*rr7
                  usr3 = dmpe(3) - usc3*rr3
                  usr5 = dmpe(5) - usc5*rr5
                  do j = 1, 3
                     prc3(j) = rc3(j) * pscale(k)
                     prc5(j) = rc5(j) * pscale(k)
                     prc7(j) = rc7(j) * pscale(k)
                     drc3(j) = rc3(j) * dscale(k)
                     drc5(j) = rc5(j) * dscale(k)
                     drc7(j) = rc7(j) * dscale(k)
                     urc3(j) = rc3(j) * uscale(k)
                     urc5(j) = rc5(j) * uscale(k)
                  end do
c
c     apply charge penetration damping to scale factors
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call damppole (r,9,alphai,alphak,dmpi,dmpk,dmpik)
                  rr3core = dmpe(3) - (1.0d0-dscale(k))*rr3
                  rr5core = dmpe(5) - (1.0d0-dscale(k))*rr5
                  rr3i = dmpe(3) - (1.0d0-dscale(k)*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-dscale(k)*dmpi(5))*rr5
                  rr7i = dmpe(7) - (1.0d0-dscale(k)*dmpi(7))*rr7
                  rr9i = dmpe(9) - (1.0d0-dscale(k)*dmpi(9))*rr9
                  rr3k = dmpe(3) - (1.0d0-dscale(k)*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-dscale(k)*dmpk(5))*rr5
                  rr7k = dmpe(7) - (1.0d0-dscale(k)*dmpk(7))*rr7
                  rr9k = dmpe(9) - (1.0d0-dscale(k)*dmpk(9))*rr9
                  rr5ik = dmpe(5) - (1.0d0-wscale(k)*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-wscale(k)*dmpik(7))*rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_thole) then
                     poti = -ukr*psr3 - ukrp*dsr3
                     potk = uir*psr3 + uirp*dsr3
                  else if (use_chgpen) then
                     poti = -2.0d0 * ukr * rr3i
                     potk = 2.0d0 * uir * rr3k
                  end if
                  pot(i) = pot(i) + poti 
                  pot(k) = pot(k) + potk 
               end if
c
c     get the induced dipole field used for dipole torques
c
               if (use_thole) then
                  tix3 = psr3*ukx + dsr3*ukxp
                  tiy3 = psr3*uky + dsr3*ukyp
                  tiz3 = psr3*ukz + dsr3*ukzp
                  tkx3 = psr3*uix + dsr3*uixp
                  tky3 = psr3*uiy + dsr3*uiyp
                  tkz3 = psr3*uiz + dsr3*uizp
                  tuir = -psr5*ukr - dsr5*ukrp
                  tukr = -psr5*uir - dsr5*uirp
               else if (use_chgpen) then
                  tix3 = 2.0d0*rr3i*ukx
                  tiy3 = 2.0d0*rr3i*uky
                  tiz3 = 2.0d0*rr3i*ukz
                  tkx3 = 2.0d0*rr3k*uix
                  tky3 = 2.0d0*rr3k*uiy
                  tkz3 = 2.0d0*rr3k*uiz
                  tuir = -2.0d0*rr5i*ukr
                  tukr = -2.0d0*rr5k*uir
               end if
               ufld(1,i) = ufld(1,i) + tix3 + xr*tuir
               ufld(2,i) = ufld(2,i) + tiy3 + yr*tuir
               ufld(3,i) = ufld(3,i) + tiz3 + zr*tuir
               ufld(1,k) = ufld(1,k) + tkx3 + xr*tukr
               ufld(2,k) = ufld(2,k) + tky3 + yr*tukr
               ufld(3,k) = ufld(3,k) + tkz3 + zr*tukr
c
c     get induced dipole field gradient used for quadrupole torques
c
               if (use_thole) then
                  tix5 = 2.0d0 * (psr5*ukx+dsr5*ukxp)
                  tiy5 = 2.0d0 * (psr5*uky+dsr5*ukyp)
                  tiz5 = 2.0d0 * (psr5*ukz+dsr5*ukzp)
                  tkx5 = 2.0d0 * (psr5*uix+dsr5*uixp)
                  tky5 = 2.0d0 * (psr5*uiy+dsr5*uiyp)
                  tkz5 = 2.0d0 * (psr5*uiz+dsr5*uizp)
                  tuir = -psr7*ukr - dsr7*ukrp
                  tukr = -psr7*uir - dsr7*uirp
               else if (use_chgpen) then
                  tix5 = 4.0d0 * (rr5i*ukx)
                  tiy5 = 4.0d0 * (rr5i*uky)
                  tiz5 = 4.0d0 * (rr5i*ukz)
                  tkx5 = 4.0d0 * (rr5k*uix)
                  tky5 = 4.0d0 * (rr5k*uiy)
                  tkz5 = 4.0d0 * (rr5k*uiz)
                  tuir = -2.0d0*rr7i*ukr 
                  tukr = -2.0d0*rr7k*uir 
               end if
               dufld(1,i) = dufld(1,i) + xr*tix5 + xr*xr*tuir
               dufld(2,i) = dufld(2,i) + xr*tiy5 + yr*tix5
     &                         + 2.0d0*xr*yr*tuir
               dufld(3,i) = dufld(3,i) + yr*tiy5 + yr*yr*tuir
               dufld(4,i) = dufld(4,i) + xr*tiz5 + zr*tix5
     &                         + 2.0d0*xr*zr*tuir
               dufld(5,i) = dufld(5,i) + yr*tiz5 + zr*tiy5
     &                         + 2.0d0*yr*zr*tuir
               dufld(6,i) = dufld(6,i) + zr*tiz5 + zr*zr*tuir
               dufld(1,k) = dufld(1,k) - xr*tkx5 - xr*xr*tukr
               dufld(2,k) = dufld(2,k) - xr*tky5 - yr*tkx5
     &                         - 2.0d0*xr*yr*tukr
               dufld(3,k) = dufld(3,k) - yr*tky5 - yr*yr*tukr
               dufld(4,k) = dufld(4,k) - xr*tkz5 - zr*tkx5
     &                         - 2.0d0*xr*zr*tukr
               dufld(5,k) = dufld(5,k) - yr*tkz5 - zr*tky5
     &                         - 2.0d0*yr*zr*tukr
               dufld(6,k) = dufld(6,k) - zr*tkz5 - zr*zr*tukr
c
c     get the dEd/dR terms used for direct polarization force
c
               if (use_thole) then
                  term1 = dmpe(5) - dsc3*rr5
                  term2 = dmpe(7) - dsc5*rr7
                  term3 = -dsr3 + term1*xr*xr - rr3*xr*drc3(1)
                  term4 = rr3*drc3(1) - term1*xr - dsr5*xr
                  term5 = term2*xr*xr - dsr5 - rr5*xr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*xr - dmpe(7)
     &                       - rr7*xr*drc7(1)
                  term7 = rr5*drc5(1) - 2.0d0*dmpe(7)*xr
     &                       + (dsc5+1.5d0*dsc7)*rr7*xr
                  tixx = ci*term3 + dix*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixx + (qiy*yr+qiz*zr)*dsc7*rr7
     &                      + 2.0d0*qix*term7 + qir*term6
                  tkxx = ck*term3 - dkx*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxx + (qky*yr+qkz*zr)*dsc7*rr7
     &                      + 2.0d0*qkx*term7 + qkr*term6
                  term3 = -dsr3 + term1*yr*yr - rr3*yr*drc3(2)
                  term4 = rr3*drc3(2) - term1*yr - dsr5*yr
                  term5 = term2*yr*yr - dsr5 - rr5*yr*drc5(2)
                  term6 = (dmpe(9)-dsc7*rr9)*yr*yr - dmpe(7)
     &                       - rr7*yr*drc7(2)
                  term7 = rr5*drc5(2) - 2.0d0*dmpe(7)*yr
     &                       + (dsc5+1.5d0*dsc7)*rr7*yr
                  tiyy = ci*term3 + diy*term4 + dir*term5
     &                      + 2.0d0*dsr5*qiyy + (qix*xr+qiz*zr)*dsc7*rr7
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkyy = ck*term3 - dky*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkyy + (qkx*xr+qkz*zr)*dsc7*rr7
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = -dsr3 + term1*zr*zr - rr3*zr*drc3(3)
                  term4 = rr3*drc3(3) - term1*zr - dsr5*zr
                  term5 = term2*zr*zr - dsr5 - rr5*zr*drc5(3)
                  term6 = (dmpe(9)-dsc7*rr9)*zr*zr - dmpe(7)
     &                       - rr7*zr*drc7(3)
                  term7 = rr5*drc5(3) - 2.0d0*dmpe(7)*zr
     &                       + (dsc5+1.5d0*dsc7)*rr7*zr
                  tizz = ci*term3 + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qizz + (qix*xr+qiy*yr)*dsc7*rr7
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkzz = ck*term3 - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkzz + (qkx*xr+qky*yr)*dsc7*rr7
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*xr*yr - rr3*yr*drc3(1)
                  term4 = rr3*drc3(1) - term1*xr
                  term5 = term2*xr*yr - rr5*yr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*yr - rr7*yr*drc7(1)
                  term7 = rr5*drc5(1) - term2*xr
                  tixy = ci*term3 - dsr5*dix*yr + diy*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixy - 2.0d0*dsr7*yr*qix
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkxy = ck*term3 + dsr5*dkx*yr - dky*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxy - 2.0d0*dsr7*yr*qkx
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = term1*xr*zr - rr3*zr*drc3(1)
                  term5 = term2*xr*zr - rr5*zr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*zr - rr7*zr*drc7(1)
                  tixz = ci*term3 - dsr5*dix*zr + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixz - 2.0d0*dsr7*zr*qix
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkxz = ck*term3 + dsr5*dkx*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxz - 2.0d0*dsr7*zr*qkx
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*yr*zr - rr3*zr*drc3(2)
                  term4 = rr3*drc3(2) - term1*yr
                  term5 = term2*yr*zr - rr5*zr*drc5(2)
                  term6 = (dmpe(9)-dsc7*rr9)*yr*zr - rr7*zr*drc7(2)
                  term7 = rr5*drc5(2) - term2*yr
                  tiyz = ci*term3 - dsr5*diy*zr + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qiyz - 2.0d0*dsr7*zr*qiy
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkyz = ck*term3 + dsr5*dky*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkyz - 2.0d0*dsr7*zr*qky
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      - tkxx*uixp - tkxy*uiyp - tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      - tkxy*uixp - tkyy*uiyp - tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      - tkxz*uixp - tkyz*uiyp - tkzz*uizp
                  frcx = depx
                  frcy = depy
                  frcz = depz
c
c     get the dEp/dR terms used for direct polarization force
c
                  term1 = dmpe(5) - psc3*rr5
                  term2 = dmpe(7) - psc5*rr7
                  term3 = -psr3 + term1*xr*xr - rr3*xr*prc3(1)
                  term4 = rr3*prc3(1) - term1*xr - psr5*xr
                  term5 = term2*xr*xr - psr5 - rr5*xr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*xr - dmpe(7)
     &                       - rr7*xr*prc7(1)
                  term7 = rr5*prc5(1) - 2.0d0*dmpe(7)*xr
     &                       + (psc5+1.5d0*psc7)*rr7*xr
                  tixx = ci*term3 + dix*term4 + dir*term5
     &                      + 2.0d0*psr5*qixx + (qiy*yr+qiz*zr)*psc7*rr7
     &                      + 2.0d0*qix*term7 + qir*term6
                  tkxx = ck*term3 - dkx*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxx + (qky*yr+qkz*zr)*psc7*rr7
     &                      + 2.0d0*qkx*term7 + qkr*term6
                  term3 = -psr3 + term1*yr*yr - rr3*yr*prc3(2)
                  term4 = rr3*prc3(2) - term1*yr - psr5*yr
                  term5 = term2*yr*yr - psr5 - rr5*yr*prc5(2)
                  term6 = (dmpe(9)-psc7*rr9)*yr*yr - dmpe(7)
     &                       - rr7*yr*prc7(2)
                  term7 = rr5*prc5(2) - 2.0d0*dmpe(7)*yr
     &                       + (psc5+1.5d0*psc7)*rr7*yr
                  tiyy = ci*term3 + diy*term4 + dir*term5
     &                      + 2.0d0*psr5*qiyy + (qix*xr+qiz*zr)*psc7*rr7
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkyy = ck*term3 - dky*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkyy + (qkx*xr+qkz*zr)*psc7*rr7
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = -psr3 + term1*zr*zr - rr3*zr*prc3(3)
                  term4 = rr3*prc3(3) - term1*zr - psr5*zr
                  term5 = term2*zr*zr - psr5 - rr5*zr*prc5(3)
                  term6 = (dmpe(9)-psc7*rr9)*zr*zr - dmpe(7)
     &                       - rr7*zr*prc7(3)
                  term7 = rr5*prc5(3) - 2.0d0*dmpe(7)*zr
     &                       + (psc5+1.5d0*psc7)*rr7*zr
                  tizz = ci*term3 + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qizz + (qix*xr+qiy*yr)*psc7*rr7
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkzz = ck*term3 - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkzz + (qkx*xr+qky*yr)*psc7*rr7
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*xr*yr - rr3*yr*prc3(1)
                  term4 = rr3*prc3(1) - term1*xr
                  term5 = term2*xr*yr - rr5*yr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*yr - rr7*yr*prc7(1)
                  term7 = rr5*prc5(1) - term2*xr
                  tixy = ci*term3 - psr5*dix*yr + diy*term4 + dir*term5
     &                      + 2.0d0*psr5*qixy - 2.0d0*psr7*yr*qix
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkxy = ck*term3 + psr5*dkx*yr - dky*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxy - 2.0d0*psr7*yr*qkx
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = term1*xr*zr - rr3*zr*prc3(1)
                  term5 = term2*xr*zr - rr5*zr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*zr - rr7*zr*prc7(1)
                  tixz = ci*term3 - psr5*dix*zr + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qixz - 2.0d0*psr7*zr*qix
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkxz = ck*term3 + psr5*dkx*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxz - 2.0d0*psr7*zr*qkx
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*yr*zr - rr3*zr*prc3(2)
                  term4 = rr3*prc3(2) - term1*yr
                  term5 = term2*yr*zr - rr5*zr*prc5(2)
                  term6 = (dmpe(9)-psc7*rr9)*yr*zr - rr7*zr*prc7(2)
                  term7 = rr5*prc5(2) - term2*yr
                  tiyz = ci*term3 - psr5*diy*zr + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qiyz - 2.0d0*psr7*zr*qiy
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkyz = ck*term3 + psr5*dky*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkyz - 2.0d0*psr7*zr*qky
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
c
c     get the field gradient for direct polarization force
c
               else if (use_chgpen) then
                  term1i = rr3i - rr5i*xr*xr
                  term1core = rr3core - rr5core*xr*xr
                  term2i = 2.0d0*rr5i*xr 
                  term3i = rr7i*xr*xr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*xr
                  term6i = rr9i*xr*xr
                  term1k = rr3k - rr5k*xr*xr
                  term2k = 2.0d0*rr5k*xr
                  term3k = rr7k*xr*xr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*xr
                  term6k = rr9k*xr*xr
                  tixx = vali*term1i + corei*term1core  
     &                      + dix*term2i - dir*term3i
     &                      - qixx*term4i + qix*term5i - qir*term6i
     &                      + (qiy*yr+qiz*zr)*rr7i
                  tkxx = valk*term1k + corek*term1core
     &                      - dkx*term2k + dkr*term3k
     &                      - qkxx*term4k + qkx*term5k - qkr*term6k
     &                      + (qky*yr+qkz*zr)*rr7k
                  term1i = rr3i - rr5i*yr*yr
                  term1core = rr3core - rr5core*yr*yr
                  term2i = 2.0d0*rr5i*yr
                  term3i = rr7i*yr*yr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*yr
                  term6i = rr9i*yr*yr
                  term1k = rr3k - rr5k*yr*yr
                  term2k = 2.0d0*rr5k*yr
                  term3k = rr7k*yr*yr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*yr
                  term6k = rr9k*yr*yr
                  tiyy = vali*term1i + corei*term1core
     &                      + diy*term2i - dir*term3i
     &                      - qiyy*term4i + qiy*term5i - qir*term6i
     &                      + (qix*xr+qiz*zr)*rr7i
                  tkyy = valk*term1k + corek*term1core
     &                      - dky*term2k + dkr*term3k
     &                      - qkyy*term4k + qky*term5k - qkr*term6k
     &                      + (qkx*xr+qkz*zr)*rr7k
                  term1i = rr3i - rr5i*zr*zr
                  term1core = rr3core - rr5core*zr*zr
                  term2i = 2.0d0*rr5i*zr
                  term3i = rr7i*zr*zr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*zr
                  term6i = rr9i*zr*zr
                  term1k = rr3k - rr5k*zr*zr
                  term2k = 2.0d0*rr5k*zr
                  term3k = rr7k*zr*zr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*zr
                  term6k = rr9k*zr*zr
                  tizz = vali*term1i + corei*term1core
     &                      + diz*term2i - dir*term3i
     &                      - qizz*term4i + qiz*term5i - qir*term6i
     &                      + (qix*xr+qiy*yr)*rr7i
                  tkzz = valk*term1k + corek*term1core
     &                      - dkz*term2k + dkr*term3k
     &                      - qkzz*term4k + qkz*term5k - qkr*term6k
     &                      + (qkx*xr+qky*yr)*rr7k
                  term2i = rr5i*xr 
                  term1i = yr * term2i
                  term1core = rr5core*xr*yr
                  term3i = rr5i*yr
                  term4i = yr * (rr7i*xr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*xr
                  term7i = 2.0d0*rr7i*yr
                  term8i = yr*rr9i*xr
                  term2k = rr5k*xr
                  term1k = yr * term2k
                  term3k = rr5k*yr
                  term4k = yr * (rr7k*xr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*xr
                  term7k = 2.0d0*rr7k*yr
                  term8k = yr*rr9k*xr
                  tixy = -vali*term1i - corei*term1core 
     &                      + diy*term2i + dix*term3i
     &                      - dir*term4i - qixy*term5i + qiy*term6i
     &                      + qix*term7i - qir*term8i
                  tkxy = -valk*term1k - corek*term1core 
     &                      - dky*term2k - dkx*term3k
     &                      + dkr*term4k - qkxy*term5k + qky*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5i*xr
                  term1i = zr * term2i
                  term1core = rr5core*xr*zr
                  term3i = rr5i*zr
                  term4i = zr * (rr7i*xr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*xr
                  term7i = 2.0d0*rr7i*zr
                  term8i = zr*rr9i*xr
                  term2k = rr5k*xr
                  term1k = zr * term2k
                  term3k = rr5k*zr
                  term4k = zr * (rr7k*xr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*xr
                  term7k = 2.0d0*rr7k*zr
                  term8k = zr*rr9k*xr
                  tixz = -vali*term1i - corei*term1core
     &                      + diz*term2i + dix*term3i
     &                      - dir*term4i - qixz*term5i + qiz*term6i
     &                      + qix*term7i - qir*term8i
                  tkxz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dkx*term3k
     &                      + dkr*term4k - qkxz*term5k + qkz*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5i*yr
                  term1i = zr * term2i
                  term1core = rr5core*yr*zr
                  term3i = rr5i*zr
                  term4i = zr * (rr7i*yr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*yr
                  term7i = 2.0d0*rr7i*zr
                  term8i = zr*rr9i*yr
                  term2k = rr5k*yr
                  term1k = zr * term2k
                  term3k = rr5k*zr
                  term4k = zr * (rr7k*yr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*yr
                  term7k = 2.0d0*rr7k*zr
                  term8k = zr*rr9k*yr
                  tiyz = -vali*term1i - corei*term1core
     &                      + diz*term2i + diy*term3i
     &                      - dir*term4i - qiyz*term5i + qiz*term6i
     &                      + qiy*term7i - qir*term8i
                  tkyz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dky*term3k
     &                      + dkr*term4k - qkyz*term5k + qkz*term6k
     &                      + qky*term7k - qkr*term8k
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = -2.0d0 * depx
                  frcy = -2.0d0 * depy
                  frcz = -2.0d0 * depz
               end if
c
c     reset Thole values if alternate direct damping was used
c
               if (use_dirdamp) then
                  sc3 = 1.0d0
                  sc5 = 1.0d0
                  do j = 1, 3
                     rc3(j) = 0.0d0
                     rc5(j) = 0.0d0
                  end do
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        sc3 = 1.0d0 - expdamp
                        sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                        temp3 = 3.0d0 * damp * expdamp / r2
                        temp5 = damp
                        rc3(1) = xr * temp3
                        rc3(2) = yr * temp3
                        rc3(3) = zr * temp3
                        rc5(1) = rc3(1) * temp5
                        rc5(2) = rc3(2) * temp5
                        rc5(3) = rc3(3) * temp5
                     end if
                  end if
                  usc3 = 1.0d0 - sc3*uscale(k)
                  usc5 = 1.0d0 - sc5*uscale(k)
                  usr3 = dmpe(3) - usc3*rr3
                  usr5 = dmpe(5) - usc5*rr5
                  do j = 1, 3
                     urc3(j) = rc3(j) * uscale(k)
                     urc5(j) = rc5(j) * uscale(k)
                  end do
               end if
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp.eq.'MUTUAL' .and. use_thole) then
                  term1 = dmpe(5) - usc3*rr5
                  term2 = dmpe(7) - usc5*rr7
                  term3 = usr5 + term1
                  term4 = rr3 * uscale(k)
                  term5 = -xr*term3 + rc3(1)*term4
                  term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                  tixx = uix*term5 + uir*term6
                  tkxx = ukx*term5 + ukr*term6
                  term5 = -yr*term3 + rc3(2)*term4
                  term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                  tiyy = uiy*term5 + uir*term6
                  tkyy = uky*term5 + ukr*term6
                  term5 = -zr*term3 + rc3(3)*term4
                  term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                  tizz = uiz*term5 + uir*term6
                  tkzz = ukz*term5 + ukr*term6
                  term4 = -usr5 * yr
                  term5 = -xr*term1 + rr3*urc3(1)
                  term6 = xr*yr*term2 - rr5*yr*urc5(1)
                  tixy = uix*term4 + uiy*term5 + uir*term6
                  tkxy = ukx*term4 + uky*term5 + ukr*term6
                  term4 = -usr5 * zr
                  term6 = xr*zr*term2 - rr5*zr*urc5(1)
                  tixz = uix*term4 + uiz*term5 + uir*term6
                  tkxz = ukx*term4 + ukz*term5 + ukr*term6
                  term5 = -yr*term1 + rr3*urc3(2)
                  term6 = yr*zr*term2 - rr5*zr*urc5(2)
                  tiyz = uiy*term4 + uiz*term5 + uir*term6
                  tkyz = uky*term4 + ukz*term5 + ukr*term6
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               else if (poltyp.eq.'MUTUAL' .and. use_chgpen) then
                  term1 = 2.0d0 * rr5ik
                  term2 = term1*xr
                  term3 = rr5ik - rr7ik*xr*xr
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr
                  term3 = rr5ik - rr7ik*yr*yr
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr
                  term3 = rr5ik - rr7ik*zr*zr
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = rr5ik*yr
                  term2 = rr5ik*xr
                  term3 = yr * (rr7ik*xr)
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = rr5ik * zr
                  term3 = zr * (rr7ik*xr)
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = rr5ik*yr
                  term3 = zr * (rr7ik*yr)
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx - depx
                  frcy = frcy - depy
                  frcz = frcz - depz
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_thole) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,ii)*xr + uopt(j,2,ii)*yr
     &                          + uopt(j,3,ii)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,kk)*xr + uopt(m,2,kk)*yr
     &                             + uopt(m,3,kk)*zr
                        term1 = dmpe(5) - usc3*rr5
                        term2 = dmpe(7) - usc5*rr7
                        term3 = usr5 + term1
                        term4 = rr3 * uscale(k)
                        term5 = -xr*term3 + rc3(1)*term4
                        term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                        tixx = uopt(j,1,ii)*term5 + uirm*term6
                        tkxx = uopt(m,1,kk)*term5 + ukrm*term6
                        term5 = -yr*term3 + rc3(2)*term4
                        term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                        tiyy = uopt(j,2,ii)*term5 + uirm*term6
                        tkyy = uopt(m,2,kk)*term5 + ukrm*term6
                        term5 = -zr*term3 + rc3(3)*term4
                        term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                        tizz = uopt(j,3,ii)*term5 + uirm*term6
                        tkzz = uopt(m,3,kk)*term5 + ukrm*term6
                        term4 = -usr5 * yr
                        term5 = -xr*term1 + rr3*urc3(1)
                        term6 = xr*yr*term2 - rr5*yr*urc5(1)
                        tixy = uopt(j,1,ii)*term4 + uopt(j,2,ii)*term5
     &                            + uirm*term6
                        tkxy = uopt(m,1,kk)*term4 + uopt(m,2,kk)*term5
     &                            + ukrm*term6
                        term4 = -usr5 * zr
                        term6 = xr*zr*term2 - rr5*zr*urc5(1)
                        tixz = uopt(j,1,ii)*term4 + uopt(j,3,ii)*term5
     &                            + uirm*term6
                        tkxz = uopt(m,1,kk)*term4 + uopt(m,3,kk)*term5
     &                            + ukrm*term6
                        term5 = -yr*term1 + rr3*urc3(2)
                        term6 = yr*zr*term2 - rr5*zr*urc5(2)
                        tiyz = uopt(j,2,ii)*term4 + uopt(j,3,ii)*term5
     &                            + uirm*term6
                        tkyz = uopt(m,2,kk)*term4 + uopt(m,3,kk)*term5
     &                            + ukrm*term6
                        depx = tixx*uoptp(m,1,kk) + tkxx*uoptp(j,1,ii)
     &                       + tixy*uoptp(m,2,kk) + tkxy*uoptp(j,2,ii)
     &                       + tixz*uoptp(m,3,kk) + tkxz*uoptp(j,3,ii)
                        depy = tixy*uoptp(m,1,kk) + tkxy*uoptp(j,1,ii)
     &                       + tiyy*uoptp(m,2,kk) + tkyy*uoptp(j,2,ii)
     &                       + tiyz*uoptp(m,3,kk) + tkyz*uoptp(j,3,ii)
                        depz = tixz*uoptp(m,1,kk) + tkxz*uoptp(j,1,ii)
     &                       + tiyz*uoptp(m,2,kk) + tkyz*uoptp(j,2,ii)
     &                       + tizz*uoptp(m,3,kk) + tkzz*uoptp(j,3,ii)
                        frcx = frcx + copm(j+m+1)*depx
                        frcy = frcy + copm(j+m+1)*depy
                        frcz = frcz + copm(j+m+1)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_chgpen) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,i)*xr + uopt(j,2,i)*yr
     &                          + uopt(j,3,i)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,k)*xr + uopt(m,2,k)*yr
     &                             + uopt(m,3,k)*zr
                        term1 = 2.0d0 * rr5ik
                        term2 = term1*xr
                        term3 = rr5ik - rr7ik*xr*xr
                        tixx = uopt(j,1,i)*term2 + uirm*term3
                        tkxx = uopt(m,1,k)*term2 + ukrm*term3
                        term2 = term1*yr
                        term3 = rr5ik - rr7ik*yr*yr
                        tiyy = uopt(j,2,i)*term2 + uirm*term3
                        tkyy = uopt(m,2,k)*term2 + ukrm*term3
                        term2 = term1*zr
                        term3 = rr5ik - rr7ik*zr*zr
                        tizz = uopt(j,3,i)*term2 + uirm*term3
                        tkzz = uopt(m,3,k)*term2 + ukrm*term3
                        term1 = rr5ik*yr
                        term2 = rr5ik*xr
                        term3 = yr * (rr7ik*xr)
                        tixy = uopt(j,1,i)*term1 + uopt(j,2,i)*term2 
     &                       - uirm*term3
                        tkxy = uopt(m,1,k)*term1 + uopt(m,2,k)*term2 
     &                       - ukrm*term3
                        term1 = rr5ik * zr
                        term3 = zr * (rr7ik*xr)
                        tixz = uopt(j,1,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        term2 = rr5ik*yr
                        term3 = zr * (rr7ik*yr)
                        tiyz = uopt(j,2,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,k) + tkxx*uoptp(j,1,i)
     &                       + tixy*uoptp(m,2,k) + tkxy*uoptp(j,2,i)
     &                       + tixz*uoptp(m,3,k) + tkxz*uoptp(j,3,i)
                        depy = tixy*uoptp(m,1,k) + tkxy*uoptp(j,1,i)
     &                       + tiyy*uoptp(m,2,k) + tkyy*uoptp(j,2,i)
     &                       + tiyz*uoptp(m,3,k) + tkyz*uoptp(j,3,i)
                        depz = tixz*uoptp(m,1,k) + tkxz*uoptp(j,1,i)
     &                       + tiyz*uoptp(m,2,k) + tkyz*uoptp(j,2,i)
     &                       + tizz*uoptp(m,3,k) + tkzz*uoptp(j,3,i)
                        frcx = frcx - copm(j+m+1)*depx
                        frcy = frcy - copm(j+m+1)*depy
                        frcz = frcz - copm(j+m+1)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for TCG polarization force
c
               else if (poltyp.eq.'TCG' .and. use_thole) then
                  do j = 1, tcgnab
                     ukx = ubd(1,kk,j)
                     uky = ubd(2,kk,j)
                     ukz = ubd(3,kk,j)
                     ukxp = ubp(1,kk,j)
                     ukyp = ubp(2,kk,j)
                     ukzp = ubp(3,kk,j)
                     uirt = uax(j)*xr + uay(j)*yr + uaz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = dmpe(5) - usc3*rr5
                     term2 = dmpe(7) - usc5*rr7
                     term3 = usr5 + term1
                     term4 = rr3 * uscale(k)
                     term5 = -xr*term3 + rc3(1)*term4
                     term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                     tixx = uax(j)*term5 + uirt*term6
                     tkxx = ukx*term5 + ukrt*term6
                     term5 = -yr*term3 + rc3(2)*term4
                     term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                     tiyy = uay(j)*term5 + uirt*term6
                     tkyy = uky*term5 + ukrt*term6
                     term5 = -zr*term3 + rc3(3)*term4
                     term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                     tizz = uaz(j)*term5 + uirt*term6
                     tkzz = ukz*term5 + ukrt*term6
                     term4 = -usr5 * yr
                     term5 = -xr*term1 + rr3*urc3(1)
                     term6 = xr*yr*term2 - rr5*yr*urc5(1)
                     tixy = uax(j)*term4 + uay(j)*term5 + uirt*term6
                     tkxy = ukx*term4 + uky*term5 + ukrt*term6
                     term4 = -usr5 * zr
                     term6 = xr*zr*term2 - rr5*zr*urc5(1)
                     tixz = uax(j)*term4 + uaz(j)*term5 + uirt*term6
                     tkxz = ukx*term4 + ukz*term5 + ukrt*term6
                     term5 = -yr*term1 + rr3*urc3(2)
                     term6 = yr*zr*term2 - rr5*zr*urc5(2)
                     tiyz = uay(j)*term4 + uaz(j)*term5 + uirt*term6
                     tkyz = uky*term4 + ukz*term5 + ukrt*term6
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*uaxp(j) + tkxy*uayp(j)
     &                         + tkxz*uazp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*uaxp(j) + tkyy*uayp(j)
     &                         + tkyz*uazp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*uaxp(j) + tkyz*uayp(j)
     &                         + tkzz*uazp(j)
                     frcx = frcx + depx
                     frcy = frcy + depy
                     frcz = frcz + depz
                     ukx = uad(1,kk,j)
                     uky = uad(2,kk,j)
                     ukz = uad(3,kk,j)
                     ukxp = uap(1,kk,j)
                     ukyp = uap(2,kk,j)
                     ukzp = uap(3,kk,j)
                     uirt = ubx(j)*xr + uby(j)*yr + ubz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = dmpe(5) - usc3*rr5
                     term2 = dmpe(7) - usc5*rr7
                     term3 = usr5 + term1
                     term4 = rr3 * uscale(k)
                     term5 = -xr*term3 + rc3(1)*term4
                     term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                     tixx = ubx(j)*term5 + uirt*term6
                     tkxx = ukx*term5 + ukrt*term6
                     term5 = -yr*term3 + rc3(2)*term4
                     term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                     tiyy = uby(j)*term5 + uirt*term6
                     tkyy = uky*term5 + ukrt*term6
                     term5 = -zr*term3 + rc3(3)*term4
                     term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                     tizz = ubz(j)*term5 + uirt*term6
                     tkzz = ukz*term5 + ukrt*term6
                     term4 = -usr5 * yr
                     term5 = -xr*term1 + rr3*urc3(1)
                     term6 = xr*yr*term2 - rr5*yr*urc5(1)
                     tixy = ubx(j)*term4 + uby(j)*term5 + uirt*term6
                     tkxy = ukx*term4 + uky*term5 + ukrt*term6
                     term4 = -usr5 * zr
                     term6 = xr*zr*term2 - rr5*zr*urc5(1)
                     tixz = ubx(j)*term4 + ubz(j)*term5 + uirt*term6
                     tkxz = ukx*term4 + ukz*term5 + ukrt*term6
                     term5 = -yr*term1 + rr3*urc3(2)
                     term6 = yr*zr*term2 - rr5*zr*urc5(2)
                     tiyz = uby(j)*term4 + ubz(j)*term5 + uirt*term6
                     tkyz = uky*term4 + ukz*term5 + ukrt*term6
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*ubxp(j) + tkxy*ubyp(j)
     &                         + tkxz*ubzp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*ubxp(j) + tkyy*ubyp(j)
     &                         + tkyz*ubzp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*ubxp(j) + tkyz*ubyp(j)
     &                         + tkzz*ubzp(j)
                     frcx = frcx + depx
                     frcy = frcy + depy
                     frcz = frcz + depz
                  end do
               end if
c
c     force and torque components scaled for self-interactions
c
               if (i .eq. k) then
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
c     increment force-based gradient on the interaction sites
c
               dep(1,i) = dep(1,i) - frcx
               dep(2,i) = dep(2,i) - frcy
               dep(3,i) = dep(3,i) - frcz
               dep(1,k) = dep(1,k) + frcx
               dep(2,k) = dep(2,k) + frcy
               dep(3,k) = dep(3,k) + frcz
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = xr * frcx
               vxy = 0.5d0 * (yr*frcx+xr*frcy)
               vxz = 0.5d0 * (zr*frcx+xr*frcz)
               vyy = yr * frcy
               vyz = 0.5d0 * (zr*frcy+yr*frcz)
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
            end if
            end do
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = 1.0d0
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
               uscale(ip14(j,i)) = 1.0d0
            end do
         end if
      end do
      end if
c
c     torque is induced field and gradient cross permanent moments
c
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         tep(1) = diz*ufld(2,i) - diy*ufld(3,i)
     &               + qixz*dufld(2,i) - qixy*dufld(4,i)
     &               + 2.0d0*qiyz*(dufld(3,i)-dufld(6,i))
     &               + (qizz-qiyy)*dufld(5,i)
         tep(2) = dix*ufld(3,i) - diz*ufld(1,i)
     &               - qiyz*dufld(2,i) + qixy*dufld(5,i)
     &               + 2.0d0*qixz*(dufld(6,i)-dufld(1,i))
     &               + (qixx-qizz)*dufld(4,i)
         tep(3) = diy*ufld(1,i) - dix*ufld(2,i)
     &               + qiyz*dufld(4,i) - qixz*dufld(5,i)
     &               + 2.0d0*qixy*(dufld(1,i)-dufld(3,i))
     &               + (qiyy-qixx)*dufld(2,i)
         call torque (ii,tep,fix,fiy,fiz,dep)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = abs(yaxis(ii))
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
            dep(1,i) = dep(1,i) + frcx
            dep(2,i) = dep(2,i) + frcy
            dep(3,i) = dep(3,i) + frcz
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
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (wscale)
      deallocate (ufld)
      deallocate (dufld)
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
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
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use poltcg
      use potent
      use virial
      implicit none
      integer i,j,ii
      integer ix,iy,iz
      real*8 f,term
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 xd,yd,zd
      real*8 xq,yq,zq
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xv,yv,zv,vterm
      real*8 xufield,yufield
      real*8 zufield
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 fix(3),fiy(3),fiz(3)
      real*8 tep(3)
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
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bsporder
      aewald = apewald
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
c     compute the total induced dipole polarization energy
c
      call epolar1e
c
c     compute the real space part of the Ewald summation
c
      call epreal1d
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip1
c
c     compute the Ewald self-energy torque and virial terms
c
      term = (4.0d0/3.0d0) * f * aewald**3 / rootpi
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         uix = 0.5d0 * (uind(1,ii)+uinp(1,ii))
         uiy = 0.5d0 * (uind(2,ii)+uinp(2,ii))
         uiz = 0.5d0 * (uind(3,ii)+uinp(3,ii))
         tep(1) = term * (diy*uiz-diz*uiy)
         tep(2) = term * (diz*uix-dix*uiz)
         tep(3) = term * (dix*uiy-diy*uix)
         call torque (ii,tep,fix,fiy,fiz,dep)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = abs(yaxis(ii))
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
     &                     + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                     + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3))
         vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                     + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
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
         do ii = 1, npole
            i = ipole(ii)
            xd = xd + rpole(2,ii) + rpole(1,ii)*x(ii)
            yd = yd + rpole(3,ii) + rpole(1,ii)*y(ii)
            zd = zd + rpole(4,ii) + rpole(1,ii)*z(ii)
            xu = xu + uind(1,ii)
            yu = yu + uind(2,ii)
            zu = zu + uind(3,ii)
            xup = xup + uinp(1,ii)
            yup = yup + uinp(2,ii)
            zup = zup + uinp(3,ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         do ii = 1, npole
            i = ipole(ii)
            dep(1,i) = dep(1,i) + term*rpole(1,ii)*(xu+xup)
            dep(2,i) = dep(2,i) + term*rpole(1,ii)*(yu+yup)
            dep(3,i) = dep(3,i) + term*rpole(1,ii)*(zu+zup)
         end do
         xufield = -term * (xu+xup)
         yufield = -term * (yu+yup)
         zufield = -term * (zu+zup)
         do ii = 1, npole
            tep(1) = rpole(3,ii)*zufield - rpole(4,ii)*yufield
            tep(2) = rpole(4,ii)*xufield - rpole(2,ii)*zufield
            tep(3) = rpole(2,ii)*yufield - rpole(3,ii)*xufield
            call torque (ii,tep,fix,fiy,fiz,dep)
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
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use deriv
      use ewald
      use math
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polopt
      use polpot
      use poltcg
      use potent
      use shunt
      use virial
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      integer ix,iy,iz
      real*8 f,pgamma
      real*8 pdi,pti,ddi
      real*8 damp,expdamp
      real*8 temp3,temp5,temp7
      real*8 sc3,sc5,sc7
      real*8 psc3,psc5,psc7
      real*8 dsc3,dsc5,dsc7
      real*8 usc3,usc5
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 usr3,usr5
      real*8 rr3core,rr5core
      real*8 rr3i,rr5i
      real*8 rr7i,rr9i
      real*8 rr3k,rr5k
      real*8 rr7k,rr9k
      real*8 rr5ik,rr7ik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 uixp,uiyp,uizp
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 ukxp,ukyp,ukzp
      real*8 dir,uir,uirp
      real*8 dkr,ukr,ukrp
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 uirm,ukrm
      real*8 uirt,ukrt
      real*8 tuir,tukr
      real*8 tixx,tiyy,tizz
      real*8 tixy,tixz,tiyz
      real*8 tkxx,tkyy,tkzz
      real*8 tkxy,tkxz,tkyz
      real*8 tix3,tiy3,tiz3
      real*8 tix5,tiy5,tiz5
      real*8 tkx3,tky3,tkz3
      real*8 tkx5,tky5,tkz5
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term6,term7
      real*8 term1core
      real*8 term1i,term2i,term3i
      real*8 term4i,term5i,term6i
      real*8 term7i,term8i
      real*8 term1k,term2k,term3k
      real*8 term4k,term5k,term6k
      real*8 term7k,term8k
      real*8 poti,potk
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
      real*8 urc3(3),urc5(3),tep(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 uax(3),uay(3),uaz(3)
      real*8 ubx(3),uby(3),ubz(3)
      real*8 uaxp(3),uayp(3),uazp(3)
      real*8 ubxp(3),ubyp(3),ubzp(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9),dmpe(9)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      real*8, allocatable :: pot(:)
      real*8, allocatable :: decfx(:)
      real*8, allocatable :: decfy(:)
      real*8, allocatable :: decfz(:)
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (wscale(n))
      allocate (ufld(3,n))
      allocate (dufld(6,n))
      allocate (pot(n))
      allocate (decfx(n))
      allocate (decfy(n))
      allocate (decfz(n))
c
c     set exclusion coefficients and arrays to store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
         wscale(i) = 1.0d0
         do j = 1, 3
            ufld(j,i) = 0.0d0
         end do
         do j = 1, 6
            dufld(j,i) = 0.0d0
         end do
         pot(i) = 0.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,x,y,z,rpole,uind,
!$OMP& uinp,pdamp,thole,dirdamp,pcore,pval,palpha,n12,i12,n13,i13,n14,
!$OMP& i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,p2scale,
!$OMP& p3scale,p4scale,p5scale,p2iscale,p3iscale,p4iscale,p5iscale,
!$OMP& d1scale,d2scale,d3scale,d4scale,u1scale,u2scale,u3scale,u4scale,
!$OMP& w2scale,w3scale,w4scale,w5scale,nelst,elst,dpequal,use_thole,
!$OMP& use_chgpen,use_chgflx,use_dirdamp,use_bounds,off2,f,aewald,
!$OMP& optorder,copm,uopt,uoptp,poltyp,tcgnab,uad,uap,ubd,ubp,xaxis,
!$OMP& yaxis,zaxis)
!$OMP& shared (dep,ufld,dufld,pot,vir)
!$OMP& firstprivate(pscale,dscale,uscale,wscale)
!$OMP DO reduction(+:dep,ufld,dufld,pot,vir) schedule(guided)
c
c     compute the dipole polarization gradient components
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
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         uixp = uinp(1,ii)
         uiyp = uinp(2,ii)
         uizp = uinp(3,ii)
         do j = 1, tcgnab
            uax(j) = uad(1,ii,j)
            uay(j) = uad(2,ii,j)
            uaz(j) = uad(3,ii,j)
            uaxp(j) = uap(1,ii,j)
            uayp(j) = uap(2,ii,j)
            uazp(j) = uap(3,ii,j)
            ubx(j) = ubd(1,ii,j)
            uby(j) = ubd(2,ii,j)
            ubz(j) = ubd(3,ii,j)
            ubxp(j) = ubp(1,ii,j)
            ubyp(j) = ubp(2,ii,j)
            ubzp(j) = ubp(3,ii,j)
         end do
         if (use_thole) then
            pdi = pdamp(ii)
            pti = thole(ii)
            ddi = dirdamp(ii)
         else if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               dscale(i12(j,i)) = pscale(i12(j,i))
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = u4scale
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
               uscale(ip14(j,i)) = u4scale
            end do
         end if
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
               ukx = uind(1,kk)
               uky = uind(2,kk)
               ukz = uind(3,kk)
               ukxp = uinp(1,kk)
               ukyp = uinp(2,kk)
               ukzp = uinp(3,kk)
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
               uir = uix*xr + uiy*yr + uiz*zr
               uirp = uixp*xr + uiyp*yr + uizp*zr
               ukr = ukx*xr + uky*yr + ukz*zr
               ukrp = ukxp*xr + ukyp*yr + ukzp*zr
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate real space Ewald error function damping
c
               call dampewald (9,r,r2,f,dmpe)
c
c     set initial values for tha damping scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
c
c     apply Thole polarization damping to scale factors
c
               if (use_thole) then
                  damp = pdi * pdamp(kk)
                  if (use_dirdamp) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(ddi,dirdamp(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           sc3 = 1.0d0 - expdamp 
                           sc5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                      +0.15d0*damp**2)
                           temp3 = 1.5d0 * damp * expdamp / r2
                           temp5 = 0.5d0 * (1.0d0+damp)
                           temp7 = 0.7d0 + 0.15d0*damp**2/temp5
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
                  else
                     pgamma = min(pti,thole(kk))
                     if (pgamma .eq. 0.0d0) then
                        pgamma = max(pti,thole(kk))
                     end if
                     if (damp.ne.0.0d0 .and. pgamma.ne.0.0d0) then
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           sc3 = 1.0d0 - expdamp
                           sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                           sc7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                      +0.6d0*damp**2)
                           temp3 = 3.0d0 * damp * expdamp / r2
                           temp5 = damp
                           temp7 = -0.2d0 + 0.6d0*damp
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
                  end if
                  psc3 = 1.0d0 - sc3*pscale(k)
                  psc5 = 1.0d0 - sc5*pscale(k)
                  psc7 = 1.0d0 - sc7*pscale(k)
                  dsc3 = 1.0d0 - sc3*dscale(k)
                  dsc5 = 1.0d0 - sc5*dscale(k)
                  dsc7 = 1.0d0 - sc7*dscale(k)
                  usc3 = 1.0d0 - sc3*uscale(k)
                  usc5 = 1.0d0 - sc5*uscale(k)
                  psr3 = dmpe(3) - psc3*rr3
                  psr5 = dmpe(5) - psc5*rr5
                  psr7 = dmpe(7) - psc7*rr7
                  dsr3 = dmpe(3) - dsc3*rr3
                  dsr5 = dmpe(5) - dsc5*rr5
                  dsr7 = dmpe(7) - dsc7*rr7
                  usr3 = dmpe(3) - usc3*rr3
                  usr5 = dmpe(5) - usc5*rr5
                  do j = 1, 3
                     prc3(j) = rc3(j) * pscale(k)
                     prc5(j) = rc5(j) * pscale(k)
                     prc7(j) = rc7(j) * pscale(k)
                     drc3(j) = rc3(j) * dscale(k)
                     drc5(j) = rc5(j) * dscale(k)
                     drc7(j) = rc7(j) * dscale(k)
                     urc3(j) = rc3(j) * uscale(k)
                     urc5(j) = rc5(j) * uscale(k)
                  end do
c
c     apply charge penetration damping to scale factors
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call damppole (r,9,alphai,alphak,dmpi,dmpk,dmpik)
                  rr3core = dmpe(3) - (1.0d0-dscale(k))*rr3
                  rr5core = dmpe(5) - (1.0d0-dscale(k))*rr5
                  rr3i = dmpe(3) - (1.0d0-dscale(k)*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-dscale(k)*dmpi(5))*rr5
                  rr7i = dmpe(7) - (1.0d0-dscale(k)*dmpi(7))*rr7
                  rr9i = dmpe(9) - (1.0d0-dscale(k)*dmpi(9))*rr9
                  rr3k = dmpe(3) - (1.0d0-dscale(k)*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-dscale(k)*dmpk(5))*rr5
                  rr7k = dmpe(7) - (1.0d0-dscale(k)*dmpk(7))*rr7
                  rr9k = dmpe(9) - (1.0d0-dscale(k)*dmpk(9))*rr9
                  rr5ik = dmpe(5) - (1.0d0-wscale(k)*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-wscale(k)*dmpik(7))*rr7
               end if
c
c     store the potential at each site for use in charge flux
c
               if (use_chgflx) then
                  if (use_thole) then
                     poti = -ukr*psr3 - ukrp*dsr3
                     potk = uir*psr3 + uirp*dsr3
                  else if (use_chgpen) then
                     poti = -2.0d0 * ukr * rr3i
                     potk = 2.0d0 * uir * rr3k
                  end if
                  pot(i) = pot(i) + poti 
                  pot(k) = pot(k) + potk 
               end if
c
c     get the induced dipole field used for dipole torques
c
               if (use_thole) then
                  tix3 = psr3*ukx + dsr3*ukxp
                  tiy3 = psr3*uky + dsr3*ukyp
                  tiz3 = psr3*ukz + dsr3*ukzp
                  tkx3 = psr3*uix + dsr3*uixp
                  tky3 = psr3*uiy + dsr3*uiyp
                  tkz3 = psr3*uiz + dsr3*uizp
                  tuir = -psr5*ukr - dsr5*ukrp
                  tukr = -psr5*uir - dsr5*uirp
               else if (use_chgpen) then
                  tix3 = 2.0d0*rr3i*ukx
                  tiy3 = 2.0d0*rr3i*uky
                  tiz3 = 2.0d0*rr3i*ukz
                  tkx3 = 2.0d0*rr3k*uix
                  tky3 = 2.0d0*rr3k*uiy
                  tkz3 = 2.0d0*rr3k*uiz
                  tuir = -2.0d0*rr5i*ukr
                  tukr = -2.0d0*rr5k*uir
               end if
               ufld(1,i) = ufld(1,i) + tix3 + xr*tuir
               ufld(2,i) = ufld(2,i) + tiy3 + yr*tuir
               ufld(3,i) = ufld(3,i) + tiz3 + zr*tuir
               ufld(1,k) = ufld(1,k) + tkx3 + xr*tukr
               ufld(2,k) = ufld(2,k) + tky3 + yr*tukr
               ufld(3,k) = ufld(3,k) + tkz3 + zr*tukr
c
c     get induced dipole field gradient used for quadrupole torques
c
               if (use_thole) then
                  tix5 = 2.0d0 * (psr5*ukx+dsr5*ukxp)
                  tiy5 = 2.0d0 * (psr5*uky+dsr5*ukyp)
                  tiz5 = 2.0d0 * (psr5*ukz+dsr5*ukzp)
                  tkx5 = 2.0d0 * (psr5*uix+dsr5*uixp)
                  tky5 = 2.0d0 * (psr5*uiy+dsr5*uiyp)
                  tkz5 = 2.0d0 * (psr5*uiz+dsr5*uizp)
                  tuir = -psr7*ukr - dsr7*ukrp
                  tukr = -psr7*uir - dsr7*uirp
               else if (use_chgpen) then
                  tix5 = 4.0d0 * (rr5i*ukx)
                  tiy5 = 4.0d0 * (rr5i*uky)
                  tiz5 = 4.0d0 * (rr5i*ukz)
                  tkx5 = 4.0d0 * (rr5k*uix)
                  tky5 = 4.0d0 * (rr5k*uiy)
                  tkz5 = 4.0d0 * (rr5k*uiz)
                  tuir = -2.0d0*rr7i*ukr 
                  tukr = -2.0d0*rr7k*uir 
               end if
               dufld(1,i) = dufld(1,i) + xr*tix5 + xr*xr*tuir
               dufld(2,i) = dufld(2,i) + xr*tiy5 + yr*tix5
     &                         + 2.0d0*xr*yr*tuir
               dufld(3,i) = dufld(3,i) + yr*tiy5 + yr*yr*tuir
               dufld(4,i) = dufld(4,i) + xr*tiz5 + zr*tix5
     &                         + 2.0d0*xr*zr*tuir
               dufld(5,i) = dufld(5,i) + yr*tiz5 + zr*tiy5
     &                         + 2.0d0*yr*zr*tuir
               dufld(6,i) = dufld(6,i) + zr*tiz5 + zr*zr*tuir
               dufld(1,k) = dufld(1,k) - xr*tkx5 - xr*xr*tukr
               dufld(2,k) = dufld(2,k) - xr*tky5 - yr*tkx5
     &                         - 2.0d0*xr*yr*tukr
               dufld(3,k) = dufld(3,k) - yr*tky5 - yr*yr*tukr
               dufld(4,k) = dufld(4,k) - xr*tkz5 - zr*tkx5
     &                         - 2.0d0*xr*zr*tukr
               dufld(5,k) = dufld(5,k) - yr*tkz5 - zr*tky5
     &                         - 2.0d0*yr*zr*tukr
               dufld(6,k) = dufld(6,k) - zr*tkz5 - zr*zr*tukr
c
c     get the dEd/dR terms used for direct polarization force
c
               if (use_thole) then
                  term1 = dmpe(5) - dsc3*rr5
                  term2 = dmpe(7) - dsc5*rr7
                  term3 = -dsr3 + term1*xr*xr - rr3*xr*drc3(1)
                  term4 = rr3*drc3(1) - term1*xr - dsr5*xr
                  term5 = term2*xr*xr - dsr5 - rr5*xr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*xr - dmpe(7)
     &                       - rr7*xr*drc7(1)
                  term7 = rr5*drc5(1) - 2.0d0*dmpe(7)*xr
     &                       + (dsc5+1.5d0*dsc7)*rr7*xr
                  tixx = ci*term3 + dix*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixx + (qiy*yr+qiz*zr)*dsc7*rr7
     &                      + 2.0d0*qix*term7 + qir*term6
                  tkxx = ck*term3 - dkx*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxx + (qky*yr+qkz*zr)*dsc7*rr7
     &                      + 2.0d0*qkx*term7 + qkr*term6
                  term3 = -dsr3 + term1*yr*yr - rr3*yr*drc3(2)
                  term4 = rr3*drc3(2) - term1*yr - dsr5*yr
                  term5 = term2*yr*yr - dsr5 - rr5*yr*drc5(2)
                  term6 = (dmpe(9)-dsc7*rr9)*yr*yr - dmpe(7)
     &                       - rr7*yr*drc7(2)
                  term7 = rr5*drc5(2) - 2.0d0*dmpe(7)*yr
     &                       + (dsc5+1.5d0*dsc7)*rr7*yr
                  tiyy = ci*term3 + diy*term4 + dir*term5
     &                      + 2.0d0*dsr5*qiyy + (qix*xr+qiz*zr)*dsc7*rr7
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkyy = ck*term3 - dky*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkyy + (qkx*xr+qkz*zr)*dsc7*rr7
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = -dsr3 + term1*zr*zr - rr3*zr*drc3(3)
                  term4 = rr3*drc3(3) - term1*zr - dsr5*zr
                  term5 = term2*zr*zr - dsr5 - rr5*zr*drc5(3)
                  term6 = (dmpe(9)-dsc7*rr9)*zr*zr - dmpe(7)
     &                       - rr7*zr*drc7(3)
                  term7 = rr5*drc5(3) - 2.0d0*dmpe(7)*zr
     &                       + (dsc5+1.5d0*dsc7)*rr7*zr
                  tizz = ci*term3 + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qizz + (qix*xr+qiy*yr)*dsc7*rr7
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkzz = ck*term3 - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkzz + (qkx*xr+qky*yr)*dsc7*rr7
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*xr*yr - rr3*yr*drc3(1)
                  term4 = rr3*drc3(1) - term1*xr
                  term5 = term2*xr*yr - rr5*yr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*yr - rr7*yr*drc7(1)
                  term7 = rr5*drc5(1) - term2*xr
                  tixy = ci*term3 - dsr5*dix*yr + diy*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixy - 2.0d0*dsr7*yr*qix
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkxy = ck*term3 + dsr5*dkx*yr - dky*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxy - 2.0d0*dsr7*yr*qkx
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = term1*xr*zr - rr3*zr*drc3(1)
                  term5 = term2*xr*zr - rr5*zr*drc5(1)
                  term6 = (dmpe(9)-dsc7*rr9)*xr*zr - rr7*zr*drc7(1)
                  tixz = ci*term3 - dsr5*dix*zr + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qixz - 2.0d0*dsr7*zr*qix
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkxz = ck*term3 + dsr5*dkx*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkxz - 2.0d0*dsr7*zr*qkx
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*yr*zr - rr3*zr*drc3(2)
                  term4 = rr3*drc3(2) - term1*yr
                  term5 = term2*yr*zr - rr5*zr*drc5(2)
                  term6 = (dmpe(9)-dsc7*rr9)*yr*zr - rr7*zr*drc7(2)
                  term7 = rr5*drc5(2) - term2*yr
                  tiyz = ci*term3 - dsr5*diy*zr + diz*term4 + dir*term5
     &                      + 2.0d0*dsr5*qiyz - 2.0d0*dsr7*zr*qiy
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkyz = ck*term3 + dsr5*dky*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*dsr5*qkyz - 2.0d0*dsr7*zr*qky
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      - tkxx*uixp - tkxy*uiyp - tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      - tkxy*uixp - tkyy*uiyp - tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      - tkxz*uixp - tkyz*uiyp - tkzz*uizp
                  frcx = depx
                  frcy = depy
                  frcz = depz
c
c     get the dEp/dR terms used for direct polarization force
c
                  term1 = dmpe(5) - psc3*rr5
                  term2 = dmpe(7) - psc5*rr7
                  term3 = -psr3 + term1*xr*xr - rr3*xr*prc3(1)
                  term4 = rr3*prc3(1) - term1*xr - psr5*xr
                  term5 = term2*xr*xr - psr5 - rr5*xr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*xr - dmpe(7)
     &                       - rr7*xr*prc7(1)
                  term7 = rr5*prc5(1) - 2.0d0*dmpe(7)*xr
     &                       + (psc5+1.5d0*psc7)*rr7*xr
                  tixx = ci*term3 + dix*term4 + dir*term5
     &                      + 2.0d0*psr5*qixx + (qiy*yr+qiz*zr)*psc7*rr7
     &                      + 2.0d0*qix*term7 + qir*term6
                  tkxx = ck*term3 - dkx*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxx + (qky*yr+qkz*zr)*psc7*rr7
     &                      + 2.0d0*qkx*term7 + qkr*term6
                  term3 = -psr3 + term1*yr*yr - rr3*yr*prc3(2)
                  term4 = rr3*prc3(2) - term1*yr - psr5*yr
                  term5 = term2*yr*yr - psr5 - rr5*yr*prc5(2)
                  term6 = (dmpe(9)-psc7*rr9)*yr*yr - dmpe(7)
     &                       - rr7*yr*prc7(2)
                  term7 = rr5*prc5(2) - 2.0d0*dmpe(7)*yr
     &                       + (psc5+1.5d0*psc7)*rr7*yr
                  tiyy = ci*term3 + diy*term4 + dir*term5
     &                      + 2.0d0*psr5*qiyy + (qix*xr+qiz*zr)*psc7*rr7
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkyy = ck*term3 - dky*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkyy + (qkx*xr+qkz*zr)*psc7*rr7
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = -psr3 + term1*zr*zr - rr3*zr*prc3(3)
                  term4 = rr3*prc3(3) - term1*zr - psr5*zr
                  term5 = term2*zr*zr - psr5 - rr5*zr*prc5(3)
                  term6 = (dmpe(9)-psc7*rr9)*zr*zr - dmpe(7)
     &                       - rr7*zr*prc7(3)
                  term7 = rr5*prc5(3) - 2.0d0*dmpe(7)*zr
     &                       + (psc5+1.5d0*psc7)*rr7*zr
                  tizz = ci*term3 + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qizz + (qix*xr+qiy*yr)*psc7*rr7
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkzz = ck*term3 - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkzz + (qkx*xr+qky*yr)*psc7*rr7
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*xr*yr - rr3*yr*prc3(1)
                  term4 = rr3*prc3(1) - term1*xr
                  term5 = term2*xr*yr - rr5*yr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*yr - rr7*yr*prc7(1)
                  term7 = rr5*prc5(1) - term2*xr
                  tixy = ci*term3 - psr5*dix*yr + diy*term4 + dir*term5
     &                      + 2.0d0*psr5*qixy - 2.0d0*psr7*yr*qix
     &                      + 2.0d0*qiy*term7 + qir*term6
                  tkxy = ck*term3 + psr5*dkx*yr - dky*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxy - 2.0d0*psr7*yr*qkx
     &                      + 2.0d0*qky*term7 + qkr*term6
                  term3 = term1*xr*zr - rr3*zr*prc3(1)
                  term5 = term2*xr*zr - rr5*zr*prc5(1)
                  term6 = (dmpe(9)-psc7*rr9)*xr*zr - rr7*zr*prc7(1)
                  tixz = ci*term3 - psr5*dix*zr + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qixz - 2.0d0*psr7*zr*qix
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkxz = ck*term3 + psr5*dkx*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkxz - 2.0d0*psr7*zr*qkx
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  term3 = term1*yr*zr - rr3*zr*prc3(2)
                  term4 = rr3*prc3(2) - term1*yr
                  term5 = term2*yr*zr - rr5*zr*prc5(2)
                  term6 = (dmpe(9)-psc7*rr9)*yr*zr - rr7*zr*prc7(2)
                  term7 = rr5*prc5(2) - term2*yr
                  tiyz = ci*term3 - psr5*diy*zr + diz*term4 + dir*term5
     &                      + 2.0d0*psr5*qiyz - 2.0d0*psr7*zr*qiy
     &                      + 2.0d0*qiz*term7 + qir*term6
                  tkyz = ck*term3 + psr5*dky*zr - dkz*term4 - dkr*term5
     &                      + 2.0d0*psr5*qkyz - 2.0d0*psr7*zr*qky
     &                      + 2.0d0*qkz*term7 + qkr*term6
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
c
c     get the field gradient for direct polarization force
c
               else if (use_chgpen) then
                  term1i = rr3i - rr5i*xr*xr
                  term1core = rr3core - rr5core*xr*xr
                  term2i = 2.0d0*rr5i*xr 
                  term3i = rr7i*xr*xr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*xr
                  term6i = rr9i*xr*xr
                  term1k = rr3k - rr5k*xr*xr
                  term2k = 2.0d0*rr5k*xr
                  term3k = rr7k*xr*xr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*xr
                  term6k = rr9k*xr*xr
                  tixx = vali*term1i + corei*term1core  
     &                      + dix*term2i - dir*term3i
     &                      - qixx*term4i + qix*term5i - qir*term6i
     &                      + (qiy*yr+qiz*zr)*rr7i
                  tkxx = valk*term1k + corek*term1core
     &                      - dkx*term2k + dkr*term3k
     &                      - qkxx*term4k + qkx*term5k - qkr*term6k
     &                      + (qky*yr+qkz*zr)*rr7k
                  term1i = rr3i - rr5i*yr*yr
                  term1core = rr3core - rr5core*yr*yr
                  term2i = 2.0d0*rr5i*yr
                  term3i = rr7i*yr*yr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*yr
                  term6i = rr9i*yr*yr
                  term1k = rr3k - rr5k*yr*yr
                  term2k = 2.0d0*rr5k*yr
                  term3k = rr7k*yr*yr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*yr
                  term6k = rr9k*yr*yr
                  tiyy = vali*term1i + corei*term1core
     &                      + diy*term2i - dir*term3i
     &                      - qiyy*term4i + qiy*term5i - qir*term6i
     &                      + (qix*xr+qiz*zr)*rr7i
                  tkyy = valk*term1k + corek*term1core
     &                      - dky*term2k + dkr*term3k
     &                      - qkyy*term4k + qky*term5k - qkr*term6k
     &                      + (qkx*xr+qkz*zr)*rr7k
                  term1i = rr3i - rr5i*zr*zr
                  term1core = rr3core - rr5core*zr*zr
                  term2i = 2.0d0*rr5i*zr
                  term3i = rr7i*zr*zr - rr5i
                  term4i = 2.0d0*rr5i
                  term5i = 5.0d0*rr7i*zr
                  term6i = rr9i*zr*zr
                  term1k = rr3k - rr5k*zr*zr
                  term2k = 2.0d0*rr5k*zr
                  term3k = rr7k*zr*zr - rr5k
                  term4k = 2.0d0*rr5k
                  term5k = 5.0d0*rr7k*zr
                  term6k = rr9k*zr*zr
                  tizz = vali*term1i + corei*term1core
     &                      + diz*term2i - dir*term3i
     &                      - qizz*term4i + qiz*term5i - qir*term6i
     &                      + (qix*xr+qiy*yr)*rr7i
                  tkzz = valk*term1k + corek*term1core
     &                      - dkz*term2k + dkr*term3k
     &                      - qkzz*term4k + qkz*term5k - qkr*term6k
     &                      + (qkx*xr+qky*yr)*rr7k
                  term2i = rr5i*xr 
                  term1i = yr * term2i
                  term1core = rr5core*xr*yr
                  term3i = rr5i*yr
                  term4i = yr * (rr7i*xr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*xr
                  term7i = 2.0d0*rr7i*yr
                  term8i = yr*rr9i*xr
                  term2k = rr5k*xr
                  term1k = yr * term2k
                  term3k = rr5k*yr
                  term4k = yr * (rr7k*xr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*xr
                  term7k = 2.0d0*rr7k*yr
                  term8k = yr*rr9k*xr
                  tixy = -vali*term1i - corei*term1core 
     &                      + diy*term2i + dix*term3i
     &                      - dir*term4i - qixy*term5i + qiy*term6i
     &                      + qix*term7i - qir*term8i
                  tkxy = -valk*term1k - corek*term1core 
     &                      - dky*term2k - dkx*term3k
     &                      + dkr*term4k - qkxy*term5k + qky*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5i*xr
                  term1i = zr * term2i
                  term1core = rr5core*xr*zr
                  term3i = rr5i*zr
                  term4i = zr * (rr7i*xr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*xr
                  term7i = 2.0d0*rr7i*zr
                  term8i = zr*rr9i*xr
                  term2k = rr5k*xr
                  term1k = zr * term2k
                  term3k = rr5k*zr
                  term4k = zr * (rr7k*xr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*xr
                  term7k = 2.0d0*rr7k*zr
                  term8k = zr*rr9k*xr
                  tixz = -vali*term1i - corei*term1core
     &                      + diz*term2i + dix*term3i
     &                      - dir*term4i - qixz*term5i + qiz*term6i
     &                      + qix*term7i - qir*term8i
                  tkxz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dkx*term3k
     &                      + dkr*term4k - qkxz*term5k + qkz*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5i*yr
                  term1i = zr * term2i
                  term1core = rr5core*yr*zr
                  term3i = rr5i*zr
                  term4i = zr * (rr7i*yr)
                  term5i = 2.0d0*rr5i
                  term6i = 2.0d0*rr7i*yr
                  term7i = 2.0d0*rr7i*zr
                  term8i = zr*rr9i*yr
                  term2k = rr5k*yr
                  term1k = zr * term2k
                  term3k = rr5k*zr
                  term4k = zr * (rr7k*yr)
                  term5k = 2.0d0*rr5k
                  term6k = 2.0d0*rr7k*yr
                  term7k = 2.0d0*rr7k*zr
                  term8k = zr*rr9k*yr
                  tiyz = -vali*term1i - corei*term1core
     &                      + diz*term2i + diy*term3i
     &                      - dir*term4i - qiyz*term5i + qiz*term6i
     &                      + qiy*term7i - qir*term8i
                  tkyz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dky*term3k
     &                      + dkr*term4k - qkyz*term5k + qkz*term6k
     &                      + qky*term7k - qkr*term8k
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = -2.0d0 * depx
                  frcy = -2.0d0 * depy
                  frcz = -2.0d0 * depz
               end if
c
c     reset Thole values if alternate direct damping was used
c
               if (use_dirdamp) then
                  sc3 = 1.0d0
                  sc5 = 1.0d0
                  do j = 1, 3
                     rc3(j) = 0.0d0
                     rc5(j) = 0.0d0
                  end do
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        sc3 = 1.0d0 - expdamp
                        sc5 = 1.0d0 - expdamp*(1.0d0+damp)
                        temp3 = 3.0d0 * damp * expdamp / r2
                        temp5 = damp
                        rc3(1) = xr * temp3
                        rc3(2) = yr * temp3
                        rc3(3) = zr * temp3
                        rc5(1) = rc3(1) * temp5
                        rc5(2) = rc3(2) * temp5
                        rc5(3) = rc3(3) * temp5
                     end if
                  end if
                  usc3 = 1.0d0 - sc3*uscale(k)
                  usc5 = 1.0d0 - sc5*uscale(k)
                  usr3 = dmpe(3) - usc3*rr3
                  usr5 = dmpe(5) - usc5*rr5
                  do j = 1, 3
                     urc3(j) = rc3(j) * uscale(k)
                     urc5(j) = rc5(j) * uscale(k)
                  end do
               end if
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp.eq.'MUTUAL' .and. use_thole) then
                  term1 = dmpe(5) - usc3*rr5
                  term2 = dmpe(7) - usc5*rr7
                  term3 = usr5 + term1
                  term4 = rr3 * uscale(k)
                  term5 = -xr*term3 + rc3(1)*term4
                  term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                  tixx = uix*term5 + uir*term6
                  tkxx = ukx*term5 + ukr*term6
                  term5 = -yr*term3 + rc3(2)*term4
                  term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                  tiyy = uiy*term5 + uir*term6
                  tkyy = uky*term5 + ukr*term6
                  term5 = -zr*term3 + rc3(3)*term4
                  term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                  tizz = uiz*term5 + uir*term6
                  tkzz = ukz*term5 + ukr*term6
                  term4 = -usr5 * yr
                  term5 = -xr*term1 + rr3*urc3(1)
                  term6 = xr*yr*term2 - rr5*yr*urc5(1)
                  tixy = uix*term4 + uiy*term5 + uir*term6
                  tkxy = ukx*term4 + uky*term5 + ukr*term6
                  term4 = -usr5 * zr
                  term6 = xr*zr*term2 - rr5*zr*urc5(1)
                  tixz = uix*term4 + uiz*term5 + uir*term6
                  tkxz = ukx*term4 + ukz*term5 + ukr*term6
                  term5 = -yr*term1 + rr3*urc3(2)
                  term6 = yr*zr*term2 - rr5*zr*urc5(2)
                  tiyz = uiy*term4 + uiz*term5 + uir*term6
                  tkyz = uky*term4 + ukz*term5 + ukr*term6
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + depx
                  frcy = frcy + depy
                  frcz = frcz + depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               else if (poltyp.eq.'MUTUAL' .and. use_chgpen) then
                  term1 = 2.0d0 * rr5ik
                  term2 = term1*xr
                  term3 = rr5ik - rr7ik*xr*xr
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr
                  term3 = rr5ik - rr7ik*yr*yr
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr
                  term3 = rr5ik - rr7ik*zr*zr
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = rr5ik*yr
                  term2 = rr5ik*xr
                  term3 = yr * (rr7ik*xr)
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = rr5ik * zr
                  term3 = zr * (rr7ik*xr)
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = rr5ik*yr
                  term3 = zr * (rr7ik*yr)
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx - depx
                  frcy = frcy - depy
                  frcz = frcz - depz
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_thole) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,ii)*xr + uopt(j,2,ii)*yr
     &                          + uopt(j,3,ii)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,kk)*xr + uopt(m,2,kk)*yr
     &                             + uopt(m,3,kk)*zr
                        term1 = dmpe(5) - usc3*rr5
                        term2 = dmpe(7) - usc5*rr7
                        term3 = usr5 + term1
                        term4 = rr3 * uscale(k)
                        term5 = -xr*term3 + rc3(1)*term4
                        term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                        tixx = uopt(j,1,ii)*term5 + uirm*term6
                        tkxx = uopt(m,1,kk)*term5 + ukrm*term6
                        term5 = -yr*term3 + rc3(2)*term4
                        term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                        tiyy = uopt(j,2,ii)*term5 + uirm*term6
                        tkyy = uopt(m,2,kk)*term5 + ukrm*term6
                        term5 = -zr*term3 + rc3(3)*term4
                        term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                        tizz = uopt(j,3,ii)*term5 + uirm*term6
                        tkzz = uopt(m,3,kk)*term5 + ukrm*term6
                        term4 = -usr5 * yr
                        term5 = -xr*term1 + rr3*urc3(1)
                        term6 = xr*yr*term2 - rr5*yr*urc5(1)
                        tixy = uopt(j,1,ii)*term4 + uopt(j,2,ii)*term5
     &                            + uirm*term6
                        tkxy = uopt(m,1,kk)*term4 + uopt(m,2,kk)*term5
     &                            + ukrm*term6
                        term4 = -usr5 * zr
                        term6 = xr*zr*term2 - rr5*zr*urc5(1)
                        tixz = uopt(j,1,ii)*term4 + uopt(j,3,ii)*term5
     &                            + uirm*term6
                        tkxz = uopt(m,1,kk)*term4 + uopt(m,3,kk)*term5
     &                            + ukrm*term6
                        term5 = -yr*term1 + rr3*urc3(2)
                        term6 = yr*zr*term2 - rr5*zr*urc5(2)
                        tiyz = uopt(j,2,ii)*term4 + uopt(j,3,ii)*term5
     &                            + uirm*term6
                        tkyz = uopt(m,2,kk)*term4 + uopt(m,3,kk)*term5
     &                            + ukrm*term6
                        depx = tixx*uoptp(m,1,kk) + tkxx*uoptp(j,1,ii)
     &                       + tixy*uoptp(m,2,kk) + tkxy*uoptp(j,2,ii)
     &                       + tixz*uoptp(m,3,kk) + tkxz*uoptp(j,3,ii)
                        depy = tixy*uoptp(m,1,kk) + tkxy*uoptp(j,1,ii)
     &                       + tiyy*uoptp(m,2,kk) + tkyy*uoptp(j,2,ii)
     &                       + tiyz*uoptp(m,3,kk) + tkyz*uoptp(j,3,ii)
                        depz = tixz*uoptp(m,1,kk) + tkxz*uoptp(j,1,ii)
     &                       + tiyz*uoptp(m,2,kk) + tkyz*uoptp(j,2,ii)
     &                       + tizz*uoptp(m,3,kk) + tkzz*uoptp(j,3,ii)
                        frcx = frcx + copm(j+m+1)*depx
                        frcy = frcy + copm(j+m+1)*depy
                        frcz = frcz + copm(j+m+1)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp.eq.'OPT' .and. use_chgpen) then
                  do j = 0, optorder-1
                     uirm = uopt(j,1,i)*xr + uopt(j,2,i)*yr
     &                          + uopt(j,3,i)*zr
                     do m = 0, optorder-j-1
                        ukrm = uopt(m,1,k)*xr + uopt(m,2,k)*yr
     &                             + uopt(m,3,k)*zr
                        term1 = 2.0d0 * rr5ik
                        term2 = term1*xr
                        term3 = rr5ik - rr7ik*xr*xr
                        tixx = uopt(j,1,i)*term2 + uirm*term3
                        tkxx = uopt(m,1,k)*term2 + ukrm*term3
                        term2 = term1*yr
                        term3 = rr5ik - rr7ik*yr*yr
                        tiyy = uopt(j,2,i)*term2 + uirm*term3
                        tkyy = uopt(m,2,k)*term2 + ukrm*term3
                        term2 = term1*zr
                        term3 = rr5ik - rr7ik*zr*zr
                        tizz = uopt(j,3,i)*term2 + uirm*term3
                        tkzz = uopt(m,3,k)*term2 + ukrm*term3
                        term1 = rr5ik*yr
                        term2 = rr5ik*xr
                        term3 = yr * (rr7ik*xr)
                        tixy = uopt(j,1,i)*term1 + uopt(j,2,i)*term2 
     &                       - uirm*term3
                        tkxy = uopt(m,1,k)*term1 + uopt(m,2,k)*term2 
     &                       - ukrm*term3
                        term1 = rr5ik * zr
                        term3 = zr * (rr7ik*xr)
                        tixz = uopt(j,1,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        term2 = rr5ik*yr
                        term3 = zr * (rr7ik*yr)
                        tiyz = uopt(j,2,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,k) + tkxx*uoptp(j,1,i)
     &                       + tixy*uoptp(m,2,k) + tkxy*uoptp(j,2,i)
     &                       + tixz*uoptp(m,3,k) + tkxz*uoptp(j,3,i)
                        depy = tixy*uoptp(m,1,k) + tkxy*uoptp(j,1,i)
     &                       + tiyy*uoptp(m,2,k) + tkyy*uoptp(j,2,i)
     &                       + tiyz*uoptp(m,3,k) + tkyz*uoptp(j,3,i)
                        depz = tixz*uoptp(m,1,k) + tkxz*uoptp(j,1,i)
     &                       + tiyz*uoptp(m,2,k) + tkyz*uoptp(j,2,i)
     &                       + tizz*uoptp(m,3,k) + tkzz*uoptp(j,3,i)
                        frcx = frcx - copm(j+m+1)*depx
                        frcy = frcy - copm(j+m+1)*depy
                        frcz = frcz - copm(j+m+1)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for TCG polarization force
c
               else if (poltyp.eq.'TCG' .and. use_thole) then
                  do j = 1, tcgnab
                     ukx = ubd(1,kk,j)
                     uky = ubd(2,kk,j)
                     ukz = ubd(3,kk,j)
                     ukxp = ubp(1,kk,j)
                     ukyp = ubp(2,kk,j)
                     ukzp = ubp(3,kk,j)
                     uirt = uax(j)*xr + uay(j)*yr + uaz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = dmpe(5) - usc3*rr5
                     term2 = dmpe(7) - usc5*rr7
                     term3 = usr5 + term1
                     term4 = rr3 * uscale(k)
                     term5 = -xr*term3 + rc3(1)*term4
                     term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                     tixx = uax(j)*term5 + uirt*term6
                     tkxx = ukx*term5 + ukrt*term6
                     term5 = -yr*term3 + rc3(2)*term4
                     term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                     tiyy = uay(j)*term5 + uirt*term6
                     tkyy = uky*term5 + ukrt*term6
                     term5 = -zr*term3 + rc3(3)*term4
                     term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                     tizz = uaz(j)*term5 + uirt*term6
                     tkzz = ukz*term5 + ukrt*term6
                     term4 = -usr5 * yr
                     term5 = -xr*term1 + rr3*urc3(1)
                     term6 = xr*yr*term2 - rr5*yr*urc5(1)
                     tixy = uax(j)*term4 + uay(j)*term5 + uirt*term6
                     tkxy = ukx*term4 + uky*term5 + ukrt*term6
                     term4 = -usr5 * zr
                     term6 = xr*zr*term2 - rr5*zr*urc5(1)
                     tixz = uax(j)*term4 + uaz(j)*term5 + uirt*term6
                     tkxz = ukx*term4 + ukz*term5 + ukrt*term6
                     term5 = -yr*term1 + rr3*urc3(2)
                     term6 = yr*zr*term2 - rr5*zr*urc5(2)
                     tiyz = uay(j)*term4 + uaz(j)*term5 + uirt*term6
                     tkyz = uky*term4 + ukz*term5 + ukrt*term6
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*uaxp(j) + tkxy*uayp(j)
     &                         + tkxz*uazp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*uaxp(j) + tkyy*uayp(j)
     &                         + tkyz*uazp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*uaxp(j) + tkyz*uayp(j)
     &                         + tkzz*uazp(j)
                     frcx = frcx + depx
                     frcy = frcy + depy
                     frcz = frcz + depz
                     ukx = uad(1,kk,j)
                     uky = uad(2,kk,j)
                     ukz = uad(3,kk,j)
                     ukxp = uap(1,kk,j)
                     ukyp = uap(2,kk,j)
                     ukzp = uap(3,kk,j)
                     uirt = ubx(j)*xr + uby(j)*yr + ubz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = dmpe(5) - usc3*rr5
                     term2 = dmpe(7) - usc5*rr7
                     term3 = usr5 + term1
                     term4 = rr3 * uscale(k)
                     term5 = -xr*term3 + rc3(1)*term4
                     term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5(1)
                     tixx = ubx(j)*term5 + uirt*term6
                     tkxx = ukx*term5 + ukrt*term6
                     term5 = -yr*term3 + rc3(2)*term4
                     term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5(2)
                     tiyy = uby(j)*term5 + uirt*term6
                     tkyy = uky*term5 + ukrt*term6
                     term5 = -zr*term3 + rc3(3)*term4
                     term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5(3)
                     tizz = ubz(j)*term5 + uirt*term6
                     tkzz = ukz*term5 + ukrt*term6
                     term4 = -usr5 * yr
                     term5 = -xr*term1 + rr3*urc3(1)
                     term6 = xr*yr*term2 - rr5*yr*urc5(1)
                     tixy = ubx(j)*term4 + uby(j)*term5 + uirt*term6
                     tkxy = ukx*term4 + uky*term5 + ukrt*term6
                     term4 = -usr5 * zr
                     term6 = xr*zr*term2 - rr5*zr*urc5(1)
                     tixz = ubx(j)*term4 + ubz(j)*term5 + uirt*term6
                     tkxz = ukx*term4 + ukz*term5 + ukrt*term6
                     term5 = -yr*term1 + rr3*urc3(2)
                     term6 = yr*zr*term2 - rr5*zr*urc5(2)
                     tiyz = uby(j)*term4 + ubz(j)*term5 + uirt*term6
                     tkyz = uky*term4 + ukz*term5 + ukrt*term6
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*ubxp(j) + tkxy*ubyp(j)
     &                         + tkxz*ubzp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*ubxp(j) + tkyy*ubyp(j)
     &                         + tkyz*ubzp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*ubxp(j) + tkyz*ubyp(j)
     &                         + tkzz*ubzp(j)
                     frcx = frcx + depx
                     frcy = frcy + depy
                     frcz = frcz + depz
                  end do
               end if
c
c     increment force-based gradient on the interaction sites
c
               dep(1,i) = dep(1,i) - frcx
               dep(2,i) = dep(2,i) - frcy
               dep(3,i) = dep(3,i) - frcz
               dep(1,k) = dep(1,k) + frcx
               dep(2,k) = dep(2,k) + frcy
               dep(3,k) = dep(3,k) + frcz
c
c     increment the virial due to pairwise Cartesian forces
c
               vxx = xr * frcx
               vxy = 0.5d0 * (yr*frcx+xr*frcy)
               vxz = 0.5d0 * (zr*frcx+xr*frcz)
               vyy = yr * frcy
               vyz = 0.5d0 * (zr*frcy+yr*frcz)
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
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               uscale(ip14(j,i)) = 1.0d0
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
               uscale(ip14(j,i)) = 1.0d0
            end do
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP DO reduction(+:dep,vir) schedule(guided)
c
c     torque is induced field and gradient cross permanent moments
c
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         tep(1) = diz*ufld(2,i) - diy*ufld(3,i)
     &               + qixz*dufld(2,i) - qixy*dufld(4,i)
     &               + 2.0d0*qiyz*(dufld(3,i)-dufld(6,i))
     &               + (qizz-qiyy)*dufld(5,i)
         tep(2) = dix*ufld(3,i) - diz*ufld(1,i)
     &               - qiyz*dufld(2,i) + qixy*dufld(5,i)
     &               + 2.0d0*qixz*(dufld(6,i)-dufld(1,i))
     &               + (qixx-qizz)*dufld(4,i)
         tep(3) = diy*ufld(1,i) - dix*ufld(2,i)
     &               + qiyz*dufld(4,i) - qixz*dufld(5,i)
     &               + 2.0d0*qixy*(dufld(1,i)-dufld(3,i))
     &               + (qiyy-qixx)*dufld(2,i)
         call torque (ii,tep,fix,fiy,fiz,dep)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = abs(yaxis(ii))
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
!$OMP    DO reduction(+:dep,vir) schedule(guided)
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            frcx = decfx(i)
            frcy = decfy(i)
            frcz = decfz(i)
            dep(1,i) = dep(1,i) + frcx
            dep(2,i) = dep(2,i) + frcy
            dep(3,i) = dep(3,i) + frcz
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
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (wscale)
      deallocate (ufld)
      deallocate (dufld)
      deallocate (pot)
      deallocate (decfx)
      deallocate (decfy)
      deallocate (decfz)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar1e  --  single-loop polarization energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epreal1e" calculates the induced dipole polarization energy
c     from the induced dipoles times the electric field
c
c
      subroutine epolar1e
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use limits
      use math
      use mpole
      use polar
      use polpot
      implicit none
      integer i,j,ii
      real*8 e,f,fi,term
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
c
c
c     set the energy unit conversion factor
c
      f = -0.5d0 * electric / dielec
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(ii,j,fi,e)
!$OMP DO reduction(+:ep) schedule(guided)
c
c     get polarization energy via induced dipoles times field
c
      do ii = 1, npole
         if (douind(ipole(ii))) then
            fi = f / polarity(ii)
            e = 0.0d0
            do j = 1, 3
               e = e + fi*uind(j,ii)*udirp(j,ii)
            end do
            ep = ep + e
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     compute the cell dipole boundary correction term
c
      if (use_ewald) then
         if (boundary .eq. 'VACUUM') then
            f = electric / dielec
            xd = 0.0d0
            yd = 0.0d0
            zd = 0.0d0
            xu = 0.0d0
            yu = 0.0d0
            zu = 0.0d0
            do ii = 1, npole
               i = ipole(ii)
               dix = rpole(2,ii)
               diy = rpole(3,ii)
               diz = rpole(4,ii)
               uix = uind(1,ii)
               uiy = uind(2,ii)
               uiz = uind(3,ii)
               xd = xd + dix + rpole(1,ii)*x(i)
               yd = yd + diy + rpole(1,ii)*y(i)
               zd = zd + diz + rpole(1,ii)*z(i)
               xu = xu + uix
               yu = yu + uiy
               zu = zu + uiz
            end do
            term = (2.0d0/3.0d0) * f * (pi/volbox)
            e = term * (xd*xu+yd*yu+zd*zu)
            ep = ep + e
         end if
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine eprecip1  --  PME recip polarize energy & derivs  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "eprecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to dipole polarization
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
      use atoms
      use bound
      use boxes
      use chgpot
      use deriv
      use ewald
      use math
      use mpole
      use mrecip
      use pme
      use polar
      use polopt
      use polpot
      use poltcg
      use potent
      use virial
      implicit none
      integer i,j,k,m,ii
      integer j1,j2,j3
      integer k1,k2,k3
      integer m1,m2,m3
      integer ix,iy,iz
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 eterm,f
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
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
      real*8 tep(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8 cphid(4),cphip(4)
      real*8 a(3,3),ftc(10,10)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: fphid(:,:)
      real*8, allocatable :: fphip(:,:)
      real*8, allocatable :: fphidp(:,:)
      real*8, allocatable :: cphidp(:,:)
      real*8, allocatable :: qgrip(:,:,:,:)
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
c     initialize variables required for the scalar summation
c
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      nff = nfft1 * nfft2
      ntot = nff * nfft3
c
c     remove scalar sum virial from prior multipole FFT
c
      if (use_mpole .and. aewald.eq.aeewald) then
         vxx = -vmxx
         vxy = -vmxy
         vxz = -vmxz
         vyy = -vmyy
         vyz = -vmyz
         vzz = -vmzz
c
c     perform dynamic allocation of some global arrays
c
      else
         if (allocated(cmp)) then
            if (size(cmp) .lt. 10*npole)  deallocate (cmp)
         end if
         if (allocated(fmp)) then
            if (size(fmp) .lt. 10*npole)  deallocate (fmp)
         end if
         if (allocated(cphi)) then
            if (size(cphi) .lt. 10*npole) deallocate (cphi)
         end if
         if (allocated(fphi)) then
            if (size(fphi) .lt. 20*npole)  deallocate (fphi)
         end if
         if (.not. allocated(cmp))  allocate (cmp(10,npole))
         if (.not. allocated(fmp))  allocate (fmp(10,npole))
         if (.not. allocated(cphi))  allocate (cphi(10,npole))
         if (.not. allocated(fphi))  allocate (fphi(20,npole))
c
c     perform dynamic allocation of some global arrays
c
         ntot = nfft1 * nfft2 * nfft3
         if (allocated(qgrid)) then
            if (size(qgrid) .ne. 2*ntot)  call fftclose
         end if
         if (allocated(qfac)) then
            if (size(qfac) .ne. ntot)  deallocate (qfac)
         end if
         if (.not. allocated(qgrid))  call fftsetup
         if (.not. allocated(qfac))  allocate (qfac(nfft1,nfft2,nfft3))
c
c     setup spatial decomposition and B-spline coefficients
c
         call getchunk
         call moduli
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
         qfac(1,1,1) = 0.0d0
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
               vxx = vxx - h1*h1*vterm + eterm
               vxy = vxy - h1*h2*vterm
               vxz = vxz - h1*h3*vterm
               vyy = vyy - h2*h2*vterm + eterm
               vyz = vyz - h2*h3*vterm
               vzz = vzz - h3*h3*vterm + eterm
            end if
            qfac(k1,k2,k3) = expterm
         end do
c
c     account for zeroth grid point for nonperiodic system
c
         if (.not. use_bounds) then
            expterm = 0.5d0 * pi / xbox
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
               fphi(j,i) = f * fphi(j,i)
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
      allocate (cphidp(10,npole))
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
         do j = 2, 10
            fphid(j,i) = f * fphid(j,i)
            fphip(j,i) = f * fphip(j,i)
         end do
         do j = 1, 20
            fphidp(j,i) = f * fphidp(j,i)
         end do
      end do
c
c     increment the dipole polarization gradient contributions
c
      do i = 1, npole
         ii = ipole(i)
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 3
            j1 = deriv1(k+1)
            j2 = deriv2(k+1)
            j3 = deriv3(k+1)
            f1 = f1 + (fuind(k,i)+fuinp(k,i))*fphi(j1,i)
            f2 = f2 + (fuind(k,i)+fuinp(k,i))*fphi(j2,i)
            f3 = f3 + (fuind(k,i)+fuinp(k,i))*fphi(j3,i)
            if (poltyp .eq. 'MUTUAL') then
               f1 = f1 + fuind(k,i)*fphip(j1,i) + fuinp(k,i)*fphid(j1,i)
               f2 = f2 + fuind(k,i)*fphip(j2,i) + fuinp(k,i)*fphid(j2,i)
               f3 = f3 + fuind(k,i)*fphip(j3,i) + fuinp(k,i)*fphid(j3,i)
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
         dep(1,ii) = dep(1,ii) + h1
         dep(2,ii) = dep(2,ii) + h2
         dep(3,ii) = dep(3,ii) + h3
      end do
c
c     set the potential to be the induced dipole average
c
      do i = 1, npole
         do j = 1, 10
            fphidp(j,i) = 0.5d0 * fphidp(j,i)
         end do
      end do
      call fphi_to_cphi (fphidp,cphidp)
c
c     get the fractional to Cartesian transformation matrix
c
      call frac_to_cart (ftc)
c
c     increment the dipole polarization virial contributions
c
      do i = 1, npole
         do j = 2, 4
            cphid(j) = 0.0d0
            cphip(j) = 0.0d0
            do k = 2, 4
               cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
               cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
            end do
         end do
         vxx = vxx - cmp(2,i)*cphidp(2,i)
     &            - 0.5d0*((uind(1,i)+uinp(1,i))*cphi(2,i))
         vxy = vxy - 0.5d0*(cphidp(2,i)*cmp(3,i)+cphidp(3,i)*cmp(2,i))
     &            - 0.25d0*((uind(2,i)+uinp(2,i))*cphi(2,i)
     &                     +(uind(1,i)+uinp(1,i))*cphi(3,i))
         vxz = vxz - 0.5d0*(cphidp(2,i)*cmp(4,i)+cphidp(4,i)*cmp(2,i))
     &            - 0.25d0*((uind(3,i)+uinp(3,i))*cphi(2,i)
     &                     +(uind(1,i)+uinp(1,i))*cphi(4,i))
         vyy = vyy - cmp(3,i)*cphidp(3,i)
     &            - 0.5d0*((uind(2,i)+uinp(2,i))*cphi(3,i))
         vyz = vyz - 0.5d0*(cphidp(3,i)*cmp(4,i)+cphidp(4,i)*cmp(3,i))
     &            - 0.25d0*((uind(3,i)+uinp(3,i))*cphi(3,i)
     &                     +(uind(2,i)+uinp(2,i))*cphi(4,i))
         vzz = vzz - cmp(4,i)*cphidp(4,i)
     &            - 0.5d0*((uind(3,i)+uinp(3,i))*cphi(4,i))
         vxx = vxx - 2.0d0*cmp(5,i)*cphidp(5,i)
     &            - cmp(8,i)*cphidp(8,i) - cmp(9,i)*cphidp(9,i)
         vxy = vxy - (cmp(5,i)+cmp(6,i))*cphidp(8,i)
     &            - 0.5d0*(cmp(8,i)*(cphidp(6,i)+cphidp(5,i))
     &                 +cmp(9,i)*cphidp(10,i)+cmp(10,i)*cphidp(9,i))
         vxz = vxz - (cmp(5,i)+cmp(7,i))*cphidp(9,i)
     &            - 0.5d0*(cmp(9,i)*(cphidp(5,i)+cphidp(7,i))
     &                 +cmp(8,i)*cphidp(10,i)+cmp(10,i)*cphidp(8,i))
         vyy = vyy - 2.0d0*cmp(6,i)*cphidp(6,i)
     &            - cmp(8,i)*cphidp(8,i) - cmp(10,i)*cphidp(10,i)
         vyz = vyz - (cmp(6,i)+cmp(7,i))*cphidp(10,i)
     &            - 0.5d0*(cmp(10,i)*(cphidp(6,i)+cphidp(7,i))
     &                 +cmp(8,i)*cphidp(9,i)+cmp(9,i)*cphidp(8,i))
         vzz = vzz - 2.0d0*cmp(7,i)*cphidp(7,i)
     &            - cmp(9,i)*cphidp(9,i) - cmp(10,i)*cphidp(10,i)
         if (poltyp .eq. 'MUTUAL') then
            vxx = vxx - 0.5d0*(cphid(2)*uinp(1,i)+cphip(2)*uind(1,i))
            vxy = vxy - 0.25d0*(cphid(2)*uinp(2,i)+cphip(2)*uind(2,i)
     &                         +cphid(3)*uinp(1,i)+cphip(3)*uind(1,i))
            vxz = vxz - 0.25d0*(cphid(2)*uinp(3,i)+cphip(2)*uind(3,i)
     &                         +cphid(4)*uinp(1,i)+cphip(4)*uind(1,i))
            vyy = vyy - 0.5d0*(cphid(3)*uinp(2,i)+cphip(3)*uind(2,i))
            vyz = vyz - 0.25d0*(cphid(3)*uinp(3,i)+cphip(3)*uind(3,i)
     &                         +cphid(4)*uinp(2,i)+cphip(4)*uind(2,i))
            vzz = vzz - 0.5d0*(cphid(4)*uinp(3,i)+cphip(4)*uind(3,i))
         end if
      end do
c
c     resolve site torques then increment forces and virial
c
      do i = 1, npole
         tep(1) = cmp(4,i)*cphidp(3,i) - cmp(3,i)*cphidp(4,i)
     &               + 2.0d0*(cmp(7,i)-cmp(6,i))*cphidp(10,i)
     &               + cmp(9,i)*cphidp(8,i) + cmp(10,i)*cphidp(6,i)
     &               - cmp(8,i)*cphidp(9,i) - cmp(10,i)*cphidp(7,i)
         tep(2) = cmp(2,i)*cphidp(4,i) - cmp(4,i)*cphidp(2,i)
     &               + 2.0d0*(cmp(5,i)-cmp(7,i))*cphidp(9,i)
     &               + cmp(8,i)*cphidp(10,i) + cmp(9,i)*cphidp(7,i)
     &               - cmp(9,i)*cphidp(5,i) - cmp(10,i)*cphidp(8,i)
         tep(3) = cmp(3,i)*cphidp(2,i) - cmp(2,i)*cphidp(3,i)
     &               + 2.0d0*(cmp(6,i)-cmp(5,i))*cphidp(8,i)
     &               + cmp(8,i)*cphidp(5,i) + cmp(10,i)*cphidp(9,i)
     &               - cmp(8,i)*cphidp(6,i) - cmp(9,i)*cphidp(10,i)
         call torque (i,tep,fix,fiy,fiz,dep)
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = abs(yaxis(i))
         if (iz .eq. 0)  iz = ii
         if (ix .eq. 0)  ix = ii
         if (iy .eq. 0)  iy = ii
         xiz = x(iz) - x(ii)
         yiz = y(iz) - y(ii)
         ziz = z(iz) - z(ii)
         xix = x(ix) - x(ii)
         yix = y(ix) - y(ii)
         zix = z(ix) - z(ii)
         xiy = x(iy) - x(ii)
         yiy = y(iy) - y(ii)
         ziy = z(iy) - z(ii)
         vxx = vxx + xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
         vxy = vxy + 0.5d0*(yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                    + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
         vxz = vxz + 0.5d0*(zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                    + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3))
         vyy = vyy + yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
         vyz = vyz + 0.5d0*(zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                    + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
         vzz = vzz + zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
      end do
c
c     account for dipole response terms in the OPT method
c
      if (poltyp .eq. 'OPT') then
         do i = 1, npole
            ii = ipole(i)
            do k = 0, optorder-1
               do j = 2, 10
                  fphid(j,i) = f * fopt(k,j,i)
                  fphip(j,i) = f * foptp(k,j,i)
               end do
               do m = 0, optorder-k-1
                  do j = 1, 3
                     fuind(j,i) = a(j,1)*uopt(m,1,i)
     &                               + a(j,2)*uopt(m,2,i)
     &                               + a(j,3)*uopt(m,3,i)
                     fuinp(j,i) = a(j,1)*uoptp(m,1,i)
     &                               + a(j,2)*uoptp(m,2,i)
     &                               + a(j,3)*uoptp(m,3,i)
                  end do
                  f1 = 0.0d0
                  f2 = 0.0d0
                  f3 = 0.0d0
                  do j = 1, 3
                     j1 = deriv1(j+1)
                     j2 = deriv2(j+1)
                     j3 = deriv3(j+1)
                     f1 = f1 + fuind(j,i)*fphip(j1,i)
     &                       + fuinp(j,i)*fphid(j1,i)
                     f2 = f2 + fuind(j,i)*fphip(j2,i)
     &                       + fuinp(j,i)*fphid(j2,i)
                     f3 = f3 + fuind(j,i)*fphip(j3,i)
     &                       + fuinp(j,i)*fphid(j3,i)
                  end do
                  f1 = 0.5d0 * dble(nfft1) * f1
                  f2 = 0.5d0 * dble(nfft2) * f2
                  f3 = 0.5d0 * dble(nfft3) * f3
                  h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
                  h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
                  h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
                  dep(1,ii) = dep(1,ii) + copm(k+m+1)*h1
                  dep(2,ii) = dep(2,ii) + copm(k+m+1)*h2
                  dep(3,ii) = dep(3,ii) + copm(k+m+1)*h3
                  do j = 2, 4
                     cphid(j) = 0.0d0
                     cphip(j) = 0.0d0
                     do j1 = 2, 4
                        cphid(j) = cphid(j) + ftc(j,j1)*fphid(j1,i)
                        cphip(j) = cphip(j) + ftc(j,j1)*fphip(j1,i)
                     end do
                  end do
                  vxx = vxx - 0.5d0*copm(k+m+1)
     &                           *(cphid(2)*uoptp(m,1,i)
     &                            +cphip(2)*uopt(m,1,i))
                  vxy = vxy - 0.25d0*copm(k+m+1)
     &                           *(cphid(2)*uoptp(m,2,i)
     &                            +cphip(2)*uopt(m,2,i)
     &                            +cphid(3)*uoptp(m,1,i)
     &                            +cphip(3)*uopt(m,1,i))
                  vxz = vxz - 0.25d0*copm(k+m+1)
     &                           *(cphid(2)*uoptp(m,3,i)
     &                            +cphip(2)*uopt(m,3,i)
     &                            +cphid(4)*uoptp(m,1,i)
     &                            +cphip(4)*uopt(m,1,i))
                  vyy = vyy - 0.5d0*copm(k+m+1)
     &                           *(cphid(3)*uoptp(m,2,i)
     &                            +cphip(3)*uopt(m,2,i))
                  vyz = vyz - 0.25d0*copm(k+m+1)
     &                           *(cphid(3)*uoptp(m,3,i)
     &                            +cphip(3)*uopt(m,3,i)
     &                            +cphid(4)*uoptp(m,2,i)
     &                            +cphip(4)*uopt(m,2,i))
                  vzz = vzz - 0.5d0*copm(k+m+1)
     &                           *(cphid(4)*uoptp(m,3,i)
     &                            +cphip(4)*uopt(m,3,i))
               end do
            end do
         end do
      end if
c
c     account for dipole response terms in the TCG method
c
      if (poltyp .eq. 'TCG') then
         do m = 1, tcgnab
            do i = 1, npole
               do j = 1, 3
                  fuind(j,i) = a(j,1)*uad(1,i,m) + a(j,2)*uad(2,i,m)
     &                            + a(j,3)*uad(3,i,m)
                  fuinp(j,i) = a(j,1)*ubp(1,i,m) + a(j,2)*ubp(2,i,m)
     &                            + a(j,3)*ubp(3,i,m)
               end do
            end do
            call grid_uind (fuind,fuinp)
            call fftfront
            do k = 1, nfft3
               do j = 1, nfft2
                  do i = 1, nfft1
                     term = qfac(i,j,k)
                     qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
                     qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
                  end do
               end do
            end do
            call fftback
            call fphi_uind (fphid,fphip,fphidp)
            do i = 1, npole
               do j = 2, 10
                  fphid(j,i) = f * fphid(j,i)
                  fphip(j,i) = f * fphip(j,i)
               end do
            end do
            do i = 1, npole
               ii = ipole(i)
               f1 = 0.0d0
               f2 = 0.0d0
               f3 = 0.0d0
               do k = 1, 3
                  j1 = deriv1(k+1)
                  j2 = deriv2(k+1)
                  j3 = deriv3(k+1)
                  f1 = f1+fuind(k,i)*fphip(j1,i)+fuinp(k,i)*fphid(j1,i)
                  f2 = f2+fuind(k,i)*fphip(j2,i)+fuinp(k,i)*fphid(j2,i)
                  f3 = f3+fuind(k,i)*fphip(j3,i)+fuinp(k,i)*fphid(j3,i)
               end do
               f1 = 0.5d0 * dble(nfft1) * f1
               f2 = 0.5d0 * dble(nfft2) * f2
               f3 = 0.5d0 * dble(nfft3) * f3
               h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
               h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
               h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
               dep(1,ii) = dep(1,ii) + h1
               dep(2,ii) = dep(2,ii) + h2
               dep(3,ii) = dep(3,ii) + h3
               do j = 2, 4
                  cphid(j) = 0.0d0
                  cphip(j) = 0.0d0
                  do k = 2, 4
                     cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
                     cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
                  end do
               end do
               vxx = vxx - 0.5d0*(cphid(2)*ubp(1,i,m)
     &                              +cphip(2)*uad(1,i,m))
               vxy = vxy - 0.25d0*(cphid(2)*ubp(2,i,m)
     &                               +cphip(2)*uad(2,i,m)
     &                               +cphid(3)*ubp(1,i,m)
     &                               +cphip(3)*uad(1,i,m))
               vxz = vxz - 0.25d0*(cphid(2)*ubp(3,i,m)
     &                               +cphip(2)*uad(3,i,m)
     &                               +cphid(4)*ubp(1,i,m)
     &                               +cphip(4)*uad(1,i,m))
               vyy = vyy - 0.5d0*(cphid(3)*ubp(2,i,m)
     &                              +cphip(3)*uad(2,i,m))
               vyz = vyz - 0.25d0*(cphid(3)*ubp(3,i,m)
     &                               +cphip(3)*uad(3,i,m)
     &                               +cphid(4)*ubp(2,i,m)
     &                               +cphip(4)*uad(2,i,m))
               vzz = vzz - 0.5d0*(cphid(4)*ubp(3,i,m)
     &                              +cphip(4)*uad(3,i,m))
            end do
            do i = 1, npole
               do j = 1, 3
                  fuind(j,i) = a(j,1)*ubd(1,i,m) + a(j,2)*ubd(2,i,m)
     &                            + a(j,3)*ubd(3,i,m)
                  fuinp(j,i) = a(j,1)*uap(1,i,m) + a(j,2)*uap(2,i,m)
     &                            + a(j,3)*uap(3,i,m)
               end do
            end do
            call grid_uind (fuind,fuinp)
            call fftfront
            do k = 1, nfft3
               do j = 1, nfft2
                  do i = 1, nfft1
                     term = qfac(i,j,k)
                     qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
                     qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
                  end do
               end do
            end do
            call fftback
            call fphi_uind (fphid,fphip,fphidp)
            do i = 1, npole
               do j = 2, 10
                  fphid(j,i) = f * fphid(j,i)
                  fphip(j,i) = f * fphip(j,i)
               end do
            end do
            do i = 1, npole
               ii = ipole(i)
               f1 = 0.0d0
               f2 = 0.0d0
               f3 = 0.0d0
               do k = 1, 3
                  j1 = deriv1(k+1)
                  j2 = deriv2(k+1)
                  j3 = deriv3(k+1)
                  f1 = f1+fuind(k,i)*fphip(j1,i)+fuinp(k,i)*fphid(j1,i)
                  f2 = f2+fuind(k,i)*fphip(j2,i)+fuinp(k,i)*fphid(j2,i)
                  f3 = f3+fuind(k,i)*fphip(j3,i)+fuinp(k,i)*fphid(j3,i)
               end do
               f1 = 0.5d0 * dble(nfft1) * f1
               f2 = 0.5d0 * dble(nfft2) * f2
               f3 = 0.5d0 * dble(nfft3) * f3
               h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
               h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
               h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
               dep(1,ii) = dep(1,ii) + h1
               dep(2,ii) = dep(2,ii) + h2
               dep(3,ii) = dep(3,ii) + h3
               do j = 2, 4
                  cphid(j) = 0.0d0
                  cphip(j) = 0.0d0
                  do k = 2, 4
                     cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
                     cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
                  end do
               end do
               vxx = vxx - 0.5d0*(cphid(2)*uap(1,i,m)
     &                              +cphip(2)*ubd(1,i,m))
               vxy = vxy - 0.25d0*(cphid(2)*uap(2,i,m)
     &                               +cphip(2)*ubd(2,i,m)
     &                               +cphid(3)*uap(1,i,m)
     &                               +cphip(3)*ubd(1,i,m))
               vxz = vxz - 0.25d0*(cphid(2)*uap(3,i,m)
     &                               +cphip(2)*ubd(3,i,m)
     &                               +cphid(4)*uap(1,i,m)
     &                               +cphip(4)*ubd(1,i,m))
               vyy = vyy - 0.5d0*(cphid(3)*uap(2,i,m)
     &                              +cphip(3)*ubd(2,i,m))
               vyz = vyz - 0.25d0*(cphid(3)*uap(3,i,m)
     &                               +cphip(3)*ubd(3,i,m)
     &                               +cphid(4)*uap(2,i,m)
     &                               +cphip(4)*ubd(2,i,m))
               vzz = vzz - 0.5d0*(cphid(4)*uap(3,i,m)
     &                              +cphip(4)*ubd(3,i,m))
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fphid)
      deallocate (fphip)
      deallocate (fphidp)
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
            else if (nonprism) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
            struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                  + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
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
c     assign only the induced dipoles to the PME grid
c     and perform the 3-D FFT forward transformation
c
      if (poltyp.eq.'DIRECT' .or. poltyp.eq.'TCG') then
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
               else if (nonprism) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               end if
               struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                     + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
               eterm = 0.5d0 * f * expterm * struc2
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
c     add back missing terms for the TCG polarization method;
c     first do the term for "UAD" dotted with "UBP"
c
      if (poltyp .eq. 'TCG') then
         do m = 1, tcgnab
            do i = 1, npole
               do j = 1, 10
                  cmp(j,i) = 0.0d0
               end do
               do j = 2, 4
                  cmp(j,i) = ubp(j-1,i,m)
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
                  cmp(j,i) = uad(j-1,i,m)
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
                  else if (nonprism) then
                     if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
                  end if
                  struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                        + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
                  eterm = 0.5d0 * f * expterm * struc2
                  vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
                  vxx = vxx + h1*h1*vterm - eterm
                  vxy = vxy + h1*h2*vterm
                  vxz = vxz + h1*h3*vterm
                  vyy = vyy + h2*h2*vterm - eterm
                  vyz = vyz + h2*h3*vterm
                  vzz = vzz + h3*h3*vterm - eterm
               end if
            end do
c
c     now do the TCG terms with "UBD" dotted with "UAP"
c
            do i = 1, npole
               do j = 1, 10
                  cmp(j,i) = 0.0d0
               end do
               do j = 2, 4
                  cmp(j,i) = uap(j-1,i,m)
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
                  cmp(j,i) = ubd(j-1,i,m)
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
                  else if (nonprism) then
                     if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
                  end if
                  struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                        + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
                  eterm = 0.5d0 * f * expterm * struc2
                  vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
                  vxx = vxx + h1*h1*vterm - eterm
                  vxy = vxy + h1*h2*vterm
                  vxz = vxz + h1*h3*vterm
                  vyy = vyy + h2*h2*vterm - eterm
                  vyz = vyz + h2*h3*vterm
                  vzz = vzz + h3*h3*vterm - eterm
               end if
            end do
         end do
      end if
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
         do i = 1, npole
            ii = ipole(i)
            pot(ii) = cphidp(1,i)
         end do
         call dcflux (pot,decfx,decfy,decfz)
         do i = 1, npole
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            frcx = decfx(ii)
            frcy = decfy(ii)
            frcz = decfz(ii)
            dep(1,ii) = dep(1,ii) + frcx
            dep(2,ii) = dep(2,ii) + frcy
            dep(3,ii) = dep(3,ii) + frcz
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
c
c     perform deallocation of some local arrays
c
      deallocate (cphidp)
      deallocate (qgrip)
      return
      end

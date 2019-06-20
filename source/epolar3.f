c
c
c     ##################################################
c     ##  COPYRIGHT (C) 2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved             ##
c     ##################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar3  --  induced dipole energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar3" calculates the induced dipole polarization energy,
c     and partitions the energy among atoms
c
c
      subroutine epolar3
      use limits
      implicit none
      logical pairwise
c
c
c     choose the method for summing over polarization interactions
c
      pairwise = .true.
      if (pairwise) then
         if (use_ewald) then
            if (use_mlist) then
               call epolar3d
            else
               call epolar3c
            end if
         else
            if (use_mlist) then
               call epolar3b
            else
               call epolar3a
            end if
         end if
      else
         call epolar3e
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar3a  --  double loop polarization analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar3a" calculates the induced dipole polarization energy
c     using a double loop, and partitions the energy among atoms
c
c
      subroutine epolar3a
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use energi
      use inform
      use inter
      use iounit
      use mplpot
      use molcul
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      real*8 e,f,damp,expdamp
      real*8 pdi,pti,ddi
      real*8 pgamma
      real*8 scale3,scale5
      real*8 scale7,scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r,r2
      real*8 rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dir,diu,qiu,uir
      real*8 dkr,dku,qku,ukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 dmpi(7),dmpk(7)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      logical header,huge
      character*6 mode
c
c
c     zero out the total polarization energy and partitioning
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
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
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dipole Polarization Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the dipole polarization energy component
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
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
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
               diu = dix*ukx + diy*uky + diz*ukz
               qiu = qix*ukx + qiy*uky + qiz*ukz
               uir = uix*xr + uiy*yr + uiz*zr
               dku = dkx*uix + dky*uiy + dkz*uiz
               qku = qkx*uix + qky*uiy + qkz*uiz
               ukr = ukx*xr + uky*yr + ukz*zr
c
c     find the energy value for Thole polarization damping
c
               if (use_thole) then
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     if (use_dirdamp) then
                        pgamma = min(ddi,dirdamp(kk))
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           scale3 = 1.0d0 - expdamp 
                           scale5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                         +0.15d0*damp**2)
                        end if
                     else
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                         +0.6d0*damp**2)
                        end if
                     end if
                  end if
                  scalek = pscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3 = scale3 * rr3
                  rr5 = scale5 * rr5
                  rr7 = scale7 * rr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
                  e = term1*rr3 + term2*rr5 + term3*rr7
c
c     find the energy value for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = dscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpi(3) * rr3
                  rr5i = dmpi(5) * rr5
                  rr7i = dmpi(7) * rr7
                  rr3k = dmpk(3) * rr3
                  rr5k = dmpk(5) * rr5
                  rr7k = dmpk(7) * rr7
                  e = uir*(corek*rr3+valk*rr3k)
     &                   - ukr*(corei*rr3+vali*rr3i)
     &                   + diu*rr3i + dku*rr3k
     &                   + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                   - dkr*uir*rr5k - dir*ukr*rr5i
     &                   + qkr*uir*rr7k - qir*ukr*rr7i
               end if
c
c     increment the overall polarization energy components
c
               if (e .ne. 0.0d0) then
                  ep = ep + e
                  nep = nep + 1
                  aep(i) = aep(i) + 0.5d0*e
                  aep(k) = aep(k) + 0.5d0*e
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 10.0d0)
               if ((debug.and.e.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),k,name(k),r,e
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
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
     &                  pscale(i12(j,i)) = p2iscale
                  end do
                  dscale(i12(j,i)) = pscale(i12(j,i))
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = p3scale
                  do k = 1, np11(i)
                     if (i13(j,i) .eq. ip11(k,i))
     &                  pscale(i13(j,i)) = p3iscale
                  end do
                  dscale(i13(j,i)) = pscale(i13(j,i))
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = p4scale
                  do k = 1, np11(i)
                     if (i14(j,i) .eq. ip11(k,i))
     &                  pscale(i14(j,i)) = p4iscale
                  end do
                  dscale(i14(j,i)) = pscale(i14(j,i))
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = p5scale
                  do k = 1, np11(i)
                     if (i15(j,i) .eq. ip11(k,i))
     &                  pscale(i15(j,i)) = p5iscale
                  end do
                  dscale(i15(j,i)) = pscale(i15(j,i))
               end do
            else
               do j = 1, n12(i)
                  pscale(i12(j,i)) = p2scale
                  do k = 1, np11(i)
                     if (i12(j,i) .eq. ip11(k,i))
     &                  pscale(i12(j,i)) = p2iscale
                  end do
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = p3scale
                  do k = 1, np11(i)
                     if (i13(j,i) .eq. ip11(k,i))
     &                  pscale(i13(j,i)) = p3iscale
                  end do
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = p4scale
                  do k = 1, np11(i)
                     if (i14(j,i) .eq. ip11(k,i))
     &                  pscale(i14(j,i)) = p4iscale
                  end do
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = p5scale
                  do k = 1, np11(i)
                     if (i15(j,i) .eq. ip11(k,i))
     &                  pscale(i15(j,i)) = p5iscale
                  end do
               end do
               do j = 1, np11(i)
                  dscale(ip11(j,i)) = d1scale
               end do
               do j = 1, np12(i)
                  dscale(ip12(j,i)) = d2scale
               end do
               do j = 1, np13(i)
                  dscale(ip13(j,i)) = d3scale
               end do
               do j = 1, np14(i)
                  dscale(ip14(j,i)) = d4scale
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
                     diu = dix*ukx + diy*uky + diz*ukz
                     qiu = qix*ukx + qiy*uky + qiz*ukz
                     uir = uix*xr + uiy*yr + uiz*zr
                     dku = dkx*uix + dky*uiy + dkz*uiz
                     qku = qkx*uix + qky*uiy + qkz*uiz
                     ukr = ukx*xr + uky*yr + ukz*zr
c
c     find the energy value for Thole polarization damping
c
                     if (use_thole) then
                        scale3 = 1.0d0
                        scale5 = 1.0d0
                        scale7 = 1.0d0
                        damp = pdi * pdamp(kk)
                        if (damp .ne. 0.0d0) then
                           if (use_dirdamp) then
                              pgamma = min(ddi,dirdamp(kk))
                              damp = pgamma * (r/damp)**(1.5d0)
                              if (damp .lt. 50.0d0) then
                                 expdamp = exp(-damp) 
                                 scale3 = 1.0d0 - expdamp 
                                 scale5 = 1.0d0 - expdamp
     &                                               *(1.0d0+0.5d0*damp)
                                 scale7 = 1.0d0 - expdamp
     &                                               *(1.0d0+0.65d0*damp
     &                                                  +0.15d0*damp**2)
                              end if
                           else
                              pgamma = min(pti,thole(kk))
                              damp = pgamma * (r/damp)**3
                              if (damp .lt. 50.0d0) then
                                 expdamp = exp(-damp)
                                 scale3 = 1.0d0 - expdamp
                                 scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                                 scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                               +0.6d0*damp**2)
                              end if
                           end if
                        end if
                        scalek = pscale(k)
                        rr3 = f * scalek / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr3 = scale3 * rr3
                        rr5 = scale5 * rr5
                        rr7 = scale7 * rr7
                        term1 = ck*uir - ci*ukr + diu + dku
                        term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                        term3 = uir*qkr - ukr*qir
                        e = term1*rr3 + term2*rr5 + term3*rr7
c
c     find the energy value for charge penetration damping
c
                     else if (use_chgpen) then
                        corek = pcore(kk)
                        valk = pval(kk)
                        alphak = palpha(kk)
                        call dampdir (r,alphai,alphak,dmpi,dmpk)
                        scalek = dscale(k)
                        rr3 = f * scalek / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr3i = dmpi(3) * rr3
                        rr5i = dmpi(5) * rr5
                        rr7i = dmpi(7) * rr7
                        rr3k = dmpk(3) * rr3
                        rr5k = dmpk(5) * rr5
                        rr7k = dmpk(7) * rr7
                        e = uir*(corek*rr3+valk*rr3k)
     &                         - ukr*(corei*rr3+vali*rr3i)
     &                         + diu*rr3i + dku*rr3k
     &                         + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                         - dkr*uir*rr5k - dir*ukr*rr5i
     &                         + qkr*uir*rr7k - qir*ukr*rr7i
                     end if
c
c     increment the overall polarization energy components
c
                     if (i .eq. k)  e = 0.5d0 * e
                     if (e .ne. 0.0d0) then
                        ep = ep + e
                        nep = nep + 1
                        aep(i) = aep(i) + 0.5d0*e
                        aep(k) = aep(k) + 0.5d0*e
                        einter = einter + e
                     end if
c
c     print a message if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 10.0d0)
                     if ((debug.and.e.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Polarization',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',8x,'Energy',/)
                        end if
                        write (iout,50)  i,name(i),k,name(k),r,e
   50                   format (' Polar',5x,2(i7,'-',a3),9x,
     &                             f10.4,2x,f12.4)
                     end if
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
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = 1.0d0
                  dscale(i13(j,i)) = 1.0d0
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = 1.0d0
                  dscale(i14(j,i)) = 1.0d0
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = 1.0d0
                  dscale(i15(j,i)) = 1.0d0
               end do
            else
               do j = 1, n12(i)
                  pscale(i12(j,i)) = 1.0d0
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = 1.0d0
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = 1.0d0
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = 1.0d0
               end do
               do j = 1, np11(i)
                  dscale(ip11(j,i)) = 1.0d0
               end do
               do j = 1, np12(i)
                  dscale(ip12(j,i)) = 1.0d0
               end do
               do j = 1, np13(i)
                  dscale(ip13(j,i)) = 1.0d0
               end do
               do j = 1, np14(i)
                  dscale(ip14(j,i)) = 1.0d0
               end do
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine epolar3b  --  polarization analysis via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "epolar3b" calculates the induced dipole polarization energy
c     using a neighbor list, and partitions the energy among atoms
c
c
      subroutine epolar3b
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use energi
      use inform
      use inter
      use iounit
      use molcul
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,f,damp,expdamp
      real*8 pdi,pti,ddi
      real*8 pgamma
      real*8 scale3,scale5
      real*8 scale7,scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r,r2
      real*8 rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dir,diu,qiu,uir
      real*8 dkr,dku,qku,ukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 dmpi(7),dmpk(7)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      logical header,huge
      character*6 mode
c
c
c     zero out the total polarization energy and partitioning
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
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
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dipole Polarization Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
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
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,rpole,x,y,z,pdamp,thole,dirdamp,pcore,pval,
!$OMP& palpha,uind,n12,i12,n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,
!$OMP& np13,ip13,np14,ip14,p2scale,p3scale,p4scale,p5scale,p2iscale,
!$OMP& p3iscale,p4iscale,p5iscale,d1scale,d2scale,d3scale,d4scale,
!$OMP& nelst,elst,dpequal,use_thole,use_dirdamp,use_chgpen,use_bounds,
!$OMP& off2,f,molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(pscale,dscale) shared (ep,nep,aep,einter)
!$OMP DO reduction(+:ep,nep,aep,einter) schedule(guided)
c
c     compute the dipole polarization energy component
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
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
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
               diu = dix*ukx + diy*uky + diz*ukz
               qiu = qix*ukx + qiy*uky + qiz*ukz
               uir = uix*xr + uiy*yr + uiz*zr
               dku = dkx*uix + dky*uiy + dkz*uiz
               qku = qkx*uix + qky*uiy + qkz*uiz
               ukr = ukx*xr + uky*yr + ukz*zr
c
c     find the energy value for Thole polarization damping
c
               if (use_thole) then
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     if (use_dirdamp) then
                        pgamma = min(ddi,dirdamp(kk))
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           scale3 = 1.0d0 - expdamp 
                           scale5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                         +0.15d0*damp**2)
                        end if
                     else
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                         +0.6d0*damp**2)
                        end if
                     end if
                  end if
                  scalek = pscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3 = scale3 * rr3
                  rr5 = scale5 * rr5
                  rr7 = scale7 * rr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
                  e = term1*rr3 + term2*rr5 + term3*rr7
c
c     find the energy value for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = dscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpi(3) * rr3
                  rr5i = dmpi(5) * rr5
                  rr7i = dmpi(7) * rr7
                  rr3k = dmpk(3) * rr3
                  rr5k = dmpk(5) * rr5
                  rr7k = dmpk(7) * rr7
                  e = uir*(corek*rr3+valk*rr3k)
     &                   - ukr*(corei*rr3+vali*rr3i)
     &                   + diu*rr3i + dku*rr3k
     &                   + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                   - dkr*uir*rr5k - dir*ukr*rr5i
     &                   + qkr*uir*rr7k - qir*ukr*rr7i
               end if
c
c     increment the overall polarization energy components
c
               if (e .ne. 0.0d0) then
                  ep = ep + e
                  nep = nep + 1
                  aep(i) = aep(i) + 0.5d0*e
                  aep(k) = aep(k) + 0.5d0*e
                  if (molcule(i) .ne. molcule(k))
     &               einter = einter + e
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 10.0d0)
               if ((debug.and.e.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Dipole Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),k,name(k),r,e
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
            end do
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar3c  --  Ewald polarization analysis; loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar3c" calculates the polarization energy and analysis with
c     respect to Cartesian coordinates using particle mesh Ewald and
c     a double loop
c
c
      subroutine epolar3c
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
c
c
c     zero out the dipole polarization energy and components
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
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
c     compute the real space part of the Ewald summation
c
      call epreal3c
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * term * uii / 3.0d0
         ep = ep + e
         nep = nep + 1
         aep(i) = aep(i) + e
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
         do ii = 1, npole
            i = ipole(ii)
            xd = xd + rpole(2,ii) + rpole(1,ii)*x(i)
            yd = yd + rpole(3,ii) + rpole(1,ii)*y(i)
            zd = zd + rpole(4,ii) + rpole(1,ii)*z(i)
            xu = xu + uind(1,ii)
            yu = yu + uind(2,ii)
            zu = zu + uind(3,ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
         nep = nep + 1
         do ii = 1, npole
            i = ipole(ii)
            aep(i) = aep(i) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epreal3c  --  real space polar analysis via loop  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epreal3c" calculates the induced dipole polarization energy and
c     analysis using particle mesh Ewald summation and a double loop
c
c
      subroutine epreal3c
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use energi
      use ewald
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      real*8 e,efull,f
      real*8 damp,expdamp
      real*8 erfc,bfac
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 pdi,pti,ddi
      real*8 pgamma
      real*8 scale3,scale5
      real*8 scale7,scalek
      real*8 sr3,sr5,sr7
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r,r2
      real*8 rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dir,diu,qiu,uir
      real*8 dkr,dku,qku,ukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 dmpi(7),dmpk(7)
      real*8 bn(0:3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      logical header,huge
      character*6 mode
      external erfc
c
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dipole Polarization Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the dipole polarization energy component
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
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
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
               diu = dix*ukx + diy*uky + diz*ukz
               qiu = qix*ukx + qiy*uky + qiz*ukz
               uir = uix*xr + uiy*yr + uiz*zr
               dku = dkx*uix + dky*uiy + dkz*uiz
               qku = qkx*uix + qky*uiy + qkz*uiz
               ukr = ukx*xr + uky*yr + ukz*zr
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 3
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 3
                  bn(j) = f * bn(j)
               end do
c
c     find the energy value for Thole polarization damping
c
               if (use_thole) then
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     if (use_dirdamp) then
                        pgamma = min(ddi,dirdamp(kk))
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           scale3 = 1.0d0 - expdamp 
                           scale5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                         +0.15d0*damp**2)
                        end if
                     else
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                         +0.6d0*damp**2)
                        end if
                     end if
                  end if
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  scalek = pscale(k)
                  sr3 = scalek * scale3 * rr3
                  sr5 = scalek * scale5 * rr5
                  sr7 = scalek * scale7 * rr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
                  efull = term1*sr3 + term2*sr5 + term3*sr7
                  sr3 = bn(1) - rr3 + sr3
                  sr5 = bn(2) - rr5 + sr5
                  sr7 = bn(3) - rr7 + sr7
                  e = term1*sr3 + term2*sr5 + term3*sr7
c
c     find the energy value for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = dscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpi(3) * rr3
                  rr5i = dmpi(5) * rr5
                  rr7i = dmpi(7) * rr7
                  rr3k = dmpk(3) * rr3
                  rr5k = dmpk(5) * rr5
                  rr7k = dmpk(7) * rr7
                  efull = uir*(corek*rr3+valk*rr3k)
     &                       - ukr*(corei*rr3+vali*rr3i)
     &                       + diu*rr3i + dku*rr3k
     &                       + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                       - dkr*uir*rr5k - dir*ukr*rr5i
     &                       + qkr*uir*rr7k - qir*ukr*rr7i
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = bn(1) - rr3 + rr3i
                  rr5i = bn(2) - rr5 + rr5i
                  rr7i = bn(3) - rr7 + rr7i
                  rr3k = bn(1) - rr3 + rr3k
                  rr5k = bn(2) - rr5 + rr5k
                  rr7k = bn(3) - rr7 + rr7k
                  rr3 = bn(1) - (1.0d0-scalek)*rr3
                  e = uir*(corek*rr3+valk*rr3k)
     &                   - ukr*(corei*rr3+vali*rr3i)
     &                   + diu*rr3i + dku*rr3k
     &                   + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                   - dkr*uir*rr5k - dir*ukr*rr5i
     &                   + qkr*uir*rr7k - qir*ukr*rr7i
               end if
c
c     compute the energy contribution for this interaction
c
               ep = ep + e
               if (efull .ne. 0.0d0) then
                  nep = nep + 1
                  aep(i) = aep(i) + 0.5d0*efull
                  aep(k) = aep(k) + 0.5d0*efull
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + efull
                  end if
               end if
c
c     print message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 10.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Dipole Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),k,name(k),r,efull
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
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
     &                  pscale(i12(j,i)) = p2iscale
                  end do
                  dscale(i12(j,i)) = pscale(i12(j,i))
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = p3scale
                  do k = 1, np11(i)
                     if (i13(j,i) .eq. ip11(k,i))
     &                  pscale(i13(j,i)) = p3iscale
                  end do
                  dscale(i13(j,i)) = pscale(i13(j,i))
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = p4scale
                  do k = 1, np11(i)
                     if (i14(j,i) .eq. ip11(k,i))
     &                  pscale(i14(j,i)) = p4iscale
                  end do
                  dscale(i14(j,i)) = pscale(i14(j,i))
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = p5scale
                  do k = 1, np11(i)
                     if (i15(j,i) .eq. ip11(k,i))
     &                  pscale(i15(j,i)) = p5iscale
                  end do
                  dscale(i15(j,i)) = pscale(i15(j,i))
               end do
            else
               do j = 1, n12(i)
                  pscale(i12(j,i)) = p2scale
                  do k = 1, np11(i)
                     if (i12(j,i) .eq. ip11(k,i))
     &                  pscale(i12(j,i)) = p2iscale
                  end do
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = p3scale
                  do k = 1, np11(i)
                     if (i13(j,i) .eq. ip11(k,i))
     &                  pscale(i13(j,i)) = p3iscale
                  end do
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = p4scale
                  do k = 1, np11(i)
                     if (i14(j,i) .eq. ip11(k,i))
     &                  pscale(i14(j,i)) = p4iscale
                  end do
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = p5scale
                  do k = 1, np11(i)
                     if (i15(j,i) .eq. ip11(k,i))
     &                  pscale(i15(j,i)) = p5iscale
                  end do
               end do
               do j = 1, np11(i)
                  dscale(ip11(j,i)) = d1scale
               end do
               do j = 1, np12(i)
                  dscale(ip12(j,i)) = d2scale
               end do
               do j = 1, np13(i)
                  dscale(ip13(j,i)) = d3scale
               end do
               do j = 1, np14(i)
                  dscale(ip14(j,i)) = d4scale
               end do
            end if
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
                     pscale(k) = 1.0d0
                     dscale(k) = 1.0d0
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
                     diu = dix*ukx + diy*uky + diz*ukz
                     qiu = qix*ukx + qiy*uky + qiz*ukz
                     uir = uix*xr + uiy*yr + uiz*zr
                     dku = dkx*uix + dky*uiy + dkz*uiz
                     qku = qkx*uix + qky*uiy + qkz*uiz
                     ukr = ukx*xr + uky*yr + ukz*zr
c
c     calculate the real space Ewald error function terms
c
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) / r
                     alsq2 = 2.0d0 * aewald**2
                     alsq2n = 0.0d0
                     if (aewald .gt. 0.0d0) then
                        alsq2n = 1.0d0 / (sqrtpi*aewald)
                     end if
                     exp2a = exp(-ralpha**2)
                     do j = 1, 3
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
                     do j = 0, 3
                        bn(j) = f * bn(j)
                     end do
c
c     find the energy value for Thole polarization damping
c
                     if (use_thole) then
                        scale3 = 1.0d0
                        scale5 = 1.0d0
                        scale7 = 1.0d0
                        damp = pdi * pdamp(kk)
                        if (damp .ne. 0.0d0) then
                           if (use_dirdamp) then
                              pgamma = min(ddi,dirdamp(kk))
                              damp = pgamma * (r/damp)**(1.5d0)
                              if (damp .lt. 50.0d0) then
                                 expdamp = exp(-damp) 
                                 scale3 = 1.0d0 - expdamp 
                                 scale5 = 1.0d0 - expdamp
     &                                               *(1.0d0+0.5d0*damp)
                                 scale7 = 1.0d0 - expdamp
     &                                               *(1.0d0+0.65d0*damp
     &                                                  +0.15d0*damp**2)
                              end if
                           else
                              pgamma = min(pti,thole(kk))
                              damp = pgamma * (r/damp)**3
                              if (damp .lt. 50.0d0) then
                                 expdamp = exp(damp)
                                 scale3 = 1.0d0 - expdamp
                                 scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                                 scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                               +0.6d0*damp**2)
                              end if
                           end if
                        end if
                        rr3 = f / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        scalek = pscale(k)
                        sr3 = scalek * scale3 * rr3
                        sr5 = scalek * scale5 * rr5
                        sr7 = scalek * scale7 * rr7
                        term1 = ck*uir - ci*ukr + diu + dku
                        term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                        term3 = uir*qkr - ukr*qir
                        efull = term1*sr3 + term2*sr5 + term3*sr7
                        sr3 = bn(1) - rr3 + sr3
                        sr5 = bn(2) - rr5 + sr5
                        sr7 = bn(3) - rr7 + sr7
                        e = term1*sr3 + term2*sr5 + term3*sr7
c
c     find the energy value for charge penetration damping
c
                     else if (use_chgpen) then
                        corek = pcore(kk)
                        valk = pval(kk)
                        alphak = palpha(kk)
                        call dampdir (r,alphai,alphak,dmpi,dmpk)
                        scalek = dscale(k)
                        rr3 = f * scalek / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr3i = dmpi(3) * rr3
                        rr5i = dmpi(5) * rr5
                        rr7i = dmpi(7) * rr7
                        rr3k = dmpk(3) * rr3
                        rr5k = dmpk(5) * rr5
                        rr7k = dmpk(7) * rr7
                        efull = uir*(corek*rr3+valk*rr3k)
     &                             - ukr*(corei*rr3+vali*rr3i)
     &                             + diu*rr3i + dku*rr3k
     &                             + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                             - dkr*uir*rr5k - dir*ukr*rr5i
     &                             + qkr*uir*rr7k - qir*ukr*rr7i
                        rr3 = f / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr3i = bn(1) - rr3 + rr3i
                        rr5i = bn(2) - rr5 + rr5i
                        rr7i = bn(3) - rr7 + rr7i
                        rr3k = bn(1) - rr3 + rr3k
                        rr5k = bn(2) - rr5 + rr5k
                        rr7k = bn(3) - rr7 + rr7k
                        rr3 = bn(1) - (1.0d0-scalek)*rr3
                        e = uir*(corek*rr3+valk*rr3k)
     &                         - ukr*(corei*rr3+vali*rr3i)
     &                         + diu*rr3i + dku*rr3k
     &                         + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                         - dkr*uir*rr5k - dir*ukr*rr5i
     &                         + qkr*uir*rr7k - qir*ukr*rr7i
                     end if
c
c     compute the energy contribution for this interaction
c
                     if (i .eq. k) then
                        e = 0.5d0 * e
                        efull = 0.5d0 * efull
                     end if
                     ep = ep + e
                     if (efull .ne. 0.0d0) then
                        nep = nep + 1
                        aep(i) = aep(i) + 0.5d0*efull
                        aep(k) = aep(k) + 0.5d0*efull
                        if (molcule(i) .ne. molcule(k)) then
                           einter = einter + efull
                        end if
                     end if
c
c     print message if the energy of this interaction is large
c
                     huge = (abs(efull) .gt. 10.0d0)
                     if ((debug.and.efull.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Dipole Polarization',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',8x,'Energy',/)
                        end if
                        write (iout,50)  i,name(i),k,name(k),r,efull
   50                   format (' Polar',5x,2(i7,'-',a3),9x,
     &                             f10.4,2x,f12.4)
                     end if
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
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = 1.0d0
                  dscale(i13(j,i)) = 1.0d0
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = 1.0d0
                  dscale(i14(j,i)) = 1.0d0
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = 1.0d0
                  dscale(i15(j,i)) = 1.0d0
               end do
            else
               do j = 1, n12(i)
                  pscale(i12(j,i)) = 1.0d0
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = 1.0d0
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = 1.0d0
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = 1.0d0
               end do
               do j = 1, np11(i)
                  dscale(ip11(j,i)) = 1.0d0
               end do
               do j = 1, np12(i)
                  dscale(ip12(j,i)) = 1.0d0
               end do
               do j = 1, np13(i)
                  dscale(ip13(j,i)) = 1.0d0
               end do
               do j = 1, np14(i)
                  dscale(ip14(j,i)) = 1.0d0
               end do
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar3d  --  Ewald polarization analysis; list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar3d" calculates the polarization energy and analysis with
c     respect to Cartesian coordinates using particle mesh Ewald and
c     a neighbor list
c
c
      subroutine epolar3d
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
c
c
c     zero out the dipole polarization energy and components
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
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
c     compute the real space part of the Ewald summation
c
      call epreal3d
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * term * uii / 3.0d0
         ep = ep + e
         nep = nep + 1
         aep(i) = aep(i) + e
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
         do ii = 1, npole
            i = ipole(ii)
            xd = xd + rpole(2,ii) + rpole(1,ii)*x(i)
            yd = yd + rpole(3,ii) + rpole(1,ii)*y(i)
            zd = zd + rpole(4,ii) + rpole(1,ii)*z(i)
            xu = xu + uind(1,ii)
            yu = yu + uind(2,ii)
            zu = zu + uind(3,ii)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
         nep = nep + 1
         do ii = 1, npole
            i = ipole(ii)
            aep(i) = aep(i) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epreal3d  --  real space polar analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epreal3d" calculates the induced dipole polarization energy
c     and analysis using particle mesh Ewald and a neighbor list
c
c
      subroutine epreal3d
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use energi
      use ewald
      use inform
      use inter
      use iounit
      use math
      use mplpot
      use molcul
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,efull,f
      real*8 damp,expdamp
      real*8 erfc,bfac
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 pdi,pti,ddi
      real*8 pgamma
      real*8 scale3,scale5
      real*8 scale7,scalek
      real*8 sr3,sr5,sr7
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r,r2
      real*8 rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dir,diu,qiu,uir
      real*8 dkr,dku,qku,ukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 dmpi(7),dmpk(7)
      real*8 bn(0:3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      logical header,huge
      character*6 mode
      external erfc
c
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dipole Polarization Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
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
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,rpole,uind,x,y,z,pdamp,thole,dirdamp,pcore,
!$OMP& pval,palpha,n12,i12,n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,
!$OMP& np13,ip13,np14,ip14,p2scale,p3scale,p4scale,p5scale,p2iscale,
!$OMP& p3iscale,p4iscale,p5iscale,d1scale,d2scale,d3scale,d4scale,
!$OMP& nelst,elst,dpequal,use_thole,use_dirdamp,use_chgpen,use_bounds,
!$OMP& off2,f,aewald,molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(pscale,dscale) shared (ep,nep,aep,einter)
!$OMP DO reduction(+:ep,nep,aep,einter) schedule(guided)
c
c     compute the dipole polarization energy component
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
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
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
               diu = dix*ukx + diy*uky + diz*ukz
               qiu = qix*ukx + qiy*uky + qiz*ukz
               uir = uix*xr + uiy*yr + uiz*zr
               dku = dkx*uix + dky*uiy + dkz*uiz
               qku = qkx*uix + qky*uiy + qkz*uiz
               ukr = ukx*xr + uky*yr + ukz*zr
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 3
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
               do j = 0, 3
                  bn(j) = f * bn(j)
               end do
c
c     find the energy value for Thole polarization damping
c
               if (use_thole) then
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     if (use_dirdamp) then
                        pgamma = min(ddi,dirdamp(kk))
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           scale3 = 1.0d0 - expdamp 
                           scale5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                         +0.15d0*damp**2)
                        end if
                     else
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                         +0.6d0*damp**2)
                        end if
                     end if
                  end if
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  scalek = pscale(k)
                  sr3 = scalek * scale3 * rr3
                  sr5 = scalek * scale5 * rr5
                  sr7 = scalek * scale7 * rr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
                  efull = term1*sr3 + term2*sr5 + term3*sr7
                  sr3 = bn(1) - rr3 + sr3
                  sr5 = bn(2) - rr5 + sr5
                  sr7 = bn(3) - rr7 + sr7
                  e = term1*sr3 + term2*sr5 + term3*sr7
c
c     find the energy value for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = dscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpi(3) * rr3
                  rr5i = dmpi(5) * rr5
                  rr7i = dmpi(7) * rr7
                  rr3k = dmpk(3) * rr3
                  rr5k = dmpk(5) * rr5
                  rr7k = dmpk(7) * rr7
                  efull = uir*(corek*rr3+valk*rr3k)
     &                       - ukr*(corei*rr3+vali*rr3i)
     &                       + diu*rr3i + dku*rr3k
     &                       + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                       - dkr*uir*rr5k - dir*ukr*rr5i
     &                       + qkr*uir*rr7k - qir*ukr*rr7i
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = bn(1) - rr3 + rr3i
                  rr5i = bn(2) - rr5 + rr5i
                  rr7i = bn(3) - rr7 + rr7i
                  rr3k = bn(1) - rr3 + rr3k
                  rr5k = bn(2) - rr5 + rr5k
                  rr7k = bn(3) - rr7 + rr7k
                  rr3 = bn(1) - (1.0d0-scalek)*rr3
                  e = uir*(corek*rr3+valk*rr3k)
     &                   - ukr*(corei*rr3+vali*rr3i)
     &                   + diu*rr3i + dku*rr3k
     &                   + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                   - dkr*uir*rr5k - dir*ukr*rr5i
     &                   + qkr*uir*rr7k - qir*ukr*rr7i
               end if
c
c     compute the energy contribution for this interaction
c
               ep = ep + e
               if (efull .ne. 0.0d0) then
                  nep = nep + 1
                  aep(i) = aep(i) + 0.5d0*efull
                  aep(k) = aep(k) + 0.5d0*efull
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + efull
                  end if
               end if
c
c     print message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 10.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Dipole Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),k,name(k),r,efull
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
            end do
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
            end do
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar3e  --  single-loop polarization analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar3e" calculates the dipole polarizability interaction
c     from the induced dipoles times the electric field
c
c
      subroutine epolar3e
      use action
      use analyz
      use atoms
      use atomid
      use boxes
      use chgpot
      use energi
      use ewald
      use inform
      use iounit
      use limits
      use math
      use mpole
      use polar
      use polpot
      use potent
      use units
      implicit none
      integer i,j,ii
      real*8 e,f,fi,term
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
      logical header,huge
c
c
c     zero out the total polarization energy and partitioning
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
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
c     print header information if debug output was requested
c
      header = .true.
      if (debug) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dipole Polarization Interactions :',
     &           //,' Type',9x,'Atom Name',24x,'Alpha',8x,'Energy',/)
      end if
c
c     set the energy unit conversion factor
c
      f = -0.5d0 * electric / dielec
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private (i,j,fi,e,ii,huge)
!$OMP DO reduction(+:ep,nep,aep) schedule(guided)
c
c     get polarization energy via induced dipoles times field
c
      do ii = 1, npole
         if (douind(ipole(ii))) then
            i = ipole(ii)
            fi = f / polarity(ii)
            e = 0.0d0
            do j = 1, 3
               e = e + fi*uind(j,ii)*udirp(j,ii)
            end do
            nep = nep + 1
            ep = ep + e
            aep(i) = aep(i) + e
c
c     print a message if the energy for this site is large
c
            huge = (abs(e) .gt. 10.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Individual Dipole Polarization',
     &                       ' Interactions :',
     &                    //,' Type',9x,'Atom Name',24x,'Alpha',
     &                       8x,'Energy',/)
               end if
               write (iout,30)  i,name(i),polarity(ii),e
   30          format (' Polar',5x,i7,'-',a3,16x,2f14.4)
            end if
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
            nep = nep + 1
            ep = ep + e
            do ii = 1, npole
               i = ipole(ii)
               aep(i) = aep(i) + e/dble(npole)
            end do
         end if
      end if
      return
      end

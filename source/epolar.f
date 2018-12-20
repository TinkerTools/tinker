c
c
c     ##################################################
c     ##  COPYRIGHT (C) 2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved             ##
c     ##################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epolar  --  induced dipole polarization energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epolar" calculates the polarization energy due to induced
c     dipole interactions
c
c
      subroutine epolar
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
               call epolar0d
            else
               call epolar0c
            end if
         else
            if (use_mlist) then
               call epolar0b
            else
               call epolar0a
            end if
         end if
      else
         call epolar0e
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar0a  --  double loop polarization energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar0a" calculates the induced dipole polarization energy
c     using a double loop, and partitions the energy among atoms
c
c
      subroutine epolar0a
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use energi
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 e,f,damp,expdamp
      real*8 pdi,pti,pgamma
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
      character*6 mode
c
c
c     zero out the total induced dipole polarization energy
c
      ep = 0.0d0
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
         pdi = pdamp(ii)
         pti = thole(ii)
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
                if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
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
                     pgamma = min(pti,thole(kk))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
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
               ep = ep + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
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
            pdi = pdamp(ii)
            pti = thole(ii)
            if (use_chgpen) then
               corei = pcore(ii)
               vali = pval(ii)
               alphai = palpha(ii)
            end if
c
c     set exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4scale * p41scale
               end do
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
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
                  if (.not. (use_polymer .and. r2.le.polycut2))
     &               pscale(k) = 1.0d0
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
                           pgamma = min(pti,thole(kk))
                           damp = -pgamma * (r/damp)**3
                           if (damp .gt. -50.0d0) then
                              expdamp = exp(damp)
                              scale3 = 1.0d0 - expdamp
                              scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                              scale7 = 1.0d0 - (1.0d0-damp
     &                                    +0.6d0*damp**2)*expdamp
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
                     ep = ep + e
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
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
c     ##  subroutine epolar0b  --  neighbor list polarization energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar0b" calculates the induced dipole polarization energy
c     using a neighbor list
c
c
      subroutine epolar0b
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use energi
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
      real*8 pdi,pti,pgamma
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
      character*6 mode
c
c
c     zero out the total polarization energy and partitioning
c
      ep = 0.0d0
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
!$OMP& shared(npole,ipole,rpole,x,y,z,pdamp,thole,pcore,pval,palpha,
!$OMP& uind,n12,i12,n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,
!$OMP& np13,ip13,np14,ip14,p2scale,p3scale,p4scale,p5scale,p41scale,
!$OMP& d1scale,d2scale,d3scale,d4scale,nelst,elst,use_thole,use_chgpen,
!$OMP& use_bounds,f,off2)
!$OMP& firstprivate(pscale,dscale) shared (ep)
!$OMP DO reduction(+:ep) schedule(guided)
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
         pdi = pdamp(ii)
         pti = thole(ii)
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
                if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
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
                     pgamma = min(pti,thole(kk))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
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
               ep = ep + e
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
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
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar0c  --  Ewald polarization derivs via loop  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar0c" calculates the dipole polarization energy with respect
c     to Cartesian coordinates using particle mesh Ewald summation and
c     a double loop
c
c
      subroutine epolar0c
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
      real*8 e,f,term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
c
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0d0
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
      call epreal0c
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip
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
         ep = ep + e
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
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal0c  --  real space polar energy via loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal0c" calculates the induced dipole polarization energy
c     using particle mesh Ewald summation and a double loop
c
c
      subroutine epreal0c
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use energi
      use ewald
      use math
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
      real*8 e,f,damp
      real*8 expdamp
      real*8 erfc,bfac
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 pdi,pti,pgamma
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
      character*6 mode
      external erfc
c
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
         pdi = pdamp(ii)
         pti = thole(ii)
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
                if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
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
                     pgamma = min(pti,thole(kk))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     end if
                  end if
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  scalek = pscale(k)
                  sr3 = scalek * scale3 * rr3
                  sr5 = scalek * scale5 * rr5
                  sr7 = scalek * scale7 * rr7
                  sr3 = bn(1) - rr3 + sr3
                  sr5 = bn(2) - rr5 + sr5
                  sr7 = bn(3) - rr7 + sr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
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
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
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
            pdi = pdamp(ii)
            pti = thole(ii)
            if (use_chgpen) then
               corei = pcore(ii)
               vali = pval(ii)
               alphai = palpha(ii)
            end if
c
c     set exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4scale * p41scale
               end do
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
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
                           pgamma = min(pti,thole(kk))
                           damp = -pgamma * (r/damp)**3
                           if (damp .gt. -50.0d0) then
                              expdamp = exp(damp)
                              scale3 = 1.0d0 - expdamp
                              scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                              scale7 = 1.0d0 - (1.0d0-damp
     &                                    +0.6d0*damp**2)*expdamp
                           end if
                        end if
                        scalek = pscale(k)
                        rr3 = f / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        sr3 = scalek * scale3 * rr3
                        sr5 = scalek * scale5 * rr5
                        sr7 = scalek * scale7 * rr7
                        sr3 = bn(1) - rr3 + sr3
                        sr5 = bn(2) - rr5 + sr5
                        sr7 = bn(3) - rr7 + sr7
                        term1 = ck*uir - ci*ukr + diu + dku
                        term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                        term3 = uir*qkr - ukr*qir
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
                     if (i .eq. k)  e = 0.5d0 * e
                     ep = ep + e
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
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
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar0d  --  Ewald polarization derivs via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar0d" calculates the dipole polarization energy with respect
c     to Cartesian coordinates using particle mesh Ewald summation and
c     a neighbor list
c
c
      subroutine epolar0d
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
      real*8 e,f,term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
c
c
c     zero out the polarization energy and derivatives
c
      ep = 0.0d0
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
      call epreal0d
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip
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
         ep = ep + e
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
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epreal0d  --  real space polar energy via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epreal0d" calculates the induced dipole polarization energy
c     using particle mesh Ewald summation and a neighbor list
c
c
      subroutine epreal0d
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use energi
      use ewald
      use math
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
      real*8 e,f,damp
      real*8 expdamp
      real*8 erfc,bfac
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 pdi,pti,pgamma
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
      character*6 mode
      external erfc
c
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
!$OMP& shared(npole,ipole,rpole,uind,x,y,z,pdamp,thole,pcore,pval,
!$OMP& palpha,n12,i12,n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,
!$OMP& np13,ip13,np14,ip14,p2scale,p3scale,p4scale,p5scale,p41scale,
!$OMP& d1scale,d2scale,d3scale,d4scale,nelst,elst,use_thole,use_chgpen,
!$OMP& use_bounds,off2,f,aewald)
!$OMP& firstprivate(pscale,dscale) shared (ep)
!$OMP DO reduction(+:ep) schedule(guided)
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
         pdi = pdamp(ii)
         pti = thole(ii)
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
                if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
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
                     pgamma = min(pti,thole(kk))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     end if
                  end if
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  scalek = pscale(k)
                  sr3 = scalek * scale3 * rr3
                  sr5 = scalek * scale5 * rr5
                  sr7 = scalek * scale7 * rr7
                  sr3 = bn(1) - rr3 + sr3
                  sr5 = bn(2) - rr5 + sr5
                  sr7 = bn(3) - rr7 + sr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
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
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar0e  --  single-loop polarization energy  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar0e" calculates the dipole polarizability interaction
c     from the induced dipoles times the electric field
c
c
      subroutine epolar0e
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
      use potent
      implicit none
      integer i,j,ii
      real*8 e,f,fi,term
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
c
c
c     zero out the total polarization energy
c
      ep = 0.0d0
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
c     ##  subroutine eprecip  --  PME recip space polarization energy  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "eprecip" evaluates the reciprocal space portion of particle
c     mesh Ewald summation energy due to dipole polarization
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
      subroutine eprecip
      use atoms
      use bound
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use mrecip
      use pme
      use polar
      use polpot
      use potent
      implicit none
      integer i,j,k
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      real*8 e,r1,r2,r3
      real*8 f,h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 struc2
      real*8 a(3,3),ftc(10,10)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
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
      if (allocated(fphi) .and. size(fphi).lt.20*npole)
     &   deallocate (fphi)
      if (.not. allocated(cmp))  allocate (cmp(10,npole))
      if (.not. allocated(fmp))  allocate (fmp(10,npole))
      if (.not. allocated(fphi))  allocate (fphi(20,npole))
c
c     get the fractional to Cartesian transformation matrix
c
      call frac_to_cart (ftc)
c
c     perform dynamic allocation of some global arrays
c
      if (.not. use_mpole) then
         ntot = nfft1 * nfft2 * nfft3
         if (allocated(qgrid) .and. size(qgrid).ne.2*ntot)
     &      deallocate(qgrid)
         if (allocated(qfac) .and. size(qfac).ne.ntot)
     &      deallocate(qfac)
         if (.not. allocated(qgrid))
     &      allocate (qgrid(2,nfft1,nfft2,nfft3))
         if (.not. allocated(qfac))
     &      allocate (qfac(nfft1,nfft2,nfft3))
c
c     setup spatial decomposition, B-splines and PME arrays
c
         call getchunk
         call moduli
         call fftsetup

c     compute B-spline coefficients and spatial decomposition
c
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
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole))
      allocate (fuinp(3,npole))
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
c     account for zeroth grid point for nonperiodic system
c
      if (.not. use_bounds) then
         call grid_uind (fuind,fuinp)
         call fftfront
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = f * expterm * struc2
         ep = ep + e
      end if
c
c     increment the induced dipole polarization energy
c
      e = 0.0d0
      do i = 1, npole
         do k = 1, 3
            e = e + fuind(k,i)*fphi(k+1,i)
         end do
      end do
      e = 0.5d0 * f * e
      ep = ep + e
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      return
      end

c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2022 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dexpol  --  variable polarizability chain rule  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dexpol" computes the chain rule terms to account for variable
c     polarizability in exchange polarization
c
c     literature reference:
c
c     M. K. J. Chung, Z. Wang, J. A. Rackers and J. W. Ponder,
c     "Classical Exchange Polarization: An Anisotropic Variable
c     Polarizability Model", Journal of Physical Chemistry B,
c     submitted, June 2022
c
c
      subroutine dexpol
      use limits
      implicit none
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_mlist) then
         call dexpol1b
      else
         call dexpol1a
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dexpol1a  --  exch-polar chain rule via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dexpol1a" finds variable polarizability chain rule gradient
c     components due to exchange polarization using a double loop
c
c
      subroutine dexpol1a
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use deriv
      use expol
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use units
      use virial
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
      real*8 sizi,sizk,sizik
      real*8 alphai,alphak
      real*8 springi,springk
      real*8 f,s2,ds2
      real*8 s2i,s2k
      real*8 ds2i,ds2k
      real*8 taper,dtaper
      real*8 uixl,ukxl
      real*8 uiyl,ukyl
      real*8 uizl,ukzl
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 frcil(3)
      real*8 frckl(3)
      real*8 frcxi,frcyi,frczi
      real*8 frcxk,frcyk,frczk
      real*8 frcx,frcy,frcz
      real*8 tqxil,tqyil
      real*8 tqxkl,tqykl
      real*8 ai(3,3)
      real*8 ak(3,3)
      real*8, allocatable :: pscale(:)
      logical epli,eplk,do_g
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'REPULS'
      call switch (mode)
c
c     do_gradient set to true
c
      do_g = .true.
c
c     set array needed to scale atom and group interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     find the exchange polarization gradient
c
      do ii = 1, npole-1
         i = ipole(ii)
         springi = kpep(ii) / polarity(ii)
         sizi = prepep(ii)
         alphai = dmppep(ii)
         epli = lpep(ii)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &            pscale(i12(j,i)) = p2iscale
            end do
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &            pscale(i13(j,i)) = p3iscale
            end do
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4iscale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &            pscale(i15(j,i)) = p5iscale
            end do
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            eplk = lpep(kk)
            if (epli .or. eplk) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  springk = kpep(kk)/polarity(kk)
                  sizk = prepep(kk)
                  alphak = dmppep(kk)
                  sizik = sizi * sizk
                  call dampexpl (r,sizik,alphai,alphak,s2,ds2,do_g)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     ds2 = ds2*taper + s2*dtaper
                     s2 = s2 * taper
                  end if
                  s2i = springi * s2 * pscale(k)
                  s2k = springk * s2 * pscale(k)
                  ds2i = springi * ds2 * pscale(k)
                  ds2k = springk * ds2 * pscale(k)
                  call rotdexpl (r,xr,yr,zr,ai,ak)
                  uixl = 0.0d0
                  ukxl = 0.0d0
                  uiyl = 0.0d0
                  ukyl = 0.0d0
                  uizl = 0.0d0
                  ukzl = 0.0d0
                  do j = 1, 3
                     uixl = uixl + uind(j,ii)*ai(1,j)
                     ukxl = ukxl - uind(j,kk)*ak(1,j)
                     uiyl = uiyl + uind(j,ii)*ai(2,j)
                     ukyl = ukyl - uind(j,kk)*ak(2,j)
                     uizl = uizl + uind(j,ii)*ai(3,j)
                     ukzl = ukzl - uind(j,kk)*ak(3,j)
                  end do
                  frcil(3) = uizl**2 * ds2i
                  frckl(3) = ukzl**2 * ds2k
c
c     compute torque in local frame
c
                  tqxil = 2.0d0 * uiyl * uizl * s2i
                  tqyil = -2.0d0 * uixl * uizl * s2i
                  tqxkl = 2.0d0 * ukyl * ukzl * s2k
                  tqykl = -2.0d0 * ukxl * ukzl * s2k
c
c     convert torque to forces
c
                  frcil(1) = -tqyil / r
                  frcil(2) = tqxil / r
                  frckl(1) = -tqykl / r
                  frckl(2) = tqxkl / r
c
c     rotate force to global frame
c
                  frcxi = 0.0d0
                  frcyi = 0.0d0
                  frczi = 0.0d0
                  frcxk = 0.0d0
                  frcyk = 0.0d0
                  frczk = 0.0d0
                  do j = 1, 3
                     frcxi = frcxi + ai(j,1)*frcil(j)
                     frcyi = frcyi + ai(j,2)*frcil(j)
                     frczi = frczi + ai(j,3)*frcil(j)
                     frcxk = frcxk + ak(j,1)*frckl(j)
                     frcyk = frcyk + ak(j,2)*frckl(j)
                     frczk = frczk + ak(j,3)*frckl(j)
                  end do
                  frcx = f * (frcxk-frcxi)
                  frcy = f * (frcyk-frcyi)
                  frcz = f * (frczk-frczi)
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
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction components with other unit cells
c
         do ii = 1, npole
            i = ipole(ii)
            springi = kpep(ii) / polarity(ii)
            sizi = prepep(ii)
            alphai = dmppep(ii)
            epli = lpep(ii)
c
c     set exclusion coefficients for connected atoms
c
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
c
c     evaluate all sites within the cutoff distance
c
            do kk = ii, npole
               k = ipole(kk)
               eplk = lpep(kk)
               if (epli .or. eplk) then
                  do jcell = 2, ncell
                     xr = x(k) - x(i)
                     yr = y(k) - y(i)
                     zr = z(k) - z(i)
                     call imager (xr,yr,zr,jcell)
                     r2 = xr*xr + yr*yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        springk = kpep(kk) / polarity(kk)
                        sizk = prepep(kk)
                        alphak = dmppep(kk)
                        sizik = sizi * sizk
                        call dampexpl (r,sizik,alphai,alphak,s2,ds2,
     &                                    do_g)
c
c     use energy switching if near the cutoff distance
c
                        if (r2 .gt. cut2) then
                           r3 = r2 * r
                           r4 = r2 * r2
                           r5 = r2 * r3
                           taper = c5*r5 + c4*r4 + c3*r3
     &                                + c2*r2 + c1*r + c0
                           dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                                 + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                           ds2 = ds2*taper + s2*dtaper
                           s2 = s2 * taper
                        end if
c
c     interaction of an atom with its own image counts half
c
                        if (i .eq. k) then
                           s2 = 0.5d0 * s2
                           ds2 = 0.5d0 * ds2
                        end if
                        s2i = springi * s2 * pscale(k)
                        s2k = springk * s2 * pscale(k)
                        ds2i = springi * ds2 * pscale(k)
                        ds2k = springk * ds2 * pscale(k)
                        call rotdexpl (r,xr,yr,zr,ai,ak)
                        uixl = 0.0d0
                        ukxl = 0.0d0
                        uiyl = 0.0d0
                        ukyl = 0.0d0
                        uizl = 0.0d0
                        ukzl = 0.0d0
                        do j = 1, 3
                           uixl = uixl + uind(j,ii)*ai(1,j)
                           ukxl = ukxl - uind(j,kk)*ak(1,j)
                           uiyl = uiyl + uind(j,ii)*ai(2,j)
                           ukyl = ukyl - uind(j,kk)*ak(2,j)
                           uizl = uizl + uind(j,ii)*ai(3,j)
                           ukzl = ukzl - uind(j,kk)*ak(3,j)
                        end do
                        frcil(3) = uizl**2 * ds2i
                        frckl(3) = ukzl**2 * ds2k
c
c     compute torque in local frame
c
                        tqxil = 2.0d0 * uiyl * uizl * s2i
                        tqyil = -2.0d0 * uixl * uizl * s2i
                        tqxkl = 2.0d0 * ukyl * ukzl * s2k
                        tqykl = -2.0d0 * ukxl * ukzl * s2k
c
c     convert torque to forces
c
                        frcil(1) = -tqyil / r
                        frcil(2) = tqxil / r
                        frckl(1) = -tqykl / r
                        frckl(2) = tqxkl / r
c
c     rotate force to global frame
c
                        frcxi = 0.0d0
                        frcyi = 0.0d0
                        frczi = 0.0d0
                        frcxk = 0.0d0
                        frcyk = 0.0d0
                        frczk = 0.0d0
                        do j = 1, 3
                           frcxi = frcxi + ai(j,1)*frcil(j)
                           frcyi = frcyi + ai(j,2)*frcil(j)
                           frczi = frczi + ai(j,3)*frcil(j)
                           frcxk = frcxk + ak(j,1)*frckl(j)
                           frcyk = frcyk + ak(j,2)*frckl(j)
                           frczk = frczk + ak(j,3)*frckl(j)
                        end do
                        frcx = f * (frcxk-frcxi)
                        frcy = f * (frcyk-frcyi)
                        frcz = f * (frczk-frczi)
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
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dexpol1b  --  exch-polar chain rule via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dexpol1b" finds variable polarizability chain rule gradient
c     components due to exchange polarization using a neighbor list
c
c
      subroutine dexpol1b
      use atoms
      use bound
      use chgpot
      use couple
      use deriv
      use expol
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use units
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
      real*8 sizi,sizk,sizik
      real*8 f,alphai,alphak
      real*8 springi,springk
      real*8 s2,ds2
      real*8 s2i, s2k
      real*8 ds2i, ds2k
      real*8 taper,dtaper
      real*8 uixl,ukxl
      real*8 uiyl,ukyl
      real*8 uizl,ukzl
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 frcil(3)
      real*8 frckl(3)
      real*8 frcxi,frcyi,frczi
      real*8 frcxk,frcyk,frczk
      real*8 frcx,frcy,frcz
      real*8 tqxil,tqyil
      real*8 tqxkl,tqykl
      real*8 ai(3,3)
      real*8 ak(3,3)
      real*8, allocatable :: pscale(:)
      logical epli,eplk,do_g
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'REPULS'
      call switch (mode)
c
c     do_gradient set to true
c
      do_g = .true.
c
c     set array needed to scale atom and group interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,kpep,prepep,dmppep,lpep,np11,ip11,n12,
!$OMP& i12,n13,i13,n14,i14,n15,i15,p2scale,p3scale,p4scale,p5scale,
!$OMP& p2iscale,p3iscale,p4iscale,p5iscale,nelst,elst,use_bounds,
!$OMP& cut2,off2,c0,c1,c2,c3,c4,c5,polarity,f,uind,do_g)
!$OMP& firstprivate(pscale)
!$OMP& shared (dep,vir)
!$OMP DO reduction(+:dep,vir) schedule(guided)
c
c     find the exchange polarization gradient
c
      do ii = 1, npole
         i = ipole(ii)
         springi = kpep(ii) / polarity(ii)
         sizi = prepep(ii)
         alphai = dmppep(ii)
         epli = lpep(ii)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &            pscale(i12(j,i)) = p2iscale
            end do
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &            pscale(i13(j,i)) = p3iscale
            end do
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4iscale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &            pscale(i15(j,i)) = p5iscale
            end do
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            eplk = lpep(kk)
            if (epli .or. eplk) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  springk = kpep(kk) / polarity(kk)
                  sizk = prepep(kk)
                  alphak = dmppep(kk)
                  sizik = sizi * sizk
                  call dampexpl (r,sizik,alphai,alphak,s2,ds2,do_g)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     ds2 = ds2*taper + s2*dtaper
                     s2 = s2 * taper
                  end if
                  s2i = springi * s2 * pscale(k)
                  s2k = springk * s2 * pscale(k)
                  ds2i = springi * ds2 * pscale(k)
                  ds2k = springk * ds2 * pscale(k)
                  call rotdexpl (r,xr,yr,zr,ai,ak)
                  uixl = 0.0d0
                  ukxl = 0.0d0
                  uiyl = 0.0d0
                  ukyl = 0.0d0
                  uizl = 0.0d0
                  ukzl = 0.0d0
                  do j = 1, 3
                     uixl = uixl + uind(j,ii)*ai(1,j)
                     ukxl = ukxl - uind(j,kk)*ak(1,j)
                     uiyl = uiyl + uind(j,ii)*ai(2,j)
                     ukyl = ukyl - uind(j,kk)*ak(2,j)
                     uizl = uizl + uind(j,ii)*ai(3,j)
                     ukzl = ukzl - uind(j,kk)*ak(3,j)
                  end do
                  frcil(3) = uizl**2 * ds2i
                  frckl(3) = ukzl**2 * ds2k
c
c     compute torque in local frame
c
                  tqxil = 2.0d0 * uiyl * uizl * s2i
                  tqyil = -2.0d0 * uixl * uizl * s2i
                  tqxkl = 2.0d0 * ukyl * ukzl * s2k
                  tqykl = -2.0d0 * ukxl * ukzl * s2k
c
c     convert torque to forces
c
                  frcil(1) = -tqyil / r
                  frcil(2) = tqxil / r
                  frckl(1) = -tqykl / r
                  frckl(2) = tqxkl / r
c
c     rotate force to global frame
c
                  frcxi = 0.0d0
                  frcyi = 0.0d0
                  frczi = 0.0d0
                  frcxk = 0.0d0
                  frcyk = 0.0d0
                  frczk = 0.0d0
                  do j = 1, 3
                     frcxi = frcxi + ai(j,1)*frcil(j)
                     frcyi = frcyi + ai(j,2)*frcil(j)
                     frczi = frczi + ai(j,3)*frcil(j)
                     frcxk = frcxk + ak(j,1)*frckl(j)
                     frcyk = frcyk + ak(j,2)*frckl(j)
                     frczk = frczk + ak(j,3)*frckl(j)
                  end do
                  frcx = f * (frcxk-frcxi)
                  frcy = f * (frcyk-frcyi)
                  frcz = f * (frczk-frczi)
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
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotdexpl  --  rotation matrix for expol derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotdexpl" finds rotation matrices for variable polarizability
c     in the exchange polarization gradient
c
c
      subroutine rotdexpl (r,xr,yr,zr,ai,ak)
      use atoms
      use math
      use mpole
      use polpot
      implicit none
      integer i,j
      real*8 xr,yr,zr
      real*8 r,dot,eps
      real*8 dx,dy,dz
      real*8 dr
      real*8 ai(3,3)
      real*8 ak(3,3)
c
c
c     compute the rotation matrix elements
c
      ai(3,1) = xr / r
      ai(3,2) = yr / r
      ai(3,3) = zr / r
      dx = 1.0d0
      dy = 0.0d0
      dz = 0.0d0
      dot = ai(3,1)
      eps = 1.0d0 / root2
      if (abs(dot) .gt. eps) then
         dx = 0.0d0
         dy = 1.0d0
         dot = ai(3,2)
      end if
      dx = dx - dot*ai(3,1)
      dy = dy - dot*ai(3,2)
      dz = dz - dot*ai(3,3)
      dr = sqrt(dx*dx + dy*dy + dz*dz)
c
c     matrix "ai" rotates a vector from global to local frame
c
      ai(1,1) = dx / dr
      ai(1,2) = dy / dr
      ai(1,3) = dz / dr
      ai(2,1) = ai(1,3)*ai(3,2) - ai(1,2)*ai(3,3)
      ai(2,2) = ai(1,1)*ai(3,3) - ai(1,3)*ai(3,1)
      ai(2,3) = ai(1,2)*ai(3,1) - ai(1,1)*ai(3,2)
      ak(1,1) = ai(1,1)
      ak(1,2) = ai(1,2)
      ak(1,3) = ai(1,3)
      ak(2,1) = -ai(2,1)
      ak(2,2) = -ai(2,2)
      ak(2,3) = -ai(2,3)
      ak(3,1) = -ai(3,1)
      ak(3,2) = -ai(3,2)
      ak(3,3) = -ai(3,3)
      return
      end

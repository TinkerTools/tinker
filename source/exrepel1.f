c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine exrepel1  --  exch repulsion energy & derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "exrepel1" calculates the exchange repulsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c     literature reference:
c
c     TBD
c
c
      subroutine exrepel1
      use limits
      implicit none
c
c
c     choose the method for summing over pairwise interactions
c
c      if (use_mlist) then
c         call exrepel1b
c      else
c         call exrepel1a
c      end if
      call exrepel1a
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine exrepel1a  --  exch repulsion analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "exrepel1a" calculates the exchange repulsion energy and first
c     derivatives with respect to Cartesian coordinates using a
c     pairwise double loop
c
c
      subroutine exrepel1a
      use atoms
      use bound
      use cell
      use couple
      use deriv
      use energi
      use group
      use mpole
      use mutant
      use repel
      use reppot
      use shunt
      use units
      use usage
      use virial
      use xrepel
      implicit none
      integer i,j,k
      integer ii,kk
      integer ind1,ind2,ind3
      real*8 e,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3
      real*8 normi
      real*8 zxri,zxrk
      real*8 vali,valk
      real*8 dmpi,dmpk
      real*8 dis,dks
      real*8 dmpip,dmpkp
      real*8 cis,cks
      real*8 cix,ckx
      real*8 ciy,cky
      real*8 ciz,ckz
      real*8 rcix,rckx
      real*8 rciy,rcky
      real*8 rciz,rckz
      real*8 cscs,cxcx,cycy
      real*8 czcz,cscz,czcs
      real*8 SS,SPz,PzS
      real*8 PxPx,PyPy,PzPz
      real*8 dSS,dSPz,dPzS
      real*8 dPxPx,dPyPy,dPzPz
      real*8 intS,intS2
      real*8 dintS
      real*8 dintSx,dintSy,dintSz
      real*8 intSR,preintSR
      real*8 pre
      real*8 term1
      real*8 term2x,term2y,term2z
      real*8 frcx,frcy,frcz
      real*8 ncix,nciy,nciz
      real*8 nckx,ncky,nckz
      real*8 nrcix,nrciy,nrciz
      real*8 nrckx,nrcky,nrckz
      real*8 drcixdx,drcixdy,drcixdz
      real*8 drciydx,drciydy,drciydz
      real*8 drcizdx,drcizdy,drcizdz
      real*8 drckxdx,drckxdy,drckxdz
      real*8 drckydx,drckydy,drckydz
      real*8 drckzdx,drckzdy,drckzdz
      real*8 tixintS,tiyintS,tizintS
      real*8 tkxintS,tkyintS,tkzintS
      real*8 bi(3)
      real*8 bj(3)
      real*8 bk(3)
      real*8 ttri(3),ttrk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8, allocatable :: rscale(:)
      real*8, allocatable :: ter(:,:)
      logical proceed,usei
      logical muti,mutk,mutik
      logical header,huge
      logical grad
      character*6 mode
c
c
c     zero out the repulsion energy and derivatives
c
      er = 0.0d0
      do i = 1, n
         der(1,i) = 0.0d0
         der(2,i) = 0.0d0
         der(3,i) = 0.0d0
      end do
      if (nrep .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     determine pseudo orbital coefficients
c
      call solvcoeff
c
c     rotate the coefficient components into the global frame
c
      call rotcoeff
c
c     perform dynamic allocation of some local arrays
c
      allocate (rscale(n))
      allocate (ter(3,n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         rscale(i) = 1.0d0
         do j = 1, 3
            ter(j,i) = 0.0d0
         end do
      end do
c
c     set cutoff distances and switching coefficients
c
      mode = 'REPULS'
      call switch (mode)
c
c     set gradient mode to true
c
      grad = .true.
c
c     calculate the exchange repulsion energy and derivative
c
      do ii = 1, nrep-1
         i = irep(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         zxri = zpxr(i)
         dmpi = dmppxr(i)
         vali = zxri
         dis = 1.0d0
         dmpip = dis * dmpi
         cis = rcpxr(1,i)
         cix = rcpxr(2,i)
         ciy = rcpxr(3,i)
         ciz = rcpxr(4,i)
         usei = use(i)
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            rscale(i12(j,i)) = r2scale
         end do
         do j = 1, n13(i)
            rscale(i13(j,i)) = r3scale
         end do
         do j = 1, n14(i)
            rscale(i14(j,i)) = r4scale
         end do
         do j = 1, n15(i)
            rscale(i15(j,i)) = r5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, nrep
            k = irep(kk)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr * xr + yr * yr + zr * zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  r3 = r2 * r
                  zxrk = zpxr(k)
                  dmpk = dmppxr(k)
                  valk = zxrk
                  dks = 1.0d0
                  dmpkp = dks * dmpk
                  cks = rcpxr(1,k)
                  ckx = rcpxr(2,k)
                  cky = rcpxr(3,k)
                  ckz = rcpxr(4,k)
c
c     choose orthogonal 2-body coordinates / solve rotation matrix
c
                  bk(1) = xr / r
                  bk(2) = yr / r
                  bk(3) = zr / r
                  ind1 = maxloc(abs(bk), dim=1)
                  ind2 = mod(ind1,3) + 1
                  ind3 = mod(ind1+1,3) + 1
                  bi(ind1) = -bk(ind2)
                  bi(ind2) = bk(ind1)
                  bi(ind3) = 0.0d0
                  normi = sqrt(bi(1)**2 + bi(2)**2 + bi(3)**2)
                  bi(1) = bi(1) / normi
                  bi(2) = bi(2) / normi
                  bi(3) = bi(3) / normi
                  bj(1) = bk(2)*bi(3) - bk(3)*bi(2)
                  bj(2) = bk(3)*bi(1) - bk(1)*bi(3)
                  bj(3) = bk(1)*bi(2) - bk(2)*bi(1)
c
c     rotate p orbital cofficients to 2-body (prolate spheroid) frame
c
                  rcix = bi(1)*cix + bi(2)*ciy + bi(3)*ciz
                  rciy = bj(1)*cix + bj(2)*ciy + bj(3)*ciz
                  rciz = bk(1)*cix + bk(2)*ciy + bk(3)*ciz
                  rckx = bi(1)*ckx + bi(2)*cky + bi(3)*ckz
                  rcky = bj(1)*ckx + bj(2)*cky + bj(3)*ckz
                  rckz = bk(1)*ckx + bk(2)*cky + bk(3)*ckz
                  cscs = cis * cks
                  cxcx = rcix * rckx
                  cycy = rciy * rcky
                  czcz = rciz * rckz
                  cscz = cis * rckz
                  czcs = rciz * cks
                  call computeOverlap (dmpi, dmpk, dmpip, dmpkp, 0.0d0,
     &                           r, grad, SS, dSS, SPz, dSPz, PzS, dPzS,
     &                           PxPx, dPxPx, PyPy, dPyPy, PzPz, dPzPz)
                  intS = cscs * SS + cxcx * PxPx + cycy * PyPy
     &                           + czcz * PzPz + cscz * SPz + czcs * PzS
                  intS2 = intS * intS
                  dintS = cscs * dSS + cxcx * dPxPx + cycy * dPyPy
     &                        + czcz * dPzPz + cscz * dSPz + czcs * dPzS
                  dintSx = dintS * bk(1)
                  dintSy = dintS * bk(2)
                  dintSz = dintS * bk(3)
                  drcixdx = bi(1)*(-rciz/r)
                  drcixdy = bi(2)*(-rciz/r)
                  drcixdz = bi(3)*(-rciz/r)
                  drciydx = bj(1)*(-rciz/r)
                  drciydy = bj(2)*(-rciz/r)
                  drciydz = bj(3)*(-rciz/r)
                  drcizdx = bi(1)*( rcix/r) + bj(1)*( rciy/r)
                  drcizdy = bi(2)*( rcix/r) + bj(2)*( rciy/r)
                  drcizdz = bi(3)*( rcix/r) + bj(3)*( rciy/r)
                  drckxdx = bi(1)*(-rckz/r)
                  drckxdy = bi(2)*(-rckz/r)
                  drckxdz = bi(3)*(-rckz/r)
                  drckydx = bj(1)*(-rckz/r)
                  drckydy = bj(2)*(-rckz/r)
                  drckydz = bj(3)*(-rckz/r)
                  drckzdx = bi(1)*( rckx/r) + bj(1)*( rcky/r)
                  drckzdy = bi(2)*( rckx/r) + bj(2)*( rcky/r)
                  drckzdz = bi(3)*( rckx/r) + bj(3)*( rcky/r)
                  dintSx = dintSx + drcizdx*cks*pzs + drcixdx*rckx*pxpx 
     &                        + drciydx*rcky*pypy + drcizdx*rckz*pzpz
                  dintSy = dintSy + drcizdy*cks*pzs + drcixdy*rckx*pxpx
     &                        + drciydy*rcky*pypy + drcizdy*rckz*pzpz
                  dintSz = dintSz + drcizdz*cks*pzs + drcixdz*rckx*pxpx
     &                        + drciydz*rcky*pypy + drcizdz*rckz*pzpz
                  dintSx = dintSx + cis*drckzdx*spz + rcix*drckxdx*pxpx
     &                        + rciy*drckydx*pypy + rciz*drckzdx*pzpz
                  dintSy = dintSy + cis*drckzdy*spz + rcix*drckxdy*pxpx
     &                        + rciy*drckydy*pypy + rciz*drckzdy*pzpz
                  dintSz = dintSz + cis*drckzdz*spz + rcix*drckxdz*pxpx
     &                        + rciy*drckydz*pypy + rciz*drckzdz*pzpz
                  pre = hartree * (zxri*valk + zxrk*vali) * rscale(k)
                  e = pre * intS2 / r
                  term1 = -intS2 / r3
                  intSR = 2.0d0 * intS / r
                  term2x = intSR * dintSx
                  term2y = intSR * dintSy
                  term2z = intSR * dintSz
                  
c
c     compute the force components for this interaction
c
                  frcx = pre * (xr * term1 + term2x)
                  frcy = pre * (yr * term1 + term2y)
                  frcz = pre * (zr * term1 + term2z)
c
c     compute the torque components for this interaction
c
                  ncix = 0.0d0
                  nciy = -ciz
                  nciz = ciy
                  nrcix = bi(1)*ncix + bi(2)*nciy + bi(3)*nciz
                  nrciy = bj(1)*ncix + bj(2)*nciy + bj(3)*nciz
                  nrciz = bk(1)*ncix + bk(2)*nciy + bk(3)*nciz
                  cscs = 0.0d0 * cks
                  cxcx = nrcix * rckx
                  cycy = nrciy * rcky
                  czcz = nrciz * rckz
                  cscz = 0.0d0 * rckz
                  czcs = nrciz * cks
                  tixintS = cscs * SS + cxcx * PxPx + cycy * PyPy
     &                           + czcz * PzPz + cscz * SPz + czcs * PzS
                  ncix = ciz
                  nciy = 0.0d0
                  nciz = -cix
                  nrcix = bi(1)*ncix + bi(2)*nciy + bi(3)*nciz
                  nrciy = bj(1)*ncix + bj(2)*nciy + bj(3)*nciz
                  nrciz = bk(1)*ncix + bk(2)*nciy + bk(3)*nciz
                  cscs = 0.0d0 * cks
                  cxcx = nrcix * rckx
                  cycy = nrciy * rcky
                  czcz = nrciz * rckz
                  cscz = 0.0d0 * rckz
                  czcs = nrciz * cks
                  tiyintS = cscs * SS + cxcx * PxPx + cycy * PyPy
     &                           + czcz * PzPz + cscz * SPz + czcs * PzS
                  ncix = -ciy
                  nciy = cix
                  nciz = 0.0d0
                  nrcix = bi(1)*ncix + bi(2)*nciy + bi(3)*nciz
                  nrciy = bj(1)*ncix + bj(2)*nciy + bj(3)*nciz
                  nrciz = bk(1)*ncix + bk(2)*nciy + bk(3)*nciz
                  cscs = 0.0d0 * cks
                  cxcx = nrcix * rckx
                  cycy = nrciy * rcky
                  czcz = nrciz * rckz
                  cscz = 0.0d0 * rckz
                  czcs = nrciz * cks
                  tizintS = cscs * SS + cxcx * PxPx + cycy * PyPy
     &                           + czcz * PzPz + cscz * SPz + czcs * PzS
                  nckx = 0.0d0
                  ncky = -ckz
                  nckz = cky
                  nrckx = bi(1)*nckx + bi(2)*ncky + bi(3)*nckz
                  nrcky = bj(1)*nckx + bj(2)*ncky + bj(3)*nckz
                  nrckz = bk(1)*nckx + bk(2)*ncky + bk(3)*nckz
                  cscs = cis * 0.0d0
                  cxcx = rcix * nrckx
                  cycy = rciy * nrcky
                  czcz = rciz * nrckz
                  cscz = cis * nrckz
                  czcs = rciz * 0.0d0
                  tkxintS = cscs * SS + cxcx * PxPx + cycy * PyPy
     &                           + czcz * PzPz + cscz * SPz + czcs * PzS
                  nckx = ckz
                  ncky = 0.0d0
                  nckz = -ckx
                  nrckx = bi(1)*nckx + bi(2)*ncky + bi(3)*nckz
                  nrcky = bj(1)*nckx + bj(2)*ncky + bj(3)*nckz
                  nrckz = bk(1)*nckx + bk(2)*ncky + bk(3)*nckz
                  cscs = cis * 0.0d0
                  cxcx = rcix * nrckx
                  cycy = rciy * nrcky
                  czcz = rciz * nrckz
                  cscz = cis * nrckz
                  czcs = rciz * 0.0d0
                  tkyintS = cscs * SS + cxcx * PxPx + cycy * PyPy
     &                           + czcz * PzPz + cscz * SPz + czcs * PzS
                  nckx = -cky
                  ncky = ckx
                  nckz = 0.0d0
                  nrckx = bi(1)*nckx + bi(2)*ncky + bi(3)*nckz
                  nrcky = bj(1)*nckx + bj(2)*ncky + bj(3)*nckz
                  nrckz = bk(1)*nckx + bk(2)*ncky + bk(3)*nckz
                  cscs = cis * 0.0d0
                  cxcx = rcix * nrckx
                  cycy = rciy * nrcky
                  czcz = rciz * nrckz
                  cscz = cis * nrckz
                  czcs = rciz * 0.0d0
                  tkzintS = cscs * SS + cxcx * PxPx + cycy * PyPy
     &                           + czcz * PzPz + cscz * SPz + czcs * PzS
                  preintSR = -pre * intSR
                  ttri(1) = preintSR * tixintS
                  ttri(2) = preintSR * tiyintS
                  ttri(3) = preintSR * tizintS
                  ttrk(1) = preintSR * tkxintS
                  ttrk(2) = preintSR * tkyintS
                  ttrk(3) = preintSR * tkzintS
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = fgrp * e
                     frcx = fgrp * frcx
                     frcy = fgrp * frcy
                     frcz = fgrp * frcz
                     do j = 1, 3
                        ttri(j) = fgrp * ttri(j)
                        ttrk(j) = fgrp * ttrk(j)
                     end do
                  end if
c
c     increment the overall exchange repulsion energy component
c
                  er = er + e
c
c     increment force-based gradient on first site
c
                  der(1,i) = der(1,i) - frcx
                  der(2,i) = der(2,i) - frcy
                  der(3,i) = der(3,i) - frcz
                  ter(1,i) = ter(1,i) + ttri(1)
                  ter(2,i) = ter(2,i) + ttri(2)
                  ter(3,i) = ter(3,i) + ttri(3)
c
c     increment force-based gradient and torque on second site
c
                  der(1,k) = der(1,k) + frcx
                  der(2,k) = der(2,k) + frcy
                  der(3,k) = der(3,k) + frcz
                  ter(1,k) = ter(1,k) + ttrk(1)
                  ter(2,k) = ter(2,k) + ttrk(2)
                  ter(3,k) = ter(3,k) + ttrk(3)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            rscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            rscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            rscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            rscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     resolve site torques then increment forces and virial
c
      do ii = 1, nrep
         i = irep(ii)
         call torque (i,ter(1,i),fix,fiy,fiz,der)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (rscale)
      deallocate (ter)
      return
      end
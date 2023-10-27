c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine exrepel3  --  exchange repulsion energy  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "exrepel" calculates the exchange repulsion energy
c
c     literature reference:
c
c     TBD
c
c
      subroutine exrepel
      use limits
      implicit none
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_mlist) then
         call exrepel0b
      else
         call exrepel0a
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine exrepel0a  --  exch repulsion energy via loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "exrepel0a" calculates the exchange repulsion energy using a
c     double loop
c
c
      subroutine exrepel0a
      use atoms
      use bound
      use cell
      use couple
      use energi
      use group
      use mutant
      use polar
      use reppot
      use shunt
      use units
      use usage
      use xrepel
      implicit none
      integer i,j,k
      integer ii,kk
      integer ind1,ind2,ind3
      integer jcell
      real*8 e,fgrp,taper
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
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
      real*8 bi(3)
      real*8 bj(3)
      real*8 bk(3)
      real*8, allocatable :: rscale(:)
      logical proceed,usei
      logical muti,mutk,mutik
      logical grad
      character*6 mode
c
c
c     zero out the repulsion energy contribution
c
      er = 0.0d0
      if (nxrep .eq. 0)  return
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
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         rscale(i) = 1.0d0
      end do
c
c     set cutoff distances and switching coefficients
c
      mode = 'REPULS'
      call switch (mode)
c
c     set gradient mode to false
c
      grad = .false.
c
c     calculate the exchange repulsion interaction energy term
c
      do ii = 1, nxrep-1
         i = ixrep(ii)
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
         do kk = ii+1, nxrep
            k = ixrep(kk)
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
                  e = hartree*(zxri*valk+zxrk*vali)*intS2/r*rscale(k)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall exchange repulsion energy component
c
                  er = er + e
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
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do ii = 1, nxrep
            i = ixrep(ii)
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
            do kk = ii, nxrep
               k = ixrep(kk)
               mutk = mut(k)
               proceed = .true.
               if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
               if (.not. use_intra)  proceed = .true.
               if (proceed)  proceed = (usei .or. use(k))
               if (proceed) then
                  do jcell = 2, ncell
                     xr = x(k) - xi
                     yr = y(k) - yi
                     zr = z(k) - zi
                     call imager (xr,yr,zr,jcell)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
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
                        call computeOverlap (dmpi, dmpk, dmpip, dmpkp,
     &                     0.0d0, r, grad, SS, dSS, SPz, dSPz, PzS,
     &                     dPzS, PxPx, dPxPx, PyPy, dPyPy, PzPz, dPzPz)
                        intS = cscs * SS + cxcx * PxPx + cycy * PyPy
     &                          + czcz * PzPz + cscz * SPz + czcs * PzS
                        intS2 = intS * intS
                        e = hartree*(zxri*valk+zxrk*vali)*intS2/r*
     &                              rscale(k)
c
c     use energy switching if near the cutoff distance
c
                        if (r2 .gt. cut2) then
                           r3 = r2 * r
                           r4 = r2 * r2
                           r5 = r2 * r3
                           taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                           e = e * taper
                        end if
c
c     scale the interaction based on its group membership
c
                        if (use_group)  e = e * fgrp
c
c     increment the overall exchange repulsion energy component;
c     interaction of an atom with its own image counts half
c
                        if (i .eq. k)  e = 0.5d0 * e
                        er = er + e
                     end if
                  end do
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
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (rscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine exrepel0B  --  exch repulsion energy via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "exrepel0b" calculates the exchange repulsion energy using a
c     pairwise neightbor list
c
c
      subroutine exrepel0b
      use atoms
      use bound
      use couple
      use energi
      use group
      use mutant
      use neigh
      use polar
      use reppot
      use shunt
      use units
      use usage
      use xrepel
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ind1,ind2,ind3
      real*8 e,fgrp,taper
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
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
      real*8 bi(3)
      real*8 bj(3)
      real*8 bk(3)
      real*8, allocatable :: rscale(:)
      logical proceed,usei
      logical muti,mutk,mutik
      logical grad
      character*6 mode
c
c
c     zero out the repulsion energy contribution
c
      er = 0.0d0
      if (nxrep .eq. 0)  return
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
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         rscale(i) = 1.0d0
      end do
c
c     set cutoff distances and switching coefficients
c
      mode = 'REPULS'
      call switch (mode)
c
c     set gradient mode to false
c
      grad = .false.
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(nxrep,ixrep,x,y,z,zpxr,dmppxr,rcpxr,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,r2scale,r3scale,r4scale,r5scale,
!$OMP& nelst,elst,use,use_group,use_intra,use_bounds,grad,
!$OMP& mut,cut2,off2,c0,c1,c2,c3,c4,c5)
!$OMP& firstprivate(rscale)
!$OMP& shared (er)
!$OMP DO reduction(+:er) schedule(guided)
c
c     calculate the exchange repulsion interaction energy term
c
      do ii = 1, nxrep
         i = ixrep(ii)
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
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ixrep(kk)
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
                  e = hartree*(zxri*valk+zxrk*vali)*intS2/r*rscale(k)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall exchange repulsion energy component
c
                  er = er + e
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (rscale)
      return
      end

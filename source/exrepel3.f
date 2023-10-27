c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine exrepel3  --  exch repulsion energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "exrepel3" calculates the exchange repulsion energy and partitions
c     the energy among the atoms
c
c     literature reference:
c
c     TBD
c
c
      subroutine exrepel3
      use limits
      implicit none
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_mlist) then
         call exrepel3b
      else
         call exrepel3a
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine exrepel3a  --  exch repulsion analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "exrepel3a" calculates the exchange repulsion energy and also
c     partitions the energy among the atoms using a double loop
c
c
      subroutine exrepel3a
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mutant
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
      logical header,huge
      logical grad
      character*6 mode
c
c
c     zero out the repulsion energy and partitioning terms
c
      ner = 0
      er = 0.0d0
      do i = 1, n
         aer(i) = 0.0d0
      end do
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
                  if (e .ne. 0.0d0) then
                     ner = ner + 1
                     er = er + e
                     aer(i) = aer(i) + 0.5d0*e
                     aer(k) = aer(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
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
                        if (e .ne. 0.0d0) then
                           ner = ner + 1
                           if (i .eq. k) then
                              er = er + 0.5d0*e
                              aer(i) = aer(i) + 0.5d0*e
                           else
                              er = er + e
                              aer(i) = aer(i) + 0.5d0*e
                              aer(k) = aer(k) + 0.5d0*e
                           end if
                        end if
c
c     increment the total intermolecular energy
c
                        einter = einter + e
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine exrepel3b  --  exch repulsion analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "exrepel3b" calculates the exchange repulsion energy and also
c     partitions the energy among the atoms using a neighbor list
c
c
      subroutine exrepel3b
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mutant
      use neigh
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
      logical header,huge
      logical grad
      character*6 mode
c
c
c     zero out the repulsion energy and partitioning terms
c
      ner = 0
      er = 0.0d0
      do i = 1, n
         aer(i) = 0.0d0
      end do
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
!$OMP& mut,cut2,off2,c0,c1,c2,c3,c4,c5,molcule,name,
!$OMP& verbose,debug,header,iout)
!$OMP& firstprivate(rscale)
!$OMP& shared (er,ner,aer,einter)
!$OMP DO reduction(+:er,ner,aer,einter) schedule(guided)
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
                  if (e .ne. 0.0d0) then
                     ner = ner + 1
                     er = er + e
                     aer(i) = aer(i) + 0.5d0*e
                     aer(k) = aer(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
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
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine solvcoeff  --  solve for orbital coefficients  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "solvcoeff" finds the coefficients for the pseudo orbital
c
c
      subroutine solvcoeff
      use repel
      use xrepel
      implicit none
      integer ii,k
      integer ind1,ind2,ind3
      real*8 cr,cs
      real*8 pcoeff(3)
      real*8 ppole(3)
      real*8 p2p1,p3p1
      logical l1,l2,l3
c
c
c     determine pseudo orbital coefficients
c
      do ii = 1, nxrep
         do k = 1, 3
            pcoeff(k) = 0.0d0
         end do
         ppole(1) = xrepole(2,ii)
         ppole(2) = xrepole(3,ii)
         ppole(3) = xrepole(4,ii)
         cr = crpxr(ii)
         l1 = (abs(ppole(1)) < 1.0d-10)
         l2 = (abs(ppole(2)) < 1.0d-10)
         l3 = (abs(ppole(3)) < 1.0d-10)
c
c     case for no dipole
c
         if (l1.and.l2.and.l3) then
            cs = 1.0d0
            ind1 = 1
            ind2 = 2
            ind3 = 3
c
c     case for p orbital coefficients set to 0
c
         else if (cr < 1.0d-10) then
            cs = 1.0d0
            ind1 = 1
            ind2 = 2
            ind3 = 3
c
c     determine normalized coefficients
c
         else
            cs = 1.0d0 / sqrt(1.0d0 + cr)
c
c     determine index for largest absolute dipole component
c
            ind1 = maxloc(abs(ppole), dim=1)
            ind2 = mod(ind1,3) + 1
            ind3 = mod(ind1+1,3) + 1
            p2p1 = ppole(ind2) / ppole(ind1)
            p3p1 = ppole(ind3) / ppole(ind1)
            pcoeff(ind1) = cs * sqrt(cr / (1.0d0 + p2p1**2 + p3p1**2))
            if (ppole(ind1) < 0.0d0) then
               pcoeff(ind1) = -pcoeff(ind1)
            end if
            pcoeff(ind2) = pcoeff(ind1) * p2p1
            pcoeff(ind3) = pcoeff(ind1) * p3p1
         end if
         cpxr(1,ii) = cs
         cpxr(ind1+1,ii) = pcoeff(ind1)
         cpxr(ind2+1,ii) = pcoeff(ind2)
         cpxr(ind3+1,ii) = pcoeff(ind3)
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine rotcoeff  --  rotate orbital coefficients  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "rotcoeff" rotates the coefficients for the pseudo orbital
c     to global coordinates
c
c
      subroutine rotcoeff
      use mpole
      use repel
      use xrepel
      implicit none
      integer isite
      integer i,j
      real*8 a(3,3)
      real*8 cpole(4)
      logical planar
      character*8 axetyp
c
c
c     rotate pseudo orbital coefficients
c
      do isite = 1, nxrep
c
c     determine rotation matrix
c
         call rotmat (isite,a,planar)
c
c     copy local frame coefficients
c
         do i = 1, 4
            cpole(i) = cpxr(i,isite)
         end do
         if (planar) then
            axetyp = polaxe(isite)
            if (axetyp .eq. 'Z-Bisect') then
               cpole(2) = 0.0d0
            else if (axetyp .eq. '3-Fold') then
               do i = 2, 4
                  cpole(i) = 0.0d0
               end do
            end if
         end if
c
c     s orbital coefficients have the same value in any coordinate frame
c
         rcpxr(1,isite) = cpole(1)
c
c     rotate the p orbital coefficients to the global coordinate frame
c
         do i = 2, 4
            rcpxr(i,isite) = 0.0d0
            do j = 2, 4
               rcpxr(i,isite) = rcpxr(i,isite) + cpole(j)*a(i-1,j-1)
            end do
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine computeOverlap  --  computes overlap integrals  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "computeOverlap" computes overlap integrals
c
c
      subroutine computeOverlap (a, b, ap, bp, z1, z2, grad, SS, dSS,
     &      SPz, dSPz, PzS, dPzS, PxPx, dPxPx, PyPy, dPyPy, PzPz, dPzPz)
      implicit none
      real*8 a,b,ap,bp
      real*8 z1,z2
      real*8 SS,SPz,PzS
      real*8 PxPx,PyPy,PzPz
      real*8 dSS,dSPz,dPzS
      real*8 dPxPx,dPyPy,dPzPz
      real*8 eps, diffa, diffb
      logical grad
c
c
      eps = 1.0d-10
      diffa = abs(a - ap)
      diffb = abs(b - bp)
      if (diffa.lt.eps .and. diffb.lt.eps) then
         call overlapAll(a, b, z1, z2, grad, SS, dSS, SPz, dSPz,
     &                              PzS, dPzS, PxPx, dPxPx, PzPz, dPzPz)
      else
         call overlapSS(a, b, z1, z2, grad, SS, dSS)
         call overlapSPz(a, bp, z1, z2, grad, SPz, dSPz)
         call overlapPzS(ap, b, z1, z2, grad, PzS, dPzS)
         call overlapPxPx(ap, bp, z1, z2, grad, PxPx, dPxPx)
         call overlapPzPz(ap, bp, z1, z2, grad, PzPz, dPzPz)
      end if
      PyPy = PxPx
      dPyPy = dPxPx
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine overlapAll  --  all overlap integrals  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "overlapAll" computes all the overlap integrals
c
c
      subroutine overlapAll (a, b, z1, z2, grad, SS, dSS, SPz, dSPz,
     &                              PzS, dPzS, PxPx, dPxPx, PzPz, dPzPz)
      implicit none
      real*8 SS,SPz,PzS
      real*8 PxPx,PzPz
      real*8 dSS,dSPz,dPzS
      real*8 dPxPx,dPzPz
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4
      real*8 exp1
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4,rhoA5
      real*8 rhoB2,rhoB3,rhoB4,rhoB5
      real*8 kappam,kappam2
      real*8 kappap,kappap2
      real*8 expA,expB
      logical grad
c
c
      diff = abs(a - b)
      eps = 0.001d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         exp1 = exp(-rho)
         SS = (1.0d0 + rho + rho2 / 3.0d0) * exp1
         SPz = -0.5d0 * rho * (1.0d0 + rho + rho2 / 3.0d0) * exp1
         PzS = -SPz
         PxPx = (1.0d0 + rho + 2.0d0/5.0d0 * rho2 + rho3/15.0d0) * exp1
         PzPz = -(-1.0d0 - rho - rho2 / 5.0d0 + 2.0d0/15.0d0 * rho3
     &                                           + rho4 / 15.0d0) * exp1
         if (grad) then
         dSS = -1.0d0/3.0d0 * a * rho * (1.0d0 + rho) * exp1
         dSPz = -0.5d0 * a * (1.0d0 + rho - rho3 / 3.0d0) * exp1
         dPzS = -dSPz
         dPxPx = -0.2d0 * a * rho * (1.0d0 + rho + rho2 / 3.0d0)
     &                                                            * exp1
         dPzPz = -0.6d0 * a * rho * (1.0d0 + rho + 2.0d0/9.0d0 * rho2
     &                                            - rho3 / 9.0d0) * exp1
         end if
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho2 = rho * rho
         rho3 = rho2 * rho
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappam2 = kappam**2
         kappap2 = kappap**2
         expA = exp(-rhoA)
         expB = exp(-rhoB)
         pre = sqrt(1.0d0 - tau**2) / (tau * rho)
         term1 =-kappam * (2.0d0 * kappap + rhoA) * expA
         term2 = kappap * (2.0d0 * kappam + rhoB) * expB
         SS = pre * (term1 + term2)
         pre = sqrt((1.0d0 + tau) / (1.0d0 - tau)) / (tau * rho2)
         term1 =-kappam2 * (6.0d0 * kappap * (1.0d0 + rhoA)
     &      + 2.0d0 * rhoA2) * expA
         term2 = kappap * (6.0d0 * kappam2 * (1.0d0 + rhoB)
     &      + 4.0d0 * kappam * rhoB2 + rhoB3) * expB
         SPz = -pre * (term1 + term2)
         pre = sqrt((1.0d0 - tau) / (1.0d0 + tau)) / (-tau * rho2)
         term1 =-kappap2 * (6.0d0 * kappam * (1.0d0 + rhoB)
     &      + 2.0d0 * rhoB2) * expB
         term2 = kappam * (6.0d0 * kappap2 * (1.0d0 + rhoA)
     &      + 4.0d0 * kappap * rhoA2 + rhoA3) * expA
         PzS = pre * (term1 + term2)
         pre = 1.0d0 / (sqrt(1.0d0 - tau**2) * tau * rho3)
         term1 =-kappam2 * (24.0d0 * kappap2 * (1.0d0 + rhoA)
     &      + 12.0d0 * kappap * rhoA2 + 2.0d0 * rhoA3) * expA
         term2 = kappap2 * (24.0d0 * kappam2 * (1.0d0 + rhoB)
     &      + 12.0d0 * kappam * rhoB2 + 2.0d0 * rhoB3) * expB
         PxPx = pre * (term1 + term2)
         term1 =-kappam2 * (48.0d0 * kappap2 * (1.0d0 + rhoA
     &      + 0.5d0 * rhoA2) + 2.0d0 * (5.0d0 + 6.0d0 * kappa) * rhoA3
     &      + 2.0d0 * rhoA4) * expA
         term2 = kappap2 * (48.0d0 * kappam2 * (1.0d0 + rhoB
     &      + 0.5d0 * rhoB2) + 2.0d0 * (5.0d0 - 6.0d0 * kappa) * rhoB3
     &      + 2.0d0 * rhoB4) * expB
         PzPz = -pre * (term1 + term2)
         if (grad) then
            rhoA5 = rhoA4 * rhoA
            rhoB5 = rhoB4 * rhoB
            pre = sqrt(1.0d0 - tau**2) / (tau * rho)
            term1 = kappam * (2.0d0 * kappap * (1.0d0 + rhoA) + rhoA2)
     &                                                            * expA
            term2 =-kappap * (2.0d0 * kappam * (1.0d0 + rhoB) + rhoB2)
     &                                                            * expB
            dSS = pre / r * (term1 + term2)
            pre = sqrt((1.0d0 + tau) / (1.0d0 - tau)) / (tau * rho2)
            term1 = 2.0d0 * kappam2 * (6.0d0 * kappap * (1.0d0 + rhoA
     &         + 0.5d0 * rhoA2) + rhoA3) * expA
            term2 = kappap * (-12.0d0 * kappam2 * (1.0d0 + rhoB + 0.5d0
     &         * rhoB2) + (1.0d0 - 4.0d0 * kappam) * rhoB3 - rhoB4)
     &         * expB
            dSPz = -pre / r * (term1 + term2)
            pre = sqrt((1.0d0 - tau) / (1.0d0 + tau)) / (-tau * rho2)
            term1 = 2.0d0 * kappap2 * (6.0d0 * kappam * (1.0d0 + rhoB
     &         + 0.5d0 * rhoB2) + rhoB3) * expB
            term2 = kappam * (-12.0d0 * kappap2 * (1.0d0 + rhoA + 0.5d0
     &         * rhoA2) + (1.0d0 - 4.0d0 * kappap) * rhoA3 - rhoA4)
     &         * expA
            dPzS = pre / r * (term1 + term2)
            pre = 1.0d0 / (sqrt(1.0d0 - tau**2) * tau * rho3)
            term1 = 2.0d0 * kappam2 * (36.0d0 * kappap2 * (1.0d0 + rhoA)
     &         + 6.0d0 * kappap * (1.0d0 + 2.0d0 * kappap) * rhoA2
     &         + 6.0d0 * kappap * rhoA3 + rhoA4) * expA
            term2 =-2.0d0 * kappap2 * (36.0d0 * kappam2 * (1.0d0 + rhoB)
     &         + 6.0d0 * kappam * (1.0d0 + 2.0d0 * kappam) * rhoB2
     &         + 6.0d0 * kappam * rhoB3 + rhoB4) * expB
            dPxPx = pre / r * (term1 + term2)
            term1 = kappam2 * (72.0d0 * kappap2 * (1.0d0 + rhoA
     &         + 0.5d0 * rhoA2 + rhoA3 / 6.0d0) + 2.0d0 * (2.0d0
     &         + 3.0d0 * kappa) * rhoA4 + rhoA5) * expA
            term2 =-kappap2 * (72.0d0 * kappam2 * (1.0d0 + rhoB
     &         + 0.5d0 * rhoB2 + rhoB3 / 6.0d0) + 2.0d0 * (2.0d0
     &         - 3.0d0 * kappa) * rhoB4 + rhoB5) * expB
            dPzPz = -2.0d0 * pre / r * (term1 + term2)
         end if
      end if
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine overlapSS  --  SS overlap integral  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "overlapSS" computes the overlap integral
c
c
      subroutine overlapSS (a, b, z1, z2, grad, SS, dSS)
      implicit none
      real*8 over,dover
      real*8 SS,dSS
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2
      real*8 exp1
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoB2
      real*8 kappam,kappap
      real*8 expA,expB
      logical grad
c
c
      diff = abs(a - b)
      eps = 0.001d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         exp1 = exp(-rho)
         over = (1.0d0 + rho + rho2 / 3.0d0) * exp1
         if (grad) then
            dover = -1.0d0/3.0d0 * a * rho * (1.0d0 + rho) * exp1
         end if
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         expA = exp(-rhoA)
         expB = exp(-rhoB)
         pre = sqrt(1.0d0 - tau**2) / (tau * rho)
         term1 =-kappam * (2.0d0 * kappap + rhoA) * expA
         term2 = kappap * (2.0d0 * kappam + rhoB) * expB
         over = pre * (term1 + term2)
         if (grad) then
            rhoA2 = rhoA * rhoA
            rhoB2 = rhoB * rhoB
            term1 = kappam * (2.0d0 * kappap * (1.0d0 + rhoA) + rhoA2)
     &                                                            * expA
            term2 =-kappap * (2.0d0 * kappam * (1.0d0 + rhoB) + rhoB2)
     &                                                            * expB
            dover = pre / r * (term1 + term2)
         end if
      end if
      SS = over
      dSS = dover
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine overlapSPz  --  SPz overlap integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "overlapSPz" computes the overlap integral
c
c
      subroutine overlapSPz (a, b, z1, z2, grad, SPz, dSPz)
      implicit none
      real*8 over,dover
      real*8 SPz,dSPz
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3
      real*8 exp1
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3
      real*8 rhoB2,rhoB3,rhoB4
      real*8 kappam,kappam2
      real*8 kappap,kappap2
      real*8 expA,expB
      logical grad
c
c
      diff = abs(a - b)
      eps = 0.001d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         exp1 = exp(-rho)
         over = 0.5d0 * rho * (1.0d0 + rho + rho2 / 3.0d0) * exp1
         if (grad) then
            rho3 = rho2 * rho
            dover = 0.5d0 * a * (1.0d0 + rho - rho3 / 3.0d0) * exp1
         end if
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho2 = rho**2
         rhoA2 = rhoA * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappam2 = kappam**2
         kappap2 = kappap**2
         expA = exp(-rhoA)
         expB = exp(-rhoB)
         pre = sqrt((1.0d0 + tau) / (1.0d0 - tau)) / (tau * rho2)
         term1 =-kappam2 * (6.0d0 * kappap * (1.0d0 + rhoA)
     &      + 2.0d0 * rhoA2) * expA
         term2 = kappap * (6.0d0 * kappam2 * (1.0d0 + rhoB)
     &      + 4.0d0 * kappam * rhoB2 + rhoB3) * expB
         over = pre * (term1 + term2)
         if (grad) then
            rhoA3 = rhoA2 * rhoA
            rhoB4 = rhoB3 * rhoB
            term1 = 2.0d0 * kappam2 * (6.0d0 * kappap * (1.0d0 + rhoA
     &         + 0.5d0 * rhoA2) + rhoA3) * expA
            term2 = kappap * (-12.0d0 * kappam2 * (1.0d0 + rhoB + 0.5d0
     &         * rhoB2) + (1.0d0 - 4.0d0 * kappam) * rhoB3 - rhoB4)
     &         * expB
            dover = pre / r * (term1 + term2)
         end if
      end if
      SPz = -over
      dSPz = -dover
      if (z1 > z2) then
         SPz = -SPz
         dSPz = -dSPz
      end if
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine overlapPzS  --  PzS overlap integral  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "overlapPzS" computes the overlap integral
c
c
      subroutine overlapPzS (a, b, z1, z2, grad, PzS, dPzS)
      implicit none
      real*8 PzS,dPzS
      real*8 a,b,z1,z2
      logical grad
c
c
      call overlapSPz(b, a, z2, z1, grad, PzS, dPzS)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine overlapPxPx  --  PxPx overlap integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "overlapPxPx" computes the overlap integral
c
c
      subroutine overlapPxPx (a, b, z1, z2, grad, PxPx, dPxPx)
      implicit none
      real*8 over,dover
      real*8 PxPx,dPxPx
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3
      real*8 exp1
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4
      real*8 rhoB2,rhoB3,rhoB4
      real*8 kappam,kappam2
      real*8 kappap,kappap2
      real*8 expA,expB
      logical grad
c
c
      diff = abs(a - b)
      eps = 0.001d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         exp1 = exp(-rho)
         over = (1.0d0 + rho + 2.0d0/5.0d0 * rho2 + rho3/15.0d0) * exp1
         if (grad) then
            dover = -0.2d0 * a * rho * (1.0d0 + rho + rho2 / 3.0d0)
     &                                                            * exp1
         end if
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappam2 = kappam**2
         kappap2 = kappap**2
         expA = exp(-rhoA)
         expB = exp(-rhoB)
         pre = 1.0d0 / (sqrt(1.0d0 - tau**2) * tau * rho3)
         term1 =-kappam2 * (24.0d0 * kappap2 * (1.0d0 + rhoA)
     &      + 12.0d0 * kappap * rhoA2 + 2.0d0 * rhoA3) * expA
         term2 = kappap2 * (24.0d0 * kappam2 * (1.0d0 + rhoB)
     &      + 12.0d0 * kappam * rhoB2 + 2.0d0 * rhoB3) * expB
         over = pre * (term1 + term2)
         if (grad) then
            rhoA4 = rhoA3 * rhoA
            rhoB4 = rhoB3 * rhoB
            term1 = 2.0d0 * kappam2 * (36.0d0 * kappap2 * (1.0d0 + rhoA)
     &         + 6.0d0 * kappap * (1.0d0 + 2.0d0 * kappap) * rhoA2
     &         + 6.0d0 * kappap * rhoA3 + rhoA4) * expA
            term2 =-2.0d0 * kappap2 * (36.0d0 * kappam2 * (1.0d0 + rhoB)
     &         + 6.0d0 * kappam * (1.0d0 + 2.0d0 * kappam) * rhoB2
     &         + 6.0d0 * kappam * rhoB3 + rhoB4) * expB
            dover = pre / r * (term1 + term2)
         end if
      end if
      PxPx = over
      dPxPx = dover
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine overlapPzPz  --  PzPz overlap integral  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "overlapPzPz" computes the overlap integral
c
c
      subroutine overlapPzPz (a, b, z1, z2, grad, PzPz, dPzPz)
      implicit none
      real*8 over,dover
      real*8 PzPz,dPzPz
      real*8 a,b,z1,z2
      real*8 diff,eps,zr,r
      real*8 rho,rho2,rho3,rho4
      real*8 exp1
      real*8 alpha,tau,rhoA,rhoB
      real*8 a2,b2,kappa
      real*8 pre,term1,term2
      real*8 rhoA2,rhoA3,rhoA4,rhoA5
      real*8 rhoB2,rhoB3,rhoB4,rhoB5
      real*8 kappam,kappam2
      real*8 kappap,kappap2
      real*8 expA,expB
      logical grad
c
c
      diff = abs(a - b)
      eps = 0.001d0
      zr = z2 - z1
      r = abs(zr)
      if (diff .lt. eps) then
         rho = a * r
         rho2 = rho * rho
         rho3 = rho2 * rho
         rho4 = rho3 * rho
         exp1 = exp(-rho)
         over = (-1.0d0 - rho - rho2 / 5.0d0 + 2.0d0/15.0d0 * rho3
     &                                           + rho4 / 15.0d0) * exp1
         if (grad) then
         dover = 0.6d0 * a * rho * (1.0d0 + rho + 2.0d0/9.0d0 * rho2
     &                                            - rho3 / 9.0d0) * exp1
         end if
      else
         alpha = 1.0d0 / 2.0d0 * (a + b)
         tau = (a - b) / (a + b)
         rho = alpha * r
         rhoA = a * r
         rhoB = b * r
         a2 = a * a
         b2 = b * b
         kappa = (a2 + b2) / (a2 - b2)
         rho3 = rho**3
         rhoA2 = rhoA * rhoA
         rhoA3 = rhoA2 * rhoA
         rhoA4 = rhoA3 * rhoA
         rhoB2 = rhoB * rhoB
         rhoB3 = rhoB2 * rhoB
         rhoB4 = rhoB3 * rhoB
         kappam = 1.0d0 - kappa
         kappap = 1.0d0 + kappa
         kappam2 = kappam**2
         kappap2 = kappap**2
         expA = exp(-rhoA)
         expB = exp(-rhoB)
         pre = 1.0d0 / (sqrt(1.0d0 - tau**2) * tau * rho3)
         term1 =-kappam2 * (48.0d0 * kappap2 * (1.0d0 + rhoA
     &      + 0.5d0 * rhoA2) + 2.0d0 * (5.0d0 + 6.0d0 * kappa) * rhoA3
     &      + 2.0d0 * rhoA4) * expA
         term2 = kappap2 * (48.0d0 * kappam2 * (1.0d0 + rhoB
     &      + 0.5d0 * rhoB2) + 2.0d0 * (5.0d0 - 6.0d0 * kappa) * rhoB3
     &      + 2.0d0 * rhoB4) * expB
         over = pre * (term1 + term2)
         if (grad) then
            rhoA5 = rhoA4 * rhoA
            rhoB5 = rhoB4 * rhoB
            term1 = kappam2 * (72.0d0 * kappap2 * (1.0d0 + rhoA
     &         + 0.5d0 * rhoA2 + rhoA3 / 6.0d0) + 2.0d0 * (2.0d0
     &         + 3.0d0 * kappa) * rhoA4 + rhoA5) * expA
            term2 =-kappap2 * (72.0d0 * kappam2 * (1.0d0 + rhoB
     &         + 0.5d0 * rhoB2 + rhoB3 / 6.0d0) + 2.0d0 * (2.0d0
     &         - 3.0d0 * kappa) * rhoB4 + rhoB5) * expB
            dover = 2.0d0 * pre / r * (term1 + term2)
         end if
      end if
      PzPz = -over
      dPzPz = -dover
      return
      end

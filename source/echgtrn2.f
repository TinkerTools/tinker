c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine echgtrn2  --  atomwise charge transfer Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "echgtrn2" calculates the second derivatives of the charge
c     transfer energy using a double loop over relevant atom pairs
c
c
      subroutine echgtrn2 (iatom)
      use atoms
      use bound
      use cell
      use chgpot
      use chgtrn
      use couple
      use ctrpot
      use group
      use hessn
      use mplpot
      use mpole
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,iii
      integer iatom,jcell
      integer nlist,list(5)
      real*8 e,dedr,d2edr2
      real*8 term,f,fgrp
      real*8 termx,termy,termz
      real*8 rr1,r,r2
      real*8 r3,r4,r5
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 chgi,chgk
      real*8 chgik
      real*8 alphai,alphak
      real*8 alphaik
      real*8 alphai2,alphak2
      real*8 expi,expk
      real*8 expik
      real*8 taper,dtaper
      real*8 d2taper
      real*8 d2e(3,3)
      real*8, allocatable :: mscale(:)
      logical proceed,usei
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHGTRN'
      call switch (mode)
c
c     check to see if the atom of interest is a vdw site
c
      nlist = 0
      do k = 1, npole
         if (ipole(k) .eq. iatom) then
            nlist = nlist + 1
            list(nlist) = k
            goto 10
         end if
      end do
      return
   10 continue
c
c     calculate the charge transfer energy term
c
      do iii = 1, nlist
         ii = list(iii)
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         chgi = chgct(ii)
         alphai = dmpct(ii)
         if (alphai .eq. 0.0d0)  alphai = 1000.0d0
         alphai2 = alphai * alphai
         usei = use(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = m2scale
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = m3scale
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = m4scale
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = 1, npole
            k = ipole(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rr1 = 1.0d0 / r
                  chgk = chgct(kk)
                  alphak = dmpct(kk)
                  if (alphak .eq. 0.0d0)  alphak = 1000.0d0
                  if (ctrntyp .eq. 'SEPARATE') then
                     alphak2 = alphak * alphak
                     expi = exp(-alphai*r)
                     expk = exp(-alphak*r)
                     e = -chgi*expk - chgk*expi
                     dedr = chgi*expk*alphak + chgk*expi*alphai
                     d2edr2 = -chgi*expk*alphak2 - chgk*expi*alphai2
                  else
                     chgik = sqrt(abs(chgi*chgk))
                     alphaik = 0.5d0 * (alphai+alphak)
                     expik = exp(-alphaik*r)
                     e = -chgik * expik
                     dedr = -e * alphaik
                     d2edr2 = -dedr * alphaik
                  end if
                  e = f * e * mscale(k)
                  dedr = f * dedr * mscale(k)
                  d2edr2 = f * d2edr2 * mscale(k)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r3 * r2
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                            + 6.0d0*c3*r + 2.0d0*c2
                     d2e = e*d2taper + 2.0d0*dedr*dtaper + d2edr2*taper
                     dedr = e*dtaper + dedr*taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     dedr = dedr * fgrp
                     d2edr2 = d2edr2 * fgrp
                  end if
c
c     set the chain rule terms for the Hessian elements
c
                  if (r2 .eq. 0.0d0) then
                     dedr = 0.0d0
                     term = 0.0d0
                  else
                     dedr = dedr / r
                     term = (d2edr2-dedr) / r2
                  end if
                  termx = term * xr
                  termy = term * yr
                  termz = term * zr
                  d2e(1,1) = termx*xr + dedr
                  d2e(1,2) = termx*yr
                  d2e(1,3) = termx*zr
                  d2e(2,1) = d2e(1,2)
                  d2e(2,2) = termy*yr + dedr
                  d2e(2,3) = termy*zr
                  d2e(3,1) = d2e(1,3)
                  d2e(3,2) = d2e(2,3)
                  d2e(3,3) = termz*zr + dedr
c
c     increment diagonal and non-diagonal Hessian elements
c
                  do j = 1, 3
                     hessx(j,i) = hessx(j,i) + d2e(1,j)
                     hessy(j,i) = hessy(j,i) + d2e(2,j)
                     hessz(j,i) = hessz(j,i) + d2e(3,j)
                     hessx(j,k) = hessx(j,k) - d2e(1,j)
                     hessy(j,k) = hessy(j,k) - d2e(2,j)
                     hessz(j,k) = hessz(j,k) - d2e(3,j)
                  end do
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = 1.0d0
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
         do iii = 1, nlist
            ii = list(iii)
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            chgi = chgct(ii)
            alphai = dmpct(ii)
            if (alphai .eq. 0.0d0)  alphai = 1000.0d0
            alphai2 = alphai * alphai
            usei = use(i)
c
c     set exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               mscale(i12(j,i)) = m2scale
            end do
            do j = 1, n13(i)
               mscale(i13(j,i)) = m3scale
            end do
            do j = 1, n14(i)
               mscale(i14(j,i)) = m4scale
            end do
            do j = 1, n15(i)
               mscale(i15(j,i)) = m5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do kk = 1, npole
               k = ipole(kk)
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
                        rr1 = 1.0d0 / r
                        chgk = chgct(kk)
                        alphak = dmpct(kk)
                        if (alphak .eq. 0.0d0)  alphak = 1000.0d0
                        if (ctrntyp .eq. 'SEPARATE') then
                           alphak2 = alphak * alphak
                           expi = exp(-alphai*r)
                           expk = exp(-alphak*r)
                           e = -chgi*expk - chgk*expi
                           dedr = chgi*expk*alphak + chgk*expi*alphai
                           d2edr2 = -chgi*expk*alphak2
     &                                 - chgk*expi*alphai2
                        else
                           chgik = sqrt(abs(chgi*chgk))
                           alphaik = 0.5d0 * (alphai+alphak)
                           expik = exp(-alphaik*r)
                           e = -chgik * expik
                           dedr = -e * alphaik
                           d2edr2 = -dedr * alphaik
                        end if
                        e = f * e * mscale(k)
                        dedr = f * dedr * mscale(k)
                        d2edr2 = f * d2edr2 * mscale(k)
c
c     use energy switching if near the cutoff distance
c
                        if (r2 .gt. cut2) then
                           r3 = r2 * r
                           r4 = r2 * r2
                           r5 = r3 * r2
                           taper = c5*r5 + c4*r4 + c3*r3
     &                                + c2*r2 + c1*r + c0
                           dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                                 + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                           d2taper = 20.0d0*c5*r3 + 12.0d0*c4*r2
     &                                  + 6.0d0*c3*r + 2.0d0*c2
                           d2e = e*d2taper + 2.0d0*dedr*dtaper
     &                              + d2edr2*taper
                           dedr = e*dtaper + dedr*taper
                        end if
c
c     scale the interaction based on its group membership
c
                        if (use_group) then
                           dedr = dedr * fgrp
                           d2edr2 = d2edr2 * fgrp
                        end if
c
c     set the chain rule terms for the Hessian elements
c
                        if (r2 .eq. 0.0d0) then
                           dedr = 0.0d0
                           term = 0.0d0
                        else
                           dedr = dedr / r
                           term = (d2edr2-dedr) / r2
                        end if
                        termx = term * xr
                        termy = term * yr
                        termz = term * zr
                        d2e(1,1) = termx*xr + dedr
                        d2e(1,2) = termx*yr
                        d2e(1,3) = termx*zr
                        d2e(2,1) = d2e(1,2)
                        d2e(2,2) = termy*yr + dedr
                        d2e(2,3) = termy*zr
                        d2e(3,1) = d2e(1,3)
                        d2e(3,2) = d2e(2,3)
                        d2e(3,3) = termz*zr + dedr
c
c     increment diagonal and non-diagonal Hessian elements
c
                        do j = 1, 3
                           hessx(j,i) = hessx(j,i) + d2e(1,j)
                           hessy(j,i) = hessy(j,i) + d2e(2,j)
                           hessz(j,i) = hessz(j,i) + d2e(3,j)
                           hessx(j,k) = hessx(j,k) - d2e(1,j)
                           hessy(j,k) = hessy(j,k) - d2e(2,j)
                           hessz(j,k) = hessz(j,k) - d2e(3,j)
                        end do
                     end if
                  end do
               end if
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               mscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               mscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               mscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               mscale(i15(j,i)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end

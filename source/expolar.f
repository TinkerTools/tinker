c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine expolar  --  ExchPol variable polarizability  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "expolar" calculates the variable polarizability due to exchange
c     polarization
c
c     literature reference:
c
c
      subroutine expolar (polscale,invpolscale)
      use limits
      use mpole
      implicit none
      real*8 polscale(3,3,npole)
      real*8 invpolscale(3,3,npole)
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_mlist) then
         call expolar0b (polscale,invpolscale)
      else
         call expolar0a (polscale,invpolscale)
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine expolar0a  --  variable polarizability via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "expolar0a" calculates the variable polarizability due to exchange
c     polarization using a double loop
c
c
      subroutine expolar0a (polscale,invpolscale)
      use atoms
      use bound
      use cell
      use couple
      use expol
      use mpole
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k,m
      integer ii,kk
      integer jcell
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
      real*8 det
      real*8 sizi,sizk,sizik
      real*8 alphai,alphak
      real*8 springi,springk
      real*8 s2,ds2
      real*8 p33i, p33k
      real*8 kS2i(3,3)
      real*8 kS2k(3,3)
      real*8 taper
      real*8 ps(3,3)
      real*8, allocatable :: pscale(:)
      real*8 polscale(3,3,npole)
      real*8 invpolscale(3,3,npole)
      logical epli,eplk
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     set the switching function coefficients
c
      mode = 'REPULS'
      call switch (mode)
c
c     set polscale and invpolscale to the identity matrix
c
      do i = 1, npole
         polscale(1,1,i) = 1.0d0
         polscale(2,1,i) = 0.0d0
         polscale(3,1,i) = 0.0d0
         polscale(1,2,i) = 0.0d0
         polscale(2,2,i) = 1.0d0
         polscale(3,2,i) = 0.0d0
         polscale(1,3,i) = 0.0d0
         polscale(2,3,i) = 0.0d0
         polscale(3,3,i) = 1.0d0
         invpolscale(1,1,i) = 1.0d0
         invpolscale(2,1,i) = 0.0d0
         invpolscale(3,1,i) = 0.0d0
         invpolscale(1,2,i) = 0.0d0
         invpolscale(2,2,i) = 1.0d0
         invpolscale(3,2,i) = 0.0d0
         invpolscale(1,3,i) = 0.0d0
         invpolscale(2,3,i) = 0.0d0
         invpolscale(3,3,i) = 1.0d0
      end do
c
c     set array needed to scale atom and group interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     find the variable polarizability
c
      do ii = 1, npole-1
         i = ipole(ii)
         springi = kpep(ii)
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
         do kk = ii+1, npole
            k = ipole(kk)
            eplk = lpep(kk)
            if (epli.or.eplk) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  springk = kpep(kk)
                  sizk = prepep(kk)
                  alphak = dmppep(kk)
                  sizik = sizi*sizk
                  call dampep (r,sizik,alphai,alphak,s2,ds2)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     s2 = s2 * taper
                  end if
                  p33i = springi*s2*pscale(k)
                  p33k = springk*s2*pscale(k)
                  call exrotate(i,k,xr,yr,zr,p33i,p33k,kS2i,kS2k)
                  if (epli) then
                     do j = 1, 3
                        do m = 1, 3
                           polscale(j,m,ii) = polscale(j,m,ii)+kS2i(j,m)
                        end do
                     end do
                  end if
                  if (eplk) then
                     do j = 1, 3
                        do m = 1, 3
                           polscale(j,m,kk) = polscale(j,m,kk)+kS2k(j,m)
                        end do
                     end do
                  end if
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
c     calculate interaction energy with other unit cells
c
         do ii = 1, npole
            i = ipole(ii)
            springi = kpep(ii)
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
               if (epli.or.eplk) then
                  do jcell = 2, ncell
                     xr = x(k) - x(i)
                     yr = y(k) - y(i)
                     zr = z(k) - z(i)
                     call imager (xr,yr,zr,jcell)
                     r2 = xr*xr + yr*yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        springk = kpep(kk)
                        sizk = prepep(kk)
                        alphak = dmppep(kk)
                        sizik = sizi*sizk
                        call dampep (r,sizik,alphai,alphak,s2,ds2)
c
c     use energy switching if near the cutoff distance
c
                        if (r2 .gt. cut2) then
                           r3 = r2 * r
                           r4 = r2 * r2
                           r5 = r2 * r3
                           taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                           s2 = s2 * taper
                        end if
c
c     interaction of an atom with its own image counts half
c
                        if (i .eq. k)  s2 = 0.5d0 * s2
                        p33i = springi*s2*pscale(k)
                        p33k = springk*s2*pscale(k)
                        call exrotate(i,k,xr,yr,zr,p33i,p33k,kS2i,kS2k)
                        if (epli) then
                           do j = 1, 3
                              do m = 1, 3
                                 polscale(j,m,ii) = polscale(j,m,ii)
     &                                            + kS2i(j,m)
                              end do
                           end do
                        end if
                        if (eplk) then
                           do j = 1, 3
                              do m = 1, 3
                                 polscale(j,m,kk) = polscale(j,m,kk)
     &                                            + kS2k(j,m)
                              end do
                           end do
                        end if
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
c     invert polscale matrix
c
      do ii = 1, npole
         do j = 1, 3
            do m = 1, 3
               ps(j,m) = polscale(j,m,ii)
            end do
         end do
         det = ps(1,1)*(ps(2,2)*ps(3,3) - ps(3,2)*ps(2,3))
     &       - ps(1,2)*(ps(2,1)*ps(3,3) - ps(2,3)*ps(3,1))
     &       + ps(1,3)*(ps(2,1)*ps(3,2) - ps(2,2)*ps(3,1))
         invpolscale(1,1,ii) = (ps(2,2)*ps(3,3)-ps(3,2)*ps(2,3))/det
         invpolscale(1,2,ii) = (ps(1,3)*ps(3,2)-ps(1,2)*ps(3,3))/det
         invpolscale(1,3,ii) = (ps(1,2)*ps(2,3)-ps(1,3)*ps(2,2))/det
         invpolscale(2,1,ii) = (ps(2,3)*ps(3,1)-ps(2,1)*ps(3,3))/det
         invpolscale(2,2,ii) = (ps(1,1)*ps(3,3)-ps(1,3)*ps(3,1))/det
         invpolscale(2,3,ii) = (ps(2,1)*ps(1,3)-ps(1,1)*ps(2,3))/det
         invpolscale(3,1,ii) = (ps(2,1)*ps(3,2)-ps(3,1)*ps(2,2))/det
         invpolscale(3,2,ii) = (ps(3,1)*ps(1,2)-ps(1,1)*ps(3,2))/det
         invpolscale(3,3,ii) = (ps(1,1)*ps(2,2)-ps(2,1)*ps(1,2))/det
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine expolar0b  --  variable polarizability via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "expolar0b" calculates the variable polarizability due to exchange
c     polarization using a neighbor list
c
c
      subroutine expolar0b (polscale,invpolscale)
      use atoms
      use bound
      use couple
      use expol
      use mpole
      use neigh
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
      real*8 det
      real*8 sizi,sizk,sizik
      real*8 alphai,alphak
      real*8 springi,springk
      real*8 s2,ds2
      real*8 p33i, p33k
      real*8 kS2i(3,3)
      real*8 kS2k(3,3)
      real*8 taper
      real*8 ps(3,3)
      real*8, allocatable :: pscale(:)
      real*8 polscale(3,3,npole)
      real*8 invpolscale(3,3,npole)
      logical epli,eplk
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     set the switching function coefficients
c
      mode = 'REPULS'
      call switch (mode)
c
c     set polscale and invpolscale to the identity matrix
c
      do i = 1, npole
         polscale(1,1,i) = 1.0d0
         polscale(2,1,i) = 0.0d0
         polscale(3,1,i) = 0.0d0
         polscale(1,2,i) = 0.0d0
         polscale(2,2,i) = 1.0d0
         polscale(3,2,i) = 0.0d0
         polscale(1,3,i) = 0.0d0
         polscale(2,3,i) = 0.0d0
         polscale(3,3,i) = 1.0d0
         invpolscale(1,1,i) = 1.0d0
         invpolscale(2,1,i) = 0.0d0
         invpolscale(3,1,i) = 0.0d0
         invpolscale(1,2,i) = 0.0d0
         invpolscale(2,2,i) = 1.0d0
         invpolscale(3,2,i) = 0.0d0
         invpolscale(1,3,i) = 0.0d0
         invpolscale(2,3,i) = 0.0d0
         invpolscale(3,3,i) = 1.0d0
      end do
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
!$OMP& cut2,off2,c0,c1,c2,c3,c4,c5)
!$OMP& firstprivate(pscale)
!$OMP& shared (polscale,invpolscale)
!$OMP DO reduction(+:polscale,invpolscale) schedule(guided)
c
c     find the variable polarizability
c
      do ii = 1, npole
         i = ipole(ii)
         springi = kpep(ii)
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
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            eplk = lpep(kk)
            if (epli.or.eplk) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  springk = kpep(kk)
                  sizk = prepep(kk)
                  alphak = dmppep(kk)
                  sizik = sizi*sizk
                  call dampep (r,sizik,alphai,alphak,s2,ds2)
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     s2 = s2 * taper
                  end if
                  p33i = springi*s2*pscale(k)
                  p33k = springk*s2*pscale(k)
                  call exrotate(i,k,xr,yr,zr,p33i,p33k,kS2i,kS2k)
                  if (epli) then
                     do j = 1, 3
                        do m = 1, 3
                           polscale(j,m,ii) = polscale(j,m,ii)+kS2i(j,m)
                        end do
                     end do
                  end if
                  if (eplk) then
                     do j = 1, 3
                        do m = 1, 3
                           polscale(j,m,kk) = polscale(j,m,kk)+kS2k(j,m)
                        end do
                     end do
                  end if
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
c
c     invert polscale matrix
c
!$OMP DO schedule(guided)
      do ii = 1, npole
         do j = 1, 3
            do m = 1, 3
               ps(j,m) = polscale(j,m,ii)
            end do
         end do
         det = ps(1,1)*(ps(2,2)*ps(3,3) - ps(3,2)*ps(2,3))
     &       - ps(1,2)*(ps(2,1)*ps(3,3) - ps(2,3)*ps(3,1))
     &       + ps(1,3)*(ps(2,1)*ps(3,2) - ps(2,2)*ps(3,1))
         invpolscale(1,1,ii) = (ps(2,2)*ps(3,3)-ps(3,2)*ps(2,3))/det
         invpolscale(1,2,ii) = (ps(1,3)*ps(3,2)-ps(1,2)*ps(3,3))/det
         invpolscale(1,3,ii) = (ps(1,2)*ps(2,3)-ps(1,3)*ps(2,2))/det
         invpolscale(2,1,ii) = (ps(2,3)*ps(3,1)-ps(2,1)*ps(3,3))/det
         invpolscale(2,2,ii) = (ps(1,1)*ps(3,3)-ps(1,3)*ps(3,1))/det
         invpolscale(2,3,ii) = (ps(2,1)*ps(1,3)-ps(1,1)*ps(2,3))/det
         invpolscale(3,1,ii) = (ps(2,1)*ps(3,2)-ps(3,1)*ps(2,2))/det
         invpolscale(3,2,ii) = (ps(3,1)*ps(1,2)-ps(1,1)*ps(3,2))/det
         invpolscale(3,3,ii) = (ps(1,1)*ps(2,2)-ps(2,1)*ps(1,2))/det
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine exrotate  --  rotate polarizability matrix  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "exrotate" rotates and inverts the variable polarizability tensor
c     due to exchange polarization
c
c
      subroutine exrotate (i,k,xr,yr,zr,p33i,p33k,kS2i,kS2k)
      use atoms
      use math
      use mpole
      use polpot
      implicit none
      integer i,k
      integer j,l
      real*8 r,dot
      real*8 eps
      real*8 xr,yr,zr
      real*8 dx,dy,dz
      real*8 p33i,p33k
      real*8 ai(3,3)
      real*8 ak(3,3)
      real*8 kS2i(3,3)
      real*8 kS2k(3,3)
c
c
c     use the identity matrix as the default rotation matrix
c
      ai(1,1) = 1.0d0
      ai(2,1) = 0.0d0
      ai(3,1) = 0.0d0
      ai(1,2) = 0.0d0
      ai(2,2) = 1.0d0
      ai(3,2) = 0.0d0
      ai(1,3) = 0.0d0
      ai(2,3) = 0.0d0
      ai(3,3) = 1.0d0
c
c     get the rotation matrix elements
c
      dx = xr
      dy = yr
      dz = zr
      r = sqrt(dx*dx + dy*dy + dz*dz)
      ai(3,1) = dx / r
      ai(3,2) = dy / r
      ai(3,3) = dz / r
      dx = 1.0d0
      dy = 0.0d0
      dz = 0.0d0
      dot = ai(3,1)
      eps = 0.707d0
      if (abs(dot) .gt. eps) then
         dx = 0.0d0
         dy = 1.0d0
         dot = ai(3,2)
      end if
      dx = dx - dot*ai(3,1)
      dy = dy - dot*ai(3,2)
      dz = dz - dot*ai(3,3)
      r = sqrt(dx*dx + dy*dy + dz*dz)
c
c     ai is the rotation matrix R that takes a vector in
c     global frame to local frame
c
      ai(1,1) = dx / r
      ai(1,2) = dy / r
      ai(1,3) = dz / r
      ai(2,1) = ai(1,3)*ai(3,2) - ai(1,2)*ai(3,3)
      ai(2,2) = ai(1,1)*ai(3,3) - ai(1,3)*ai(3,1)
      ai(2,3) = ai(1,2)*ai(3,1) - ai(1,1)*ai(3,2)
      ak = ai
      ak(2,1) = -ak(2,1)
      ak(2,2) = -ak(2,2)
      ak(2,3) = -ak(2,3)
      ak(3,1) = -ak(3,1)
      ak(3,2) = -ak(3,2)
      ak(3,3) = -ak(3,3)
c
c     apply R^T from left and R from right to rotate kS2 matrix
c
      do j = 1, 3
         do l = 1, 3
            kS2i(j,l) = p33i*ai(3,j)*ai(3,l)
            kS2k(j,l) = p33k*ak(3,j)*ak(3,l)
         end do
      end do
      return
      end
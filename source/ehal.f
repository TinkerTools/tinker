c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ehal  --  buffered 14-7 van der Waals energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ehal" calculates the buffered 14-7 van der Waals energy
c
c
      subroutine ehal
      implicit none
      real*8 elrc
      include 'cutoff.i'
      include 'energi.i'
      include 'vdwpot.i'
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_lights) then
         call ehal0b
      else if (use_vlist) then
         call ehal0c
      else
         call ehal0a
      end if
c
c     apply long range van der Waals correction if desired
c
      if (use_vcorr) then
         call evcorr (elrc)
         ev = ev + elrc
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal0a  --  buffered 14-7 vdw via double loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal0a" calculates the buffered 14-7 van der Waals energy
c     using a pairwise double loop
c
c
      subroutine ehal0a
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer iv14(maxatm)
      real*8 e,rv,rv7
      real*8 eps,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      logical proceed,usei
c
c
c     zero out the van der Waals energy contribution
c
      ev = 0.0d0
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     apply any reduction factor to the atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find the van der Waals energy via double loop search
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         iv = ired(i)
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + ghal*rv7
                  tau = (dhal+1.0d0) / (rik + dhal*rv)
                  e = eps * rv7 * tau**7
     &                   * ((ghal+1.0d0)*rv7/rho-2.0d0)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  ev = ev + e
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               do j = 1, ncell
                  xr = xi - xred(k)
                  yr = yi - yred(k)
                  zr = zi - zred(k)
                  call imager (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
                  if (rik2 .le. off2) then
                     rik = sqrt(rik2)
                     rv = radmin(kt,it)
                     eps = epsilon(kt,it)
                     if (use_polymer) then
                        if (rik2 .le. polycut2) then
                           if (iv14(k) .eq. i) then
                              rv = radmin4(kt,it)
                              eps = epsilon4(kt,it)
                           end if
                           eps = eps * vscale(k)
                        end if
                     end if
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0d0)*rv7/rho-2.0d0)
c
c     use energy switching if near the cutoff distance
c
                     if (rik2 .gt. cut2) then
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        e = e * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy component;
c     interaction of an atom with its own image counts half
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ev = ev + e
                  end if
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ehal0b  --  buffered 14-7 vdw energy via lights  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ehal0b" calculates the buffered 14-7 van der Waals energy
c     using the method of lights
c
c
      subroutine ehal0b
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'light.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer kgy,kgz
      integer start,stop
      integer iv14(maxatm)
      real*8 e,rv,rv7
      real*8 eps,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      real*8 xsort(maxlight)
      real*8 ysort(maxlight)
      real*8 zsort(maxlight)
      logical proceed,usei
      logical prime,repeat
c
c
c     zero out the van der Waals energy contribution
c
      ev = 0.0d0
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     apply any reduction factor to the atomic coordinates
c
      do j = 1, nvdw
         i = ivdw(j)
         iv = ired(i)
         rdn = kred(i)
         xred(j) = rdn*(x(i)-x(iv)) + x(iv)
         yred(j) = rdn*(y(i)-y(iv)) + y(iv)
         zred(j) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      call lights (off,nvdw,xsort,ysort,zsort)
c
c     loop over all atoms computing the interactions
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         it = jvdw(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         usei = (use(i) .or. use(iv))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     loop over method of lights neighbors of current atom
c
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            k = ivdw(kk-((kk-1)/nvdw)*nvdw)
            kv = ired(k)
            prime = (kk .le. nvdw)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               if (use_bounds) then
                  if (abs(xr) .gt. xcell2)  xr = xr - sign(xcell,xr)
                  if (abs(yr) .gt. ycell2)  yr = yr - sign(ycell,yr)
                  if (abs(zr) .gt. zcell2)  zr = zr - sign(zcell,zr)
                  if (monoclinic) then
                     xr = xr + zr*beta_cos
                     zr = zr * beta_sin
                  else if (triclinic) then
                     xr = xr + yr*gamma_cos + zr*beta_cos
                     yr = yr*gamma_sin + zr*beta_term
                     zr = zr * gamma_term
                  end if
               end if
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (prime) then
                     if (iv14(k) .eq. i) then
                        rv = radmin4(kt,it)
                        eps = epsilon4(kt,it)
                     end if
                     eps = eps * vscale(k)
                  end if
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + ghal*rv7
                  tau = (dhal+1.0d0) / (rik + dhal*rv)
                  e = eps * rv7 * tau**7
     &                   * ((ghal+1.0d0)*rv7/rho-2.0d0)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  ev = ev + e
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
         end if
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ehal0c  --  buffered 14-7 vdw energy via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ehal0c" calculates the buffered 14-7 van der Waals energy
c     using a pairwise neighbor list
c
c
      subroutine ehal0c
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'neigh.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer iv14(maxatm)
      real*8 e,eps,rdn,evt
      real*8 fgrp,rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8 xred(maxatm)
      real*8 yred(maxatm)
      real*8 zred(maxatm)
      real*8 vscale(maxatm)
      logical proceed,usei
c
c
c     zero out the van der Waals energy contribution
c
      ev = 0.0d0
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c
      call switch ('VDW')
c
c     apply any reduction factor to the atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     transfer global to local copies for OpenMP calculation
c
      evt = ev
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nvdw,ivdw,ired,kred,
!$OMP& jvdw,xred,yred,zred,use,nvlst,vlst,n12,n13,n14,n15,
!$OMP& i12,i13,i14,i15,v2scale,v3scale,v4scale,v5scale,
!$OMP& use_group,fgrp,off2,radmin,epsilon,radmin4,epsilon4,
!$OMP& ghal,dhal,cut2,c0,c1,c2,c3,c4,c5)
!$OMP& firstprivate(vscale)
!$OMP& shared(evt)
!$OMP DO reduction(+:evt) schedule(dynamic)
c
c     find the van der Waals energy via neighbor list search
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = 1, nvlst(ii)
            k = ivdw(vlst(kk,ii))
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rik = sqrt(rik2)
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  rv7 = rv**7
                  rik7 = rik**7
                  rho = rik7 + ghal*rv7
                  tau = (dhal+1.0d0) / (rik + dhal*rv)
                  e = eps * rv7 * tau**7
     &                   * ((ghal+1.0d0)*rv7/rho-2.0d0)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c
c     increment the overall van der Waals energy components
c
                  evt = evt + e
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     transfer local to global copies for OpenMP calculation
c
      ev = evt
      return
      end

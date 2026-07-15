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
      use dlmda
      use mutant
      use energi
      use limits
      use vdwpot
      implicit none
      real*8 elrc
      character*6 mode
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_evdt) then
         if (use_rel) then
            call ehal0dr
         else
            call ehal0d
         end if
      else
         if (use_lights) then
            call ehal0b
         else if (use_vlist) then
            call ehal0c
         else
            call ehal0a
         end if
c
c     apply the long range van der Waals correction if used
c
         if (use_vcorr) then
            mode = 'VDW'
            call evcorr (mode,elrc)
            ev = ev + elrc
         end if
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
      use atomid
      use atoms
      use bound
      use cell
      use couple
      use energi
      use group
      use mutant
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,k
      integer ii,it,iv
      integer kk,kt,kv
      integer, allocatable :: iv14(:)
      real*8 e,rv,rv7
      real*8 eps,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 scal,t1,t2
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical muti,mutk,mutik
      character*6 mode
c
c
c     zero out the van der Waals energy contribution
c
      ev = 0.0d0
      if (nvdw .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         iv14(i) = 0
         vscale(i) = 1.0d0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
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
         it = jvdw(i)
         iv = ired(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
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
            kt = jvdw(k)
            kv = ired(k)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
            if (use_subsys .and. proceed) then
               if (.not.subon(i) .or. .not.subon(k))  proceed = .false.
            end if
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
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
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = .false.
                  if ((muti .or. mutk) .and. .not.use_subsys) then
                     if (vcouple .eq. 1) then
                        mutik = .true.
                     else if (.not.muti .or. .not.mutk) then
                        mutik = .true.
                     end if
                  end if
c
c     get interaction energy, via soft core lambda scaling as needed
c
                  if (mutik) then
                     rho = rik / rv
                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0d0-vlambda)**2
                     t1 = (1.0d0+dhal)**7 / (scal+(rho+dhal)**7)
                     t2 = (1.0d0+ghal) / (scal+rho**7+ghal)
                     e = eps * t1 * (t2-2.0d0)
                  else
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0d0)*rv7/rho-2.0d0)
                  end if
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
c     reset exclusion coefficients for connected atoms
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
         it = jvdw(i)
         iv = ired(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
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
            kt = jvdw(k)
            kv = ired(k)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
            if (use_subsys .and. proceed) then
               if (.not.subon(i) .or. .not.subon(k))  proceed = .false.
            end if
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 2, ncell
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
c
c     set use of lambda scaling for decoupling or annihilation
c
                     mutik = .false.
                     if ((muti .or. mutk) .and. .not.use_subsys) then
                        if (vcouple .eq. 1) then
                           mutik = .true.
                        else if (.not.muti .or. .not.mutk) then
                           mutik = .true.
                        end if
                     end if
c
c     get interaction energy, via soft core lambda scaling as needed
c
                     if (mutik) then
                        rho = rik / rv
                        eps = eps * vlambda**scexp
                        scal = scalpha * (1.0d0-vlambda)**2
                        t1 = (1.0d0+dhal)**7 / (scal+(rho+dhal)**7)
                        t2 = (1.0d0+ghal) / (scal+rho**7+ghal)
                        e = eps * t1 * (t2-2.0d0)
                     else
                        rv7 = rv**7
                        rik7 = rik**7
                        rho = rik7 + ghal*rv7
                        tau = (dhal+1.0d0) / (rik + dhal*rv)
                        e = eps * rv7 * tau**7
     &                         * ((ghal+1.0d0)*rv7/rho-2.0d0)
                     end if
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
c     reset exclusion coefficients for connected atoms
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
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (vscale)
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
      use atomid
      use atoms
      use bound
      use boxes
      use cell
      use couple
      use energi
      use group
      use light
      use mutant
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,k
      integer ii,it,iv
      integer kk,kt,kv
      integer kgy,kgz
      integer start,stop
      integer, allocatable :: iv14(:)
      real*8 e,rv,rv7
      real*8 eps,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 scal,t1,t2
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical unique,repeat
      logical muti,mutk,mutik
      character*6 mode
c
c
c     zero out the van der Waals energy contribution
c
      ev = 0.0d0
      if (nvdw .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (vscale(n))
      allocate (xsort(8*n))
      allocate (ysort(8*n))
      allocate (zsort(8*n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         iv14(i) = 0
         vscale(i) = 1.0d0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
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
      unique = .true.
      call lights (off,nvdw,xsort,ysort,zsort,unique)
c
c     loop over all atoms computing the interactions
c
      do ii = 1, nvdw
         i = ivdw(ii)
         it = jvdw(i)
         iv = ired(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         usei = (use(i) .or. use(iv))
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
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
            kt = jvdw(k)
            kv = ired(k)
            mutk = mut(k)
            prime = (kk .le. nvdw)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
            if (use_subsys .and. proceed) then
               if (.not.subon(i) .or. .not.subon(k))  proceed = .false.
            end if
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
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
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = .false.
                  if ((muti .or. mutk) .and. .not.use_subsys) then
                     if (vcouple .eq. 1) then
                        mutik = .true.
                     else if (.not.muti .or. .not.mutk) then
                        mutik = .true.
                     end if
                  end if
c
c     get interaction energy, via soft core lambda scaling as needed
c
                  if (mutik) then
                     rho = rik / rv
                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0d0-vlambda)**2
                     t1 = (1.0d0+dhal)**7 / (scal+(rho+dhal)**7)
                     t2 = (1.0d0+ghal) / (scal+rho**7+ghal)
                     e = eps * t1 * (t2-2.0d0)
                  else
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0d0)*rv7/rho-2.0d0)
                  end if
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
c     reset exclusion coefficients for connected atoms
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
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (vscale)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
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
      use atomid
      use atoms
      use bound
      use couple
      use energi
      use group
      use mutant
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      implicit none
      integer i,j,k
      integer ii,it,iv
      integer kk,kt,kv
      integer, allocatable :: iv14(:)
      real*8 e,eps,rdn
      real*8 fgrp,rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rho,tau,taper
      real*8 scal,t1,t2
      real*8 rik,rik2,rik3
      real*8 rik4,rik5,rik7
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical muti,mutk,mutik
      character*6 mode
c
c
c     zero out the van der Waals energy contribution
c
      ev = 0.0d0
      if (nvdw .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         iv14(i) = 0
         vscale(i) = 1.0d0
      end do
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
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
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nvdw,ivdw,jvdw,ired,
!$OMP& subon,use_subsys,
!$OMP& xred,yred,zred,use,nvlst,vlst,n12,n13,n14,n15,i12,i13,
!$OMP& i14,i15,v2scale,v3scale,v4scale,v5scale,use_group,
!$OMP& off2,radmin,epsilon,radmin4,epsilon4,ghal,dhal,vcouple,
!$OMP& vlambda,mut,scexp,scalpha,cut2,c0,c1,c2,c3,c4,c5)
!$OMP& firstprivate(vscale,iv14) shared(ev)
!$OMP DO reduction(+:ev)
c
c     find the van der Waals energy via neighbor list search
c
      do ii = 1, nvdw
         i = ivdw(ii)
         it = jvdw(i)
         iv = ired(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
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
         do kk = 1, nvlst(i)
            k = vlst(kk,i)
            kt = jvdw(k)
            kv = ired(k)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
            if (use_subsys .and. proceed) then
               if (.not.subon(i) .or. .not.subon(k))  proceed = .false.
            end if
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
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
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = .false.
                  if ((muti .or. mutk) .and. .not.use_subsys) then
                     if (vcouple .eq. 1) then
                        mutik = .true.
                     else if (.not.muti .or. .not.mutk) then
                        mutik = .true.
                     end if
                  end if
c
c     get interaction energy, via soft core lambda scaling as needed
c
                  if (mutik) then
                     rho = rik / rv
                     eps = eps * vlambda**scexp
                     scal = scalpha * (1.0d0-vlambda)**2
                     t1 = (1.0d0+dhal)**7 / (scal+(rho+dhal)**7)
                     t2 = (1.0d0+ghal) / (scal+rho**7+ghal)
                     e = eps * t1 * (t2-2.0d0)
                  else
                     rv7 = rv**7
                     rik7 = rik**7
                     rho = rik7 + ghal*rv7
                     tau = (dhal+1.0d0) / (rik + dhal*rv)
                     e = eps * rv7 * tau**7
     &                      * ((ghal+1.0d0)*rv7/rho-2.0d0)
                  end if
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
c     reset exclusion coefficients for connected atoms
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (vscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ehal0d  --  dual topology buffered 14-7 energy   ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ehal0d" calculates the buffered 14-7 van der Waals energy with
c     the dual topology method, in which the fully coupled (vlambda=1)
c     and fully decoupled (vlambda=0) states are each evaluated in full
c     and combined by a power law interpolation in vlambda
c
c     note the long range correction depends on vlambda, so it must be
c     evaluated separately for each end state and interpolated along
c     with the pairwise energy
c
c
      subroutine ehal0d
      use dlmda
      use energi
      use mutant
      implicit none
      real*8 ev1,ev0
      real*8 vlambdaorig
      real*8 weight1
c
c
c     compute energy of the fully coupled vlambda = 1 state
c
      vlambdaorig = vlambda
      vlambda = 1.0d0
      call ehal0sub
      ev1 = ev
c
c     compute energy of the fully decoupled vlambda = 0 state
c
      vlambda = 0.0d0
      call ehal0sub
      ev0 = ev
c
c     restore the original vlambda value
c
      vlambda = vlambdaorig
c
c     interpolate the dual topology energy
c
      weight1 = vlambda**evdtexp
      ev = weight1*ev1 + (1.0d0-weight1)*ev0
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ehal0sub  --  subsystem buffered 14-7 energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ehal0sub" evaluates the buffered 14-7 van der Waals energy for
c     the atom subsystem currently flagged by "subon", using the same
c     standard routine selection as "ehal0d"; the long range correction
c     follows the active subsystem when enabled
c
c
      subroutine ehal0sub
      use energi
      use limits
      use vdwpot
      implicit none
      real*8 elrc
      character*6 mode
c
c
      if (use_lights) then
         call ehal0b
      else if (use_vlist) then
         call ehal0c
      else
         call ehal0a
      end if
      if (use_vcorr) then
         mode = 'VDW'
         call evcorr (mode,elrc)
         ev = ev + elrc
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ehal0dr  --  relative dual topology 14-7 energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ehal0dr" calculates the buffered 14-7 van der Waals energy for a
c     two-ligand relative dual topology calculation by combining four
c     subsystem energies, E1 = E(A+env) + E(B) and E0 = E(B+env) + E(A),
c     so that intramolecular van der Waals is preserved, only the ligand
c     environment coupling is scaled, and the two ligands never interact
c
c
      subroutine ehal0dr
      use dlmda
      use energi
      use mutant
      implicit none
      real*8 evae,evbe
      real*8 eva,evb
      real*8 ev1,ev0
      real*8 weight1,weight0
c
c
      call submask (.true.,.false.,.true.)
      call ehal0sub
      evae = ev
      call submask (.false.,.true.,.true.)
      call ehal0sub
      evbe = ev
      call submask (.true.,.false.,.false.)
      call ehal0sub
      eva = ev
      call submask (.false.,.true.,.false.)
      call ehal0sub
      evb = ev
      call submask (.true.,.true.,.true.)
      weight1 = vlambda**evdtexp
      weight0 = 1.0d0 - weight1
      ev1 = evae + evb
      ev0 = evbe + eva
      ev = weight1*ev1 + weight0*ev0
      return
      end

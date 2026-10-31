c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine elj1  --  Lennard-Jones energy & derivatives  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj1" calculates the Lennard-Jones 6-12 van der Waals energy
c     and its first derivatives with respect to Cartesian coordinates
c
c
      subroutine elj1
      use dlmda
      use energi
      use limits
      use vdwpot
      use virial
      use warp
      implicit none
      integer i,j
      real*8 elrc,vlrc
      real*8 delrc,dvlrc
      real*8 d2elrc
      logical dodl,dod2
      character*6 mode
c
c
c     take lambda derivatives as "vdlflags" sets for the van der
c     Waals lambda
c
      call vdlflags (dodl,dod2)
c
c     zero out the van der Waals virial, which is kept apart from
c     the main virial for the lambda derivatives
c
      do i = 1, 3
         do j = 1, 3
            evvir(j,i) = 0.0d0
         end do
      end do
c
c     choose the method for summing over pairwise interactions
c
      if (use_stophat) then
         call elj1e
      else if (use_smooth) then
         call elj1d
      else if (use_vlist) then
         call elj1c
      else if (use_lights) then
         call elj1b
      else
         call elj1a
      end if
c
c     apply the long range van der Waals correction if used,
c     along with its lambda derivative
c
      if (use_vcorr) then
         mode = 'VDW'
         call evcorr1 (mode,elrc,vlrc,delrc,dvlrc,d2elrc)
         ev = ev + elrc
         if (dodl)  devdl = devdl + delrc
         if (dod2)  d2evdl2 = d2evdl2 + d2elrc
         do i = 1, 3
            evvir(i,i) = evvir(i,i) + vlrc
            if (dod2)  devvirdl(i,i) = devvirdl(i,i) + dvlrc
         end do
      end if
c
c     add the van der Waals virial to main virial
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + evvir(j,i)
         end do
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine elj1a  --  double loop Lennard-Jones vdw derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "elj1a" calculates the Lennard-Jones 6-12 van der Waals energy
c     and its first derivatives using a pairwise double loop
c
c
      subroutine elj1a
      use atomid
      use atoms
      use bound
      use cell
      use couple
      use deriv
      use dlmda
      use energi
      use group
      use mutant
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      implicit none
      integer i,j,k
      integer ii,it,iv
      integer kk,kt,kv
      integer ig
      integer, allocatable :: iv14(:)
      real*8 e,de,p6,p12
      real*8 eps
      real*8 vlsc(2),dvlsc(2),d2vlsc(2)
      real*8 vshf(2),dvshf(2),d2vshf(2)
      real*8 rv,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 taper,dtaper
      real*8 dlambda,dlambda2
      real*8 dlde
      real*8 dldedx,dldedy,dldedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      logical mutik
      logical dodl,dod2,mutdl,mutd2
      character*6 mode
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
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
c     take lambda derivatives as "vdlflags" sets for the van der
c     Waals lambda
c
      call vdlflags (dodl,dod2)
c
c     set the soft core lambda terms for each ligand group, along
c     with their lambda derivatives as needed
c
      if (dodl) then
         call grpsoftd (scexplj,scalphalj,scexpljm,scexpljp,vlsc,dvlsc,
     &                  d2vlsc,vshf,dvshf,d2vshf)
      else
         call grpsoft (scexplj,scalphalj,scexpljm,scexpljp,vlsc,vshf)
      end if
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
c     find van der Waals energy and derivatives via double loop
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         it = jvdw(i)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
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
            ig = mutcls(mutg(k),mutg(i))
c
c     the two ligands of a relative dual topology never interact
c
            proceed = (ig .ge. 0)
            if (use_group .and. proceed)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
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
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  rik = sqrt(rik2)
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (ig .gt. 0)
                  mutdl = mutik .and. dodl
                  mutd2 = mutik .and. dod2
c
c     get interaction energy, via soft core lambda scaling as needed
c
                  if (mutdl) then
                     call eljsc4 (rik,rik2,rv,eps,ig,vlsc,dvlsc,
     &                            d2vlsc,vshf,dvshf,d2vshf,e,de,
     &                            dlambda,dlambda2,dlde)
                  else if (mutik) then
                     call eljsc (rik,rik2,rv,eps,ig,vlsc,vshf,
     &                           .true.,e,de)
                  else
                     p6 = rv**6 / rik2**3
                     p12 = p6 * p6
                     e = eps * (p12 - 2.0d0*p6)
                     de = eps * (p12-p6) * (-12.0d0/rik)
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
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                     if (mutd2) then
                        dlde = dlambda * dtaper + dlde * taper
                        dlambda2 = dlambda2 * taper
                     end if
                     if (mutdl)  dlambda = dlambda * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     if (mutdl)  dlambda = dlambda * fgrp
                     if (mutd2) then
                        dlambda2 = dlambda2 * fgrp
                        dlde = dlde * fgrp
                     end if
                  end if
c
c     find the chain rule terms for derivative components; atoms at
c     the same point have no direction and take zero force
c
                  if (rik .ne. 0.0d0) then
                     de = de / rik
                  else
                     de = 0.0d0
                  end if
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  ev = ev + e
                  if (mutd2) then
                     if (rik .ne. 0.0d0) then
                        dlde = dlde / rik
                     else
                        dlde = 0.0d0
                     end if
                     dldedx = dlde * xr
                     dldedy = dlde * yr
                     dldedz = dlde * zr
                  end if
                  if (i .eq. iv) then
                     dev(1,i) = dev(1,i) + dedx
                     dev(2,i) = dev(2,i) + dedy
                     dev(3,i) = dev(3,i) + dedz
                     if (mutd2) then
                        dfvdl(1,i) = dfvdl(1,i) + dldedx
                        dfvdl(2,i) = dfvdl(2,i) + dldedy
                        dfvdl(3,i) = dfvdl(3,i) + dldedz
                     end if
                  else
                     dev(1,i) = dev(1,i) + dedx*redi
                     dev(2,i) = dev(2,i) + dedy*redi
                     dev(3,i) = dev(3,i) + dedz*redi
                     dev(1,iv) = dev(1,iv) + dedx*rediv
                     dev(2,iv) = dev(2,iv) + dedy*rediv
                     dev(3,iv) = dev(3,iv) + dedz*rediv
                     if (mutd2) then
                        dfvdl(1,i) = dfvdl(1,i) + dldedx*redi
                        dfvdl(2,i) = dfvdl(2,i) + dldedy*redi
                        dfvdl(3,i) = dfvdl(3,i) + dldedz*redi
                        dfvdl(1,iv) = dfvdl(1,iv) + dldedx*rediv
                        dfvdl(2,iv) = dfvdl(2,iv) + dldedy*rediv
                        dfvdl(3,iv) = dfvdl(3,iv) + dldedz*rediv
                     end if
                  end if
                  if (k .eq. kv) then
                     dev(1,k) = dev(1,k) - dedx
                     dev(2,k) = dev(2,k) - dedy
                     dev(3,k) = dev(3,k) - dedz
                     if (mutd2) then
                        dfvdl(1,k) = dfvdl(1,k) - dldedx
                        dfvdl(2,k) = dfvdl(2,k) - dldedy
                        dfvdl(3,k) = dfvdl(3,k) - dldedz
                     end if
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     dev(1,k) = dev(1,k) - dedx*redk
                     dev(2,k) = dev(2,k) - dedy*redk
                     dev(3,k) = dev(3,k) - dedz*redk
                     dev(1,kv) = dev(1,kv) - dedx*redkv
                     dev(2,kv) = dev(2,kv) - dedy*redkv
                     dev(3,kv) = dev(3,kv) - dedz*redkv
                     if (mutd2) then
                        dfvdl(1,k) = dfvdl(1,k) - dldedx*redk
                        dfvdl(2,k) = dfvdl(2,k) - dldedy*redk
                        dfvdl(3,k) = dfvdl(3,k) - dldedz*redk
                        dfvdl(1,kv) = dfvdl(1,kv) - dldedx*redkv
                        dfvdl(2,kv) = dfvdl(2,kv) - dldedy*redkv
                        dfvdl(3,kv) = dfvdl(3,kv) - dldedz*redkv
                     end if
                  end if
                  if (mutdl)  devdl = devdl + dlambda
                  if (mutd2)  d2evdl2 = d2evdl2 + dlambda2
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  evvir(1,1) = evvir(1,1) + vxx
                  evvir(2,1) = evvir(2,1) + vyx
                  evvir(3,1) = evvir(3,1) + vzx
                  evvir(1,2) = evvir(1,2) + vyx
                  evvir(2,2) = evvir(2,2) + vyy
                  evvir(3,2) = evvir(3,2) + vzy
                  evvir(1,3) = evvir(1,3) + vzx
                  evvir(2,3) = evvir(2,3) + vzy
                  evvir(3,3) = evvir(3,3) + vzz
                  if (mutd2) then
                     call dlcvir (xr,yr,zr,dldedx,dldedy,dldedz,
     &                            devvirdl)
                  end if
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
         redi = kred(i)
         rediv = 1.0d0 - redi
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
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
            ig = mutcls(mutg(k),mutg(i))
c
c     the two ligands of a relative dual topology never interact
c
            proceed = (ig .ge. 0)
            if (use_group .and. proceed)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
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
                     rik = sqrt(rik2)
c
c     set use of lambda scaling for decoupling or annihilation
c
                     mutik = (ig .gt. 0)
                     mutdl = mutik .and. dodl
                     mutd2 = mutik .and. dod2
c
c     get interaction energy, via soft core lambda scaling as needed
c
                     if (mutdl) then
                        call eljsc4 (rik,rik2,rv,eps,ig,vlsc,dvlsc,
     &                               d2vlsc,vshf,dvshf,d2vshf,e,de,
     &                               dlambda,dlambda2,dlde)
                     else if (mutik) then
                        call eljsc (rik,rik2,rv,eps,ig,vlsc,vshf,
     &                              .true.,e,de)
                     else
                        p6 = rv**6 / rik2**3
                        p12 = p6 * p6
                        e = eps * (p12 - 2.0d0*p6)
                        de = eps * (p12-p6) * (-12.0d0/rik)
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
                        dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                              + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                        de = e*dtaper + de*taper
                        e = e * taper
                        if (mutd2) then
                           dlde = dlambda * dtaper + dlde * taper
                           dlambda2 = dlambda2 * taper
                        end if
                        if (mutdl)  dlambda = dlambda * taper
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                        if (mutdl)  dlambda = dlambda * fgrp
                        if (mutd2) then
                           dlambda2 = dlambda2 * fgrp
                           dlde = dlde * fgrp
                        end if
                     end if
c
c     find the chain rule terms for derivative components; atoms at
c     the same point have no direction and take zero force
c
                     if (rik .ne. 0.0d0) then
                        de = de / rik
                     else
                        de = 0.0d0
                     end if
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                     if (i .eq. k) then
                        e = 0.5d0 * e
                        if (mutdl)  dlambda = 0.5d0 * dlambda
                        if (mutd2)  dlambda2 = 0.5d0 * dlambda2
                     end if
                     ev = ev + e
                     if (mutd2) then
                        if (rik .ne. 0.0d0) then
                           dlde = dlde / rik
                        else
                           dlde = 0.0d0
                        end if
                        dldedx = dlde * xr
                        dldedy = dlde * yr
                        dldedz = dlde * zr
                     end if
                     if (i .eq. iv) then
                        dev(1,i) = dev(1,i) + dedx
                        dev(2,i) = dev(2,i) + dedy
                        dev(3,i) = dev(3,i) + dedz
                        if (mutd2) then
                           dfvdl(1,i) = dfvdl(1,i) + dldedx
                           dfvdl(2,i) = dfvdl(2,i) + dldedy
                           dfvdl(3,i) = dfvdl(3,i) + dldedz
                        end if
                     else
                        dev(1,i) = dev(1,i) + dedx*redi
                        dev(2,i) = dev(2,i) + dedy*redi
                        dev(3,i) = dev(3,i) + dedz*redi
                        dev(1,iv) = dev(1,iv) + dedx*rediv
                        dev(2,iv) = dev(2,iv) + dedy*rediv
                        dev(3,iv) = dev(3,iv) + dedz*rediv
                        if (mutd2) then
                           dfvdl(1,i) = dfvdl(1,i) + dldedx*redi
                           dfvdl(2,i) = dfvdl(2,i) + dldedy*redi
                           dfvdl(3,i) = dfvdl(3,i) + dldedz*redi
                           dfvdl(1,iv) = dfvdl(1,iv) + dldedx*rediv
                           dfvdl(2,iv) = dfvdl(2,iv) + dldedy*rediv
                           dfvdl(3,iv) = dfvdl(3,iv) + dldedz*rediv
                        end if
                     end if
                     if (i .ne. k) then
                        if (k .eq. kv) then
                           dev(1,k) = dev(1,k) - dedx
                           dev(2,k) = dev(2,k) - dedy
                           dev(3,k) = dev(3,k) - dedz
                           if (mutd2) then
                              dfvdl(1,k) = dfvdl(1,k) - dldedx
                              dfvdl(2,k) = dfvdl(2,k) - dldedy
                              dfvdl(3,k) = dfvdl(3,k) - dldedz
                           end if
                        else
                           redk = kred(k)
                           redkv = 1.0d0 - redk
                           dev(1,k) = dev(1,k) - dedx*redk
                           dev(2,k) = dev(2,k) - dedy*redk
                           dev(3,k) = dev(3,k) - dedz*redk
                           dev(1,kv) = dev(1,kv) - dedx*redkv
                           dev(2,kv) = dev(2,kv) - dedy*redkv
                           dev(3,kv) = dev(3,kv) - dedz*redkv
                           if (mutd2) then
                              dfvdl(1,k) = dfvdl(1,k) - dldedx*redk
                              dfvdl(2,k) = dfvdl(2,k) - dldedy*redk
                              dfvdl(3,k) = dfvdl(3,k) - dldedz*redk
                              dfvdl(1,kv) = dfvdl(1,kv) - dldedx*redkv
                              dfvdl(2,kv) = dfvdl(2,kv) - dldedy*redkv
                              dfvdl(3,kv) = dfvdl(3,kv) - dldedz*redkv
                           end if
                        end if
                     end if
                     if (mutdl)  devdl = devdl + dlambda
                     if (mutd2)  d2evdl2 = d2evdl2 + dlambda2
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     evvir(1,1) = evvir(1,1) + vxx
                     evvir(2,1) = evvir(2,1) + vyx
                     evvir(3,1) = evvir(3,1) + vzx
                     evvir(1,2) = evvir(1,2) + vyx
                     evvir(2,2) = evvir(2,2) + vyy
                     evvir(3,2) = evvir(3,2) + vzy
                     evvir(1,3) = evvir(1,3) + vzx
                     evvir(2,3) = evvir(2,3) + vzy
                     evvir(3,3) = evvir(3,3) + vzz
                     if (mutd2) then
                        call dlcvir (xr,yr,zr,dldedx,dldedy,dldedz,
     &                               devvirdl)
                     end if
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine elj1b  --  Lennard-Jones vdw derivs via lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "elj1b" calculates the Lennard-Jones 6-12 van der Waals energy
c     and its first derivatives using the method of lights
c
c
      subroutine elj1b
      use atomid
      use atoms
      use bound
      use boxes
      use cell
      use couple
      use deriv
      use dlmda
      use energi
      use group
      use light
      use mutant
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      implicit none
      integer i,j,k
      integer ii,it,iv
      integer kk,kt,kv
      integer kgy,kgz
      integer start,stop
      integer ig
      integer, allocatable :: iv14(:)
      real*8 e,de,p6,p12
      real*8 eps
      real*8 vlsc(2),dvlsc(2),d2vlsc(2)
      real*8 vshf(2),dvshf(2),d2vshf(2)
      real*8 rv,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 taper,dtaper
      real*8 dlambda,dlambda2
      real*8 dlde
      real*8 dldedx,dldedy,dldedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical unique,repeat
      logical mutik
      logical dodl,dod2,mutdl,mutd2
      character*6 mode
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
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
c     take lambda derivatives as "vdlflags" sets for the van der
c     Waals lambda
c
      call vdlflags (dodl,dod2)
c
c     set the soft core lambda terms for each ligand group, along
c     with their lambda derivatives as needed
c
      if (dodl) then
         call grpsoftd (scexplj,scalphalj,scexpljm,scexpljp,vlsc,dvlsc,
     &                  d2vlsc,vshf,dvshf,d2vshf)
      else
         call grpsoft (scexplj,scalphalj,scexpljm,scexpljp,vlsc,vshf)
      end if
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
         redi = kred(i)
         rediv = 1.0d0 - redi
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         usei = (use(i) .or. use(iv))
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
            ig = mutcls(mutg(k),mutg(i))
            prime = (kk .le. nvdw)
c
c     decide whether to compute the current interaction
c
c
c     the two ligands of a relative dual topology never interact
c
            proceed = (ig .ge. 0)
            if (use_group .and. proceed)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
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
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (prime) then
                     if (iv14(k) .eq. i) then
                        rv = radmin4(kt,it)
                        eps = epsilon4(kt,it)
                     end if
                     eps = eps * vscale(k)
                  end if
                  rik = sqrt(rik2)
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (ig .gt. 0)
                  mutdl = mutik .and. dodl
                  mutd2 = mutik .and. dod2
c
c     get interaction energy, via soft core lambda scaling as needed
c
                  if (mutdl) then
                     call eljsc4 (rik,rik2,rv,eps,ig,vlsc,dvlsc,
     &                            d2vlsc,vshf,dvshf,d2vshf,e,de,
     &                            dlambda,dlambda2,dlde)
                  else if (mutik) then
                     call eljsc (rik,rik2,rv,eps,ig,vlsc,vshf,
     &                           .true.,e,de)
                  else
                     p6 = rv**6 / rik2**3
                     p12 = p6 * p6
                     e = eps * (p12 - 2.0d0*p6)
                     de = eps * (p12-p6) * (-12.0d0/rik)
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
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                     if (mutd2) then
                        dlde = dlambda * dtaper + dlde * taper
                        dlambda2 = dlambda2 * taper
                     end if
                     if (mutdl)  dlambda = dlambda * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     if (mutdl)  dlambda = dlambda * fgrp
                     if (mutd2) then
                        dlambda2 = dlambda2 * fgrp
                        dlde = dlde * fgrp
                     end if
                  end if
c
c     find the chain rule terms for derivative components; atoms at
c     the same point have no direction and take zero force
c
                  if (rik .ne. 0.0d0) then
                     de = de / rik
                  else
                     de = 0.0d0
                  end if
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  ev = ev + e
                  if (mutd2) then
                     if (rik .ne. 0.0d0) then
                        dlde = dlde / rik
                     else
                        dlde = 0.0d0
                     end if
                     dldedx = dlde * xr
                     dldedy = dlde * yr
                     dldedz = dlde * zr
                  end if
                  if (i .eq. iv) then
                     dev(1,i) = dev(1,i) + dedx
                     dev(2,i) = dev(2,i) + dedy
                     dev(3,i) = dev(3,i) + dedz
                     if (mutd2) then
                        dfvdl(1,i) = dfvdl(1,i) + dldedx
                        dfvdl(2,i) = dfvdl(2,i) + dldedy
                        dfvdl(3,i) = dfvdl(3,i) + dldedz
                     end if
                  else
                     dev(1,i) = dev(1,i) + dedx*redi
                     dev(2,i) = dev(2,i) + dedy*redi
                     dev(3,i) = dev(3,i) + dedz*redi
                     dev(1,iv) = dev(1,iv) + dedx*rediv
                     dev(2,iv) = dev(2,iv) + dedy*rediv
                     dev(3,iv) = dev(3,iv) + dedz*rediv
                     if (mutd2) then
                        dfvdl(1,i) = dfvdl(1,i) + dldedx*redi
                        dfvdl(2,i) = dfvdl(2,i) + dldedy*redi
                        dfvdl(3,i) = dfvdl(3,i) + dldedz*redi
                        dfvdl(1,iv) = dfvdl(1,iv) + dldedx*rediv
                        dfvdl(2,iv) = dfvdl(2,iv) + dldedy*rediv
                        dfvdl(3,iv) = dfvdl(3,iv) + dldedz*rediv
                     end if
                  end if
                  if (k .eq. kv) then
                     dev(1,k) = dev(1,k) - dedx
                     dev(2,k) = dev(2,k) - dedy
                     dev(3,k) = dev(3,k) - dedz
                     if (mutd2) then
                        dfvdl(1,k) = dfvdl(1,k) - dldedx
                        dfvdl(2,k) = dfvdl(2,k) - dldedy
                        dfvdl(3,k) = dfvdl(3,k) - dldedz
                     end if
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     dev(1,k) = dev(1,k) - dedx*redk
                     dev(2,k) = dev(2,k) - dedy*redk
                     dev(3,k) = dev(3,k) - dedz*redk
                     dev(1,kv) = dev(1,kv) - dedx*redkv
                     dev(2,kv) = dev(2,kv) - dedy*redkv
                     dev(3,kv) = dev(3,kv) - dedz*redkv
                     if (mutd2) then
                        dfvdl(1,k) = dfvdl(1,k) - dldedx*redk
                        dfvdl(2,k) = dfvdl(2,k) - dldedy*redk
                        dfvdl(3,k) = dfvdl(3,k) - dldedz*redk
                        dfvdl(1,kv) = dfvdl(1,kv) - dldedx*redkv
                        dfvdl(2,kv) = dfvdl(2,kv) - dldedy*redkv
                        dfvdl(3,kv) = dfvdl(3,kv) - dldedz*redkv
                     end if
                  end if
                  if (mutdl)  devdl = devdl + dlambda
                  if (mutd2)  d2evdl2 = d2evdl2 + dlambda2
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  evvir(1,1) = evvir(1,1) + vxx
                  evvir(2,1) = evvir(2,1) + vyx
                  evvir(3,1) = evvir(3,1) + vzx
                  evvir(1,2) = evvir(1,2) + vyx
                  evvir(2,2) = evvir(2,2) + vyy
                  evvir(3,2) = evvir(3,2) + vzy
                  evvir(1,3) = evvir(1,3) + vzx
                  evvir(2,3) = evvir(2,3) + vzy
                  evvir(3,3) = evvir(3,3) + vzz
                  if (mutd2) then
                     call dlcvir (xr,yr,zr,dldedx,dldedy,dldedz,
     &                            devvirdl)
                  end if
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine elj1c  --  Lennard-Jones vdw derivs via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "elj1c" calculates the Lennard-Jones 12-6 van der Waals energy
c     and its first derivatives using a pairwise neighbor list
c
c
      subroutine elj1c
      use atomid
      use atoms
      use bound
      use couple
      use deriv
      use dlmda
      use energi
      use group
      use mutant
      use neigh
      use shunt
      use usage
      use vdw
      use vdwpot
      use virial
      implicit none
      integer i,j,k
      integer ii,it,iv
      integer kk,kt,kv
      integer ig
      integer, allocatable :: iv14(:)
      real*8 e,de,p6,p12
      real*8 eps
      real*8 vlsc(2),dvlsc(2),d2vlsc(2)
      real*8 vshf(2),dvshf(2),d2vshf(2)
      real*8 rv,rdn,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 taper,dtaper
      real*8 dlambda,dlambda2
      real*8 dlde
      real*8 dldedx,dldedy,dldedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: vscale(:)
      real*8, allocatable :: dfl(:,:)
      logical proceed,usei
      logical mutik
      logical dodl,dod2,mutdl,mutd2
      character*6 mode
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
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
c     take lambda derivatives as "vdlflags" sets for the van der
c     Waals lambda
c
      call vdlflags (dodl,dod2)
c
c     set the soft core lambda terms for each ligand group, along
c     with their lambda derivatives as needed
c
      if (dodl) then
         call grpsoftd (scexplj,scalphalj,scexpljm,scexpljp,vlsc,dvlsc,
     &                  d2vlsc,vshf,dvshf,d2vshf)
      else
         call grpsoft (scexplj,scalphalj,scexpljm,scexpljp,vlsc,vshf)
      end if
c
c     the lambda derivative of the force is reduced over threads in
c     a local array that only has full size when it is needed
c
      if (dod2) then
         allocate (dfl(3,n))
      else
         allocate (dfl(3,1))
      end if
      dfl = 0.0d0
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
!$OMP PARALLEL default(private) shared(nvdw,ivdw,jvdw,ired,kred,
!$OMP& xred,yred,zred,use,nvlst,vlst,n12,n13,n14,n15,i12,i13,i14,
!$OMP& i15,v2scale,v3scale,v4scale,v5scale,use_group,off2,radmin,
!$OMP& epsilon,radmin4,epsilon4,mutcls,mutg,cut2,c0,c1,
!$OMP& c2,c3,c4,c5,dodl,dod2,vlsc,dvlsc,d2vlsc,vshf,dvshf,d2vshf)
!$OMP& firstprivate(vscale,iv14)
!$OMP& shared(ev,dev,evvir,devdl,d2evdl2,dfl,devvirdl)
!$OMP DO reduction(+:ev,dev,evvir,devdl,d2evdl2,dfl,devvirdl)
c
c     find van der Waals energy and derivatives via neighbor list
c
      do ii = 1, nvdw
         i = ivdw(ii)
         it = jvdw(i)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
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
            ig = mutcls(mutg(k),mutg(i))
c
c     the two ligands of a relative dual topology never interact
c
            proceed = (ig .ge. 0)
            if (use_group .and. proceed)
     &         call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
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
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  rik = sqrt(rik2)
c
c     set use of lambda scaling for decoupling or annihilation
c
                  mutik = (ig .gt. 0)
                  mutdl = mutik .and. dodl
                  mutd2 = mutik .and. dod2
c
c     get interaction energy, via soft core lambda scaling as needed
c
                  if (mutdl) then
                     call eljsc4 (rik,rik2,rv,eps,ig,vlsc,dvlsc,
     &                            d2vlsc,vshf,dvshf,d2vshf,e,de,
     &                            dlambda,dlambda2,dlde)
                  else if (mutik) then
                     call eljsc (rik,rik2,rv,eps,ig,vlsc,vshf,
     &                           .true.,e,de)
                  else
                     p6 = rv**6 / rik2**3
                     p12 = p6 * p6
                     e = eps * (p12 - 2.0d0*p6)
                     de = eps * (p12-p6) * (-12.0d0/rik)
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
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                     if (mutd2) then
                        dlde = dlambda * dtaper + dlde * taper
                        dlambda2 = dlambda2 * taper
                     end if
                     if (mutdl)  dlambda = dlambda * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     if (mutdl)  dlambda = dlambda * fgrp
                     if (mutd2) then
                        dlambda2 = dlambda2 * fgrp
                        dlde = dlde * fgrp
                     end if
                  end if
c
c     find the chain rule terms for derivative components; atoms at
c     the same point have no direction and take zero force
c
                  if (rik .ne. 0.0d0) then
                     de = de / rik
                  else
                     de = 0.0d0
                  end if
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                  ev = ev + e
                  if (mutd2) then
                     if (rik .ne. 0.0d0) then
                        dlde = dlde / rik
                     else
                        dlde = 0.0d0
                     end if
                     dldedx = dlde * xr
                     dldedy = dlde * yr
                     dldedz = dlde * zr
                  end if
                  if (i .eq. iv) then
                     dev(1,i) = dev(1,i) + dedx
                     dev(2,i) = dev(2,i) + dedy
                     dev(3,i) = dev(3,i) + dedz
                     if (mutd2) then
                        dfl(1,i) = dfl(1,i) + dldedx
                        dfl(2,i) = dfl(2,i) + dldedy
                        dfl(3,i) = dfl(3,i) + dldedz
                     end if
                  else
                     dev(1,i) = dev(1,i) + dedx*redi
                     dev(2,i) = dev(2,i) + dedy*redi
                     dev(3,i) = dev(3,i) + dedz*redi
                     dev(1,iv) = dev(1,iv) + dedx*rediv
                     dev(2,iv) = dev(2,iv) + dedy*rediv
                     dev(3,iv) = dev(3,iv) + dedz*rediv
                     if (mutd2) then
                        dfl(1,i) = dfl(1,i) + dldedx*redi
                        dfl(2,i) = dfl(2,i) + dldedy*redi
                        dfl(3,i) = dfl(3,i) + dldedz*redi
                        dfl(1,iv) = dfl(1,iv) + dldedx*rediv
                        dfl(2,iv) = dfl(2,iv) + dldedy*rediv
                        dfl(3,iv) = dfl(3,iv) + dldedz*rediv
                     end if
                  end if
                  if (k .eq. kv) then
                     dev(1,k) = dev(1,k) - dedx
                     dev(2,k) = dev(2,k) - dedy
                     dev(3,k) = dev(3,k) - dedz
                     if (mutd2) then
                        dfl(1,k) = dfl(1,k) - dldedx
                        dfl(2,k) = dfl(2,k) - dldedy
                        dfl(3,k) = dfl(3,k) - dldedz
                     end if
                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
                     dev(1,k) = dev(1,k) - dedx*redk
                     dev(2,k) = dev(2,k) - dedy*redk
                     dev(3,k) = dev(3,k) - dedz*redk
                     dev(1,kv) = dev(1,kv) - dedx*redkv
                     dev(2,kv) = dev(2,kv) - dedy*redkv
                     dev(3,kv) = dev(3,kv) - dedz*redkv
                     if (mutd2) then
                        dfl(1,k) = dfl(1,k) - dldedx*redk
                        dfl(2,k) = dfl(2,k) - dldedy*redk
                        dfl(3,k) = dfl(3,k) - dldedz*redk
                        dfl(1,kv) = dfl(1,kv) - dldedx*redkv
                        dfl(2,kv) = dfl(2,kv) - dldedy*redkv
                        dfl(3,kv) = dfl(3,kv) - dldedz*redkv
                     end if
                  end if
                  if (mutdl)  devdl = devdl + dlambda
                  if (mutd2)  d2evdl2 = d2evdl2 + dlambda2
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  evvir(1,1) = evvir(1,1) + vxx
                  evvir(2,1) = evvir(2,1) + vyx
                  evvir(3,1) = evvir(3,1) + vzx
                  evvir(1,2) = evvir(1,2) + vyx
                  evvir(2,2) = evvir(2,2) + vyy
                  evvir(3,2) = evvir(3,2) + vzy
                  evvir(1,3) = evvir(1,3) + vzx
                  evvir(2,3) = evvir(2,3) + vzy
                  evvir(3,3) = evvir(3,3) + vzz
                  if (mutd2) then
                     call dlcvir (xr,yr,zr,dldedx,dldedy,dldedz,
     &                            devvirdl)
                  end if
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
c     add the lambda derivative of the force to its global array
c
      if (dod2) then
         do i = 1, n
            dfvdl(1,i) = dfvdl(1,i) + dfl(1,i)
            dfvdl(2,i) = dfvdl(2,i) + dfl(2,i)
            dfvdl(3,i) = dfvdl(3,i) + dfl(3,i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (vscale)
      deallocate (dfl)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine elj1d  --  Lennard-Jones derivs for smoothing  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "elj1d" calculates the Lennard-Jones 6-12 van der Waals energy
c      and its first derivatives via a Gaussian approximation for
c      potential energy smoothing
c
c
      subroutine elj1d
      use math
      use vdwpot
      implicit none
c
c
c     set coefficients for a two-Gaussian fit to Lennard-Jones
c
      ngauss = 2
      igauss(1,1) = 14487.1d0
      igauss(2,1) = 9.05148d0 * twosix**2
      igauss(1,2) = -5.55338d0
      igauss(2,2) = 1.22536d0 * twosix**2
c
c     compute Gaussian approximation to Lennard-Jones potential
c
      call egauss1
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine elj1e  --  Lennard-Jones derivs for stophat  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "elj1e" calculates the van der Waals interaction energy and its
c     first derivatives for use with stophat potential energy smoothing
c
c
      subroutine elj1e
      use atomid
      use atoms
      use couple
      use deriv
      use energi
      use group
      use usage
      use vdw
      use vdwpot
      use warp
      implicit none
      integer i,j,k
      integer ii,it,iv
      integer kk,kt,kv
      integer, allocatable :: iv14(:)
      real*8 e,de,rdn
      real*8 p6,denom
      real*8 eps,rv,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rik,rik2
      real*8 rik3,rik4
      real*8 rik5,rik6
      real*8 rik7,rik8
      real*8 width,width2
      real*8 width3,width4
      real*8 width5,width6
      real*8 width7
      real*8, allocatable :: vscale(:)
      logical proceed,usei
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
      end do
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
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the extent of smoothing to be performed
c
      width = deform * diffv
      width2 = width * width
      width3 = width2 * width
      width4 = width2 * width2
      width5 = width2 * width3
      width6 = width3 * width3
      width7 = width3 * width4
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
c     find van der Waals energy and derivatives via double loop
c
      do ii = 1, nvdw-1
         i = ivdw(ii)
         it = jvdw(i)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
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
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               rik2 = xr*xr + yr*yr + zr*zr
               eps = epsilon(kt,it)
               rv = radmin(kt,it)
               if (iv14(k) .eq. i) then
                  eps = epsilon4(kt,it)
                  rv = radmin4(kt,it)
               end if
               eps = eps * vscale(k)
               p6 = rv**6
               rik = sqrt(rik2)
               rik3 = rik2 * rik
               rik4 = rik2 * rik2
               rik5 = rik2 * rik3
               rik6 = rik3 * rik3
               rik7 = rik3 * rik4
               rik8 = rik4 * rik4
               denom = rik*(rik+2.0d0*width)
               denom = denom**9
c
c     transform the potential function via smoothing
c
               e = rik5 * (30.0d0*rik7 + 360.0d0*rik6*width
     &                + 1800.0d0*rik5*width2 + 4800.0d0*rik4*width3
     &                + 7200.0d0*rik3*width4 + 5760.0d0*rik2*width5
     &                + 1920.0d0*rik*width6)
               e = -e + p6 * (15.0d0*rik6 + 90.0d0*rik5*width
     &                + 288.0d0*rik4*width2 + 552.0d0*rik3*width3
     &                + 648.0d0*rik2*width4 + 432.0d0*rik*width5
     &                + 128.0d0*width6)
               e = e*eps*p6 / (15.0d0*denom)
               de = rik5 * (5.0d0*rik8 + 65.0d0*rik7*width
     &                 + 360.0d0*rik6*width2 + 1100.0d0*rik5*width3
     &                 + 2000.0d0*rik4*width4 + 2160.0d0*rik3*width5
     &                 + 1280.0d0*rik2*width6 + 320.0d0*rik*width7)
               de = de - p6 * (5.0d0*rik7 + 35.0d0*rik6*width
     &                 + 132.0d0*rik5*width2 + 310.0d0*rik4*width3
     &                 + 472.0d0*rik3*width4 + 456.0d0*rik2*width5
     &                 + 256.0d0*rik*width6 + 64.0d0*width7)
               de = 12.0d0*de*eps*p6
     &                 / (5.0d0*denom*rik*(rik+2.0d0*width))
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  de = de * fgrp
               end if
c
c     find the chain rule terms for derivative components
c
               de = de / rik
               dedx = de * xr
               dedy = de * yr
               dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
               ev = ev + e
               if (i .eq. iv) then
                  dev(1,i) = dev(1,i) + dedx
                  dev(2,i) = dev(2,i) + dedy
                  dev(3,i) = dev(3,i) + dedz
               else
                  dev(1,i) = dev(1,i) + dedx*redi
                  dev(2,i) = dev(2,i) + dedy*redi
                  dev(3,i) = dev(3,i) + dedz*redi
                  dev(1,iv) = dev(1,iv) + dedx*rediv
                  dev(2,iv) = dev(2,iv) + dedy*rediv
                  dev(3,iv) = dev(3,iv) + dedz*rediv
               end if
               if (k .eq. kv) then
                  dev(1,k) = dev(1,k) - dedx
                  dev(2,k) = dev(2,k) - dedy
                  dev(3,k) = dev(3,k) - dedz
               else
                  redk = kred(k)
                  redkv = 1.0d0 - redk
                  dev(1,k) = dev(1,k) - dedx*redk
                  dev(2,k) = dev(2,k) - dedy*redk
                  dev(3,k) = dev(3,k) - dedz*redk
                  dev(1,kv) = dev(1,kv) - dedx*redkv
                  dev(2,kv) = dev(2,kv) - dedy*redkv
                  dev(3,kv) = dev(3,kv) - dedz*redkv
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
c     perform deallocation of some local arrays
c
      deallocate (iv14)
      deallocate (vscale)
      return
      end

c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echgtrn1  --  charge transfer energy & derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echgtrn1" calculates the charge transfer energy and first
c     derivatives with respect to Cartesian coordinates, along with
c     the lambda derivatives when polarization carries one
c
c     charge transfer is a pairwise sum that is linear in its lambda
c     scale from "ctlmdasc", so the lambda derivatives are the unscaled
c     pair energy, force and virial of the mutated pairs times the
c     first or second derivative of that scale; they are added to the
c     polarization lambda derivatives, which "gradient" zeroes before
c     it calls the polarization routines and then this one
c
c
      subroutine echgtrn1
      use limits
      use virial
      implicit none
      integer i,j
c
c
c     choose method for summing over charge transfer interactions
c
      if (use_mlist) then
         call echgtrn1b
      else
         call echgtrn1a
      end if
c
c     add the charge transfer virial to main virial
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + ectvir(j,i)
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine echgtrn1a  --  charge transfer derivs via loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "echgtrn1a" calculates the charge transfer interaction energy,
c     first derivatives and lambda derivatives using a double loop
c
c
      subroutine echgtrn1a
      use atoms
      use bound
      use chgpot
      use chgtrn
      use cell
      use couple
      use ctrpot
      use deriv
      use dlmda
      use energi
      use group
      use mplpot
      use mpole
      use mutant
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk
      integer mutgi,mutgk
      integer jcell
      real*8 e,de,f,fgrp
      real*8 rr1,r,r2
      real*8 r3,r4,r5
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 chgi,chgk
      real*8 chgik
      real*8 alphai,alphak
      real*8 alphaik
      real*8 expi,expk
      real*8 expik
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 taper,dtaper
      real*8 dlfrcx,dlfrcy,dlfrcz
      real*8 dldvxx,dldvyy,dldvzz
      real*8 dldvxy,dldvxz,dldvyz
      real*8 scik,dscik
      real*8 ctsc(0:2,0:2)
      real*8 dctsc(0:2,0:2)
      real*8 d2ctsc(0:2,0:2)
      real*8, allocatable :: mscale(:)
      logical proceed,usei
      logical live(0:2,0:2)
      logical dodl,dod2
      character*6 mode
c
c
c     zero out the charge transfer energy, derivatives and virial
c
      ect = 0.0d0
      do i = 1, n
         dect(1,i) = 0.0d0
         dect(2,i) = 0.0d0
         dect(3,i) = 0.0d0
      end do
      do i = 1, 3
         do j = 1, 3
            ectvir(j,i) = 0.0d0
         end do
      end do
      if (nct .eq. 0)  return
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
c     set the polarization lambda scale of each pair of atom groups,
c     taking lambda derivatives whenever polarization has one, since
c     "lmdachain" zeroes them while the polarization map is flat; a
c     group pair is skipped once it carries neither energy nor a
c     needed lambda derivative
c
      dodl = use_pdlmda
      dod2 = (dodl .and. use_d2lmda)
      call ctlmdasc (dodl,ctsc,dctsc,d2ctsc,live)
c
c     calculate the charge transfer energy and derivatives
c
      do ii = 1, npole-1
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         chgi = chgct(i)
         alphai = dmpct(i)
         if (alphai .eq. 0.0d0)  alphai = 1000.0d0
         usei = use(i)
         mutgi = mutg(i)
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
         do kk = ii+1, npole
            k = ipole(kk)
            mutgk = mutg(k)
            proceed = live(mutgk,mutgi)
            if (proceed .and. use_group) then
               call groups (proceed,fgrp,i,k,0,0,0,0)
               if (.not. use_intra)  proceed = .true.
            end if
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
                  chgk = chgct(k)
                  alphak = dmpct(k)
                  if (alphak .eq. 0.0d0)  alphak = 1000.0d0
                  if (ctrntyp .eq. 'SEPARATE') then
                     expi = exp(-alphai*r)
                     expk = exp(-alphak*r)
                     e = -chgi*expk - chgk*expi
                     de = chgi*expk*alphak + chgk*expi*alphai
                  else
                     chgik = sqrt(abs(chgi*chgk))
                     alphaik = 0.5d0 * (alphai+alphak)
                     expik = exp(-alphaik*r)
                     e = -chgik * expik
                     de = -e * alphaik
                  end if
                  e = f * e * mscale(k)
                  de = f * de * mscale(k)
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
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     increment the polarization lambda derivatives from the
c     unscaled pair terms, then apply polarization lambda scaling
c     for interaction annihilation
c
                  if (mutgi+mutgk .ne. 0) then
                     scik = ctsc(mutgk,mutgi)
                     dscik = dctsc(mutgk,mutgi)
                     if (dodl)  depdl = depdl + dscik*e
                     if (dod2) then
                        d2epdl2 = d2epdl2 + d2ctsc(mutgk,mutgi)*e
                        dlfrcx = dscik * de * xr * rr1
                        dlfrcy = dscik * de * yr * rr1
                        dlfrcz = dscik * de * zr * rr1
                        dfpdl(1,i) = dfpdl(1,i) - dlfrcx
                        dfpdl(2,i) = dfpdl(2,i) - dlfrcy
                        dfpdl(3,i) = dfpdl(3,i) - dlfrcz
                        dfpdl(1,k) = dfpdl(1,k) + dlfrcx
                        dfpdl(2,k) = dfpdl(2,k) + dlfrcy
                        dfpdl(3,k) = dfpdl(3,k) + dlfrcz
                        dldvxx = xr * dlfrcx
                        dldvxy = yr * dlfrcx
                        dldvxz = zr * dlfrcx
                        dldvyy = yr * dlfrcy
                        dldvyz = zr * dlfrcy
                        dldvzz = zr * dlfrcz
                        depvirdl(1,1) = depvirdl(1,1) + dldvxx
                        depvirdl(2,1) = depvirdl(2,1) + dldvxy
                        depvirdl(3,1) = depvirdl(3,1) + dldvxz
                        depvirdl(1,2) = depvirdl(1,2) + dldvxy
                        depvirdl(2,2) = depvirdl(2,2) + dldvyy
                        depvirdl(3,2) = depvirdl(3,2) + dldvyz
                        depvirdl(1,3) = depvirdl(1,3) + dldvxz
                        depvirdl(2,3) = depvirdl(2,3) + dldvyz
                        depvirdl(3,3) = depvirdl(3,3) + dldvzz
                     end if
                     e = e * scik
                     de = de * scik
                  end if
c
c     compute the force components for this interaction
c
                  frcx = de * xr * rr1
                  frcy = de * yr * rr1
                  frcz = de * zr * rr1
c     
c     increment the total charge transfer energy and derivatives
c
                  ect = ect + e
                  dect(1,i) = dect(1,i) - frcx
                  dect(2,i) = dect(2,i) - frcy
                  dect(3,i) = dect(3,i) - frcz
                  dect(1,k) = dect(1,k) + frcx
                  dect(2,k) = dect(2,k) + frcy
                  dect(3,k) = dect(3,k) + frcz
c
c     increment the internal virial tensor components
c
                  vxx = xr * frcx
                  vxy = yr * frcx
                  vxz = zr * frcx
                  vyy = yr * frcy
                  vyz = zr * frcy
                  vzz = zr * frcz
                  ectvir(1,1) = ectvir(1,1) + vxx
                  ectvir(2,1) = ectvir(2,1) + vxy
                  ectvir(3,1) = ectvir(3,1) + vxz
                  ectvir(1,2) = ectvir(1,2) + vxy
                  ectvir(2,2) = ectvir(2,2) + vyy
                  ectvir(3,2) = ectvir(3,2) + vyz
                  ectvir(1,3) = ectvir(1,3) + vxz
                  ectvir(2,3) = ectvir(2,3) + vyz
                  ectvir(3,3) = ectvir(3,3) + vzz
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
c     calculate interaction components with other unit cells
c
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            chgi = chgct(i)
            alphai = dmpct(i)
            if (alphai .eq. 0.0d0)  alphai = 1000.0d0
            usei = use(i)
            mutgi = mutg(i)
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
            do kk = ii, npole
               k = ipole(kk)
               mutgk = mutg(k)
               proceed = live(mutgk,mutgi)
               if (proceed .and. use_group) then
                  call groups (proceed,fgrp,i,k,0,0,0,0)
                  if (.not. use_intra)  proceed = .true.
               end if
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
                        chgk = chgct(k)
                        alphak = dmpct(k)
                        if (alphak .eq. 0.0d0)  alphak = 1000.0d0
                        if (ctrntyp .eq. 'SEPARATE') then
                           expi = exp(-alphai*r)
                           expk = exp(-alphak*r)
                           e = -chgi*expk - chgk*expi
                           de = chgi*expk*alphak + chgk*expi*alphai
                        else
                           chgik = sqrt(abs(chgi*chgk))
                           alphaik = 0.5d0 * (alphai+alphak)
                           expik = exp(-alphaik*r)
                           e = -chgik * expik
                           de = -e * alphaik
                        end if
                        e = f * e * mscale(k)
                        de = f * de * mscale(k)
                        if (i .eq. k) then
                           e = 0.5d0 * e
                           de = 0.5d0 * de
                        end if
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
                           de = e*dtaper + de*taper
                           e = e * taper
                        end if
c
c     scale the interaction based on its group membership
c
                        if (use_group) then
                           e = e * fgrp
                           de = de * fgrp
                        end if
c
c     increment the polarization lambda derivatives from the
c     unscaled pair terms, then apply polarization lambda scaling
c     for interaction annihilation
c
                        if (mutgi+mutgk .ne. 0) then
                           scik = ctsc(mutgk,mutgi)
                           dscik = dctsc(mutgk,mutgi)
                           if (dodl)  depdl = depdl + dscik*e
                           if (dod2) then
                              d2epdl2 = d2epdl2 + d2ctsc(mutgk,mutgi)*e
                              dlfrcx = dscik * de * xr * rr1
                              dlfrcy = dscik * de * yr * rr1
                              dlfrcz = dscik * de * zr * rr1
                              dfpdl(1,i) = dfpdl(1,i) - dlfrcx
                              dfpdl(2,i) = dfpdl(2,i) - dlfrcy
                              dfpdl(3,i) = dfpdl(3,i) - dlfrcz
                              dfpdl(1,k) = dfpdl(1,k) + dlfrcx
                              dfpdl(2,k) = dfpdl(2,k) + dlfrcy
                              dfpdl(3,k) = dfpdl(3,k) + dlfrcz
                              dldvxx = xr * dlfrcx
                              dldvxy = yr * dlfrcx
                              dldvxz = zr * dlfrcx
                              dldvyy = yr * dlfrcy
                              dldvyz = zr * dlfrcy
                              dldvzz = zr * dlfrcz
                              depvirdl(1,1) = depvirdl(1,1) + dldvxx
                              depvirdl(2,1) = depvirdl(2,1) + dldvxy
                              depvirdl(3,1) = depvirdl(3,1) + dldvxz
                              depvirdl(1,2) = depvirdl(1,2) + dldvxy
                              depvirdl(2,2) = depvirdl(2,2) + dldvyy
                              depvirdl(3,2) = depvirdl(3,2) + dldvyz
                              depvirdl(1,3) = depvirdl(1,3) + dldvxz
                              depvirdl(2,3) = depvirdl(2,3) + dldvyz
                              depvirdl(3,3) = depvirdl(3,3) + dldvzz
                           end if
                           e = e * scik
                           de = de * scik
                        end if
c
c     compute the force components for this interaction
c
                        frcx = de * xr * rr1
                        frcy = de * yr * rr1
                        frcz = de * zr * rr1
c     
c     increment the total charge transfer energy and derivatives
c
                        ect = ect + e
                        dect(1,i) = dect(1,i) - frcx
                        dect(2,i) = dect(2,i) - frcy
                        dect(3,i) = dect(3,i) - frcz
                        dect(1,k) = dect(1,k) + frcx
                        dect(2,k) = dect(2,k) + frcy
                        dect(3,k) = dect(3,k) + frcz
c
c     increment the internal virial tensor components
c
                        vxx = xr * frcx
                        vxy = yr * frcx
                        vxz = zr * frcx
                        vyy = yr * frcy
                        vyz = zr * frcy
                        vzz = zr * frcz
                        ectvir(1,1) = ectvir(1,1) + vxx
                        ectvir(2,1) = ectvir(2,1) + vxy
                        ectvir(3,1) = ectvir(3,1) + vxz
                        ectvir(1,2) = ectvir(1,2) + vxy
                        ectvir(2,2) = ectvir(2,2) + vyy
                        ectvir(3,2) = ectvir(3,2) + vyz
                        ectvir(1,3) = ectvir(1,3) + vxz
                        ectvir(2,3) = ectvir(2,3) + vyz
                        ectvir(3,3) = ectvir(3,3) + vzz
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
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine echgtrn1b  --  charge transfer derivs via list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "echgtrn1b" calculates the charge transfer energy, first
c     derivatives and lambda derivatives using a pairwise neighbor list
c
c
      subroutine echgtrn1b
      use atoms
      use bound
      use chgpot
      use chgtrn
      use cell
      use couple
      use ctrpot
      use deriv
      use dlmda
      use energi
      use group
      use mplpot
      use mpole
      use mutant
      use neigh
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer mutgi,mutgk
      real*8 e,de,f,fgrp
      real*8 rr1,r,r2
      real*8 r3,r4,r5
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 chgi,chgk
      real*8 chgik
      real*8 alphai,alphak
      real*8 alphaik
      real*8 expi,expk
      real*8 expik
      real*8 frcx,frcy,frcz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 taper,dtaper
      real*8 dlfrcx,dlfrcy,dlfrcz
      real*8 dldvxx,dldvyy,dldvzz
      real*8 dldvxy,dldvxz,dldvyz
      real*8 scik,dscik
      real*8 ctsc(0:2,0:2)
      real*8 dctsc(0:2,0:2)
      real*8 d2ctsc(0:2,0:2)
      real*8, allocatable :: mscale(:)
      logical proceed,usei
      logical live(0:2,0:2)
      logical dodl,dod2
      character*6 mode
c
c
c     zero out the charge transfer energy, derivatives and virial
c
      ect = 0.0d0
      do i = 1, n
         dect(1,i) = 0.0d0
         dect(2,i) = 0.0d0
         dect(3,i) = 0.0d0
      end do
      do i = 1, 3
         do j = 1, 3
            ectvir(j,i) = 0.0d0
         end do
      end do
      if (nct .eq. 0)  return
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
c     set the polarization lambda scale of each pair of atom groups,
c     taking lambda derivatives whenever polarization has one, since
c     "lmdachain" zeroes them while the polarization map is flat; a
c     group pair is skipped once it carries neither energy nor a
c     needed lambda derivative
c
      dodl = use_pdlmda
      dod2 = (dodl .and. use_d2lmda)
      call ctlmdasc (dodl,ctsc,dctsc,d2ctsc,live)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,chgct,dmpct,n12,i12,n13,i13,
!$OMP& n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,nelst,
!$OMP& elst,use,use_group,use_intra,use_bounds,ctrntyp,f,off2,
!$OMP& ctsc,dctsc,d2ctsc,live,dodl,dod2,mutg,cut2,c0,c1,c2,c3,
!$OMP& c4,c5)
!$OMP& firstprivate(mscale) shared(ect,dect,ectvir,depdl,d2epdl2,
!$OMP& depvirdl,dfpdl)
!$OMP DO reduction(+:ect,dect,ectvir,depdl,d2epdl2,depvirdl,dfpdl)
c
c     compute the charge transfer energy and derivatives
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         chgi = chgct(i)
         alphai = dmpct(i)
         if (alphai .eq. 0.0d0)  alphai = 1000.0d0
         usei = use(i)
         mutgi = mutg(i)
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
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            mutgk = mutg(k)
            proceed = live(mutgk,mutgi)
            if (proceed .and. use_group) then
               call groups (proceed,fgrp,i,k,0,0,0,0)
               if (.not. use_intra)  proceed = .true.
            end if
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
                  chgk = chgct(k)
                  alphak = dmpct(k)
                  if (alphak .eq. 0.0d0)  alphak = 1000.0d0
                  if (ctrntyp .eq. 'SEPARATE') then
                     expi = exp(-alphai*r)
                     expk = exp(-alphak*r)
                     e = -chgi*expk - chgk*expi
                     de = chgi*expk*alphak + chgk*expi*alphai
                  else
                     chgik = sqrt(abs(chgi*chgk))
                     alphaik = 0.5d0 * (alphai+alphak)
                     expik = exp(-alphaik*r)
                     e = -chgik * expik
                     de = -e * alphaik
                  end if
                  e = f * e * mscale(k)
                  de = f * de * mscale(k)
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
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     increment the polarization lambda derivatives from the
c     unscaled pair terms, then apply polarization lambda scaling
c     for interaction annihilation
c
                  if (mutgi+mutgk .ne. 0) then
                     scik = ctsc(mutgk,mutgi)
                     dscik = dctsc(mutgk,mutgi)
                     if (dodl)  depdl = depdl + dscik*e
                     if (dod2) then
                        d2epdl2 = d2epdl2 + d2ctsc(mutgk,mutgi)*e
                        dlfrcx = dscik * de * xr * rr1
                        dlfrcy = dscik * de * yr * rr1
                        dlfrcz = dscik * de * zr * rr1
                        dfpdl(1,i) = dfpdl(1,i) - dlfrcx
                        dfpdl(2,i) = dfpdl(2,i) - dlfrcy
                        dfpdl(3,i) = dfpdl(3,i) - dlfrcz
                        dfpdl(1,k) = dfpdl(1,k) + dlfrcx
                        dfpdl(2,k) = dfpdl(2,k) + dlfrcy
                        dfpdl(3,k) = dfpdl(3,k) + dlfrcz
                        dldvxx = xr * dlfrcx
                        dldvxy = yr * dlfrcx
                        dldvxz = zr * dlfrcx
                        dldvyy = yr * dlfrcy
                        dldvyz = zr * dlfrcy
                        dldvzz = zr * dlfrcz
                        depvirdl(1,1) = depvirdl(1,1) + dldvxx
                        depvirdl(2,1) = depvirdl(2,1) + dldvxy
                        depvirdl(3,1) = depvirdl(3,1) + dldvxz
                        depvirdl(1,2) = depvirdl(1,2) + dldvxy
                        depvirdl(2,2) = depvirdl(2,2) + dldvyy
                        depvirdl(3,2) = depvirdl(3,2) + dldvyz
                        depvirdl(1,3) = depvirdl(1,3) + dldvxz
                        depvirdl(2,3) = depvirdl(2,3) + dldvyz
                        depvirdl(3,3) = depvirdl(3,3) + dldvzz
                     end if
                     e = e * scik
                     de = de * scik
                  end if
c
c     compute the force components for this interaction
c
                  frcx = de * xr * rr1
                  frcy = de * yr * rr1
                  frcz = de * zr * rr1
c     
c     increment the total charge transfer energy and derivatives
c
                  ect = ect + e
                  dect(1,i) = dect(1,i) - frcx
                  dect(2,i) = dect(2,i) - frcy
                  dect(3,i) = dect(3,i) - frcz
                  dect(1,k) = dect(1,k) + frcx
                  dect(2,k) = dect(2,k) + frcy
                  dect(3,k) = dect(3,k) + frcz
c
c     increment the internal virial tensor components
c
                  vxx = xr * frcx
                  vxy = yr * frcx
                  vxz = zr * frcx
                  vyy = yr * frcy
                  vyz = zr * frcy
                  vzz = zr * frcz
                  ectvir(1,1) = ectvir(1,1) + vxx
                  ectvir(2,1) = ectvir(2,1) + vxy
                  ectvir(3,1) = ectvir(3,1) + vxz
                  ectvir(1,2) = ectvir(1,2) + vxy
                  ectvir(2,2) = ectvir(2,2) + vyy
                  ectvir(3,2) = ectvir(3,2) + vyz
                  ectvir(1,3) = ectvir(1,3) + vxz
                  ectvir(2,3) = ectvir(2,3) + vyz
                  ectvir(3,3) = ectvir(3,3) + vzz
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end

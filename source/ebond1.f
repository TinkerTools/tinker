c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ebond1  --  bond stretch energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ebond1" calculates the bond stretching energy and
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine ebond1
      use sizes
      use atoms
      use bndpot
      use bndstr
      use bound
      use deriv
      use energi
      use group
      use usage
      use virial
      implicit none
      integer i,ia,ib
      real*8 e,ebo
      real*8 ideal,force
      real*8 expterm,bde,fgrp
      real*8 dt,dt2,deddt
      real*8 de,dedx,dedy,dedz
      real*8 xab,yab,zab,rab
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 viro(3,3)
      real*8, allocatable :: debo(:,:)
      logical proceed
c
c
c     zero out the bond energy and first derivatives
c
      eb = 0.0d0
      do i = 1, n
         deb(1,i) = 0.0d0
         deb(2,i) = 0.0d0
         deb(3,i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (debo(3,n))
c
c     transfer global to local copies for OpenMP calculation
c
      ebo = eb
      do i = 1, n
         debo(1,i) = deb(1,i)
         debo(2,i) = deb(2,i)
         debo(3,i) = deb(3,i)
      end do
      do i = 1, 3
         viro(1,i) = vir(1,i)
         viro(2,i) = vir(2,i)
         viro(3,i) = vir(3,i)
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nbond,ibnd,bl,bk,use,
!$OMP& x,y,z,cbnd,qbnd,bndtyp,bndunit,use_group,use_polymer)
!$OMP& shared(ebo,debo,viro)
!$OMP DO reduction(+:ebo,debo,viro) schedule(guided)
c
c     calculate the bond stretch energy and first derivatives
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ideal = bl(i)
         force = bk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,0,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib))
c
c     compute the value of the bond length deviation
c
         if (proceed) then
            xab = x(ia) - x(ib)
            yab = y(ia) - y(ib)
            zab = z(ia) - z(ib)
            if (use_polymer)  call image (xab,yab,zab)
            rab = sqrt(xab*xab + yab*yab + zab*zab)
            dt = rab - ideal
c
c     harmonic potential uses Taylor expansion of Morse potential
c     through the fourth power of the bond length deviation
c
            if (bndtyp .eq. 'HARMONIC') then
               dt2 = dt * dt
               e = bndunit * force * dt2 * (1.0d0+cbnd*dt+qbnd*dt2)
               deddt = 2.0d0 * bndunit * force * dt
     &                    * (1.0d0+1.5d0*cbnd*dt+2.0d0*qbnd*dt2)
c
c     Morse potential uses energy = BDE * (1 - e**(-alpha*dt))**2)
c     with the approximations alpha = sqrt(ForceConst/BDE) = -2
c     and BDE = Bond Dissociation Energy = ForceConst/alpha**2
c
            else if (bndtyp .eq. 'MORSE') then
               expterm = exp(-2.0d0*dt)
               bde = 0.25d0 * bndunit * force
               e = bde * (1.0d0-expterm)**2
               deddt = 4.0d0 * bde * (1.0d0-expterm) * expterm
            end if
c
c     scale the interaction based on its group membership
c
            if (use_group) then
               e = e * fgrp
               deddt = deddt * fgrp
            end if
c
c     compute chain rule terms needed for derivatives
c
            if (rab .eq. 0.0d0) then
               de = 0.0d0
            else
               de = deddt / rab
            end if
            dedx = de * xab
            dedy = de * yab
            dedz = de * zab
c
c     increment the total bond energy and first derivatives
c
            ebo = ebo + e
            debo(1,ia) = debo(1,ia) + dedx
            debo(2,ia) = debo(2,ia) + dedy
            debo(3,ia) = debo(3,ia) + dedz
            debo(1,ib) = debo(1,ib) - dedx
            debo(2,ib) = debo(2,ib) - dedy
            debo(3,ib) = debo(3,ib) - dedz
c
c     increment the internal virial tensor components
c
            vxx = xab * dedx
            vyx = yab * dedx
            vzx = zab * dedx
            vyy = yab * dedy
            vzy = zab * dedy
            vzz = zab * dedz
            viro(1,1) = viro(1,1) + vxx
            viro(2,1) = viro(2,1) + vyx
            viro(3,1) = viro(3,1) + vzx
            viro(1,2) = viro(1,2) + vyx
            viro(2,2) = viro(2,2) + vyy
            viro(3,2) = viro(3,2) + vzy
            viro(1,3) = viro(1,3) + vzx
            viro(2,3) = viro(2,3) + vzy
            viro(3,3) = viro(3,3) + vzz
         end if
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     transfer local to global copies for OpenMP calculation
c
      eb = ebo
      do i = 1, n
         deb(1,i) = debo(1,i)
         deb(2,i) = debo(2,i)
         deb(3,i) = debo(3,i)
      end do
      do i = 1, 3
         vir(1,i) = viro(1,i)
         vir(2,i) = viro(2,i)
         vir(3,i) = viro(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (debo)
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eurey1  --  bond stretch energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eurey1" calculates the Urey-Bradley interaction energy and
c     its first derivatives with respect to Cartesian coordinates
c
c
      subroutine eurey1
      use sizes
      use atoms
      use bound
      use deriv
      use energi
      use group
      use urey
      use urypot
      use usage
      use virial
      implicit none
      integer i,ia,ic
      real*8 e,eubo
      real*8 de,ideal,force
      real*8 dt,dt2,deddt,fgrp
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 xac,yac,zac,rac
      real*8 viro(3,3)
      real*8, allocatable :: deubo(:,:)
      logical proceed
c
c
c     zero out the Urey-Bradley energy and first derivatives
c
      eub = 0.0d0
      do i = 1, n
         deub(1,i) = 0.0d0
         deub(2,i) = 0.0d0
         deub(3,i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (deubo(3,n))
c
c     transfer global to local copies for OpenMP calculation
c
      eubo = eub
      do i = 1, n
         deubo(1,i) = deub(1,i)
         deubo(2,i) = deub(2,i)
         deubo(3,i) = deub(3,i)
      end do
      do i = 1, 3
         viro(1,i) = vir(1,i)
         viro(2,i) = vir(2,i)
         viro(3,i) = vir(3,i)
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nurey,iury,ul,uk,
!$OMP& use,x,y,z,cury,qury,ureyunit,use_group,use_polymer)
!$OMP& shared(eubo,deubo,viro)
!$OMP DO reduction(+:eubo,deubo,viro) schedule(guided)
c
c     calculate the Urey-Bradley 1-3 energy and first derivatives
c
      do i = 1, nurey
         ia = iury(1,i)
         ic = iury(3,i)
         ideal = ul(i)
         force = uk(i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ic,0,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ic))
c
c     compute the value of the 1-3 distance deviation
c
         if (proceed) then
            xac = x(ia) - x(ic)
            yac = y(ia) - y(ic)
            zac = z(ia) - z(ic)
            if (use_polymer)  call image (xac,yac,zac)
            rac = sqrt(xac*xac + yac*yac + zac*zac)
            dt = rac - ideal
            dt2 = dt * dt
            e = ureyunit * force * dt2 * (1.0d0+cury*dt+qury*dt2)
            deddt = 2.0d0 * ureyunit * force * dt
     &                 * (1.0d0+1.5d0*cury*dt+2.0d0*qury*dt2)
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
            de = deddt / rac
            dedx = de * xac
            dedy = de * yac
            dedz = de * zac
c
c     increment the total Urey-Bradley energy and first derivatives
c
            eubo = eubo + e
            deubo(1,ia) = deubo(1,ia) + dedx
            deubo(2,ia) = deubo(2,ia) + dedy
            deubo(3,ia) = deubo(3,ia) + dedz
            deubo(1,ic) = deubo(1,ic) - dedx
            deubo(2,ic) = deubo(2,ic) - dedy
            deubo(3,ic) = deubo(3,ic) - dedz
c
c     increment the internal virial tensor components
c
            vxx = xac * dedx
            vyx = yac * dedx
            vzx = zac * dedx
            vyy = yac * dedy
            vzy = zac * dedy
            vzz = zac * dedz
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
      eub = eubo
      do i = 1, n
         deub(1,i) = deubo(1,i)
         deub(2,i) = deubo(2,i)
         deub(3,i) = deubo(3,i)
      end do
      do i = 1, 3
         vir(1,i) = viro(1,i)
         vir(2,i) = viro(2,i)
         vir(3,i) = viro(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (deubo)
      return
      end

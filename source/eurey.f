c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine eurey  --  Urey-Bradley potential energy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "eurey" calculates the Urey-Bradley 1-3 interaction energy
c
c
      subroutine eurey
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'energi.i'
      include 'group.i'
      include 'urey.i'
      include 'urypot.i'
      include 'usage.i'
      integer i,ia,ic
      real*8 e,ideal,force
      real*8 dt,dt2,fgrp
      real*8 xac,yac,zac,rac
      logical proceed
c
c
c     zero out the Urey-Bradley interaction energy
c
      eub = 0.0d0
c
c     calculate the Urey-Bradley 1-3 energy term
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
c
c     calculate the Urey-Bradley energy for this interaction
c
            e = ureyunit * force * dt2 * (1.0d0+cury*dt+qury*dt2)
c
c     scale the interaction based on its group membership
c
            if (use_group)  e = e * fgrp
c
c     increment the total Urey-Bradley energy
c
            eub = eub + e
         end if
      end do
      return
      end

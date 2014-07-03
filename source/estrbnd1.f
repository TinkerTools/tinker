c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine estrbnd1   --  stretch-bend energy and derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "estrbnd1" calculates the stretch-bend potential energy and
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine estrbnd1
      use sizes
      use angbnd
      use angpot
      use atoms
      use bndstr
      use bound
      use deriv
      use energi
      use group
      use math
      use strbnd
      use usage
      use virial
      implicit none
      integer i,j,k,istrbnd
      integer ia,ib,ic
      real*8 e,ebao
      real*8 dr1,dr2,dt
      real*8 fgrp,angle
      real*8 force1,force2
      real*8 dot,cosine
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 rab,rab2
      real*8 rcb,rcb2
      real*8 xp,yp,zp,rp
      real*8 term1,term2
      real*8 termr,term1t,term2t
      real*8 ddtdxia,ddtdyia,ddtdzia
      real*8 ddtdxic,ddtdyic,ddtdzic
      real*8 ddrdxia,ddrdyia,ddrdzia
      real*8 ddrdxic,ddrdyic,ddrdzic
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 viro(3,3)
      real*8, allocatable :: debao(:,:)
      logical proceed
c
c
c     zero out the energy and first derivative components
c
      eba = 0.0d0
      do i = 1, n
         deba(1,i) = 0.0d0
         deba(2,i) = 0.0d0
         deba(3,i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (debao(3,n))
c
c     transfer global to local copies for OpenMP calculation
c
      ebao = eba
      do i = 1, n
         debao(1,i) = deba(1,i)
         debao(2,i) = deba(2,i)
         debao(3,i) = deba(3,i)
      end do
      do i = 1, 3
         viro(1,i) = vir(1,i)
         viro(2,i) = vir(2,i)
         viro(3,i) = vir(3,i)
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nstrbnd,isb,iang,sbk,
!$OMP& anat,bl,bk,use,x,y,z,stbnunit,use_group,use_polymer)
!$OMP& shared(ebao,debao,viro)
!$OMP DO reduction(+:ebao,debao,viro) schedule(guided)
c
c     calculate the stretch-bend energy and first derivatives
c
      do istrbnd = 1, nstrbnd
         i = isb(1,istrbnd)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         force1 = sbk(1,istrbnd)
         force2 = sbk(2,istrbnd)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,0,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic))
c
c     get the coordinates of the atoms in the angle
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
c
c     compute the value of the bond angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            if (use_polymer) then
               call image (xab,yab,zab)
               call image (xcb,ycb,zcb)
            end if
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            if (rab2.ne.0.0d0 .and. rcb2.ne.0.0d0) then
               rab = sqrt(rab2)
               rcb = sqrt(rcb2)
               xp = ycb*zab - zcb*yab
               yp = zcb*xab - xcb*zab
               zp = xcb*yab - ycb*xab
               rp = sqrt(xp*xp + yp*yp + zp*zp)
               rp = max(rp,0.001d0)
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / (rab*rcb)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
c
c     find chain rule terms for the bond angle deviation
c
               dt = angle - anat(i)
               term1 = -radian / (rab2*rp)
               term2 = radian / (rcb2*rp)
               ddtdxia = term1 * (yab*zp-zab*yp)
               ddtdyia = term1 * (zab*xp-xab*zp)
               ddtdzia = term1 * (xab*yp-yab*xp)
               ddtdxic = term2 * (ycb*zp-zcb*yp)
               ddtdyic = term2 * (zcb*xp-xcb*zp)
               ddtdzic = term2 * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the bond length deviations
c
               j = isb(2,istrbnd)
               k = isb(3,istrbnd)
               dr1 = rab - bl(j)
               term1 = 1.0d0 / rab
               dr2 = rcb - bl(k)
               term2 = 1.0d0 / rcb
               ddrdxia = term1 * xab
               ddrdyia = term1 * yab
               ddrdzia = term1 * zab
               ddrdxic = term2 * xcb
               ddrdyic = term2 * ycb
               ddrdzic = term2 * zcb
c
c     abbreviations used in defining chain rule terms
c
               term1 = stbnunit * force1
               term2 = stbnunit * force2
               termr = term1*dr1 + term2*dr2
               term1t = term1 * dt
               term2t = term2 * dt
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  termr = termr * fgrp
                  term1t = term1t * fgrp
                  term2t = term2t * fgrp
               end if
c
c     get the energy and master chain rule terms for derivatives
c
               e = termr * dt
               dedxia = term1t*ddrdxia + termr*ddtdxia
               dedyia = term1t*ddrdyia + termr*ddtdyia
               dedzia = term1t*ddrdzia + termr*ddtdzia
               dedxic = term2t*ddrdxic + termr*ddtdxic
               dedyic = term2t*ddrdyic + termr*ddtdyic
               dedzic = term2t*ddrdzic + termr*ddtdzic
               dedxib = -dedxia - dedxic
               dedyib = -dedyia - dedyic
               dedzib = -dedzia - dedzic
c
c     increment the total stretch-bend energy and derivatives
c
               ebao = ebao + e
               debao(1,ia) = debao(1,ia) + dedxia
               debao(2,ia) = debao(2,ia) + dedyia
               debao(3,ia) = debao(3,ia) + dedzia
               debao(1,ib) = debao(1,ib) + dedxib
               debao(2,ib) = debao(2,ib) + dedyib
               debao(3,ib) = debao(3,ib) + dedzib
               debao(1,ic) = debao(1,ic) + dedxic
               debao(2,ic) = debao(2,ic) + dedyic
               debao(3,ic) = debao(3,ic) + dedzic
c
c     increment the internal virial tensor components
c
               vxx = xab*dedxia + xcb*dedxic
               vyx = yab*dedxia + ycb*dedxic
               vzx = zab*dedxia + zcb*dedxic
               vyy = yab*dedyia + ycb*dedyic
               vzy = zab*dedyia + zcb*dedyic
               vzz = zab*dedzia + zcb*dedzic
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
      eba = ebao
      do i = 1, n
         deba(1,i) = debao(1,i)
         deba(2,i) = debao(2,i)
         deba(3,i) = debao(3,i)
      end do
      do i = 1, 3
         vir(1,i) = viro(1,i)
         vir(2,i) = viro(2,i)
         vir(3,i) = viro(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (debao)
      return
      end

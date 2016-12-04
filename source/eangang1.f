c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine eangang1  --  angle-angle energy & derivatives  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "eangang1" calculates the angle-angle potential energy and
c     first derivatives with respect to Cartesian coordinates
c
c
      subroutine eangang1
      use sizes
      use angang
      use angbnd
      use angpot
      use atoms
      use bound
      use deriv
      use energi
      use group
      use math
      use usage
      use virial
      implicit none
      integer i,k,iangang
      integer ia,ib,ic,id,ie
      real*8 e,eaao
      real*8 angle,fgrp
      real*8 dot,cosine
      real*8 dt1,deddt1
      real*8 dt2,deddt2
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xie,yie,zie
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xdb,ydb,zdb
      real*8 xeb,yeb,zeb
      real*8 rab2,rcb2
      real*8 rdb2,reb2
      real*8 xp,yp,zp,rp
      real*8 xq,yq,zq,rq
      real*8 terma,termc
      real*8 termd,terme
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 dedxie,dedyie,dedzie
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 viro(3,3)
      real*8, allocatable :: deaao(:,:)
      logical proceed
c
c
c     zero out the angle-angle energy and first derivatives
c
      eaa = 0.0d0
      do i = 1, n
         deaa(1,i) = 0.0d0
         deaa(2,i) = 0.0d0
         deaa(3,i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (deaao(3,n))
c
c     transfer global to local copies for OpenMP calculation
c
      eaao = eaa
      do i = 1, n
         deaao(1,i) = deaa(1,i)
         deaao(2,i) = deaa(2,i)
         deaao(3,i) = deaa(3,i)
      end do
      do i = 1, 3
         viro(1,i) = vir(1,i)
         viro(2,i) = vir(2,i)
         viro(3,i) = vir(3,i)
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nangang,iaa,iang,
!$OMP& use,x,y,z,anat,kaa,aaunit,use_group,use_polymer)
!$OMP& shared(eaao,deaao,viro)
!$OMP DO reduction(+:eaao,deaao,viro) schedule(guided)
c
c     find the energy of each angle-angle interaction
c
      do iangang = 1, nangang
         i = iaa(1,iangang)
         k = iaa(2,iangang)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(1,k)
         ie = iang(3,k)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,ie,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or. use(ic)
     &                               .or. use(id) .or. use(ie))
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
            xid = x(id)
            yid = y(id)
            zid = z(id)
            xie = x(ie)
            yie = y(ie)
            zie = z(ie)
c
c     compute the values of the two bond angles
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            xeb = xie - xib
            yeb = yie - yib
            zeb = zie - zib
            if (use_polymer) then
               call image (xab,yab,zab)
               call image (xcb,ycb,zcb)
               call image (xdb,ydb,zdb)
               call image (xeb,yeb,zeb)
            end if
            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            reb2 = xeb*xeb + yeb*yeb + zeb*zeb
            xp = ycb*zab - zcb*yab
            yp = zcb*xab - xcb*zab
            zp = xcb*yab - ycb*xab
            xq = yeb*zdb - zeb*ydb
            yq = zeb*xdb - xeb*zdb
            zq = xeb*ydb - yeb*xdb
            rp = sqrt(xp*xp + yp*yp + zp*zp)
            rq = sqrt(xq*xq + yq*yq + zq*zq)
            if (rp*rq .ne. 0.0d0) then
               dot = xab*xcb + yab*ycb + zab*zcb
               cosine = dot / sqrt(rab2*rcb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt1 = angle - anat(i)
               dot = xdb*xeb + ydb*yeb + zdb*zeb
               cosine = dot / sqrt(rdb2*reb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)
               dt2 = angle - anat(k)
c
c     get the energy and master chain rule terms for derivatives
c
               e = aaunit * kaa(iangang) * dt1 * dt2
               deddt1 = radian * e / dt1
               deddt2 = radian * e / dt2
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  deddt1 = deddt1 * fgrp
                  deddt2 = deddt2 * fgrp
               end if
c
c     find chain rule terms for the first bond angle deviation
c
               terma = -deddt1 / (rab2*rp)
               termc = deddt1 / (rcb2*rp)
               dedxia = terma * (yab*zp-zab*yp)
               dedyia = terma * (zab*xp-xab*zp)
               dedzia = terma * (xab*yp-yab*xp)
               dedxic = termc * (ycb*zp-zcb*yp)
               dedyic = termc * (zcb*xp-xcb*zp)
               dedzic = termc * (xcb*yp-ycb*xp)
c
c     find chain rule terms for the second bond angle deviation
c
               termd = -deddt2 / (rdb2*rq)
               terme = deddt2 / (reb2*rq)
               dedxid = termd * (ydb*zq-zdb*yq)
               dedyid = termd * (zdb*xq-xdb*zq)
               dedzid = termd * (xdb*yq-ydb*xq)
               dedxie = terme * (yeb*zq-zeb*yq)
               dedyie = terme * (zeb*xq-xeb*zq)
               dedzie = terme * (xeb*yq-yeb*xq)
c
c     get the central atom derivative terms by difference
c
               dedxib = -dedxia - dedxic - dedxid - dedxie
               dedyib = -dedyia - dedyic - dedyid - dedyie
               dedzib = -dedzia - dedzic - dedzid - dedzie
c
c     increment the total angle-angle energy and derivatives
c
               eaao = eaao + e
               deaao(1,ia) = deaao(1,ia) + dedxia
               deaao(2,ia) = deaao(2,ia) + dedyia
               deaao(3,ia) = deaao(3,ia) + dedzia
               deaao(1,ib) = deaao(1,ib) + dedxib
               deaao(2,ib) = deaao(2,ib) + dedyib
               deaao(3,ib) = deaao(3,ib) + dedzib
               deaao(1,ic) = deaao(1,ic) + dedxic
               deaao(2,ic) = deaao(2,ic) + dedyic
               deaao(3,ic) = deaao(3,ic) + dedzic
               deaao(1,id) = deaao(1,id) + dedxid
               deaao(2,id) = deaao(2,id) + dedyid
               deaao(3,id) = deaao(3,id) + dedzid
               deaao(1,ie) = deaao(1,ie) + dedxie
               deaao(2,ie) = deaao(2,ie) + dedyie
               deaao(3,ie) = deaao(3,ie) + dedzie
c
c     increment the internal virial tensor components
c
               vxx = xab*dedxia + xcb*dedxic + xdb*dedxid + xeb*dedxie
               vyx = yab*dedxia + ycb*dedxic + ydb*dedxid + yeb*dedxie
               vzx = zab*dedxia + zcb*dedxic + zdb*dedxid + zeb*dedxie
               vyy = yab*dedyia + ycb*dedyic + ydb*dedyid + yeb*dedyie
               vzy = zab*dedyia + zcb*dedyic + zdb*dedyid + zeb*dedyie
               vzz = zab*dedzia + zcb*dedzic + zdb*dedzid + zeb*dedzie
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
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     transfer local to global copies for OpenMP calculation
c
      eaa = eaao
      do i = 1, n
         deaa(1,i) = deaao(1,i)
         deaa(2,i) = deaao(2,i)
         deaa(3,i) = deaao(3,i)
      end do
      do i = 1, 3
         vir(1,i) = viro(1,i)
         vir(2,i) = viro(2,i)
         vir(3,i) = viro(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (deaao)
      return
      end

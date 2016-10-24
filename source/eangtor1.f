c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2014 by Chao Lu & Jay William Ponder  ##
c     ##                  All Rights Reserved                 ##
c     ##########################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine eangtor1  --  angle-torsion energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "eangtor1" calculates the angle-torsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c
      subroutine eangtor1
      use sizes
      use angbnd
      use angtor
      use atoms
      use bound
      use deriv
      use energi
      use group
      use math
      use torpot
      use tors
      use usage
      use virial
      implicit none
      integer i,k,iangtor
      integer ia,ib,ic,id
      real*8 e,e1,e2,eato
      real*8 rcb,fgrp
      real*8 ddt,dedphi
      real*8 rt2,ru2,rtru
      real*8 rba2,rcb2,rdc2
      real*8 dot,dt
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 v1,v2,v3
      real*8 c1,c2,c3
      real*8 s1,s2,s3
      real*8 sine,cosine
      real*8 sine2,cosine2
      real*8 sine3,cosine3
      real*8 phi1,phi2,phi3
      real*8 dphi1,dphi2,dphi3
      real*8 angle,cosang
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 terma,termb
      real*8 termc,termd
      real*8 dedxt,dedyt,dedzt
      real*8 dedxu,dedyu,dedzu
      real*8 dedxia,dedyia,dedzia
      real*8 dedxib,dedyib,dedzib
      real*8 dedxic,dedyic,dedzic
      real*8 dedxid,dedyid,dedzid
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 viro(3,3)
      real*8, allocatable :: deato(:,:)
      logical proceed
c
c
c     zero out the angle-torsion energy and first derivatives
c
      eat = 0.0d0
      do i = 1, n
         deat(1,i) = 0.0d0
         deat(2,i) = 0.0d0
         deat(3,i) = 0.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (deato(3,n))
c
c     transfer global to local copies for OpenMP calculation
c
      eato = eat
      do i = 1, n
         deato(1,i) = deat(1,i)
         deato(2,i) = deat(2,i)
         deato(3,i) = deat(3,i)
      end do
      do i = 1, 3
         viro(1,i) = vir(1,i)
         viro(2,i) = vir(2,i)
         viro(3,i) = vir(3,i)
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nangtor,iat,itors,kant,anat,
!$OMP& tors1,tors2,tors3,use,x,y,z,atorunit,use_group,use_polymer)
!$OMP& shared(eato,deato,viro)
!$OMP DO reduction(+:eato,deato,viro) schedule(guided)
c
c     calculate the angle-torsion energy and first derviatives
c
      do iangtor = 1, nangtor
         i = iat(1,iangtor)
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     compute the value of the torsional angle
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
            xba = xib - xia
            yba = yib - yia
            zba = zib - zia
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdc = xid - xic
            ydc = yid - yic
            zdc = zid - zic
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
            end if
            rba2 = xba*xba + yba*yba + zba*zba
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rdc2 = xdc*xdc + ydc*ydc + zdc*zdc
            if (min(rba2,rcb2,rdc2) .ne. 0.0d0) then
               xt = yba*zcb - ycb*zba
               yt = zba*xcb - zcb*xba
               zt = xba*ycb - xcb*yba
               xu = ycb*zdc - ydc*zcb
               yu = zcb*xdc - zdc*xcb
               zu = xcb*ydc - xdc*ycb
               xtu = yt*zu - yu*zt
               ytu = zt*xu - zu*xt
               ztu = xt*yu - xu*yt
               rt2 = xt*xt + yt*yt + zt*zt
               rt2 = max(rt2,0.000001d0)
               ru2 = xu*xu + yu*yu + zu*zu
               ru2 = max(ru2,0.000001d0)
               rtru = sqrt(rt2*ru2)
               xca = xic - xia
               yca = yic - yia
               zca = zic - zia
               xdb = xid - xib
               ydb = yid - yib
               zdb = zid - zib
               if (use_polymer) then
                  call image (xca,yca,zca)
                  call image (xdb,ydb,zdb)
               end if
               rcb = sqrt(rcb2)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     compute multiple angle trigonometry and phase terms
c
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
               cosine2 = cosine*cosine - sine*sine
               sine2 = 2.0d0 * cosine * sine
               cosine3 = cosine*cosine2 - sine*sine2
               sine3 = cosine*sine2 + sine*cosine2
               phi1 = 1.0d0 + (cosine*c1 + sine*s1)
               phi2 = 1.0d0 + (cosine2*c2 + sine2*s2)
               phi3 = 1.0d0 + (cosine3*c3 + sine3*s3)
               dphi1 = (cosine*s1 - sine*c1)
               dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0d0 * (cosine3*s3 - sine3*c3)
c
c     set the angle-torsion parameters for the first angle
c
               v1 = kant(1,iangtor)
               v2 = kant(2,iangtor)
               v3 = kant(3,iangtor)
               k = iat(2,iangtor)
               dot = xba*xcb + yba*ycb + zba*zcb
               cosang = -dot / sqrt(rba2*rcb2)
               angle = radian * acos(cosang)
               dt = angle - anat(k)
               e1 = atorunit * dt * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphi = atorunit * dt * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               ddt = atorunit * radian * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e1 = e1 * fgrp
                  dedphi = dedphi * fgrp
                  ddt = ddt * fgrp
               end if
c
c     compute derivative components for this interaction
c
               dedxt = dedphi * (zcb*yt-ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (xcb*zt-zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (ycb*xt-xcb*yt) / (rt2*rcb)
               dedxu = dedphi * (ycb*zu-zcb*yu) / (ru2*rcb)
               dedyu = dedphi * (zcb*xu-xcb*zu) / (ru2*rcb)
               dedzu = dedphi * (xcb*yu-ycb*xu) / (ru2*rcb)
c
c     increment chain rule components for the first angle
c
               terma = -ddt / (rba2*sqrt(rt2))
               termc = ddt / (rcb2*sqrt(rt2))
               dedxia = terma*(zba*yt-yba*zt) + zcb*dedyt - ycb*dedzt
               dedyia = terma*(xba*zt-zba*xt) + xcb*dedzt - zcb*dedxt
               dedzia = terma*(yba*xt-xba*yt) + ycb*dedxt - xcb*dedyt
               dedxib = terma*(yba*zt-zba*yt) + termc*(zcb*yt-ycb*zt)
     &                     + yca*dedzt - zca*dedyt
     &                     + zdc*dedyu - ydc*dedzu
               dedyib = terma*(zba*xt-xba*zt) + termc*(xcb*zt-zcb*xt)
     &                     + zca*dedxt - xca*dedzt
     &                     + xdc*dedzu - zdc*dedxu
               dedzib = terma*(xba*yt-yba*xt) + termc*(ycb*xt-xcb*yt)
     &                     + xca*dedyt - yca*dedxt
     &                     + ydc*dedxu - xdc*dedyu
               dedxic = termc*(ycb*zt-zcb*yt) + zba*dedyt
     &                     - yba*dedzt + ydb*dedzu - zdb*dedyu
               dedyic = termc*(zcb*xt-xcb*zt) + xba*dedzt
     &                     - zba*dedxt + zdb*dedxu - xdb*dedzu
               dedzic = termc*(xcb*yt-ycb*xt) + yba*dedxt
     &                     - xba*dedyt + xdb*dedyu - ydb*dedxu
               dedxid = zcb*dedyu - ycb*dedzu
               dedyid = xcb*dedzu - zcb*dedxu
               dedzid = ycb*dedxu - xcb*dedyu
c
c     get the angle-torsion values for the second angle
c
               v1 = kant(4,iangtor)
               v2 = kant(5,iangtor)
               v3 = kant(6,iangtor)
               k = iat(3,iangtor)
               dot = xcb*xdc + ycb*ydc + zcb*zdc
               cosang = -dot / sqrt(rcb2*rdc2)
               angle = radian * acos(cosang)
               dt = angle - anat(k)
               e2 = atorunit * dt * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphi = atorunit * dt * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               ddt = atorunit * radian * (v1*phi1 + v2*phi2 + v3*phi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e2 = e2 * fgrp
                  dedphi = dedphi * fgrp
                  ddt = ddt * fgrp
               end if
c
c     compute derivative components for this interaction
c
               dedxt = dedphi * (zcb*yt-ycb*zt) / (rt2*rcb)
               dedyt = dedphi * (xcb*zt-zcb*xt) / (rt2*rcb)
               dedzt = dedphi * (ycb*xt-xcb*yt) / (rt2*rcb)
               dedxu = dedphi * (ycb*zu-zcb*yu) / (ru2*rcb)
               dedyu = dedphi * (zcb*xu-xcb*zu) / (ru2*rcb)
               dedzu = dedphi * (xcb*yu-ycb*xu) / (ru2*rcb)
c
c     increment chain rule components for the second angle
c
               termb = -ddt / (rcb2*sqrt(ru2))
               termd = ddt / (rdc2*sqrt(ru2))
               dedxia = dedxia + zcb*dedyt - ycb*dedzt
               dedyia = dedyia + xcb*dedzt - zcb*dedxt
               dedzia = dedzia + ycb*dedxt - xcb*dedyt
               dedxib = dedxib + termb*(zcb*yu-ycb*zu) + yca*dedzt
     &                     - zca*dedyt + zdc*dedyu - ydc*dedzu
               dedyib = dedyib + termb*(xcb*zu-zcb*xu) + zca*dedxt
     &                     - xca*dedzt + xdc*dedzu - zdc*dedxu
               dedzib = dedzib + termb*(ycb*xu-xcb*yu) + xca*dedyt
     &                     - yca*dedxt + ydc*dedxu - xdc*dedyu
               dedxic = dedxic + termb*(ycb*zu-zcb*yu)
     &                     + termd*(zdc*yu-ydc*zu) + zba*dedyt
     &                     - yba*dedzt + ydb*dedzu - zdb*dedyu
               dedyic = dedyic + termb*(zcb*xu-xcb*zu)
     &                     + termd*(xdc*zu-zdc*xu) + xba*dedzt
     &                     - zba*dedxt + zdb*dedxu - xdb*dedzu
               dedzic = dedzic + termb*(xcb*yu-ycb*xu)
     &                     + termd*(ydc*xu-xdc*yu) + yba*dedxt
     &                     - xba*dedyt + xdb*dedyu - ydb*dedxu
               dedxid = dedxid + termd*(ydc*zu-zdc*yu)
     &                     + zcb*dedyu - ycb*dedzu
               dedyid = dedyid + termd*(zdc*xu-xdc*zu)
     &                     + xcb*dedzu - zcb*dedxu
               dedzid = dedzid + termd*(xdc*yu-ydc*xu)
     &                     + ycb*dedxu - xcb*dedyu
c
c     increment the angle-torsion energy and gradient
c
               e = e1 + e2
               eato = eato + e
               deato(1,ia) = deato(1,ia) + dedxia
               deato(2,ia) = deato(2,ia) + dedyia
               deato(3,ia) = deato(3,ia) + dedzia
               deato(1,ib) = deato(1,ib) + dedxib
               deato(2,ib) = deato(2,ib) + dedyib
               deato(3,ib) = deato(3,ib) + dedzib
               deato(1,ic) = deato(1,ic) + dedxic
               deato(2,ic) = deato(2,ic) + dedyic
               deato(3,ic) = deato(3,ic) + dedzic
               deato(1,id) = deato(1,id) + dedxid
               deato(2,id) = deato(2,id) + dedyid
               deato(3,id) = deato(3,id) + dedzid
c
c     increment the internal virial tensor components
c
               vxx = xcb*(dedxic+dedxid) - xba*dedxia + xdc*dedxid
               vyx = ycb*(dedxic+dedxid) - yba*dedxia + ydc*dedxid
               vzx = zcb*(dedxic+dedxid) - zba*dedxia + zdc*dedxid
               vyy = ycb*(dedyic+dedyid) - yba*dedyia + ydc*dedyid
               vzy = zcb*(dedyic+dedyid) - zba*dedyia + zdc*dedyid
               vzz = zcb*(dedzic+dedzid) - zba*dedzia + zdc*dedzid
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
      eat = eato
      do i = 1, n
         deat(1,i) = deato(1,i)
         deat(2,i) = deato(2,i)
         deat(3,i) = deato(3,i)
      end do
      do i = 1, 3
         vir(1,i) = viro(1,i)
         vir(2,i) = viro(2,i)
         vir(3,i) = viro(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (deato)
      return
      end

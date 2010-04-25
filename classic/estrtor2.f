c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine estrtor2  --  atomwise stretch-torsion Hessian  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "estrtor2" calculates the stretch-torsion potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine estrtor2 (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bond.i'
      include 'bound.i'
      include 'group.i'
      include 'hessn.i'
      include 'strtor.i'
      include 'torpot.i'
      include 'tors.i'
      integer i,j,k,istrtor
      integer ia,ib,ic,id
      real*8 dedphi,d2edphi2,fgrp
      real*8 rt2,ru2,rtru,rcb
      real*8 xt,yt,zt,xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 v1,v2,v3
      real*8 c1,c2,c3
      real*8 s1,s2,s3
      real*8 sine,cosine
      real*8 sine2,cosine2
      real*8 sine3,cosine3
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xba,yba,zba
      real*8 xcb,ycb,zcb
      real*8 xdc,ydc,zdc
      real*8 xca,yca,zca
      real*8 xdb,ydb,zdb
      real*8 phi1,phi2,phi3
      real*8 dphi1,dphi2,dphi3
      real*8 d2phi1,d2phi2,d2phi3
      real*8 dr,ddr,d2dr
      real*8 ddrdx,ddrdy,ddrdz
      real*8 d2drdxx,d2drdyy,d2drdzz
      real*8 d2drdxy,d2drdxz,d2drdyz
      real*8 dphidxt,dphidyt,dphidzt
      real*8 dphidxu,dphidyu,dphidzu
      real*8 dphidxia,dphidyia,dphidzia
      real*8 dphidxib,dphidyib,dphidzib
      real*8 dphidxic,dphidyic,dphidzic
      real*8 dphidxid,dphidyid,dphidzid
      real*8 xycb2,xzcb2,yzcb2
      real*8 rcbxt,rcbyt,rcbzt,rcbt2
      real*8 rcbxu,rcbyu,rcbzu,rcbu2
      real*8 dphidxibt,dphidyibt,dphidzibt
      real*8 dphidxibu,dphidyibu,dphidzibu
      real*8 dphidxict,dphidyict,dphidzict
      real*8 dphidxicu,dphidyicu,dphidzicu
      real*8 dxia,dyia,dzia
      real*8 dxib,dyib,dzib
      real*8 dxic,dyic,dzic
      real*8 dxid,dyid,dzid
      real*8 dxiaxia,dyiayia,dziazia
      real*8 dxibxib,dyibyib,dzibzib
      real*8 dxicxic,dyicyic,dziczic
      real*8 dxidxid,dyidyid,dzidzid
      real*8 dxiayia,dxiazia,dyiazia
      real*8 dxibyib,dxibzib,dyibzib
      real*8 dxicyic,dxiczic,dyiczic
      real*8 dxidyid,dxidzid,dyidzid
      real*8 dxiaxib,dxiayib,dxiazib
      real*8 dyiaxib,dyiayib,dyiazib
      real*8 dziaxib,dziayib,dziazib
      real*8 dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic
      real*8 dziaxic,dziayic,dziazic
      real*8 dxiaxid,dxiayid,dxiazid
      real*8 dyiaxid,dyiayid,dyiazid
      real*8 dziaxid,dziayid,dziazid
      real*8 dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic
      real*8 dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid
      real*8 dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid
      real*8 dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid
      real*8 dzicxid,dzicyid,dziczid
      logical proceed
c
c
c     compute the Hessian elements of the stretch-torsions
c
      do istrtor = 1, nstrtor
         j = ist(1,istrtor)
         ia = itors(1,j)
         ib = itors(2,j)
         ic = itors(3,j)
         id = itors(4,j)
c
c     decide whether to compute the current interaction
c
         proceed = (i.eq.ia .or. i.eq.ib .or. i.eq.ic .or. i.eq.id)
         if (proceed .and. use_group)
     &      call groups (proceed,fgrp,ia,ib,ic,id,0,0)
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
            ru2 = xu*xu + yu*yu + zu*zu
            rtru = sqrt(rt2 * ru2)
            if (rtru .ne. 0.0d0) then
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     set the stretch-torsional parameters for this angle
c
               v1 = kst(1,istrtor)
               c1 = tors1(3,i)
               s1 = tors1(4,i)
               v2 = kst(2,istrtor)
               c2 = tors2(3,i)
               s2 = tors2(4,i)
               v3 = kst(3,istrtor)
               c3 = tors3(3,i)
               s3 = tors3(4,i)
c
c     compute the multiple angle trigonometry and the phase terms
c
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
               d2phi1 = -(cosine*c1 + sine*s1)
               d2phi2 = -4.0d0 * (cosine2*c2 + sine2*s2)
               d2phi3 = -9.0d0 * (cosine3*c3 + sine3*s3)
c
c     calculate the bond-stretch for the central bond
c
               k = ist(2,istrtor)
               dr = rcb - bl(k)
c
c     calculate the stretch-torsion master chain rule terms
c
               dedphi = storunit * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2 = storunit * dr
     &                       * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3)
               ddr = 1.0d0 / rcb
               d2dr = -storunit * (v1*phi1 + v2*phi2 + v3*phi3) / rcb**3
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
                  d2dr = d2dr * fgrp
               end if
c
c     first and second derivative components for the bond stretch
c
               ddrdx = xcb * ddr
               ddrdy = ycb * ddr
               ddrdz = zcb * ddr
               d2drdxx = (xcb**2-rcb**2) * d2dr
               d2drdyy = (ycb**2-rcb**2) * d2dr
               d2drdzz = (zcb**2-rcb**2) * d2dr
               d2drdxy = xcb * ycb * d2dr
               d2drdxz = xcb * zcb * d2dr
               d2drdyz = ycb * zcb * d2dr
c
c     abbreviations for first derivative chain rule terms
c
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
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb)
c
c     abbreviations for second derivative chain rule terms
c
               xycb2 = xcb*xcb + ycb*ycb
               xzcb2 = xcb*xcb + zcb*zcb
               yzcb2 = ycb*ycb + zcb*zcb
               rcbxt = -2.0d0 * rcb * dphidxt
               rcbyt = -2.0d0 * rcb * dphidyt
               rcbzt = -2.0d0 * rcb * dphidzt
               rcbt2 = rcb * rt2
               rcbxu = 2.0d0 * rcb * dphidxu
               rcbyu = 2.0d0 * rcb * dphidyu
               rcbzu = 2.0d0 * rcb * dphidzu
               rcbu2 = rcb * ru2
               dphidxibt = yca*dphidzt - zca*dphidyt
               dphidxibu = zdc*dphidyu - ydc*dphidzu
               dphidyibt = zca*dphidxt - xca*dphidzt
               dphidyibu = xdc*dphidzu - zdc*dphidxu
               dphidzibt = xca*dphidyt - yca*dphidxt
               dphidzibu = ydc*dphidxu - xdc*dphidyu
               dphidxict = zba*dphidyt - yba*dphidzt
               dphidxicu = ydb*dphidzu - zdb*dphidyu
               dphidyict = xba*dphidzt - zba*dphidxt
               dphidyicu = zdb*dphidxu - xdb*dphidzu
               dphidzict = yba*dphidxt - xba*dphidyt
               dphidzicu = xdb*dphidyu - ydb*dphidxu
c
c     intermediate terms for first derivative components
c
               dphidxia = zcb*dphidyt - ycb*dphidzt
               dphidyia = xcb*dphidzt - zcb*dphidxt
               dphidzia = ycb*dphidxt - xcb*dphidyt
               dphidxib = dphidxibt + dphidxibu
               dphidyib = dphidyibt + dphidyibu
               dphidzib = dphidzibt + dphidzibu
               dphidxic = dphidxict + dphidxicu
               dphidyic = dphidyict + dphidyicu
               dphidzic = dphidzict + dphidzicu
               dphidxid = zcb*dphidyu - ycb*dphidzu
               dphidyid = xcb*dphidzu - zcb*dphidxu
               dphidzid = ycb*dphidxu - xcb*dphidyu
c
c     chain rule terms for first derivative components
c
               dxia = dedphi * dphidxia
               dyia = dedphi * dphidyia
               dzia = dedphi * dphidzia
               dxib = dedphi * dphidxib
               dyib = dedphi * dphidyib
               dzib = dedphi * dphidzib
               dxic = dedphi * dphidxic
               dyic = dedphi * dphidyic
               dzic = dedphi * dphidzic
               dxid = dedphi * dphidxid
               dyid = dedphi * dphidyid
               dzid = dedphi * dphidzid
               dedphi = dedphi * dr
c
c     chain rule terms for second derivative components
c
               dxiaxia = rcbxt*dphidxia
               dxiayia = rcbxt*dphidyia - zcb*rcb/rt2
               dxiazia = rcbxt*dphidzia + ycb*rcb/rt2
               dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2
               dxiayic = rcbxt*dphidyict - dphidzt
     &                      - (xba*zcb*xcb+zba*yzcb2)/rcbt2
               dxiazic = rcbxt*dphidzict + dphidyt
     &                      + (xba*ycb*xcb+yba*yzcb2)/rcbt2
               dxiaxid = 0.0d0
               dxiayid = 0.0d0
               dxiazid = 0.0d0
               dyiayia = rcbyt*dphidyia
               dyiazia = rcbyt*dphidzia - xcb*rcb/rt2
               dyiaxib = rcbyt*dphidxibt - dphidzt
     &                      - (yca*zcb*ycb+zca*xzcb2)/rcbt2
               dyiaxic = rcbyt*dphidxict + dphidzt
     &                      + (yba*zcb*ycb+zba*xzcb2)/rcbt2
               dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2
               dyiazic = rcbyt*dphidzict - dphidxt
     &                      - (yba*xcb*ycb+xba*xzcb2)/rcbt2
               dyiaxid = 0.0d0
               dyiayid = 0.0d0
               dyiazid = 0.0d0
               dziazia = rcbzt*dphidzia
               dziaxib = rcbzt*dphidxibt + dphidyt
     &                      + (zca*ycb*zcb+yca*xycb2)/rcbt2
               dziayib = rcbzt*dphidyibt - dphidxt
     &                      - (zca*xcb*zcb+xca*xycb2)/rcbt2
               dziaxic = rcbzt*dphidxict - dphidyt
     &                      - (zba*ycb*zcb+yba*xycb2)/rcbt2
               dziayic = rcbzt*dphidyict + dphidxt
     &                      + (zba*xcb*zcb+xba*xycb2)/rcbt2
               dziazic = rcbzt*dphidzict + zcb*zt/rcbt2
               dziaxid = 0.0d0
               dziayid = 0.0d0
               dziazid = 0.0d0
               dxibxic = -xcb*dphidxib/(rcb*rcb)
     &             - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidxibt/rt2
     &             - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidxibu/ru2
               dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu
     &             - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidxibt/rt2
     &             + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidxibu/ru2
               dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2
               dxibyid = rcbyu*dphidxibu - dphidzu
     &                      - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2
               dxibzid = rcbzu*dphidxibu + dphidyu
     &                      + (zdc*ycb*zcb+ydc*xycb2)/rcbu2
               dyibzib = ycb*dphidzib/(rcb*rcb)
     &             - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
     &             - 2.0d0*(xt*zca-xca*zt)*dphidzibt/rt2
     &             + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
     &             + 2.0d0*(xu*zdc-xdc*zu)*dphidzibu/ru2
               dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu
     &             + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidyibt/rt2
     &             - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidyibu/ru2
               dyibyic = -ycb*dphidyib/(rcb*rcb)
     &             - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
     &             - 2.0d0*(zt*xba-zba*xt)*dphidyibt/rt2
     &             - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
     &             + 2.0d0*(zu*xdb-zdb*xu)*dphidyibu/ru2
               dyibxid = rcbxu*dphidyibu + dphidzu
     &                      + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2
               dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2
               dyibzid = rcbzu*dphidyibu - dphidxu
     &                      - (zdc*xcb*zcb+xdc*xycb2)/rcbu2
               dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu
     &             - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
     &             - 2.0d0*(yt*zba-yba*zt)*dphidzibt/rt2
     &             + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
     &             + 2.0d0*(yu*zdb-ydb*zu)*dphidzibu/ru2
               dzibzic = -zcb*dphidzib/(rcb*rcb)
     &             - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
     &             - 2.0d0*(xt*yba-xba*yt)*dphidzibt/rt2
     &             - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
     &             + 2.0d0*(xu*ydb-xdb*yu)*dphidzibu/ru2
               dzibxid = rcbxu*dphidzibu - dphidyu
     &                      - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2
               dzibyid = rcbyu*dphidzibu + dphidxu
     &                      + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2
               dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2
               dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2
               dxicyid = rcbyu*dphidxicu + dphidzu
     &                      + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2
               dxiczid = rcbzu*dphidxicu - dphidyu
     &                      - (zdb*ycb*zcb+ydb*xycb2)/rcbu2
               dyicxid = rcbxu*dphidyicu - dphidzu
     &                      - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2
               dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2
               dyiczid = rcbzu*dphidyicu + dphidxu
     &                      + (zdb*xcb*zcb+xdb*xycb2)/rcbu2
               dzicxid = rcbxu*dphidzicu + dphidyu
     &                      + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2
               dzicyid = rcbyu*dphidzicu - dphidxu
     &                      - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2
               dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2
               dxidxid = rcbxu*dphidxid
               dxidyid = rcbxu*dphidyid + zcb*rcb/ru2
               dxidzid = rcbxu*dphidzid - ycb*rcb/ru2
               dyidyid = rcbyu*dphidyid
               dyidzid = rcbyu*dphidzid + xcb*rcb/ru2
               dzidzid = rcbzu*dphidzid
c
c     get some second derivative chain rule terms by difference
c
               dxiaxib = -dxiaxia - dxiaxic - dxiaxid
               dxiayib = -dxiayia - dxiayic - dxiayid
               dxiazib = -dxiazia - dxiazic - dxiazid
               dyiayib = -dyiayia - dyiayic - dyiayid
               dyiazib = -dyiazia - dyiazic - dyiazid
               dziazib = -dziazia - dziazic - dziazid
               dxibxib = -dxiaxib - dxibxic - dxibxid
               dxibyib = -dyiaxib - dxibyic - dxibyid
               dxibzib = -dxiazib - dzibxic - dzibxid
               dxibzic = -dziaxib - dxibzib - dxibzid
               dyibyib = -dyiayib - dyibyic - dyibyid
               dyibzic = -dziayib - dyibzib - dyibzid
               dzibzib = -dziazib - dzibzic - dzibzid
               dzibyic = -dyiazib - dyibzib - dzibyid
               dxicxic = -dxiaxic - dxibxic - dxicxid
               dxicyic = -dyiaxic - dyibxic - dxicyid
               dxiczic = -dziaxic - dzibxic - dxiczid
               dyicyic = -dyiayic - dyibyic - dyicyid
               dyiczic = -dziayic - dzibyic - dyiczid
               dziczic = -dziazic - dzibzic - dziczid
c
c     increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                             + d2edphi2*dphidxia*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                             + d2edphi2*dphidyia*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                             + d2edphi2*dphidzia*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxia*dphidxib
     &                             - dxia*ddrdx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                             + d2edphi2*dphidyia*dphidxib
     &                             - dyia*ddrdx
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                             + d2edphi2*dphidzia*dphidxib
     &                             - dzia*ddrdx
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                             + d2edphi2*dphidxia*dphidyib
     &                             - dxia*ddrdy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                             + d2edphi2*dphidyia*dphidyib
     &                             - dyia*ddrdy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                             + d2edphi2*dphidzia*dphidyib
     &                             - dzia*ddrdy
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                             + d2edphi2*dphidxia*dphidzib
     &                             - dxia*ddrdz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                             + d2edphi2*dphidyia*dphidzib
     &                             - dyia*ddrdz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                             + d2edphi2*dphidzia*dphidzib
     &                             - dzia*ddrdz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxia*dphidxic
     &                             + dxia*ddrdx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                             + d2edphi2*dphidyia*dphidxic
     &                             + dyia*ddrdx
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                             + d2edphi2*dphidzia*dphidxic
     &                             + dzia*ddrdx
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                             + d2edphi2*dphidxia*dphidyic
     &                             + dxia*ddrdy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                             + d2edphi2*dphidyia*dphidyic
     &                             + dyia*ddrdy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                             + d2edphi2*dphidzia*dphidyic
     &                             + dzia*ddrdy
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                             + d2edphi2*dphidxia*dphidzic
     &                             + dxia*ddrdz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                             + d2edphi2*dphidyia*dphidzic
     &                             + dyia*ddrdz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                             + d2edphi2*dphidzia*dphidzic
     &                             + dzia*ddrdz
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxia*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                             + d2edphi2*dphidyia*dphidxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                             + d2edphi2*dphidzia*dphidxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                             + d2edphi2*dphidxia*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                             + d2edphi2*dphidyia*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                             + d2edphi2*dphidzia*dphidyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                             + d2edphi2*dphidxia*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                             + d2edphi2*dphidyia*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                             + d2edphi2*dphidzia*dphidzid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                             + d2edphi2*dphidxib*dphidxib
     &                             - 2.0d0*dxib*ddrdx + d2drdxx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             - dyib*ddrdx - dxib*ddrdy + d2drdxy
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             - dzib*ddrdx - dxib*ddrdz + d2drdxz
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             - dxib*ddrdy - dyib*ddrdx + d2drdxy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                             + d2edphi2*dphidyib*dphidyib
     &                             - 2.0d0*dyib*ddrdy + d2drdyy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             - dzib*ddrdy - dyib*ddrdz + d2drdyz
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             - dxib*ddrdz - dzib*ddrdx + d2drdxz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             - dyib*ddrdz - dzib*ddrdy + d2drdyz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                             + d2edphi2*dphidzib*dphidzib
     &                             - 2.0d0*dzib*ddrdz + d2drdzz
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxib*dphidxia
     &                             - dxia*ddrdx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                             + d2edphi2*dphidyib*dphidxia
     &                             - dxia*ddrdy
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                             + d2edphi2*dphidzib*dphidxia
     &                             - dxia*ddrdz
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                             + d2edphi2*dphidxib*dphidyia
     &                             - dyia*ddrdx
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                             + d2edphi2*dphidyib*dphidyia
     &                             - dyia*ddrdy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                             + d2edphi2*dphidzib*dphidyia
     &                             - dyia*ddrdz
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                             + d2edphi2*dphidxib*dphidzia
     &                             - dzia*ddrdx
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                             + d2edphi2*dphidyib*dphidzia
     &                             - dzia*ddrdy
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                             + d2edphi2*dphidzib*dphidzia
     &                             - dzia*ddrdz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                             + d2edphi2*dphidxib*dphidxic
     &                             + (dxib-dxic)*ddrdx - d2drdxx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                             + d2edphi2*dphidyib*dphidxic
     &                             + dyib*ddrdx - dxic*ddrdy - d2drdxy
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                             + d2edphi2*dphidzib*dphidxic
     &                             + dzib*ddrdx - dxic*ddrdz - d2drdxz
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                             + d2edphi2*dphidxib*dphidyic
     &                             + dxib*ddrdy - dyic*ddrdx - d2drdxy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                             + d2edphi2*dphidyib*dphidyic
     &                             + (dyib-dyic)*ddrdy - d2drdyy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                             + d2edphi2*dphidzib*dphidyic
     &                             + dzib*ddrdy - dyic*ddrdz - d2drdyz
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                             + d2edphi2*dphidxib*dphidzic
     &                             + dxib*ddrdz - dzic*ddrdx - d2drdxz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                             + d2edphi2*dphidyib*dphidzic
     &                             + dyib*ddrdz - dzic*ddrdy - d2drdyz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                             + d2edphi2*dphidzib*dphidzic
     &                             + (dzib-dzic)*ddrdz - d2drdzz
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                             + d2edphi2*dphidxib*dphidxid
     &                             - dxid*ddrdx
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                             + d2edphi2*dphidyib*dphidxid
     &                             - dxid*ddrdy
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                             + d2edphi2*dphidzib*dphidxid
     &                             - dxid*ddrdz
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                             + d2edphi2*dphidxib*dphidyid
     &                             - dyid*ddrdx
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                             + d2edphi2*dphidyib*dphidyid
     &                             - dyid*ddrdy
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                             + d2edphi2*dphidzib*dphidyid
     &                             - dyid*ddrdz
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                             + d2edphi2*dphidxib*dphidzid
     &                             - dzid*ddrdx
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                             + d2edphi2*dphidyib*dphidzid
     &                             - dzid*ddrdy
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                             + d2edphi2*dphidzib*dphidzid
     &                             - dzid*ddrdz
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                             + d2edphi2*dphidxic*dphidxic
     &                             + 2.0d0*dxic*ddrdx + d2drdxx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + dyic*ddrdx + dxic*ddrdy + d2drdxy
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + dzic*ddrdx + dxic*ddrdz + d2drdxz
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + dxic*ddrdy + dyic*ddrdx + d2drdxy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                             + d2edphi2*dphidyic*dphidyic
     &                             + 2.0d0*dyic*ddrdy + d2drdyy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + dzic*ddrdy + dyic*ddrdz + d2drdyz
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + dxic*ddrdz + dzic*ddrdx + d2drdxz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + dyic*ddrdz + dzic*ddrdy + d2drdyz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                             + d2edphi2*dphidzic*dphidzic
     &                             + 2.0d0*dzic*ddrdz + d2drdzz
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxic*dphidxia
     &                             + dxia*ddrdx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                             + d2edphi2*dphidyic*dphidxia
     &                             + dxia*ddrdy
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                             + d2edphi2*dphidzic*dphidxia
     &                             + dxia*ddrdz
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                             + d2edphi2*dphidxic*dphidyia
     &                             + dyia*ddrdx
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                             + d2edphi2*dphidyic*dphidyia
     &                             + dyia*ddrdy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                             + d2edphi2*dphidzic*dphidyia
     &                             + dyia*ddrdz
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                             + d2edphi2*dphidxic*dphidzia
     &                             + dzia*ddrdx
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                             + d2edphi2*dphidyic*dphidzia
     &                             + dzia*ddrdy
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                             + d2edphi2*dphidzic*dphidzia
     &                             + dzia*ddrdz
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                             + d2edphi2*dphidxic*dphidxib
     &                             - (dxic-dxib)*ddrdx - d2drdxx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                             + d2edphi2*dphidyic*dphidxib
     &                             - dyic*ddrdx + dxib*ddrdy - d2drdxy
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                             + d2edphi2*dphidzic*dphidxib
     &                             - dzic*ddrdx + dxib*ddrdz - d2drdxz
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                             + d2edphi2*dphidxic*dphidyib
     &                             - dxic*ddrdy + dyib*ddrdx - d2drdxy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                             + d2edphi2*dphidyic*dphidyib
     &                             - (dyic-dyib)*ddrdy - d2drdyy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                             + d2edphi2*dphidzic*dphidyib
     &                             - dzic*ddrdy + dyib*ddrdz - d2drdyz
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                             + d2edphi2*dphidxic*dphidzib
     &                             - dxic*ddrdz + dzib*ddrdx - d2drdxz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                             + d2edphi2*dphidyic*dphidzib
     &                             - dyic*ddrdz + dzib*ddrdy - d2drdyz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                             + d2edphi2*dphidzic*dphidzib
     &                             - (dzic-dzib)*ddrdz - d2drdzz
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                             + d2edphi2*dphidxic*dphidxid
     &                             + dxid*ddrdx
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                             + d2edphi2*dphidyic*dphidxid
     &                             + dxid*ddrdy
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                             + d2edphi2*dphidzic*dphidxid
     &                             + dxid*ddrdz
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                             + d2edphi2*dphidxic*dphidyid
     &                             + dyid*ddrdx
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                             + d2edphi2*dphidyic*dphidyid
     &                             + dyid*ddrdy
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                             + d2edphi2*dphidzic*dphidyid
     &                             + dyid*ddrdz
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                             + d2edphi2*dphidxic*dphidzid
     &                             + dzid*ddrdx
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                             + d2edphi2*dphidyic*dphidzid
     &                             + dzid*ddrdy
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                             + d2edphi2*dphidzic*dphidzid
     &                             + dzid*ddrdz
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                             + d2edphi2*dphidxid*dphidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                             + d2edphi2*dphidyid*dphidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                             + d2edphi2*dphidzid*dphidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxid*dphidxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                             + d2edphi2*dphidyid*dphidxia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                             + d2edphi2*dphidzid*dphidxia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                             + d2edphi2*dphidxid*dphidyia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                             + d2edphi2*dphidyid*dphidyia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                             + d2edphi2*dphidzid*dphidyia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                             + d2edphi2*dphidxid*dphidzia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                             + d2edphi2*dphidyid*dphidzia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                             + d2edphi2*dphidzid*dphidzia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                             + d2edphi2*dphidxid*dphidxib
     &                             - dxid*ddrdx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                             + d2edphi2*dphidyid*dphidxib
     &                             - dyid*ddrdx
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                             + d2edphi2*dphidzid*dphidxib
     &                             - dzid*ddrdx
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                             + d2edphi2*dphidxid*dphidyib
     &                             - dxid*ddrdy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                             + d2edphi2*dphidyid*dphidyib
     &                             - dyid*ddrdy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                             + d2edphi2*dphidzid*dphidyib
     &                             - dzid*ddrdy
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                             + d2edphi2*dphidxid*dphidzib
     &                             - dxid*ddrdz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                             + d2edphi2*dphidyid*dphidzib
     &                             - dyid*ddrdz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                             + d2edphi2*dphidzid*dphidzib
     &                             - dzid*ddrdz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                             + d2edphi2*dphidxid*dphidxic
     &                             + dxid*ddrdx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                             + d2edphi2*dphidyid*dphidxic
     &                             + dyid*ddrdx
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                             + d2edphi2*dphidzid*dphidxic
     &                             + dzid*ddrdx
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                             + d2edphi2*dphidxid*dphidyic
     &                             + dxid*ddrdy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                             + d2edphi2*dphidyid*dphidyic
     &                             + dyid*ddrdy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                             + d2edphi2*dphidzid*dphidyic
     &                             + dzid*ddrdy
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                             + d2edphi2*dphidxid*dphidzic
     &                             + dxid*ddrdz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                             + d2edphi2*dphidyid*dphidzic
     &                             + dyid*ddrdz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                             + d2edphi2*dphidzid*dphidzic
     &                             + dzid*ddrdz
               end if
            end if
         end if
      end do
      return
      end

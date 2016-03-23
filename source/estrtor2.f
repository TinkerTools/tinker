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
      use sizes
      use atoms
      use bndstr
      use bound
      use group
      use hessn
      use strtor
      use torpot
      use tors
      implicit none
      integer i,j,k,istrtor
      integer ia,ib,ic,id
      real*8 fgrp,dedphi
      real*8 d2edphi2
      real*8 rt2,ru2,rtru
      real*8 rba,rcb,rdc
      real*8 xt,yt,zt
      real*8 xu,yu,zu
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
      real*8 force,d2drprime
      real*8 dedphiprime
      real*8 d2edphi2prime
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
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb)
               rba = sqrt(xba*xba + yba*yba + zba*zba)
               rdc = sqrt(xdc*xdc + ydc*ydc + zdc*zdc)
               cosine = (xt*xu + yt*yu + zt*zu) / rtru
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru)
c
c     abbreviations for first derivative chain rule terms
c
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
c     chain rule terms for torsion second derivative components
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
               dphi1 = cosine*s1 - sine*c1
               dphi2 = 2.0d0 * (cosine2*s2 - sine2*c2)
               dphi3 = 3.0d0 * (cosine3*s3 - sine3*c3)
               d2phi1 = -cosine*c1 - sine*s1
               d2phi2 = -4.0d0 * (cosine2*c2 + sine2*s2)
               d2phi3 = -9.0d0 * (cosine3*c3 + sine3*s3)
c
c     get the stretch-torsion values for the first bond
c
               v1 = kst(1,istrtor)
               v2 = kst(2,istrtor)
               v3 = kst(3,istrtor)
               k = ist(2,istrtor)
               dr = rba - bl(k)
               force = bk(k)
               dedphi = storunit * 2.0d0 * force 
     &                     * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2 = storunit * 2.0d0 * force * dr
     &                       * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3)
               ddr = 1.0d0 / rba
               d2dr = -storunit * 2.0d0 * force 
     &                    * (v1*phi1 + v2*phi2 + v3*phi3) / rba**3
               dedphiprime = storunit * 2.0d0 * force
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2prime = storunit * 2.0d0 * force 
     &                            * (v1*phi1 + v2*phi2 + v3*phi3)
     &                            * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3)
               d2drprime = storunit * 2.0d0 * force
     &                        * (v1*phi1 + v2*phi2 + v3*phi3)
     &                        * (v1*dphi1 + v2*dphi2 + v3*dphi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
                  d2dr = d2dr * fgrp
                  dedphiprime = dedphiprime * fgrp
                  d2edphi2prime = d2edphi2prime * fgrp
                  d2drprime = d2drprime * fgrp
               end if
c
c     compute derivative components for this interaction
c
               ddrdx = xba * ddr
               ddrdy = yba * ddr
               ddrdz = zba * ddr
               d2drdxx = (xba*xba-rba*rba) * d2dr
               d2drdyy = (yba*yba-rba*rba) * d2dr
               d2drdzz = (zba*zba-rba*rba) * d2dr
               d2drdxy = xba * yba * d2dr
               d2drdxz = xba * zba * d2dr
               d2drdyz = yba * zba * d2dr
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
c     increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                             + d2edphi2*dphidxia*dphidxia
     &                             - 2.0d0*dxia*ddrdx
     &                             + dedphiprime*dphidxia*dphidxia
     &                             + d2edphi2prime*dphidxia*dphidxia
     &                             + d2drprime*dxiaxia + d2drdxx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             - dxia*ddrdy - dyia*ddrdx
     &                             + dedphiprime*dphidxia*dphidyia
     &                             + d2edphi2prime*dphidxia*dphidyia
     &                             + d2drprime*dxiayia + d2drdxy
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             - dxia*ddrdz - dzia*ddrdx
     &                             + dedphiprime*dphidxia*dphidzia
     &                             + d2edphi2prime*dphidxia*dphidzia
     &                             + d2drprime*dxiazia + d2drdxz
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             - dyia*ddrdx - dxia*ddrdy
     &                             + dedphiprime*dphidxia*dphidyia
     &                             + d2edphi2prime*dphidxia*dphidyia
     &                             + d2drprime*dxiayia + d2drdxy
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                             + d2edphi2*dphidyia*dphidyia
     &                             - 2.0d0*dyia*ddrdy
     &                             + dedphiprime*dphidyia*dphidyia
     &                             + d2edphi2prime*dphidyia*dphidyia
     &                             + d2drprime*dyiayia + d2drdyy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             - dyia*ddrdz - dzia*ddrdy
     &                             + dedphiprime*dphidyia*dphidzia
     &                             + d2edphi2prime*dphidyia*dphidzia
     &                             + d2drprime*dyiazia + d2drdyz
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             - dxia*ddrdz - dzia*ddrdx
     &                             + dedphiprime*dphidxia*dphidzia
     &                             + d2edphi2prime*dphidxia*dphidzia
     &                             + d2drprime*dxiazia + d2drdxz
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             - dyia*ddrdz - dzia*ddrdy
     &                             + dedphiprime*dphidyia*dphidzia
     &                             + d2edphi2prime*dphidyia*dphidzia
     &                             + d2drprime*dyiazia + d2drdyz
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                             + d2edphi2*dphidzia*dphidzia
     &                             - 2.0d0*dzia*ddrdz
     &                             + dedphiprime*dphidzia*dphidzia
     &                             + d2edphi2prime*dphidzia*dphidzia
     &                             + d2drprime*dziazia + d2drdzz
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxia*dphidxib
     &                             + dxia*ddrdx - dxib*ddrdx
     &                             + dedphiprime*dphidxia*dphidxib
     &                             + d2edphi2prime*dphidxia*dphidxib
     &                             + d2drprime*dxiaxib - d2drdxx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                             + d2edphi2*dphidyia*dphidxib
     &                             + dyia*ddrdx - dxib*ddrdy
     &                             + dedphiprime*dphidyia*dphidxib
     &                             + d2edphi2prime*dphidyia*dphidxib
     &                             + d2drprime*dyiaxib - d2drdxy
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                             + d2edphi2*dphidzia*dphidxib
     &                             + dzia*ddrdx - dxib*ddrdz
     &                             + dedphiprime*dphidzia*dphidxib
     &                             + d2edphi2prime*dphidzia*dphidxib
     &                             + d2drprime*dziaxib - d2drdxz
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                             + d2edphi2*dphidxia*dphidyib
     &                             + dxia*ddrdy - dyib*ddrdx
     &                             + dedphiprime*dphidxia*dphidyib
     &                             + d2edphi2prime*dphidxia*dphidyib
     &                             + d2drprime*dxiayib - d2drdxy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                             + d2edphi2*dphidyia*dphidyib
     &                             + dyia*ddrdy - dyib*ddrdy
     &                             + dedphiprime*dphidyia*dphidyib
     &                             + d2edphi2prime*dphidyia*dphidyib
     &                             + d2drprime*dyiayib - d2drdyy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                             + d2edphi2*dphidzia*dphidyib
     &                             + dzia*ddrdy - dyib*ddrdz
     &                             + dedphiprime*dphidzia*dphidyib
     &                             + d2edphi2prime*dphidzia*dphidyib
     &                             + d2drprime*dziayib - d2drdyz
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                             + d2edphi2*dphidxia*dphidzib
     &                             + dxia*ddrdz - dzib*ddrdx
     &                             + dedphiprime*dphidxia*dphidzib
     &                             + d2edphi2prime*dphidxia*dphidzib
     &                             + d2drprime*dxiazib - d2drdxz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                             + d2edphi2*dphidyia*dphidzib
     &                             + dyia*ddrdz - dzib*ddrdy
     &                             + dedphiprime*dphidyia*dphidzib
     &                             + d2edphi2prime*dphidyia*dphidzib
     &                             + d2drprime*dyiazib - d2drdyz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                             + d2edphi2*dphidzia*dphidzib
     &                             + dzia*ddrdz - dzib*ddrdz
     &                             + dedphiprime*dphidzia*dphidzib
     &                             + d2edphi2prime*dphidzia*dphidzib
     &                             + d2drprime*dziazib - d2drdzz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxia*dphidxic
     &                             + dedphiprime*dphidxia*dphidxic
     &                             + d2edphi2prime*dphidxia*dphidxic
     &                             + d2drprime*dxiaxic - dxic*ddrdx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                             + d2edphi2*dphidyia*dphidxic
     &                             + dedphiprime*dphidyia*dphidxic
     &                             + d2edphi2prime*dphidyia*dphidxic
     &                             + d2drprime*dyiaxic - dxic*ddrdy
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                             + d2edphi2*dphidzia*dphidxic
     &                             + dedphiprime*dphidzia*dphidxic
     &                             + d2edphi2prime*dphidzia*dphidxic
     &                             + d2drprime*dziaxic - dxic*ddrdz
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                             + d2edphi2*dphidxia*dphidyic
     &                             + dedphiprime*dphidxia*dphidyic
     &                             + d2edphi2prime*dphidxia*dphidyic
     &                             + d2drprime*dxiayic - dyic*ddrdx
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                             + d2edphi2*dphidyia*dphidyic
     &                             + dedphiprime*dphidyia*dphidyic
     &                             + d2edphi2prime*dphidyia*dphidyic
     &                             + d2drprime*dyiayic - dyic*ddrdy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                             + d2edphi2*dphidzia*dphidyic
     &                             + dedphiprime*dphidzia*dphidyic
     &                             + d2edphi2prime*dphidzia*dphidyic
     &                             + d2drprime*dziayic - dyic*ddrdz
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                             + d2edphi2*dphidxia*dphidzic
     &                             + dedphiprime*dphidxia*dphidzic
     &                             + d2edphi2prime*dphidxia*dphidzic
     &                             + d2drprime*dxiazic - dzic*ddrdx
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                             + d2edphi2*dphidyia*dphidzic
     &                             + dedphiprime*dphidyia*dphidzic
     &                             + d2edphi2prime*dphidyia*dphidzic
     &                             + d2drprime*dyiazic - dzic*ddrdy
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                             + d2edphi2*dphidzia*dphidzic
     &                             + dedphiprime*dphidzia*dphidzic
     &                             + d2edphi2prime*dphidzia*dphidzic
     &                             + d2drprime*dziazic - dzic*ddrdz
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxia*dphidxid
     &                             + dedphiprime*dphidxia*dphidxid
     &                             + d2edphi2prime*dphidxia*dphidxid
     &                             + d2drprime*dxiaxid - dxid*ddrdx
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                             + d2edphi2*dphidyia*dphidxid
     &                             + dedphiprime*dphidyia*dphidxid
     &                             + d2edphi2prime*dphidyia*dphidxid
     &                             + d2drprime*dyiaxid - dxid*ddrdy
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                             + d2edphi2*dphidzia*dphidxid
     &                             + dedphiprime*dphidzia*dphidxid
     &                             + d2edphi2prime*dphidzia*dphidxid
     &                             + d2drprime*dziaxid - dxid*ddrdz
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                             + d2edphi2*dphidxia*dphidyid
     &                             + dedphiprime*dphidxia*dphidyid
     &                             + d2edphi2prime*dphidxia*dphidyid
     &                             + d2drprime*dxiayid - dyid*ddrdx
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                             + d2edphi2*dphidyia*dphidyid
     &                             + dedphiprime*dphidyia*dphidyid
     &                             + d2edphi2prime*dphidyia*dphidyid
     &                             + d2drprime*dyiayid - dyid*ddrdy
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                             + d2edphi2*dphidzia*dphidyid
     &                             + dedphiprime*dphidzia*dphidyid
     &                             + d2edphi2prime*dphidzia*dphidyid
     &                             + d2drprime*dziayid - dyid*ddrdz
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                             + d2edphi2*dphidxia*dphidzid
     &                             + dedphiprime*dphidxia*dphidzid
     &                             + d2edphi2prime*dphidxia*dphidzid
     &                             + d2drprime*dxiazid - dzid*ddrdx
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                             + d2edphi2*dphidyia*dphidzid
     &                             + dedphiprime*dphidyia*dphidzid
     &                             + d2edphi2prime*dphidyia*dphidzid
     &                             + d2drprime*dyiazid - dzid*ddrdy
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                             + d2edphi2*dphidzia*dphidzid
     &                             + dedphiprime*dphidzia*dphidzid
     &                             + d2edphi2prime*dphidzia*dphidzid
     &                             + d2drprime*dziazid - dzid*ddrdz
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                             + d2edphi2*dphidxib*dphidxib
     &                             + 2.0d0*dxib*ddrdx + d2drdxx
     &                             + dedphiprime*dphidxib*dphidxib
     &                             + d2edphi2prime*dphidxib*dphidxib
     &                             + d2drprime*dxibxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + dyib*ddrdx + dxib*ddrdy
     &                             + dedphiprime*dphidxib*dphidyib
     &                             + d2edphi2prime*dphidxib*dphidyib
     &                             + d2drprime*dxibyib + d2drdxy
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + dzib*ddrdx + dxib*ddrdz
     &                             + dedphiprime*dphidxib*dphidzib
     &                             + d2edphi2prime*dphidxib*dphidzib
     &                             + d2drprime*dxibzib + d2drdxz
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + dxib*ddrdy + dyib*ddrdx
     &                             + dedphiprime*dphidxib*dphidyib
     &                             + d2edphi2prime*dphidxib*dphidyib
     &                             + d2drprime*dxibyib + d2drdxy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                             + d2edphi2*dphidyib*dphidyib
     &                             + 2.0d0*dyib*ddrdy + d2drdyy
     &                             + dedphiprime*dphidyib*dphidyib
     &                             + d2edphi2prime*dphidyib*dphidyib
     &                             + d2drprime*dyibyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + dzib*ddrdy + dyib*ddrdz
     &                             + dedphiprime*dphidyib*dphidzib
     &                             + d2edphi2prime*dphidyib*dphidzib
     &                             + d2drprime*dyibzib + d2drdyz
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + dxib*ddrdz + dzib*ddrdx
     &                             + dedphiprime*dphidxib*dphidzib
     &                             + d2edphi2prime*dphidxib*dphidzib
     &                             + d2drprime*dxibzib + d2drdxz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + dyib*ddrdz + dzib*ddrdy
     &                             + dedphiprime*dphidyib*dphidzib
     &                             + d2edphi2prime*dphidyib*dphidzib
     &                             + d2drprime*dyibzib + d2drdyz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                             + d2edphi2*dphidzib*dphidzib
     &                             + 2.0d0*dzib*ddrdz + d2drdzz
     &                             + dedphiprime*dphidzib*dphidzib
     &                             + d2edphi2prime*dphidzib*dphidzib
     &                             + d2drprime*dzibzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxib*dphidxia
     &                             + dxia*ddrdx - dxib*ddrdx
     &                             + dedphiprime*dphidxib*dphidxia
     &                             + d2edphi2prime*dphidxib*dphidxia
     &                             + d2drprime*dxiaxib - d2drdxx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                             + d2edphi2*dphidyib*dphidxia
     &                             + dxia*ddrdy - dyib*ddrdx
     &                             + dedphiprime*dphidyib*dphidxia
     &                             + d2edphi2prime*dphidyib*dphidxia
     &                             + d2drprime*dxiayib - d2drdxy
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                             + d2edphi2*dphidzib*dphidxia
     &                             + dxia*ddrdz - dzib*ddrdx
     &                             + dedphiprime*dphidzib*dphidxia
     &                             + d2edphi2prime*dphidzib*dphidxia
     &                             + d2drprime*dxiazib - d2drdxz
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                             + d2edphi2*dphidxib*dphidyia
     &                             + dyia*ddrdx - dxib*ddrdy
     &                             + dedphiprime*dphidxib*dphidyia
     &                             + d2edphi2prime*dphidxib*dphidyia
     &                             + d2drprime*dyiaxib - d2drdxy
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                             + d2edphi2*dphidyib*dphidyia
     &                             + dyia*ddrdy - dyib*ddrdy
     &                             + dedphiprime*dphidyib*dphidyia
     &                             + d2edphi2prime*dphidyib*dphidyia
     &                             + d2drprime*dyiayib - d2drdyy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                             + d2edphi2*dphidzib*dphidyia
     &                             + dyia*ddrdz - dzib*ddrdy
     &                             + dedphiprime*dphidzib*dphidyia
     &                             + d2edphi2prime*dphidzib*dphidyia
     &                             + d2drprime*dyiazib - d2drdyz
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                             + d2edphi2*dphidxib*dphidzia
     &                             + dzia*ddrdx - dxib*ddrdz
     &                             + dedphiprime*dphidxib*dphidzia
     &                             + d2edphi2prime*dphidxib*dphidzia
     &                             + d2drprime*dziaxib - d2drdxz
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                             + d2edphi2*dphidyib*dphidzia
     &                             + dzia*ddrdy - dyib*ddrdz
     &                             + dedphiprime*dphidyib*dphidzia
     &                             + d2edphi2prime*dphidyib*dphidzia
     &                             + d2drprime*dziayib - d2drdyz
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                             + d2edphi2*dphidzib*dphidzia
     &                             + dzia*ddrdz - dzib*ddrdz
     &                             + dedphiprime*dphidzib*dphidzia
     &                             + d2edphi2prime*dphidzib*dphidzia
     &                             + d2drprime*dziazib - d2drdzz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                             + d2edphi2*dphidxib*dphidxic
     &                             + dedphiprime*dphidxib*dphidxic
     &                             + d2edphi2prime*dphidxib*dphidxic
     &                             + d2drprime*dxibxic + dxic*ddrdx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                             + d2edphi2*dphidyib*dphidxic
     &                             + dedphiprime*dphidyib*dphidxic
     &                             + d2edphi2prime*dphidyib*dphidxic
     &                             + d2drprime*dyibxic + dxic*ddrdy
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                             + d2edphi2*dphidzib*dphidxic
     &                             + dedphiprime*dphidzib*dphidxic
     &                             + d2edphi2prime*dphidzib*dphidxic
     &                             + d2drprime*dzibxic + dxic*ddrdz
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                             + d2edphi2*dphidxib*dphidyic
     &                             + dedphiprime*dphidxib*dphidyic
     &                             + d2edphi2prime*dphidxib*dphidyic
     &                             + d2drprime*dxibyic + dyic*ddrdx
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                             + d2edphi2*dphidyib*dphidyic
     &                             + dedphiprime*dphidyib*dphidyic
     &                             + d2edphi2prime*dphidyib*dphidyic
     &                             + d2drprime*dyibyic + dyic*ddrdy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                             + d2edphi2*dphidzib*dphidyic
     &                             + dedphiprime*dphidzib*dphidyic
     &                             + d2edphi2prime*dphidzib*dphidyic
     &                             + d2drprime*dzibyic + dyic*ddrdz
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                             + d2edphi2*dphidxib*dphidzic
     &                             + dedphiprime*dphidxib*dphidzic
     &                             + d2edphi2prime*dphidxib*dphidzic
     &                             + d2drprime*dxibzic + dzic*ddrdx
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                             + d2edphi2*dphidyib*dphidzic
     &                             + dedphiprime*dphidyib*dphidzic
     &                             + d2edphi2prime*dphidyib*dphidzic
     &                             + d2drprime*dyibzic + dzic*ddrdy
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                             + d2edphi2*dphidzib*dphidzic
     &                             + dedphiprime*dphidzib*dphidzic
     &                             + d2edphi2prime*dphidzib*dphidzic
     &                             + d2drprime*dzibzic + dzic*ddrdz
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                             + d2edphi2*dphidxib*dphidxid
     &                             + dedphiprime*dphidxib*dphidxid
     &                             + d2edphi2prime*dphidxib*dphidxid
     &                             + d2drprime*dxibxid + dxid*ddrdx
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                             + d2edphi2*dphidyib*dphidxid
     &                             + dedphiprime*dphidyib*dphidxid
     &                             + d2edphi2prime*dphidyib*dphidxid
     &                             + d2drprime*dyibxid + dxid*ddrdy
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                             + d2edphi2*dphidzib*dphidxid
     &                             + dedphiprime*dphidzib*dphidxid
     &                             + d2edphi2prime*dphidzib*dphidxid
     &                             + d2drprime*dzibxid + dxid*ddrdz
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                             + d2edphi2*dphidxib*dphidyid
     &                             + dedphiprime*dphidxib*dphidyid
     &                             + d2edphi2prime*dphidxib*dphidyid
     &                             + d2drprime*dxibyid + dyid*ddrdx
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                             + d2edphi2*dphidyib*dphidyid
     &                             + dedphiprime*dphidyib*dphidyid
     &                             + d2edphi2prime*dphidyib*dphidyid
     &                             + d2drprime*dyibyid + dyid*ddrdy
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                             + d2edphi2*dphidzib*dphidyid
     &                             + dedphiprime*dphidzib*dphidyid
     &                             + d2edphi2prime*dphidzib*dphidyid
     &                             + d2drprime*dzibyid + dyid*ddrdz
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                             + d2edphi2*dphidxib*dphidzid
     &                             + dedphiprime*dphidxib*dphidzid
     &                             + d2edphi2prime*dphidxib*dphidzid
     &                             + d2drprime*dxibzid + dzid*ddrdx
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                             + d2edphi2*dphidyib*dphidzid
     &                             + dedphiprime*dphidyib*dphidzid
     &                             + d2edphi2prime*dphidyib*dphidzid
     &                             + d2drprime*dyibzid + dzid*ddrdy
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                             + d2edphi2*dphidzib*dphidzid
     &                             + dedphiprime*dphidzib*dphidzid
     &                             + d2edphi2prime*dphidzib*dphidzid
     &                             + d2drprime*dzibzid + dzid*ddrdz
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                             + d2edphi2*dphidxic*dphidxic
     &                             + dedphiprime*dphidxic*dphidxic
     &                             + d2edphi2prime*dphidxic*dphidxic
     &                             + d2drprime*dxicxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + dedphiprime*dphidxic*dphidyic
     &                             + d2edphi2prime*dphidxic*dphidyic
     &                             + d2drprime*dxicyic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + dedphiprime*dphidxic*dphidzic
     &                             + d2edphi2prime*dphidxic*dphidzic
     &                             + d2drprime*dxiczic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + dedphiprime*dphidxic*dphidyic
     &                             + d2edphi2prime*dphidxic*dphidyic
     &                             + d2drprime*dxicyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                             + d2edphi2*dphidyic*dphidyic
     &                             + dedphiprime*dphidyic*dphidyic
     &                             + d2edphi2prime*dphidyic*dphidyic
     &                             + d2drprime*dyicyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + dedphiprime*dphidyic*dphidzic
     &                             + d2edphi2prime*dphidyic*dphidzic
     &                             + d2drprime*dyiczic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + dedphiprime*dphidxic*dphidzic
     &                             + d2edphi2prime*dphidxic*dphidzic
     &                             + d2drprime*dxiczic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + dedphiprime*dphidyic*dphidzic
     &                             + d2edphi2prime*dphidyic*dphidzic
     &                             + d2drprime*dyiczic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                             + d2edphi2*dphidzic*dphidzic
     &                             + dedphiprime*dphidzic*dphidzic
     &                             + d2edphi2prime*dphidzic*dphidzic
     &                             + d2drprime*dziczic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxic*dphidxia
     &                             + dedphiprime*dphidxic*dphidxia
     &                             + d2edphi2prime*dphidxic*dphidxia
     &                             + d2drprime*dxiaxic - dxic*ddrdx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                             + d2edphi2*dphidyic*dphidxia
     &                             + dedphiprime*dphidyic*dphidxia
     &                             + d2edphi2prime*dphidyic*dphidxia
     &                             + d2drprime*dxiayic - dyic*ddrdx
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                             + d2edphi2*dphidzic*dphidxia
     &                             + dedphiprime*dphidzic*dphidxia
     &                             + d2edphi2prime*dphidzic*dphidxia
     &                             + d2drprime*dxiazic - dzic*ddrdx
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                             + d2edphi2*dphidxic*dphidyia
     &                             + dedphiprime*dphidxic*dphidyia
     &                             + d2edphi2prime*dphidxic*dphidyia
     &                             + d2drprime*dyiaxic - dxic*ddrdy
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                             + d2edphi2*dphidyic*dphidyia
     &                             + dedphiprime*dphidyic*dphidyia
     &                             + d2edphi2prime*dphidyic*dphidyia
     &                             + d2drprime*dyiayic - dyic*ddrdy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                             + d2edphi2*dphidzic*dphidyia
     &                             + dedphiprime*dphidzic*dphidyia
     &                             + d2edphi2prime*dphidzic*dphidyia
     &                             + d2drprime*dyiazic - dzic*ddrdy
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                             + d2edphi2*dphidxic*dphidzia
     &                             + dedphiprime*dphidxic*dphidzia
     &                             + d2edphi2prime*dphidxic*dphidzia
     &                             + d2drprime*dziaxic - dxic*ddrdz
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                             + d2edphi2*dphidyic*dphidzia
     &                             + dedphiprime*dphidyic*dphidzia
     &                             + d2edphi2prime*dphidyic*dphidzia
     &                             + d2drprime*dziayic - dyic*ddrdz
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                             + d2edphi2*dphidzic*dphidzia
     &                             + dedphiprime*dphidzic*dphidzia
     &                             + d2edphi2prime*dphidzic*dphidzia
     &                             + d2drprime*dziazic - dzic*ddrdz
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                             + d2edphi2*dphidxic*dphidxib
     &                             + dedphiprime*dphidxic*dphidxib
     &                             + d2edphi2prime*dphidxic*dphidxib
     &                             + d2drprime*dxibxic + dxic*ddrdx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                             + d2edphi2*dphidyic*dphidxib
     &                             + dedphiprime*dphidyic*dphidxib
     &                             + d2edphi2prime*dphidyic*dphidxib
     &                             + d2drprime*dxibyic + dyic*ddrdx
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                             + d2edphi2*dphidzic*dphidxib
     &                             + dedphiprime*dphidzic*dphidxib
     &                             + d2edphi2prime*dphidzic*dphidxib
     &                             + d2drprime*dxibzic + dzic*ddrdx
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                             + d2edphi2*dphidxic*dphidyib
     &                             + dedphiprime*dphidxic*dphidyib
     &                             + d2edphi2prime*dphidxic*dphidyib
     &                             + d2drprime*dyibxic + dxic*ddrdy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                             + d2edphi2*dphidyic*dphidyib
     &                             + dedphiprime*dphidyic*dphidyib
     &                             + d2edphi2prime*dphidyic*dphidyib
     &                             + d2drprime*dyibyic + dyic*ddrdy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                             + d2edphi2*dphidzic*dphidyib
     &                             + dedphiprime*dphidzic*dphidyib
     &                             + d2edphi2prime*dphidzic*dphidyib
     &                             + d2drprime*dyibzic + dzic*ddrdy
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                             + d2edphi2*dphidxic*dphidzib
     &                             + dedphiprime*dphidxic*dphidzib
     &                             + d2edphi2prime*dphidxic*dphidzib
     &                             + d2drprime*dzibxic + dxic*ddrdz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                             + d2edphi2*dphidyic*dphidzib
     &                             + dedphiprime*dphidyic*dphidzib
     &                             + d2edphi2prime*dphidyic*dphidzib
     &                             + d2drprime*dzibyic + dyic*ddrdz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                             + d2edphi2*dphidzic*dphidzib
     &                             + dedphiprime*dphidzic*dphidzib
     &                             + d2edphi2prime*dphidzic*dphidzib
     &                             + d2drprime*dzibzic + dzic*ddrdz
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                             + d2edphi2*dphidxic*dphidxid
     &                             + dedphiprime*dphidxic*dphidxid
     &                             + d2edphi2prime*dphidxic*dphidxid
     &                             + d2drprime*dxicxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                             + d2edphi2*dphidyic*dphidxid
     &                             + dedphiprime*dphidyic*dphidxid
     &                             + d2edphi2prime*dphidyic*dphidxid
     &                             + d2drprime*dyicxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                             + d2edphi2*dphidzic*dphidxid
     &                             + dedphiprime*dphidzic*dphidxid
     &                             + d2edphi2prime*dphidzic*dphidxid
     &                             + d2drprime*dzicxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                             + d2edphi2*dphidxic*dphidyid
     &                             + dedphiprime*dphidxic*dphidyid
     &                             + d2edphi2prime*dphidxic*dphidyid
     &                             + d2drprime*dxicyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                             + d2edphi2*dphidyic*dphidyid
     &                             + dedphiprime*dphidyic*dphidyid
     &                             + d2edphi2prime*dphidyic*dphidyid
     &                             + d2drprime*dyicyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                             + d2edphi2*dphidzic*dphidyid
     &                             + dedphiprime*dphidzic*dphidyid
     &                             + d2edphi2prime*dphidzic*dphidyid
     &                             + d2drprime*dzicyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                             + d2edphi2*dphidxic*dphidzid
     &                             + dedphiprime*dphidxic*dphidzid
     &                             + d2edphi2prime*dphidxic*dphidzid
     &                             + d2drprime*dxiczid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                             + d2edphi2*dphidyic*dphidzid
     &                             + dedphiprime*dphidyic*dphidzid
     &                             + d2edphi2prime*dphidyic*dphidzid
     &                             + d2drprime*dyiczid
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                             + d2edphi2*dphidzic*dphidzid
     &                             + dedphiprime*dphidzic*dphidzid
     &                             + d2edphi2prime*dphidzic*dphidzid
     &                             + d2drprime*dziczid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                             + d2edphi2*dphidxid*dphidxid
     &                             + dedphiprime*dphidxid*dphidxid
     &                             + d2edphi2prime*dphidxid*dphidxid
     &                             + d2drprime*dxidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + dedphiprime*dphidxid*dphidyid
     &                             + d2edphi2prime*dphidxid*dphidyid
     &                             + d2drprime*dxidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + dedphiprime*dphidxid*dphidzid
     &                             + d2edphi2prime*dphidxid*dphidzid
     &                             + d2drprime*dxidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + dedphiprime*dphidxid*dphidyid
     &                             + d2edphi2prime*dphidxid*dphidyid
     &                             + d2drprime*dxidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                             + d2edphi2*dphidyid*dphidyid
     &                             + dedphiprime*dphidyid*dphidyid
     &                             + d2edphi2prime*dphidyid*dphidyid
     &                             + d2drprime*dyidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + dedphiprime*dphidyid*dphidzid
     &                             + d2edphi2prime*dphidyid*dphidzid
     &                             + d2drprime*dyidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + dedphiprime*dphidxid*dphidzid
     &                             + d2edphi2prime*dphidxid*dphidzid
     &                             + d2drprime*dxidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + dedphiprime*dphidyid*dphidzid
     &                             + d2edphi2prime*dphidyid*dphidzid
     &                             + d2drprime*dyidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                             + d2edphi2*dphidzid*dphidzid
     &                             + dedphiprime*dphidzid*dphidzid
     &                             + d2edphi2prime*dphidzid*dphidzid
     &                             + d2drprime*dzidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxid*dphidxia
     &                             + dedphiprime*dphidxid*dphidxia
     &                             + d2edphi2prime*dphidxid*dphidxia
     &                             + d2drprime*dxiaxid - dxid*ddrdx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                             + d2edphi2*dphidyid*dphidxia
     &                             + dedphiprime*dphidyid*dphidxia
     &                             + d2edphi2prime*dphidyid*dphidxia
     &                             + d2drprime*dxiayid - dyid*ddrdx
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                             + d2edphi2*dphidzid*dphidxia
     &                             + dedphiprime*dphidzid*dphidxia
     &                             + d2edphi2prime*dphidzid*dphidxia
     &                             + d2drprime*dxiazid - dzid*ddrdx
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                             + d2edphi2*dphidxid*dphidyia
     &                             + dedphiprime*dphidxid*dphidyia
     &                             + d2edphi2prime*dphidxid*dphidyia
     &                             + d2drprime*dyiaxid - dxid*ddrdy
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                             + d2edphi2*dphidyid*dphidyia
     &                             + dedphiprime*dphidyid*dphidyia
     &                             + d2edphi2prime*dphidyid*dphidyia
     &                             + d2drprime*dyiayid - dyid*ddrdy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                             + d2edphi2*dphidzid*dphidyia
     &                             + dedphiprime*dphidzid*dphidyia
     &                             + d2edphi2prime*dphidzid*dphidyia
     &                             + d2drprime*dyiazid - dzid*ddrdy
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                             + d2edphi2*dphidxid*dphidzia
     &                             + dedphiprime*dphidxid*dphidzia
     &                             + d2edphi2prime*dphidxid*dphidzia
     &                             + d2drprime*dziaxid - dxid*ddrdz
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                             + d2edphi2*dphidyid*dphidzia
     &                             + dedphiprime*dphidyid*dphidzia
     &                             + d2edphi2prime*dphidyid*dphidzia
     &                             + d2drprime*dziayid - dyid*ddrdz
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                             + d2edphi2*dphidzid*dphidzia
     &                             + dedphiprime*dphidzid*dphidzia
     &                             + d2edphi2prime*dphidzid*dphidzia
     &                             + d2drprime*dziazid - dzid*ddrdz
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                             + d2edphi2*dphidxid*dphidxib
     &                             + dedphiprime*dphidxid*dphidxib
     &                             + d2edphi2prime*dphidxid*dphidxib
     &                             + d2drprime*dxibxid + dxid*ddrdx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                             + d2edphi2*dphidyid*dphidxib
     &                             + dedphiprime*dphidyid*dphidxib
     &                             + d2edphi2prime*dphidyid*dphidxib
     &                             + d2drprime*dxibyid + dyid*ddrdx
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                             + d2edphi2*dphidzid*dphidxib
     &                             + dedphiprime*dphidzid*dphidxib
     &                             + d2edphi2prime*dphidzid*dphidxib
     &                             + d2drprime*dxibzid + dzid*ddrdx
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                             + d2edphi2*dphidxid*dphidyib
     &                             + dedphiprime*dphidxid*dphidyib
     &                             + d2edphi2prime*dphidxid*dphidyib
     &                             + d2drprime*dyibxid + dxid*ddrdy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                             + d2edphi2*dphidyid*dphidyib
     &                             + dedphiprime*dphidyid*dphidyib
     &                             + d2edphi2prime*dphidyid*dphidyib
     &                             + d2drprime*dyibyid + dyid*ddrdy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                             + d2edphi2*dphidzid*dphidyib
     &                             + dedphiprime*dphidzid*dphidyib
     &                             + d2edphi2prime*dphidzid*dphidyib
     &                             + d2drprime*dyibzid + dzid*ddrdy
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                             + d2edphi2*dphidxid*dphidzib
     &                             + dedphiprime*dphidxid*dphidzib
     &                             + d2edphi2prime*dphidxid*dphidzib
     &                             + d2drprime*dzibxid + dxid*ddrdz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                             + d2edphi2*dphidyid*dphidzib
     &                             + dedphiprime*dphidyid*dphidzib
     &                             + d2edphi2prime*dphidyid*dphidzib
     &                             + d2drprime*dzibyid + dyid*ddrdz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                             + d2edphi2*dphidzid*dphidzib
     &                             + dedphiprime*dphidzid*dphidzib
     &                             + d2edphi2prime*dphidzid*dphidzib
     &                             + d2drprime*dzibzid + dzid*ddrdz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                             + d2edphi2*dphidxid*dphidxic
     &                             + dedphiprime*dphidxid*dphidxic
     &                             + d2edphi2prime*dphidxid*dphidxic
     &                             + d2drprime*dxicxid
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                             + d2edphi2*dphidyid*dphidxic
     &                             + dedphiprime*dphidyid*dphidxic
     &                             + d2edphi2prime*dphidyid*dphidxic
     &                             + d2drprime*dxicyid
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                             + d2edphi2*dphidzid*dphidxic
     &                             + dedphiprime*dphidzid*dphidxic
     &                             + d2edphi2prime*dphidzid*dphidxic
     &                             + d2drprime*dxiczid
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                             + d2edphi2*dphidxid*dphidyic
     &                             + dedphiprime*dphidxid*dphidyic
     &                             + d2edphi2prime*dphidxid*dphidyic
     &                             + d2drprime*dyicxid
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                             + d2edphi2*dphidyid*dphidyic
     &                             + dedphiprime*dphidyid*dphidyic
     &                             + d2edphi2prime*dphidyid*dphidyic
     &                             + d2drprime*dyicyid
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                             + d2edphi2*dphidzid*dphidyic
     &                             + dedphiprime*dphidzid*dphidyic
     &                             + d2edphi2prime*dphidzid*dphidyic
     &                             + d2drprime*dyiczid
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                             + d2edphi2*dphidxid*dphidzic
     &                             + dedphiprime*dphidxid*dphidzic
     &                             + d2edphi2prime*dphidxid*dphidzic
     &                             + d2drprime*dzicxid
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                             + d2edphi2*dphidyid*dphidzic
     &                             + dedphiprime*dphidyid*dphidzic
     &                             + d2edphi2prime*dphidyid*dphidzic
     &                             + d2drprime*dzicyid
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                             + d2edphi2*dphidzid*dphidzic
     &                             + dedphiprime*dphidzid*dphidzic
     &                             + d2edphi2prime*dphidzid*dphidzic
     &                             + d2drprime*dziczid
               end if
c
c     get the stretch-torsion values for the second bond
c
               v1 = kst(4,istrtor)
               v2 = kst(5,istrtor)
               v3 = kst(6,istrtor)
               k = ist(3,istrtor)
               dr = rcb - bl(k)
               force = bk(k)
               dedphi = storunit * 2.0d0 * force 
     &                       * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2 = storunit * 2.0d0 * force * dr
     &                       * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3)
               ddr = 1.0d0 / rcb
               d2dr = -storunit * 2.0d0 * force 
     &                    * (v1*phi1 + v2*phi2 + v3*phi3) / rcb**3
               dedphiprime = storunit * 2.0d0 * force
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2prime = storunit * 2.0d0 * force
     &                            * (v1*phi1 + v2*phi2 + v3*phi3)
     &                            * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3)
               d2drprime = storunit * 2.0d0 * force
     &                        * (v1*phi1 + v2*phi2 + v3*phi3)
     &                        * (v1*dphi1 + v2*dphi2 + v3*dphi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
                  d2dr = d2dr * fgrp
                  dedphiprime = dedphiprime * fgrp
                  d2edphi2prime = d2edphi2prime * fgrp
                  d2drprime = d2drprime * fgrp
               end if
c
c     compute derivative components for this interaction
c
               ddrdx = xcb * ddr
               ddrdy = ycb * ddr
               ddrdz = zcb * ddr
               d2drdxx = (xcb*xcb-rcb*rcb) * d2dr
               d2drdyy = (ycb*ycb-rcb*rcb) * d2dr
               d2drdzz = (zcb*zcb-rcb*rcb) * d2dr
               d2drdxy = xcb * ycb * d2dr
               d2drdxz = xcb * zcb * d2dr
               d2drdyz = ycb * zcb * d2dr
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
c     increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                             + d2edphi2*dphidxia*dphidxia
     &                             + dedphiprime*dphidxia*dphidxia
     &                             + d2edphi2prime*dphidxia*dphidxia
     &                             + d2drprime*dxiaxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             + dedphiprime*dphidxia*dphidyia
     &                             + d2edphi2prime*dphidxia*dphidyia
     &                             + d2drprime*dxiayia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             + dedphiprime*dphidxia*dphidzia
     &                             + d2edphi2prime*dphidxia*dphidzia
     &                             + d2drprime*dxiazia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             + dedphiprime*dphidxia*dphidyia
     &                             + d2edphi2prime*dphidxia*dphidyia
     &                             + d2drprime*dxiayia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                             + d2edphi2*dphidyia*dphidyia
     &                             + dedphiprime*dphidyia*dphidyia
     &                             + d2edphi2prime*dphidyia*dphidyia
     &                             + d2drprime*dyiayia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             + dedphiprime*dphidyia*dphidzia
     &                             + d2edphi2prime*dphidyia*dphidzia
     &                             + d2drprime*dyiazia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             + dedphiprime*dphidxia*dphidzia
     &                             + d2edphi2prime*dphidxia*dphidzia
     &                             + d2drprime*dxiazia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             + dedphiprime*dphidyia*dphidzia
     &                             + d2edphi2prime*dphidyia*dphidzia
     &                             + d2drprime*dyiazia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                             + d2edphi2*dphidzia*dphidzia
     &                             + dedphiprime*dphidzia*dphidzia
     &                             + d2edphi2prime*dphidzia*dphidzia
     &                             + d2drprime*dziazia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxia*dphidxib
     &                             + dedphiprime*dphidxia*dphidxib
     &                             + d2edphi2prime*dphidxia*dphidxib
     &                             + d2drprime*dxiaxib - dxia*ddrdx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                             + d2edphi2*dphidyia*dphidxib
     &                             + dedphiprime*dphidyia*dphidxib
     &                             + d2edphi2prime*dphidyia*dphidxib
     &                             + d2drprime*dyiaxib - dyia*ddrdx
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                             + d2edphi2*dphidzia*dphidxib
     &                             + dedphiprime*dphidzia*dphidxib
     &                             + d2edphi2prime*dphidzia*dphidxib
     &                             + d2drprime*dziaxib - dzia*ddrdx
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                             + d2edphi2*dphidxia*dphidyib
     &                             + dedphiprime*dphidxia*dphidyib
     &                             + d2edphi2prime*dphidxia*dphidyib
     &                             + d2drprime*dxiayib - dxia*ddrdy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                             + d2edphi2*dphidyia*dphidyib
     &                             + dedphiprime*dphidyia*dphidyib
     &                             + d2edphi2prime*dphidyia*dphidyib
     &                             + d2drprime*dyiayib - dyia*ddrdy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                             + d2edphi2*dphidzia*dphidyib
     &                             + dedphiprime*dphidzia*dphidyib
     &                             + d2edphi2prime*dphidzia*dphidyib
     &                             + d2drprime*dziayib - dzia*ddrdy
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                             + d2edphi2*dphidxia*dphidzib
     &                             + dedphiprime*dphidxia*dphidzib
     &                             + d2edphi2prime*dphidxia*dphidzib
     &                             + d2drprime*dxiazib - dxia*ddrdz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                             + d2edphi2*dphidyia*dphidzib
     &                             + dedphiprime*dphidyia*dphidzib
     &                             + d2edphi2prime*dphidyia*dphidzib
     &                             + d2drprime*dyiazib - dyia*ddrdz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                             + d2edphi2*dphidzia*dphidzib
     &                             + dedphiprime*dphidzia*dphidzib
     &                             + d2edphi2prime*dphidzia*dphidzib
     &                             + d2drprime*dziazib - dzia*ddrdz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxia*dphidxic
     &                             + dedphiprime*dphidxia*dphidxic
     &                             + d2edphi2prime*dphidxia*dphidxic
     &                             + d2drprime*dxiaxic + dxia*ddrdx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                             + d2edphi2*dphidyia*dphidxic
     &                             + dedphiprime*dphidyia*dphidxic
     &                             + d2edphi2prime*dphidyia*dphidxic
     &                             + d2drprime*dyiaxic + dyia*ddrdx
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                             + d2edphi2*dphidzia*dphidxic
     &                             + dedphiprime*dphidzia*dphidxic
     &                             + d2edphi2prime*dphidzia*dphidxic
     &                             + d2drprime*dziaxic + dzia*ddrdx
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                             + d2edphi2*dphidxia*dphidyic
     &                             + dedphiprime*dphidxia*dphidyic
     &                             + d2edphi2prime*dphidxia*dphidyic
     &                             + d2drprime*dxiayic + dxia*ddrdy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                             + d2edphi2*dphidyia*dphidyic
     &                             + dedphiprime*dphidyia*dphidyic
     &                             + d2edphi2prime*dphidyia*dphidyic
     &                             + d2drprime*dyiayic + dyia*ddrdy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                             + d2edphi2*dphidzia*dphidyic
     &                             + dedphiprime*dphidzia*dphidyic
     &                             + d2edphi2prime*dphidzia*dphidyic
     &                             + d2drprime*dziayic + dzia*ddrdy
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                             + d2edphi2*dphidxia*dphidzic
     &                             + dedphiprime*dphidxia*dphidzic
     &                             + d2edphi2prime*dphidxia*dphidzic
     &                             + d2drprime*dxiazic + dxia*ddrdz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                             + d2edphi2*dphidyia*dphidzic
     &                             + dedphiprime*dphidyia*dphidzic
     &                             + d2edphi2prime*dphidyia*dphidzic
     &                             + d2drprime*dyiazic + dyia*ddrdz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                             + d2edphi2*dphidzia*dphidzic
     &                             + dedphiprime*dphidzia*dphidzic
     &                             + d2edphi2prime*dphidzia*dphidzic
     &                             + d2drprime*dziazic + dzia*ddrdz
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxia*dphidxid
     &                             + dedphiprime*dphidxia*dphidxid
     &                             + d2edphi2prime*dphidxia*dphidxid
     &                             + d2drprime*dxiaxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                             + d2edphi2*dphidyia*dphidxid
     &                             + dedphiprime*dphidyia*dphidxid
     &                             + d2edphi2prime*dphidyia*dphidxid
     &                             + d2drprime*dyiaxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                             + d2edphi2*dphidzia*dphidxid
     &                             + dedphiprime*dphidzia*dphidxid
     &                             + d2edphi2prime*dphidzia*dphidxid
     &                             + d2drprime*dziaxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                             + d2edphi2*dphidxia*dphidyid
     &                             + dedphiprime*dphidxia*dphidyid
     &                             + d2edphi2prime*dphidxia*dphidyid
     &                             + d2drprime*dxiayid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                             + d2edphi2*dphidyia*dphidyid
     &                             + dedphiprime*dphidyia*dphidyid
     &                             + d2edphi2prime*dphidyia*dphidyid
     &                             + d2drprime*dyiayid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                             + d2edphi2*dphidzia*dphidyid
     &                             + dedphiprime*dphidzia*dphidyid
     &                             + d2edphi2prime*dphidzia*dphidyid
     &                             + d2drprime*dziayid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                             + d2edphi2*dphidxia*dphidzid
     &                             + dedphiprime*dphidxia*dphidzid
     &                             + d2edphi2prime*dphidxia*dphidzid
     &                             + d2drprime*dxiazid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                             + d2edphi2*dphidyia*dphidzid
     &                             + dedphiprime*dphidyia*dphidzid
     &                             + d2edphi2prime*dphidyia*dphidzid
     &                             + d2drprime*dyiazid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                             + d2edphi2*dphidzia*dphidzid
     &                             + dedphiprime*dphidzia*dphidzid
     &                             + d2edphi2prime*dphidzia*dphidzid
     &                             + d2drprime*dziazid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                             + d2edphi2*dphidxib*dphidxib
     &                             - 2.0d0*dxib*ddrdx + d2drdxx
     &                             + dedphiprime*dphidxib*dphidxib
     &                             + d2edphi2prime*dphidxib*dphidxib
     &                             + d2drprime*dxibxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             - dyib*ddrdx - dxib*ddrdy
     &                             + dedphiprime*dphidxib*dphidyib
     &                             + d2edphi2prime*dphidxib*dphidyib
     &                             + d2drprime*dxibyib + d2drdxy
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             - dzib*ddrdx - dxib*ddrdz
     &                             + dedphiprime*dphidxib*dphidzib
     &                             + d2edphi2prime*dphidxib*dphidzib
     &                             + d2drprime*dxibzib + d2drdxz
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             - dxib*ddrdy - dyib*ddrdx
     &                             + dedphiprime*dphidxib*dphidyib
     &                             + d2edphi2prime*dphidxib*dphidyib
     &                             + d2drprime*dxibyib + d2drdxy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                             + d2edphi2*dphidyib*dphidyib
     &                             - 2.0d0*dyib*ddrdy + d2drdyy
     &                             + dedphiprime*dphidyib*dphidyib
     &                             + d2edphi2prime*dphidyib*dphidyib
     &                             + d2drprime*dyibyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             - dzib*ddrdy - dyib*ddrdz
     &                             + dedphiprime*dphidyib*dphidzib
     &                             + d2edphi2prime*dphidyib*dphidzib
     &                             + d2drprime*dyibzib + d2drdyz
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             - dxib*ddrdz - dzib*ddrdx
     &                             + dedphiprime*dphidxib*dphidzib
     &                             + d2edphi2prime*dphidxib*dphidzib
     &                             + d2drprime*dxibzib + d2drdxz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             - dyib*ddrdz - dzib*ddrdy
     &                             + dedphiprime*dphidyib*dphidzib
     &                             + d2edphi2prime*dphidyib*dphidzib
     &                             + d2drprime*dyibzib + d2drdyz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                             + d2edphi2*dphidzib*dphidzib
     &                             - 2.0d0*dzib*ddrdz + d2drdzz
     &                             + dedphiprime*dphidzib*dphidzib
     &                             + d2edphi2prime*dphidzib*dphidzib
     &                             + d2drprime*dzibzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxib*dphidxia
     &                             + dedphiprime*dphidxib*dphidxia
     &                             + d2edphi2prime*dphidxib*dphidxia
     &                             + d2drprime*dxiaxib - dxia*ddrdx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                             + d2edphi2*dphidyib*dphidxia
     &                             + dedphiprime*dphidyib*dphidxia
     &                             + d2edphi2prime*dphidyib*dphidxia
     &                             + d2drprime*dxiayib - dxia*ddrdy
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                             + d2edphi2*dphidzib*dphidxia
     &                             + dedphiprime*dphidzib*dphidxia
     &                             + d2edphi2prime*dphidzib*dphidxia
     &                             + d2drprime*dxiazib - dxia*ddrdz
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                             + d2edphi2*dphidxib*dphidyia
     &                             + dedphiprime*dphidxib*dphidyia
     &                             + d2edphi2prime*dphidxib*dphidyia
     &                             + d2drprime*dyiaxib - dyia*ddrdx
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                             + d2edphi2*dphidyib*dphidyia
     &                             + dedphiprime*dphidyib*dphidyia
     &                             + d2edphi2prime*dphidyib*dphidyia
     &                             + d2drprime*dyiayib - dyia*ddrdy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                             + d2edphi2*dphidzib*dphidyia
     &                             + dedphiprime*dphidzib*dphidyia
     &                             + d2edphi2prime*dphidzib*dphidyia
     &                             + d2drprime*dyiazib - dyia*ddrdz
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                             + d2edphi2*dphidxib*dphidzia
     &                             - dzia*ddrdx
     &                             + dedphiprime*dphidxib*dphidzia
     &                             + d2edphi2prime*dphidxib*dphidzia
     &                             + d2drprime*dziaxib
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                             + d2edphi2*dphidyib*dphidzia
     &                             - dzia*ddrdy
     &                             + dedphiprime*dphidyib*dphidzia
     &                             + d2edphi2prime*dphidyib*dphidzia
     &                             + d2drprime*dziayib
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                             + d2edphi2*dphidzib*dphidzia
     &                             - dzia*ddrdz
     &                             + dedphiprime*dphidzib*dphidzia
     &                             + d2edphi2prime*dphidzib*dphidzia
     &                             + d2drprime*dziazib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                             + d2edphi2*dphidxib*dphidxic
     &                             + (dxib-dxic)*ddrdx - d2drdxx
     &                             + dedphiprime*dphidxib*dphidxic
     &                             + d2edphi2prime*dphidxib*dphidxic
     &                             + d2drprime*dxibxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                             + d2edphi2*dphidyib*dphidxic
     &                             + dyib*ddrdx - dxic*ddrdy
     &                             + dedphiprime*dphidyib*dphidxic
     &                             + d2edphi2prime*dphidyib*dphidxic
     &                             + d2drprime*dyibxic - d2drdxy
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                             + d2edphi2*dphidzib*dphidxic
     &                             + dzib*ddrdx - dxic*ddrdz
     &                             + dedphiprime*dphidzib*dphidxic
     &                             + d2edphi2prime*dphidzib*dphidxic
     &                             + d2drprime*dzibxic - d2drdxz
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                             + d2edphi2*dphidxib*dphidyic
     &                             + dxib*ddrdy - dyic*ddrdx
     &                             + dedphiprime*dphidxib*dphidyic
     &                             + d2edphi2prime*dphidxib*dphidyic
     &                             + d2drprime*dxibyic - d2drdxy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                             + d2edphi2*dphidyib*dphidyic
     &                             + (dyib-dyic)*ddrdy - d2drdyy
     &                             + dedphiprime*dphidyib*dphidyic
     &                             + d2edphi2prime*dphidyib*dphidyic
     &                             + d2drprime*dyibyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                             + d2edphi2*dphidzib*dphidyic
     &                             + dzib*ddrdy - dyic*ddrdz
     &                             + dedphiprime*dphidzib*dphidyic
     &                             + d2edphi2prime*dphidzib*dphidyic
     &                             + d2drprime*dzibyic - d2drdyz
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                             + d2edphi2*dphidxib*dphidzic
     &                             + dxib*ddrdz - dzic*ddrdx
     &                             + dedphiprime*dphidxib*dphidzic
     &                             + d2edphi2prime*dphidxib*dphidzic
     &                             + d2drprime*dxibzic - d2drdxz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                             + d2edphi2*dphidyib*dphidzic
     &                             + dyib*ddrdz - dzic*ddrdy
     &                             + dedphiprime*dphidyib*dphidzic
     &                             + d2edphi2prime*dphidyib*dphidzic
     &                             + d2drprime*dyibzic - d2drdyz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                             + d2edphi2*dphidzib*dphidzic
     &                             + (dzib-dzic)*ddrdz - d2drdzz
     &                             + dedphiprime*dphidzib*dphidzic
     &                             + d2edphi2prime*dphidzib*dphidzic
     &                             + d2drprime*dzibzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                             + d2edphi2*dphidxib*dphidxid
     &                             + dedphiprime*dphidxib*dphidxid
     &                             + d2edphi2prime*dphidxib*dphidxid
     &                             + d2drprime*dxibxid - dxid*ddrdx
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                             + d2edphi2*dphidyib*dphidxid
     &                             + dedphiprime*dphidyib*dphidxid
     &                             + d2edphi2prime*dphidyib*dphidxid
     &                             + d2drprime*dyibxid - dxid*ddrdy
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                             + d2edphi2*dphidzib*dphidxid
     &                             + dedphiprime*dphidzib*dphidxid
     &                             + d2edphi2prime*dphidzib*dphidxid
     &                             + d2drprime*dzibxid - dxid*ddrdz
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                             + d2edphi2*dphidxib*dphidyid
     &                             + dedphiprime*dphidxib*dphidyid
     &                             + d2edphi2prime*dphidxib*dphidyid
     &                             + d2drprime*dxibyid - dyid*ddrdx
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                             + d2edphi2*dphidyib*dphidyid
     &                             + dedphiprime*dphidyib*dphidyid
     &                             + d2edphi2prime*dphidyib*dphidyid
     &                             + d2drprime*dyibyid - dyid*ddrdy
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                             + d2edphi2*dphidzib*dphidyid
     &                             + dedphiprime*dphidzib*dphidyid
     &                             + d2edphi2prime*dphidzib*dphidyid
     &                             + d2drprime*dzibyid - dyid*ddrdz
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                             + d2edphi2*dphidxib*dphidzid
     &                             + dedphiprime*dphidxib*dphidzid
     &                             + d2edphi2prime*dphidxib*dphidzid
     &                             + d2drprime*dxibzid - dzid*ddrdx
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                             + d2edphi2*dphidyib*dphidzid
     &                             + dedphiprime*dphidyib*dphidzid
     &                             + d2edphi2prime*dphidyib*dphidzid
     &                             + d2drprime*dyibzid - dzid*ddrdy
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                             + d2edphi2*dphidzib*dphidzid
     &                             + dedphiprime*dphidzib*dphidzid
     &                             + d2edphi2prime*dphidzib*dphidzid
     &                             + d2drprime*dzibzid - dzid*ddrdz
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                             + d2edphi2*dphidxic*dphidxic
     &                             + 2.0d0*dxic*ddrdx + d2drdxx
     &                             + dedphiprime*dphidxic*dphidxic
     &                             + d2edphi2prime*dphidxic*dphidxic
     &                             + d2drprime*dxicxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + dyic*ddrdx + dxic*ddrdy
     &                             + dedphiprime*dphidxic*dphidyic
     &                             + d2edphi2prime*dphidxic*dphidyic
     &                             + d2drprime*dxicyic + d2drdxy
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + dzic*ddrdx + dxic*ddrdz
     &                             + dedphiprime*dphidxic*dphidzic
     &                             + d2edphi2prime*dphidxic*dphidzic
     &                             + d2drprime*dxiczic + d2drdxz
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + dxic*ddrdy + dyic*ddrdx
     &                             + dedphiprime*dphidxic*dphidyic
     &                             + d2edphi2prime*dphidxic*dphidyic
     &                             + d2drprime*dxicyic + d2drdxy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                             + d2edphi2*dphidyic*dphidyic
     &                             + 2.0d0*dyic*ddrdy + d2drdyy
     &                             + dedphiprime*dphidyic*dphidyic
     &                             + d2edphi2prime*dphidyic*dphidyic
     &                             + d2drprime*dyicyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + dzic*ddrdy + dyic*ddrdz
     &                             + dedphiprime*dphidyic*dphidzic
     &                             + d2edphi2prime*dphidyic*dphidzic
     &                             + d2drprime*dyiczic + d2drdyz
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + dxic*ddrdz + dzic*ddrdx
     &                             + dedphiprime*dphidxic*dphidzic
     &                             + d2edphi2prime*dphidxic*dphidzic
     &                             + d2drprime*dxiczic + d2drdxz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + dyic*ddrdz + dzic*ddrdy
     &                             + dedphiprime*dphidyic*dphidzic
     &                             + d2edphi2prime*dphidyic*dphidzic
     &                             + d2drprime*dyiczic + d2drdyz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                             + d2edphi2*dphidzic*dphidzic
     &                             + 2.0d0*dzic*ddrdz + d2drdzz
     &                             + dedphiprime*dphidzic*dphidzic
     &                             + d2edphi2prime*dphidzic*dphidzic
     &                             + d2drprime*dziczic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxic*dphidxia
     &                             + dedphiprime*dphidxic*dphidxia
     &                             + d2edphi2prime*dphidxic*dphidxia
     &                             + d2drprime*dxiaxic + dxia*ddrdx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                             + d2edphi2*dphidyic*dphidxia
     &                             + dedphiprime*dphidyic*dphidxia
     &                             + d2edphi2prime*dphidyic*dphidxia
     &                             + d2drprime*dxiayic + dxia*ddrdy
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                             + d2edphi2*dphidzic*dphidxia
     &                             + dedphiprime*dphidzic*dphidxia
     &                             + d2edphi2prime*dphidzic*dphidxia
     &                             + d2drprime*dxiazic + dxia*ddrdz
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                             + d2edphi2*dphidxic*dphidyia
     &                             + dedphiprime*dphidxic*dphidyia
     &                             + d2edphi2prime*dphidxic*dphidyia
     &                             + d2drprime*dyiaxic + dyia*ddrdx
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                             + d2edphi2*dphidyic*dphidyia
     &                             + dedphiprime*dphidyic*dphidyia
     &                             + d2edphi2prime*dphidyic*dphidyia
     &                             + d2drprime*dyiayic + dyia*ddrdy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                             + d2edphi2*dphidzic*dphidyia
     &                             + dedphiprime*dphidzic*dphidyia
     &                             + d2edphi2prime*dphidzic*dphidyia
     &                             + d2drprime*dyiazic + dyia*ddrdz
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                             + d2edphi2*dphidxic*dphidzia
     &                             + dedphiprime*dphidxic*dphidzia
     &                             + d2edphi2prime*dphidxic*dphidzia
     &                             + d2drprime*dziaxic + dzia*ddrdx
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                             + d2edphi2*dphidyic*dphidzia
     &                             + dedphiprime*dphidyic*dphidzia
     &                             + d2edphi2prime*dphidyic*dphidzia
     &                             + d2drprime*dziayic + dzia*ddrdy
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                             + d2edphi2*dphidzic*dphidzia
     &                             + dedphiprime*dphidzic*dphidzia
     &                             + d2edphi2prime*dphidzic*dphidzia
     &                             + d2drprime*dziazic + dzia*ddrdz
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                             + d2edphi2*dphidxic*dphidxib
     &                             - (dxic-dxib)*ddrdx - d2drdxx
     &                             + dedphiprime*dphidxic*dphidxib
     &                             + d2edphi2prime*dphidxic*dphidxib
     &                             + d2drprime*dxibxic
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                             + d2edphi2*dphidyic*dphidxib
     &                             - dyic*ddrdx + dxib*ddrdy
     &                             + dedphiprime*dphidyic*dphidxib
     &                             + d2edphi2prime*dphidyic*dphidxib
     &                             + d2drprime*dxibyic - d2drdxy
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                             + d2edphi2*dphidzic*dphidxib
     &                             - dzic*ddrdx + dxib*ddrdz
     &                             + dedphiprime*dphidzic*dphidxib
     &                             + d2edphi2prime*dphidzic*dphidxib
     &                             + d2drprime*dxibzic - d2drdxz
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                             + d2edphi2*dphidxic*dphidyib
     &                             - dxic*ddrdy + dyib*ddrdx
     &                             + dedphiprime*dphidxic*dphidyib
     &                             + d2edphi2prime*dphidxic*dphidyib
     &                             + d2drprime*dyibxic - d2drdxy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                             + d2edphi2*dphidyic*dphidyib
     &                             - (dyic-dyib)*ddrdy - d2drdyy
     &                             + dedphiprime*dphidyic*dphidyib
     &                             + d2edphi2prime*dphidyic*dphidyib
     &                             + d2drprime*dyibyic
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                             + d2edphi2*dphidzic*dphidyib
     &                             - dzic*ddrdy + dyib*ddrdz
     &                             + dedphiprime*dphidzic*dphidyib
     &                             + d2edphi2prime*dphidzic*dphidyib
     &                             + d2drprime*dyibzic - d2drdyz
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                             + d2edphi2*dphidxic*dphidzib
     &                             - dxic*ddrdz + dzib*ddrdx
     &                             + dedphiprime*dphidxic*dphidzib
     &                             + d2edphi2prime*dphidxic*dphidzib
     &                             + d2drprime*dzibxic - d2drdxz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                             + d2edphi2*dphidyic*dphidzib
     &                             - dyic*ddrdz + dzib*ddrdy
     &                             + dedphiprime*dphidyic*dphidzib
     &                             + d2edphi2prime*dphidyic*dphidzib
     &                             + d2drprime*dzibyic - d2drdyz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                             + d2edphi2*dphidzic*dphidzib
     &                             - (dzic-dzib)*ddrdz - d2drdzz
     &                             + dedphiprime*dphidzic*dphidzib
     &                             + d2edphi2prime*dphidzic*dphidzib 
     &                             + d2drprime*dzibzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                             + d2edphi2*dphidxic*dphidxid
     &                             + dedphiprime*dphidxic*dphidxid
     &                             + d2edphi2prime*dphidxic*dphidxid
     &                             + d2drprime*dxicxid + dxid*ddrdx
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                             + d2edphi2*dphidyic*dphidxid
     &                             + dedphiprime*dphidyic*dphidxid
     &                             + d2edphi2prime*dphidyic*dphidxid
     &                             + d2drprime*dyicxid + dxid*ddrdy
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                             + d2edphi2*dphidzic*dphidxid
     &                             + dedphiprime*dphidzic*dphidxid
     &                             + d2edphi2prime*dphidzic*dphidxid
     &                             + d2drprime*dzicxid + dxid*ddrdz
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                             + d2edphi2*dphidxic*dphidyid
     &                             + dedphiprime*dphidxic*dphidyid
     &                             + d2edphi2prime*dphidxic*dphidyid
     &                             + d2drprime*dxicyid + dyid*ddrdx
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                             + d2edphi2*dphidyic*dphidyid
     &                             + dedphiprime*dphidyic*dphidyid
     &                             + d2edphi2prime*dphidyic*dphidyid
     &                             + d2drprime*dyicyid + dyid*ddrdy
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                             + d2edphi2*dphidzic*dphidyid
     &                             + dedphiprime*dphidzic*dphidyid
     &                             + d2edphi2prime*dphidzic*dphidyid
     &                             + d2drprime*dzicyid + dyid*ddrdz
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                             + d2edphi2*dphidxic*dphidzid
     &                             + dedphiprime*dphidxic*dphidzid
     &                             + d2edphi2prime*dphidxic*dphidzid
     &                             + d2drprime*dxiczid + dzid*ddrdx
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                             + d2edphi2*dphidyic*dphidzid
     &                             + dedphiprime*dphidyic*dphidzid
     &                             + d2edphi2prime*dphidyic*dphidzid
     &                             + d2drprime*dyiczid + dzid*ddrdy
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                             + d2edphi2*dphidzic*dphidzid
     &                             + dedphiprime*dphidzic*dphidzid
     &                             + d2edphi2prime*dphidzic*dphidzid
     &                             + d2drprime*dziczid + dzid*ddrdz
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                             + d2edphi2*dphidxid*dphidxid
     &                             + dedphiprime*dphidxid*dphidxid
     &                             + d2edphi2prime*dphidxid*dphidxid
     &                             + d2drprime*dxidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + dedphiprime*dphidxid*dphidyid
     &                             + d2edphi2prime*dphidxid*dphidyid
     &                             + d2drprime*dxidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + dedphiprime*dphidxid*dphidzid
     &                             + d2edphi2prime*dphidxid*dphidzid
     &                             + d2drprime*dxidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + dedphiprime*dphidxid*dphidyid
     &                             + d2edphi2prime*dphidxid*dphidyid
     &                             + d2drprime*dxidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                             + d2edphi2*dphidyid*dphidyid
     &                             + dedphiprime*dphidyid*dphidyid
     &                             + d2edphi2prime*dphidyid*dphidyid
     &                             + d2drprime*dyidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + dedphiprime*dphidyid*dphidzid
     &                             + d2edphi2prime*dphidyid*dphidzid
     &                             + d2drprime*dyidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + dedphiprime*dphidxid*dphidzid
     &                             + d2edphi2prime*dphidxid*dphidzid
     &                             + d2drprime*dxidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + dedphiprime*dphidyid*dphidzid
     &                             + d2edphi2prime*dphidyid*dphidzid
     &                             + d2drprime*dyidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                             + d2edphi2*dphidzid*dphidzid
     &                             + dedphiprime*dphidzid*dphidzid
     &                             + d2edphi2prime*dphidzid*dphidzid
     &                             + d2drprime*dzidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxid*dphidxia
     &                             + dedphiprime*dphidxid*dphidxia
     &                             + d2edphi2prime*dphidxid*dphidxia
     &                             + d2drprime*dxiaxid
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                             + d2edphi2*dphidyid*dphidxia
     &                             + dedphiprime*dphidyid*dphidxia
     &                             + d2edphi2prime*dphidyid*dphidxia
     &                             + d2drprime*dxiayid
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                             + d2edphi2*dphidzid*dphidxia
     &                             + dedphiprime*dphidzid*dphidxia
     &                             + d2edphi2prime*dphidzid*dphidxia
     &                             + d2drprime*dxiazid
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                             + d2edphi2*dphidxid*dphidyia
     &                             + dedphiprime*dphidxid*dphidyia
     &                             + d2edphi2prime*dphidxid*dphidyia
     &                             + d2drprime*dyiaxid
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                             + d2edphi2*dphidyid*dphidyia
     &                             + dedphiprime*dphidyid*dphidyia
     &                             + d2edphi2prime*dphidyid*dphidyia
     &                             + d2drprime*dyiayid
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                             + d2edphi2*dphidzid*dphidyia
     &                             + dedphiprime*dphidzid*dphidyia
     &                             + d2edphi2prime*dphidzid*dphidyia
     &                             + d2drprime*dyiazid
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                             + d2edphi2*dphidxid*dphidzia
     &                             + dedphiprime*dphidxid*dphidzia
     &                             + d2edphi2prime*dphidxid*dphidzia
     &                             + d2drprime*dziaxid
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                             + d2edphi2*dphidyid*dphidzia
     &                             + dedphiprime*dphidyid*dphidzia
     &                             + d2edphi2prime*dphidyid*dphidzia
     &                             + d2drprime*dziayid
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                             + d2edphi2*dphidzid*dphidzia
     &                             + dedphiprime*dphidzid*dphidzia
     &                             + d2edphi2prime*dphidzid*dphidzia
     &                             + d2drprime*dziazid
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                             + d2edphi2*dphidxid*dphidxib
     &                             + dedphiprime*dphidxid*dphidxib
     &                             + d2edphi2prime*dphidxid*dphidxib
     &                             + d2drprime*dxibxid - dxid*ddrdx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                             + d2edphi2*dphidyid*dphidxib
     &                             + dedphiprime*dphidyid*dphidxib
     &                             + d2edphi2prime*dphidyid*dphidxib
     &                             + d2drprime*dxibyid - dyid*ddrdx
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                             + d2edphi2*dphidzid*dphidxib
     &                             + dedphiprime*dphidzid*dphidxib
     &                             + d2edphi2prime*dphidzid*dphidxib
     &                             + d2drprime*dxibzid - dzid*ddrdx
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                             + d2edphi2*dphidxid*dphidyib
     &                             + dedphiprime*dphidxid*dphidyib
     &                             + d2edphi2prime*dphidxid*dphidyib
     &                             + d2drprime*dyibxid - dxid*ddrdy
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                             + d2edphi2*dphidyid*dphidyib
     &                             + dedphiprime*dphidyid*dphidyib
     &                             + d2edphi2prime*dphidyid*dphidyib
     &                             + d2drprime*dyibyid - dyid*ddrdy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                             + d2edphi2*dphidzid*dphidyib
     &                             + dedphiprime*dphidzid*dphidyib
     &                             + d2edphi2prime*dphidzid*dphidyib 
     &                             + d2drprime*dyibzid - dzid*ddrdy
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                             + d2edphi2*dphidxid*dphidzib
     &                             + dedphiprime*dphidxid*dphidzib
     &                             + d2edphi2prime*dphidxid*dphidzib
     &                             + d2drprime*dzibxid - dxid*ddrdz
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                             + d2edphi2*dphidyid*dphidzib
     &                             + dedphiprime*dphidyid*dphidzib
     &                             + d2edphi2prime*dphidyid*dphidzib
     &                             + d2drprime*dzibyid - dyid*ddrdz
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                             + d2edphi2*dphidzid*dphidzib
     &                             + dedphiprime*dphidzid*dphidzib
     &                             + d2edphi2prime*dphidzid*dphidzib
     &                             + d2drprime*dzibzid - dzid*ddrdz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                             + d2edphi2*dphidxid*dphidxic
     &                             + dedphiprime*dphidxid*dphidxic
     &                             + d2edphi2prime*dphidxid*dphidxic
     &                             + d2drprime*dxicxid + dxid*ddrdx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                             + d2edphi2*dphidyid*dphidxic
     &                             + dedphiprime*dphidyid*dphidxic
     &                             + d2edphi2prime*dphidyid*dphidxic
     &                             + d2drprime*dxicyid + dyid*ddrdx
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                             + d2edphi2*dphidzid*dphidxic
     &                             + dedphiprime*dphidzid*dphidxic
     &                             + d2edphi2prime*dphidzid*dphidxic
     &                             + d2drprime*dxiczid + dzid*ddrdx
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                             + d2edphi2*dphidxid*dphidyic
     &                             + dedphiprime*dphidxid*dphidyic
     &                             + d2edphi2prime*dphidxid*dphidyic
     &                             + d2drprime*dyicxid + dxid*ddrdy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                             + d2edphi2*dphidyid*dphidyic
     &                             + dedphiprime*dphidyid*dphidyic
     &                             + d2edphi2prime*dphidyid*dphidyic
     &                             + d2drprime*dyicyid + dyid*ddrdy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                             + d2edphi2*dphidzid*dphidyic
     &                             + dedphiprime*dphidzid*dphidyic
     &                             + d2edphi2prime*dphidzid*dphidyic
     &                             + d2drprime*dyiczid + dzid*ddrdy
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                             + d2edphi2*dphidxid*dphidzic
     &                             + dedphiprime*dphidxid*dphidzic
     &                             + d2edphi2prime*dphidxid*dphidzic
     &                             + d2drprime*dzicxid + dxid*ddrdz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                             + d2edphi2*dphidyid*dphidzic
     &                             + dedphiprime*dphidyid*dphidzic
     &                             + d2edphi2prime*dphidyid*dphidzic
     &                             + d2drprime*dzicyid + dyid*ddrdz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                             + d2edphi2*dphidzid*dphidzic
     &                             + dedphiprime*dphidzid*dphidzic
     &                             + d2edphi2prime*dphidzid*dphidzic
     &                             + d2drprime*dziczid + dzid*ddrdz
               end if
c
c     get the stretch-torsion values for the third bond
c
               v1 = kst(7,istrtor)
               v2 = kst(8,istrtor)
               v3 = kst(9,istrtor)
               k = ist(4,istrtor)
               dr = rdc - bl(k)
               force = bk(k)
               dedphi = storunit * 2.0d0 * force 
     &                     * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2 = storunit * 2.0d0 * force * dr
     &                       * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3)
               ddr = 1.0d0 / rdc
               d2dr = -storunit * 2.0d0 * force 
     &                    * (v1*phi1 + v2*phi2 + v3*phi3) / rdc**3
               dedphiprime = storunit * 2.0d0 * force
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2prime = storunit * 2.0d0 * force
     &                            * (v1*phi1 + v2*phi2 + v3*phi3)
     &                            * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3)
               d2drprime = storunit * 2.0d0 * force
     &                        * (v1*phi1 + v2*phi2 + v3*phi3)
     &                        * (v1*dphi1 + v2*dphi2 + v3*dphi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
                  d2dr = d2dr * fgrp
               end if
c
c     compute derivative components for this interaction
c
               ddrdx = xdc * ddr
               ddrdy = ydc * ddr
               ddrdz = zdc * ddr
               d2drdxx = (xdc*xdc-rdc*rdc) * d2dr
               d2drdyy = (ydc*ydc-rdc*rdc) * d2dr
               d2drdzz = (zdc*zdc-rdc*rdc) * d2dr
               d2drdxy = xdc * ydc * d2dr
               d2drdxz = xdc * zdc * d2dr
               d2drdyz = ydc * zdc * d2dr
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
c     increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                             + d2edphi2*dphidxia*dphidxia
     &                             + dedphiprime*dphidxia*dphidxia
     &                             + d2edphi2prime*dphidxia*dphidxia
     &                             + d2drprime*dxiaxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             + dedphiprime*dphidxia*dphidyia
     &                             + d2edphi2prime*dphidxia*dphidyia
     &                             + d2drprime*dxiayia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             + dedphiprime*dphidxia*dphidzia
     &                             + d2edphi2prime*dphidxia*dphidzia
     &                             + d2drprime*dxiazia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             + dedphiprime*dphidxia*dphidyia
     &                             + d2edphi2prime*dphidxia*dphidyia
     &                             + d2drprime*dxiayia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                             + d2edphi2*dphidyia*dphidyia
     &                             + dedphiprime*dphidyia*dphidyia
     &                             + d2edphi2prime*dphidyia*dphidyia 
     &                             + d2drprime*dyiayia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             + dedphiprime*dphidyia*dphidzia
     &                             + d2edphi2prime*dphidyia*dphidzia
     &                             + d2drprime*dyiazia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             + dedphiprime*dphidxia*dphidzia
     &                             + d2edphi2prime*dphidxia*dphidzia
     &                             + d2drprime*dxiazia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             + dedphiprime*dphidyia*dphidzia
     &                             + d2edphi2prime*dphidyia*dphidzia
     &                             + d2drprime*dyiazia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                             + d2edphi2*dphidzia*dphidzia
     &                             + dedphiprime*dphidzia*dphidzia
     &                             + d2edphi2prime*dphidzia*dphidzia
     &                             + d2drprime*dziazia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxia*dphidxib
     &                             + dedphiprime*dphidxia*dphidxib
     &                             + d2edphi2prime*dphidxia*dphidxib
     &                             + d2drprime*dxiaxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                             + d2edphi2*dphidyia*dphidxib
     &                             + dedphiprime*dphidyia*dphidxib
     &                             + d2edphi2prime*dphidyia*dphidxib
     &                             + d2drprime*dyiaxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                             + d2edphi2*dphidzia*dphidxib
     &                             + dedphiprime*dphidzia*dphidxib
     &                             + d2edphi2prime*dphidzia*dphidxib
     &                             + d2drprime*dziaxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                             + d2edphi2*dphidxia*dphidyib
     &                             + dedphiprime*dphidxia*dphidyib
     &                             + d2edphi2prime*dphidxia*dphidyib
     &                             + d2drprime*dxiayib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                             + d2edphi2*dphidyia*dphidyib
     &                             + dedphiprime*dphidyia*dphidyib
     &                             + d2edphi2prime*dphidyia*dphidyib
     &                             + d2drprime*dyiayib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                             + d2edphi2*dphidzia*dphidyib
     &                             + dedphiprime*dphidzia*dphidyib
     &                             + d2edphi2prime*dphidzia*dphidyib
     &                             + d2drprime*dziayib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                             + d2edphi2*dphidxia*dphidzib
     &                             + dedphiprime*dphidxia*dphidzib
     &                             + d2edphi2prime*dphidxia*dphidzib
     &                             + d2drprime*dxiazib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                             + d2edphi2*dphidyia*dphidzib
     &                             + dedphiprime*dphidyia*dphidzib
     &                             + d2edphi2prime*dphidyia*dphidzib
     &                             + d2drprime*dyiazib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                             + d2edphi2*dphidzia*dphidzib
     &                             + dedphiprime*dphidzia*dphidzib
     &                             + d2edphi2prime*dphidzia*dphidzib
     &                             + d2drprime*dziazib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxia*dphidxic
     &                             + dedphiprime*dphidxia*dphidxic
     &                             + d2edphi2prime*dphidxia*dphidxic
     &                             + d2drprime*dxiaxic - dxia*ddrdx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                             + d2edphi2*dphidyia*dphidxic
     &                             + dedphiprime*dphidyia*dphidxic
     &                             + d2edphi2prime*dphidyia*dphidxic
     &                             + d2drprime*dyiaxic - dyia*ddrdx
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                             + d2edphi2*dphidzia*dphidxic
     &                             + dedphiprime*dphidzia*dphidxic
     &                             + d2edphi2prime*dphidzia*dphidxic
     &                             + d2drprime*dziaxic - dzia*ddrdx
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                             + d2edphi2*dphidxia*dphidyic
     &                             + dedphiprime*dphidxia*dphidyic
     &                             + d2edphi2prime*dphidxia*dphidyic
     &                             + d2drprime*dxiayic - dxia*ddrdy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                             + d2edphi2*dphidyia*dphidyic
     &                             + dedphiprime*dphidyia*dphidyic
     &                             + d2edphi2prime*dphidyia*dphidyic
     &                             + d2drprime*dyiayic - dyia*ddrdy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                             + d2edphi2*dphidzia*dphidyic
     &                             + dedphiprime*dphidzia*dphidyic
     &                             + d2edphi2prime*dphidzia*dphidyic
     &                             + d2drprime*dziayic - dzia*ddrdy
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                             + d2edphi2*dphidxia*dphidzic
     &                             + dedphiprime*dphidxia*dphidzic
     &                             + d2edphi2prime*dphidxia*dphidzic
     &                             + d2drprime*dxiazic - dxia*ddrdz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                             + d2edphi2*dphidyia*dphidzic
     &                             + dedphiprime*dphidyia*dphidzic
     &                             + d2edphi2prime*dphidyia*dphidzic
     &                             + d2drprime*dyiazic - dyia*ddrdz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                             + d2edphi2*dphidzia*dphidzic
     &                             + dedphiprime*dphidzia*dphidzic
     &                             + d2edphi2prime*dphidzia*dphidzic
     &                             + d2drprime*dziazic - dzia*ddrdz
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxia*dphidxid
     &                             + dedphiprime*dphidxia*dphidxid
     &                             + d2edphi2prime*dphidxia*dphidxid
     &                             + d2drprime*dxiaxid + dxia*ddrdx
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                             + d2edphi2*dphidyia*dphidxid
     &                             + dedphiprime*dphidyia*dphidxid
     &                             + d2edphi2prime*dphidyia*dphidxid
     &                             + d2drprime*dyiaxid + dyia*ddrdx
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                             + d2edphi2*dphidzia*dphidxid
     &                             + dedphiprime*dphidzia*dphidxid
     &                             + d2edphi2prime*dphidzia*dphidxid
     &                             + d2drprime*dziaxid + dzia*ddrdx
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                             + d2edphi2*dphidxia*dphidyid
     &                             + dedphiprime*dphidxia*dphidyid
     &                             + d2edphi2prime*dphidxia*dphidyid
     &                             + d2drprime*dxiayid + dxia*ddrdy
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                             + d2edphi2*dphidyia*dphidyid
     &                             + dedphiprime*dphidyia*dphidyid
     &                             + d2edphi2prime*dphidyia*dphidyid
     &                             + d2drprime*dyiayid + dyia*ddrdy
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                             + d2edphi2*dphidzia*dphidyid
     &                             + dedphiprime*dphidzia*dphidyid
     &                             + d2edphi2prime*dphidzia*dphidyid
     &                             + d2drprime*dziayid + dzia*ddrdy
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                             + d2edphi2*dphidxia*dphidzid
     &                             + dedphiprime*dphidxia*dphidzid
     &                             + d2edphi2prime*dphidxia*dphidzid
     &                             + d2drprime*dxiazid + dxia*ddrdz
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                             + d2edphi2*dphidyia*dphidzid
     &                             + dedphiprime*dphidyia*dphidzid
     &                             + d2edphi2prime*dphidyia*dphidzid
     &                             + d2drprime*dyiazid + dyia*ddrdz
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                             + d2edphi2*dphidzia*dphidzid
     &                             + dedphiprime*dphidzia*dphidzid
     &                             + d2edphi2prime*dphidzia*dphidzid
     &                             + d2drprime*dziazid + dzia*ddrdz
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib
     &                             + d2edphi2*dphidxib*dphidxib
     &                             + dedphiprime*dphidxib*dphidxib
     &                             + d2edphi2prime*dphidxib*dphidxib
     &                             + d2drprime*dxibxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + dedphiprime*dphidxib*dphidyib
     &                             + d2edphi2prime*dphidxib*dphidyib
     &                             + d2drprime*dxibyib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + dedphiprime*dphidxib*dphidzib
     &                             + d2edphi2prime*dphidxib*dphidzib
     &                             + d2drprime*dxibzib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + dedphiprime*dphidxib*dphidyib
     &                             + d2edphi2prime*dphidxib*dphidyib
     &                             + d2drprime*dxibyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib
     &                             + d2edphi2*dphidyib*dphidyib
     &                             + dedphiprime*dphidyib*dphidyib
     &                             + d2edphi2prime*dphidyib*dphidyib
     &                             + d2drprime*dyibyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + dedphiprime*dphidyib*dphidzib
     &                             + d2edphi2prime*dphidyib*dphidzib
     &                             + d2drprime*dyibzib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + dedphiprime*dphidxib*dphidzib
     &                             + d2edphi2prime*dphidxib*dphidzib
     &                             + d2drprime*dxibzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + dedphiprime*dphidyib*dphidzib
     &                             + d2edphi2prime*dphidyib*dphidzib
     &                             + d2drprime*dyibzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib
     &                             + d2edphi2*dphidzib*dphidzib
     &                             + dedphiprime*dphidzib*dphidzib
     &                             + d2edphi2prime*dphidzib*dphidzib
     &                             + d2drprime*dzibzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxib*dphidxia
     &                             + dedphiprime*dphidxib*dphidxia
     &                             + d2edphi2prime*dphidxib*dphidxia
     &                             + d2drprime*dxiaxib
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                             + d2edphi2*dphidyib*dphidxia
     &                             + dedphiprime*dphidyib*dphidxia
     &                             + d2edphi2prime*dphidyib*dphidxia
     &                             + d2drprime*dxiayib
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                             + d2edphi2*dphidzib*dphidxia
     &                             + dedphiprime*dphidzib*dphidxia
     &                             + d2edphi2prime*dphidzib*dphidxia
     &                             + d2drprime*dxiazib
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                             + d2edphi2*dphidxib*dphidyia
     &                             + dedphiprime*dphidxib*dphidyia
     &                             + d2edphi2prime*dphidxib*dphidyia
     &                             + d2drprime*dyiaxib
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                             + d2edphi2*dphidyib*dphidyia
     &                             + dedphiprime*dphidyib*dphidyia
     &                             + d2edphi2prime*dphidyib*dphidyia
     &                             + d2drprime*dyiayib
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                             + d2edphi2*dphidzib*dphidyia
     &                             + dedphiprime*dphidzib*dphidyia
     &                             + d2edphi2prime*dphidzib*dphidyia
     &                             + d2drprime*dyiazib
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                             + d2edphi2*dphidxib*dphidzia
     &                             + dedphiprime*dphidxib*dphidzia
     &                             + d2edphi2prime*dphidxib*dphidzia
     &                             + d2drprime*dziaxib
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                             + d2edphi2*dphidyib*dphidzia
     &                             + dedphiprime*dphidyib*dphidzia
     &                             + d2edphi2prime*dphidyib*dphidzia
     &                             + d2drprime*dziayib
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                             + d2edphi2*dphidzib*dphidzia
     &                             + dedphiprime*dphidzib*dphidzia
     &                             + d2edphi2prime*dphidzib*dphidzia
     &                             + d2drprime*dziazib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic
     &                             + d2edphi2*dphidxib*dphidxic
     &                             + dedphiprime*dphidxib*dphidxic
     &                             + d2edphi2prime*dphidxib*dphidxic
     &                             + d2drprime*dxibxic - dxib*ddrdx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic
     &                             + d2edphi2*dphidyib*dphidxic
     &                             + dedphiprime*dphidyib*dphidxic
     &                             + d2edphi2prime*dphidyib*dphidxic
     &                             + d2drprime*dyibxic - dyib*ddrdx
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic
     &                             + d2edphi2*dphidzib*dphidxic
     &                             + dedphiprime*dphidzib*dphidxic
     &                             + d2edphi2prime*dphidzib*dphidxic
     &                             + d2drprime*dzibxic - dzib*ddrdx
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic
     &                             + d2edphi2*dphidxib*dphidyic
     &                             + dedphiprime*dphidxib*dphidyic
     &                             + d2edphi2prime*dphidxib*dphidyic
     &                             + d2drprime*dxibyic - dxib*ddrdy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic
     &                             + d2edphi2*dphidyib*dphidyic
     &                             + dedphiprime*dphidyib*dphidyic
     &                             + d2edphi2prime*dphidyib*dphidyic
     &                             + d2drprime*dyibyic - dyib*ddrdy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic
     &                             + d2edphi2*dphidzib*dphidyic
     &                             + dedphiprime*dphidzib*dphidyic
     &                             + d2edphi2prime*dphidzib*dphidyic
     &                             + d2drprime*dzibyic - dzib*ddrdy
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic
     &                             + d2edphi2*dphidxib*dphidzic
     &                             + dedphiprime*dphidxib*dphidzic
     &                             + d2edphi2prime*dphidxib*dphidzic
     &                             + d2drprime*dxibzic - dxib*ddrdz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic
     &                             + d2edphi2*dphidyib*dphidzic
     &                             + dedphiprime*dphidyib*dphidzic
     &                             + d2edphi2prime*dphidyib*dphidzic
     &                             + d2drprime*dyibzic - dyib*ddrdz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic
     &                             + d2edphi2*dphidzib*dphidzic
     &                             + dedphiprime*dphidzib*dphidzic
     &                             + d2edphi2prime*dphidzib*dphidzic
     &                             + d2drprime*dzibzic - dzib*ddrdz
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                             + d2edphi2*dphidxib*dphidxid
     &                             + dedphiprime*dphidxib*dphidxid
     &                             + d2edphi2prime*dphidxib*dphidxid
     &                             + d2drprime*dxibxid + dxib*ddrdx
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                             + d2edphi2*dphidyib*dphidxid
     &                             + dedphiprime*dphidyib*dphidxid
     &                             + d2edphi2prime*dphidyib*dphidxid
     &                             + d2drprime*dyibxid + dyib*ddrdx
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                             + d2edphi2*dphidzib*dphidxid
     &                             + dedphiprime*dphidzib*dphidxid
     &                             + d2edphi2prime*dphidzib*dphidxid
     &                             + d2drprime*dzibxid + dzib*ddrdx
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                             + d2edphi2*dphidxib*dphidyid
     &                             + dedphiprime*dphidxib*dphidyid
     &                             + d2edphi2prime*dphidxib*dphidyid
     &                             + d2drprime*dxibyid + dxib*ddrdy
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                             + d2edphi2*dphidyib*dphidyid
     &                             + dedphiprime*dphidyib*dphidyid
     &                             + d2edphi2prime*dphidyib*dphidyid
     &                             + d2drprime*dyibyid + dyib*ddrdy
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                             + d2edphi2*dphidzib*dphidyid
     &                             + dedphiprime*dphidzib*dphidyid
     &                             + d2edphi2prime*dphidzib*dphidyid
     &                             + d2drprime*dzibyid + dzib*ddrdy
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                             + d2edphi2*dphidxib*dphidzid
     &                             + dedphiprime*dphidxib*dphidzid
     &                             + d2edphi2prime*dphidxib*dphidzid
     &                             + d2drprime*dxibzid + dxib*ddrdz
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                             + d2edphi2*dphidyib*dphidzid
     &                             + dedphiprime*dphidyib*dphidzid
     &                             + d2edphi2prime*dphidyib*dphidzid
     &                             + d2drprime*dyibzid + dyib*ddrdz
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                             + d2edphi2*dphidzib*dphidzid
     &                             + dedphiprime*dphidzib*dphidzid
     &                             + d2edphi2prime*dphidzib*dphidzid
     &                             + d2drprime*dzibzid + dzib*ddrdz
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic
     &                             + d2edphi2*dphidxic*dphidxic
     &                             - 2.0d0*dxic*ddrdx + d2drdxx
     &                             + dedphiprime*dphidxic*dphidxic
     &                             + d2edphi2prime*dphidxic*dphidxic
     &                             + d2drprime*dxicxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             - dyic*ddrdx - dxic*ddrdy
     &                             + dedphiprime*dphidxic*dphidyic
     &                             + d2edphi2prime*dphidxic*dphidyic
     &                             + d2drprime*dxicyic + d2drdxy
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             - dzic*ddrdx - dxic*ddrdz
     &                             + dedphiprime*dphidxic*dphidzic
     &                             + d2edphi2prime*dphidxic*dphidzic
     &                             + d2drprime*dxiczic + d2drdxz
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             - dxic*ddrdy - dyic*ddrdx
     &                             + dedphiprime*dphidxic*dphidyic
     &                             + d2edphi2prime*dphidxic*dphidyic
     &                             + d2drprime*dxicyic + d2drdxy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic
     &                             + d2edphi2*dphidyic*dphidyic
     &                             - 2.0d0*dyic*ddrdy + d2drdyy
     &                             + dedphiprime*dphidyic*dphidyic
     &                             + d2edphi2prime*dphidyic*dphidyic
     &                             + d2drprime*dyicyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             - dzic*ddrdy - dyic*ddrdz
     &                             + dedphiprime*dphidyic*dphidzic
     &                             + d2edphi2prime*dphidyic*dphidzic
     &                             + d2drprime*dyiczic + d2drdyz
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             - dxic*ddrdz - dzic*ddrdx
     &                             + dedphiprime*dphidxic*dphidzic
     &                             + d2edphi2prime*dphidxic*dphidzic
     &                             + d2drprime*dxiczic + d2drdxz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             - dyic*ddrdz - dzic*ddrdy
     &                             + dedphiprime*dphidyic*dphidzic
     &                             + d2edphi2prime*dphidyic*dphidzic
     &                             + d2drprime*dyiczic + d2drdyz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic
     &                             + d2edphi2*dphidzic*dphidzic
     &                             - 2.0d0*dzic*ddrdz + d2drdzz
     &                             + dedphiprime*dphidzic*dphidzic
     &                             + d2edphi2prime*dphidzic*dphidzic
     &                             + d2drprime*dziczic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxic*dphidxia
     &                             + dedphiprime*dphidxic*dphidxia
     &                             + d2edphi2prime*dphidxic*dphidxia
     &                             + d2drprime*dxiaxic - dxia*ddrdx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                             + d2edphi2*dphidyic*dphidxia
     &                             + dedphiprime*dphidyic*dphidxia
     &                             + d2edphi2prime*dphidyic*dphidxia
     &                             + d2drprime*dxiayic - dxia*ddrdy
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                             + d2edphi2*dphidzic*dphidxia
     &                             + dedphiprime*dphidzic*dphidxia
     &                             + d2edphi2prime*dphidzic*dphidxia
     &                             + d2drprime*dxiazic - dxia*ddrdz
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                             + d2edphi2*dphidxic*dphidyia
     &                             + dedphiprime*dphidxic*dphidyia
     &                             + d2edphi2prime*dphidxic*dphidyia
     &                             + d2drprime*dyiaxic - dyia*ddrdx
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                             + d2edphi2*dphidyic*dphidyia
     &                             + dedphiprime*dphidyic*dphidyia
     &                             + d2edphi2prime*dphidyic*dphidyia
     &                             + d2drprime*dyiayic - dyia*ddrdy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                             + d2edphi2*dphidzic*dphidyia
     &                             + dedphiprime*dphidzic*dphidyia
     &                             + d2edphi2prime*dphidzic*dphidyia
     &                             + d2drprime*dyiazic - dyia*ddrdz
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                             + d2edphi2*dphidxic*dphidzia
     &                             + dedphiprime*dphidxic*dphidzia
     &                             + d2edphi2prime*dphidxic*dphidzia
     &                             + d2drprime*dziaxic - dzia*ddrdx
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                             + d2edphi2*dphidyic*dphidzia
     &                             + dedphiprime*dphidyic*dphidzia
     &                             + d2edphi2prime*dphidyic*dphidzia
     &                             + d2drprime*dziayic - dzia*ddrdy
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                             + d2edphi2*dphidzic*dphidzia
     &                             + dedphiprime*dphidzic*dphidzia
     &                             + d2edphi2prime*dphidzic*dphidzia
     &                             + d2drprime*dziazic - dzia*ddrdz
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                             + d2edphi2*dphidxic*dphidxib
     &                             + dedphiprime*dphidxic*dphidxib
     &                             + d2edphi2prime*dphidxic*dphidxib
     &                             + d2drprime*dxibxic - dxib*ddrdx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic
     &                             + d2edphi2*dphidyic*dphidxib
     &                             + dedphiprime*dphidyic*dphidxib
     &                             + d2edphi2prime*dphidyic*dphidxib 
     &                             + d2drprime*dxibyic - dxib*ddrdy
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic
     &                             + d2edphi2*dphidzic*dphidxib
     &                             + dedphiprime*dphidzic*dphidxib
     &                             + d2edphi2prime*dphidzic*dphidxib
     &                             + d2drprime*dxibzic - dxib*ddrdz
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic
     &                             + d2edphi2*dphidxic*dphidyib
     &                             + dedphiprime*dphidxic*dphidyib
     &                             + d2edphi2prime*dphidxic*dphidyib
     &                             + d2drprime*dyibxic - dyib*ddrdx
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic
     &                             + d2edphi2*dphidyic*dphidyib
     &                             + dedphiprime*dphidyic*dphidyib
     &                             + d2edphi2prime*dphidyic*dphidyib
     &                             + d2drprime*dyibyic - dyib*ddrdy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic
     &                             + d2edphi2*dphidzic*dphidyib
     &                             + dedphiprime*dphidzic*dphidyib
     &                             + d2edphi2prime*dphidzic*dphidyib
     &                             + d2drprime*dyibzic - dyib*ddrdz
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic
     &                             + d2edphi2*dphidxic*dphidzib
     &                             + dedphiprime*dphidxic*dphidzib
     &                             + d2edphi2prime*dphidxic*dphidzib
     &                             + d2drprime*dzibxic - dzib*ddrdx
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic
     &                             + d2edphi2*dphidyic*dphidzib
     &                             + dedphiprime*dphidyic*dphidzib
     &                             + d2edphi2prime*dphidyic*dphidzib
     &                             + d2drprime*dzibyic - dzib*ddrdy
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic
     &                             + d2edphi2*dphidzic*dphidzib
     &                             + dedphiprime*dphidzic*dphidzib
     &                             + d2edphi2prime*dphidzic*dphidzib
     &                             + d2drprime*dzibzic - dzib*ddrdz
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                             + d2edphi2*dphidxic*dphidxid
     &                             - dxid*ddrdx + dxic*ddrdx
     &                             + dedphiprime*dphidxic*dphidxid
     &                             + d2edphi2prime*dphidxic*dphidxid
     &                             + d2drprime*dxicxid - d2drdxx
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                             + d2edphi2*dphidyic*dphidxid
     &                             - dxid*ddrdy + dyic*ddrdx
     &                             + dedphiprime*dphidyic*dphidxid
     &                             + d2edphi2prime*dphidyic*dphidxid
     &                             + d2drprime*dyicxid - d2drdxy
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                             + d2edphi2*dphidzic*dphidxid
     &                             - dxid*ddrdz + dzic*ddrdx
     &                             + dedphiprime*dphidzic*dphidxid
     &                             + d2edphi2prime*dphidzic*dphidxid
     &                             + d2drprime*dzicxid - d2drdxz
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                             + d2edphi2*dphidxic*dphidyid
     &                             - dyid*ddrdx + dxic*ddrdy
     &                             + dedphiprime*dphidxic*dphidyid
     &                             + d2edphi2prime*dphidxic*dphidyid
     &                             + d2drprime*dxicyid - d2drdxy
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                             + d2edphi2*dphidyic*dphidyid
     &                             - dyid*ddrdy + dyic*ddrdy
     &                             + dedphiprime*dphidyic*dphidyid
     &                             + d2edphi2prime*dphidyic*dphidyid
     &                             + d2drprime*dyicyid - d2drdyy
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                             + d2edphi2*dphidzic*dphidyid
     &                             - dyid*ddrdz + dzic*ddrdy
     &                             + dedphiprime*dphidzic*dphidyid
     &                             + d2edphi2prime*dphidzic*dphidyid
     &                             + d2drprime*dzicyid - d2drdyz
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                             + d2edphi2*dphidxic*dphidzid
     &                             - dzid*ddrdx + dxic*ddrdz
     &                             + dedphiprime*dphidxic*dphidzid
     &                             + d2edphi2prime*dphidxic*dphidzid
     &                             + d2drprime*dxiczid - d2drdxz
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                             + d2edphi2*dphidyic*dphidzid
     &                             - dzid*ddrdy + dyic*ddrdz
     &                             + dedphiprime*dphidyic*dphidzid
     &                             + d2edphi2prime*dphidyic*dphidzid
     &                             + d2drprime*dyiczid - d2drdyz
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                             + d2edphi2*dphidzic*dphidzid
     &                             - dzid*ddrdz + dzic*ddrdz
     &                             + dedphiprime*dphidzic*dphidzid
     &                             + d2edphi2prime*dphidzic*dphidzid
     &                             + d2drprime*dziczid - d2drdzz
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                             + d2edphi2*dphidxid*dphidxid
     &                             + 2.0d0*dxid*ddrdx
     &                             + dedphiprime*dphidxid*dphidxid
     &                             + d2edphi2prime*dphidxid*dphidxid
     &                             + d2drprime*dxidxid + d2drdxx
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + dxid*ddrdy + dyid*ddrdx
     &                             + dedphiprime*dphidxid*dphidyid
     &                             + d2edphi2prime*dphidxid*dphidyid
     &                             + d2drprime*dxidyid + d2drdxy
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + dxid*ddrdz + dzid*ddrdx
     &                             + dedphiprime*dphidxid*dphidzid
     &                             + d2edphi2prime*dphidxid*dphidzid
     &                             + d2drprime*dxidzid + d2drdxz
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + dxid*ddrdy + dyid*ddrdx
     &                             + dedphiprime*dphidxid*dphidyid
     &                             + d2edphi2prime*dphidxid*dphidyid
     &                             + d2drprime*dxidyid + d2drdxy
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                             + d2edphi2*dphidyid*dphidyid
     &                             + 2.0d0*dyid*ddrdy
     &                             + dedphiprime*dphidyid*dphidyid
     &                             + d2edphi2prime*dphidyid*dphidyid
     &                             + d2drprime*dyidyid + d2drdyy
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + dyid*ddrdz + dzid*ddrdy
     &                             + dedphiprime*dphidyid*dphidzid
     &                             + d2edphi2prime*dphidyid*dphidzid
     &                             + d2drprime*dyidzid + d2drdyz
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + dxid*ddrdz + dzid*ddrdx
     &                             + dedphiprime*dphidxid*dphidzid
     &                             + d2edphi2prime*dphidxid*dphidzid
     &                             + d2drprime*dxidzid + d2drdxz
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + dyid*ddrdz + dzid*ddrdy
     &                             + dedphiprime*dphidyid*dphidzid
     &                             + d2edphi2prime*dphidyid*dphidzid
     &                             + d2drprime*dyidzid + d2drdyz
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                             + d2edphi2*dphidzid*dphidzid
     &                             + 2.0d0*dzid*ddrdz
     &                             + dedphiprime*dphidzid*dphidzid
     &                             + d2edphi2prime*dphidzid*dphidzid
     &                             + d2drprime*dzidzid + d2drdzz
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxid*dphidxia
     &                             + dedphiprime*dphidxid*dphidxia
     &                             + d2edphi2prime*dphidxid*dphidxia
     &                             + d2drprime*dxiaxid + dxia*ddrdx
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                             + d2edphi2*dphidyid*dphidxia
     &                             + dedphiprime*dphidyid*dphidxia
     &                             + d2edphi2prime*dphidyid*dphidxia
     &                             + d2drprime*dxiayid + dxia*ddrdy
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                             + d2edphi2*dphidzid*dphidxia
     &                             + dedphiprime*dphidzid*dphidxia
     &                             + d2edphi2prime*dphidzid*dphidxia
     &                             + d2drprime*dxiazid + dxia*ddrdz
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                             + d2edphi2*dphidxid*dphidyia
     &                             + dedphiprime*dphidxid*dphidyia
     &                             + d2edphi2prime*dphidxid*dphidyia
     &                             + d2drprime*dyiaxid + dyia*ddrdx
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                             + d2edphi2*dphidyid*dphidyia
     &                             + dedphiprime*dphidyid*dphidyia
     &                             + d2edphi2prime*dphidyid*dphidyia
     &                             + d2drprime*dyiayid + dyia*ddrdy
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                             + d2edphi2*dphidzid*dphidyia
     &                             + dedphiprime*dphidzid*dphidyia
     &                             + d2edphi2prime*dphidzid*dphidyia
     &                             + d2drprime*dyiazid + dyia*ddrdz
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                             + d2edphi2*dphidxid*dphidzia
     &                             + dedphiprime*dphidxid*dphidzia
     &                             + d2edphi2prime*dphidxid*dphidzia
     &                             + d2drprime*dziaxid + dzia*ddrdx
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                             + d2edphi2*dphidyid*dphidzia
     &                             + dedphiprime*dphidyid*dphidzia
     &                             + d2edphi2prime*dphidyid*dphidzia
     &                             + d2drprime*dziayid + dzia*ddrdy
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                             + d2edphi2*dphidzid*dphidzia
     &                             + dedphiprime*dphidzid*dphidzia
     &                             + d2edphi2prime*dphidzid*dphidzia
     &                             + d2drprime*dziazid + dzia*ddrdz
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                             + d2edphi2*dphidxid*dphidxib
     &                             + dedphiprime*dphidxid*dphidxib
     &                             + d2edphi2prime*dphidxid*dphidxib
     &                             + d2drprime*dxibxid + dxib*ddrdx
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                             + d2edphi2*dphidyid*dphidxib
     &                             + dedphiprime*dphidyid*dphidxib
     &                             + d2edphi2prime*dphidyid*dphidxib
     &                             + d2drprime*dxibyid + dxib*ddrdy
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                             + d2edphi2*dphidzid*dphidxib
     &                             + dedphiprime*dphidzid*dphidxib
     &                             + d2edphi2prime*dphidzid*dphidxib
     &                             + d2drprime*dxibzid + dxib*ddrdz
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                             + d2edphi2*dphidxid*dphidyib
     &                             + dedphiprime*dphidxid*dphidyib
     &                             + d2edphi2prime*dphidxid*dphidyib
     &                             + d2drprime*dyibxid + dyib*ddrdx
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                             + d2edphi2*dphidyid*dphidyib
     &                             + dedphiprime*dphidyid*dphidyib
     &                             + d2edphi2prime*dphidyid*dphidyib
     &                             + d2drprime*dyibyid + dyib*ddrdy
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                             + d2edphi2*dphidzid*dphidyib
     &                             + dedphiprime*dphidzid*dphidyib
     &                             + d2edphi2prime*dphidzid*dphidyib
     &                             + d2drprime*dyibzid + dyib*ddrdz
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                             + d2edphi2*dphidxid*dphidzib
     &                             + dedphiprime*dphidxid*dphidzib
     &                             + d2edphi2prime*dphidxid*dphidzib
     &                             + d2drprime*dzibxid + dzib*ddrdx
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                             + d2edphi2*dphidyid*dphidzib
     &                             + dedphiprime*dphidyid*dphidzib
     &                             + d2edphi2prime*dphidyid*dphidzib
     &                             + d2drprime*dzibyid + dzib*ddrdy
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                             + d2edphi2*dphidzid*dphidzib
     &                             + dedphiprime*dphidzid*dphidzib
     &                             + d2edphi2prime*dphidzid*dphidzib
     &                             + d2drprime*dzibzid + dzib*ddrdz
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                             + d2edphi2*dphidxid*dphidxic
     &                             + dxic*ddrdx - dxid*ddrdx
     &                             + dedphiprime*dphidxid*dphidxic
     &                             + d2edphi2prime*dphidxid*dphidxic
     &                             + d2drprime*dxicxid - d2drdxx
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                             + d2edphi2*dphidyid*dphidxic
     &                             + dxic*ddrdy - dyid*ddrdx
     &                             + dedphiprime*dphidyid*dphidxic
     &                             + d2edphi2prime*dphidyid*dphidxic
     &                             + d2drprime*dxicyid - d2drdxy
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                             + d2edphi2*dphidzid*dphidxic
     &                             + dxic*ddrdz - dzid*ddrdx
     &                             + dedphiprime*dphidzid*dphidxic
     &                             + d2edphi2prime*dphidzid*dphidxic
     &                             + d2drprime*dxiczid - d2drdxz
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                             + d2edphi2*dphidxid*dphidyic
     &                             + dyic*ddrdx - dxid*ddrdy
     &                             + dedphiprime*dphidxid*dphidyic
     &                             + d2edphi2prime*dphidxid*dphidyic
     &                             + d2drprime*dyicxid - d2drdxy
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                             + d2edphi2*dphidyid*dphidyic
     &                             + dyic*ddrdy - dyid*ddrdy
     &                             + dedphiprime*dphidyid*dphidyic
     &                             + d2edphi2prime*dphidyid*dphidyic
     &                             + d2drprime*dyicyid - d2drdyy
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                             + d2edphi2*dphidzid*dphidyic
     &                             + dyic*ddrdz - dzid*ddrdy
     &                             + dedphiprime*dphidzid*dphidyic
     &                             + d2edphi2prime*dphidzid*dphidyic
     &                             + d2drprime*dyiczid - d2drdyz
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                             + d2edphi2*dphidxid*dphidzic
     &                             + dzic*ddrdx - dxid*ddrdz
     &                             + dedphiprime*dphidxid*dphidzic
     &                             + d2edphi2prime*dphidxid*dphidzic
     &                             + d2drprime*dzicxid - d2drdxz
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                             + d2edphi2*dphidyid*dphidzic
     &                             + dzic*ddrdy - dyid*ddrdz
     &                             + dedphiprime*dphidyid*dphidzic
     &                             + d2edphi2prime*dphidyid*dphidzic
     &                             + d2drprime*dzicyid - d2drdyz
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                             + d2edphi2*dphidzid*dphidzic
     &                             + dzic*ddrdz - dzid*ddrdz
     &                             + dedphiprime*dphidzid*dphidzic
     &                             + d2edphi2prime*dphidzid*dphidzic
     &                             + d2drprime*dziczid - d2drdzz
               end if
            end if
         end if
      end do
      return
      end

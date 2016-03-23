c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2014 by Chao Lv & Jay William Ponder  ##
c     ##                  All Rights Reserved                 ##
c     ##########################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine eangtor2  --  atomwise angle-torsion Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "eangtor2" calculates the angle-torsion potential energy
c     second derivatives with respect to Cartesian coordinates
c
c
      subroutine eangtor2 (i)
      use sizes
      use angbnd
      use angtor
      use atoms
      use bound
      use group
      use hessn
      use math
      use torpot
      use tors
      implicit none
      integer i,j,k,iangtor
      integer ia,ib,ic,id
      real*8 dedphi,d2edphi2,fgrp
      real*8 rt2,ru2,rtru,rcb
      real*8 rba2,rcb2,rdc2
      real*8 dot,dt,d2dt
      real*8 xt,yt,zt
      real*8 xu,yu,zu
      real*8 xtu,ytu,ztu
      real*8 v1,v2,v3
      real*8 c1,c2,c3
      real*8 s1,s2,s3
      real*8 terma,termb
      real*8 termc,termd
      real*8 sine,cosine
      real*8 sine2,cosine2
      real*8 sine3,cosine3
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
      real*8 xab,yab,zab
      real*8 xbc,ybc,zbc
      real*8 xrab,yrab,zrab
      real*8 xrcb,yrcb,zrcb
      real*8 xabp,yabp,zabp
      real*8 xcbp,ycbp,zcbp
      real*8 xrbc,yrbc,zrbc
      real*8 xrdc,yrdc,zrdc
      real*8 xbcp,ybcp,zbcp
      real*8 xdcp,ydcp,zdcp
      real*8 phi1,phi2,phi3
      real*8 dphi1,dphi2,dphi3
      real*8 d2phi1,d2phi2,d2phi3
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
      real*8 domegadxia,domegadyia,domegadzia
      real*8 domegadxib,domegadyib,domegadzib
      real*8 domegadxic,domegadyic,domegadzic
      real*8 domegadxid,domegadyid,domegadzid
      real*8 doxiaxia,doyiayia,doziazia
      real*8 doxibxib,doyibyib,dozibzib
      real*8 doxicxic,doyicyic,doziczic
      real*8 doxidxid,doyidyid,dozidzid
      real*8 doxiayia,doxiazia,doyiazia
      real*8 doxibyib,doxibzib,doyibzib
      real*8 doxicyic,doxiczic,doyiczic
      real*8 doxidyid,doxidzid,doyidzid
      real*8 doxiaxic,doxiayic,doxiazic
      real*8 doyiaxic,doyiayic,doyiazic
      real*8 doziaxic,doziayic,doziazic
      real*8 doxibxic,doxibyic,doxibzic
      real*8 doyibxic,doyibyic,doyibzic
      real*8 dozibxic,dozibyic,dozibzic
      real*8 doxibxid,doxibyid,doxibzid
      real*8 doyibxid,doyibyid,doyibzid
      real*8 dozibxid,dozibyid,dozibzid
      real*8 doxicxid,doxicyid,doxiczid
      real*8 doyicxid,doyicyid,doyiczid
      real*8 dozicxid,dozicyid,doziczid
      real*8 doxibxia,doxibyia,doxibzia
      real*8 doyibxia,doyibyia,doyibzia
      real*8 dozibxia,dozibyia,dozibzia
      real*8 doxicxib,doxicyib,doxiczib
      real*8 doyicxib,doyicyib,doyiczib
      real*8 dozicxib,dozicyib,doziczib
      real*8 force,dedphiprime
      real*8 d2edphi2prime
      real*8 dtdphiprime
      logical proceed
c
c
c     compute the Hessian elements of the stretch-torsions
c
      do iangtor = 1, nangtor
         j = iat(1,iangtor)
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
            xab = -xba
            yab = -yba
            zab = -zba
            xbc = -xcb
            ybc = -ycb
            zbc = -zcb
            rba2 = xba*xba + yba*yba + zba*zba
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
            rdc2 = xdc*xdc + ydc*ydc + zdc*zdc
            if (use_polymer) then
               call image (xba,yba,zba)
               call image (xcb,ycb,zcb)
               call image (xdc,ydc,zdc)
               call image (xab,yab,zab)
               call image (xbc,ybc,zbc)
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
               d2phi1 = -(cosine*c1 + sine*s1)
               d2phi2 = -4.0d0 * (cosine2*c2 + sine2*s2)
               d2phi3 = -9.0d0 * (cosine3*c3 + sine3*s3)
c
c     set the angle-torsion parameters for the first angle
c
               v1 = kant(1,iangtor)
               v2 = kant(2,iangtor)
               v3 = kant(3,iangtor)
               k = iat(2,iangtor)
               dot = xab*xcb + yab*ycb + zab*zcb
               cosang = dot / sqrt(rba2*rcb2)
               angle = radian * acos(cosang)
               dt = angle - anat(k)
               force = ak(k)
               dedphi = atorunit * 2.0d0 * force 
     &                     * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2 = atorunit * dt * 2.0d0 * force
     &                       * (v1*d2phi1 + v2* d2phi2 + v3*d2phi3)
               d2dt = atorunit * radian * 2.0d0 * force
     &                   * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphiprime = atorunit * 2.0d0 * force
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2prime = atorunit * 2.0d0 * force
     &                            * (v1*phi1 + v2*phi2 + v3*phi3)
     &                            * (v1*d2phi1 + v2* d2phi2 + v3*d2phi3)
               dtdphiprime = atorunit * 2.0d0 * force
     &                          * (v1*phi1 + v2*phi2 + v3*phi3)
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
                  d2dt = d2dt * fgrp
                  dedphiprime = dedphiprime * fgrp
                  d2edphi2prime = d2edphi2prime * fgrp
                  dtdphiprime = dtdphiprime * fgrp
               end if
c
c     first and second derivative components for the first angle
c
               terma = -1.0d0 / (rba2*sqrt(rt2))
               termc = 1.0d0 / (rcb2*sqrt(rt2))
               domegadxia = terma * (zba*yt-yba*zt)
               domegadyia = terma * (xba*zt-zba*xt)
               domegadzia = terma * (yba*xt-xba*yt)
               domegadxic = termc * (ycb*zt-zcb*yt)
               domegadyic = termc * (zcb*xt-xcb*zt)
               domegadzic = termc * (xcb*yt-ycb*xt)
               domegadxib = -domegadxia - domegadxic
               domegadyib = -domegadyia - domegadyic
               domegadzib = -domegadzia - domegadzic
c
c     abbreviations used in defining chain rule terms
c
               xrab = 2.0d0 * xab / rba2
               yrab = 2.0d0 * yab / rba2
               zrab = 2.0d0 * zab / rba2
               xrcb = 2.0d0 * xcb / rcb2
               yrcb = 2.0d0 * ycb / rcb2
               zrcb = 2.0d0 * zcb / rcb2
               xabp = (yab*zt-zab*yt) / rt2
               yabp = (zab*xt-xab*zt) / rt2
               zabp = (xab*yt-yab*xt) / rt2
               xcbp = (ycb*zt-zcb*yt) / rt2
               ycbp = (zcb*xt-xcb*zt) / rt2
               zcbp = (xcb*yt-ycb*xt) / rt2
c
c     chain rule terms for second derivative components
c
               doxiaxia = terma*(xab*xcb-dot) + domegadxia*(xcbp-xrab)
               doxiayia = terma*(zt+yab*xcb) + domegadxia*(ycbp-yrab)
               doxiazia = terma*(zab*xcb-yt) + domegadxia*(zcbp-zrab)
               doyiayia = terma*(yab*ycb-dot) + domegadyia*(ycbp-yrab)
               doyiazia = terma*(xt+zab*ycb) + domegadyia*(zcbp-zrab)
               doziazia = terma*(zab*zcb-dot) + domegadzia*(zcbp-zrab)
               doxicxic = termc*(dot-xab*xcb) - domegadxic*(xabp+xrcb)
               doxicyic = termc*(zt-ycb*xab) - domegadxic*(yabp+yrcb)
               doxiczic = -termc*(yt+zcb*xab) - domegadxic*(zabp+zrcb)
               doyicyic = termc*(dot-yab*ycb) - domegadyic*(yabp+yrcb)
               doyiczic = termc*(xt-zcb*yab) - domegadyic*(zabp+zrcb)
               doziczic = termc*(dot-zab*zcb) - domegadzic*(zabp+zrcb)
               doxiaxic = terma*(yab*yab+zab*zab) - domegadxia*xabp
               doxiayic = -terma*xab*yab - domegadxia*yabp
               doxiazic = -terma*xab*zab - domegadxia*zabp
               doyiaxic = -terma*xab*yab - domegadyia*xabp
               doyiayic = terma*(xab*xab+zab*zab) - domegadyia*yabp
               doyiazic = -terma*yab*zab - domegadyia*zabp
               doziaxic = -terma*xab*zab - domegadzia*xabp
               doziayic = -terma*yab*zab - domegadzia*yabp
               doziazic = terma*(xab*xab+yab*yab) - domegadzia*zabp
c
c     get some second derivative chain rule terms by difference
c
               doxibxia = -doxiaxia - doxiaxic
               doxibyia = -doxiayia - doyiaxic
               doxibzia = -doxiazia - doziaxic
               doyibxia = -doxiayia - doxiayic
               doyibyia = -doyiayia - doyiayic
               doyibzia = -doyiazia - doziayic
               dozibxia = -doxiazia - doxiazic
               dozibyia = -doyiazia - doyiazic
               dozibzia = -doziazia - doziazic
               doxibxic = -doxicxic - doxiaxic
               doxibyic = -doxicyic - doxiayic
               doxibzic = -doxiczic - doxiazic
               doyibxic = -doxicyic - doyiaxic
               doyibyic = -doyicyic - doyiayic
               doyibzic = -doyiczic - doyiazic
               dozibxic = -doxiczic - doziaxic
               dozibyic = -doyiczic - doziayic
               dozibzic = -doziczic - doziazic
               doxibxib = -doxibxia - doxibxic
               doxibyib = -doxibyia - doxibyic
               doxibzib = -doxibzia - doxibzic
               doyibyib = -doyibyia - doyibyic
               doyibzib = -doyibzia - doyibzic
               dozibzib = -dozibzia - dozibzic
c
c     scale the first derivatives of the first angle
c
               domegadxia = domegadxia * radian
               domegadyia = domegadyia * radian
               domegadzia = domegadzia * radian
               domegadxic = domegadxic * radian
               domegadyic = domegadyic * radian
               domegadzic = domegadzic * radian
               domegadxib = domegadxib * radian
               domegadyib = domegadyib * radian
               domegadzib = domegadzib * radian
c
c     scale the second derivatives of the first angle
c
               doxiaxia = doxiaxia * d2dt
               doxiayia = doxiayia * d2dt
               doxiazia = doxiazia * d2dt
               doyiayia = doyiayia * d2dt
               doyiazia = doyiazia * d2dt
               doziazia = doziazia * d2dt
               doxicxic = doxicxic * d2dt
               doxicyic = doxicyic * d2dt
               doxiczic = doxiczic * d2dt
               doyicyic = doyicyic * d2dt
               doyiczic = doyiczic * d2dt
               doziczic = doziczic * d2dt
               doxiaxic = doxiaxic * d2dt
               doxiayic = doxiayic * d2dt
               doxiazic = doxiazic * d2dt
               doyiaxic = doyiaxic * d2dt
               doyiayic = doyiayic * d2dt
               doyiazic = doyiazic * d2dt
               doziaxic = doziaxic * d2dt
               doziayic = doziayic * d2dt
               doziazic = doziazic * d2dt
               doxibxia = doxibxia * d2dt
               doxibyia = doxibyia * d2dt
               doxibzia = doxibzia * d2dt
               doyibxia = doyibxia * d2dt
               doyibyia = doyibyia * d2dt
               doyibzia = doyibzia * d2dt
               dozibxia = dozibxia * d2dt
               dozibyia = dozibyia * d2dt
               dozibzia = dozibzia * d2dt
               doxibxic = doxibxic * d2dt
               doxibyic = doxibyic * d2dt
               doxibzic = doxibzic * d2dt
               doyibxic = doyibxic * d2dt
               doyibyic = doyibyic * d2dt
               doyibzic = doyibzic * d2dt
               dozibxic = dozibxic * d2dt
               dozibyic = dozibyic * d2dt
               dozibzic = dozibzic * d2dt
               doxibxib = doxibxib * d2dt
               doxibyib = doxibyib * d2dt
               doxibzib = doxibzib * d2dt
               doyibyib = doyibyib * d2dt
               doyibzib = doyibzib * d2dt
               dozibzib = dozibzib * d2dt
c
c     abbreviations for first derivative chain rule terms
c
               dphidxt = (yt*zcb-ycb*zt) / (rt2*rcb)
               dphidyt = (zt*xcb-zcb*xt) / (rt2*rcb)
               dphidzt = (xt*ycb-xcb*yt) / (rt2*rcb)
               dphidxu = -(yu*zcb-ycb*zu) / (ru2*rcb)
               dphidyu = -(zu*xcb-zcb*xu) / (ru2*rcb)
               dphidzu = -(xu*ycb-xcb*yu) / (ru2*rcb)
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
               dedphi = dedphi * dt
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
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia + doxiaxia
     &                             + d2edphi2*dphidxia*dphidxia
     &                             + 2.0d0*domegadxia*dxia
     &                             + dedphiprime*dphidxia*dphidxia
     &                             + d2edphi2prime*dphidxia*dphidxia
     &                             + dtdphiprime*dxiaxia               
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia + doxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             + domegadxia*dyia + domegadyia*dxia
     &                             + dedphiprime*dphidxia*dphidyia
     &                             + d2edphi2prime*dphidxia*dphidyia
     &                             + dtdphiprime*dxiayia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia + doxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             + domegadxia*dzia + domegadzia*dxia
     &                             + dedphiprime*dphidxia*dphidzia
     &                             + d2edphi2prime*dphidxia*dphidzia
     &                             + dtdphiprime*dxiazia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia + doxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             + domegadyia*dxia + domegadxia*dyia
     &                             + dedphiprime*dphidxia*dphidyia
     &                             + d2edphi2prime*dphidxia*dphidyia
     &                             + dtdphiprime*dxiayia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia + doyiayia
     &                             + d2edphi2*dphidyia*dphidyia
     &                             + 2.0d0*domegadyia*dyia
     &                             + dedphiprime*dphidyia*dphidyia
     &                             + d2edphi2prime*dphidyia*dphidyia
     &                             + dtdphiprime*dyiayia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia + doyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             + domegadyia*dzia + domegadzia*dyia
     &                             + dedphiprime*dphidyia*dphidzia
     &                             + d2edphi2prime*dphidyia*dphidzia
     &                             + dtdphiprime*dyiazia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia + doxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             + domegadxia*dzia + domegadzia*dxia
     &                             + dedphiprime*dphidxia*dphidzia
     &                             + d2edphi2prime*dphidxia*dphidzia
     &                             + dtdphiprime*dxiazia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia + doyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             + domegadyia*dzia + domegadzia*dyia
     &                             + dedphiprime*dphidyia*dphidzia
     &                             + d2edphi2prime*dphidyia*dphidzia
     &                             + dtdphiprime*dyiazia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia + doziazia
     &                             + d2edphi2*dphidzia*dphidzia
     &                             + 2.0d0*domegadzia*dzia
     &                             + dedphiprime*dphidzia*dphidzia
     &                             + d2edphi2prime*dphidzia*dphidzia
     &                             + dtdphiprime*dziazia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib + doxibxia
     &                             + d2edphi2*dphidxia*dphidxib
     &                             + domegadxia*dxib + domegadxib*dxia
     &                             + dedphiprime*dphidxia*dphidxib
     &                             + d2edphi2prime*dphidxia*dphidxib
     &                             + dtdphiprime*dxiaxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib + doxibyia
     &                             + d2edphi2*dphidyia*dphidxib
     &                             + domegadyia*dxib + domegadxib*dyia
     &                             + dedphiprime*dphidyia*dphidxib
     &                             + d2edphi2prime*dphidyia*dphidxib
     &                             + dtdphiprime*dyiaxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib + doxibzia
     &                             + d2edphi2*dphidzia*dphidxib
     &                             + domegadzia*dxib + domegadxib*dzia
     &                             + dedphiprime*dphidzia*dphidxib
     &                             + d2edphi2prime*dphidzia*dphidxib
     &                             + dtdphiprime*dziaxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib + doyibxia
     &                             + d2edphi2*dphidxia*dphidyib
     &                             + domegadxia*dyib + domegadyib*dxia
     &                             + dedphiprime*dphidxia*dphidyib
     &                             + d2edphi2prime*dphidxia*dphidyib
     &                             + dtdphiprime*dxiayib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib + doyibyia
     &                             + d2edphi2*dphidyia*dphidyib
     &                             + domegadyia*dyib + domegadyib*dyia
     &                             + dedphiprime*dphidyia*dphidyib
     &                             + d2edphi2prime*dphidyia*dphidyib
     &                             + dtdphiprime*dyiayib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib + doyibzia
     &                             + d2edphi2*dphidzia*dphidyib
     &                             + domegadzia*dyib + domegadyib*dzia
     &                             + dedphiprime*dphidzia*dphidyib
     &                             + d2edphi2prime*dphidzia*dphidyib
     &                             + dtdphiprime*dziayib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib + dozibxia
     &                             + d2edphi2*dphidxia*dphidzib
     &                             + domegadxia*dzib + domegadzib*dxia
     &                             + dedphiprime*dphidxia*dphidzib
     &                             + d2edphi2prime*dphidxia*dphidzib
     &                             + dtdphiprime*dxiazib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib + dozibyia
     &                             + d2edphi2*dphidyia*dphidzib
     &                             + domegadyia*dzib + domegadzib*dyia
     &                             + dedphiprime*dphidyia*dphidzib
     &                             + d2edphi2prime*dphidyia*dphidzib
     &                             + dtdphiprime*dyiazib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib + dozibzia
     &                             + d2edphi2*dphidzia*dphidzib
     &                             + domegadzia*dzib + domegadzib*dzia
     &                             + dedphiprime*dphidzia*dphidzib
     &                             + d2edphi2prime*dphidzia*dphidzib
     &                             + dtdphiprime*dziazib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic + doxiaxic
     &                             + d2edphi2*dphidxia*dphidxic
     &                             + domegadxia*dxic + domegadxic*dxia
     &                             + dedphiprime*dphidxia*dphidxic
     &                             + d2edphi2prime*dphidxia*dphidxic
     &                             + dtdphiprime*dxiaxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic + doyiaxic
     &                             + d2edphi2*dphidyia*dphidxic
     &                             + domegadyia*dxic + domegadxic*dyia
     &                             + dedphiprime*dphidyia*dphidxic
     &                             + d2edphi2prime*dphidyia*dphidxic
     &                             + dtdphiprime*dyiaxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic + doziaxic
     &                             + d2edphi2*dphidzia*dphidxic
     &                             + domegadzia*dxic + domegadxic*dzia
     &                             + dedphiprime*dphidzia*dphidxic
     &                             + d2edphi2prime*dphidzia*dphidxic
     &                             + dtdphiprime*dziaxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic + doxiayic
     &                             + d2edphi2*dphidxia*dphidyic
     &                             + domegadxia*dyic + domegadyic*dxia
     &                             + dedphiprime*dphidxia*dphidyic
     &                             + d2edphi2prime*dphidxia*dphidyic
     &                             + dtdphiprime*dxiayic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic + doyiayic
     &                             + d2edphi2*dphidyia*dphidyic
     &                             + domegadyia*dyic + domegadyic*dyia
     &                             + dedphiprime*dphidyia*dphidyic
     &                             + d2edphi2prime*dphidyia*dphidyic
     &                             + dtdphiprime*dyiayic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic + doziayic
     &                             + d2edphi2*dphidzia*dphidyic
     &                             + domegadzia*dyic + domegadyic*dzia
     &                             + dedphiprime*dphidzia*dphidyic
     &                             + d2edphi2prime*dphidzia*dphidyic
     &                             + dtdphiprime*dziayic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                             + d2edphi2*dphidxia*dphidzic
     &                             + domegadxia*dzic + domegadzic*dxia
     &                             + doxiazic
     &                             + dedphiprime*dphidxia*dphidzic
     &                             + d2edphi2prime*dphidxia*dphidzic
     &                             + dtdphiprime*dxiazic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic + doyiazic
     &                             + d2edphi2*dphidyia*dphidzic
     &                             + domegadyia*dzic + domegadzic*dyia
     &                             + dedphiprime*dphidyia*dphidzic
     &                             + d2edphi2prime*dphidyia*dphidzic
     &                             + dtdphiprime*dyiazic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic + doziazic
     &                             + d2edphi2*dphidzia*dphidzic
     &                             + domegadzia*dzic + domegadzic*dzia
     &                             + dedphiprime*dphidzia*dphidzic
     &                             + d2edphi2prime*dphidzia*dphidzic
     &                             + dtdphiprime*dziazic
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxia*dphidxid
     &                             + domegadxia*dxid
     &                             + dedphiprime*dphidxia*dphidxid
     &                             + d2edphi2prime*dphidxia*dphidxid
     &                             + dtdphiprime*dxiaxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                             + d2edphi2*dphidyia*dphidxid
     &                             + domegadyia*dxid
     &                             + dedphiprime*dphidyia*dphidxid
     &                             + d2edphi2prime*dphidyia*dphidxid
     &                             + dtdphiprime*dyiaxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                             + d2edphi2*dphidzia*dphidxid
     &                             + domegadzia*dxid
     &                             + dedphiprime*dphidzia*dphidxid
     &                             + d2edphi2prime*dphidzia*dphidxid
     &                             + dtdphiprime*dziaxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                             + d2edphi2*dphidxia*dphidyid
     &                             + domegadxia*dyid
     &                             + dedphiprime*dphidxia*dphidyid
     &                             + d2edphi2prime*dphidxia*dphidyid
     &                             + dtdphiprime*dxiayid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                             + d2edphi2*dphidyia*dphidyid
     &                             + domegadyia*dyid
     &                             + dedphiprime*dphidyia*dphidyid
     &                             + d2edphi2prime*dphidyia*dphidyid
     &                             + dtdphiprime*dyiayid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                             + d2edphi2*dphidzia*dphidyid
     &                             + domegadzia*dyid
     &                             + dedphiprime*dphidzia*dphidyid
     &                             + d2edphi2prime*dphidzia*dphidyid
     &                             + dtdphiprime*dziayid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                             + d2edphi2*dphidxia*dphidzid
     &                             + domegadxia*dzid
     &                             + dedphiprime*dphidxia*dphidzid
     &                             + d2edphi2prime*dphidxia*dphidzid
     &                             + dtdphiprime*dxiazid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                             + d2edphi2*dphidyia*dphidzid
     &                             + domegadyia*dzid
     &                             + dedphiprime*dphidyia*dphidzid
     &                             + d2edphi2prime*dphidyia*dphidzid
     &                             + dtdphiprime*dyiazid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                             + d2edphi2*dphidzia*dphidzid
     &                             + domegadzia*dzid
     &                             + dedphiprime*dphidzia*dphidzid
     &                             + d2edphi2prime*dphidzia*dphidzid
     &                             + dtdphiprime*dziazid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib + doxibxib
     &                             + d2edphi2*dphidxib*dphidxib
     &                             + 2.0d0*domegadxib*dxib
     &                             + dedphiprime*dphidxib*dphidxib
     &                             + d2edphi2prime*dphidxib*dphidxib
     &                             + dtdphiprime*dxibxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib + doxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + domegadxib*dyib + domegadyib*dxib
     &                             + dedphiprime*dphidxib*dphidyib
     &                             + d2edphi2prime*dphidxib*dphidyib
     &                             + dtdphiprime*dxibyib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib + doxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + domegadxib*dzib + domegadzib*dxib
     &                             + dedphiprime*dphidxib*dphidzib
     &                             + d2edphi2prime*dphidxib*dphidzib
     &                             + dtdphiprime*dxibzib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib + doxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + domegadxib*dyib + domegadyib*dxib
     &                             + dedphiprime*dphidxib*dphidyib
     &                             + d2edphi2prime*dphidxib*dphidyib
     &                             + dtdphiprime*dxibyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib + doyibyib
     &                             + d2edphi2*dphidyib*dphidyib
     &                             + 2.0d0*domegadyib*dyib
     &                             + dedphiprime*dphidyib*dphidyib
     &                             + d2edphi2prime*dphidyib*dphidyib
     &                             + dtdphiprime*dyibyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib + doyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + domegadyib*dzib + domegadzib*dyib
     &                             + dedphiprime*dphidyib*dphidzib
     &                             + d2edphi2prime*dphidyib*dphidzib
     &                             + dtdphiprime*dyibzib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib + doxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + domegadxib*dzib + domegadzib*dxib
     &                             + dedphiprime*dphidxib*dphidzib
     &                             + d2edphi2prime*dphidxib*dphidzib
     &                             + dtdphiprime*dxibzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib + doyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + domegadyib*dzib + domegadzib*dyib
     &                             + dedphiprime*dphidyib*dphidzib
     &                             + d2edphi2prime*dphidyib*dphidzib
     &                             + dtdphiprime*dyibzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib + dozibzib
     &                             + d2edphi2*dphidzib*dphidzib
     &                             + 2.0d0*domegadzib*dzib
     &                             + dedphiprime*dphidzib*dphidzib
     &                             + d2edphi2prime*dphidzib*dphidzib
     &                             + dtdphiprime*dzibzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib + doxibxia
     &                             + d2edphi2*dphidxib*dphidxia
     &                             + domegadxia*dxib + domegadxib*dxia
     &                             + dedphiprime*dphidxib*dphidxia
     &                             + d2edphi2prime*dphidxib*dphidxia
     &                             + dtdphiprime*dxiaxib
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib + doyibxia
     &                             + d2edphi2*dphidyib*dphidxia
     &                             + domegadxia*dyib + domegadyib*dxia
     &                             + dedphiprime*dphidyib*dphidxia
     &                             + d2edphi2prime*dphidyib*dphidxia
     &                             + dtdphiprime*dxiayib
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib + dozibxia
     &                             + d2edphi2*dphidzib*dphidxia
     &                             + domegadxia*dzib + domegadzib*dxia
     &                             + dedphiprime*dphidzib*dphidxia
     &                             + d2edphi2prime*dphidzib*dphidxia
     &                             + dtdphiprime*dxiazib
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib + doxibyia
     &                             + d2edphi2*dphidxib*dphidyia
     &                             + domegadyia*dxib + domegadxib*dyia
     &                             + dedphiprime*dphidxib*dphidyia
     &                             + d2edphi2prime*dphidxib*dphidyia
     &                             + dtdphiprime*dyiaxib
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib + doyibyia
     &                             + d2edphi2*dphidyib*dphidyia
     &                             + domegadyia*dyib + domegadyib*dyia
     &                             + dedphiprime*dphidyib*dphidyia
     &                             + d2edphi2prime*dphidyib*dphidyia
     &                             + dtdphiprime*dyiayib
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib + dozibyia
     &                             + d2edphi2*dphidzib*dphidyia
     &                             + domegadyia*dzib + domegadzib*dyia
     &                             + dedphiprime*dphidzib*dphidyia
     &                             + d2edphi2prime*dphidzib*dphidyia
     &                             + dtdphiprime*dyiazib
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib + doxibzia
     &                             + d2edphi2*dphidxib*dphidzia
     &                             + domegadzia*dxib + domegadxib*dzia
     &                             + dedphiprime*dphidxib*dphidzia
     &                             + d2edphi2prime*dphidxib*dphidzia
     &                             + dtdphiprime*dziaxib
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib + doyibzia
     &                             + d2edphi2*dphidyib*dphidzia
     &                             + domegadzia*dyib + domegadyib*dzia
     &                             + dedphiprime*dphidyib*dphidzia
     &                             + d2edphi2prime*dphidyib*dphidzia
     &                             + dtdphiprime*dziayib
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib + dozibzia
     &                             + d2edphi2*dphidzib*dphidzia
     &                             + domegadzia*dzib + domegadzib*dzia
     &                             + dedphiprime*dphidzib*dphidzia
     &                             + d2edphi2prime*dphidzib*dphidzia
     &                             + dtdphiprime*dziazib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic + doxibxic
     &                             + d2edphi2*dphidxib*dphidxic
     &                             + domegadxib*dxic + domegadxic*dxib
     &                             + dedphiprime*dphidxib*dphidxic
     &                             + d2edphi2prime*dphidxib*dphidxic
     &                             + dtdphiprime*dxibxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic + doyibxic
     &                             + d2edphi2*dphidyib*dphidxic
     &                             + domegadyib*dxic + domegadxic*dyib
     &                             + dedphiprime*dphidyib*dphidxic
     &                             + d2edphi2prime*dphidyib*dphidxic
     &                             + dtdphiprime*dyibxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic + dozibxic
     &                             + d2edphi2*dphidzib*dphidxic
     &                             + domegadzib*dxic + domegadxic*dzib
     &                             + dedphiprime*dphidzib*dphidxic
     &                             + d2edphi2prime*dphidzib*dphidxic
     &                             + dtdphiprime*dzibxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic + doxibyic
     &                             + d2edphi2*dphidxib*dphidyic
     &                             + domegadxib*dyic + domegadyic*dxib
     &                             + dedphiprime*dphidxib*dphidyic
     &                             + d2edphi2prime*dphidxib*dphidyic
     &                             + dtdphiprime*dxibyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic + doyibyic
     &                             + d2edphi2*dphidyib*dphidyic
     &                             + domegadyib*dyic + domegadyic*dyib
     &                             + dedphiprime*dphidyib*dphidyic
     &                             + d2edphi2prime*dphidyib*dphidyic
     &                             + dtdphiprime*dyibyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic + dozibyic
     &                             + d2edphi2*dphidzib*dphidyic
     &                             + domegadzib*dyic + domegadyic*dzib
     &                             + dedphiprime*dphidzib*dphidyic
     &                             + d2edphi2prime*dphidzib*dphidyic
     &                             + dtdphiprime*dzibyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic + doxibzic
     &                             + d2edphi2*dphidxib*dphidzic
     &                             + domegadxib*dzic + domegadzic*dxib
     &                             + dedphiprime*dphidxib*dphidzic
     &                             + d2edphi2prime*dphidxib*dphidzic
     &                             + dtdphiprime*dxibzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic + doyibzic
     &                             + d2edphi2*dphidyib*dphidzic
     &                             + domegadyib*dzic + domegadzic*dyib
     &                             + dedphiprime*dphidyib*dphidzic
     &                             + d2edphi2prime*dphidyib*dphidzic
     &                             + dtdphiprime*dyibzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic + dozibzic
     &                             + d2edphi2*dphidzib*dphidzic
     &                             + domegadzib*dzic + domegadzic*dzib
     &                             + dedphiprime*dphidzib*dphidzic
     &                             + d2edphi2prime*dphidzib*dphidzic
     &                             + dtdphiprime*dzibzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid
     &                             + d2edphi2*dphidxib*dphidxid
     &                             + domegadxib*dxid
     &                             + dedphiprime*dphidxib*dphidxid
     &                             + d2edphi2prime*dphidxib*dphidxid
     &                             + dtdphiprime*dxibxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid
     &                             + d2edphi2*dphidyib*dphidxid
     &                             + domegadyib*dxid
     &                             + dedphiprime*dphidyib*dphidxid
     &                             + d2edphi2prime*dphidyib*dphidxid
     &                             + dtdphiprime*dyibxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid
     &                             + d2edphi2*dphidzib*dphidxid
     &                             + domegadzib*dxid
     &                             + dedphiprime*dphidzib*dphidxid
     &                             + d2edphi2prime*dphidzib*dphidxid
     &                             + dtdphiprime*dzibxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid
     &                             + d2edphi2*dphidxib*dphidyid
     &                             + domegadxib*dyid
     &                             + dedphiprime*dphidxib*dphidyid
     &                             + d2edphi2prime*dphidxib*dphidyid
     &                             + dtdphiprime*dxibyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid
     &                             + d2edphi2*dphidyib*dphidyid
     &                             + domegadyib*dyid
     &                             + dedphiprime*dphidyib*dphidyid
     &                             + d2edphi2prime*dphidyib*dphidyid
     &                             + dtdphiprime*dyibyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid
     &                             + d2edphi2*dphidzib*dphidyid
     &                             + domegadzib*dyid
     &                             + dedphiprime*dphidzib*dphidyid
     &                             + d2edphi2prime*dphidzib*dphidyid
     &                             + dtdphiprime*dzibyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid
     &                             + d2edphi2*dphidxib*dphidzid
     &                             + domegadxib*dzid
     &                             + dedphiprime*dphidxib*dphidzid
     &                             + d2edphi2prime*dphidxib*dphidzid
     &                             + dtdphiprime*dxibzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid
     &                             + d2edphi2*dphidyib*dphidzid
     &                             + domegadyib*dzid
     &                             + dedphiprime*dphidyib*dphidzid
     &                             + d2edphi2prime*dphidyib*dphidzid
     &                             + dtdphiprime*dyibzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid
     &                             + d2edphi2*dphidzib*dphidzid
     &                             + domegadzib*dzid
     &                             + dedphiprime*dphidzib*dphidzid
     &                             + d2edphi2prime*dphidzib*dphidzid
     &                             + dtdphiprime*dzibzid
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic + doxicxic
     &                             + d2edphi2*dphidxic*dphidxic
     &                             + 2.0d0*domegadxic*dxic
     &                             + dedphiprime*dphidxic*dphidxic
     &                             + d2edphi2prime*dphidxic*dphidxic
     &                             + dtdphiprime*dxicxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic + doxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + domegadxic*dyic + domegadyic*dxic
     &                             + dedphiprime*dphidxic*dphidyic
     &                             + d2edphi2prime*dphidxic*dphidyic
     &                             + dtdphiprime*dxicyic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic + doxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + domegadxic*dzic + domegadzic*dxic
     &                             + dedphiprime*dphidxic*dphidzic
     &                             + d2edphi2prime*dphidxic*dphidzic
     &                             + dtdphiprime*dxiczic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic + doxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + domegadxic*dyic + domegadyic*dxic
     &                             + dedphiprime*dphidxic*dphidyic
     &                             + d2edphi2prime*dphidxic*dphidyic
     &                             + dtdphiprime*dxicyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic + doyicyic
     &                             + d2edphi2*dphidyic*dphidyic
     &                             + 2.0d0*domegadyic*dyic
     &                             + dedphiprime*dphidyic*dphidyic
     &                             + d2edphi2prime*dphidyic*dphidyic
     &                             + dtdphiprime*dyicyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic + doyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + domegadyic*dzic + domegadzic*dyic
     &                             + dedphiprime*dphidyic*dphidzic
     &                             + d2edphi2prime*dphidyic*dphidzic
     &                             + dtdphiprime*dyiczic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic + doxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + domegadxic*dzic + domegadzic*dxic
     &                             + dedphiprime*dphidxic*dphidzic
     &                             + d2edphi2prime*dphidxic*dphidzic
     &                             + dtdphiprime*dxiczic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic + doyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + domegadyic*dzic + domegadzic*dyic
     &                             + dedphiprime*dphidyic*dphidzic
     &                             + d2edphi2prime*dphidyic*dphidzic
     &                             + dtdphiprime*dyiczic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic + doziczic
     &                             + d2edphi2*dphidzic*dphidzic
     &                             + 2.0d0*domegadzic*dzic
     &                             + dedphiprime*dphidzic*dphidzic
     &                             + d2edphi2prime*dphidzic*dphidzic
     &                             + dtdphiprime*dziczic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic + doxiaxic
     &                             + d2edphi2*dphidxic*dphidxia
     &                             + domegadxia*dxic + domegadxic*dxia
     &                             + dedphiprime*dphidxic*dphidxia
     &                             + d2edphi2prime*dphidxic*dphidxia
     &                             + dtdphiprime*dxiaxic
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic + doxiayic
     &                             + d2edphi2*dphidyic*dphidxia
     &                             + domegadxia*dyic + domegadyic*dxia
     &                             + dedphiprime*dphidyic*dphidxia
     &                             + d2edphi2prime*dphidyic*dphidxia
     &                             + dtdphiprime*dxiayic
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic + doxiazic
     &                             + d2edphi2*dphidzic*dphidxia
     &                             + domegadxia*dzic + domegadzic*dxia
     &                             + dedphiprime*dphidzic*dphidxia
     &                             + d2edphi2prime*dphidzic*dphidxia
     &                             + dtdphiprime*dxiazic
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic + doyiaxic
     &                             + d2edphi2*dphidxic*dphidyia
     &                             + domegadyia*dxic + domegadxic*dyia
     &                             + dedphiprime*dphidxic*dphidyia
     &                             + d2edphi2prime*dphidxic*dphidyia
     &                             + dtdphiprime*dyiaxic
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic + doyiayic
     &                             + d2edphi2*dphidyic*dphidyia
     &                             + domegadyia*dyic + domegadyic*dyia
     &                             + dedphiprime*dphidyic*dphidyia
     &                             + d2edphi2prime*dphidyic*dphidyia
     &                             + dtdphiprime*dyiayic
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic + doyiazic
     &                             + d2edphi2*dphidzic*dphidyia
     &                             + domegadyia*dzic + domegadzic*dyia
     &                             + dedphiprime*dphidzic*dphidyia
     &                             + d2edphi2prime*dphidzic*dphidyia
     &                             + dtdphiprime*dyiazic
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic + doziaxic
     &                             + d2edphi2*dphidxic*dphidzia
     &                             + domegadzia*dxic + domegadxic*dzia
     &                             + dedphiprime*dphidxic*dphidzia
     &                             + d2edphi2prime*dphidxic*dphidzia
     &                             + dtdphiprime*dziaxic
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic + doziayic
     &                             + d2edphi2*dphidyic*dphidzia
     &                             + domegadzia*dyic + domegadyic*dzia
     &                             + dedphiprime*dphidyic*dphidzia
     &                             + d2edphi2prime*dphidyic*dphidzia
     &                             + dtdphiprime*dziayic
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic + doziazic
     &                             + d2edphi2*dphidzic*dphidzia
     &                             + domegadzia*dzic + domegadzic*dzia
     &                             + dedphiprime*dphidzic*dphidzia
     &                             + d2edphi2prime*dphidzic*dphidzia
     &                             + dtdphiprime*dziazic
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic + doxibxic
     &                             + d2edphi2*dphidxic*dphidxib
     &                             + domegadxib*dxic + domegadxic*dxib
     &                             + dedphiprime*dphidxic*dphidxib
     &                             + d2edphi2prime*dphidxic*dphidxib
     &                             + dtdphiprime*dxibxic
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic + doxibyic
     &                             + d2edphi2*dphidyic*dphidxib
     &                             + domegadxib*dyic + domegadyic*dxib
     &                             + dedphiprime*dphidyic*dphidxib
     &                             + d2edphi2prime*dphidyic*dphidxib
     &                             + dtdphiprime*dxibyic
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic + doxibzic
     &                             + d2edphi2*dphidzic*dphidxib
     &                             + domegadxib*dzic + domegadzic*dxib
     &                             + dedphiprime*dphidzic*dphidxib
     &                             + d2edphi2prime*dphidzic*dphidxib
     &                             + dtdphiprime*dxibzic
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic + doyibxic
     &                             + d2edphi2*dphidxic*dphidyib
     &                             + domegadyib*dxic + domegadxic*dyib
     &                             + dedphiprime*dphidxic*dphidyib
     &                             + d2edphi2prime*dphidxic*dphidyib
     &                             + dtdphiprime*dyibxic
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic + doyibyic
     &                             + d2edphi2*dphidyic*dphidyib
     &                             + domegadyib*dyic + domegadyic*dyib
     &                             + dedphiprime*dphidyic*dphidyib
     &                             + d2edphi2prime*dphidyic*dphidyib
     &                             + dtdphiprime*dyibyic
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic + doyibzic
     &                             + d2edphi2*dphidzic*dphidyib
     &                             + domegadyib*dzic + domegadzic*dyib
     &                             + dedphiprime*dphidzic*dphidyib
     &                             + d2edphi2prime*dphidzic*dphidyib
     &                             + dtdphiprime*dyibzic
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic + dozibxic
     &                             + d2edphi2*dphidxic*dphidzib
     &                             + domegadzib*dxic + domegadxic*dzib
     &                             + dedphiprime*dphidxic*dphidzib
     &                             + d2edphi2prime*dphidxic*dphidzib
     &                             + dtdphiprime*dzibxic
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic + dozibyic
     &                             + d2edphi2*dphidyic*dphidzib
     &                             + domegadzib*dyic + domegadyic*dzib
     &                             + dedphiprime*dphidyic*dphidzib
     &                             + d2edphi2prime*dphidyic*dphidzib
     &                             + dtdphiprime*dzibyic
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic + dozibzic
     &                             + d2edphi2*dphidzic*dphidzib
     &                             + domegadzib*dzic + domegadzic*dzib
     &                             + dedphiprime*dphidzic*dphidzib
     &                             + d2edphi2prime*dphidzic*dphidzib
     &                             + dtdphiprime*dzibzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid
     &                             + d2edphi2*dphidxic*dphidxid
     &                             + domegadxic*dxid
     &                             + dedphiprime*dphidxic*dphidxid
     &                             + d2edphi2prime*dphidxic*dphidxid
     &                             + dtdphiprime*dxicxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid
     &                             + d2edphi2*dphidyic*dphidxid
     &                             + domegadyic*dxid
     &                             + dedphiprime*dphidyic*dphidxid
     &                             + d2edphi2prime*dphidyic*dphidxid
     &                             + dtdphiprime*dyicxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid
     &                             + d2edphi2*dphidzic*dphidxid
     &                             + domegadzic*dxid
     &                             + dedphiprime*dphidzic*dphidxid
     &                             + d2edphi2prime*dphidzic*dphidxid
     &                             + dtdphiprime*dzicxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid
     &                             + d2edphi2*dphidxic*dphidyid
     &                             + domegadxic*dyid
     &                             + dedphiprime*dphidxic*dphidyid
     &                             + d2edphi2prime*dphidxic*dphidyid
     &                             + dtdphiprime*dxicyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid
     &                             + d2edphi2*dphidyic*dphidyid
     &                             + domegadyic*dyid
     &                             + dedphiprime*dphidyic*dphidyid
     &                             + d2edphi2prime*dphidyic*dphidyid
     &                             + dtdphiprime*dyicyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid
     &                             + d2edphi2*dphidzic*dphidyid
     &                             + domegadzic*dyid
     &                             + dedphiprime*dphidzic*dphidyid
     &                             + d2edphi2prime*dphidzic*dphidyid
     &                             + dtdphiprime*dzicyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid
     &                             + d2edphi2*dphidxic*dphidzid
     &                             + domegadxic*dzid
     &                             + dedphiprime*dphidxic*dphidzid
     &                             + d2edphi2prime*dphidxic*dphidzid
     &                             + dtdphiprime*dxiczid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid
     &                             + d2edphi2*dphidyic*dphidzid
     &                             + domegadyic*dzid
     &                             + dedphiprime*dphidyic*dphidzid
     &                             + d2edphi2prime*dphidyic*dphidzid
     &                             + dtdphiprime*dyiczid
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid
     &                             + d2edphi2*dphidzic*dphidzid
     &                             + domegadzic*dzid
     &                             + dedphiprime*dphidzic*dphidzid
     &                             + d2edphi2prime*dphidzic*dphidzid
     &                             + dtdphiprime*dziczid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid
     &                             + d2edphi2*dphidxid*dphidxid
     &                             + dedphiprime*dphidxid*dphidxid
     &                             + d2edphi2prime*dphidxid*dphidxid
     &                             + dtdphiprime*dxidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + dedphiprime*dphidxid*dphidyid
     &                             + d2edphi2prime*dphidxid*dphidyid
     &                             + dtdphiprime*dxidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + dedphiprime*dphidxid*dphidzid
     &                             + d2edphi2prime*dphidxid*dphidzid
     &                             + dtdphiprime*dxidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + dedphiprime*dphidxid*dphidyid
     &                             + d2edphi2prime*dphidxid*dphidyid
     &                             + dtdphiprime*dxidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid
     &                             + d2edphi2*dphidyid*dphidyid
     &                             + dedphiprime*dphidyid*dphidyid
     &                             + d2edphi2prime*dphidyid*dphidyid
     &                             + dtdphiprime*dyidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + dedphiprime*dphidyid*dphidzid
     &                             + d2edphi2prime*dphidyid*dphidzid
     &                             + dtdphiprime*dyidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + dedphiprime*dphidxid*dphidzid
     &                             + d2edphi2prime*dphidxid*dphidzid
     &                             + dtdphiprime*dxidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + dedphiprime*dphidyid*dphidzid
     &                             + d2edphi2prime*dphidyid*dphidzid
     &                             + dtdphiprime*dyidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid
     &                             + d2edphi2*dphidzid*dphidzid
     &                             + dedphiprime*dphidzid*dphidzid
     &                             + d2edphi2prime*dphidzid*dphidzid
     &                             + dtdphiprime*dzidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxid*dphidxia
     &                             + domegadxia*dxid
     &                             + dedphiprime*dphidxid*dphidxia
     &                             + d2edphi2prime*dphidxid*dphidxia
     &                             + dtdphiprime*dxiaxid
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                             + d2edphi2*dphidyid*dphidxia
     &                             + domegadxia*dyid
     &                             + dedphiprime*dphidyid*dphidxia
     &                             + d2edphi2prime*dphidyid*dphidxia
     &                             + dtdphiprime*dxiayid
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                             + d2edphi2*dphidzid*dphidxia
     &                             + domegadxia*dzid
     &                             + dedphiprime*dphidzid*dphidxia
     &                             + d2edphi2prime*dphidzid*dphidxia
     &                             + dtdphiprime*dxiazid
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                             + d2edphi2*dphidxid*dphidyia
     &                             + domegadyia*dxid
     &                             + dedphiprime*dphidxid*dphidyia
     &                             + d2edphi2prime*dphidxid*dphidyia
     &                             + dtdphiprime*dyiaxid
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                             + d2edphi2*dphidyid*dphidyia
     &                             + domegadyia*dyid
     &                             + dedphiprime*dphidyid*dphidyia
     &                             + d2edphi2prime*dphidyid*dphidyia
     &                             + dtdphiprime*dyiayid
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                             + d2edphi2*dphidzid*dphidyia
     &                             + domegadyia*dzid
     &                             + dedphiprime*dphidzid*dphidyia
     &                             + d2edphi2prime*dphidzid*dphidyia
     &                             + dtdphiprime*dyiazid
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                             + d2edphi2*dphidxid*dphidzia
     &                             + domegadzia*dxid
     &                             + dedphiprime*dphidxid*dphidzia
     &                             + d2edphi2prime*dphidxid*dphidzia
     &                             + dtdphiprime*dziaxid
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                             + d2edphi2*dphidyid*dphidzia
     &                             + domegadzia*dyid
     &                             + dedphiprime*dphidyid*dphidzia
     &                             + d2edphi2prime*dphidyid*dphidzia
     &                             + dtdphiprime*dziayid
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                             + d2edphi2*dphidzid*dphidzia
     &                             + domegadzia*dzid
     &                             + dedphiprime*dphidzid*dphidzia
     &                             + d2edphi2prime*dphidzid*dphidzia
     &                             + dtdphiprime*dziazid
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid
     &                             + d2edphi2*dphidxid*dphidxib
     &                             + domegadxib*dxid
     &                             + dedphiprime*dphidxid*dphidxib
     &                             + d2edphi2prime*dphidxid*dphidxib
     &                             + dtdphiprime*dxibxid
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid
     &                             + d2edphi2*dphidyid*dphidxib
     &                             + domegadxib*dyid
     &                             + dedphiprime*dphidyid*dphidxib
     &                             + d2edphi2prime*dphidyid*dphidxib
     &                             + dtdphiprime*dxibyid
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid
     &                             + d2edphi2*dphidzid*dphidxib
     &                             + domegadxib*dzid
     &                             + dedphiprime*dphidzid*dphidxib
     &                             + d2edphi2prime*dphidzid*dphidxib
     &                             + dtdphiprime*dxibzid
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid
     &                             + d2edphi2*dphidxid*dphidyib
     &                             + domegadyib*dxid
     &                             + dedphiprime*dphidxid*dphidyib
     &                             + d2edphi2prime*dphidxid*dphidyib
     &                             + dtdphiprime*dyibxid
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid
     &                             + d2edphi2*dphidyid*dphidyib
     &                             + domegadyib*dyid
     &                             + dedphiprime*dphidyid*dphidyib
     &                             + d2edphi2prime*dphidyid*dphidyib
     &                             + dtdphiprime*dyibyid
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid
     &                             + d2edphi2*dphidzid*dphidyib
     &                             + domegadyib*dzid
     &                             + dedphiprime*dphidzid*dphidyib
     &                             + d2edphi2prime*dphidzid*dphidyib
     &                             + dtdphiprime*dyibzid
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid
     &                             + d2edphi2*dphidxid*dphidzib
     &                             + domegadzib*dxid
     &                             + dedphiprime*dphidxid*dphidzib
     &                             + d2edphi2prime*dphidxid*dphidzib
     &                             + dtdphiprime*dzibxid
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid
     &                             + d2edphi2*dphidyid*dphidzib
     &                             + domegadzib*dyid
     &                             + dedphiprime*dphidyid*dphidzib
     &                             + d2edphi2prime*dphidyid*dphidzib
     &                             + dtdphiprime*dzibyid
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid
     &                             + d2edphi2*dphidzid*dphidzib
     &                             + domegadzib*dzid
     &                             + dedphiprime*dphidzid*dphidzib
     &                             + d2edphi2prime*dphidzid*dphidzib
     &                             + dtdphiprime*dzibzid
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid
     &                             + d2edphi2*dphidxid*dphidxic
     &                             + domegadxic*dxid
     &                             + dedphiprime*dphidxid*dphidxic
     &                             + d2edphi2prime*dphidxid*dphidxic
     &                             + dtdphiprime*dxicxid
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid
     &                             + d2edphi2*dphidyid*dphidxic
     &                             + domegadxic*dyid
     &                             + dedphiprime*dphidyid*dphidxic
     &                             + d2edphi2prime*dphidyid*dphidxic
     &                             + dtdphiprime*dxicyid
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid
     &                             + d2edphi2*dphidzid*dphidxic
     &                             + domegadxic*dzid
     &                             + dedphiprime*dphidzid*dphidxic
     &                             + d2edphi2prime*dphidzid*dphidxic
     &                             + dtdphiprime*dxiczid
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid
     &                             + d2edphi2*dphidxid*dphidyic
     &                             + domegadyic*dxid
     &                             + dedphiprime*dphidxid*dphidyic
     &                             + d2edphi2prime*dphidxid*dphidyic
     &                             + dtdphiprime*dyicxid
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid
     &                             + d2edphi2*dphidyid*dphidyic
     &                             + domegadyic*dyid
     &                             + dedphiprime*dphidyid*dphidyic
     &                             + d2edphi2prime*dphidyid*dphidyic
     &                             + dtdphiprime*dyicyid
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid
     &                             + d2edphi2*dphidzid*dphidyic
     &                             + domegadyic*dzid
     &                             + dedphiprime*dphidzid*dphidyic
     &                             + d2edphi2prime*dphidzid*dphidyic
     &                             + dtdphiprime*dyiczid
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid
     &                             + d2edphi2*dphidxid*dphidzic
     &                             + domegadzic*dxid
     &                             + dedphiprime*dphidxid*dphidzic
     &                             + d2edphi2prime*dphidxid*dphidzic
     &                             + dtdphiprime*dzicxid
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid
     &                             + d2edphi2*dphidyid*dphidzic
     &                             + domegadzic*dyid
     &                             + dedphiprime*dphidyid*dphidzic
     &                             + d2edphi2prime*dphidyid*dphidzic
     &                             + dtdphiprime*dzicyid
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid
     &                             + d2edphi2*dphidzid*dphidzic
     &                             + domegadzic*dzid
     &                             + dedphiprime*dphidzid*dphidzic
     &                             + d2edphi2prime*dphidzid*dphidzic
     &                             + dtdphiprime*dziczid
               end if
c
c     get the angle-torsion values for the second angle
c
               v1 = kant(4,iangtor)
               v2 = kant(5,iangtor)
               v3 = kant(6,iangtor)
               k = iat(3,iangtor)
               dot = xbc*xdc + ybc*ydc + zbc*zdc
               cosang = dot / sqrt(rcb2*rdc2)
               angle = radian * acos(cosang)
               dt = angle - anat(k)
               force = ak(k)
               dedphi = atorunit * 2.0d0 * force
     &                     * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2 = atorunit * dt * 2.0d0 * force
     &                       * (v1*d2phi1 + v2* d2phi2 + v3*d2phi3)
               d2dt = atorunit * radian * 2.0d0 * force 
     &                   * (v1*phi1 + v2*phi2 + v3*phi3)
               dedphiprime = atorunit * 2.0d0 * force
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
               d2edphi2prime = atorunit * 2.0d0 * force
     &                            * (v1*phi1 + v2*phi2 + v3*phi3)
     &                            * (v1*d2phi1 + v2* d2phi2 + v3*d2phi3)
               dtdphiprime = atorunit * 2.0d0 * force
     &                          * (v1*phi1 + v2*phi2 + v3*phi3)
     &                          * (v1*dphi1 + v2*dphi2 + v3*dphi3)
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  dedphi = dedphi * fgrp
                  d2edphi2 = d2edphi2 * fgrp
                  d2dt = d2dt * fgrp   
                  dedphiprime = dedphiprime * fgrp
                  d2edphi2prime = d2edphi2prime * fgrp
                  dtdphiprime = dedphiprime * fgrp
               end if
c
c     first and second derivative components for the second angle
c
               termb = -1.0d0 / (rcb2*sqrt(ru2))
               termd = 1.0d0 / (rdc2*sqrt(ru2))
               domegadxib = termb * (ybc*zu - zbc*yu)
               domegadyib = termb * (zbc*xu - xbc*zu)
               domegadzib = termb * (xbc*yu - ybc*xu)
               domegadxid = termd * (ydc*zu - zdc*yu)
               domegadyid = termd * (zdc*xu - xdc*zu)
               domegadzid = termd * (xdc*yu - ydc*xu)
               domegadxic = -domegadxib - domegadxid
               domegadyic = -domegadyib - domegadyid
               domegadzic = -domegadzib - domegadzid
c
c     abbreviations used in defining chain rule terms
c
               xrbc = 2.0d0 * xbc / rcb2
               yrbc = 2.0d0 * ybc / rcb2
               zrbc = 2.0d0 * zbc / rcb2
               xrdc = 2.0d0 * xdc / rdc2
               yrdc = 2.0d0 * ydc / rdc2
               zrdc = 2.0d0 * zdc / rdc2
               xbcp = (ybc*zu-zbc*yu) / ru2
               ybcp = (zbc*xu-xbc*zu) / ru2
               zbcp = (xbc*yu-ybc*xu) / ru2
               xdcp = (ydc*zu-zdc*yu) / ru2
               ydcp = (zdc*xu-xdc*zu) / ru2
               zdcp = (xdc*yu-ydc*xu) / ru2
c
c     chain rule terms for second derivative components
c
               doxibxib = termb*(xbc*xdc-dot) + domegadxib*(xdcp-xrbc)
               doxibyib = termb*(zu+ybc*xdc) + domegadxib*(ydcp-yrbc)
               doxibzib = termb*(zbc*xdc-yu) + domegadxib*(zdcp-zrbc)
               doyibyib = termb*(ybc*ydc-dot) + domegadyib*(ydcp-yrbc)
               doyibzib = termb*(xu+zbc*ydc) + domegadyib*(zdcp-zrbc)
               dozibzib = termb*(zbc*zdc-dot) + domegadzib*(zdcp-zrbc)
               doxidxid = termd*(dot-xbc*xdc) - domegadxid*(xbcp+xrdc)
               doxidyid = termd*(zu-ydc*xbc) - domegadxid*(ybcp+yrdc)
               doxidzid = -termd*(yu+zdc*xbc) - domegadxid*(zbcp+zrdc)
               doyidyid = termd*(dot-ybc*ydc) - domegadyid*(ybcp+yrdc)
               doyidzid = termd*(xu-zdc*ybc) - domegadyid*(zbcp+zrdc)
               dozidzid = termd*(dot-zbc*zdc) - domegadzid*(zbcp+zrdc)
               doxibxid = termb*(ybc*ybc+zbc*zbc) - domegadxib*xbcp
               doxibyid = -termb*xbc*ybc - domegadxib*ybcp
               doxibzid = -termb*xbc*zbc - domegadxib*zbcp
               doyibxid = -termb*xbc*ybc - domegadyib*xbcp
               doyibyid = termb*(xbc*xbc+zbc*zbc) - domegadyib*ybcp
               doyibzid = -termb*ybc*zbc - domegadyib*zbcp
               dozibxid = -termb*xbc*zbc - domegadzib*xbcp
               dozibyid = -termb*ybc*zbc - domegadzib*ybcp
               dozibzid = termb*(xbc*xbc+ybc*ybc) - domegadzib*zbcp
c
c     get some second derivative chain rule terms by difference
c
               doxicxib = -doxibxib - doxibxid
               doxicyib = -doxibyib - doyibxid
               doxiczib = -doxibzib - dozibxid
               doyicxib = -doxibyib - doxibyid
               doyicyib = -doyibyib - doyibyid
               doyiczib = -doyibzib - dozibyid
               dozicxib = -doxibzib - doxibzid
               dozicyib = -doyibzib - doyibzid
               doziczib = -dozibzib - dozibzid
               doxicxid = -doxidxid - doxibxid
               doxicyid = -doxidyid - doxibyid
               doxiczid = -doxidzid - doxibzid
               doyicxid = -doxidyid - doyibxid
               doyicyid = -doyidyid - doyibyid
               doyiczid = -doyidzid - doyibzid
               dozicxid = -doxidzid - dozibxid
               dozicyid = -doyidzid - dozibyid
               doziczid = -dozidzid - dozibzid
               doxicxic = -doxicxib - doxicxid
               doxicyic = -doxicyib - doxicyid
               doxiczic = -doxiczib - doxiczid
               doyicyic = -doyicyib - doyicyid
               doyiczic = -doyiczib - doyiczid
               doziczic = -doziczib - doziczid
c
c     scale the first-derivatives of the second angle
c
               domegadxib = domegadxib * radian
               domegadyib = domegadyib * radian
               domegadzib = domegadzib * radian
               domegadxid = domegadxid * radian
               domegadyid = domegadyid * radian
               domegadzid = domegadzid * radian
               domegadxic = domegadxic * radian
               domegadyic = domegadyic * radian
               domegadzic = domegadzic * radian
c
c     scale the second-derivatives of the second angle
c
               doxibxib = doxibxib * d2dt
               doxibyib = doxibyib * d2dt
               doxibzib = doxibzib * d2dt
               doyibyib = doyibyib * d2dt
               doyibzib = doyibzib * d2dt
               dozibzib = dozibzib * d2dt
               doxidxid = doxidxid * d2dt
               doxidyid = doxidyid * d2dt
               doxidzid = doxidzid * d2dt
               doyidyid = doyidyid * d2dt
               doyidzid = doyidzid * d2dt
               dozidzid = dozidzid * d2dt
               doxibxid = doxibxid * d2dt
               doxibyid = doxibyid * d2dt
               doxibzid = doxibzid * d2dt
               doyibxid = doyibxid * d2dt
               doyibyid = doyibyid * d2dt
               doyibzid = doyibzid * d2dt
               dozibxid = dozibxid * d2dt
               dozibyid = dozibyid * d2dt
               dozibzid = dozibzid * d2dt
               doxicxib = doxicxib * d2dt
               doxicyib = doxicyib * d2dt
               doxiczib = doxiczib * d2dt
               doyicxib = doyicxib * d2dt
               doyicyib = doyicyib * d2dt
               doyiczib = doyiczib * d2dt
               dozicxib = dozicxib * d2dt
               dozicyib = dozicyib * d2dt
               doziczib = doziczib * d2dt
               doxicxid = doxicxid * d2dt
               doxicyid = doxicyid * d2dt
               doxiczid = doxiczid * d2dt
               doyicxid = doyicxid * d2dt
               doyicyid = doyicyid * d2dt
               doyiczid = doyiczid * d2dt
               dozicxid = dozicxid * d2dt
               dozicyid = dozicyid * d2dt
               doziczid = doziczid * d2dt
               doxicxic = doxicxic * d2dt
               doxicyic = doxicyic * d2dt
               doxiczic = doxiczic * d2dt
               doyicyic = doyicyic * d2dt
               doyiczic = doyiczic * d2dt
               doziczic = doziczic * d2dt
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
               dedphi = dedphi * dt
c
c     increment diagonal and off-diagonal Hessian elements
c
               if (i .eq. ia) then
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxia
     &                             + d2edphi2*dphidxia*dphidxia
     &                             + d2edphi2prime*dphidxia*dphidxia
     &                             + dedphiprime*dphidxia*dphidxia
     &                             + dtdphiprime*dxiaxia
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             + d2edphi2prime*dphidxia*dphidyia
     &                             + dedphiprime*dphidxia*dphidyia
     &                             + dtdphiprime*dxiayia
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             + d2edphi2prime*dphidxia*dphidzia
     &                             + dedphiprime*dphidxia*dphidzia
     &                             + dtdphiprime*dxiazia
                  hessx(2,ia) = hessx(2,ia) + dedphi*dxiayia
     &                             + d2edphi2*dphidxia*dphidyia
     &                             + d2edphi2prime*dphidxia*dphidyia
     &                             + dedphiprime*dphidxia*dphidyia
     &                             + dtdphiprime*dxiayia
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayia
     &                             + d2edphi2*dphidyia*dphidyia
     &                             + d2edphi2prime*dphidyia*dphidyia
     &                             + dedphiprime*dphidyia*dphidyia
     &                             + dtdphiprime*dyiayia
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             + d2edphi2prime*dphidyia*dphidzia
     &                             + dedphiprime*dphidyia*dphidzia
     &                             + dtdphiprime*dyiazia
                  hessx(3,ia) = hessx(3,ia) + dedphi*dxiazia
     &                             + d2edphi2*dphidxia*dphidzia
     &                             + d2edphi2prime*dphidxia*dphidzia
     &                             + dedphiprime*dphidxia*dphidzia
     &                             + dtdphiprime*dxiazia
                  hessy(3,ia) = hessy(3,ia) + dedphi*dyiazia
     &                             + d2edphi2*dphidyia*dphidzia
     &                             + d2edphi2prime*dphidyia*dphidzia
     &                             + dedphiprime*dphidyia*dphidzia
     &                             + dtdphiprime*dyiazia
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazia
     &                             + d2edphi2*dphidzia*dphidzia
     &                             + d2edphi2prime*dphidzia*dphidzia
     &                             + dedphiprime*dphidzia*dphidzia
     &                             + dtdphiprime*dziazia
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxia*dphidxib
     &                             + domegadxib*dxia
     &                             + d2edphi2prime*dphidxia*dphidxib
     &                             + dedphiprime*dphidxia*dphidxib
     &                             + dtdphiprime*dxiaxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dyiaxib
     &                             + d2edphi2*dphidyia*dphidxib
     &                             + domegadxib*dyia
     &                             + d2edphi2prime*dphidyia*dphidxib
     &                             + dedphiprime*dphidyia*dphidxib
     &                             + dtdphiprime*dyiaxib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dziaxib
     &                             + d2edphi2*dphidzia*dphidxib
     &                             + domegadxib*dzia
     &                             + d2edphi2prime*dphidzia*dphidxib
     &                             + dedphiprime*dphidzia*dphidxib 
     &                             + dtdphiprime*dziaxib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxiayib
     &                             + d2edphi2*dphidxia*dphidyib
     &                             + domegadyib*dxia
     &                             + d2edphi2prime*dphidxia*dphidyib
     &                             + dedphiprime*dphidxia*dphidyib
     &                             + dtdphiprime*dxiayib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyiayib
     &                             + d2edphi2*dphidyia*dphidyib
     &                             + domegadyib*dyia
     &                             + d2edphi2prime*dphidyia*dphidyib
     &                             + dedphiprime*dphidyia*dphidyib
     &                             + dtdphiprime*dyiayib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dziayib
     &                             + d2edphi2*dphidzia*dphidyib
     &                             + domegadyib*dzia
     &                             + d2edphi2prime*dphidzia*dphidyib
     &                             + dedphiprime*dphidzia*dphidyib
     &                             + dtdphiprime*dziayib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxiazib
     &                             + d2edphi2*dphidxia*dphidzib
     &                             + domegadzib*dxia
     &                             + d2edphi2prime*dphidxia*dphidzib
     &                             + dedphiprime*dphidxia*dphidzib
     &                             + dtdphiprime*dxiazib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyiazib
     &                             + d2edphi2*dphidyia*dphidzib
     &                             + domegadzib*dyia
     &                             + d2edphi2prime*dphidyia*dphidzib
     &                             + dedphiprime*dphidyia*dphidzib
     &                             + dtdphiprime*dyiazib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dziazib
     &                             + d2edphi2*dphidzia*dphidzib
     &                             + domegadzib*dzia
     &                             + d2edphi2prime*dphidzia*dphidzib 
     &                             + dedphiprime*dphidzia*dphidzib
     &                             + dtdphiprime*dziazib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxia*dphidxic
     &                             + domegadxic*dxia
     &                             + d2edphi2prime*dphidxia*dphidxic
     &                             + dedphiprime*dphidxia*dphidxic
     &                             + dtdphiprime*dxiaxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyiaxic
     &                             + d2edphi2*dphidyia*dphidxic
     &                             + domegadxic*dyia
     &                             + d2edphi2prime*dphidyia*dphidxic
     &                             + dedphiprime*dphidyia*dphidxic
     &                             + dtdphiprime*dyiaxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dziaxic
     &                             + d2edphi2*dphidzia*dphidxic
     &                             + domegadxic*dzia
     &                             + d2edphi2prime*dphidzia*dphidxic
     &                             + dedphiprime*dphidzia*dphidxic
     &                             + dtdphiprime*dziaxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxiayic
     &                             + d2edphi2*dphidxia*dphidyic
     &                             + domegadyic*dxia
     &                             + d2edphi2prime*dphidxia*dphidyic
     &                             + dedphiprime*dphidxia*dphidyic
     &                             + dtdphiprime*dxiayic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyiayic
     &                             + d2edphi2*dphidyia*dphidyic
     &                             + domegadyic*dyia
     &                             + d2edphi2prime*dphidyia*dphidyic
     &                             + dedphiprime*dphidyia*dphidyic
     &                             + dtdphiprime*dyiayic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dziayic
     &                             + d2edphi2*dphidzia*dphidyic
     &                             + domegadyic*dzia
     &                             + d2edphi2prime*dphidzia*dphidyic
     &                             + dedphiprime*dphidzia*dphidyic
     &                             + dtdphiprime*dziayic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiazic
     &                             + d2edphi2*dphidxia*dphidzic
     &                             + domegadzic*dxia
     &                             + d2edphi2prime*dphidxia*dphidzic
     &                             + dedphiprime*dphidxia*dphidzic
     &                             + dtdphiprime*dxiazic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiazic
     &                             + d2edphi2*dphidyia*dphidzic
     &                             + domegadzic*dyia
     &                             + d2edphi2prime*dphidyia*dphidzic
     &                             + dedphiprime*dphidyia*dphidzic
     &                             + dtdphiprime*dyiazic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziazic
     &                             + d2edphi2*dphidzia*dphidzic
     &                             + domegadzic*dzia
     &                             + d2edphi2prime*dphidzia*dphidzic
     &                             + dedphiprime*dphidzia*dphidzic
     &                             + dtdphiprime*dziazic
                  hessx(1,id) = hessx(1,id) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxia*dphidxid
     &                             + domegadxid*dxia
     &                             + d2edphi2prime*dphidxia*dphidxid
     &                             + dedphiprime*dphidxia*dphidxid
     &                             + dtdphiprime*dxiaxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyiaxid
     &                             + d2edphi2*dphidyia*dphidxid
     &                             + domegadxid*dyia
     &                             + d2edphi2prime*dphidyia*dphidxid
     &                             + dedphiprime*dphidyia*dphidxid
     &                             + dtdphiprime*dyiaxid
                  hessz(1,id) = hessz(1,id) + dedphi*dziaxid
     &                             + d2edphi2*dphidzia*dphidxid
     &                             + domegadxid*dzia
     &                             + d2edphi2prime*dphidzia*dphidxid
     &                             + dedphiprime*dphidzia*dphidxid
     &                             + dtdphiprime*dziaxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxiayid
     &                             + d2edphi2*dphidxia*dphidyid
     &                             + domegadyid*dxia
     &                             + d2edphi2prime*dphidxia*dphidyid
     &                             + dedphiprime*dphidxia*dphidyid
     &                             + dtdphiprime*dxiayid
                  hessy(2,id) = hessy(2,id) + dedphi*dyiayid
     &                             + d2edphi2*dphidyia*dphidyid
     &                             + domegadyid*dyia
     &                             + d2edphi2prime*dphidyia*dphidyid
     &                             + dedphiprime*dphidyia*dphidyid
     &                             + dtdphiprime*dyiayid
                  hessz(2,id) = hessz(2,id) + dedphi*dziayid
     &                             + d2edphi2*dphidzia*dphidyid
     &                             + domegadyid*dzia
     &                             + d2edphi2prime*dphidzia*dphidyid
     &                             + dedphiprime*dphidzia*dphidyid
     &                             + dtdphiprime*dziayid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiazid
     &                             + d2edphi2*dphidxia*dphidzid
     &                             + domegadzid*dxia
     &                             + d2edphi2prime*dphidxia*dphidzid
     &                             + dedphiprime*dphidxia*dphidzid
     &                             + dtdphiprime*dxiazid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiazid
     &                             + d2edphi2*dphidyia*dphidzid
     &                             + domegadzid*dyia
     &                             + d2edphi2prime*dphidyia*dphidzid
     &                             + dedphiprime*dphidyia*dphidzid
     &                             + dtdphiprime*dyiazid
                  hessz(3,id) = hessz(3,id) + dedphi*dziazid
     &                             + d2edphi2*dphidzia*dphidzid
     &                             + domegadzid*dzia
     &                             + d2edphi2prime*dphidzia*dphidzid
     &                             + dedphiprime*dphidzia*dphidzid
     &                             + dtdphiprime*dziazid
               else if (i .eq. ib) then
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxib + doxibxib
     &                             + d2edphi2*dphidxib*dphidxib
     &                             + 2.0d0*domegadxib*dxib
     &                             + d2edphi2prime*dphidxib*dphidxib
     &                             + dedphiprime*dphidxib*dphidxib
     &                             + dtdphiprime*dxibxib
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyib + doxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + domegadxib*dyib + domegadyib*dxib
     &                             + d2edphi2prime*dphidxib*dphidyib
     &                             + dedphiprime*dphidxib*dphidyib
     &                             + dtdphiprime*dxibyib
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzib + doxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + domegadxib*dzib + domegadzib*dxib
     &                             + d2edphi2prime*dphidxib*dphidzib
     &                             + dedphiprime*dphidxib*dphidzib
     &                             + dtdphiprime*dxibzib
                  hessx(2,ib) = hessx(2,ib) + dedphi*dxibyib
     &                             + d2edphi2*dphidxib*dphidyib
     &                             + domegadxib*dyib + domegadyib*dxib
     &                             + doxibyib
     &                             + d2edphi2prime*dphidxib*dphidyib
     &                             + dedphiprime*dphidxib*dphidyib
     &                             + dtdphiprime*dxibyib
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyib + doyibyib
     &                             + d2edphi2*dphidyib*dphidyib
     &                             + 2.0d0*domegadyib*dyib
     &                             + d2edphi2prime*dphidyib*dphidyib
     &                             + dedphiprime*dphidyib*dphidyib
     &                             + dtdphiprime*dyibyib
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzib + doyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + domegadyib*dzib + domegadzib*dyib
     &                             + d2edphi2prime*dphidyib*dphidzib
     &                             + dedphiprime*dphidyib*dphidzib
     &                             + dtdphiprime*dyibzib
                  hessx(3,ib) = hessx(3,ib) + dedphi*dxibzib + doxibzib
     &                             + d2edphi2*dphidxib*dphidzib
     &                             + domegadxib*dzib + domegadzib*dxib
     &                             + d2edphi2prime*dphidxib*dphidzib
     &                             + dedphiprime*dphidxib*dphidzib
     &                             + dtdphiprime*dxibzib
                  hessy(3,ib) = hessy(3,ib) + dedphi*dyibzib + doyibzib
     &                             + d2edphi2*dphidyib*dphidzib
     &                             + domegadyib*dzib + domegadzib*dyib
     &                             + d2edphi2prime*dphidyib*dphidzib
     &                             + dedphiprime*dphidyib*dphidzib
     &                             + dtdphiprime*dyibzib
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzib + dozibzib
     &                             + d2edphi2*dphidzib*dphidzib
     &                             + 2.0d0*domegadzib*dzib
     &                             + d2edphi2prime*dphidzib*dphidzib
     &                             + dedphiprime*dphidzib*dphidzib
     &                             + dtdphiprime*dzibzib
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxib
     &                             + d2edphi2*dphidxib*dphidxia
     &                             + domegadxib*dxia
     &                             + d2edphi2prime*dphidxib*dphidxia
     &                             + dedphiprime*dphidxib*dphidxia
     &                             + dtdphiprime*dxiaxib
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayib
     &                             + d2edphi2*dphidyib*dphidxia
     &                             + domegadyib*dxia
     &                             + d2edphi2prime*dphidyib*dphidxia
     &                             + dedphiprime*dphidyib*dphidxia
     &                             + dtdphiprime*dxiayib
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazib
     &                             + d2edphi2*dphidzib*dphidxia
     &                             + domegadzib*dxia
     &                             + d2edphi2prime*dphidzib*dphidxia
     &                             + dedphiprime*dphidzib*dphidxia
     &                             + dtdphiprime*dxiazib
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxib
     &                             + d2edphi2*dphidxib*dphidyia
     &                             + domegadxib*dyia
     &                             + d2edphi2prime*dphidxib*dphidyia
     &                             + dedphiprime*dphidxib*dphidyia
     &                             + dtdphiprime*dyiaxib
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayib
     &                             + d2edphi2*dphidyib*dphidyia
     &                             + domegadyib*dyia
     &                             + d2edphi2prime*dphidyib*dphidyia
     &                             + dedphiprime*dphidyib*dphidyia
     &                             + dtdphiprime*dyiayib
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazib
     &                             + d2edphi2*dphidzib*dphidyia
     &                             + domegadzib*dyia
     &                             + d2edphi2prime*dphidzib*dphidyia
     &                             + dedphiprime*dphidzib*dphidyia
     &                             + dtdphiprime*dyiazib
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxib
     &                             + d2edphi2*dphidxib*dphidzia
     &                             + domegadxib*dzia
     &                             + d2edphi2prime*dphidxib*dphidzia
     &                             + dedphiprime*dphidxib*dphidzia
     &                             + dtdphiprime*dziaxib
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayib
     &                             + d2edphi2*dphidyib*dphidzia
     &                             + domegadyib*dzia
     &                             + d2edphi2prime*dphidyib*dphidzia
     &                             + dedphiprime*dphidyib*dphidzia
     &                             + dtdphiprime*dziayib
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazib
     &                             + d2edphi2*dphidzib*dphidzia
     &                             + domegadzib*dzia
     &                             + d2edphi2prime*dphidzib*dphidzia
     &                             + dedphiprime*dphidzib*dphidzia
     &                             + dtdphiprime*dziazib
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxibxic + doxicxib
     &                             + d2edphi2*dphidxib*dphidxic
     &                             + domegadxib*dxic + domegadxic*dxib
     &                             + d2edphi2prime*dphidxib*dphidxic
     &                             + dedphiprime*dphidxib*dphidxic
     &                             + dtdphiprime*dxibxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dyibxic + doxicyib
     &                             + d2edphi2*dphidyib*dphidxic
     &                             + domegadyib*dxic + domegadxic*dyib
     &                             + d2edphi2prime*dphidyib*dphidxic
     &                             + dedphiprime*dphidyib*dphidxic
     &                             + dtdphiprime*dyibxic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dzibxic + doxiczib
     &                             + d2edphi2*dphidzib*dphidxic
     &                             + domegadzib*dxic + domegadxic*dzib
     &                             + d2edphi2prime*dphidzib*dphidxic
     &                             + dedphiprime*dphidzib*dphidxic
     &                             + dtdphiprime*dzibxic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxibyic + doyicxib
     &                             + d2edphi2*dphidxib*dphidyic
     &                             + domegadxib*dyic + domegadyic*dxib
     &                             + d2edphi2prime*dphidxib*dphidyic
     &                             + dedphiprime*dphidxib*dphidyic
     &                             + dtdphiprime*dxibyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyibyic + doyicyib
     &                             + d2edphi2*dphidyib*dphidyic
     &                             + domegadyib*dyic + domegadyic*dyib
     &                             + d2edphi2prime*dphidyib*dphidyic
     &                             + dedphiprime*dphidyib*dphidyic
     &                             + dtdphiprime*dyibyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dzibyic + doyiczib
     &                             + d2edphi2*dphidzib*dphidyic
     &                             + domegadzib*dyic + domegadyic*dzib
     &                             + d2edphi2prime*dphidzib*dphidyic
     &                             + dedphiprime*dphidzib*dphidyic
     &                             + dtdphiprime*dzibyic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxibzic + dozicxib
     &                             + d2edphi2*dphidxib*dphidzic
     &                             + domegadxib*dzic + domegadzic*dxib
     &                             + d2edphi2prime*dphidxib*dphidzic
     &                             + dedphiprime*dphidxib*dphidzic
     &                             + dtdphiprime*dxibzic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyibzic + dozicyib
     &                             + d2edphi2*dphidyib*dphidzic
     &                             + domegadyib*dzic + domegadzic*dyib
     &                             + d2edphi2prime*dphidyib*dphidzic
     &                             + dedphiprime*dphidyib*dphidzic
     &                             + dtdphiprime*dyibzic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dzibzic + doziczib
     &                             + d2edphi2*dphidzib*dphidzic
     &                             + domegadzib*dzic + domegadzic*dzib
     &                             + d2edphi2prime*dphidzib*dphidzic
     &                             + dedphiprime*dphidzib*dphidzic
     &                             + dtdphiprime*dzibzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxibxid + doxibxid
     &                             + d2edphi2*dphidxib*dphidxid
     &                             + domegadxib*dxid + domegadxid*dxib
     &                             + d2edphi2prime*dphidxib*dphidxid
     &                             + dedphiprime*dphidxib*dphidxid
     &                             + dtdphiprime*dxibxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyibxid + doyibxid
     &                             + d2edphi2*dphidyib*dphidxid
     &                             + domegadyib*dxid + domegadxid*dyib
     &                             + d2edphi2prime*dphidyib*dphidxid
     &                             + dedphiprime*dphidyib*dphidxid
     &                             + dtdphiprime*dyibxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzibxid + dozibxid
     &                             + d2edphi2*dphidzib*dphidxid
     &                             + domegadzib*dxid + domegadxid*dzib
     &                             + d2edphi2prime*dphidzib*dphidxid
     &                             + dedphiprime*dphidzib*dphidxid
     &                             + dtdphiprime*dzibxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxibyid + doxibyid
     &                             + d2edphi2*dphidxib*dphidyid
     &                             + domegadxib*dyid + domegadyid*dxib
     &                             + d2edphi2prime*dphidxib*dphidyid
     &                             + dedphiprime*dphidxib*dphidyid
     &                             + dtdphiprime*dxibyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyibyid + doyibyid
     &                             + d2edphi2*dphidyib*dphidyid
     &                             + domegadyib*dyid + domegadyid*dyib
     &                             + d2edphi2prime*dphidyib*dphidyid
     &                             + dedphiprime*dphidyib*dphidyid
     &                             + dtdphiprime*dyibyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzibyid + dozibyid
     &                             + d2edphi2*dphidzib*dphidyid
     &                             + domegadzib*dyid + domegadyid*dzib
     &                             + d2edphi2prime*dphidzib*dphidyid
     &                             + dedphiprime*dphidzib*dphidyid
     &                             + dtdphiprime*dzibyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxibzid + doxibzid
     &                             + d2edphi2*dphidxib*dphidzid
     &                             + domegadxib*dzid + domegadzid*dxib
     &                             + d2edphi2prime*dphidxib*dphidzid
     &                             + dedphiprime*dphidxib*dphidzid
     &                             + dtdphiprime*dxibzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyibzid + doyibzid
     &                             + d2edphi2*dphidyib*dphidzid
     &                             + domegadyib*dzid + domegadzid*dyib
     &                             + d2edphi2prime*dphidyib*dphidzid
     &                             + dedphiprime*dphidyib*dphidzid
     &                             + dtdphiprime*dyibzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzibzid + dozibzid
     &                             + d2edphi2*dphidzib*dphidzid
     &                             + domegadzib*dzid + domegadzid*dzib
     &                             + d2edphi2prime*dphidzib*dphidzid
     &                             + dedphiprime*dphidzib*dphidzid
     &                             + dtdphiprime*dzibzid
               else if (i .eq. ic) then
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxic + doxicxic
     &                             + d2edphi2*dphidxic*dphidxic
     &                             + 2.0d0*domegadxic*dxic
     &                             + d2edphi2prime*dphidxic*dphidxic
     &                             + dedphiprime*dphidxic*dphidxic
     &                             + dtdphiprime*dxicxic
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyic + doxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + domegadxic*dyic + domegadyic*dxic
     &                             + d2edphi2prime*dphidxic*dphidyic
     &                             + dedphiprime*dphidxic*dphidyic
     &                             + dtdphiprime*dxicyic
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczic + doxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + domegadxic*dzic + domegadzic*dxic
     &                             + d2edphi2prime*dphidxic*dphidzic
     &                             + dedphiprime*dphidxic*dphidzic
     &                             + dtdphiprime*dxiczic
                  hessx(2,ic) = hessx(2,ic) + dedphi*dxicyic + doxicyic
     &                             + d2edphi2*dphidxic*dphidyic
     &                             + domegadxic*dyic + domegadyic*dxic
     &                             + d2edphi2prime*dphidxic*dphidyic
     &                             + dedphiprime*dphidxic*dphidyic
     &                             + dtdphiprime*dxicyic
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyic + doyicyic
     &                             + d2edphi2*dphidyic*dphidyic
     &                             + 2.0d0*domegadyic*dyic
     &                             + d2edphi2prime*dphidyic*dphidyic
     &                             + dedphiprime*dphidyic*dphidyic
     &                             + dtdphiprime*dyicyic
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczic + doyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + domegadyic*dzic + domegadzic*dyic
     &                             + d2edphi2prime*dphidyic*dphidzic
     &                             + dedphiprime*dphidyic*dphidzic
     &                             + dtdphiprime*dyiczic
                  hessx(3,ic) = hessx(3,ic) + dedphi*dxiczic + doxiczic
     &                             + d2edphi2*dphidxic*dphidzic
     &                             + domegadxic*dzic + domegadzic*dxic
     &                             + d2edphi2prime*dphidxic*dphidzic
     &                             + dedphiprime*dphidxic*dphidzic
     &                             + dtdphiprime*dxiczic
                  hessy(3,ic) = hessy(3,ic) + dedphi*dyiczic + doyiczic
     &                             + d2edphi2*dphidyic*dphidzic
     &                             + domegadyic*dzic + domegadzic*dyic
     &                             + d2edphi2prime*dphidyic*dphidzic
     &                             + dedphiprime*dphidyic*dphidzic
     &                             + dtdphiprime*dyiczic
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczic + doziczic
     &                             + d2edphi2*dphidzic*dphidzic
     &                             + 2.0d0*domegadzic*dzic
     &                             + d2edphi2prime*dphidzic*dphidzic
     &                             + dedphiprime*dphidzic*dphidzic
     &                             + dtdphiprime*dziczic
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxic
     &                             + d2edphi2*dphidxic*dphidxia
     &                             + domegadxic*dxia
     &                             + d2edphi2prime*dphidxic*dphidxia
     &                             + dedphiprime*dphidxic*dphidxia
     &                             + dtdphiprime*dxiaxic
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayic
     &                             + d2edphi2*dphidyic*dphidxia
     &                             + domegadyic*dxia
     &                             + d2edphi2prime*dphidyic*dphidxia
     &                             + dedphiprime*dphidyic*dphidxia
     &                             + dtdphiprime*dxiayic
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazic
     &                             + d2edphi2*dphidzic*dphidxia
     &                             + domegadzic*dxia
     &                             + d2edphi2prime*dphidzic*dphidxia
     &                             + dedphiprime*dphidzic*dphidxia
     &                             + dtdphiprime*dxiazic
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxic
     &                             + d2edphi2*dphidxic*dphidyia
     &                             + domegadxic*dyia
     &                             + d2edphi2prime*dphidxic*dphidyia
     &                             + dedphiprime*dphidxic*dphidyia
     &                             + dtdphiprime*dyiaxic
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayic
     &                             + d2edphi2*dphidyic*dphidyia
     &                             + domegadyic*dyia
     &                             + d2edphi2prime*dphidyic*dphidyia
     &                             + dedphiprime*dphidyic*dphidyia
     &                             + dtdphiprime*dyiayic
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazic
     &                             + d2edphi2*dphidzic*dphidyia
     &                             + domegadzic*dyia
     &                             + d2edphi2prime*dphidzic*dphidyia
     &                             + dedphiprime*dphidzic*dphidyia
     &                             + dtdphiprime*dyiazic
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxic
     &                             + d2edphi2*dphidxic*dphidzia
     &                             + domegadxic*dzia
     &                             + d2edphi2prime*dphidxic*dphidzia
     &                             + dedphiprime*dphidxic*dphidzia
     &                             + dtdphiprime*dziaxic
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayic
     &                             + d2edphi2*dphidyic*dphidzia
     &                             + domegadyic*dzia
     &                             + d2edphi2prime*dphidyic*dphidzia
     &                             + dedphiprime*dphidyic*dphidzia
     &                             + dtdphiprime*dziayic
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazic
     &                             + d2edphi2*dphidzic*dphidzia
     &                             + domegadzic*dzia
     &                             + d2edphi2prime*dphidzic*dphidzia
     &                             + dedphiprime*dphidzic*dphidzia
     &                             + dtdphiprime*dziazic
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxic
     &                             + d2edphi2*dphidxic*dphidxib
     &                             + domegadxib*dxic + domegadxic*dxib
     &                             + doxicxib
     &                             + d2edphi2prime*dphidxic*dphidxib
     &                             + dedphiprime*dphidxic*dphidxib
     &                             + dtdphiprime*dxibxic
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyic + doyicxib
     &                             + d2edphi2*dphidyic*dphidxib
     &                             + domegadxib*dyic + domegadyic*dxib
     &                             + d2edphi2prime*dphidyic*dphidxib
     &                             + dedphiprime*dphidyic*dphidxib
     &                             + dtdphiprime*dxibyic
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzic + dozicxib
     &                             + d2edphi2*dphidzic*dphidxib
     &                             + domegadxib*dzic + domegadzic*dxib
     &                             + d2edphi2prime*dphidzic*dphidxib
     &                             + dedphiprime*dphidzic*dphidxib
     &                             + dtdphiprime*dxibzic
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxic + doxicyib
     &                             + d2edphi2*dphidxic*dphidyib
     &                             + domegadyib*dxic + domegadxic*dyib
     &                             + d2edphi2prime*dphidxic*dphidyib
     &                             + dedphiprime*dphidxic*dphidyib
     &                             + dtdphiprime*dyibxic
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyic + doyicyib
     &                             + d2edphi2*dphidyic*dphidyib
     &                             + domegadyib*dyic + domegadyic*dyib
     &                             + d2edphi2prime*dphidyic*dphidyib
     &                             + dedphiprime*dphidyic*dphidyib
     &                             + dtdphiprime*dyibyic
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzic + dozicyib
     &                             + d2edphi2*dphidzic*dphidyib
     &                             + domegadyib*dzic + domegadzic*dyib
     &                             + d2edphi2prime*dphidzic*dphidyib
     &                             + dedphiprime*dphidzic*dphidyib
     &                             + dtdphiprime*dyibzic
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxic + doxiczib
     &                             + d2edphi2*dphidxic*dphidzib
     &                             + domegadzib*dxic + domegadxic*dzib
     &                             + d2edphi2prime*dphidxic*dphidzib
     &                             + dedphiprime*dphidxic*dphidzib
     &                             + dtdphiprime*dzibxic
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyic + doyiczib
     &                             + d2edphi2*dphidyic*dphidzib
     &                             + domegadzib*dyic + domegadyic*dzib
     &                             + d2edphi2prime*dphidyic*dphidzib
     &                             + dedphiprime*dphidyic*dphidzib
     &                             + dtdphiprime*dzibyic
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzic + doziczib
     &                             + d2edphi2*dphidzic*dphidzib
     &                             + domegadzib*dzic + domegadzic*dzib
     &                             + d2edphi2prime*dphidzic*dphidzib
     &                             + dedphiprime*dphidzic*dphidzib
     &                             + dtdphiprime*dzibzic
                  hessx(1,id) = hessx(1,id) + dedphi*dxicxid + doxicxid
     &                             + d2edphi2*dphidxic*dphidxid
     &                             + domegadxic*dxid + domegadxid*dxic
     &                             + d2edphi2prime*dphidxic*dphidxid
     &                             + dedphiprime*dphidxic*dphidxid
     &                             + dtdphiprime*dxicxid
                  hessy(1,id) = hessy(1,id) + dedphi*dyicxid + doyicxid
     &                             + d2edphi2*dphidyic*dphidxid
     &                             + domegadyic*dxid + domegadxid*dyic
     &                             + d2edphi2prime*dphidyic*dphidxid
     &                             + dedphiprime*dphidyic*dphidxid
     &                             + dtdphiprime*dyicxid
                  hessz(1,id) = hessz(1,id) + dedphi*dzicxid + dozicxid
     &                             + d2edphi2*dphidzic*dphidxid
     &                             + domegadzic*dxid + domegadxid*dzic
     &                             + d2edphi2prime*dphidzic*dphidxid
     &                             + dedphiprime*dphidzic*dphidxid
     &                             + dtdphiprime*dzicxid
                  hessx(2,id) = hessx(2,id) + dedphi*dxicyid + doxicyid
     &                             + d2edphi2*dphidxic*dphidyid
     &                             + domegadxic*dyid + domegadyid*dxic
     &                             + d2edphi2prime*dphidxic*dphidyid
     &                             + dedphiprime*dphidxic*dphidyid
     &                             + dtdphiprime*dxicyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyicyid + doyicyid
     &                             + d2edphi2*dphidyic*dphidyid
     &                             + domegadyic*dyid + domegadyid*dyic
     &                             + d2edphi2prime*dphidyic*dphidyid
     &                             + dedphiprime*dphidyic*dphidyid
     &                             + dtdphiprime*dyicyid
                  hessz(2,id) = hessz(2,id) + dedphi*dzicyid + dozicyid
     &                             + d2edphi2*dphidzic*dphidyid
     &                             + domegadzic*dyid + domegadyid*dzic
     &                             + d2edphi2prime*dphidzic*dphidyid
     &                             + dedphiprime*dphidzic*dphidyid
     &                             + dtdphiprime*dzicyid
                  hessx(3,id) = hessx(3,id) + dedphi*dxiczid + doxiczid
     &                             + d2edphi2*dphidxic*dphidzid
     &                             + domegadxic*dzid + domegadzid*dxic
     &                             + d2edphi2prime*dphidxic*dphidzid
     &                             + dedphiprime*dphidxic*dphidzid
     &                             + dtdphiprime*dxiczid
                  hessy(3,id) = hessy(3,id) + dedphi*dyiczid + doyiczid
     &                             + d2edphi2*dphidyic*dphidzid
     &                             + domegadyic*dzid + domegadzid*dyic
     &                             + d2edphi2prime*dphidyic*dphidzid
     &                             + dedphiprime*dphidyic*dphidzid
     &                             + dtdphiprime*dyiczid
                  hessz(3,id) = hessz(3,id) + dedphi*dziczid + doziczid
     &                             + d2edphi2*dphidzic*dphidzid
     &                             + domegadzic*dzid + domegadzid*dzic
     &                             + d2edphi2prime*dphidzic*dphidzid
     &                             + dedphiprime*dphidzic*dphidzid
     &                             + dtdphiprime*dziczid
               else if (i .eq. id) then
                  hessx(1,id) = hessx(1,id) + dedphi*dxidxid + doxidxid
     &                             + d2edphi2*dphidxid*dphidxid
     &                             + 2.0d0*domegadxid*dxid
     &                             + d2edphi2prime*dphidxid*dphidxid
     &                             + dedphiprime*dphidxid*dphidxid
     &                             + dtdphiprime*dxidxid
                  hessy(1,id) = hessy(1,id) + dedphi*dxidyid + doxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + domegadxid*dyid + domegadyid*dxid
     &                             + d2edphi2prime*dphidxid*dphidyid
     &                             + dedphiprime*dphidxid*dphidyid
     &                             + dtdphiprime*dxidyid
                  hessz(1,id) = hessz(1,id) + dedphi*dxidzid + doxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + domegadxid*dzid + domegadzid*dxid
     &                             + d2edphi2prime*dphidxid*dphidzid
     &                             + dedphiprime*dphidxid*dphidzid
     &                             + dtdphiprime*dxidzid
                  hessx(2,id) = hessx(2,id) + dedphi*dxidyid + doxidyid
     &                             + d2edphi2*dphidxid*dphidyid
     &                             + domegadxid*dyid + domegadyid*dxid
     &                             + d2edphi2prime*dphidxid*dphidyid
     &                             + dedphiprime*dphidxid*dphidyid
     &                             + dtdphiprime*dxidyid
                  hessy(2,id) = hessy(2,id) + dedphi*dyidyid + doyidyid
     &                             + d2edphi2*dphidyid*dphidyid
     &                             + 2.0d0*domegadyid*dyid
     &                             + d2edphi2prime*dphidyid*dphidyid
     &                             + dedphiprime*dphidyid*dphidyid
     &                             + dtdphiprime*dyidyid
                  hessz(2,id) = hessz(2,id) + dedphi*dyidzid + doyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + domegadyid*dzid + domegadzid*dyid
     &                             + d2edphi2prime*dphidyid*dphidzid
     &                             + dedphiprime*dphidyid*dphidzid
     &                             + dtdphiprime*dyidzid
                  hessx(3,id) = hessx(3,id) + dedphi*dxidzid + doxidzid
     &                             + d2edphi2*dphidxid*dphidzid
     &                             + domegadxid*dzid + domegadzid*dxid
     &                             + d2edphi2prime*dphidxid*dphidzid
     &                             + dedphiprime*dphidxid*dphidzid
     &                             + dtdphiprime*dxidzid
                  hessy(3,id) = hessy(3,id) + dedphi*dyidzid + doyidzid
     &                             + d2edphi2*dphidyid*dphidzid
     &                             + domegadyid*dzid + domegadzid*dyid
     &                             + d2edphi2prime*dphidyid*dphidzid
     &                             + dedphiprime*dphidyid*dphidzid
     &                             + dtdphiprime*dyidzid
                  hessz(3,id) = hessz(3,id) + dedphi*dzidzid + dozidzid
     &                             + d2edphi2*dphidzid*dphidzid
     &                             + 2.0d0*domegadzid*dzid
     &                             + d2edphi2prime*dphidzid*dphidzid
     &                             + dedphiprime*dphidzid*dphidzid
     &                             + dtdphiprime*dzidzid
                  hessx(1,ia) = hessx(1,ia) + dedphi*dxiaxid
     &                             + d2edphi2*dphidxid*dphidxia
     &                             + domegadxid*dxia
     &                             + d2edphi2prime*dphidxid*dphidxia
     &                             + dedphiprime*dphidxid*dphidxia
     &                             + dtdphiprime*dxiaxid
                  hessy(1,ia) = hessy(1,ia) + dedphi*dxiayid
     &                             + d2edphi2*dphidyid*dphidxia
     &                             + domegadyid*dxia
     &                             + d2edphi2prime*dphidyid*dphidxia
     &                             + dedphiprime*dphidyid*dphidxia
     &                             + dtdphiprime*dxiayid
                  hessz(1,ia) = hessz(1,ia) + dedphi*dxiazid
     &                             + d2edphi2*dphidzid*dphidxia
     &                             + domegadzid*dxia
     &                             + d2edphi2prime*dphidzid*dphidxia
     &                             + dedphiprime*dphidzid*dphidxia
     &                             + dtdphiprime*dxiazid
                  hessx(2,ia) = hessx(2,ia) + dedphi*dyiaxid
     &                             + d2edphi2*dphidxid*dphidyia
     &                             + domegadxid*dyia
     &                             + d2edphi2prime*dphidxid*dphidyia
     &                             + dedphiprime*dphidxid*dphidyia
     &                             + dtdphiprime*dyiaxid
                  hessy(2,ia) = hessy(2,ia) + dedphi*dyiayid
     &                             + d2edphi2*dphidyid*dphidyia
     &                             + domegadyid*dyia
     &                             + d2edphi2prime*dphidyid*dphidyia
     &                             + dedphiprime*dphidyid*dphidyia
     &                             + dtdphiprime*dyiayid
                  hessz(2,ia) = hessz(2,ia) + dedphi*dyiazid
     &                             + d2edphi2*dphidzid*dphidyia
     &                             + domegadzid*dyia
     &                             + d2edphi2prime*dphidzid*dphidyia
     &                             + dedphiprime*dphidzid*dphidyia
     &                             + dtdphiprime*dyiazid
                  hessx(3,ia) = hessx(3,ia) + dedphi*dziaxid
     &                             + d2edphi2*dphidxid*dphidzia
     &                             + domegadxid*dzia
     &                             + d2edphi2prime*dphidxid*dphidzia
     &                             + dedphiprime*dphidxid*dphidzia
     &                             + dtdphiprime*dziaxid
                  hessy(3,ia) = hessy(3,ia) + dedphi*dziayid
     &                             + d2edphi2*dphidyid*dphidzia
     &                             + domegadyid*dzia
     &                             + d2edphi2prime*dphidyid*dphidzia
     &                             + dedphiprime*dphidyid*dphidzia
     &                             + dtdphiprime*dziayid
                  hessz(3,ia) = hessz(3,ia) + dedphi*dziazid
     &                             + d2edphi2*dphidzid*dphidzia
     &                             + domegadzid*dzia
     &                             + d2edphi2prime*dphidzid*dphidzia
     &                             + dedphiprime*dphidzid*dphidzia
     &                             + dtdphiprime*dziazid
                  hessx(1,ib) = hessx(1,ib) + dedphi*dxibxid + doxibxid
     &                             + d2edphi2*dphidxid*dphidxib
     &                             + domegadxib*dxid + domegadxid*dxib
     &                             + d2edphi2prime*dphidxid*dphidxib
     &                             + dedphiprime*dphidxid*dphidxib
     &                             + dtdphiprime*dxibxid
                  hessy(1,ib) = hessy(1,ib) + dedphi*dxibyid + doxibyid
     &                             + d2edphi2*dphidyid*dphidxib
     &                             + domegadxib*dyid + domegadyid*dxib
     &                             + d2edphi2prime*dphidyid*dphidxib
     &                             + dedphiprime*dphidyid*dphidxib
     &                             + dtdphiprime*dxibyid
                  hessz(1,ib) = hessz(1,ib) + dedphi*dxibzid + doxibzid
     &                             + d2edphi2*dphidzid*dphidxib
     &                             + domegadxib*dzid + domegadzid*dxib
     &                             + d2edphi2prime*dphidzid*dphidxib
     &                             + dedphiprime*dphidzid*dphidxib
     &                             + dtdphiprime*dxibzid
                  hessx(2,ib) = hessx(2,ib) + dedphi*dyibxid + doyibxid
     &                             + d2edphi2*dphidxid*dphidyib
     &                             + domegadyib*dxid + domegadxid*dyib
     &                             + d2edphi2prime*dphidxid*dphidyib
     &                             + dedphiprime*dphidxid*dphidyib
     &                             + dtdphiprime*dyibxid
                  hessy(2,ib) = hessy(2,ib) + dedphi*dyibyid + doyibyid
     &                             + d2edphi2*dphidyid*dphidyib
     &                             + domegadyib*dyid + domegadyid*dyib
     &                             + d2edphi2prime*dphidyid*dphidyib
     &                             + dedphiprime*dphidyid*dphidyib
     &                             + dtdphiprime*dyibyid
                  hessz(2,ib) = hessz(2,ib) + dedphi*dyibzid + doyibzid
     &                             + d2edphi2*dphidzid*dphidyib
     &                             + domegadyib*dzid + domegadzid*dyib
     &                             + d2edphi2prime*dphidzid*dphidyib
     &                             + dedphiprime*dphidzid*dphidyib
     &                             + dtdphiprime*dyibzid
                  hessx(3,ib) = hessx(3,ib) + dedphi*dzibxid + dozibxid
     &                             + d2edphi2*dphidxid*dphidzib
     &                             + domegadzib*dxid + domegadxid*dzib
     &                             + d2edphi2prime*dphidxid*dphidzib
     &                             + dedphiprime*dphidxid*dphidzib
     &                             + dtdphiprime*dzibxid
                  hessy(3,ib) = hessy(3,ib) + dedphi*dzibyid + dozibyid
     &                             + d2edphi2*dphidyid*dphidzib
     &                             + domegadzib*dyid + domegadyid*dzib
     &                             + d2edphi2prime*dphidyid*dphidzib
     &                             + dedphiprime*dphidyid*dphidzib
     &                             + dtdphiprime*dzibyid
                  hessz(3,ib) = hessz(3,ib) + dedphi*dzibzid + dozibzid
     &                             + d2edphi2*dphidzid*dphidzib
     &                             + domegadzib*dzid + domegadzid*dzib
     &                             + d2edphi2prime*dphidzid*dphidzib
     &                             + dedphiprime*dphidzid*dphidzib
     &                             + dtdphiprime*dzibzid
                  hessx(1,ic) = hessx(1,ic) + dedphi*dxicxid + doxicxid
     &                             + d2edphi2*dphidxid*dphidxic
     &                             + domegadxic*dxid + domegadxid*dxic
     &                             + d2edphi2prime*dphidxid*dphidxic
     &                             + dedphiprime*dphidxid*dphidxic
     &                             + dtdphiprime*dxicxid
                  hessy(1,ic) = hessy(1,ic) + dedphi*dxicyid + doxicyid
     &                             + d2edphi2*dphidyid*dphidxic
     &                             + domegadxic*dyid + domegadyid*dxic
     &                             + d2edphi2prime*dphidyid*dphidxic
     &                             + dedphiprime*dphidyid*dphidxic
     &                             + dtdphiprime*dxicyid
                  hessz(1,ic) = hessz(1,ic) + dedphi*dxiczid + doxiczid
     &                             + d2edphi2*dphidzid*dphidxic
     &                             + domegadxic*dzid + domegadzid*dxic
     &                             + d2edphi2prime*dphidzid*dphidxic
     &                             + dedphiprime*dphidzid*dphidxic
     &                             + dtdphiprime*dxiczid
                  hessx(2,ic) = hessx(2,ic) + dedphi*dyicxid + doyicxid
     &                             + d2edphi2*dphidxid*dphidyic
     &                             + domegadyic*dxid + domegadxid*dyic
     &                             + d2edphi2prime*dphidxid*dphidyic 
     &                             + dedphiprime*dphidxid*dphidyic
     &                             + dtdphiprime*dyicxid
                  hessy(2,ic) = hessy(2,ic) + dedphi*dyicyid + doyicyid
     &                             + d2edphi2*dphidyid*dphidyic
     &                             + domegadyic*dyid + domegadyid*dyic
     &                             + d2edphi2prime*dphidyid*dphidyic
     &                             + dedphiprime*dphidyid*dphidyic
     &                             + dtdphiprime*dyicyid
                  hessz(2,ic) = hessz(2,ic) + dedphi*dyiczid + doyiczid
     &                             + d2edphi2*dphidzid*dphidyic
     &                             + domegadyic*dzid + domegadzid*dyic
     &                             + d2edphi2prime*dphidzid*dphidyic
     &                             + dedphiprime*dphidzid*dphidyic
     &                             + dtdphiprime*dyiczid
                  hessx(3,ic) = hessx(3,ic) + dedphi*dzicxid + dozicxid
     &                             + d2edphi2*dphidxid*dphidzic
     &                             + domegadzic*dxid + domegadxid*dzic
     &                             + d2edphi2prime*dphidxid*dphidzic
     &                             + dedphiprime*dphidxid*dphidzic
     &                             + dtdphiprime*dzicxid
                  hessy(3,ic) = hessy(3,ic) + dedphi*dzicyid + dozicyid
     &                             + d2edphi2*dphidyid*dphidzic
     &                             + domegadzic*dyid + domegadyid*dzic
     &                             + d2edphi2prime*dphidyid*dphidzic
     &                             + dedphiprime*dphidyid*dphidzic
     &                             + dtdphiprime*dzicyid
                  hessz(3,ic) = hessz(3,ic) + dedphi*dziczid + doziczid
     &                             + d2edphi2*dphidzid*dphidzic
     &                             + domegadzic*dzid + domegadzid*dzic
     &                             + d2edphi2prime*dphidzid*dphidzic
     &                             + dedphiprime*dphidzid*dphidzic
     &                             + dtdphiprime*dziczid
               end if
            end if
         end if
      end do
      return
      end

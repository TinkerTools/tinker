c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine hessfull  --  calculates the Hessian elements  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "hessfull" calls subroutines to calculate the Hessian matrix
c     of potential energy with respect to cartesian coordinates
c
c
      subroutine hessfull (hessn)
      implicit none
      include 'sizes.for'
      integer i,j,k,l,n,inp,iout
      real*8 hessn(3,maxhes,3,maxhes)
      real*8 heb,hea,heba,het,hev,he14,hec,hecd,hed,hex
      common /atoms / n
      common /hessn1/ heb(3,maxhes,3,maxhes),hea(3,maxhes,3,maxhes),
     &                heba(3,maxhes,3,maxhes),het(3,maxhes,3,maxhes),
     &                hev(3,maxhes,3,maxhes),he14(3,maxhes,3,maxhes),
     &                hec(3,maxhes,3,maxhes),hecd(3,maxhes,3,maxhes),
     &                hed(3,maxhes,3,maxhes),hex(3,maxhes,3,maxhes)
      common /iodev / inp,iout
c
c
c     print message and exit if too many atoms for hessian
c
      if (n .gt. maxhes) then
         write (iout,10)  n,maxhes
   10    format (/' ******* Too many HESSIAN --> Molecule '
     &              'contains',i4,', Current maximum is',i4,
     &           /20x,'******* Increase "MAXHES" *******'//)
         stop
      end if
c
c     calculate the hessian matrix elements
c
      call ebond5
      call eangle5
      call etors5
      call evdw5
      call charge5
      call chgdpl5
      call dipole5
      call extra5
c
c     now sum up to give the total hessian matrix
c        (diagonal and above diagonal elements)
c
      do i = 1, n
         do j = 1, 3
            do l = j, 3
               hessn(j,i,l,i) = heb(j,i,l,i) + hea(j,i,l,i) +
     &             heba(j,i,l,i) + het(j,i,l,i) + hev(j,i,l,i) +
     &             he14(j,i,l,i) + hec(j,i,l,i) + hecd(j,i,l,i) +
     &             hed(j,i,l,i) + hex(j,i,l,i)
            end do
         end do
         do k = i+1, n
            do j = 1, 3
               do l = 1, 3
                  hessn(j,i,l,k) = heb(j,i,l,k) + hea(j,i,l,k) +
     &                heba(j,i,l,k) + het(j,i,l,k) + hev(j,i,l,k) +
     &                he14(j,i,l,k) + hec(j,i,l,k) + hecd(j,i,l,k) +
     &                hed(j,i,l,k) + hex(j,i,l,k)
               end do
            end do
         end do
      end do
c
c     finally, set below diagonal components and total Hessian
c        (this version does not use symmetry of "hex" term)
c
      do i = 1, n
         do j = 2, 3
            do l = j-1, 1, -1
               heb(j,i,l,i) = heb(l,i,j,i)
               hea(j,i,l,i) = hea(l,i,j,i)
               heba(j,i,l,i) = heba(l,i,j,i)
               het(j,i,l,i) = het(l,i,j,i)
               hev(j,i,l,i) = hev(l,i,j,i)
               he14(j,i,l,i) = he14(l,i,j,i)
               hec(j,i,l,i) = hec(l,i,j,i)
               hecd(j,i,l,i) = hecd(l,i,j,i)
               hed(j,i,l,i) = hed(l,i,j,i)
c              hex(j,i,l,i) = hex(l,i,j,i)
               hessn(j,i,l,i) = hessn(l,i,j,i)
     &                             - hex(l,i,j,i) + hex(j,i,l,i)
            end do
         end do
         do k = 1, i-1
            do j = 1, 3
               do l = 1, 3
                  heb(j,i,l,k) = heb(l,k,j,i)
                  hea(j,i,l,k) = hea(l,k,j,i)
                  heba(j,i,l,k) = heba(l,k,j,i)
                  het(j,i,l,k) = het(l,k,j,i)
                  hev(j,i,l,k) = hev(l,k,j,i)
                  he14(j,i,l,k) = he14(l,k,j,i)
                  hec(j,i,l,k) = hec(l,k,j,i)
                  hecd(j,i,l,k) = hecd(l,k,j,i)
                  hed(j,i,l,k) = hed(l,k,j,i)
c                 hex(j,i,l,k) = hex(l,k,j,i)
                  hessn(j,i,l,k) = hessn(l,k,j,i)
     &                                - hex(l,k,j,i) + hex(j,i,l,k)
               end do
            end do
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ebond5  --  bond energy hessian; cart. version  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ebond5" calculates second derivatives (hessian) of the bond
c     compression energy with respect to the cartesian coordinates
c
c
      subroutine ebond5
      implicit none
      include 'sizes.for'
      integer i,j,k,l,n,nbond,ibnd,ibond
      real*8 x,y,z,bk,bl,cstr,dt
      real*8 rik,rik2,rikx,riky,rikz,deddt,d2eddt2
      real*8 de,term,termx,termy,termz,d2e(3,3)
      real*8 heb,hea,heba,het,hev,he14,hec,hecd,hed
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm)
      common /bond  / nbond,ibnd(2,maxbnd),bk(maxbnd),bl(maxbnd),cstr
      common /hessn1/ heb(3,maxhes,3,maxhes),hea(3,maxhes,3,maxhes),
     &                heba(3,maxhes,3,maxhes),het(3,maxhes,3,maxhes),
     &                hev(3,maxhes,3,maxhes),he14(3,maxhes,3,maxhes),
     &                hec(3,maxhes,3,maxhes),hecd(3,maxhes,3,maxhes),
     &                hed(3,maxhes,3,maxhes)
c
c
c     zero out the bond energy hessian matrix elements
c
      do l = 1, n
         do k = 1, 3
            do j = 1, n
               do i = 1, 3
                  heb(i,j,k,l) = 0.0d0
               end do
            end do
         end do
      end do
c
c     compute the hessian elements of the compression energy
c
      do ibond = 1, nbond
         i = ibnd(1,ibond)
         k = ibnd(2,ibond)
         rikx = x(i) - x(k)
         riky = y(i) - y(k)
         rikz = z(i) - z(k)
         rik2 = rikx*rikx + riky*riky + rikz*rikz
         rik = sqrt(rik2)
         dt = rik - bl(ibond)
         d2eddt2 = 143.88d0 * bk(ibond)
         deddt = d2eddt2 * dt
         if (dt .le. 0.2d0) then
            d2eddt2 = d2eddt2 * (1.0d0+3.0d0*cstr*dt)
            deddt = deddt * (1.0d0+1.5d0*cstr*dt)
         end if
c
c     set the chain rule terms for the hessian elements
c
         de = deddt / rik
         term = (d2eddt2-de) / rik2
         termx = term * rikx
         termy = term * riky
         termz = term * rikz
         d2e(1,1) = termx*rikx + de
         d2e(1,2) = termx*riky
         d2e(1,3) = termx*rikz
         d2e(2,1) = d2e(1,2)
         d2e(2,2) = termy*riky + de
         d2e(2,3) = termy*rikz
         d2e(3,1) = d2e(1,3)
         d2e(3,2) = d2e(2,3)
         d2e(3,3) = termz*rikz + de
c
c     increment diagonal and above-diagonal hessian elements
c
         heb(1,i,1,i) = heb(1,i,1,i) + d2e(1,1)
         heb(1,i,2,i) = heb(1,i,2,i) + d2e(1,2)
         heb(2,i,2,i) = heb(2,i,2,i) + d2e(2,2)
         heb(1,i,3,i) = heb(1,i,3,i) + d2e(1,3)
         heb(2,i,3,i) = heb(2,i,3,i) + d2e(2,3)
         heb(3,i,3,i) = heb(3,i,3,i) + d2e(3,3)
         heb(1,k,1,k) = heb(1,k,1,k) + d2e(1,1)
         heb(1,k,2,k) = heb(1,k,2,k) + d2e(1,2)
         heb(2,k,2,k) = heb(2,k,2,k) + d2e(2,2)
         heb(1,k,3,k) = heb(1,k,3,k) + d2e(1,3)
         heb(2,k,3,k) = heb(2,k,3,k) + d2e(2,3)
         heb(3,k,3,k) = heb(3,k,3,k) + d2e(3,3)
         if (i .lt. k) then
            heb(1,i,1,k) = heb(1,i,1,k) - d2e(1,1)
            heb(2,i,1,k) = heb(2,i,1,k) - d2e(2,1)
            heb(3,i,1,k) = heb(3,i,1,k) - d2e(3,1)
            heb(1,i,2,k) = heb(1,i,2,k) - d2e(1,2)
            heb(2,i,2,k) = heb(2,i,2,k) - d2e(2,2)
            heb(3,i,2,k) = heb(3,i,2,k) - d2e(3,2)
            heb(1,i,3,k) = heb(1,i,3,k) - d2e(1,3)
            heb(2,i,3,k) = heb(2,i,3,k) - d2e(2,3)
            heb(3,i,3,k) = heb(3,i,3,k) - d2e(3,3)
         else if (k .lt. i) then
            heb(1,k,1,i) = heb(1,k,1,i) - d2e(1,1)
            heb(2,k,1,i) = heb(2,k,1,i) - d2e(2,1)
            heb(3,k,1,i) = heb(3,k,1,i) - d2e(3,1)
            heb(1,k,2,i) = heb(1,k,2,i) - d2e(1,2)
            heb(2,k,2,i) = heb(2,k,2,i) - d2e(2,2)
            heb(3,k,2,i) = heb(3,k,2,i) - d2e(3,2)
            heb(1,k,3,i) = heb(1,k,3,i) - d2e(1,3)
            heb(2,k,3,i) = heb(2,k,3,i) - d2e(2,3)
            heb(3,k,3,i) = heb(3,k,3,i) - d2e(3,3)
         end if
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine eangle5  --  bend & str-bend hessian; cartesian  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "eangle5" calculates second derivatives (hessian) of the bend
c     and stretch bend energy with respect to the cartesian coordinates
c     by one-sided finite differences based on first derivatives
c
c
      subroutine eangle5
      implicit none
      include 'sizes.for'
      integer i,j,k,l,n
      real*8 x,y,z,old,eps,dea0(3,maxatm),deba0(3,maxatm)
      real*8 deb,dea,deba,det,dev,de14,dec,decd,ded
      real*8 heb,hea,heba,het,hev,he14,hec,hecd,hed
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm)
      common /deriv1/ deb(3,maxatm),dea(3,maxatm),deba(3,maxatm),
     &                det(3,maxatm),dev(3,maxatm),de14(3,maxatm),
     &                dec(3,maxatm),decd(3,maxatm),ded(3,maxatm)
      common /hessn1/ heb(3,maxhes,3,maxhes),hea(3,maxhes,3,maxhes),
     &                heba(3,maxhes,3,maxhes),het(3,maxhes,3,maxhes),
     &                hev(3,maxhes,3,maxhes),he14(3,maxhes,3,maxhes),
     &                hec(3,maxhes,3,maxhes),hecd(3,maxhes,3,maxhes),
     &                hed(3,maxhes,3,maxhes)
c
c
c     zero out the bend and stretch-bend energy hessian elements
c
      do l = 1, n
         do k = 1, 3
            do j = 1, n
               do i = 1, 3
                  hea(i,j,k,l) = 0.0d0
                  heba(i,j,k,l) = 0.0d0
               end do
            end do
         end do
      end do
c
c     find first derivatives for the base structure
c
      eps = 1.0d-8
      call eangle2
      do k = 1, n
         do l = 1, 3
            dea0(l,k) = dea(l,k)
            deba0(l,k) = deba(l,k)
         end do
      end do
c
c     compute numerical bend and stretch-bend Hessian elements
c
      do i = 1, n
         old = x(i)
         x(i) = x(i) + eps
         call eangle2
         x(i) = old
         do k = 1, n
            do l = 1, 3
               hea(l,k,1,i) = (dea(l,k) - dea0(l,k)) / eps
               heba(l,k,1,i) = (deba(l,k) - deba0(l,k)) / eps
            end do
         end do
         old = y(i)
         y(i) = y(i) + eps
         call eangle2
         y(i) = old
         do k = 1, n
            do l = 1, 3
               hea(l,k,2,i) = (dea(l,k) - dea0(l,k)) / eps
               heba(l,k,2,i) = (deba(l,k) - deba0(l,k)) / eps
            end do
         end do
         old = z(i)
         z(i) = z(i) + eps
         call eangle2
         z(i) = old
         do k = 1, n
            do l = 1, 3
               hea(l,k,3,i) = (dea(l,k) - dea0(l,k)) / eps
               heba(l,k,3,i) = (deba(l,k) - deba0(l,k)) / eps
            end do
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine etors5  --  torsional hessian; cart. version  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "etors5" calculates second derivatives (hessian) of the
c     torsional energy with respect to the cartesian coordinates
c
c
      subroutine etors5
      implicit none
      include 'sizes.for'
      integer i,j,k,l,ia,ib,ic,id,icon,n,ntors,itors,ktb
      real*8 dedcos,d2edcos2,d1d2,sqrtd1d2,cos,v1,v2,v3
      real*8 xia,yia,zia,xib,yib,zib,xic,yic,zic,xid,yid,zid
      real*8 a1,b1,c1,a2,b2,c2,d1,d2,d3,d3d1,d3d2
      real*8 xba,yba,zba,xcb,ycb,zcb,xdc,ydc,zdc
      real*8 xca,yca,zca,xdb,ydb,zdb
      real*8 x,y,z,tcon1,tcon2,tcon3
      real*8 dd1,dd2,dd3,dd3d1,dd3d2,dd1d2
      real*8 cb1yzcb,ac1zxcb,ba1xycb,bc1zyca,cb2yzdc,ca1xzca
      real*8 ac2zxdc,ab1yxca,ba2xydc,cb1yzba,bc2zydb,ac1zxba
      real*8 ca2xzdb,ba1xyba,ab2yxdb,cb2yzcb,ac2zxcb,ba2xycb
      real*8 dcosdxia,dcosdyia,dcosdzia,dcosdxib,dcosdyib,dcosdzib
      real*8 dcosdxic,dcosdyic,dcosdzic,dcosdxid,dcosdyid,dcosdzid
      real*8 dxiaxia,dyiayia,dziazia,dxibxib,dyibyib,dzibzib
      real*8 dxicxic,dyicyic,dziczic,dxidxid,dyidyid,dzidzid
      real*8 dxiayia,dxiazia,dyiazia,dxibyib,dxibzib,dyibzib
      real*8 dxicyic,dxiczic,dyiczic,dxidyid,dxidzid,dyidzid
      real*8 dxiaxib,dxiayib,dxiazib,dyiaxib,dyiayib,dyiazib
      real*8 dziaxib,dziayib,dziazib,dxiaxic,dxiayic,dxiazic
      real*8 dyiaxic,dyiayic,dyiazic,dziaxic,dziayic,dziazic
      real*8 dxiaxid,dxiayid,dxiazid,dyiaxid,dyiayid,dyiazid
      real*8 dziaxid,dziayid,dziazid,dxibxic,dxibyic,dxibzic
      real*8 dyibxic,dyibyic,dyibzic,dzibxic,dzibyic,dzibzic
      real*8 dxibxid,dxibyid,dxibzid,dyibxid,dyibyid,dyibzid
      real*8 dzibxid,dzibyid,dzibzid,dxicxid,dxicyid,dxiczid
      real*8 dyicxid,dyicyid,dyiczid,dzicxid,dzicyid,dziczid
      real*8 heb,hea,heba,het,hev,he14,hec,hecd,hed
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm)
      common /hessn1/ heb(3,maxhes,3,maxhes),hea(3,maxhes,3,maxhes),
     &                heba(3,maxhes,3,maxhes),het(3,maxhes,3,maxhes),
     &                hev(3,maxhes,3,maxhes),he14(3,maxhes,3,maxhes),
     &                hec(3,maxhes,3,maxhes),hecd(3,maxhes,3,maxhes),
     &                hed(3,maxhes,3,maxhes)
      common /tors  / ntors,itors(4,maxtors),ktb(maxtors),
     &                tcon1(maxtcon),tcon2(maxtcon),tcon3(maxtcon)
c
c
c     zero out the torsional energy hessian elements
c
      do l = 1, n
         do k = 1, 3
            do j = 1, n
               do i = 1, 3
                  het(i,j,k,l) = 0.0d0
               end do
            end do
         end do
      end do
c
c     for each torsional angle we first calculate
c     the cosine of the dihedral between ia-ib-ic-id
c
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
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
         a1 = yba*zcb - ycb*zba
         b1 = xcb*zba - xba*zcb
         c1 = xba*ycb - xcb*yba
         a2 = ycb*zdc - ydc*zcb
         b2 = xdc*zcb - xcb*zdc
         c2 = xcb*ydc - xdc*ycb
         d1 = a1*a1 + b1*b1 + c1*c1
         d2 = a2*a2 + b2*b2 + c2*c2
         d3 = a1*a2 + b1*b2 + c1*c2
         d1d2 = d1 * d2
         sqrtd1d2 = sqrt(d1d2)
         cos = d3 / sqrtd1d2
c
c     zero out the chain rule term for this angle, and
c     set the appropriate torsional constants
c
         dedcos = 0.0d0
         d2edcos2 = 0.0d0
         icon = ktb(i)
         v1 = tcon1(icon)
         v2 = tcon2(icon)
         v3 = tcon3(icon)
c
c     calculate the chain rule term for this dihedral angle
c
         if (v1 .ne. 0.0d0) then
            dedcos = 0.5d0*v1
         end if
         if (v2 .ne. 0.0d0) then
            dedcos = dedcos - 2.0d0*v2*cos
            d2edcos2 = -2.0d0 * v2
         end if
         if (v3 .ne. 0.0d0) then
            dedcos = dedcos + 3.0d0*v3*(2.0d0*cos*cos-0.5d0)
            d2edcos2 = d2edcos2 + 12.0d0*v3*cos
         end if
c
c     prepare to increment the torsional energy and derivatives
c
         dedcos = dedcos / sqrtd1d2
         d2edcos2 = d2edcos2 / d1d2
         xca = xic - xia
         yca = yic - yia
         zca = zic - zia
         xdb = xid - xib
         ydb = yid - yib
         zdb = zid - zib
         d3d1 = d3 / d1
         d3d2 = d3 / d2
c
c     abbreviations for the chain rule terms
c
         cb1yzcb = c1*ycb - b1*zcb
         ac1zxcb = a1*zcb - c1*xcb
         ba1xycb = b1*xcb - a1*ycb
         bc1zyca = b1*zca - c1*yca
         cb2yzdc = c2*ydc - b2*zdc
         ca1xzca = c1*xca - a1*zca
         ac2zxdc = a2*zdc - c2*xdc
         ab1yxca = a1*yca - b1*xca
         ba2xydc = b2*xdc - a2*ydc
         cb1yzba = c1*yba - b1*zba
         bc2zydb = b2*zdb - c2*ydb
         ac1zxba = a1*zba - c1*xba
         ca2xzdb = c2*xdb - a2*zdb
         ba1xyba = b1*xba - a1*yba
         ab2yxdb = a2*ydb - b2*xdb
         cb2yzcb = c2*ycb - b2*zcb
         ac2zxcb = a2*zcb - c2*xcb
         ba2xycb = b2*xcb - a2*ycb
c
c     chain rule terms for the first derivatives
c
         dcosdxia = b2*zcb - c2*ycb + cb1yzcb*d3d1
         dcosdyia = c2*xcb - a2*zcb + ac1zxcb*d3d1
         dcosdzia = a2*ycb - b2*xcb + ba1xycb*d3d1
         dcosdxib = c2*yca - b2*zca + bc1zyca*d3d1 +
     &              b1*zdc - c1*ydc + cb2yzdc*d3d2
         dcosdyib = a2*zca - c2*xca + ca1xzca*d3d1 +
     &              c1*xdc - a1*zdc + ac2zxdc*d3d2
         dcosdzib = b2*xca - a2*yca + ab1yxca*d3d1 +
     &              a1*ydc - b1*xdc + ba2xydc*d3d2
         dcosdxic = b2*zba - c2*yba + cb1yzba*d3d1 +
     &              c1*ydb - b1*zdb + bc2zydb*d3d2
         dcosdyic = c2*xba - a2*zba + ac1zxba*d3d1 +
     &              a1*zdb - c1*xdb + ca2xzdb*d3d2
         dcosdzic = a2*yba - b2*xba + ba1xyba*d3d1 +
     &              b1*xdb - a1*ydb + ab2yxdb*d3d2
         dcosdxid = b1*zcb - c1*ycb + cb2yzcb*d3d2
         dcosdyid = c1*xcb - a1*zcb + ac2zxcb*d3d2
         dcosdzid = a1*ycb - b1*xcb + ba2xycb*d3d2
c
c     chain rule terms for second derivative components
c
         dd1 = -2.0d0 * cb1yzcb
         dd3 = -cb2yzcb
         dd3d1 = (dd3-d3d1*dd1)/d1
         dd3d2 = dd3/d2
         dd1d2 = 0.5d0 * (dd1*d2)/d1d2
         dxiaxia = -(ycb*ycb+zcb*zcb)*d3d1 + cb1yzcb*dd3d1
     &                - dcosdxia*dd1d2
         dxiayia = (ycb*xcb)*d3d1 + ac1zxcb*dd3d1 - dcosdyia*dd1d2
         dxiazia = (zcb*xcb)*d3d1 + ba1xycb*dd3d1 - dcosdzia*dd1d2
         dxiaxib = zcb*zdc + ycb*ydc + bc1zyca*dd3d1 + cb2yzdc*dd3d2
     &                + (zca*zcb+yca*ycb)*d3d1 - dcosdxib*dd1d2
         dxiayib = c2 - xdc*ycb + ca1xzca*dd3d1 + ac2zxdc*dd3d2
     &                - (c1+xca*ycb)*d3d1 - dcosdyib*dd1d2
         dxiazib = -b2 - xdc*zcb + ab1yxca*dd3d1 + ba2xydc*dd3d2
     &                + (b1-xca*zcb)*d3d1 - dcosdzib*dd1d2
         dxiaxic = -ycb*ydb - zcb*zdb + cb1yzba*dd3d1 + bc2zydb*dd3d2
     &                - (yba*ycb+zba*zcb)*d3d1 - dcosdxic*dd1d2
         dxiayic = -c2 + xdb*ycb + ac1zxba*dd3d1 + ca2xzdb*dd3d2
     &                + (c1+xba*ycb)*d3d1 - dcosdyic*dd1d2
         dxiazic = b2 + xdb*zcb + ba1xyba*dd3d1 + ab2yxdb*dd3d2
     &                + (xba*zcb-b1)*d3d1 - dcosdzic*dd1d2
         dxiaxid = zcb*zcb + ycb*ycb + cb2yzcb*dd3d2 - dcosdxid*dd1d2
         dxiayid = -xcb*ycb + ac2zxcb*dd3d2 - dcosdyid*dd1d2
         dxiazid = -xcb*zcb + ba2xycb*dd3d2 - dcosdzid*dd1d2
c
         dd1 = -2.0d0 * ac1zxcb
         dd3 = -ac2zxcb
         dd3d1 = (dd3-d3d1*dd1)/d1
         dd3d2 = dd3/d2
         dd1d2 = 0.5d0 * (dd1*d2)/d1d2
         dyiayia = -(zcb*zcb+xcb*xcb)*d3d1 + ac1zxcb*dd3d1
     &                - dcosdyia*dd1d2
         dyiazia = (zcb*ycb)*d3d1 + ba1xycb*dd3d1 - dcosdzia*dd1d2
         dyiaxib = -c2 - ydc*xcb + bc1zyca*dd3d1 + cb2yzdc*dd3d2
     &                + (c1-yca*xcb)*d3d1 - dcosdxib*dd1d2
         dyiayib = xcb*xdc + zcb*zdc + ca1xzca*dd3d1 + ac2zxdc*dd3d2
     &                + (xca*xcb+zca*zcb)*d3d1 - dcosdyib*dd1d2
         dyiazib = a2 - ydc*zcb + ab1yxca*dd3d1 + ba2xydc*dd3d2
     &                - (a1+yca*zcb)*d3d1 - dcosdzib*dd1d2
         dyiaxic = c2 + ydb*xcb + cb1yzba*dd3d1 + bc2zydb*dd3d2
     &                + (yba*xcb-c1)*d3d1 - dcosdxic*dd1d2
         dyiayic = -zcb*zdb - xcb*xdb + ac1zxba*dd3d1 + ca2xzdb*dd3d2
     &                - (zba*zcb+xba*xcb)*d3d1 - dcosdyic*dd1d2
         dyiazic = -a2 + ydb*zcb + ba1xyba*dd3d1 + ab2yxdb*dd3d2
     &                + (a1+yba*zcb)*d3d1 - dcosdzic*dd1d2
         dyiaxid = -ycb*xcb + cb2yzcb*dd3d2 - dcosdxid*dd1d2
         dyiayid = xcb*xcb + zcb*zcb + ac2zxcb*dd3d2 - dcosdyid*dd1d2
         dyiazid = -ycb*zcb + ba2xycb*dd3d2 - dcosdzid*dd1d2
c
         dd1 = -2.0d0 * ba1xycb
         dd3 = -ba2xycb
         dd3d1 = (dd3-d3d1*dd1)/d1
         dd3d2 = dd3/d2
         dd1d2 = 0.5d0 * (dd1*d2)/d1d2
         dziazia = -(xcb*xcb+ycb*ycb)*d3d1 + ba1xycb*dd3d1
     &                - dcosdzia*dd1d2
         dziaxib = b2 - zdc*xcb + bc1zyca*dd3d1 + cb2yzdc*dd3d2
     &                - (b1+zca*xcb)*d3d1 - dcosdxib*dd1d2
         dziayib = -a2 - zdc*ycb + ca1xzca*dd3d1 + ac2zxdc*dd3d2
     &                + (a1-zca*ycb)*d3d1 - dcosdyib*dd1d2
         dziazib = ycb*ydc + xcb*xdc + ab1yxca*dd3d1 + ba2xydc*dd3d2
     &                + (yca*ycb+xca*xcb)*d3d1 - dcosdzib*dd1d2
         dziaxic = -b2 + zdb*xcb + cb1yzba*dd3d1 + bc2zydb*dd3d2
     &                + (b1+zba*xcb)*d3d1 - dcosdxic*dd1d2
         dziayic = a2 + zdb*ycb + ac1zxba*dd3d1 + ca2xzdb*dd3d2
     &                + (zba*ycb-a1)*d3d1 - dcosdyic*dd1d2
         dziazic = -xcb*xdb - ycb*ydb + ba1xyba*dd3d1 + ab2yxdb*dd3d2
     &                - (xba*xcb+yba*ycb)*d3d1 - dcosdzic*dd1d2
         dziaxid = -zcb*xcb + cb2yzcb*dd3d2 - dcosdxid*dd1d2
         dziayid = -zcb*ycb + ac2zxcb*dd3d2 - dcosdyid*dd1d2
         dziazid = ycb*ycb + xcb*xcb + ba2xycb*dd3d2 - dcosdzid*dd1d2
c
         dd1 = 2.0d0 * (cb1yzba+cb1yzcb)
         dd2 = -2.0d0 * cb2yzdc
         dd3 = cb2yzcb + b1*zdc - b2*zba - c1*ydc + c2*yba
         dd3d1 = (dd3-d3d1*dd1)/d1
         dd3d2 = (dd3-d3d2*dd2)/d2
         dd1d2 = 0.5d0 * (d1*dd2+dd1*d2)/d1d2
         dxibxib = -ydc*(yba+yca+ycb) - zdc*(zba+zca+zcb)
     &                - (zca*(zba+zcb)+yca*(yba+ycb))*d3d1
     &                - (ydc*ydc+zdc*zdc)*d3d2 + bc1zyca*dd3d1
     &                + cb2yzdc*dd3d2 - dcosdxib*dd1d2
         dxibyib = xca*ydc + xdc*(yba+ycb) + (xca*(yba+ycb))*d3d1
     &                + (xdc*ydc)*d3d2 + ca1xzca*dd3d1
     &                + ac2zxdc*dd3d2 - dcosdyib*dd1d2
         dxibzib = xca*zdc + xdc*(zba+zcb) + (xca*(zba+zcb))*d3d1
     &                + (xdc*zdc)*d3d2 + ab1yxca*dd3d1
     &                + ba2xydc*dd3d2 - dcosdzib*dd1d2
         dxibxic = zba*zdc + yba*ydc + ydb*(yba+ycb) + zdb*(zba+zcb)
     &                + (yba*(yba+ycb)+zba*(zba+zcb))*d3d1
     &                + (zdb*zdc+ydb*ydc)*d3d2 + cb1yzba*dd3d1
     &                + bc2zydb*dd3d2 - dcosdxic*dd1d2
         dxibyic = c2 - xba*ydc + c1 - xdb*(yba+ycb)
     &                - (c1+xba*(yba+ycb))*d3d1
     &                - (c2+xdb*ydc)*d3d2 + ac1zxba*dd3d1
     &                + ca2xzdb*dd3d2 - dcosdyic*dd1d2
         dxibzic = -b2 - xba*zdc - b1 - xdb*(zba+zcb)
     &                + (b1-xba*(zba+zcb))*d3d1
     &                + (b2-xdb*zdc)*d3d2 + ba1xyba*dd3d1
     &                + ab2yxdb*dd3d2 - dcosdzic*dd1d2
c
         dd1 = 2.0d0 * (ac1zxba+ac1zxcb)
         dd2 = -2.0d0 * ac2zxdc
         dd3 = ac2zxcb - a1*zdc + a2*zba + c1*xdc - c2*xba
         dd3d1 = (dd3-d3d1*dd1)/d1
         dd3d2 = (dd3-d3d2*dd2)/d2
         dd1d2 = 0.5d0 * (d1*dd2+dd1*d2)/d1d2
         dyibyib = -zdc*(zba+zca+zcb) - xdc*(xba+xca+xcb)
     &                - (xca*(xba+xcb)+zca*(zba+zcb))*d3d1
     &                - (zdc*zdc+xdc*xdc)*d3d2 + ca1xzca*dd3d1
     &                + ac2zxdc*dd3d2 - dcosdyib*dd1d2
         dyibzib = yca*zdc + ydc*(zba+zcb) + (yca*(zba+zcb))*d3d1
     &                + (ydc*zdc)*d3d2 + ab1yxca*dd3d1
     &                + ba2xydc*dd3d2 - dcosdzib*dd1d2
         dyibxic = -c2 - yba*xdc - c1 - ydb*(xba+xcb)
     &                + (c1-yba*(xba+xcb))*d3d1
     &                + (c2-ydb*xdc)*d3d2 + cb1yzba*dd3d1
     &                + bc2zydb*dd3d2 - dcosdxic*dd1d2
         dyibyic = xba*xdc + zba*zdc + zdb*(zba+zcb) + xdb*(xba+xcb)
     &                + (zba*(zba+zcb)+xba*(xba+xcb))*d3d1
     &                + (xdb*xdc+zdb*zdc)*d3d2 + ac1zxba*dd3d1
     &                + ca2xzdb*dd3d2 - dcosdyic*dd1d2
         dyibzic = a2 - yba*zdc + a1 - ydb*(zba+zcb)
     &                - (a1+yba*(zba+zcb))*d3d1
     &                - (a2+ydb*zdc)*d3d2 + ba1xyba*dd3d1
     &                + ab2yxdb*dd3d2 - dcosdzic*dd1d2
c
         dd1 = 2.0d0 * (ba1xyba+ba1xycb)
         dd2 = -2.0d0 * ba2xydc
         dd3 = ba2xycb + a1*ydc - a2*yba - b1*xdc + b2*xba
         dd3d1 = (dd3-d3d1*dd1)/d1
         dd3d2 = (dd3-d3d2*dd2)/d2
         dd1d2 = 0.5d0 * (d1*dd2+dd1*d2)/d1d2
         dzibzib = -xdc*(xba+xca+xcb) - ydc*(yba+yca+ycb)
     &                - (yca*(yba+ycb)+xca*(xba+xcb))*d3d1
     &                - (xdc*xdc+ydc*ydc)*d3d2 + ab1yxca*dd3d1
     &                + ba2xydc*dd3d2 - dcosdzib*dd1d2
         dzibxic = b2 - zba*xdc + b1 - zdb*(xba+xcb)
     &                - (b1+zba*(xba+xcb))*d3d1
     &                - (b2+zdb*xdc)*d3d2 + cb1yzba*dd3d1
     &                + bc2zydb*dd3d2 - dcosdxic*dd1d2
         dzibyic = -a2 - zba*ydc - a1 - zdb*(yba+ycb)
     &                + (a1-zba*(yba+ycb))*d3d1
     &                + (a2-zdb*ydc)*d3d2 + ac1zxba*dd3d1
     &                + ca2xzdb*dd3d2 - dcosdyic*dd1d2
         dzibzic = yba*ydc + xba*xdc + xdb*(xba+xcb) + ydb*(yba+ycb)
     &                + (xba*(xba+xcb)+yba*(yba+ycb))*d3d1
     &                + (ydb*ydc+xdb*xdc)*d3d2 + ba1xyba*dd3d1
     &                + ab2yxdb*dd3d2 - dcosdzic*dd1d2
c
         dd1 = -2.0d0 * cb1yzba
         dd2 = 2.0d0 * (cb2yzcb+cb2yzdc)
         dd3 = cb1yzcb -b1*zdc + b2*zba + c1*ydc - c2*yba
         dd3d1 = (dd3-d3d1*dd1)/d1
         dd3d2 = (dd3-d3d2*dd2)/d2
         dd1d2 = 0.5d0 * (d1*dd2+dd1*d2)/d1d2
         dxicxic = - zba*(zcb+zdb+zdc) - yba*(ycb+ydb+ydc)
     &                - (zdb*(zcb+zdc)+ydb*(ycb+ydc))*d3d2
     &                - (yba*yba+zba*zba)*d3d1 + cb1yzba*dd3d1
     &                + bc2zydb*dd3d2 - dcosdxic*dd1d2
         dxicyic = xdb*yba + xba*(ycb+ydc) + (xba*yba)*d3d1
     &                + (xdb*(ycb+ydc))*d3d2 + ac1zxba*dd3d1
     &                + ca2xzdb*dd3d2 - dcosdyic*dd1d2
         dxiczic = xdb*zba + xba*(zcb+zdc) + (xba*zba)*d3d1
     &                + (xdb*(zcb+zdc))*d3d2 + ba1xyba*dd3d1
     &                + ab2yxdb*dd3d2 - dcosdzic*dd1d2
c
         dd1 = -2.0d0 * ac1zxba
         dd2 = 2.0d0 * (ac2zxcb+ac2zxdc)
         dd3 = ac1zxcb + a1*zdc - a2*zba - c1*xdc + c2*xba
         dd3d1 = (dd3-d3d1*dd1)/d1
         dd3d2 = (dd3-d3d2*dd2)/d2
         dd1d2 = 0.5d0 * (d1*dd2+dd1*d2)/d1d2
         dyicyic = -xba*(xcb+xdb+xdc) - zba*(zcb+zdb+zdc)
     &                - (xdb*(xcb+xdc)+zdb*(zcb+zdc))*d3d2
     &                - (zba*zba+xba*xba)*d3d1 + ac1zxba*dd3d1
     &                + ca2xzdb*dd3d2 - dcosdyic*dd1d2
         dyiczic = ydb*zba + yba*(zcb+zdc) + (yba*zba)*d3d1
     &                + (ydb*(zcb+zdc))*d3d2 + ba1xyba*dd3d1
     &                + ab2yxdb*dd3d2 - dcosdzic*dd1d2
c
         dd1 = -2.0d0 * ba1xyba
         dd2 = 2.0d0 * (ba2xycb+ba2xydc)
         dd3 = ba1xycb - a1*ydc + a2*yba + b1*xdc - b2*xba
         dd3d1 = (dd3-d3d1*dd1)/d1
         dd3d2 = (dd3-d3d2*dd2)/d2
         dd1d2 = 0.5d0 * (d1*dd2+dd1*d2)/d1d2
         dziczic = -yba*(ycb+ydb+ydc) - xba*(xcb+xdb+xdc)
     &                - (ydb*(ycb+ydc)+xdb*(xcb+xdc))*d3d2
     &                - (xba*xba+yba*yba)*d3d1 + ba1xyba*dd3d1
     &                + ab2yxdb*dd3d2 - dcosdzic*dd1d2
c
         dd2 = -2.0d0 * cb2yzcb
         dd3 = -cb1yzcb
         dd3d1 = dd3/d1
         dd3d2 = (dd3-d3d2*dd2)/d2
         dd1d2 = 0.5d0 * (d1*dd2)/d1d2
         dxibxid = -yca*ycb - zca*zcb + bc1zyca*dd3d1 + cb2yzdc*dd3d2
     &                - (ycb*ydc+zcb*zdc)*d3d2 - dcosdxib*dd1d2
         dyibxid = c1 + xca*ycb + ca1xzca*dd3d1 + ac2zxdc*dd3d2
     &                + (xdc*ycb-c2)*d3d2 - dcosdyib*dd1d2
         dzibxid = -b1 + xca*zcb + ab1yxca*dd3d1 + ba2xydc*dd3d2
     &                + (b2+xdc*zcb)*d3d2 - dcosdzib*dd1d2
         dxicxid = zba*zcb + yba*ycb + cb1yzba*dd3d1 + bc2zydb*dd3d2
     &                + (zcb*zdb+ycb*ydb)*d3d2 - dcosdxic*dd1d2
         dyicxid = -c1 - xba*ycb + ac1zxba*dd3d1 + ca2xzdb*dd3d2
     &                + (c2-xdb*ycb)*d3d2 - dcosdyic*dd1d2
         dzicxid = b1 - xba*zcb + ba1xyba*dd3d1 + ab2yxdb*dd3d2
     &                - (b2+xdb*zcb)*d3d2 - dcosdzic*dd1d2
         dxidxid = -(ycb*ycb+zcb*zcb)*d3d2 + cb2yzcb*dd3d2
     &                - dcosdxid*dd1d2
         dxidyid = (ycb*xcb)*d3d2 + ac2zxcb*dd3d2 - dcosdyid*dd1d2
         dxidzid = (zcb*xcb)*d3d2 + ba2xycb*dd3d2 - dcosdzid*dd1d2
c
         dd2 = -2.0d0 * ac2zxcb
         dd3 = -ac1zxcb
         dd3d1 = dd3/d1
         dd3d2 = (dd3-d3d2*dd2)/d2
         dd1d2 = 0.5d0 * (d1*dd2)/d1d2
         dxibyid = -c1 + yca*xcb + bc1zyca*dd3d1 + cb2yzdc*dd3d2
     &                + (c2+ydc*xcb)*d3d2 - dcosdxib*dd1d2
         dyibyid = -zca*zcb - xca*xcb + ca1xzca*dd3d1 + ac2zxdc*dd3d2
     &                - (zcb*zdc+xcb*xdc)*d3d2 - dcosdyib*dd1d2
         dzibyid = a1 + yca*zcb + ab1yxca*dd3d1 + ba2xydc*dd3d2
     &                + (ydc*zcb-a2)*d3d2 - dcosdzib*dd1d2
         dxicyid = c1 - yba*xcb + cb1yzba*dd3d1 + bc2zydb*dd3d2
     &                - (c2+ydb*xcb)*d3d2 - dcosdxic*dd1d2
         dyicyid = xba*xcb + zba*zcb + ac1zxba*dd3d1 + ca2xzdb*dd3d2
     &                + (xcb*xdb+zcb*zdb)*d3d2 - dcosdyic*dd1d2
         dzicyid = -a1 - yba*zcb + ba1xyba*dd3d1 + ab2yxdb*dd3d2
     &                + (a2-ydb*zcb)*d3d2 - dcosdzic*dd1d2
         dyidyid = -(zcb*zcb+xcb*xcb)*d3d2 + ac2zxcb*dd3d2
     &                - dcosdyid*dd1d2
         dyidzid = (zcb*ycb)*d3d2 + ba2xycb*dd3d2 - dcosdzid*dd1d2
c
         dd2 = -2.0d0 * ba2xycb
         dd3 = -ba1xycb
         dd3d1 = dd3/d1
         dd3d2 = (dd3-d3d2*dd2)/d2
         dd1d2 = 0.5d0 * (d1*dd2)/d1d2
         dxibzid = b1 + zca*xcb + bc1zyca*dd3d1 + cb2yzdc*dd3d2
     &                + (zdc*xcb-b2)*d3d2 - dcosdxib*dd1d2
         dyibzid = -a1 + zca*ycb + ca1xzca*dd3d1 + ac2zxdc*dd3d2
     &                + (a2+zdc*ycb)*d3d2 - dcosdyib*dd1d2
         dzibzid = -xca*xcb - yca*ycb + ab1yxca*dd3d1 + ba2xydc*dd3d2
     &                - (xcb*xdc+ycb*ydc)*d3d2 - dcosdzib*dd1d2
         dxiczid = -b1 - zba*xcb + cb1yzba*dd3d1 + bc2zydb*dd3d2
     &                + (b2-zdb*xcb)*d3d2 - dcosdxic*dd1d2
         dyiczid = a1 - zba*ycb + ac1zxba*dd3d1 + ca2xzdb*dd3d2
     &                - (a2+zdb*ycb)*d3d2 - dcosdyic*dd1d2
         dziczid = yba*ycb + xba*xcb + ba1xyba*dd3d1 + ab2yxdb*dd3d2
     &                + (ycb*ydb+xcb*xdb)*d3d2 - dcosdzic*dd1d2
         dzidzid = -(xcb*xcb+ycb*ycb)*d3d2 + ba2xycb*dd3d2
     &                - dcosdzid*dd1d2
c
c     now, increment the diagonal hessian block elements
c
         het(1,ia,1,ia) = het(1,ia,1,ia) + dedcos*dxiaxia
     &                     + d2edcos2*dcosdxia*dcosdxia
         het(1,ia,2,ia) = het(1,ia,2,ia) + dedcos*dxiayia
     &                     + d2edcos2*dcosdxia*dcosdyia
         het(2,ia,2,ia) = het(2,ia,2,ia) + dedcos*dyiayia
     &                     + d2edcos2*dcosdyia*dcosdyia
         het(1,ia,3,ia) = het(1,ia,3,ia) + dedcos*dxiazia
     &                     + d2edcos2*dcosdxia*dcosdzia
         het(2,ia,3,ia) = het(2,ia,3,ia) + dedcos*dyiazia
     &                     + d2edcos2*dcosdyia*dcosdzia
         het(3,ia,3,ia) = het(3,ia,3,ia) + dedcos*dziazia
     &                     + d2edcos2*dcosdzia*dcosdzia
         het(1,ib,1,ib) = het(1,ib,1,ib) + dedcos*dxibxib
     &                     + d2edcos2*dcosdxib*dcosdxib
         het(1,ib,2,ib) = het(1,ib,2,ib) + dedcos*dxibyib
     &                     + d2edcos2*dcosdxib*dcosdyib
         het(2,ib,2,ib) = het(2,ib,2,ib) + dedcos*dyibyib
     &                     + d2edcos2*dcosdyib*dcosdyib
         het(1,ib,3,ib) = het(1,ib,3,ib) + dedcos*dxibzib
     &                     + d2edcos2*dcosdxib*dcosdzib
         het(2,ib,3,ib) = het(2,ib,3,ib) + dedcos*dyibzib
     &                     + d2edcos2*dcosdyib*dcosdzib
         het(3,ib,3,ib) = het(3,ib,3,ib) + dedcos*dzibzib
     &                     + d2edcos2*dcosdzib*dcosdzib
         het(1,ic,1,ic) = het(1,ic,1,ic) + dedcos*dxicxic
     &                     + d2edcos2*dcosdxic*dcosdxic
         het(1,ic,2,ic) = het(1,ic,2,ic) + dedcos*dxicyic
     &                     + d2edcos2*dcosdxic*dcosdyic
         het(2,ic,2,ic) = het(2,ic,2,ic) + dedcos*dyicyic
     &                     + d2edcos2*dcosdyic*dcosdyic
         het(1,ic,3,ic) = het(1,ic,3,ic) + dedcos*dxiczic
     &                     + d2edcos2*dcosdxic*dcosdzic
         het(2,ic,3,ic) = het(2,ic,3,ic) + dedcos*dyiczic
     &                     + d2edcos2*dcosdyic*dcosdzic
         het(3,ic,3,ic) = het(3,ic,3,ic) + dedcos*dziczic
     &                     + d2edcos2*dcosdzic*dcosdzic
         het(1,id,1,id) = het(1,id,1,id) + dedcos*dxidxid
     &                     + d2edcos2*dcosdxid*dcosdxid
         het(1,id,2,id) = het(1,id,2,id) + dedcos*dxidyid
     &                     + d2edcos2*dcosdxid*dcosdyid
         het(2,id,2,id) = het(2,id,2,id) + dedcos*dyidyid
     &                     + d2edcos2*dcosdyid*dcosdyid
         het(1,id,3,id) = het(1,id,3,id) + dedcos*dxidzid
     &                     + d2edcos2*dcosdxid*dcosdzid
         het(2,id,3,id) = het(2,id,3,id) + dedcos*dyidzid
     &                     + d2edcos2*dcosdyid*dcosdzid
         het(3,id,3,id) = het(3,id,3,id) + dedcos*dzidzid
     &                     + d2edcos2*dcosdzid*dcosdzid
c
c     finally, increment above-diagonal hessian block elements
c
         if (ia .lt. ib) then
            het(1,ia,1,ib) = het(1,ia,1,ib) + dedcos*dxiaxib
     &                        + d2edcos2*dcosdxia*dcosdxib
            het(2,ia,1,ib) = het(2,ia,1,ib) + dedcos*dyiaxib
     &                        + d2edcos2*dcosdyia*dcosdxib
            het(3,ia,1,ib) = het(3,ia,1,ib) + dedcos*dziaxib
     &                        + d2edcos2*dcosdzia*dcosdxib
            het(1,ia,2,ib) = het(1,ia,2,ib) + dedcos*dxiayib
     &                        + d2edcos2*dcosdxia*dcosdyib
            het(2,ia,2,ib) = het(2,ia,2,ib) + dedcos*dyiayib
     &                        + d2edcos2*dcosdyia*dcosdyib
            het(3,ia,2,ib) = het(3,ia,2,ib) + dedcos*dziayib
     &                        + d2edcos2*dcosdzia*dcosdyib
            het(1,ia,3,ib) = het(1,ia,3,ib) + dedcos*dxiazib
     &                        + d2edcos2*dcosdxia*dcosdzib
            het(2,ia,3,ib) = het(2,ia,3,ib) + dedcos*dyiazib
     &                        + d2edcos2*dcosdyia*dcosdzib
            het(3,ia,3,ib) = het(3,ia,3,ib) + dedcos*dziazib
     &                        + d2edcos2*dcosdzia*dcosdzib
         else
            het(1,ib,1,ia) = het(1,ib,1,ia) + dedcos*dxiaxib
     &                        + d2edcos2*dcosdxib*dcosdxia
            het(2,ib,1,ia) = het(2,ib,1,ia) + dedcos*dxiayib
     &                        + d2edcos2*dcosdyib*dcosdxia
            het(3,ib,1,ia) = het(3,ib,1,ia) + dedcos*dxiazib
     &                        + d2edcos2*dcosdzib*dcosdxia
            het(1,ib,2,ia) = het(1,ib,2,ia) + dedcos*dyiaxib
     &                        + d2edcos2*dcosdxib*dcosdyia
            het(2,ib,2,ia) = het(2,ib,2,ia) + dedcos*dyiayib
     &                        + d2edcos2*dcosdyib*dcosdyia
            het(3,ib,2,ia) = het(3,ib,2,ia) + dedcos*dyiazib
     &                        + d2edcos2*dcosdzib*dcosdyia
            het(1,ib,3,ia) = het(1,ib,3,ia) + dedcos*dziaxib
     &                        + d2edcos2*dcosdxib*dcosdzia
            het(2,ib,3,ia) = het(2,ib,3,ia) + dedcos*dziayib
     &                        + d2edcos2*dcosdyib*dcosdzia
            het(3,ib,3,ia) = het(3,ib,3,ia) + dedcos*dziazib
     &                        + d2edcos2*dcosdzib*dcosdzia
         end if
c
         if (ia .lt. ic) then
            het(1,ia,1,ic) = het(1,ia,1,ic) + dedcos*dxiaxic
     &                        + d2edcos2*dcosdxia*dcosdxic
            het(2,ia,1,ic) = het(2,ia,1,ic) + dedcos*dyiaxic
     &                        + d2edcos2*dcosdyia*dcosdxic
            het(3,ia,1,ic) = het(3,ia,1,ic) + dedcos*dziaxic
     &                        + d2edcos2*dcosdzia*dcosdxic
            het(1,ia,2,ic) = het(1,ia,2,ic) + dedcos*dxiayic
     &                        + d2edcos2*dcosdxia*dcosdyic
            het(2,ia,2,ic) = het(2,ia,2,ic) + dedcos*dyiayic
     &                        + d2edcos2*dcosdyia*dcosdyic
            het(3,ia,2,ic) = het(3,ia,2,ic) + dedcos*dziayic
     &                        + d2edcos2*dcosdzia*dcosdyic
            het(1,ia,3,ic) = het(1,ia,3,ic) + dedcos*dxiazic
     &                        + d2edcos2*dcosdxia*dcosdzic
            het(2,ia,3,ic) = het(2,ia,3,ic) + dedcos*dyiazic
     &                        + d2edcos2*dcosdyia*dcosdzic
            het(3,ia,3,ic) = het(3,ia,3,ic) + dedcos*dziazic
     &                        + d2edcos2*dcosdzia*dcosdzic
         else
            het(1,ic,1,ia) = het(1,ic,1,ia) + dedcos*dxiaxic
     &                        + d2edcos2*dcosdxic*dcosdxia
            het(2,ic,1,ia) = het(2,ic,1,ia) + dedcos*dxiayic
     &                        + d2edcos2*dcosdyic*dcosdxia
            het(3,ic,1,ia) = het(3,ic,1,ia) + dedcos*dxiazic
     &                        + d2edcos2*dcosdzic*dcosdxia
            het(1,ic,2,ia) = het(1,ic,2,ia) + dedcos*dyiaxic
     &                        + d2edcos2*dcosdxic*dcosdyia
            het(2,ic,2,ia) = het(2,ic,2,ia) + dedcos*dyiayic
     &                        + d2edcos2*dcosdyic*dcosdyia
            het(3,ic,2,ia) = het(3,ic,2,ia) + dedcos*dyiazic
     &                        + d2edcos2*dcosdzic*dcosdyia
            het(1,ic,3,ia) = het(1,ic,3,ia) + dedcos*dziaxic
     &                        + d2edcos2*dcosdxic*dcosdzia
            het(2,ic,3,ia) = het(2,ic,3,ia) + dedcos*dziayic
     &                        + d2edcos2*dcosdyic*dcosdzia
            het(3,ic,3,ia) = het(3,ic,3,ia) + dedcos*dziazic
     &                        + d2edcos2*dcosdzic*dcosdzia
         end if
c
         if (ia .lt. id) then
            het(1,ia,1,id) = het(1,ia,1,id) + dedcos*dxiaxid
     &                        + d2edcos2*dcosdxia*dcosdxid
            het(2,ia,1,id) = het(2,ia,1,id) + dedcos*dyiaxid
     &                        + d2edcos2*dcosdyia*dcosdxid
            het(3,ia,1,id) = het(3,ia,1,id) + dedcos*dziaxid
     &                        + d2edcos2*dcosdzia*dcosdxid
            het(1,ia,2,id) = het(1,ia,2,id) + dedcos*dxiayid
     &                        + d2edcos2*dcosdxia*dcosdyid
            het(2,ia,2,id) = het(2,ia,2,id) + dedcos*dyiayid
     &                        + d2edcos2*dcosdyia*dcosdyid
            het(3,ia,2,id) = het(3,ia,2,id) + dedcos*dziayid
     &                        + d2edcos2*dcosdzia*dcosdyid
            het(1,ia,3,id) = het(1,ia,3,id) + dedcos*dxiazid
     &                        + d2edcos2*dcosdxia*dcosdzid
            het(2,ia,3,id) = het(2,ia,3,id) + dedcos*dyiazid
     &                        + d2edcos2*dcosdyia*dcosdzid
            het(3,ia,3,id) = het(3,ia,3,id) + dedcos*dziazid
     &                        + d2edcos2*dcosdzia*dcosdzid
         else
            het(1,id,1,ia) = het(1,id,1,ia) + dedcos*dxiaxid
     &                        + d2edcos2*dcosdxid*dcosdxia
            het(2,id,1,ia) = het(2,id,1,ia) + dedcos*dxiayid
     &                        + d2edcos2*dcosdyid*dcosdxia
            het(3,id,1,ia) = het(3,id,1,ia) + dedcos*dxiazid
     &                        + d2edcos2*dcosdzid*dcosdxia
            het(1,id,2,ia) = het(1,id,2,ia) + dedcos*dyiaxid
     &                        + d2edcos2*dcosdxid*dcosdyia
            het(2,id,2,ia) = het(2,id,2,ia) + dedcos*dyiayid
     &                        + d2edcos2*dcosdyid*dcosdyia
            het(3,id,2,ia) = het(3,id,2,ia) + dedcos*dyiazid
     &                        + d2edcos2*dcosdzid*dcosdyia
            het(1,id,3,ia) = het(1,id,3,ia) + dedcos*dziaxid
     &                        + d2edcos2*dcosdxid*dcosdzia
            het(2,id,3,ia) = het(2,id,3,ia) + dedcos*dziayid
     &                        + d2edcos2*dcosdyid*dcosdzia
            het(3,id,3,ia) = het(3,id,3,ia) + dedcos*dziazid
     &                        + d2edcos2*dcosdzid*dcosdzia
         end if
c
         if (ib .lt. ic) then
            het(1,ib,1,ic) = het(1,ib,1,ic) + dedcos*dxibxic
     &                        + d2edcos2*dcosdxib*dcosdxic
            het(2,ib,1,ic) = het(2,ib,1,ic) + dedcos*dyibxic
     &                        + d2edcos2*dcosdyib*dcosdxic
            het(3,ib,1,ic) = het(3,ib,1,ic) + dedcos*dzibxic
     &                        + d2edcos2*dcosdzib*dcosdxic
            het(1,ib,2,ic) = het(1,ib,2,ic) + dedcos*dxibyic
     &                        + d2edcos2*dcosdxib*dcosdyic
            het(2,ib,2,ic) = het(2,ib,2,ic) + dedcos*dyibyic
     &                        + d2edcos2*dcosdyib*dcosdyic
            het(3,ib,2,ic) = het(3,ib,2,ic) + dedcos*dzibyic
     &                        + d2edcos2*dcosdzib*dcosdyic
            het(1,ib,3,ic) = het(1,ib,3,ic) + dedcos*dxibzic
     &                        + d2edcos2*dcosdxib*dcosdzic
            het(2,ib,3,ic) = het(2,ib,3,ic) + dedcos*dyibzic
     &                        + d2edcos2*dcosdyib*dcosdzic
            het(3,ib,3,ic) = het(3,ib,3,ic) + dedcos*dzibzic
     &                        + d2edcos2*dcosdzib*dcosdzic
         else
            het(1,ic,1,ib) = het(1,ic,1,ib) + dedcos*dxibxic
     &                        + d2edcos2*dcosdxic*dcosdxib
            het(2,ic,1,ib) = het(2,ic,1,ib) + dedcos*dxibyic
     &                        + d2edcos2*dcosdyic*dcosdxib
            het(3,ic,1,ib) = het(3,ic,1,ib) + dedcos*dxibzic
     &                        + d2edcos2*dcosdzic*dcosdxib
            het(1,ic,2,ib) = het(1,ic,2,ib) + dedcos*dyibxic
     &                        + d2edcos2*dcosdxic*dcosdyib
            het(2,ic,2,ib) = het(2,ic,2,ib) + dedcos*dyibyic
     &                        + d2edcos2*dcosdyic*dcosdyib
            het(3,ic,2,ib) = het(3,ic,2,ib) + dedcos*dyibzic
     &                        + d2edcos2*dcosdzic*dcosdyib
            het(1,ic,3,ib) = het(1,ic,3,ib) + dedcos*dzibxic
     &                        + d2edcos2*dcosdxic*dcosdzib
            het(2,ic,3,ib) = het(2,ic,3,ib) + dedcos*dzibyic
     &                        + d2edcos2*dcosdyic*dcosdzib
            het(3,ic,3,ib) = het(3,ic,3,ib) + dedcos*dzibzic
     &                        + d2edcos2*dcosdzic*dcosdzib
         end if
c
         if (ib .lt. id) then
            het(1,ib,1,id) = het(1,ib,1,id) + dedcos*dxibxid
     &                        + d2edcos2*dcosdxib*dcosdxid
            het(2,ib,1,id) = het(2,ib,1,id) + dedcos*dyibxid
     &                        + d2edcos2*dcosdyib*dcosdxid
            het(3,ib,1,id) = het(3,ib,1,id) + dedcos*dzibxid
     &                        + d2edcos2*dcosdzib*dcosdxid
            het(1,ib,2,id) = het(1,ib,2,id) + dedcos*dxibyid
     &                        + d2edcos2*dcosdxib*dcosdyid
            het(2,ib,2,id) = het(2,ib,2,id) + dedcos*dyibyid
     &                        + d2edcos2*dcosdyib*dcosdyid
            het(3,ib,2,id) = het(3,ib,2,id) + dedcos*dzibyid
     &                        + d2edcos2*dcosdzib*dcosdyid
            het(1,ib,3,id) = het(1,ib,3,id) + dedcos*dxibzid
     &                        + d2edcos2*dcosdxib*dcosdzid
            het(2,ib,3,id) = het(2,ib,3,id) + dedcos*dyibzid
     &                        + d2edcos2*dcosdyib*dcosdzid
            het(3,ib,3,id) = het(3,ib,3,id) + dedcos*dzibzid
     &                        + d2edcos2*dcosdzib*dcosdzid
         else
            het(1,id,1,ib) = het(1,id,1,ib) + dedcos*dxibxid
     &                        + d2edcos2*dcosdxid*dcosdxib
            het(2,id,1,ib) = het(2,id,1,ib) + dedcos*dxibyid
     &                        + d2edcos2*dcosdyid*dcosdxib
            het(3,id,1,ib) = het(3,id,1,ib) + dedcos*dxibzid
     &                        + d2edcos2*dcosdzid*dcosdxib
            het(1,id,2,ib) = het(1,id,2,ib) + dedcos*dyibxid
     &                        + d2edcos2*dcosdxid*dcosdyib
            het(2,id,2,ib) = het(2,id,2,ib) + dedcos*dyibyid
     &                        + d2edcos2*dcosdyid*dcosdyib
            het(3,id,2,ib) = het(3,id,2,ib) + dedcos*dyibzid
     &                        + d2edcos2*dcosdzid*dcosdyib
            het(1,id,3,ib) = het(1,id,3,ib) + dedcos*dzibxid
     &                        + d2edcos2*dcosdxid*dcosdzib
            het(2,id,3,ib) = het(2,id,3,ib) + dedcos*dzibyid
     &                        + d2edcos2*dcosdyid*dcosdzib
            het(3,id,3,ib) = het(3,id,3,ib) + dedcos*dzibzid
     &                        + d2edcos2*dcosdzid*dcosdzib
         end if
c
         if (ic .lt. id) then
            het(1,ic,1,id) = het(1,ic,1,id) + dedcos*dxicxid
     &                        + d2edcos2*dcosdxic*dcosdxid
            het(2,ic,1,id) = het(2,ic,1,id) + dedcos*dyicxid
     &                        + d2edcos2*dcosdyic*dcosdxid
            het(3,ic,1,id) = het(3,ic,1,id) + dedcos*dzicxid
     &                        + d2edcos2*dcosdzic*dcosdxid
            het(1,ic,2,id) = het(1,ic,2,id) + dedcos*dxicyid
     &                        + d2edcos2*dcosdxic*dcosdyid
            het(2,ic,2,id) = het(2,ic,2,id) + dedcos*dyicyid
     &                        + d2edcos2*dcosdyic*dcosdyid
            het(3,ic,2,id) = het(3,ic,2,id) + dedcos*dzicyid
     &                        + d2edcos2*dcosdzic*dcosdyid
            het(1,ic,3,id) = het(1,ic,3,id) + dedcos*dxiczid
     &                        + d2edcos2*dcosdxic*dcosdzid
            het(2,ic,3,id) = het(2,ic,3,id) + dedcos*dyiczid
     &                        + d2edcos2*dcosdyic*dcosdzid
            het(3,ic,3,id) = het(3,ic,3,id) + dedcos*dziczid
     &                        + d2edcos2*dcosdzic*dcosdzid
         else
            het(1,id,1,ic) = het(1,id,1,ic) + dedcos*dxicxid
     &                        + d2edcos2*dcosdxid*dcosdxic
            het(2,id,1,ic) = het(2,id,1,ic) + dedcos*dxicyid
     &                        + d2edcos2*dcosdyid*dcosdxic
            het(3,id,1,ic) = het(3,id,1,ic) + dedcos*dxiczid
     &                        + d2edcos2*dcosdzid*dcosdxic
            het(1,id,2,ic) = het(1,id,2,ic) + dedcos*dyicxid
     &                        + d2edcos2*dcosdxid*dcosdyic
            het(2,id,2,ic) = het(2,id,2,ic) + dedcos*dyicyid
     &                        + d2edcos2*dcosdyid*dcosdyic
            het(3,id,2,ic) = het(3,id,2,ic) + dedcos*dyiczid
     &                        + d2edcos2*dcosdzid*dcosdyic
            het(1,id,3,ic) = het(1,id,3,ic) + dedcos*dzicxid
     &                        + d2edcos2*dcosdxid*dcosdzic
            het(2,id,3,ic) = het(2,id,3,ic) + dedcos*dzicyid
     &                        + d2edcos2*dcosdyid*dcosdzic
            het(3,id,3,ic) = het(3,id,3,ic) + dedcos*dziczid
     &                        + d2edcos2*dcosdzid*dcosdzic
         end if
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine evdw5  --  van der Waals hessian; cart. version  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "evdw5" calculates second derivatives (hessian) of the van
c     der Waals energy with respect to the cartesian coordinates
c
c
      subroutine evdw5
      implicit none
      include 'sizes.for'
      integer i,j,k,l,it,kt,ii,kk,ich,ikch,chtype,skip(maxatm)
      integer n,itype,nvdw,nv14,n12,i12,n13,i13,n14,i14
      real*8 xr(maxatm),yr(maxatm),zr(maxatm)
      real*8 reduc,vrad,veps,seps,epsch,radch,epscd,radcd,vdwcut
      real*8 p,p2,p6,p7,p8,eps,rv,epsi,vradi,cutoff
      real*8 x,y,z,xi,yi,zi,rikx,riky,rikz,rik,rik2
      real*8 redi,redii,redk,redkk
      real*8 expterm,dedrik,d2edrik2,term
      real*8 termx,termy,termz,d2e(3,3)
      real*8 heb,hea,heba,het,hev,he14,hec,hecd,hed
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm),itype(maxatm)
      common /attach/ n12(maxatm),i12(4,maxatm),n13(maxatm),
     &                i13(12,maxatm),n14(maxatm),i14(36,maxatm)
      common /hessn1/ heb(3,maxhes,3,maxhes),hea(3,maxhes,3,maxhes),
     &                heba(3,maxhes,3,maxhes),het(3,maxhes,3,maxhes),
     &                hev(3,maxhes,3,maxhes),he14(3,maxhes,3,maxhes),
     &                hec(3,maxhes,3,maxhes),hecd(3,maxhes,3,maxhes),
     &                hed(3,maxhes,3,maxhes)
      common /vdw   / nvdw,nv14,vrad(maxtyp),veps(maxtyp),seps(maxtyp),
     &                reduc(maxtyp),chtype(maxtyp),epsch,radch,epscd,
     &                radcd,vdwcut
c
c
c     zero out the van der Waals energy hessian elements
c
      do l = 1, n
         do k = 1, 3
            do j = 1, n
               do i = 1, 3
                  hev(i,j,k,l) = 0.0d0
                  he14(i,j,k,l) = 0.0d0
               end do
            end do
         end do
      end do
      cutoff = vdwcut
c
c     calculate "reduced" atomic coordinates for hydrogens
c
      do i = 1, n
         redi = reduc(itype(i))
         if (redi .eq. 0.0d0) then
            xr(i) = x(i)
            yr(i) = y(i)
            zr(i) = z(i)
         else
            ii = i12(1,i)
            xr(i) = redi*(x(i)-x(ii)) + x(ii)
            yr(i) = redi*(y(i)-y(ii)) + y(ii)
            zr(i) = redi*(z(i)-z(ii)) + z(ii)
         end if
      end do
c
c     find van der Waals hessian elements via double loop
c
      do i = 1, n-1
         do k = 1, n12(i)
            skip(i12(k,i)) = i
         end do
         do k = 1, n13(i)
            skip(i13(k,i)) = i
         end do
         do k = 1, n14(i)
            skip(i14(k,i)) = -i
         end do
         it = itype(i)
         redi = reduc(it)
         if (redi .ne. 0.0d0) then
            ii = i12(1,i)
            redii = 1.0d0 - redi
         end if
         ich = chtype(it)
         epsi = seps(it)
         vradi = vrad(it)
         xi = xr(i)
         yi = yr(i)
         zi = zr(i)
         do k = i+1, n
            if (skip(k) .ne. i) then
               rikx = xi - xr(k)
               riky = yi - yr(k)
               rikz = zi - zr(k)
               rik2 = rikx*rikx + riky*riky + rikz*rikz
               if (rik2 .le. cutoff) then
                  kt = itype(k)
                  redk = reduc(kt)
                  if (redk .ne. 0.0d0) then
                     kk = i12(1,k)
                     redkk = 1.0d0 - redk
                  end if
                  ikch = ich * chtype(kt)
                  if (ikch .ge. 0) then
                     eps = epsi * seps(kt)
                     rv = vradi + vrad(kt)
                  else if (ikch .eq. -1) then
                     eps = epsch
                     rv = radch
                  else if (ikch .eq. -2) then
                     eps = epscd
                     rv = radcd
                  end if
c
c     compute hessian elements for the "i-k" interaction
c
                  p2 = (rv*rv)/rik2
                  if (p2 .le. 0.04d0) then
                     p6 = p2 * p2 * p2
                     dedrik = eps * 13.5d0 * p6/rik2
                     d2edrik2 = eps * -94.5d0 * p6/rik2
                  else if (p2 .le. 10.963d0) then
                     rik = sqrt(rik2)
                     p = rv/rik
                     p7 = p2 * p2 * p2 * p
                     p8 = p7 * p
                     expterm = exp(-12.5d0/p)
                     dedrik = (eps/(rv*rik)) *
     &                  (-3625000.0d0*expterm+13.5d0*p7)
                     d2edrik2 = (eps/(rv*rv)) *
     &                  (45312500.0d0*expterm-94.5d0*p8)
                  else
                     dedrik = eps * -772.352d0 * p2/rik2
                     d2edrik2 = eps * 2317.056d0 * p2/rik2
                  end if
c
c     get chain rule terms for van der Waals hessian elements
c
                  term = (d2edrik2-dedrik) / rik2
                  termx = term * rikx
                  termy = term * riky
                  termz = term * rikz
                  d2e(1,1) = termx*rikx + dedrik
                  d2e(1,2) = termx*riky
                  d2e(1,3) = termx*rikz
                  d2e(2,1) = d2e(1,2)
                  d2e(2,2) = termy*riky + dedrik
                  d2e(2,3) = termy*rikz
                  d2e(3,1) = d2e(1,3)
                  d2e(3,2) = d2e(2,3)
                  d2e(3,3) = termz*rikz + dedrik
c
c     increment diagonal and above-diagonal hessian elements
c
                  if (skip(k) .eq. -i) then
                     if (redi.eq.0.0d0 .and. redk.eq.0.0d0) then
                        do j = 1, 3
                           do l = j, 3
                              he14(j,i,l,i) = he14(j,i,l,i) + d2e(j,l)
                              he14(j,k,l,k) = he14(j,k,l,k) + d2e(j,l)
                           end do
                        end do
                        do j = 1, 3
                           do l = 1, 3
                              he14(j,i,l,k) = he14(j,i,l,k) - d2e(j,l)
                           end do
                        end do
                     else if (redk .eq. 0.0d0) then
                        do j = 1, 3
                           do l = j, 3
                              he14(j,i,l,i) = he14(j,i,l,i) +
     &                           d2e(j,l)*redi*redi
                              he14(j,ii,l,ii) = he14(j,ii,l,ii) +
     &                           d2e(j,l)*redii*redii
                              he14(j,k,l,k) = he14(j,k,l,k) + d2e(j,l)
                           end do
                        end do
                        do j = 1, 3
                           do l = 1, 3
                              he14(j,i,l,k) = he14(j,i,l,k) -
     &                           d2e(j,l)*redi
                           end do
                        end do
                        if (i .lt. ii) then
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,i,l,ii) = he14(j,i,l,ii) +
     &                              d2e(j,l)*redi*redii
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,ii,l,i) = he14(j,ii,l,i) +
     &                              d2e(j,l)*redi*redii
                              end do
                           end do
                        end if
                        if (k .lt. ii) then
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,k,l,ii) = he14(j,k,l,ii) -
     &                              d2e(j,l)*redii
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,ii,l,k) = he14(j,ii,l,k) -
     &                              d2e(j,l)*redii
                              end do
                           end do
                        end if
                     else if (redi .eq. 0.0d0) then
                        do j = 1, 3
                           do l = j, 3
                              he14(j,i,l,i) = he14(j,i,l,i) + d2e(j,l)
                              he14(j,k,l,k) = he14(j,k,l,k) +
     &                           d2e(j,l)*redk*redk
                              he14(j,kk,l,kk) = he14(j,kk,l,kk) +
     &                           d2e(j,l)*redkk*redkk
                           end do
                        end do
                        do j = 1, 3
                           do l = 1, 3
                              he14(j,i,l,k) = he14(j,i,l,k) -
     &                           d2e(j,l)*redk
                           end do
                        end do
                        if (k .lt. kk) then
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,k,l,kk) = he14(j,k,l,kk) +
     &                              d2e(j,l)*redk*redkk
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,kk,l,k) = he14(j,kk,l,k) +
     &                              d2e(j,l)*redk*redkk
                              end do
                           end do
                        end if
                        if (i .lt. kk) then
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,i,l,kk) = he14(j,i,l,kk) -
     &                              d2e(j,l)*redkk
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,kk,l,i) = he14(j,kk,l,i) -
     &                              d2e(j,l)*redkk
                              end do
                           end do
                        end if
                     else
                        do j = 1, 3
                           do l = j, 3
                              he14(j,i,l,i) = he14(j,i,l,i) +
     &                           d2e(j,l)*redi*redi
                              he14(j,ii,l,ii) = he14(j,ii,l,ii) +
     &                           d2e(j,l)*redii*redii
                              he14(j,k,l,k) = he14(j,k,l,k) +
     &                           d2e(j,l)*redk*redk
                              he14(j,kk,l,kk) = he14(j,kk,l,kk) +
     &                           d2e(j,l)*redkk*redkk
                           end do
                        end do
                        do j = 1, 3
                           do l = 1, 3
                              he14(j,i,l,k) = he14(j,i,l,k) -
     &                           d2e(j,l)*redi*redk
                           end do
                        end do
                        if (i .lt. ii) then
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,i,l,ii) = he14(j,i,l,ii) +
     &                              d2e(j,l)*redi*redii
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,ii,l,i) = he14(j,ii,l,i) +
     &                              d2e(j,l)*redi*redii
                              end do
                           end do
                        end if
                        if (k .lt. kk) then
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,k,l,kk) = he14(j,k,l,kk) +
     &                              d2e(j,l)*redk*redkk
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,kk,l,k) = he14(j,kk,l,k) +
     &                              d2e(j,l)*redk*redkk
                              end do
                           end do
                        end if
                        if (i .lt. kk) then
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,i,l,kk) = he14(j,i,l,kk) -
     &                              d2e(j,l)*redi*redkk
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,kk,l,i) = he14(j,kk,l,i) -
     &                              d2e(j,l)*redi*redkk
                              end do
                           end do
                        end if
                        if (k .lt. ii) then
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,k,l,ii) = he14(j,k,l,ii) -
     &                              d2e(j,l)*redk*redii
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,ii,l,k) = he14(j,ii,l,k) -
     &                              d2e(j,l)*redk*redii
                              end do
                           end do
                        end if
                        if (ii .lt. kk) then
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,ii,l,kk) = he14(j,ii,l,kk) -
     &                              d2e(j,l)*redii*redkk
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 he14(j,kk,l,ii) = he14(j,kk,l,ii) -
     &                              d2e(j,l)*redii*redkk
                              end do
                           end do
                        end if
                     end if
                  else
                     if (redi.eq.0.0d0 .and. redk.eq.0.0d0) then
                        do j = 1, 3
                           do l = j, 3
                              hev(j,i,l,i) = hev(j,i,l,i) + d2e(j,l)
                              hev(j,k,l,k) = hev(j,k,l,k) + d2e(j,l)
                           end do
                        end do
                        do j = 1, 3
                           do l = 1, 3
                              hev(j,i,l,k) = hev(j,i,l,k) - d2e(j,l)
                           end do
                        end do
                     else if (redk .eq. 0.0d0) then
                        do j = 1, 3
                           do l = j, 3
                              hev(j,i,l,i) = hev(j,i,l,i) +
     &                           d2e(j,l)*redi*redi
                              hev(j,ii,l,ii) = hev(j,ii,l,ii) +
     &                           d2e(j,l)*redii*redii
                              hev(j,k,l,k) = hev(j,k,l,k) + d2e(j,l)
                           end do
                        end do
                        do j = 1, 3
                           do l = 1, 3
                              hev(j,i,l,k) = hev(j,i,l,k) -
     &                           d2e(j,l)*redi
                           end do
                        end do
                        if (i .lt. ii) then
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,i,l,ii) = hev(j,i,l,ii) +
     &                              d2e(j,l)*redi*redii
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,ii,l,i) = hev(j,ii,l,i) +
     &                              d2e(j,l)*redi*redii
                              end do
                           end do
                        end if
                        if (k .lt. ii) then
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,k,l,ii) = hev(j,k,l,ii) -
     &                              d2e(j,l)*redii
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,ii,l,k) = hev(j,ii,l,k) -
     &                              d2e(j,l)*redii
                              end do
                           end do
                        end if
                     else if (redi .eq. 0.0d0) then
                        do j = 1, 3
                           do l = j, 3
                              hev(j,i,l,i) = hev(j,i,l,i) + d2e(j,l)
                              hev(j,k,l,k) = hev(j,k,l,k) +
     &                           d2e(j,l)*redk*redk
                              hev(j,kk,l,kk) = hev(j,kk,l,kk) +
     &                           d2e(j,l)*redkk*redkk
                           end do
                        end do
                        do j = 1, 3
                           do l = 1, 3
                              hev(j,i,l,k) = hev(j,i,l,k) -
     &                           d2e(j,l)*redk
                           end do
                        end do
                        if (k .lt. kk) then
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,k,l,kk) = hev(j,k,l,kk) +
     &                              d2e(j,l)*redk*redkk
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,kk,l,k) = hev(j,kk,l,k) +
     &                              d2e(j,l)*redk*redkk
                              end do
                           end do
                        end if
                        if (i .lt. kk) then
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,i,l,kk) = hev(j,i,l,kk) -
     &                              d2e(j,l)*redkk
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,kk,l,i) = hev(j,kk,l,i) -
     &                              d2e(j,l)*redkk
                              end do
                           end do
                        end if
                     else
                        do j = 1, 3
                           do l = j, 3
                              hev(j,i,l,i) = hev(j,i,l,i) +
     &                           d2e(j,l)*redi*redi
                              hev(j,ii,l,ii) = hev(j,ii,l,ii) +
     &                           d2e(j,l)*redii*redii
                              hev(j,k,l,k) = hev(j,k,l,k) +
     &                           d2e(j,l)*redk*redk
                              hev(j,kk,l,kk) = hev(j,kk,l,kk) +
     &                           d2e(j,l)*redkk*redkk
                           end do
                        end do
                        do j = 1, 3
                           do l = 1, 3
                              hev(j,i,l,k) = hev(j,i,l,k) -
     &                           d2e(j,l)*redi*redk
                           end do
                        end do
                        if (i .lt. ii) then
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,i,l,ii) = hev(j,i,l,ii) +
     &                              d2e(j,l)*redi*redii
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,ii,l,i) = hev(j,ii,l,i) +
     &                              d2e(j,l)*redi*redii
                              end do
                           end do
                        end if
                        if (k .lt. kk) then
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,k,l,kk) = hev(j,k,l,kk) +
     &                              d2e(j,l)*redk*redkk
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,kk,l,k) = hev(j,kk,l,k) +
     &                              d2e(j,l)*redk*redkk
                              end do
                           end do
                        end if
                        if (i .lt. kk) then
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,i,l,kk) = hev(j,i,l,kk) -
     &                              d2e(j,l)*redi*redkk
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,kk,l,i) = hev(j,kk,l,i) -
     &                              d2e(j,l)*redi*redkk
                              end do
                           end do
                        end if
                        if (k .lt. ii) then
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,k,l,ii) = hev(j,k,l,ii) -
     &                              d2e(j,l)*redk*redii
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,ii,l,k) = hev(j,ii,l,k) -
     &                              d2e(j,l)*redk*redii
                              end do
                           end do
                        end if
                        if (ii .lt. kk) then
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,ii,l,kk) = hev(j,ii,l,kk) -
     &                              d2e(j,l)*redii*redkk
                              end do
                           end do
                        else
                           do j = 1, 3
                              do l = 1, 3
                                 hev(j,kk,l,ii) = hev(j,kk,l,ii) -
     &                              d2e(j,l)*redii*redkk
                              end do
                           end do
                        end if
                     end if
                  end if
               end if
            end if
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine charge5  --  charge hessian; cartesian version  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "charge5" calculates second derivatives (hessian) of charge
c     interaction energy with respect to the cartesian coordinates
c
c
      subroutine charge5
      implicit none
      include 'sizes.for'
      integer i,j,k,l,ii,kk,n,ndc,nion,iion
      integer n12,i12,n13,i13,skip(maxatm)
      real*8 f,fi,fik,de,term,termx,termy,termz,d2e(3,3)
      real*8 xi,yi,zi,rikx,riky,rikz,rik2,rik,cutoff
      real*8 x,y,z,pchg,dielec,dieled,chgcut,dplcut
      real*8 heb,hea,heba,het,hev,he14,hec,hecd,hed
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm)
      common /attach/ n12(maxatm),i12(4,maxatm),n13(maxatm),
     &                i13(12,maxatm)
      common /charge/ nion,iion(maxatm),pchg(maxatm)
      common /electr/ ndc,dielec,dieled,chgcut,dplcut
      common /hessn1/ heb(3,maxhes,3,maxhes),hea(3,maxhes,3,maxhes),
     &                heba(3,maxhes,3,maxhes),het(3,maxhes,3,maxhes),
     &                hev(3,maxhes,3,maxhes),he14(3,maxhes,3,maxhes),
     &                hec(3,maxhes,3,maxhes),hecd(3,maxhes,3,maxhes),
     &                hed(3,maxhes,3,maxhes)
c
c
c     zero out charge interaction energy hessian
c     elements and set the constants
c
      do l = 1, n
         do k = 1, 3
            do j = 1, n
               do i = 1, 3
                  hec(i,j,k,l) = 0.0d0
               end do
            end do
         end do
      end do
      f = 332.05382d0 / dielec
      cutoff = chgcut
c
c     calculate charge interaction energy hessian; excluding
c     atoms bound to each other or to a common atom
c
      do ii = 1, nion-1
         i = iion(ii)
         do k = 1, n12(i)
            skip(i12(k,i)) = i
         end do
         do k = 1, n13(i)
            skip(i13(k,i)) = i
         end do
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
         do kk = ii+1, nion
            k = iion(kk)
            if (skip(k) .ne. i) then
               rikx = xi - x(k)
               riky = yi - y(k)
               rikz = zi - z(k)
               rik2 = rikx*rikx + riky*riky + rikz*rikz
               if (rik2 .le. chgcut) then
                  rik = sqrt(rik2)
                  fik = fi * pchg(kk)
c
c     set chain rule terms for hessian matrix elements
c
                  de = fik / (rik*rik2)
                  term = 3.0d0 * de/rik2
                  termx = term * rikx
                  termy = term * riky
                  termz = term * rikz
                  d2e(1,1) = termx*rikx - de
                  d2e(1,2) = termx*riky
                  d2e(1,3) = termx*rikz
                  d2e(2,1) = d2e(1,2)
                  d2e(2,2) = termy*riky - de
                  d2e(2,3) = termy*rikz
                  d2e(3,1) = d2e(1,3)
                  d2e(3,2) = d2e(2,3)
                  d2e(3,3) = termz*rikz - de
c
c     increment diagonal and above-diagonal hessian elements
c
                  do j = 1, 3
                     do l = j, 3
                        hec(j,i,l,i) = hec(j,i,l,i) + d2e(j,l)
                        hec(j,k,l,k) = hec(j,k,l,k) + d2e(j,l)
                     end do
                  end do
                  do j = 1, 3
                     do l = 1, 3
                        hec(j,i,l,k) = hec(j,i,l,k) - d2e(j,l)
                     end do
                  end do
               end if
            end if
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine chgdpl5  --  charge-dipole hessian; cartesian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "chgdpl5" calculates second derivatives (hessian) of charge-
c     dipole energy with respect to the cartesian coordinates
c
c
      subroutine chgdpl5
      implicit none
      include 'sizes.for'
      integer i,j,k,l,i1,k1,k2,n,n12,i12,skip(maxatm)
      integer ndc,ndipole,mup,mun,nion,iion
      real*8 x,y,z,dielec,dieled,chgcut,dplcut,pchg,bdpl
      real*8 f,fi,fik,cutoff
      real*8 xi,yi,zi,xk,yk,zk,xr,yr,zr
      real*8 r2,rk2,rkr3,dotk,term,term2,part,part2
      real*8 termx,termy,termz,termxk,termyk,termzk
      real*8 xrr2,yrr2,zrr2,xkrk2,ykrk2,zkrk2
      real*8 dotk2,dotkr2,dotkrk2,factor,factork
      real*8 dtdxi1,dtdyi1,dtdzi1,dtdxk1,dtdyk1,dtdzk1
      real*8 dtdxk2,dtdyk2,dtdzk2
      real*8 dtxdxi1,dtxkdxi1,dtxdxk1,dtxkdxk1,dtxdxk2,dtxkdxk2
      real*8 dtydxi1,dtykdxi1,dtydxk1,dtykdxk1,dtydxk2,dtykdxk2
      real*8 dtzdxi1,dtzkdxi1,dtzdxk1,dtzkdxk1,dtzdxk2,dtzkdxk2
      real*8 dtxdyi1,dtxkdyi1,dtxdyk1,dtxkdyk1,dtxdyk2,dtxkdyk2
      real*8 dtydyi1,dtykdyi1,dtydyk1,dtykdyk1,dtydyk2,dtykdyk2
      real*8 dtzdyi1,dtzkdyi1,dtzdyk1,dtzkdyk1,dtzdyk2,dtzkdyk2
      real*8 dtxdzi1,dtxkdzi1,dtxdzk1,dtxkdzk1,dtxdzk2,dtxkdzk2
      real*8 dtydzi1,dtykdzi1,dtydzk1,dtykdzk1,dtydzk2,dtykdzk2
      real*8 dtzdzi1,dtzkdzi1,dtzdzk1,dtzkdzk1,dtzdzk2,dtzkdzk2
      real*8 heb,hea,heba,het,hev,he14,hec,hecd,hed
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm)
      common /attach/ n12(maxatm),i12(4,maxatm)
      common /charge/ nion,iion(maxatm),pchg(maxatm)
      common /dipole/ ndipole,mup(maxbnd),mun(maxbnd),bdpl(maxbnd)
      common /electr/ ndc,dielec,dieled,chgcut,dplcut
      common /hessn1/ heb(3,maxhes,3,maxhes),hea(3,maxhes,3,maxhes),
     &                heba(3,maxhes,3,maxhes),het(3,maxhes,3,maxhes),
     &                hev(3,maxhes,3,maxhes),he14(3,maxhes,3,maxhes),
     &                hec(3,maxhes,3,maxhes),hecd(3,maxhes,3,maxhes),
     &                hed(3,maxhes,3,maxhes)
c
c
c     zero out charge-dipole interaction hessian elements,
c     then set up the constants for the calculation
c
      do l = 1, n
         do k = 1, 3
            do j = 1, n
               do i = 1, 3
                  hecd(i,j,k,l) = 0.0d0
               end do
            end do
         end do
      end do
      f = 69.120d0 / sqrt(dielec*dieled)
      cutoff = sqrt(chgcut*dplcut)
c
c     calculate the total energy by suming up the
c     contributions over each charge-dipole pair
c
      do i = 1, nion
         i1 = iion(i)
         skip(i1) = i1
         do k = 1, n12(i1)
            skip(i12(k,i1)) = i1
         end do
         xi = x(i1)
         yi = y(i1)
         zi = z(i1)
         fi = f * pchg(i)
         do k = 1, ndipole
            k1 = mup(k)
            k2 = mun(k)
            if (skip(k1).ne.i1 .and. skip(k2).ne.i1) then
               xr = xi - 0.5d0*(x(k1)+x(k2))
               yr = yi - 0.5d0*(y(k1)+y(k2))
               zr = zi - 0.5d0*(z(k1)+z(k2))
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. cutoff) then
                  xk = x(k2) - x(k1)
                  yk = y(k2) - y(k1)
                  zk = z(k2) - z(k1)
                  rk2 = xk*xk + yk*yk + zk*zk
                  rkr3 = sqrt(rk2*r2) * r2
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fi * bdpl(k)
c
c     some abbreviations used in various chain rule terms
c
                  xrr2 = xr / r2
                  yrr2 = yr / r2
                  zrr2 = zr / r2
                  xkrk2 = xk / rk2
                  ykrk2 = yk / rk2
                  zkrk2 = zk / rk2
                  dotk2 = 2.0d0 * dotk
                  dotkr2 = dotk / r2
                  dotkrk2 = dotk / rk2
c
c     now, form the chain rule terms for first derivatives
c
                  term = fik / rkr3
                  term2 = -3.0d0 * dotk
                  termx = 0.5d0 * term * (xk+xrr2*term2)
                  termy = 0.5d0 * term * (yk+yrr2*term2)
                  termz = 0.5d0 * term * (zk+zrr2*term2)
                  termxk = term * (xr-dotk*xkrk2)
                  termyk = term * (yr-dotk*ykrk2)
                  termzk = term * (zr-dotk*zkrk2)
c
c     next, find the second derivative chain rule terms
c
                  dtdxi1 = -3.0d0 * xrr2
                  part = xk - dotk2*xrr2
                  factor = -1.5d0 * (dotkr2 + xrr2*part)
                  factork = 1.0d0 - xk*xkrk2
                  dtxdxi1 = dtdxi1*termx + term*factor
                  dtxkdxi1 = dtdxi1*termxk + term*factork
                  factor = -1.5d0 * yrr2 * part
                  factork = -yk * xkrk2
                  dtydxi1 = dtdxi1*termy + term*factor
                  dtykdxi1 = dtdxi1*termyk + term*factork
                  factor = -1.5d0 * zrr2 * part
                  factork = -zk * xkrk2
                  dtzdxi1 = dtdxi1*termz + term*factor
                  dtzkdxi1 = dtdxi1*termzk + term*factork
c
                  dtdyi1 = -3.0d0 * yrr2
                  part = yk - dotk2*yrr2
                  factor = -1.5d0 * xrr2 * part
                  factork = -xk * ykrk2
                  dtxdyi1 = dtdyi1*termx + term*factor
                  dtxkdyi1 = dtdyi1*termxk + term*factork
                  factor = -1.5d0 * (dotkr2 + yrr2*part)
                  factork = 1.0d0 - yk*ykrk2
                  dtydyi1 = dtdyi1*termy + term*factor
                  dtykdyi1 = dtdyi1*termyk + term*factork
                  factor = -1.5d0 * zrr2 * part
                  factork = -zk * ykrk2
                  dtzdyi1 = dtdyi1*termz + term*factor
                  dtzkdyi1 = dtdyi1*termzk + term*factork
c
                  dtdzi1 = -3.0d0 * zrr2
                  part = zk - dotk2*zrr2
                  factor = -1.5d0 * xrr2 * part
                  factork = -xk * zkrk2
                  dtxdzi1 = dtdzi1*termx + term*factor
                  dtxkdzi1 = dtdzi1*termxk + term*factork
                  factor = -1.5d0 * yrr2 * part
                  factork = -yk * zkrk2
                  dtydzi1 = dtdzi1*termy + term*factor
                  dtykdzi1 = dtdzi1*termyk + term*factork
                  factor = -1.5d0 * (dotkr2 + zrr2*part)
                  factork = 1.0d0 - zk*zkrk2
                  dtzdzi1 = dtdzi1*termz + term*factor
                  dtzkdzi1 = dtdzi1*termzk + term*factork
c
                  dtdxk1 = 1.5d0*xrr2 + xkrk2
                  part = 0.5d0*xk + xr
                  part2 = dotk*xrr2 - part
                  factor = -0.5d0 - 1.5d0*xrr2*part2 + 0.75d0*dotkr2
                  factork = -0.5d0 - dotk2*xkrk2*xkrk2
     &                         + xkrk2*part + dotkrk2
                  dtxdxk1 = dtdxk1*termx + term*factor
                  dtxkdxk1 = dtdxk1*termxk + term*factork
                  factor = -1.5d0 * yrr2 * part2
                  factork = -dotk2*ykrk2*xkrk2 + ykrk2*part
                  dtydxk1 = dtdxk1*termy + term*factor
                  dtykdxk1 = dtdxk1*termyk + term*factork
                  factor = -1.5d0 * zrr2 * part2
                  factork = -dotk2*zkrk2*xkrk2 + zkrk2*part
                  dtzdxk1 = dtdxk1*termz + term*factor
                  dtzkdxk1 = dtdxk1*termzk + term*factork
c
                  dtdyk1 = 1.5d0*yrr2 + ykrk2
                  part = 0.5d0*yk + yr
                  part2 = dotk*yrr2 - part
                  factor = -1.5d0 * xrr2 * part2
                  factork = -dotk2*xkrk2*ykrk2 + xkrk2*part
                  dtxdyk1 = dtdyk1*termx + term*factor
                  dtxkdyk1 = dtdyk1*termxk + term*factork
                  factor = -0.5d0 - 1.5d0*yrr2*part2 + 0.75d0*dotkr2
                  factork = -0.5d0 - dotk2*ykrk2*ykrk2
     &                         + ykrk2*part + dotkrk2
                  dtydyk1 = dtdyk1*termy + term*factor
                  dtykdyk1 = dtdyk1*termyk + term*factork
                  factor = -1.5d0 * zrr2 * part2
                  factork = -dotk2*zkrk2*ykrk2 + zkrk2*part
                  dtzdyk1 = dtdyk1*termz + term*factor
                  dtzkdyk1 = dtdyk1*termzk + term*factork
c
                  dtdzk1 = 1.5d0*zrr2 + zkrk2
                  part = 0.5d0*zk + zr
                  part2 = dotk*zrr2 - part
                  factor = -1.5d0 * xrr2 * part2
                  factork = -dotk2*xkrk2*zkrk2 + xkrk2*part
                  dtxdzk1 = dtdzk1*termx + term*factor
                  dtxkdzk1 = dtdzk1*termxk + term*factork
                  factor = -1.5d0 * yrr2 * part2
                  factork = -dotk2*ykrk2*zkrk2 + ykrk2*part
                  dtydzk1 = dtdzk1*termy + term*factor
                  dtykdzk1 = dtdzk1*termyk + term*factork
                  factor = -0.5d0 - 1.5d0*zrr2*part2 + 0.75d0*dotkr2
                  factork = -0.5d0 - dotk2*zkrk2*zkrk2
     &                         + zkrk2*part + dotkrk2
                  dtzdzk1 = dtdzk1*termz + term*factor
                  dtzkdzk1 = dtdzk1*termzk + term*factork
c
                  dtdxk2 = 1.5d0*xrr2 - xkrk2
                  part = 0.5d0*xk - xr
                  part2 = dotk*xrr2 - part
                  factor = 0.5d0 - 1.5d0*xrr2*part2 + 0.75d0*dotkr2
                  factork = -0.5d0 + dotk2*xkrk2*xkrk2
     &                         + xkrk2*part - dotkrk2
                  dtxdxk2 = dtdxk2*termx + term*factor
                  dtxkdxk2 = dtdxk2*termxk + term*factork
                  factor = -1.5d0 * yrr2 * part2
                  factork = dotk2*ykrk2*xkrk2 + ykrk2*part
                  dtydxk2 = dtdxk2*termy + term*factor
                  dtykdxk2 = dtdxk2*termyk + term*factork
                  factor = -1.5d0 * zrr2 * part2
                  factork = dotk2*zkrk2*xkrk2 + zkrk2*part
                  dtzdxk2 = dtdxk2*termz + term*factor
                  dtzkdxk2 = dtdxk2*termzk + term*factork
c
                  dtdyk2 = 1.5d0*yrr2 - ykrk2
                  part = 0.5d0*yk - yr
                  part2 = dotk*yrr2 - part
                  factor = -1.5d0 * xrr2 * part2
                  factork = dotk2*xkrk2*ykrk2 + xkrk2*part
                  dtxdyk2 = dtdyk2*termx + term*factor
                  dtxkdyk2 = dtdyk2*termxk + term*factork
                  factor = 0.5d0 - 1.5d0*yrr2*part2 + 0.75d0*dotkr2
                  factork = -0.5d0 + dotk2*ykrk2*ykrk2
     &                         + ykrk2*part - dotkrk2
                  dtydyk2 = dtdyk2*termy + term*factor
                  dtykdyk2 = dtdyk2*termyk + term*factork
                  factor = -1.5d0 * zrr2 * part2
                  factork = dotk2*zkrk2*ykrk2 + zkrk2*part
                  dtzdyk2 = dtdyk2*termz + term*factor
                  dtzkdyk2 = dtdyk2*termzk + term*factork
c
                  dtdzk2 = 1.5d0*zrr2 - zkrk2
                  part = 0.5d0*zk - zr
                  part2 = dotk*zrr2 - part
                  factor = -1.5d0 * xrr2 * part2
                  factork = dotk2*xkrk2*zkrk2 + xkrk2*part
                  dtxdzk2 = dtdzk2*termx + term*factor
                  dtxkdzk2 = dtdzk2*termxk + term*factork
                  factor = -1.5d0 * yrr2 * part2
                  factork = dotk2*ykrk2*zkrk2 + ykrk2*part
                  dtydzk2 = dtdzk2*termy + term*factor
                  dtykdzk2 = dtdzk2*termyk + term*factork
                  factor = 0.5d0 - 1.5d0*zrr2*part2 + 0.75d0*dotkr2
                  factork = -0.5d0 + dotk2*zkrk2*zkrk2
     &                         + zkrk2*part - dotkrk2
                  dtzdzk2 = dtdzk2*termz + term*factor
                  dtzkdzk2 = dtdzk2*termzk + term*factork
c
c     now, increment the diagonal hessian block elements
c
                  hecd(1,i1,1,i1) = hecd(1,i1,1,i1) + 2.0d0*dtxdxi1
                  hecd(1,i1,2,i1) = hecd(1,i1,2,i1) + 2.0d0*dtydxi1
                  hecd(2,i1,2,i1) = hecd(2,i1,2,i1) + 2.0d0*dtydyi1
                  hecd(1,i1,3,i1) = hecd(1,i1,3,i1) + 2.0d0*dtzdxi1
                  hecd(2,i1,3,i1) = hecd(2,i1,3,i1) + 2.0d0*dtzdyi1
                  hecd(3,i1,3,i1) = hecd(3,i1,3,i1) + 2.0d0*dtzdzi1
                  hecd(1,k1,1,k1) = hecd(1,k1,1,k1) - dtxdxk1 - dtxkdxk1
                  hecd(1,k1,2,k1) = hecd(1,k1,2,k1) - dtydxk1 - dtykdxk1
                  hecd(2,k1,2,k1) = hecd(2,k1,2,k1) - dtydyk1 - dtykdyk1
                  hecd(1,k1,3,k1) = hecd(1,k1,3,k1) - dtzdxk1 - dtzkdxk1
                  hecd(2,k1,3,k1) = hecd(2,k1,3,k1) - dtzdyk1 - dtzkdyk1
                  hecd(3,k1,3,k1) = hecd(3,k1,3,k1) - dtzdzk1 - dtzkdzk1
                  hecd(1,k2,1,k2) = hecd(1,k2,1,k2) - dtxdxk2 + dtxkdxk2
                  hecd(1,k2,2,k2) = hecd(1,k2,2,k2) - dtydxk2 + dtykdxk2
                  hecd(2,k2,2,k2) = hecd(2,k2,2,k2) - dtydyk2 + dtykdyk2
                  hecd(1,k2,3,k2) = hecd(1,k2,3,k2) - dtzdxk2 + dtzkdxk2
                  hecd(2,k2,3,k2) = hecd(2,k2,3,k2) - dtzdyk2 + dtzkdyk2
                  hecd(3,k2,3,k2) = hecd(3,k2,3,k2) - dtzdzk2 + dtzkdzk2
c
c     finally, increment above-diagonal hessian block elements
c
                  if (i1 .lt. k1) then
                     hecd(1,i1,1,k1) = hecd(1,i1,1,k1)-dtxdxi1-dtxkdxi1
                     hecd(2,i1,1,k1) = hecd(2,i1,1,k1)-dtxdyi1-dtxkdyi1
                     hecd(3,i1,1,k1) = hecd(3,i1,1,k1)-dtxdzi1-dtxkdzi1
                     hecd(1,i1,2,k1) = hecd(1,i1,2,k1)-dtydxi1-dtykdxi1
                     hecd(2,i1,2,k1) = hecd(2,i1,2,k1)-dtydyi1-dtykdyi1
                     hecd(3,i1,2,k1) = hecd(3,i1,2,k1)-dtydzi1-dtykdzi1
                     hecd(1,i1,3,k1) = hecd(1,i1,3,k1)-dtzdxi1-dtzkdxi1
                     hecd(2,i1,3,k1) = hecd(2,i1,3,k1)-dtzdyi1-dtzkdyi1
                     hecd(3,i1,3,k1) = hecd(3,i1,3,k1)-dtzdzi1-dtzkdzi1
                  else
                     hecd(1,k1,1,i1) = hecd(1,k1,1,i1) + 2.0d0*dtxdxk1
                     hecd(2,k1,1,i1) = hecd(2,k1,1,i1) + 2.0d0*dtxdyk1
                     hecd(3,k1,1,i1) = hecd(3,k1,1,i1) + 2.0d0*dtxdzk1
                     hecd(1,k1,2,i1) = hecd(1,k1,2,i1) + 2.0d0*dtydxk1
                     hecd(2,k1,2,i1) = hecd(2,k1,2,i1) + 2.0d0*dtydyk1
                     hecd(3,k1,2,i1) = hecd(3,k1,2,i1) + 2.0d0*dtydzk1
                     hecd(1,k1,3,i1) = hecd(1,k1,3,i1) + 2.0d0*dtzdxk1
                     hecd(2,k1,3,i1) = hecd(2,k1,3,i1) + 2.0d0*dtzdyk1
                     hecd(3,k1,3,i1) = hecd(3,k1,3,i1) + 2.0d0*dtzdzk1
                  end if
c
                  if (i1 .lt. k2) then
                     hecd(1,i1,1,k2) = hecd(1,i1,1,k2)-dtxdxi1+dtxkdxi1
                     hecd(2,i1,1,k2) = hecd(2,i1,1,k2)-dtxdyi1+dtxkdyi1
                     hecd(3,i1,1,k2) = hecd(3,i1,1,k2)-dtxdzi1+dtxkdzi1
                     hecd(1,i1,2,k2) = hecd(1,i1,2,k2)-dtydxi1+dtykdxi1
                     hecd(2,i1,2,k2) = hecd(2,i1,2,k2)-dtydyi1+dtykdyi1
                     hecd(3,i1,2,k2) = hecd(3,i1,2,k2)-dtydzi1+dtykdzi1
                     hecd(1,i1,3,k2) = hecd(1,i1,3,k2)-dtzdxi1+dtzkdxi1
                     hecd(2,i1,3,k2) = hecd(2,i1,3,k2)-dtzdyi1+dtzkdyi1
                     hecd(3,i1,3,k2) = hecd(3,i1,3,k2)-dtzdzi1+dtzkdzi1
                  else
                     hecd(1,k2,1,i1) = hecd(1,k2,1,i1) + 2.0d0*dtxdxk2
                     hecd(2,k2,1,i1) = hecd(2,k2,1,i1) + 2.0d0*dtxdyk2
                     hecd(3,k2,1,i1) = hecd(3,k2,1,i1) + 2.0d0*dtxdzk2
                     hecd(1,k2,2,i1) = hecd(1,k2,2,i1) + 2.0d0*dtydxk2
                     hecd(2,k2,2,i1) = hecd(2,k2,2,i1) + 2.0d0*dtydyk2
                     hecd(3,k2,2,i1) = hecd(3,k2,2,i1) + 2.0d0*dtydzk2
                     hecd(1,k2,3,i1) = hecd(1,k2,3,i1) + 2.0d0*dtzdxk2
                     hecd(2,k2,3,i1) = hecd(2,k2,3,i1) + 2.0d0*dtzdyk2
                     hecd(3,k2,3,i1) = hecd(3,k2,3,i1) + 2.0d0*dtzdzk2
                  end if
c
                  if (k1 .lt. k2) then
                     hecd(1,k1,1,k2) = hecd(1,k1,1,k2)-dtxdxk1+dtxkdxk1
                     hecd(2,k1,1,k2) = hecd(2,k1,1,k2)-dtxdyk1+dtxkdyk1
                     hecd(3,k1,1,k2) = hecd(3,k1,1,k2)-dtxdzk1+dtxkdzk1
                     hecd(1,k1,2,k2) = hecd(1,k1,2,k2)-dtydxk1+dtykdxk1
                     hecd(2,k1,2,k2) = hecd(2,k1,2,k2)-dtydyk1+dtykdyk1
                     hecd(3,k1,2,k2) = hecd(3,k1,2,k2)-dtydzk1+dtykdzk1
                     hecd(1,k1,3,k2) = hecd(1,k1,3,k2)-dtzdxk1+dtzkdxk1
                     hecd(2,k1,3,k2) = hecd(2,k1,3,k2)-dtzdyk1+dtzkdyk1
                     hecd(3,k1,3,k2) = hecd(3,k1,3,k2)-dtzdzk1+dtzkdzk1
                  else
                     hecd(1,k2,1,k1) = hecd(1,k2,1,k1)-dtxdxk2-dtxkdxk2
                     hecd(2,k2,1,k1) = hecd(2,k2,1,k1)-dtxdyk2-dtxkdyk2
                     hecd(3,k2,1,k1) = hecd(3,k2,1,k1)-dtxdzk2-dtxkdzk2
                     hecd(1,k2,2,k1) = hecd(1,k2,2,k1)-dtydxk2-dtykdxk2
                     hecd(2,k2,2,k1) = hecd(2,k2,2,k1)-dtydyk2-dtykdyk2
                     hecd(3,k2,2,k1) = hecd(3,k2,2,k1)-dtydzk2-dtykdzk2
                     hecd(1,k2,3,k1) = hecd(1,k2,3,k1)-dtzdxk2-dtzkdxk2
                     hecd(2,k2,3,k1) = hecd(2,k2,3,k1)-dtzdyk2-dtzkdyk2
                     hecd(3,k2,3,k1) = hecd(3,k2,3,k1)-dtzdzk2-dtzkdzk2
                  end if
               end if
            end if
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dipole5  --  dipole hessian; cartesian version  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dipole5" calculates second derivatives (hessian) of dipole
c     interaction energy with respect to the cartesian coordinates
c
c
      subroutine dipole5
      implicit none
      include 'sizes.for'
      integer i,j,k,l,n,ndc,ndipole,mup,mun,i1,i2,k1,k2
      real*8 x,y,z,bdpl,dielec,dieled,chgcut,dplcut
      real*8 f,fi,fik,cutoff
      real*8 xi,yi,zi,xk,yk,zk,xq,yq,zq,xr,yr,zr
      real*8 e,r2,ri2,rk2,rirkr5,doti,dotk,dotp
      real*8 term,dedr,dedrirk,deddoti,deddotk,deddotp
      real*8 termx,termy,termz
      real*8 termxi,termyi,termzi,termxk,termyk,termzk
      real*8 enum,r2inv,ri2inv,rk2inv,dotik,xrr2,yrr2,zrr2
      real*8 xiri2,yiri2,ziri2,xkrk2,ykrk2,zkrk2
      real*8 xixr,xiyr,xizr,yixr,yiyr,yizr,zixr,ziyr,zizr
      real*8 xkxr,xkyr,xkzr,ykxr,ykyr,ykzr,zkxr,zkyr,zkzr
      real*8 xixk,xiyk,xizk,yixk,yiyk,yizk,zixk,ziyk,zizk
      real*8 xrxr,xryr,xrzr,yryr,yrzr,zrzr
      real*8 xidotk,yidotk,zidotk,xkdoti,ykdoti,zkdoti
      real*8 factor,factori,factork,part,partik
      real*8 dtdxi1,dtdyi1,dtdzi1,dtdxi2,dtdyi2,dtdzi2
      real*8 dtdxk1,dtdyk1,dtdzk1,dtdxk2,dtdyk2,dtdzk2
      real*8 dtxdxi1,dtxidxi1,dtxkdxi1,dtxdxi2,dtxidxi2,dtxkdxi2
      real*8 dtydxi1,dtyidxi1,dtykdxi1,dtydxi2,dtyidxi2,dtykdxi2
      real*8 dtzdxi1,dtzidxi1,dtzkdxi1,dtzdxi2,dtzidxi2,dtzkdxi2
      real*8 dtxdyi1,dtxidyi1,dtxkdyi1,dtxdyi2,dtxidyi2,dtxkdyi2
      real*8 dtydyi1,dtyidyi1,dtykdyi1,dtydyi2,dtyidyi2,dtykdyi2
      real*8 dtzdyi1,dtzidyi1,dtzkdyi1,dtzdyi2,dtzidyi2,dtzkdyi2
      real*8 dtxdzi1,dtxidzi1,dtxkdzi1,dtxdzi2,dtxidzi2,dtxkdzi2
      real*8 dtydzi1,dtyidzi1,dtykdzi1,dtydzi2,dtyidzi2,dtykdzi2
      real*8 dtzdzi1,dtzidzi1,dtzkdzi1,dtzdzi2,dtzidzi2,dtzkdzi2
      real*8 dtxdxk1,dtxidxk1,dtxkdxk1,dtxdxk2,dtxidxk2,dtxkdxk2
      real*8 dtydxk1,dtyidxk1,dtykdxk1,dtydxk2,dtyidxk2,dtykdxk2
      real*8 dtzdxk1,dtzidxk1,dtzkdxk1,dtzdxk2,dtzidxk2,dtzkdxk2
      real*8 dtxdyk1,dtxidyk1,dtxkdyk1,dtxdyk2,dtxidyk2,dtxkdyk2
      real*8 dtydyk1,dtyidyk1,dtykdyk1,dtydyk2,dtyidyk2,dtykdyk2
      real*8 dtzdyk1,dtzidyk1,dtzkdyk1,dtzdyk2,dtzidyk2,dtzkdyk2
      real*8 dtxdzk1,dtxidzk1,dtxkdzk1,dtxdzk2,dtxidzk2,dtxkdzk2
      real*8 dtydzk1,dtyidzk1,dtykdzk1,dtydzk2,dtyidzk2,dtykdzk2
      real*8 dtzdzk1,dtzidzk1,dtzkdzk1,dtzdzk2,dtzidzk2,dtzkdzk2
      real*8 heb,hea,heba,het,hev,he14,hec,hecd,hed
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm)
      common /dipole/ ndipole,mup(maxbnd),mun(maxbnd),bdpl(maxbnd)
      common /electr/ ndc,dielec,dieled,chgcut,dplcut
      common /hessn1/ heb(3,maxhes,3,maxhes),hea(3,maxhes,3,maxhes),
     &                heba(3,maxhes,3,maxhes),het(3,maxhes,3,maxhes),
     &                hev(3,maxhes,3,maxhes),he14(3,maxhes,3,maxhes),
     &                hec(3,maxhes,3,maxhes),hecd(3,maxhes,3,maxhes),
     &                hed(3,maxhes,3,maxhes)
c
c
c     zero out the dipole interaction energy hessian elements
c
      do l = 1, n
         do k = 1, 3
            do j = 1, n
               do i = 1, 3
                  hed(i,j,k,l) = 0.0d0
               end do
            end do
         end do
      end do
      f = 14.388d0 / dieled
      cutoff = dplcut
c
c     calculate the dipole interaction energy second derivatives
c     by suming up the contributions over each pair of dipoles
c
      do i = 1, ndipole-1
         i1 = mup(i)
         i2 = mun(i)
         xi = x(i2) - x(i1)
         yi = y(i2) - y(i1)
         zi = z(i2) - z(i1)
         ri2 = xi*xi + yi*yi + zi*zi
         xq = x(i1) + x(i2)
         yq = y(i1) + y(i2)
         zq = z(i1) + z(i2)
         fi = f * bdpl(i)
         do k = i+1, ndipole
            k1 = mup(k)
            k2 = mun(k)
            if (k1.ne.i1 .and. k1.ne.i2 .and.
     &          k2.ne.i1 .and. k2.ne.i2) then
               xr = 0.5d0 * (xq-x(k1)-x(k2))
               yr = 0.5d0 * (yq-y(k1)-y(k2))
               zr = 0.5d0 * (zq-z(k1)-z(k2))
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. cutoff) then
                  xk = x(k2) - x(k1)
                  yk = y(k2) - y(k1)
                  zk = z(k2) - z(k1)
                  rk2 = xk*xk + yk*yk + zk*zk
                  rirkr5 = sqrt(ri2*rk2*r2) * r2 * r2
                  dotp = xi*xk + yi*yk + zi*zk
                  doti = xi*xr + yi*yr + zi*zr
                  dotk = xk*xr + yk*yr + zk*zr
                  fik = fi * bdpl(k)
c
c     some abbreviations used in various chain rule terms
c
                  dotik = doti * dotk
                  enum = dotp*r2 - 3.0d0*dotik
                  r2inv = 3.75d0 / r2
                  ri2inv = 1.0d0 / ri2
                  rk2inv = 1.0d0 / rk2
                  xrr2 = xr / r2
                  yrr2 = yr / r2
                  zrr2 = zr / r2
                  xiri2 = xi / ri2
                  yiri2 = yi / ri2
                  ziri2 = zi / ri2
                  xkrk2 = xk / rk2
                  ykrk2 = yk / rk2
                  zkrk2 = zk / rk2
                  xixr = xi * xr
                  xiyr = xi * yr
                  xizr = xi * zr
                  yixr = yi * xr
                  yiyr = yi * yr
                  yizr = yi * zr
                  zixr = zi * xr
                  ziyr = zi * yr
                  zizr = zi * zr
                  xkxr = xk * xr
                  xkyr = xk * yr
                  xkzr = xk * zr
                  ykxr = yk * xr
                  ykyr = yk * yr
                  ykzr = yk * zr
                  zkxr = zk * xr
                  zkyr = zk * yr
                  zkzr = zk * zr
                  xixk = xi * xk
                  xiyk = xi * yk
                  xizk = xi * zk
                  yixk = yi * xk
                  yiyk = yi * yk
                  yizk = yi * zk
                  zixk = zi * xk
                  ziyk = zi * yk
                  zizk = zi * zk
                  xrxr = 3.0d0 * xr * xr
                  xryr = 3.0d0 * xr * yr
                  xrzr = 3.0d0 * xr * zr
                  yryr = 3.0d0 * yr * yr
                  yrzr = 3.0d0 * yr * zr
                  zrzr = 3.0d0 * zr * zr
                  xidotk = xi * dotk
                  yidotk = yi * dotk
                  zidotk = zi * dotk
                  xkdoti = xk * doti
                  ykdoti = yk * doti
                  zkdoti = zk * doti
c
c     now, form the chain rule terms for first derivatives
c
                  term = -fik / rirkr5
                  deddotp = -term * r2
                  deddoti = term * 3.0d0*dotk
                  deddotk = term * 3.0d0*doti
                  dedr = term * (1.5d0*dotp-7.5d0*dotik/r2)
                  dedrirk = term * enum
c
c     more first derivative chain rule expressions
c
                  termx = dedr*xr + 0.5d0*(deddoti*xi + deddotk*xk)
                  termy = dedr*yr + 0.5d0*(deddoti*yi + deddotk*yk)
                  termz = dedr*zr + 0.5d0*(deddoti*zi + deddotk*zk)
                  termxi = dedrirk*xiri2 + deddotp*xk + deddoti*xr
                  termyi = dedrirk*yiri2 + deddotp*yk + deddoti*yr
                  termzi = dedrirk*ziri2 + deddotp*zk + deddoti*zr
                  termxk = dedrirk*xkrk2 + deddotp*xi + deddotk*xr
                  termyk = dedrirk*ykrk2 + deddotp*yi + deddotk*yr
                  termzk = dedrirk*zkrk2 + deddotp*zi + deddotk*zr
c
c     next, find the second derivative chain rule terms
c
                  dtdxi1 = -2.5d0*xrr2 + xiri2
                  part = xkdoti - 2.0d0*dotk*xr + xidotk
     &                      - 2.0d0*dotik*xrr2
                  partik = -xk*r2 + dotp*xr - 1.5d0*xkdoti
     &                        + 3.0d0*xr*dotk - 1.5d0*xidotk
                  factor = 0.75d0*dotp - 3.0d0*xkxr + 1.5d0*xixk
     &                        - 1.5d0*dotk - r2inv*(xr*part+dotik)
                  factori = 1.5d0*dotk + 0.5d0*xkxr + xiri2*partik
     &                         - enum*(ri2inv-2.0d0*xiri2*xiri2)
                  factork = r2 + 1.5d0*doti + 0.5d0*xixr
     &                         - xrxr + xkrk2*partik
                  dtxdxi1 = dtdxi1*termx + term*factor
                  dtxidxi1 = dtdxi1*termxi + term*factori
                  dtxkdxi1 = dtdxi1*termxk + term*factork
                  factor = -1.5d0*xkyr - 1.5d0*ykxr + 0.75d0*xiyk
     &                        + 0.75d0*yixk - r2inv*yr*part
                  factori = -ykxr + 1.5d0*xkyr + yiri2*partik
     &                         + enum*(2.0d0*yiri2*xiri2)
                  factork = -yixr - xryr + 1.5d0*xiyr
     &                         + ykrk2*partik
                  dtydxi1 = dtdxi1*termy + term*factor
                  dtyidxi1 = dtdxi1*termyi + term*factori
                  dtykdxi1 = dtdxi1*termyk + term*factork
                  factor = -1.5d0*xkzr - 1.5d0*zkxr + 0.75d0*xizk
     &                        + 0.75d0*zixk - r2inv*zr*part
                  factori = -zkxr + 1.5d0*xkzr + ziri2*partik
     &                         + enum*(2.0d0*ziri2*xiri2)
                  factork = -zixr - xrzr + 1.5d0*xizr
     &                         + zkrk2*partik
                  dtzdxi1 = dtdxi1*termz + term*factor
                  dtzidxi1 = dtdxi1*termzi + term*factori
                  dtzkdxi1 = dtdxi1*termzk + term*factork
c
                  dtdyi1 = -2.5d0*yrr2 + yiri2
                  part = ykdoti - 2.0d0*dotk*yr + yidotk
     &                      - 2.0d0*dotik*yrr2
                  partik = -yk*r2 + dotp*yr - 1.5d0*ykdoti
     &                        + 3.0d0*yr*dotk - 1.5d0*yidotk
                  factor = -1.5d0*ykxr - 1.5d0*xkyr + 0.75d0*yixk
     &                        + 0.75d0*xiyk - r2inv*xr*part
                  factori = -xkyr + 1.5d0*ykxr + xiri2*partik
     &                         + enum*(2.0d0*xiri2*yiri2)
                  factork = -xiyr - xryr + 1.5d0*yixr
     &                         + xkrk2*partik
                  dtxdyi1 = dtdyi1*termx + term*factor
                  dtxidyi1 = dtdyi1*termxi + term*factori
                  dtxkdyi1 = dtdyi1*termxk + term*factork
                  factor = 0.75d0*dotp - 3.0d0*ykyr + 1.5d0*yiyk
     &                        - 1.5d0*dotk - r2inv*(yr*part+dotik)
                  factori = 1.5d0*dotk + 0.5d0*ykyr + yiri2*partik
     &                         - enum*(ri2inv-2.0d0*yiri2*yiri2)
                  factork = r2 + 1.5d0*doti + 0.5d0*yiyr
     &                         - yryr + ykrk2*partik
                  dtydyi1 = dtdyi1*termy + term*factor
                  dtyidyi1 = dtdyi1*termyi + term*factori
                  dtykdyi1 = dtdyi1*termyk + term*factork
                  factor = -1.5d0*ykzr - 1.5d0*zkyr + 0.75d0*yizk
     &                        + 0.75d0*ziyk - r2inv*zr*part
                  factori = -zkyr + 1.5d0*ykzr + ziri2*partik
     &                         + enum*(2.0d0*ziri2*yiri2)
                  factork = -ziyr - yrzr + 1.5d0*yizr
     &                         + zkrk2*partik
                  dtzdyi1 = dtdyi1*termz + term*factor
                  dtzidyi1 = dtdyi1*termzi + term*factori
                  dtzkdyi1 = dtdyi1*termzk + term*factork
c
                  dtdzi1 = -2.5d0*zrr2 + ziri2
                  part = zkdoti - 2.0d0*dotk*zr + zidotk
     &                      - 2.0d0*dotik*zrr2
                  partik = -zk*r2 + dotp*zr - 1.5d0*zkdoti
     &                      + 3.0d0*zr*dotk - 1.5d0*zidotk
                  factor = -1.5d0*zkxr - 1.5d0*xkzr + 0.75d0*zixk
     &                        + 0.75d0*xizk - r2inv*xr*part
                  factori = -xkzr + 1.5d0*zkxr + xiri2*partik
     &                         + enum*(2.0d0*xiri2*ziri2)
                  factork = -xizr - xrzr + 1.5d0*zixr
     &                         + xkrk2*partik
                  dtxdzi1 = dtdzi1*termx + term*factor
                  dtxidzi1 = dtdzi1*termxi + term*factori
                  dtxkdzi1 = dtdzi1*termxk + term*factork
                  factor = -1.5d0*zkyr - 1.5d0*ykzr + 0.75d0*ziyk
     &                        + 0.75d0*yizk - r2inv*yr*part
                  factori = -ykzr + 1.5d0*zkyr + yiri2*partik
     &                         + enum*(2.0d0*yiri2*ziri2)
                  factork = -yizr - yrzr + 1.5d0*ziyr
     &                         + ykrk2*partik
                  dtydzi1 = dtdzi1*termy + term*factor
                  dtyidzi1 = dtdzi1*termyi + term*factori
                  dtykdzi1 = dtdzi1*termyk + term*factork
                  factor = 0.75d0*dotp - 3.0d0*zkzr + 1.5d0*zizk
     &                        - 1.5d0*dotk - r2inv*(zr*part+dotik)
                  factori = 1.5d0*dotk + 0.5d0*zkzr + ziri2*partik
     &                         - enum*(ri2inv-2.0d0*ziri2*ziri2)
                  factork = r2 + 1.5d0*doti + 0.5d0*zizr
     &                         - zrzr + zkrk2*partik
                  dtzdzi1 = dtdzi1*termz + term*factor
                  dtzidzi1 = dtdzi1*termzi + term*factori
                  dtzkdzi1 = dtdzi1*termzk + term*factork
c
                  dtdxi2 = -2.5d0*xrr2 - xiri2
                  part = xkdoti + 2.0d0*dotk*xr + xidotk
     &                      - 2.0d0*dotik*xrr2
                  partik = xk*r2 + dotp*xr - 1.5d0*xkdoti
     &                        -3.0d0*xr*dotk - 1.5d0*xidotk
                  factor = 0.75d0*dotp + 3.0d0*xkxr + 1.5d0*xixk
     &                        + 1.5d0*dotk - r2inv*(xr*part+dotik)
                  factori = 1.5d0*dotk + 0.5d0*xkxr + xiri2*partik
     &                         + enum*(ri2inv-2.0d0*xiri2*xiri2)
                  factork = -r2 + 1.5d0*doti + 0.5d0*xixr
     &                         + xrxr + xkrk2*partik
                  dtxdxi2 = dtdxi2*termx + term*factor
                  dtxidxi2 = dtdxi2*termxi + term*factori
                  dtxkdxi2 = dtdxi2*termxk + term*factork
                  factor = 1.5d0*xkyr + 1.5d0*ykxr + 0.75d0*xiyk
     &                        + 0.75d0*yixk - r2inv*yr*part
                  factori = -ykxr + 1.5d0*xkyr + yiri2*partik
     &                         - enum*(2.0d0*yiri2*xiri2)
                  factork = -yixr + xryr + 1.5d0*xiyr
     &                         + ykrk2*partik
                  dtydxi2 = dtdxi2*termy + term*factor
                  dtyidxi2 = dtdxi2*termyi + term*factori
                  dtykdxi2 = dtdxi2*termyk + term*factork
                  factor = 1.5d0*xkzr + 1.5d0*zkxr + 0.75d0*xizk
     &                        + 0.75d0*zixk - r2inv*zr*part
                  factori = -zkxr + 1.5d0*xkzr + ziri2*partik
     &                         - enum * (2.0d0*ziri2*xiri2)
                  factork = -zixr + xrzr + 1.5d0*xizr
     &                         + zkrk2*partik
                  dtzdxi2 = dtdxi2*termz + term*factor
                  dtzidxi2 = dtdxi2*termzi + term*factori
                  dtzkdxi2 = dtdxi2*termzk + term*factork
c
                  dtdyi2 = -2.5d0*yrr2 - yiri2
                  part = ykdoti + 2.0d0*dotk*yr + yidotk
     &                      - 2.0d0*dotik*yrr2
                  partik = yk*r2 + dotp*yr - 1.5d0*ykdoti
     &                        -3.0d0*yr*dotk - 1.5d0*yidotk
                  factor = 1.5d0*ykxr + 1.5d0*xkyr + 0.75d0*yixk
     &                        + 0.75d0*xiyk - r2inv*xr*part
                  factori = -xkyr + 1.5d0*ykxr + xiri2*partik
     &                         - enum*(2.0d0*xiri2*yiri2)
                  factork = -xiyr + xryr + 1.5d0*yixr
     &                         + xkrk2*partik
                  dtxdyi2 = dtdyi2*termx + term*factor
                  dtxidyi2 = dtdyi2*termxi + term*factori
                  dtxkdyi2 = dtdyi2*termxk + term*factork
                  factor = 0.75d0*dotp + 3.0d0*ykyr + 1.5d0*yiyk
     &                        + 1.5d0*dotk - r2inv*(yr*part+dotik)
                  factori = 1.5d0*dotk + 0.5d0*ykyr + yiri2*partik
     &                         + enum*(ri2inv-2.0d0*yiri2*yiri2)
                  factork = -r2 + 1.5d0*doti + 0.5d0*yiyr
     &                         + yryr + ykrk2*partik
                  dtydyi2 = dtdyi2*termy + term*factor
                  dtyidyi2 = dtdyi2*termyi + term*factori
                  dtykdyi2 = dtdyi2*termyk + term*factork
                  factor = 1.5d0*ykzr + 1.5d0*zkyr + 0.75d0*yizk
     &                        + 0.75d0*ziyk - r2inv*zr*part
                  factori = -zkyr + 1.5d0*ykzr + ziri2*partik
     &                         - enum * (2.0d0*ziri2*yiri2)
                  factork = -ziyr + yrzr + 1.5d0*yizr
     &                         + zkrk2*partik
                  dtzdyi2 = dtdyi2*termz + term*factor
                  dtzidyi2 = dtdyi2*termzi + term*factori
                  dtzkdyi2 = dtdyi2*termzk + term*factork
c
                  dtdzi2 = -2.5d0*zrr2 - ziri2
                  part = zkdoti + 2.0d0*dotk*zr + zidotk
     &                      - 2.0d0*dotik*zrr2
                  partik = zk*r2 + dotp*zr - 1.5d0*zkdoti
     &                        - 3.0d0*zr*dotk - 1.5d0*zidotk
                  factor = 1.5d0*zkxr + 1.5d0*xkzr + 0.75d0*zixk
     &                        + 0.75d0*xizk - r2inv*xr*part
                  factori = -xkzr + 1.5d0*zkxr + xiri2*partik
     &                         - enum*(2.0d0*xiri2*ziri2)
                  factork = -xizr + xrzr + 1.5d0*zixr
     &                         + xkrk2*partik
                  dtxdzi2 = dtdzi2*termx + term*factor
                  dtxidzi2 = dtdzi2*termxi + term*factori
                  dtxkdzi2 = dtdzi2*termxk + term*factork
                  factor = 1.5d0*zkyr + 1.5d0*ykzr + 0.75d0*ziyk
     &                        + 0.75d0*yizk - r2inv*yr*part
                  factori = -ykzr + 1.5d0*zkyr + yiri2*partik
     &                         - enum*(2.0d0*yiri2*ziri2)
                  factork = -yizr + yrzr + 1.5d0*ziyr
     &                         + ykrk2*partik
                  dtydzi2 = dtdzi2*termy + term*factor
                  dtyidzi2 = dtdzi2*termyi + term*factori
                  dtykdzi2 = dtdzi2*termyk + term*factork
                  factor = 0.75d0*dotp + 3.0d0*zkzr + 1.5d0*zizk
     &                        + 1.5d0*dotk - r2inv*(zr*part+dotik)
                  factori = 1.5d0*dotk + 0.5d0*zkzr + ziri2*partik
     &                         + enum*(ri2inv-2.0d0*ziri2*ziri2)
                  factork = -r2 + 1.5d0*doti + 0.5d0*zizr
     &                         + zrzr + zkrk2*partik
                  dtzdzi2 = dtdzi2*termz + term*factor
                  dtzidzi2 = dtdzi2*termzi + term*factori
                  dtzkdzi2 = dtdzi2*termzk + term*factork
c
                  dtdxk1 =  2.5d0*xrr2 + xkrk2
                  part = xkdoti + 2.0d0*doti*xr + xidotk
     &                      - 2.0d0*dotik*xrr2
                  partik = -xi*r2 - dotp*xr + 1.5d0*xidotk
     &                        + 3.0d0*xr*doti + 1.5d0*xkdoti
                  factor = -0.75d0*dotp - 3.0d0*xixr - 1.5d0*xixk
     &                        - 1.5d0*doti + r2inv*(xr*part+dotik)
                  factori = r2 - 1.5d0*dotk - 0.5d0*xkxr
     &                         - xrxr + xiri2*partik
                  factork = -1.5d0*doti - 0.5d0*xixr + xkrk2*partik
     &                         - enum*(rk2inv-2.0d0*xkrk2*xkrk2)
                  dtxdxk1 = dtdxk1*termx + term*factor
                  dtxidxk1 = dtdxk1*termxi + term*factori
                  dtxkdxk1 = dtdxk1*termxk + term*factork
                  factor = -1.5d0*xiyr - 1.5d0*yixr - 0.75d0*xiyk
     &                        - 0.75d0*yixk + r2inv*yr*part
                  factori = ykxr - xryr - 1.5d0*xkyr
     &                         + yiri2*partik
                  factork = yixr - 1.5d0*xiyr + ykrk2*partik
     &                         + enum*(2.0d0*ykrk2*xkrk2)
                  dtydxk1 = dtdxk1*termy + term*factor
                  dtyidxk1 = dtdxk1*termyi + term*factori
                  dtykdxk1 = dtdxk1*termyk + term*factork
                  factor = -1.5d0*xizr - 1.5d0*zixr - 0.75d0*xizk
     &                        - 0.75d0*zixk + r2inv*zr*part
                  factori = zkxr - xrzr - 1.5d0*xkzr
     &                         + ziri2*partik
                  factork = zixr - 1.5d0*xizr + zkrk2*partik
     &                         + enum*(2.0d0*zkrk2*xkrk2)
                  dtzdxk1 = dtdxk1*termz + term*factor
                  dtzidxk1 = dtdxk1*termzi + term*factori
                  dtzkdxk1 = dtdxk1*termzk + term*factork
c
                  dtdyk1 =  2.5d0*yrr2 + ykrk2
                  part = ykdoti + 2.0d0*doti*yr + yidotk
     &                      - 2.0d0*dotik*yrr2
                  partik = -yi*r2 - dotp*yr + 1.5d0*ykdoti
     &                        + 3.0d0*yr*doti + 1.5d0*yidotk
                  factor = -1.5d0*yixr - 1.5d0*xiyr - 0.75d0*yixk
     &                        - 0.75d0*xiyk + r2inv*xr*part
                  factori = xkyr - xryr - 1.5d0*ykxr
     &                         + xiri2*partik
                  factork = xiyr - 1.5d0*yixr + xkrk2*partik
     &                         + enum*(2.0d0*xkrk2*ykrk2)
                  dtxdyk1 = dtdyk1*termx + term*factor
                  dtxidyk1 = dtdyk1*termxi + term*factori
                  dtxkdyk1 = dtdyk1*termxk + term*factork
                  factor = -0.75d0*dotp - 3.0d0*yiyr - 1.5d0*yiyk
     &                        - 1.5d0*doti + r2inv*(yr*part+dotik)
                  factori = r2 - 1.5d0*dotk - 0.5d0*ykyr
     &                         - yryr + yiri2*partik
                  factork = -1.5d0*doti - 0.5d0*yiyr + ykrk2*partik
     &                         - enum*(rk2inv-2.0d0*ykrk2*ykrk2)
                  dtydyk1 = dtdyk1*termy + term*factor
                  dtyidyk1 = dtdyk1*termyi + term*factori
                  dtykdyk1 = dtdyk1*termyk + term*factork
                  factor = -1.5d0*yizr - 1.5d0*ziyr - 0.75d0*yizk
     &                        - 0.75d0*ziyk + r2inv*zr*part
                  factori = zkyr - yrzr - 1.5d0*ykzr
     &                         + ziri2*partik
                  factork = ziyr - 1.5d0*yizr + zkrk2*partik
     &                         + enum*(2.0d0*zkrk2*ykrk2)
                  dtzdyk1 = dtdyk1*termz + term*factor
                  dtzidyk1 = dtdyk1*termzi + term*factori
                  dtzkdyk1 = dtdyk1*termzk + term*factork
c
                  dtdzk1 =  2.5d0*zrr2 + zkrk2
                  part = zkdoti + 2.0d0*doti*zr + zidotk
     &                      - 2.0d0*dotik*zrr2
                  partik = -zi*r2 - dotp*zr + 1.5d0*zkdoti
     &                        + 3.0d0*zr*doti + 1.5d0*zidotk
                  factor = -1.5d0*zixr - 1.5d0*xizr - 0.75d0*zixk
     &                        - 0.75d0*xizk + r2inv*xr*part
                  factori = xkzr - xrzr - 1.5d0*zkxr
     &                         + xiri2*partik
                  factork = xizr - 1.5d0*zixr + xkrk2*partik
     &                         + enum*(2.0d0*xkrk2*zkrk2)
                  dtxdzk1 = dtdzk1*termx + term*factor
                  dtxidzk1 = dtdzk1*termxi + term*factori
                  dtxkdzk1 = dtdzk1*termxk + term*factork
                  factor = -1.5d0*ziyr - 1.5d0*yizr - 0.75d0*ziyk
     &                        - 0.75d0*yizk + r2inv*yr*part
                  factori = ykzr - yrzr - 1.5d0*zkyr
     &                         + yiri2*partik
                  factork = yizr - 1.5d0*ziyr + ykrk2*partik
     &                         + enum*(2.0d0*ykrk2*zkrk2)
                  dtydzk1 = dtdzk1*termy + term*factor
                  dtyidzk1 = dtdzk1*termyi + term*factori
                  dtykdzk1 = dtdzk1*termyk + term*factork
                  factor = -0.75d0*dotp - 3.0d0*zizr - 1.5d0*zizk
     &                        - 1.5d0*doti + r2inv*(zr*part+dotik)
                  factori = r2 - 1.5d0*dotk - 0.5d0*zkzr
     &                         - zrzr + ziri2*partik
                  factork = -1.5d0*doti - 0.5d0*zizr + zkrk2*partik
     &                         - enum*(rk2inv-2.0d0*zkrk2*zkrk2)
                  dtzdzk1 = dtdzk1*termz + term*factor
                  dtzidzk1 = dtdzk1*termzi + term*factori
                  dtzkdzk1 = dtdzk1*termzk + term*factork
c
                  dtdxk2 =  2.5d0*xrr2 - xkrk2
                  part = xkdoti - 2.0d0*doti*xr + xidotk
     &                      - 2.0d0*dotik*xrr2
                  partik = xi*r2 - dotp*xr + 1.5d0*xidotk
     &                        - 3.0d0*xr*doti + 1.5d0*xkdoti
                  factor = -0.75d0*dotp + 3.0d0*xixr - 1.5d0*xixk
     &                        + 1.5d0*doti + r2inv*(xr*part+dotik)
                  factori = -r2 - 1.5d0*dotk - 0.5d0*xkxr
     &                         + xrxr + xiri2*partik
                  factork = -1.5d0*doti - 0.5d0*xixr + xkrk2*partik
     &                         + enum*(rk2inv-2.0d0*xkrk2*xkrk2)
                  dtxdxk2 = dtdxk2*termx + term*factor
                  dtxidxk2 = dtdxk2*termxi + term*factori
                  dtxkdxk2 = dtdxk2*termxk + term*factork
                  factor = 1.5d0*xiyr + 1.5d0*yixr - 0.75d0*xiyk
     &                        - 0.75d0*yixk + r2inv*yr*part
                  factori = ykxr + xryr - 1.5d0*xkyr
     &                         + yiri2*partik
                  factork = yixr - 1.5d0*xiyr + ykrk2*partik
     &                         - enum*(2.0d0*ykrk2*xkrk2)
                  dtydxk2 = dtdxk2*termy + term*factor
                  dtyidxk2 = dtdxk2*termyi + term*factori
                  dtykdxk2 = dtdxk2*termyk + term*factork
                  factor = 1.5d0*xizr + 1.5d0*zixr - 0.75d0*xizk
     &                        - 0.75d0*zixk + r2inv*zr*part
                  factori = zkxr + xrzr - 1.5d0*xkzr
     &                         + ziri2*partik
                  factork = zixr - 1.5d0*xizr + zkrk2*partik
     &                         - enum*(2.0d0*zkrk2*xkrk2)
                  dtzdxk2 = dtdxk2*termz + term*factor
                  dtzidxk2 = dtdxk2*termzi + term*factori
                  dtzkdxk2 = dtdxk2*termzk + term*factork
c
                  dtdyk2 =  2.5d0*yrr2 - ykrk2
                  part = ykdoti - 2.0d0*doti*yr + yidotk
     &                      - 2.0d0*dotik*yrr2
                  partik = yi*r2 - dotp*yr + 1.5d0*ykdoti
     &                        - 3.0d0*yr*doti + 1.5d0*yidotk
                  factor = 1.5d0*yixr + 1.5d0*xiyr - 0.75d0*yixk
     &                        - 0.75d0*xiyk + r2inv*xr*part
                  factori = xkyr + xryr - 1.5d0*ykxr
     &                         + xiri2*partik
                  factork = xiyr - 1.5d0*yixr + xkrk2*partik
     &                         - enum*(2.0d0*xkrk2*ykrk2)
                  dtxdyk2 = dtdyk2*termx + term*factor
                  dtxidyk2 = dtdyk2*termxi + term*factori
                  dtxkdyk2 = dtdyk2*termxk + term*factork
                  factor = -0.75d0*dotp + 3.0d0*yiyr - 1.5d0*yiyk
     &                        + 1.5d0*doti + r2inv*(yr*part+dotik)
                  factori = -r2 - 1.5d0*dotk - 0.5d0*ykyr
     &                         + yryr + yiri2*partik
                  factork = -1.5d0*doti - 0.5d0*yiyr + ykrk2*partik
     &                         + enum*(rk2inv-2.0d0*ykrk2*ykrk2)
                  dtydyk2 = dtdyk2*termy + term*factor
                  dtyidyk2 = dtdyk2*termyi + term*factori
                  dtykdyk2 = dtdyk2*termyk + term*factork
                  factor = 1.5d0*yizr + 1.5d0*ziyr - 0.75d0*yizk
     &                        - 0.75d0*ziyk + r2inv*zr*part
                  factori = zkyr + yrzr - 1.5d0*ykzr
     &                         + ziri2*partik
                  factork = ziyr - 1.5d0*yizr + zkrk2*partik
     &                         - enum*(2.0d0*zkrk2*ykrk2)
                  dtzdyk2 = dtdyk2*termz + term*factor
                  dtzidyk2 = dtdyk2*termzi + term*factori
                  dtzkdyk2 = dtdyk2*termzk + term*factork
c
                  dtdzk2 =  2.5d0*zrr2 - zkrk2
                  part = zkdoti - 2.0d0*doti*zr + zidotk
     &                      - 2.0d0*dotik*zrr2
                  partik = zi*r2 - dotp*zr + 1.5d0*zkdoti
     &                        - 3.0d0*zr*doti + 1.5d0*zidotk
                  factor = 1.5d0*zixr + 1.5d0*xizr - 0.75d0*zixk
     &                        - 0.75d0*xizk + r2inv*xr*part
                  factori = xkzr + xrzr - 1.5d0*zkxr
     &                         + xiri2*partik
                  factork = xizr - 1.5d0*zixr + xkrk2*partik
     &                         - enum*(2.0d0*xkrk2*zkrk2)
                  dtxdzk2 = dtdzk2*termx + term*factor
                  dtxidzk2 = dtdzk2*termxi + term*factori
                  dtxkdzk2 = dtdzk2*termxk + term*factork
                  factor = 1.5d0*ziyr + 1.5d0*yizr - 0.75d0*ziyk
     &                        - 0.75d0*yizk + r2inv*yr*part
                  factori = ykzr + yrzr - 1.5d0*zkyr
     &                         + yiri2*partik
                  factork = yizr - 1.5d0*ziyr + ykrk2*partik
     &                         - enum*(2.0d0*ykrk2*zkrk2)
                  dtydzk2 = dtdzk2*termy + term*factor
                  dtyidzk2 = dtdzk2*termyi + term*factori
                  dtykdzk2 = dtdzk2*termyk + term*factork
                  factor = -0.75d0*dotp + 3.0d0*zizr - 1.5d0*zizk
     &                        + 1.5d0*doti + r2inv*(zr*part+dotik)
                  factori = -r2 - 1.5d0*dotk - 0.5d0*zkzr
     &                         + zrzr + ziri2*partik
                  factork = -1.5d0*doti - 0.5d0*zizr + zkrk2*partik
     &                         + enum*(rk2inv-2.0d0*zkrk2*zkrk2)
                  dtzdzk2 = dtdzk2*termz + term*factor
                  dtzidzk2 = dtdzk2*termzi + term*factori
                  dtzkdzk2 = dtdzk2*termzk + term*factork
c
c     now, increment the diagonal hessian block elements
c
                  hed(1,i1,1,i1) = hed(1,i1,1,i1) + dtxdxi1 - dtxidxi1
                  hed(1,i1,2,i1) = hed(1,i1,2,i1) + dtydxi1 - dtyidxi1
                  hed(2,i1,2,i1) = hed(2,i1,2,i1) + dtydyi1 - dtyidyi1
                  hed(1,i1,3,i1) = hed(1,i1,3,i1) + dtzdxi1 - dtzidxi1
                  hed(2,i1,3,i1) = hed(2,i1,3,i1) + dtzdyi1 - dtzidyi1
                  hed(3,i1,3,i1) = hed(3,i1,3,i1) + dtzdzi1 - dtzidzi1
                  hed(1,i2,1,i2) = hed(1,i2,1,i2) + dtxdxi2 + dtxidxi2
                  hed(1,i2,2,i2) = hed(1,i2,2,i2) + dtydxi2 + dtyidxi2
                  hed(2,i2,2,i2) = hed(2,i2,2,i2) + dtydyi2 + dtyidyi2
                  hed(1,i2,3,i2) = hed(1,i2,3,i2) + dtzdxi2 + dtzidxi2
                  hed(2,i2,3,i2) = hed(2,i2,3,i2) + dtzdyi2 + dtzidyi2
                  hed(3,i2,3,i2) = hed(3,i2,3,i2) + dtzdzi2 + dtzidzi2
                  hed(1,k1,1,k1) = hed(1,k1,1,k1) - dtxdxk1 - dtxkdxk1
                  hed(1,k1,2,k1) = hed(1,k1,2,k1) - dtydxk1 - dtykdxk1
                  hed(2,k1,2,k1) = hed(2,k1,2,k1) - dtydyk1 - dtykdyk1
                  hed(1,k1,3,k1) = hed(1,k1,3,k1) - dtzdxk1 - dtzkdxk1
                  hed(2,k1,3,k1) = hed(2,k1,3,k1) - dtzdyk1 - dtzkdyk1
                  hed(3,k1,3,k1) = hed(3,k1,3,k1) - dtzdzk1 - dtzkdzk1
                  hed(1,k2,1,k2) = hed(1,k2,1,k2) - dtxdxk2 + dtxkdxk2
                  hed(1,k2,2,k2) = hed(1,k2,2,k2) - dtydxk2 + dtykdxk2
                  hed(2,k2,2,k2) = hed(2,k2,2,k2) - dtydyk2 + dtykdyk2
                  hed(1,k2,3,k2) = hed(1,k2,3,k2) - dtzdxk2 + dtzkdxk2
                  hed(2,k2,3,k2) = hed(2,k2,3,k2) - dtzdyk2 + dtzkdyk2
                  hed(3,k2,3,k2) = hed(3,k2,3,k2) - dtzdzk2 + dtzkdzk2
c
c     finally, increment above-diagonal hessian block elements
c
                  if (i1 .lt. i2) then
                     hed(1,i1,1,i2) = hed(1,i1,1,i2) + dtxdxi1+dtxidxi1
                     hed(2,i1,1,i2) = hed(2,i1,1,i2) + dtxdyi1+dtxidyi1
                     hed(3,i1,1,i2) = hed(3,i1,1,i2) + dtxdzi1+dtxidzi1
                     hed(1,i1,2,i2) = hed(1,i1,2,i2) + dtydxi1+dtyidxi1
                     hed(2,i1,2,i2) = hed(2,i1,2,i2) + dtydyi1+dtyidyi1
                     hed(3,i1,2,i2) = hed(3,i1,2,i2) + dtydzi1+dtyidzi1
                     hed(1,i1,3,i2) = hed(1,i1,3,i2) + dtzdxi1+dtzidxi1
                     hed(2,i1,3,i2) = hed(2,i1,3,i2) + dtzdyi1+dtzidyi1
                     hed(3,i1,3,i2) = hed(3,i1,3,i2) + dtzdzi1+dtzidzi1
                  else
                     hed(1,i2,1,i1) = hed(1,i2,1,i1) + dtxdxi2-dtxidxi2
                     hed(2,i2,1,i1) = hed(2,i2,1,i1) + dtxdyi2-dtxidyi2
                     hed(3,i2,1,i1) = hed(3,i2,1,i1) + dtxdzi2-dtxidzi2
                     hed(1,i2,2,i1) = hed(1,i2,2,i1) + dtydxi2-dtyidxi2
                     hed(2,i2,2,i1) = hed(2,i2,2,i1) + dtydyi2-dtyidyi2
                     hed(3,i2,2,i1) = hed(3,i2,2,i1) + dtydzi2-dtyidzi2
                     hed(1,i2,3,i1) = hed(1,i2,3,i1) + dtzdxi2-dtzidxi2
                     hed(2,i2,3,i1) = hed(2,i2,3,i1) + dtzdyi2-dtzidyi2
                     hed(3,i2,3,i1) = hed(3,i2,3,i1) + dtzdzi2-dtzidzi2
                  end if
c
                  if (i1 .lt. k1) then
                     hed(1,i1,1,k1) = hed(1,i1,1,k1) - dtxdxi1-dtxkdxi1
                     hed(2,i1,1,k1) = hed(2,i1,1,k1) - dtxdyi1-dtxkdyi1
                     hed(3,i1,1,k1) = hed(3,i1,1,k1) - dtxdzi1-dtxkdzi1
                     hed(1,i1,2,k1) = hed(1,i1,2,k1) - dtydxi1-dtykdxi1
                     hed(2,i1,2,k1) = hed(2,i1,2,k1) - dtydyi1-dtykdyi1
                     hed(3,i1,2,k1) = hed(3,i1,2,k1) - dtydzi1-dtykdzi1
                     hed(1,i1,3,k1) = hed(1,i1,3,k1) - dtzdxi1-dtzkdxi1
                     hed(2,i1,3,k1) = hed(2,i1,3,k1) - dtzdyi1-dtzkdyi1
                     hed(3,i1,3,k1) = hed(3,i1,3,k1) - dtzdzi1-dtzkdzi1
                  else
                     hed(1,k1,1,i1) = hed(1,k1,1,i1) + dtxdxk1-dtxidxk1
                     hed(2,k1,1,i1) = hed(2,k1,1,i1) + dtxdyk1-dtxidyk1
                     hed(3,k1,1,i1) = hed(3,k1,1,i1) + dtxdzk1-dtxidzk1
                     hed(1,k1,2,i1) = hed(1,k1,2,i1) + dtydxk1-dtyidxk1
                     hed(2,k1,2,i1) = hed(2,k1,2,i1) + dtydyk1-dtyidyk1
                     hed(3,k1,2,i1) = hed(3,k1,2,i1) + dtydzk1-dtyidzk1
                     hed(1,k1,3,i1) = hed(1,k1,3,i1) + dtzdxk1-dtzidxk1
                     hed(2,k1,3,i1) = hed(2,k1,3,i1) + dtzdyk1-dtzidyk1
                     hed(3,k1,3,i1) = hed(3,k1,3,i1) + dtzdzk1-dtzidzk1
                  end if
c
                  if (i1 .lt. k2) then
                     hed(1,i1,1,k2) = hed(1,i1,1,k2) - dtxdxi1+dtxkdxi1
                     hed(2,i1,1,k2) = hed(2,i1,1,k2) - dtxdyi1+dtxkdyi1
                     hed(3,i1,1,k2) = hed(3,i1,1,k2) - dtxdzi1+dtxkdzi1
                     hed(1,i1,2,k2) = hed(1,i1,2,k2) - dtydxi1+dtykdxi1
                     hed(2,i1,2,k2) = hed(2,i1,2,k2) - dtydyi1+dtykdyi1
                     hed(3,i1,2,k2) = hed(3,i1,2,k2) - dtydzi1+dtykdzi1
                     hed(1,i1,3,k2) = hed(1,i1,3,k2) - dtzdxi1+dtzkdxi1
                     hed(2,i1,3,k2) = hed(2,i1,3,k2) - dtzdyi1+dtzkdyi1
                     hed(3,i1,3,k2) = hed(3,i1,3,k2) - dtzdzi1+dtzkdzi1
                  else
                     hed(1,k2,1,i1) = hed(1,k2,1,i1) + dtxdxk2-dtxidxk2
                     hed(2,k2,1,i1) = hed(2,k2,1,i1) + dtxdyk2-dtxidyk2
                     hed(3,k2,1,i1) = hed(3,k2,1,i1) + dtxdzk2-dtxidzk2
                     hed(1,k2,2,i1) = hed(1,k2,2,i1) + dtydxk2-dtyidxk2
                     hed(2,k2,2,i1) = hed(2,k2,2,i1) + dtydyk2-dtyidyk2
                     hed(3,k2,2,i1) = hed(3,k2,2,i1) + dtydzk2-dtyidzk2
                     hed(1,k2,3,i1) = hed(1,k2,3,i1) + dtzdxk2-dtzidxk2
                     hed(2,k2,3,i1) = hed(2,k2,3,i1) + dtzdyk2-dtzidyk2
                     hed(3,k2,3,i1) = hed(3,k2,3,i1) + dtzdzk2-dtzidzk2
                  end if
c
                  if (i2 .lt. k1) then
                     hed(1,i2,1,k1) = hed(1,i2,1,k1) - dtxdxi2-dtxkdxi2
                     hed(2,i2,1,k1) = hed(2,i2,1,k1) - dtxdyi2-dtxkdyi2
                     hed(3,i2,1,k1) = hed(3,i2,1,k1) - dtxdzi2-dtxkdzi2
                     hed(1,i2,2,k1) = hed(1,i2,2,k1) - dtydxi2-dtykdxi2
                     hed(2,i2,2,k1) = hed(2,i2,2,k1) - dtydyi2-dtykdyi2
                     hed(3,i2,2,k1) = hed(3,i2,2,k1) - dtydzi2-dtykdzi2
                     hed(1,i2,3,k1) = hed(1,i2,3,k1) - dtzdxi2-dtzkdxi2
                     hed(2,i2,3,k1) = hed(2,i2,3,k1) - dtzdyi2-dtzkdyi2
                     hed(3,i2,3,k1) = hed(3,i2,3,k1) - dtzdzi2-dtzkdzi2
                  else
                     hed(1,k1,1,i2) = hed(1,k1,1,i2) + dtxdxk1+dtxidxk1
                     hed(2,k1,1,i2) = hed(2,k1,1,i2) + dtxdyk1+dtxidyk1
                     hed(3,k1,1,i2) = hed(3,k1,1,i2) + dtxdzk1+dtxidzk1
                     hed(1,k1,2,i2) = hed(1,k1,2,i2) + dtydxk1+dtyidxk1
                     hed(2,k1,2,i2) = hed(2,k1,2,i2) + dtydyk1+dtyidyk1
                     hed(3,k1,2,i2) = hed(3,k1,2,i2) + dtydzk1+dtyidzk1
                     hed(1,k1,3,i2) = hed(1,k1,3,i2) + dtzdxk1+dtzidxk1
                     hed(2,k1,3,i2) = hed(2,k1,3,i2) + dtzdyk1+dtzidyk1
                     hed(3,k1,3,i2) = hed(3,k1,3,i2) + dtzdzk1+dtzidzk1
                  end if
c
                  if (i2 .lt. k2) then
                     hed(1,i2,1,k2) = hed(1,i2,1,k2) - dtxdxi2+dtxkdxi2
                     hed(2,i2,1,k2) = hed(2,i2,1,k2) - dtxdyi2+dtxkdyi2
                     hed(3,i2,1,k2) = hed(3,i2,1,k2) - dtxdzi2+dtxkdzi2
                     hed(1,i2,2,k2) = hed(1,i2,2,k2) - dtydxi2+dtykdxi2
                     hed(2,i2,2,k2) = hed(2,i2,2,k2) - dtydyi2+dtykdyi2
                     hed(3,i2,2,k2) = hed(3,i2,2,k2) - dtydzi2+dtykdzi2
                     hed(1,i2,3,k2) = hed(1,i2,3,k2) - dtzdxi2+dtzkdxi2
                     hed(2,i2,3,k2) = hed(2,i2,3,k2) - dtzdyi2+dtzkdyi2
                     hed(3,i2,3,k2) = hed(3,i2,3,k2) - dtzdzi2+dtzkdzi2
                  else
                     hed(1,k2,1,i2) = hed(1,k2,1,i2) + dtxdxk2+dtxidxk2
                     hed(2,k2,1,i2) = hed(2,k2,1,i2) + dtxdyk2+dtxidyk2
                     hed(3,k2,1,i2) = hed(3,k2,1,i2) + dtxdzk2+dtxidzk2
                     hed(1,k2,2,i2) = hed(1,k2,2,i2) + dtydxk2+dtyidxk2
                     hed(2,k2,2,i2) = hed(2,k2,2,i2) + dtydyk2+dtyidyk2
                     hed(3,k2,2,i2) = hed(3,k2,2,i2) + dtydzk2+dtyidzk2
                     hed(1,k2,3,i2) = hed(1,k2,3,i2) + dtzdxk2+dtzidxk2
                     hed(2,k2,3,i2) = hed(2,k2,3,i2) + dtzdyk2+dtzidyk2
                     hed(3,k2,3,i2) = hed(3,k2,3,i2) + dtzdzk2+dtzidzk2
                  end if
c
                  if (k1 .lt. k2) then
                     hed(1,k1,1,k2) = hed(1,k1,1,k2) - dtxdxk1+dtxkdxk1
                     hed(2,k1,1,k2) = hed(2,k1,1,k2) - dtxdyk1+dtxkdyk1
                     hed(3,k1,1,k2) = hed(3,k1,1,k2) - dtxdzk1+dtxkdzk1
                     hed(1,k1,2,k2) = hed(1,k1,2,k2) - dtydxk1+dtykdxk1
                     hed(2,k1,2,k2) = hed(2,k1,2,k2) - dtydyk1+dtykdyk1
                     hed(3,k1,2,k2) = hed(3,k1,2,k2) - dtydzk1+dtykdzk1
                     hed(1,k1,3,k2) = hed(1,k1,3,k2) - dtzdxk1+dtzkdxk1
                     hed(2,k1,3,k2) = hed(2,k1,3,k2) - dtzdyk1+dtzkdyk1
                     hed(3,k1,3,k2) = hed(3,k1,3,k2) - dtzdzk1+dtzkdzk1
                  else
                     hed(1,k2,1,k1) = hed(1,k2,1,k1) - dtxdxk2-dtxkdxk2
                     hed(2,k2,1,k1) = hed(2,k2,1,k1) - dtxdyk2-dtxkdyk2
                     hed(3,k2,1,k1) = hed(3,k2,1,k1) - dtxdzk2-dtxkdzk2
                     hed(1,k2,2,k1) = hed(1,k2,2,k1) - dtydxk2-dtykdxk2
                     hed(2,k2,2,k1) = hed(2,k2,2,k1) - dtydyk2-dtykdyk2
                     hed(3,k2,2,k1) = hed(3,k2,2,k1) - dtydzk2-dtykdzk2
                     hed(1,k2,3,k1) = hed(1,k2,3,k1) - dtzdxk2-dtzkdxk2
                     hed(2,k2,3,k1) = hed(2,k2,3,k1) - dtzdyk2-dtzkdyk2
                     hed(3,k2,3,k1) = hed(3,k2,3,k1) - dtzdzk2-dtzkdzk2
                  end if
               end if
            end if
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine extra5  --  user defined extra potentials  ##
c     ##                                                        ##
c     ############################################################
c
c
      subroutine extra5
      implicit none
      include 'sizes.for'
      integer i,j,k,l,n
      real*8 heb,hea,heba,het,hev,he14,hec,hecd,hed,hex
      common /atoms / n
      common /hessn1/ heb(3,maxhes,3,maxhes),hea(3,maxhes,3,maxhes),
     &                heba(3,maxhes,3,maxhes),het(3,maxhes,3,maxhes),
     &                hev(3,maxhes,3,maxhes),he14(3,maxhes,3,maxhes),
     &                hec(3,maxhes,3,maxhes),hecd(3,maxhes,3,maxhes),
     &                hed(3,maxhes,3,maxhes),hex(3,maxhes,3,maxhes)
c
c
c     zero out the extra energy term hessian matrix elements
c
      do l = 1, n
         do k = 1, 3
            do j = 1, n
               do i = 1, 3
                  hex(i,j,k,l) = 0.0d0
               end do
            end do
         end do
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine eangle5  --  bend & str-bend hessian; cartesian  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "eangle5" calculates second derivatives (hessian) of the bend
c     and stretch bend energy with respect to the cartesian coordinates
c
c
      subroutine eangle5
      implicit none
      include 'sizes.for'
      integer n,itype,n12,i12,nbond,ibnd
      integer nangle,nstbnd,iang,ksbopb,jopb
      integer i,j,k,l,m,ii,kk,kang,ia,ib,ic,id,nang,nsb,nopb
      real*8 x,y,z,xi,yi,zi,xp,yp,zp,bx,by,bz,cx,cy,cz
      real*8 bk,bl,anat,acon,sbk,opbk,sf,radian
      real*8 dt,dt2,dq,dq2,dr,a,b,c,d,s
      real*8 ang(6),angin(6),angout(6),dot(6),norm(6),cos(6)
      real*8 dx(4),dy(4),dz(4),ds(4),dxx(3),dyy(3),dzz(3),dss(3)
      real*8 doti(3),normi(3),cosi(3),doto(3),normo(3),coso(3)
      real*8 ddtdxi,ddtdyi,ddtdzi,ddrdxi,ddrdyi,ddrdzi
      real*8 ddtdxia,ddtdyia,ddtdzia,ddtdxib,ddtdyib,ddtdzib
      real*8 ddrdxia,ddrdyia,ddrdzia,ddrdxib,ddrdyib,ddrdzib
      real*8 dcosodxi,dcosodyi,dcosodzi,dcosidxia,dcosidyia
      real*8 dcosidzia,dcosidxib,dcosidyib,dcosidzib
      real*8 dcosidxic,dcosidyic,dcosidzic
      real*8 dcosidxi,dcosidyi,dcosidzi,dcosodxia,dcosodyia
      real*8 dcosodzia,dcosodxib,dcosodyib,dcosodzib
      real*8 dcosodxic,dcosodyic,dcosodzic
      real*8 ddxxdxi,ddxxdyi,ddxxdzi,ddyydxi,ddyydyi
      real*8 ddyydzi,ddzzdxi,ddzzdyi,ddzzdzi
      real*8 ddyydxia,ddzzdxia,ddxxdyia,ddzzdyia,ddxxdzia
      real*8 ddyydzia,ddyydxib,ddzzdxib,ddxxdyib,ddzzdyib
      real*8 ddxxdzib,ddyydzib,ddyydxic,ddzzdxic
      real*8 ddxxdyic,ddzzdyic,ddxxdzic,ddyydzic
      real*8 dbdxia,dcdxia,dddxia,dsdxia,dadyia,dcdyia
      real*8 dddyia,dsdyia,dadzia,dbdzia,dddzia,dsdzia
      real*8 dbdxib,dcdxib,dddxib,dsdxib,dadyib,dcdyib
      real*8 dddyib,dsdyib,dadzib,dbdzib,dddzib,dsdzib
      real*8 dbdxic,dcdxic,dddxic,dsdxic,dadyic,dcdyic
      real*8 dddyic,dsdyic,dadzic,dbdzic,dddzic,dsdzic
      real*8 ddsdxi(3),ddsdyi(3),ddsdzi(3)
      real*8 ddssdxi(3),ddssdyi(3),ddssdzi(3)
      real*8 dxxss(3),dyyss(3),dzzss(3)
      real*8 ddxdxia(3),ddxxdxia(3),ddsdxia(3),ddssdxia(3)
      real*8 ddydyia(3),ddyydyia(3),ddsdyia(3),ddssdyia(3)
      real*8 ddzdzia(3),ddzzdzia(3),ddsdzia(3),ddssdzia(3)
      real*8 ddxdxib(3),ddxxdxib(3),ddsdxib(3),ddssdxib(3)
      real*8 ddydyib(3),ddyydyib(3),ddsdyib(3),ddssdyib(3)
      real*8 ddzdzib(3),ddzzdzib(3),ddsdzib(3),ddssdzib(3)
      real*8 ddxdxic(3),ddxxdxic(3),ddsdxic(3),ddssdxic(3)
      real*8 ddydyic(3),ddyydyic(3),ddsdyic(3),ddssdyic(3)
      real*8 ddzdzic(3),ddzzdzic(3),ddsdzic(3),ddssdzic(3)
      real*8 deddt,ddtdcos,term,terma,termb,sine
      real*8 deoddq,deodcoso,ddqdcoso,deiddt,deidcosi,ddtdcosi
c
      real*8 d2dtdcos2,d2eddt2,cotan
      real*8 d2dqcoso2,d2eoddq2,d2eodcoso2
      real*8 d2dtcosi2,d2eiddt2,d2eidcosi2
      real*8 dtxixi,dtxiyi,dtxizi,dtyiyi,dtyizi,dtzizi
      real*8 dtxiaxia,dtxiayia,dtxiazia,dtyiayia,dtyiazia,dtziazia
      real*8 dtxibxib,dtxibyib,dtxibzib,dtyibyib,dtyibzib,dtzibzib
      real*8 dtxixia,dtxiyia,dtxizia,dtyixia,dtyiyia
      real*8 dtyizia,dtzixia,dtziyia,dtzizia
      real*8 dtxixib,dtxiyib,dtxizib,dtyixib,dtyiyib
      real*8 dtyizib,dtzixib,dtziyib,dtzizib
      real*8 dtxiaxib,dtxiayib,dtxiazib,dtyiaxib,dtyiayib
      real*8 dtyiazib,dtziaxib,dtziayib,dtziazib
      real*8 drxixi,drxiyi,drxizi,dryiyi,dryizi,drzizi
      real*8 drxiaxia,drxiayia,drxiazia,dryiayia,dryiazia,drziazia
      real*8 drxibxib,drxibyib,drxibzib,dryibyib,dryibzib,drzibzib
      real*8 heb,hea,heba,het,hev,he14,hec,hecd,hed
      common /angle / nangle,nstbnd,iang(4,maxang),anat(maxang),
     &                acon(maxang),ksbopb(4,maxang),sbk(4),opbk(60),
     &                sf,jopb(maxtyp)
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm),itype(maxatm)
      common /attach/ n12(maxatm),i12(4,maxatm)
      common /bond  / nbond,ibnd(2,maxbnd),bk(maxbnd),bl(maxbnd)
      common /hessn1/ heb(3,maxhes,3,maxhes),hea(3,maxhes,3,maxhes),
     &                heba(3,maxhes,3,maxhes),het(3,maxhes,3,maxhes),
     &                hev(3,maxhes,3,maxhes),he14(3,maxhes,3,maxhes),
     &                hec(3,maxhes,3,maxhes),hecd(3,maxhes,3,maxhes),
     &                hed(3,maxhes,3,maxhes)
c
c
c     zero out the bend and stretch-bend energy hessian elements
c
      do l = 1, n
         do k = 1, 3
            do j = 1, n
               do i = 1, 3
                  hea(i,j,k,l) = 0.0d0
                  heba(i,j,k,l) = 0.0d0
               end do
            end do
         end do
      end do
      if (nangle .eq. 0)  return
      radian = 180.0d0 / acos(-1.0d0)
      kang = 0
c
c     loop over all atoms, calculating the energy for angles
c     centered on atom "i" before going to the next atom;
c     "nang" is the number of angles centered on atom "i"
c
      do i = 1, n
         if (n12(i) .ge. 2) then
            ia = i12(1,i)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            dx(1) = x(ia) - xi
            dy(1) = y(ia) - yi
            dz(1) = z(ia) - zi
            ds(1) = dsqrt(dx(1)*dx(1)+dy(1)*dy(1)+dz(1)*dz(1))
            ib = i12(2,i)
            dx(2) = x(ib) - xi
            dy(2) = y(ib) - yi
            dz(2) = z(ib) - zi
            ds(2) = dsqrt(dx(2)*dx(2)+dy(2)*dy(2)+dz(2)*dz(2))
            dot(1) = dx(1)*dx(2)+dy(1)*dy(2)+dz(1)*dz(2)
            norm(1) = ds(1) * ds(2)
            cos(1) = dot(1) / norm(1)
            if (abs(cos(1)) .gt. 1.0d0)
     &         cos(1) = cos(1)/abs(cos(1))
            ang(1) = radian * acos(cos(1))
            nang = 1
            if (n12(i) .ne. 2) then
               ic = i12(3,i)
               dx(3) = x(ic) - xi
               dy(3) = y(ic) - yi
               dz(3) = z(ic) - zi
               ds(3) = dsqrt(dx(3)*dx(3)+dy(3)*dy(3)+dz(3)*dz(3))
               dot(2) = dx(1)*dx(3)+dy(1)*dy(3)+dz(1)*dz(3)
               norm(2) = ds(1) * ds(3)
               cos(2) = dot(2) / norm(2)
               if (abs(cos(2)) .gt. 1.0d0)
     &            cos(2) = cos(2)/abs(cos(2))
               ang(2) = radian * acos(cos(2))
               dot(3) = dx(2)*dx(3)+dy(2)*dy(3)+dz(2)*dz(3)
               norm(3) = ds(2) * ds(3)
               cos(3) = dot(3) / norm(3)
               if (abs(cos(3)) .gt. 1.0d0)
     &            cos(3) = cos(3)/abs(cos(3))
               ang(3) = radian * acos(cos(3))
               nang = 3
               if (n12(i) .ne. 3) then
                  id = i12(4,i)
                  dx(4) = x(id) - xi
                  dy(4) = y(id) - yi
                  dz(4) = z(id) - zi
                  ds(4) = dsqrt(dx(4)*dx(4)+dy(4)*dy(4)+dz(4)*dz(4))
                  dot(4) = dot(3)
                  norm(4) = norm(3)
                  cos(4) = cos(3)
                  ang(4) = ang(3)
                  dot(3) = dx(1)*dx(4)+dy(1)*dy(4)+dz(1)*dz(4)
                  norm(3) = ds(1) * ds(4)
                  cos(3) = dot(3) / norm(3)
                  if (abs(cos(3)) .gt. 1.0d0)
     &               cos(3) = cos(3)/abs(cos(3))
                  ang(3) = radian * acos(cos(3))
                  dot(5) = dx(2)*dx(4)+dy(2)*dy(4)+dz(2)*dz(4)
                  norm(5) = ds(2) * ds(4)
                  cos(5) = dot(5) / norm(5)
                  if (abs(cos(5)) .gt. 1.0d0)
     &               cos(5) = cos(5)/abs(cos(5))
                  ang(5) = radian * acos(cos(5))
                  dot(6) = dx(3)*dx(4)+dy(3)*dy(4)+dz(3)*dz(4)
                  norm(6) = ds(3) * ds(4)
                  cos(6) = dot(6) / norm(6)
                  if (abs(cos(6)) .gt. 1.0d0)
     &               cos(6) = cos(6)/abs(cos(6))
                  ang(6) = radian * acos(cos(6))
                  nang = 6
               else
c
c     this section calculates in-plane and out-of-plane angles;
c     types with "jopb" of 1 are used in out-of-plane bending
c
                  if (jopb(itype(i)) .eq. 1) then
                     bx = dx(2) - dx(1)
                     by = dy(2) - dy(1)
                     bz = dz(2) - dz(1)
                     cx = dx(3) - dx(1)
                     cy = dy(3) - dy(1)
                     cz = dz(3) - dz(1)
                     a = by*cz - cy*bz
                     b = cx*bz - bx*cz
                     c = bx*cy - cx*by
                     s = a*a + b*b + c*c
                     d = a*dx(1) + b*dy(1) + c*dz(1)
                     xp = a*d/s
                     yp = b*d/s
                     zp = c*d/s
                     dxx(1) = dx(1) - xp
                     dyy(1) = dy(1) - yp
                     dzz(1) = dz(1) - zp
                     dxx(2) = dx(2) - xp
                     dyy(2) = dy(2) - yp
                     dzz(2) = dz(2) - zp
                     dxx(3) = dx(3) - xp
                     dyy(3) = dy(3) - yp
                     dzz(3) = dz(3) - zp
                     dss(1) = dsqrt(dxx(1)*dxx(1) + dyy(1)*dyy(1) +
     &                              dzz(1)*dzz(1))
                     dss(2) = dsqrt(dxx(2)*dxx(2) + dyy(2)*dyy(2) +
     &                              dzz(2)*dzz(2))
                     dss(3) = dsqrt(dxx(3)*dxx(3) + dyy(3)*dyy(3) +
     &                              dzz(3)*dzz(3))
c
c     calculate the in-plane angles
c
                     doti(1) = dxx(1)*dxx(2)+dyy(1)*dyy(2)+dzz(1)*dzz(2)
                     normi(1) = dss(1) * dss(2)
                     cosi(1) = doti(1) / normi(1)
                     if (abs(cosi(1)) .gt. 1.0d0)
     &                  cosi(1) = cosi(1)/abs(cosi(1))
                     angin(1) = radian * acos(cosi(1))
                     doti(2) = dxx(1)*dxx(3)+dyy(1)*dyy(3)+dzz(1)*dzz(3)
                     normi(2) = dss(1) * dss(3)
                     cosi(2) = doti(2) / normi(2)
                     if (abs(cosi(2)) .gt. 1.0d0)
     &                  cosi(2) = cosi(2)/abs(cosi(2))
                     angin(2) = radian * acos(cosi(2))
                     doti(3) = dxx(2)*dxx(3)+dyy(2)*dyy(3)+dzz(2)*dzz(3)
                     normi(3) = dss(2) * dss(3)
                     cosi(3) = doti(3) / normi(3)
                     if (abs(cosi(3)) .gt. 1.0d0)
     &                  cosi(3) = cosi(3)/abs(cosi(3))
                     angin(3) = radian * acos(cosi(3))
c
c     calculate the out-of-plane angles
c
                     doto(1) = dx(1)*dxx(1)+dy(1)*dyy(1)+dz(1)*dzz(1)
                     normo(1) = ds(1) * dss(1)
                     coso(1) = doto(1) / normo(1)
                     if (abs(coso(1)) .gt. 1.0d0)
     &                  coso(1) = coso(1)/abs(coso(1))
                     angout(1) = radian * acos(coso(1))
                     doto(2) = dx(2)*dxx(2)+dy(2)*dyy(2)+dz(2)*dzz(2)
                     normo(2) = ds(2) * dss(2)
                     coso(2) = doto(2) / normo(2)
                     if (abs(coso(2)) .gt. 1.0d0)
     &                  coso(2) = coso(2)/abs(coso(2))
                     angout(2) = radian * acos(coso(2))
                     doto(3) = dx(3)*dxx(3)+dy(3)*dyy(3)+dz(3)*dzz(3)
                     normo(3) = ds(3) * dss(3)
                     coso(3) = doto(3) / normo(3)
                     if (abs(coso(3)) .gt. 1.0d0)
     &                  coso(3) = coso(3)/abs(coso(3))
                     angout(3) = radian * acos(coso(3))
                  end if
               end if
            end if
c
c     the next loop performs the actual energy calculations;
c     "kang" is the index into the "ksbopb" array
c
            l = 1
            m = 2
            do k = 1, nang
               kang = kang + 1
               ia = i12(l,i)
               ib = i12(m,i)
               nsb = ksbopb(1,kang)
               nopb = ksbopb(2,kang)
c
c     set up some terms used in the calculations below
c
               if (nsb.ne.0 .or. nopb.eq.0) then
                  dt = ang(k) - anat(kang)
                  sine = dsqrt(1.0d0-cos(k)*cos(k))
                  ddtdcos = -1.0d0 / (norm(k)*sine)
                  terma = -cos(k) * ds(m)/ds(l)
                  termb = -cos(k) * ds(l)/ds(m)
                  ddtdxia = ddtdcos * (dx(m) + dx(l)*terma)
                  ddtdyia = ddtdcos * (dy(m) + dy(l)*terma)
                  ddtdzia = ddtdcos * (dz(m) + dz(l)*terma)
                  ddtdxib = ddtdcos * (dx(l) + dx(m)*termb)
                  ddtdyib = ddtdcos * (dy(l) + dy(m)*termb)
                  ddtdzib = ddtdcos * (dz(l) + dz(m)*termb)
c
c     more chain rule terms used in the below
c
                  cotan = cos(k) / sine
c
                  dtxiaxia = -2.0d0*ddtdxia*dx(l)/ds(l)**2
     &      - cotan*(ddtdxia**2+dx(l)*dx(l)/ds(l)**4-1.0d0/ds(l)**2)
                  dtxiayia = -(ddtdxia*dy(l)+ddtdyia*dx(l))/ds(l)**2
     &      - cotan*(ddtdxia*ddtdyia+dx(l)*dy(l)/ds(l)**4)
                  dtxiazia = -(ddtdxia*dz(l)+ddtdzia*dx(l))/ds(l)**2
     &      - cotan*(ddtdxia*ddtdzia+dx(l)*dz(l)/ds(l)**4)
                  dtyiayia = -2.0d0*ddtdyia*dy(l)/ds(l)**2
     &      - cotan*(ddtdyia**2+dy(l)*dy(l)/ds(l)**4-1.0d0/ds(l)**2)
                  dtyiazia = -(ddtdyia*dz(l)+ddtdzia*dy(l))/ds(l)**2
     &      - cotan*(ddtdyia*ddtdzia+dy(l)*dz(l)/ds(l)**4)
                  dtziazia = -2.0d0*ddtdzia*dz(l)/ds(l)**2
     &      - cotan*(ddtdzia**2+dz(l)*dz(l)/ds(l)**4-1.0d0/ds(l)**2)
c
                  dtxibxib = -2.0d0*ddtdxib*dx(m)/ds(m)**2
     &      - cotan*(ddtdxib**2+dx(m)*dx(m)/ds(m)**4-1.0d0/ds(m)**2)
                  dtxibyib = -(ddtdxib*dy(m)+ddtdyib*dx(m))/ds(m)**2
     &      - cotan*(ddtdxib*ddtdyib+dx(m)*dy(m)/ds(m)**4)
                  dtxibzib = -(ddtdxib*dz(m)+ddtdzib*dx(m))/ds(m)**2
     &      - cotan*(ddtdxib*ddtdzib+dx(m)*dz(m)/ds(m)**4)
                  dtyibyib = -2.0d0*ddtdyib*dy(m)/ds(m)**2
     &      - cotan*(ddtdyib**2+dy(m)*dy(m)/ds(m)**4-1.0d0/ds(m)**2)
                  dtyibzib = -(ddtdyib*dz(m)+ddtdzib*dy(m))/ds(m)**2
     &      - cotan*(ddtdyib*ddtdzib+dy(m)*dz(m)/ds(m)**4)
                  dtzibzib = -2.0d0*ddtdzib*dz(m)/ds(m)**2
     &      - cotan*(ddtdzib**2+dz(m)*dz(m)/ds(m)**4-1.0d0/ds(m)**2)
c
                  dtxiaxib = ddtdxib*(-dx(l)/ds(l)**2-cotan*ddtdxia)
     &      + (dx(m)**2/ds(m)**2-1.0d0)/(ds(l)*ds(m)*sine)
                  dtxiayib = ddtdyib*(-dx(l)/ds(l)**2-cotan*ddtdxia)
     &      + (dx(m)*dy(m)/ds(m)**2)/(ds(l)*ds(m)*sine)
                  dtxiazib = ddtdzib*(-dx(l)/ds(l)**2-cotan*ddtdxia)
     &      + (dx(m)*dz(m)/ds(m)**2)/(ds(l)*ds(m)*sine)
                  dtyiaxib = ddtdxib*(-dy(l)/ds(l)**2-cotan*ddtdyia)
     &      + (dy(m)*dx(m)/ds(m)**2)/(ds(l)*ds(m)*sine)
                  dtyiayib = ddtdyib*(-dy(l)/ds(l)**2-cotan*ddtdyia)
     &      + (dy(m)**2/ds(m)**2-1.0d0)/(ds(l)*ds(m)*sine)
                  dtyiazib = ddtdzib*(-dy(l)/ds(l)**2-cotan*ddtdyia)
     &      + (dy(m)*dz(m)/ds(m)**2)/(ds(l)*ds(m)*sine)
                  dtziaxib = ddtdxib*(-dz(l)/ds(l)**2-cotan*ddtdzia)
     &      + (dz(m)*dx(m)/ds(m)**2)/(ds(l)*ds(m)*sine)
                  dtziayib = ddtdyib*(-dz(l)/ds(l)**2-cotan*ddtdzia)
     &      + (dz(m)*dy(m)/ds(m)**2)/(ds(l)*ds(m)*sine)
                  dtziazib = ddtdzib*(-dz(l)/ds(l)**2-cotan*ddtdzia)
     &      + (dz(m)**2/ds(m)**2-1.0d0)/(ds(l)*ds(m)*sine)
c
                  ddtdxia = ddtdxia * radian
                  ddtdyia = ddtdyia * radian
                  ddtdzia = ddtdzia * radian
                  ddtdxib = ddtdxib * radian
                  ddtdyib = ddtdyib * radian
                  ddtdzib = ddtdzib * radian
                  ddtdxi = -ddtdxia - ddtdxib
                  ddtdyi = -ddtdyia - ddtdyib
                  ddtdzi = -ddtdzia - ddtdzib
c
                  dtxiaxia = dtxiaxia * radian
                  dtxiayia = dtxiayia * radian
                  dtxiazia = dtxiazia * radian
                  dtyiayia = dtyiayia * radian
                  dtyiazia = dtyiazia * radian
                  dtziazia = dtziazia * radian
                  dtxibxib = dtxibxib * radian
                  dtxibyib = dtxibyib * radian
                  dtxibzib = dtxibzib * radian
                  dtyibyib = dtyibyib * radian
                  dtyibzib = dtyibzib * radian
                  dtzibzib = dtzibzib * radian
                  dtxiaxib = dtxiaxib * radian
                  dtxiayib = dtxiayib * radian
                  dtxiazib = dtxiazib * radian
                  dtyiaxib = dtyiaxib * radian
                  dtyiayib = dtyiayib * radian
                  dtyiazib = dtyiazib * radian
                  dtziaxib = dtziaxib * radian
                  dtziayib = dtziayib * radian
                  dtziazib = dtziazib * radian
c
                  dtxixia = -dtxiaxia - dtxiaxib
                  dtxiyia = -dtxiayia - dtyiaxib
                  dtxizia = -dtxiazia - dtziaxib
                  dtyixia = -dtxiayia - dtxiayib
                  dtyiyia = -dtyiayia - dtyiayib
                  dtyizia = -dtyiazia - dtziayib
                  dtzixia = -dtxiazia - dtxiazib
                  dtziyia = -dtyiazia - dtyiazib
                  dtzizia = -dtziazia - dtziazib
c
                  dtxixib = -dtxibxib - dtxiaxib
                  dtxiyib = -dtxibyib - dtxiayib
                  dtxizib = -dtxibzib - dtxiazib
                  dtyixib = -dtxibyib - dtyiaxib
                  dtyiyib = -dtyibyib - dtyiayib
                  dtyizib = -dtyibzib - dtyiazib
                  dtzixib = -dtxibzib - dtziaxib
                  dtziyib = -dtyibzib - dtziayib
                  dtzizib = -dtzibzib - dtziazib
c
                  dtxixi = -dtxixia + dtxibxib + dtxiaxib
                  dtxiyi = -dtxiyia + dtxibyib + dtyiaxib
                  dtxizi = -dtxizia + dtxibzib + dtziaxib
                  dtyiyi = -dtyiyia + dtyibyib + dtyiayib
                  dtyizi = -dtyizia + dtyibzib + dtziayib
                  dtzizi = -dtzizia + dtzibzib + dtziazib
               end if
c
c     stretch bend calculation
c
               if (nsb .ne. 0) then
                  ii = ksbopb(3,kang)
                  kk = ksbopb(4,kang)
                  term = 2.51118d0 * sbk(nsb)
                  if (ii .eq. 0) then
                     dr = 0.0d0
                     ddrdxia = 0.0d0
                     ddrdyia = 0.0d0
                     ddrdzia = 0.0d0
                     drxiaxia = 0.0d0
                     drxiayia = 0.0d0
                     drxiazia = 0.0d0
                     dryiayia = 0.0d0
                     dryiazia = 0.0d0
                     drziazia = 0.0d0
                  else
                     dr = ds(l) - bl(ii)
                     ddrdxia = dx(l) / ds(l)
                     ddrdyia = dy(l) / ds(l)
                     ddrdzia = dz(l) / ds(l)
                     drxiaxia = (1.0d0-ddrdxia*ddrdxia) / ds(l)
                     drxiayia = -ddrdxia*ddrdyia / ds(l)
                     drxiazia = -ddrdxia*ddrdzia / ds(l)
                     dryiayia = (1.0d0-ddrdyia*ddrdyia) / ds(l)
                     dryiazia = -ddrdyia*ddrdzia / ds(l)
                     drziazia = (1.0d0-ddrdzia*ddrdzia) / ds(l)
                  end if
                  if (kk .eq. 0) then
                     ddrdxib = 0.0d0
                     ddrdyib = 0.0d0
                     ddrdzib = 0.0d0
                     drxibxib = 0.0d0
                     drxibyib = 0.0d0
                     drxibzib = 0.0d0
                     dryibyib = 0.0d0
                     dryibzib = 0.0d0
                     drzibzib = 0.0d0
                  else
                     dr = dr + ds(m) - bl(kk)
                     ddrdxib = dx(m) / ds(m)
                     ddrdyib = dy(m) / ds(m)
                     ddrdzib = dz(m) / ds(m)
                     drxibxib = (1.0d0-ddrdxib*ddrdxib) / ds(m)
                     drxibyib = -ddrdxib*ddrdyib / ds(m)
                     drxibzib = -ddrdxib*ddrdzib / ds(m)
                     dryibyib = (1.0d0-ddrdyib*ddrdyib) / ds(m)
                     dryibzib = -ddrdyib*ddrdzib / ds(m)
                     drzibzib = (1.0d0-ddrdzib*ddrdzib) / ds(m)
                  end if
c
                  dr = dr * term
                  ddrdxia = ddrdxia * term
                  ddrdyia = ddrdyia * term
                  ddrdzia = ddrdzia * term
                  ddrdxib = ddrdxib * term
                  ddrdyib = ddrdyib * term
                  ddrdzib = ddrdzib * term
                  ddrdxi = -ddrdxia - ddrdxib
                  ddrdyi = -ddrdyia - ddrdyib
                  ddrdzi = -ddrdzia - ddrdzib
c
                  drxiaxia = drxiaxia * term
                  drxiayia = drxiayia * term
                  drxiazia = drxiazia * term
                  dryiayia = dryiayia * term
                  dryiazia = dryiazia * term
                  drziazia = drziazia * term
                  drxibxib = drxibxib * term
                  drxibyib = drxibyib * term
                  drxibzib = drxibzib * term
                  dryibyib = dryibyib * term
                  dryibzib = dryibzib * term
                  drzibzib = drzibzib * term
c
                  drxixi = drxiaxia + drxibxib
                  drxiyi = drxiayia + drxibyib
                  drxizi = drxiazia + drxibzib
                  dryiyi = dryiayia + dryibyib
                  dryizi = dryiazia + dryibzib
                  drzizi = drziazia + drzibzib
c
c     stretch bend -- diagonal hessian block elements
c
                  heba(1,i,1,i) = heba(1,i,1,i) + dt*drxixi
     &               + 2.0d0*ddtdxi*ddrdxi + dr*dtxixi
                  heba(1,i,2,i) = heba(1,i,2,i) + dt*drxiyi
     &               + ddtdxi*ddrdyi + ddtdyi*ddrdxi + dr*dtxiyi
                  heba(1,i,3,i) = heba(1,i,3,i) + dt*drxizi
     &               + ddtdxi*ddrdzi + ddtdzi*ddrdxi + dr*dtxizi
                  heba(2,i,2,i) = heba(2,i,2,i) + dt*dryiyi
     &               + 2.0d0*ddtdyi*ddrdyi + dr*dtyiyi
                  heba(2,i,3,i) = heba(2,i,3,i) + dt*dryizi
     &               + ddtdyi*ddrdzi + ddtdzi*ddrdyi + dr*dtyizi
                  heba(3,i,3,i) = heba(3,i,3,i) + dt*drzizi
     &               + 2.0d0*ddtdzi*ddrdzi + dr*dtzizi
c
                  heba(1,ia,1,ia) = heba(1,ia,1,ia) + dt*drxiaxia
     &               + 2.0d0*ddtdxia*ddrdxia + dr*dtxiaxia
                  heba(1,ia,2,ia) = heba(1,ia,2,ia) + dt*drxiayia
     &               + ddtdxia*ddrdyia + ddtdyia*ddrdxia + dr*dtxiayia
                  heba(1,ia,3,ia) = heba(1,ia,3,ia) + dt*drxiazia
     &               + ddtdxia*ddrdzia + ddtdzia*ddrdxia + dr*dtxiazia
                  heba(2,ia,2,ia) = heba(2,ia,2,ia) + dt*dryiayia
     &               + 2.0d0*ddtdyia*ddrdyia + dr*dtyiayia
                  heba(2,ia,3,ia) = heba(2,ia,3,ia) + dt*dryiazia
     &               + ddtdyia*ddrdzia + ddtdzia*ddrdyia + dr*dtyiazia
                  heba(3,ia,3,ia) = heba(3,ia,3,ia) + dt*drziazia
     &               + 2.0d0*ddtdzia*ddrdzia + dr*dtziazia
c
                  heba(1,ib,1,ib) = heba(1,ib,1,ib) + dt*drxibxib
     &               + 2.0d0*ddtdxib*ddrdxib + dr*dtxibxib
                  heba(1,ib,2,ib) = heba(1,ib,2,ib) + dt*drxibyib
     &               + ddtdxib*ddrdyib + ddtdyib*ddrdxib + dr*dtxibyib
                  heba(1,ib,3,ib) = heba(1,ib,3,ib) + dt*drxibzib
     &               + ddtdxib*ddrdzib + ddtdzib*ddrdxib + dr*dtxibzib
                  heba(2,ib,2,ib) = heba(2,ib,2,ib) + dt*dryibyib
     &               + 2.0d0*ddtdyib*ddrdyib + dr*dtyibyib
                  heba(2,ib,3,ib) = heba(2,ib,3,ib) + dt*dryibzib
     &               + ddtdyib*ddrdzib + ddtdzib*ddrdyib + dr*dtyibzib
                  heba(3,ib,3,ib) = heba(3,ib,3,ib) + dt*drzibzib
     &               + 2.0d0*ddtdzib*ddrdzib + dr*dtzibzib
c
c     stretch bend -- above-diagonal hessian block elements
c
                  if (i .lt. ia) then
                     heba(1,i,1,ia) = heba(1,i,1,ia) - dt*drxiaxia
     &                  + ddtdxi*ddrdxia + ddtdxia*ddrdxi + dr*dtxixia
                     heba(1,i,2,ia) = heba(1,i,2,ia) - dt*drxiayia
     &                  + ddtdxi*ddrdyia + ddtdyia*ddrdxi + dr*dtxiyia
                     heba(1,i,3,ia) = heba(1,i,3,ia) - dt*drxiazia
     &                  + ddtdxi*ddrdzia + ddtdzia*ddrdxi + dr*dtxizia
                     heba(2,i,1,ia) = heba(2,i,1,ia) - dt*drxiayia
     &                  + ddtdyi*ddrdxia + ddtdxia*ddrdyi + dr*dtyixia
                     heba(2,i,2,ia) = heba(2,i,2,ia) - dt*dryiayia
     &                  + ddtdyi*ddrdyia + ddtdyia*ddrdyi + dr*dtyiyia
                     heba(2,i,3,ia) = heba(2,i,3,ia) - dt*dryiazia
     &                  + ddtdyi*ddrdzia + ddtdzia*ddrdyi + dr*dtyizia
                     heba(3,i,1,ia) = heba(3,i,1,ia) - dt*drxiazia
     &                  + ddtdzi*ddrdxia + ddtdxia*ddrdzi + dr*dtzixia
                     heba(3,i,2,ia) = heba(3,i,2,ia) - dt*dryiazia
     &                  + ddtdzi*ddrdyia + ddtdyia*ddrdzi + dr*dtziyia
                     heba(3,i,3,ia) = heba(3,i,3,ia) - dt*drziazia
     &                  + ddtdzi*ddrdzia + ddtdzia*ddrdzi + dr*dtzizia
                  else
                     heba(1,ia,1,i) = heba(1,ia,1,i) - dt*drxiaxia
     &                  + ddtdxia*ddrdxi + ddtdxi*ddrdxia + dr*dtxixia
                     heba(1,ia,2,i) = heba(1,ia,2,i) - dt*drxiayia
     &                  + ddtdxia*ddrdyi + ddtdyi*ddrdxia + dr*dtxiyia
                     heba(1,ia,3,i) = heba(1,ia,3,i) - dt*drxiazia
     &                  + ddtdxia*ddrdzi + ddtdzi*ddrdxia + dr*dtxizia
                     heba(2,ia,1,i) = heba(2,ia,1,i) - dt*drxiayia
     &                  + ddtdyia*ddrdxi + ddtdxi*ddrdyia + dr*dtyixia
                     heba(2,ia,2,i) = heba(2,ia,2,i) - dt*dryiayia
     &                  + ddtdyia*ddrdyi + ddtdyi*ddrdyia + dr*dtyiyia
                     heba(2,ia,3,i) = heba(2,ia,3,i) - dt*dryiazia
     &                  + ddtdyia*ddrdzi + ddtdzi*ddrdyia + dr*dtyizia
                     heba(3,ia,1,i) = heba(3,ia,1,i) - dt*drxiazia
     &                  + ddtdzia*ddrdxi + ddtdxi*ddrdzia + dr*dtzixia
                     heba(3,ia,2,i) = heba(3,ia,2,i) - dt*dryiazia
     &                  + ddtdzia*ddrdyi + ddtdyi*ddrdzia + dr*dtziyia
                     heba(3,ia,3,i) = heba(3,ia,3,i) - dt*drziazia
     &                  + ddtdzia*ddrdzi + ddtdzi*ddrdzia + dr*dtzizia
                  end if
c
                  if (i .lt. ib) then
                     heba(1,i,1,ib) = heba(1,i,1,ib) - dt*drxibxib
     &                  + ddtdxi*ddrdxib + ddtdxib*ddrdxi + dr*dtxixib
                     heba(1,i,2,ib) = heba(1,i,2,ib) - dt*drxibyib
     &                  + ddtdxi*ddrdyib + ddtdyib*ddrdxi + dr*dtxiyib
                     heba(1,i,3,ib) = heba(1,i,3,ib) - dt*drxibzib
     &                  + ddtdxi*ddrdzib + ddtdzib*ddrdxi + dr*dtxizib
                     heba(2,i,1,ib) = heba(2,i,1,ib) - dt*drxibyib
     &                  + ddtdyi*ddrdxib + ddtdxib*ddrdyi + dr*dtyixib
                     heba(2,i,2,ib) = heba(2,i,2,ib) - dt*dryibyib
     &                  + ddtdyi*ddrdyib + ddtdyib*ddrdyi + dr*dtyiyib
                     heba(2,i,3,ib) = heba(2,i,3,ib) - dt*dryibzib
     &                  + ddtdyi*ddrdzib + ddtdzib*ddrdyi + dr*dtyizib
                     heba(3,i,1,ib) = heba(3,i,1,ib) - dt*drxibzib
     &                  + ddtdzi*ddrdxib + ddtdxib*ddrdzi + dr*dtzixib
                     heba(3,i,2,ib) = heba(3,i,2,ib) - dt*dryibzib
     &                  + ddtdzi*ddrdyib + ddtdyib*ddrdzi + dr*dtziyib
                     heba(3,i,3,ib) = heba(3,i,3,ib) - dt*drzibzib
     &                  + ddtdzi*ddrdzib + ddtdzib*ddrdzi + dr*dtzizib
                  else
                     heba(1,ib,1,i) = heba(1,ib,1,i) - dt*drxibxib
     &                  + ddtdxib*ddrdxi + ddtdxi*ddrdxib + dr*dtxixib
                     heba(1,ib,2,i) = heba(1,ib,2,i) - dt*drxibyib
     &                  + ddtdxib*ddrdyi + ddtdyi*ddrdxib + dr*dtxiyib
                     heba(1,ib,3,i) = heba(1,ib,3,i) - dt*drxibzib
     &                  + ddtdxib*ddrdzi + ddtdzi*ddrdxib + dr*dtxizib
                     heba(2,ib,1,i) = heba(2,ib,1,i) - dt*drxibyib
     &                  + ddtdyib*ddrdxi + ddtdxi*ddrdyib + dr*dtyixib
                     heba(2,ib,2,i) = heba(2,ib,2,i) - dt*dryibyib
     &                  + ddtdyib*ddrdyi + ddtdyi*ddrdyib + dr*dtyiyib
                     heba(2,ib,3,i) = heba(2,ib,3,i) - dt*dryibzib
     &                  + ddtdyib*ddrdzi + ddtdzi*ddrdyib + dr*dtyizib
                     heba(3,ib,1,i) = heba(3,ib,1,i) - dt*drxibzib
     &                  + ddtdzib*ddrdxi + ddtdxi*ddrdzib + dr*dtzixib
                     heba(3,ib,2,i) = heba(3,ib,2,i) - dt*dryibzib
     &                  + ddtdzib*ddrdyi + ddtdyi*ddrdzib + dr*dtziyib
                     heba(3,ib,3,i) = heba(3,ib,3,i) - dt*drzibzib
     &                  + ddtdzib*ddrdzi + ddtdzi*ddrdzib + dr*dtzizib
                  end if
c
                  if (ia .lt. ib) then
                     heba(1,ia,1,ib) = heba(1,ia,1,ib) + ddtdxia*ddrdxib
     &                  + ddtdxib*ddrdxia + dr*dtxiaxib
                     heba(1,ia,2,ib) = heba(1,ia,2,ib) + ddtdxia*ddrdyib
     &                  + ddtdyib*ddrdxia + dr*dtxiayib
                     heba(1,ia,3,ib) = heba(1,ia,3,ib) + ddtdxia*ddrdzib
     &                  + ddtdzib*ddrdxia + dr*dtxiazib
                     heba(2,ia,1,ib) = heba(2,ia,1,ib) + ddtdyia*ddrdxib
     &                  + ddtdxib*ddrdyia + dr*dtyiaxib
                     heba(2,ia,2,ib) = heba(2,ia,2,ib) + ddtdyia*ddrdyib
     &                  + ddtdyib*ddrdyia + dr*dtyiayib
                     heba(2,ia,3,ib) = heba(2,ia,3,ib) + ddtdyia*ddrdzib
     &                  + ddtdzib*ddrdyia + dr*dtyiazib
                     heba(3,ia,1,ib) = heba(3,ia,1,ib) + ddtdzia*ddrdxib
     &                  + ddtdxib*ddrdzia + dr*dtziaxib
                     heba(3,ia,2,ib) = heba(3,ia,2,ib) + ddtdzia*ddrdyib
     &                  + ddtdyib*ddrdzia + dr*dtziayib
                     heba(3,ia,3,ib) = heba(3,ia,3,ib) + ddtdzia*ddrdzib
     &                  + ddtdzib*ddrdzia + dr*dtziazib
                  else
                     heba(1,ib,1,ia) = heba(1,ib,1,ia) + ddtdxib*ddrdxia
     &                  + ddtdxia*ddrdxib + dr*dtxiaxib
                     heba(1,ib,2,ia) = heba(1,ib,2,ia) + ddtdxib*ddrdyia
     &                  + ddtdyia*ddrdxib + dr*dtyiaxib
                     heba(1,ib,3,ia) = heba(1,ib,3,ia) + ddtdxib*ddrdzia
     &                  + ddtdzia*ddrdxib + dr*dtziaxib
                     heba(2,ib,1,ia) = heba(2,ib,1,ia) + ddtdyib*ddrdxia
     &                  + ddtdxia*ddrdyib + dr*dtxiayib
                     heba(2,ib,2,ia) = heba(2,ib,2,ia) + ddtdyib*ddrdyia
     &                  + ddtdyia*ddrdyib + dr*dtyiayib
                     heba(2,ib,3,ia) = heba(2,ib,3,ia) + ddtdyib*ddrdzia
     &                  + ddtdzia*ddrdyib + dr*dtziayib
                     heba(3,ib,1,ia) = heba(3,ib,1,ia) + ddtdzib*ddrdxia
     &                  + ddtdxia*ddrdzib + dr*dtxiazib
                     heba(3,ib,2,ia) = heba(3,ib,2,ia) + ddtdzib*ddrdyia
     &                  + ddtdyia*ddrdzib + dr*dtyiazib
                     heba(3,ib,3,ia) = heba(3,ib,3,ia) + ddtdzib*ddrdzia
     &                  + ddtdzia*ddrdzib + dr*dtziazib
                  end if
               end if
c
c     standard angle bend calculation
c
               if (nopb .eq. 0) then
                  dt2 = dt * dt
                  deddt = 0.02191418d0 * acon(kang) * dt
     &                       * (2.0d0 + 0.00006d0*sf*dt2*dt2)
                  d2eddt2 = 0.02191418d0 * acon(kang)
     &                         * (2.0d0 + 0.0003d0*sf*dt2*dt2)
c
c     standard angle bend -- diagonal hessian block elements
c
                  hea(1,i,1,i)   = hea(1,i,1,i) + deddt*dtxixi
     &                                 + d2eddt2*ddtdxi*ddtdxi
                  hea(1,i,2,i)   = hea(1,i,2,i) + deddt*dtxiyi
     &                                 + d2eddt2*ddtdxi*ddtdyi
                  hea(1,i,3,i)   = hea(1,i,3,i) + deddt*dtxizi
     &                                 + d2eddt2*ddtdxi*ddtdzi
                  hea(2,i,2,i)   = hea(2,i,2,i) + deddt*dtyiyi
     &                                 + d2eddt2*ddtdyi*ddtdyi
                  hea(2,i,3,i)   = hea(2,i,3,i) + deddt*dtyizi
     &                                 + d2eddt2*ddtdyi*ddtdzi
                  hea(3,i,3,i)   = hea(3,i,3,i) + deddt*dtzizi
     &                                 + d2eddt2*ddtdzi*ddtdzi
c
                  hea(1,ia,1,ia) = hea(1,ia,1,ia) + deddt*dtxiaxia
     &                                 + d2eddt2*ddtdxia*ddtdxia
                  hea(1,ia,2,ia) = hea(1,ia,2,ia) + deddt*dtxiayia
     &                                 + d2eddt2*ddtdxia*ddtdyia
                  hea(1,ia,3,ia) = hea(1,ia,3,ia) + deddt*dtxiazia
     &                                 + d2eddt2*ddtdxia*ddtdzia
                  hea(2,ia,2,ia) = hea(2,ia,2,ia) + deddt*dtyiayia
     &                                 + d2eddt2*ddtdyia*ddtdyia
                  hea(2,ia,3,ia) = hea(2,ia,3,ia) + deddt*dtyiazia
     &                                 + d2eddt2*ddtdyia*ddtdzia
                  hea(3,ia,3,ia) = hea(3,ia,3,ia) + deddt*dtziazia
     &                                 + d2eddt2*ddtdzia*ddtdzia
c
                  hea(1,ib,1,ib) = hea(1,ib,1,ib) + deddt*dtxibxib
     &                                 + d2eddt2*ddtdxib*ddtdxib
                  hea(1,ib,2,ib) = hea(1,ib,2,ib) + deddt*dtxibyib
     &                                 + d2eddt2*ddtdxib*ddtdyib
                  hea(1,ib,3,ib) = hea(1,ib,3,ib) + deddt*dtxibzib
     &                                 + d2eddt2*ddtdxib*ddtdzib
                  hea(2,ib,2,ib) = hea(2,ib,2,ib) + deddt*dtyibyib
     &                                 + d2eddt2*ddtdyib*ddtdyib
                  hea(2,ib,3,ib) = hea(2,ib,3,ib) + deddt*dtyibzib
     &                                 + d2eddt2*ddtdyib*ddtdzib
                  hea(3,ib,3,ib) = hea(3,ib,3,ib) + deddt*dtzibzib
     &                                 + d2eddt2*ddtdzib*ddtdzib
c
c     standard angle bend -- above-diagonal hessian block elements
c
                  if (i .lt. ia) then
                     hea(1,i,1,ia) = hea(1,i,1,ia) + deddt*dtxixia
     &                                    + d2eddt2*ddtdxi*ddtdxia
                     hea(1,i,2,ia) = hea(1,i,2,ia) + deddt*dtxiyia
     &                                    + d2eddt2*ddtdxi*ddtdyia
                     hea(1,i,3,ia) = hea(1,i,3,ia) + deddt*dtxizia
     &                                    + d2eddt2*ddtdxi*ddtdzia
                     hea(2,i,1,ia) = hea(2,i,1,ia) + deddt*dtyixia
     &                                    + d2eddt2*ddtdyi*ddtdxia
                     hea(2,i,2,ia) = hea(2,i,2,ia) + deddt*dtyiyia
     &                                    + d2eddt2*ddtdyi*ddtdyia
                     hea(2,i,3,ia) = hea(2,i,3,ia) + deddt*dtyizia
     &                                    + d2eddt2*ddtdyi*ddtdzia
                     hea(3,i,1,ia) = hea(3,i,1,ia) + deddt*dtzixia
     &                                    + d2eddt2*ddtdzi*ddtdxia
                     hea(3,i,2,ia) = hea(3,i,2,ia) + deddt*dtziyia
     &                                    + d2eddt2*ddtdzi*ddtdyia
                     hea(3,i,3,ia) = hea(3,i,3,ia) + deddt*dtzizia
     &                                    + d2eddt2*ddtdzi*ddtdzia
                  else
                     hea(1,ia,1,i) = hea(1,ia,1,i) + deddt*dtxixia
     &                                    + d2eddt2*ddtdxia*ddtdxi
                     hea(1,ia,2,i) = hea(1,ia,2,i) + deddt*dtyixia
     &                                    + d2eddt2*ddtdxia*ddtdyi
                     hea(1,ia,3,i) = hea(1,ia,3,i) + deddt*dtzixia
     &                                    + d2eddt2*ddtdxia*ddtdzi
                     hea(2,ia,1,i) = hea(2,ia,1,i) + deddt*dtxiyia
     &                                    + d2eddt2*ddtdyia*ddtdxi
                     hea(2,ia,2,i) = hea(2,ia,2,i) + deddt*dtyiyia
     &                                    + d2eddt2*ddtdyia*ddtdyi
                     hea(2,ia,3,i) = hea(2,ia,3,i) + deddt*dtziyia
     &                                    + d2eddt2*ddtdyia*ddtdzi
                     hea(3,ia,1,i) = hea(3,ia,1,i) + deddt*dtxizia
     &                                    + d2eddt2*ddtdzia*ddtdxi
                     hea(3,ia,2,i) = hea(3,ia,2,i) + deddt*dtyizia
     &                                    + d2eddt2*ddtdzia*ddtdyi
                     hea(3,ia,3,i) = hea(3,ia,3,i) + deddt*dtzizia
     &                                    + d2eddt2*ddtdzia*ddtdzi
                  end if
c
                  if (i .lt. ib) then
                     hea(1,i,1,ib) = hea(1,i,1,ib) + deddt*dtxixib
     &                                    + d2eddt2*ddtdxi*ddtdxib
                     hea(1,i,2,ib) = hea(1,i,2,ib) + deddt*dtxiyib
     &                                    + d2eddt2*ddtdxi*ddtdyib
                     hea(1,i,3,ib) = hea(1,i,3,ib) + deddt*dtxizib
     &                                    + d2eddt2*ddtdxi*ddtdzib
                     hea(2,i,1,ib) = hea(2,i,1,ib) + deddt*dtyixib
     &                                    + d2eddt2*ddtdyi*ddtdxib
                     hea(2,i,2,ib) = hea(2,i,2,ib) + deddt*dtyiyib
     &                                    + d2eddt2*ddtdyi*ddtdyib
                     hea(2,i,3,ib) = hea(2,i,3,ib) + deddt*dtyizib
     &                                    + d2eddt2*ddtdyi*ddtdzib
                     hea(3,i,1,ib) = hea(3,i,1,ib) + deddt*dtzixib
     &                                    + d2eddt2*ddtdzi*ddtdxib
                     hea(3,i,2,ib) = hea(3,i,2,ib) + deddt*dtziyib
     &                                    + d2eddt2*ddtdzi*ddtdyib
                     hea(3,i,3,ib) = hea(3,i,3,ib) + deddt*dtzizib
     &                                    + d2eddt2*ddtdzi*ddtdzib
                  else
                     hea(1,ib,1,i) = hea(1,ib,1,i) + deddt*dtxixib
     &                                    + d2eddt2*ddtdxib*ddtdxi
                     hea(1,ib,2,i) = hea(1,ib,2,i) + deddt*dtyixib
     &                                    + d2eddt2*ddtdxib*ddtdyi
                     hea(1,ib,3,i) = hea(1,ib,3,i) + deddt*dtzixib
     &                                    + d2eddt2*ddtdxib*ddtdzi
                     hea(2,ib,1,i) = hea(2,ib,1,i) + deddt*dtxiyib
     &                                    + d2eddt2*ddtdyib*ddtdxi
                     hea(2,ib,2,i) = hea(2,ib,2,i) + deddt*dtyiyib
     &                                    + d2eddt2*ddtdyib*ddtdyi
                     hea(2,ib,3,i) = hea(2,ib,3,i) + deddt*dtziyib
     &                                    + d2eddt2*ddtdyib*ddtdzi
                     hea(3,ib,1,i) = hea(3,ib,1,i) + deddt*dtxizib
     &                                    + d2eddt2*ddtdzib*ddtdxi
                     hea(3,ib,2,i) = hea(3,ib,2,i) + deddt*dtyizib
     &                                    + d2eddt2*ddtdzib*ddtdyi
                     hea(3,ib,3,i) = hea(3,ib,3,i) + deddt*dtzizib
     &                                    + d2eddt2*ddtdzib*ddtdzi
                  end if
c
                  if (ia .lt. ib) then
                     hea(1,ia,1,ib) = hea(1,ia,1,ib) + deddt*dtxiaxib
     &                                    + d2eddt2*ddtdxia*ddtdxib
                     hea(1,ia,2,ib) = hea(1,ia,2,ib) + deddt*dtxiayib
     &                                    + d2eddt2*ddtdxia*ddtdyib
                     hea(1,ia,3,ib) = hea(1,ia,3,ib) + deddt*dtxiazib
     &                                    + d2eddt2*ddtdxia*ddtdzib
                     hea(2,ia,1,ib) = hea(2,ia,1,ib) + deddt*dtyiaxib
     &                                    + d2eddt2*ddtdyia*ddtdxib
                     hea(2,ia,2,ib) = hea(2,ia,2,ib) + deddt*dtyiayib
     &                                    + d2eddt2*ddtdyia*ddtdyib
                     hea(2,ia,3,ib) = hea(2,ia,3,ib) + deddt*dtyiazib
     &                                    + d2eddt2*ddtdyia*ddtdzib
                     hea(3,ia,1,ib) = hea(3,ia,1,ib) + deddt*dtziaxib
     &                                    + d2eddt2*ddtdzia*ddtdxib
                     hea(3,ia,2,ib) = hea(3,ia,2,ib) + deddt*dtziayib
     &                                    + d2eddt2*ddtdzia*ddtdyib
                     hea(3,ia,3,ib) = hea(3,ia,3,ib) + deddt*dtziazib
     &                                    + d2eddt2*ddtdzia*ddtdzib
                  else
                     hea(1,ib,1,ia) = hea(1,ib,1,ia) + deddt*dtxiaxib
     &                                    + d2eddt2*ddtdxib*ddtdxia
                     hea(1,ib,2,ia) = hea(1,ib,2,ia) + deddt*dtyiaxib
     &                                    + d2eddt2*ddtdxib*ddtdyia
                     hea(1,ib,3,ia) = hea(1,ib,3,ia) + deddt*dtziaxib
     &                                    + d2eddt2*ddtdxib*ddtdzia
                     hea(2,ib,1,ia) = hea(2,ib,1,ia) + deddt*dtxiayib
     &                                    + d2eddt2*ddtdyib*ddtdxia
                     hea(2,ib,2,ia) = hea(2,ib,2,ia) + deddt*dtyiayib
     &                                    + d2eddt2*ddtdyib*ddtdyia
                     hea(2,ib,3,ia) = hea(2,ib,3,ia) + deddt*dtziayib
     &                                    + d2eddt2*ddtdyib*ddtdzia
                     hea(3,ib,1,ia) = hea(3,ib,1,ia) + deddt*dtxiazib
     &                                    + d2eddt2*ddtdzib*ddtdxia
                     hea(3,ib,2,ia) = hea(3,ib,2,ia) + deddt*dtyiazib
     &                                    + d2eddt2*ddtdzib*ddtdyia
                     hea(3,ib,3,ia) = hea(3,ib,3,ia) + deddt*dtziazib
     &                                    + d2eddt2*ddtdzib*ddtdzia
                  end if
c
c     out-of-plane bend calculation; out-of-plane portion
c
               else
                  ia = i12(1,i)
                  ib = i12(2,i)
                  ic = i12(3,i)
                  dq = angout(k)
                  dq2 = dq * dq
                  deoddq = 0.02191418d0 * opbk(nopb) * dq
     &                        * (2.0d0 + 0.00006d0*sf*dq2*dq2)
                  d2eoddq2 = 0.02191418d0 * opbk(nopb)
     &                          * (2.0d0 + 0.0003d0*sf*dq2*dq2)
                  if (dq .eq. 0.0d0) then
                     sine = 1.0d0   !! deoddq is zero anyway
                  else if (dq .lt. 0.1d0) then
                     dq = dq / radian
                     sine = dq - dq**3/6.0d0 + dq**5/120.0d0
                  else
                     sine = dsqrt(1.0d0-coso(k)*coso(k))
                  end if
                  ddqdcoso = -radian / (normo(k)*sine)
                  d2dqcoso2 = ddqdcoso * coso(k)/(normo(k)*sine**2)
                  deodcoso = deoddq * ddqdcoso
                  d2eodcoso2 = d2eoddq2*ddqdcoso + deoddq*d2dqcoso2
                  if (k .eq. 1) then
                     dxxss(1) = dxx(1) / dss(1)
                     dxxss(2) = dxx(2) / dss(2)
                     dxxss(3) = dxx(3) / dss(3)
                     dyyss(1) = dyy(1) / dss(1)
                     dyyss(2) = dyy(2) / dss(2)
                     dyyss(3) = dyy(3) / dss(3)
                     dzzss(1) = dzz(1) / dss(1)
                     dzzss(2) = dzz(2) / dss(2)
                     dzzss(3) = dzz(3) / dss(3)
c
                     ddxxdxi = a*a/s - 1.0d0
                     ddyydxi = b*a/s
                     ddzzdxi = c*a/s
                     ddsdxi(1) = -dx(1) / ds(1)
                     ddsdxi(2) = -dx(2) / ds(2)
                     ddsdxi(3) = -dx(3) / ds(3)
                     ddssdxi(1) = ddxxdxi*dxxss(1) +
     &                    ddyydxi*dyyss(1) + ddzzdxi*dzzss(1)
                     ddssdxi(2) = ddxxdxi*dxxss(2) +
     &                    ddyydxi*dyyss(2) + ddzzdxi*dzzss(2)
                     ddssdxi(3) = ddxxdxi*dxxss(3) +
     &                    ddyydxi*dyyss(3) + ddzzdxi*dzzss(3)
                     ddxxdyi = a*b/s
                     ddyydyi = b*b/s - 1.0d0
                     ddzzdyi = c*b/s
                     ddsdyi(1) = -dy(1) / ds(1)
                     ddsdyi(2) = -dy(2) / ds(2)
                     ddsdyi(3) = -dy(3) / ds(3)
                     ddssdyi(1) = ddxxdyi*dxxss(1) +
     &                    ddyydyi*dyyss(1) + ddzzdyi*dzzss(1)
                     ddssdyi(2) = ddxxdyi*dxxss(2) +
     &                    ddyydyi*dyyss(2) + ddzzdyi*dzzss(2)
                     ddssdyi(3) = ddxxdyi*dxxss(3) +
     &                    ddyydyi*dyyss(3) + ddzzdyi*dzzss(3)
                     ddxxdzi = a*c/s
                     ddyydzi = b*c/s
                     ddzzdzi = c*c/s - 1.0d0
                     ddsdzi(1) = -dz(1) / ds(1)
                     ddsdzi(2) = -dz(2) / ds(2)
                     ddsdzi(3) = -dz(3) / ds(3)
                     ddssdzi(1) = ddxxdzi*dxxss(1) +
     &                    ddyydzi*dyyss(1) + ddzzdzi*dzzss(1)
                     ddssdzi(2) = ddxxdzi*dxxss(2) +
     &                    ddyydzi*dyyss(2) + ddzzdzi*dzzss(2)
                     ddssdzi(3) = ddxxdzi*dxxss(3) +
     &                    ddyydzi*dyyss(3) + ddzzdzi*dzzss(3)
c
                     dbdxia = dz(3) - dz(2)
                     dcdxia = dy(2) - dy(3)
                     dddxia = a + dy(1)*dbdxia + dz(1)*dcdxia
                     dsdxia = 2.0d0 * (b*dbdxia + c*dcdxia)
                     ddxdxia(1) = 1.0d0
                     ddxdxia(2) = 0.0d0
                     ddxdxia(3) = 0.0d0
                     ddxxdxia(1) = (xp*dsdxia-a*dddxia)/s + 1.0d0
                     ddxxdxia(2) = (xp*dsdxia-a*dddxia)/s
                     ddxxdxia(3) = (xp*dsdxia-a*dddxia)/s
                     ddyydxia = (yp*dsdxia-b*dddxia-d*dbdxia)/s
                     ddzzdxia = (zp*dsdxia-c*dddxia-d*dcdxia)/s
                     ddsdxia(1) = dx(1) / ds(1)
                     ddsdxia(2) = 0.0d0
                     ddsdxia(3) = 0.0d0
                     ddssdxia(1) = ddxxdxia(1)*dxxss(1) +
     &                    ddyydxia*dyyss(1) + ddzzdxia*dzzss(1)
                     ddssdxia(2) = ddxxdxia(2)*dxxss(2) +
     &                    ddyydxia*dyyss(2) + ddzzdxia*dzzss(2)
                     ddssdxia(3) = ddxxdxia(3)*dxxss(3) +
     &                    ddyydxia*dyyss(3) + ddzzdxia*dzzss(3)
                     dadyia = dz(2) - dz(3)
                     dcdyia = dx(3) - dx(2)
                     dddyia = dx(1)*dadyia + b + dz(1)*dcdyia
                     dsdyia = 2.0d0 * (a*dadyia + c*dcdyia)
                     ddydyia(1) = 1.0d0
                     ddydyia(2) = 0.0d0
                     ddydyia(3) = 0.0d0
                     ddxxdyia = (xp*dsdyia-a*dddyia-d*dadyia)/s
                     ddyydyia(1) = (yp*dsdyia-b*dddyia)/s + 1.0d0
                     ddyydyia(2) = (yp*dsdyia-b*dddyia)/s
                     ddyydyia(3) = (yp*dsdyia-b*dddyia)/s
                     ddzzdyia = (zp*dsdyia-c*dddyia-d*dcdyia)/s
                     ddsdyia(1) = dy(1) / ds(1)
                     ddsdyia(2) = 0.0d0
                     ddsdyia(3) = 0.0d0
                     ddssdyia(1) = ddxxdyia*dxxss(1) +
     &                    ddyydyia(1)*dyyss(1) + ddzzdyia*dzzss(1)
                     ddssdyia(2) = ddxxdyia*dxxss(2) +
     &                    ddyydyia(2)*dyyss(2) + ddzzdyia*dzzss(2)
                     ddssdyia(3) = ddxxdyia*dxxss(3) +
     &                    ddyydyia(3)*dyyss(3) + ddzzdyia*dzzss(3)
                     dadzia = dy(3) - dy(2)
                     dbdzia = dx(2) - dx(3)
                     dddzia = dx(1)*dadzia + dy(1)*dbdzia + c
                     dsdzia = 2.0d0 * (a*dadzia + b*dbdzia)
                     ddzdzia(1) = 1.0d0
                     ddzdzia(2) = 0.0d0
                     ddzdzia(3) = 0.0d0
                     ddxxdzia = (xp*dsdzia-a*dddzia-d*dadzia)/s
                     ddyydzia = (yp*dsdzia-b*dddzia-d*dbdzia)/s
                     ddzzdzia(1) = (zp*dsdzia-c*dddzia)/s + 1.0d0
                     ddzzdzia(2) = (zp*dsdzia-c*dddzia)/s
                     ddzzdzia(3) = (zp*dsdzia-c*dddzia)/s
                     ddsdzia(1) = dz(1) / ds(1)
                     ddsdzia(2) = 0.0d0
                     ddsdzia(3) = 0.0d0
                     ddssdzia(1) = ddxxdzia*dxxss(1) +
     &                    ddyydzia*dyyss(1) + ddzzdzia(1)*dzzss(1)
                     ddssdzia(2) = ddxxdzia*dxxss(2) +
     &                    ddyydzia*dyyss(2) + ddzzdzia(2)*dzzss(2)
                     ddssdzia(3) = ddxxdzia*dxxss(3) +
     &                    ddyydzia*dyyss(3) + ddzzdzia(3)*dzzss(3)
c
                     dbdxib = dz(1) - dz(3)
                     dcdxib = dy(3) - dy(1)
                     dddxib = dy(1)*dbdxib + dz(1)*dcdxib
                     dsdxib = 2.0d0 * (b*dbdxib + c*dcdxib)
                     ddxdxib(1) = 0.0d0
                     ddxdxib(2) = 1.0d0
                     ddxdxib(3) = 0.0d0
                     ddxxdxib(1) = (xp*dsdxib-a*dddxib)/s
                     ddxxdxib(2) = (xp*dsdxib-a*dddxib)/s + 1.0d0
                     ddxxdxib(3) = (xp*dsdxib-a*dddxib)/s
                     ddyydxib = (yp*dsdxib-b*dddxib-d*dbdxib)/s
                     ddzzdxib = (zp*dsdxib-c*dddxib-d*dcdxib)/s
                     ddsdxib(1) = 0.0d0
                     ddsdxib(2) = dx(2) / ds(2)
                     ddsdxib(3) = 0.0d0
                     ddssdxib(1) = ddxxdxib(1)*dxxss(1) +
     &                    ddyydxib*dyyss(1) + ddzzdxib*dzzss(1)
                     ddssdxib(2) = ddxxdxib(2)*dxxss(2) +
     &                    ddyydxib*dyyss(2) + ddzzdxib*dzzss(2)
                     ddssdxib(3) = ddxxdxib(3)*dxxss(3) +
     &                    ddyydxib*dyyss(3) + ddzzdxib*dzzss(3)
                     dadyib = dz(3) - dz(1)
                     dcdyib = dx(1) - dx(3)
                     dddyib = dx(1)*dadyib + dz(1)*dcdyib
                     dsdyib = 2.0d0 * (a*dadyib + c*dcdyib)
                     ddydyib(1) = 0.0d0
                     ddydyib(2) = 1.0d0
                     ddydyib(3) = 0.0d0
                     ddxxdyib = (xp*dsdyib-a*dddyib-d*dadyib)/s
                     ddyydyib(1) = (yp*dsdyib-b*dddyib)/s
                     ddyydyib(2) = (yp*dsdyib-b*dddyib)/s + 1.0d0
                     ddyydyib(3) = (yp*dsdyib-b*dddyib)/s
                     ddzzdyib = (zp*dsdyib-c*dddyib-d*dcdyib)/s
                     ddsdyib(1) = 0.0d0
                     ddsdyib(2) = dy(2) / ds(2)
                     ddsdyib(3) = 0.0d0
                     ddssdyib(1) = ddxxdyib*dxxss(1) +
     &                    ddyydyib(1)*dyyss(1) + ddzzdyib*dzzss(1)
                     ddssdyib(2) = ddxxdyib*dxxss(2) +
     &                    ddyydyib(2)*dyyss(2) + ddzzdyib*dzzss(2)
                     ddssdyib(3) = ddxxdyib*dxxss(3) +
     &                    ddyydyib(3)*dyyss(3) + ddzzdyib*dzzss(3)
                     dadzib = dy(1) - dy(3)
                     dbdzib = dx(3) - dx(1)
                     dddzib = dx(1)*dadzib + dy(1)*dbdzib
                     dsdzib = 2.0d0 * (a*dadzib + b*dbdzib)
                     ddzdzib(1) = 0.0d0
                     ddzdzib(2) = 1.0d0
                     ddzdzib(3) = 0.0d0
                     ddxxdzib = (xp*dsdzib-a*dddzib-d*dadzib)/s
                     ddyydzib = (yp*dsdzib-b*dddzib-d*dbdzib)/s
                     ddzzdzib(1) = (zp*dsdzib-c*dddzib)/s
                     ddzzdzib(2) = (zp*dsdzib-c*dddzib)/s + 1.0d0
                     ddzzdzib(3) = (zp*dsdzib-c*dddzib)/s
                     ddsdzib(1) = 0.0d0
                     ddsdzib(2) = dz(2) / ds(2)
                     ddsdzib(3) = 0.0d0
                     ddssdzib(1) = ddxxdzib*dxxss(1) +
     &                    ddyydzib*dyyss(1) + ddzzdzib(1)*dzzss(1)
                     ddssdzib(2) = ddxxdzib*dxxss(2) +
     &                    ddyydzib*dyyss(2) + ddzzdzib(2)*dzzss(2)
                     ddssdzib(3) = ddxxdzib*dxxss(3) +
     &                    ddyydzib*dyyss(3) + ddzzdzib(3)*dzzss(3)
c
c                    dbdxic = dz(2) - dz(1)
c                    dcdxic = dy(1) - dy(2)
c                    dddxic = dy(1)*dbdxic + dz(1)*dcdxic
c                    dsdxic = 2.0d0 * (b*dbdxic + c*dcdxic)
c                    ddxdxic(1) = 0.0d0
c                    ddxdxic(2) = 0.0d0
c                    ddxdxic(3) = 1.0d0
c                    ddxxdxic(1) = (xp*dsdxic-a*dddxic)/s
c                    ddxxdxic(2) = (xp*dsdxic-a*dddxic)/s
c                    ddxxdxic(3) = (xp*dsdxic-a*dddxic)/s + 1.0d0
c                    ddyydxic = (yp*dsdxic-b*dddxic-d*dbdxic)/s
c                    ddzzdxic = (zp*dsdxic-c*dddxic-d*dcdxic)/s
c                    ddsdxic(1) = 0.0d0
c                    ddsdxic(2) = 0.0d0
c                    ddsdxic(3) = dx(3) / ds(3)
c                    ddssdxic(1) = ddxxdxic(1)*dxxss(1) +
c    &                    ddyydxic*dyyss(1) + ddzzdxic*dzzss(1)
c                    ddssdxic(2) = ddxxdxic(2)*dxxss(2) +
c    &                    ddyydxic*dyyss(2) + ddzzdxic*dzzss(2)
c                    ddssdxic(3) = ddxxdxic(3)*dxxss(3) +
c    &                    ddyydxic*dyyss(3) + ddzzdxic*dzzss(3)
c                    dadyic = dz(1) - dz(2)
c                    dcdyic = dx(2) - dx(1)
c                    dddyic = dx(1)*dadyic + dz(1)*dcdyic
c                    dsdyic = 2.0d0 * (a*dadyic + c*dcdyic)
c                    ddydyic(1) = 0.0d0
c                    ddydyic(2) = 0.0d0
c                    ddydyic(3) = 1.0d0
c                    ddxxdyic = (xp*dsdyic-a*dddyic-d*dadyic)/s
c                    ddyydyic(1) = (yp*dsdyic-b*dddyic)/s
c                    ddyydyic(2) = (yp*dsdyic-b*dddyic)/s
c                    ddyydyic(3) = (yp*dsdyic-b*dddyic)/s + 1.0d0
c                    ddzzdyic = (zp*dsdyic-c*dddyic-d*dcdyic)/s
c                    ddsdyic(1) = 0.0d0
c                    ddsdyic(2) = 0.0d0
c                    ddsdyic(3) = dy(3) / ds(3)
c                    ddssdyic(1) = ddxxdyic*dxxss(1) +
c    &                    ddyydyic(1)*dyyss(1) + ddzzdyic*dzzss(1)
c                    ddssdyic(2) = ddxxdyic*dxxss(2) +
c    &                    ddyydyic(2)*dyyss(2) + ddzzdyic*dzzss(2)
c                    ddssdyic(3) = ddxxdyic*dxxss(3) +
c    &                    ddyydyic(3)*dyyss(3) + ddzzdyic*dzzss(3)
c                    dadzic = dy(2) - dy(1)
c                    dbdzic = dx(1) - dx(2)
c                    dddzic = dx(1)*dadzic + dy(1)*dbdzic
c                    dsdzic = 2.0d0 * (a*dadzic + b*dbdzic)
c                    ddzdzic(1) = 0.0d0
c                    ddzdzic(2) = 0.0d0
c                    ddzdzic(3) = 1.0d0
c                    ddxxdzic = (xp*dsdzic-a*dddzic-d*dadzic)/s
c                    ddyydzic = (yp*dsdzic-b*dddzic-d*dbdzic)/s
c                    ddzzdzic(1) = (zp*dsdzic-c*dddzic)/s
c                    ddzzdzic(2) = (zp*dsdzic-c*dddzic)/s
c                    ddzzdzic(3) = (zp*dsdzic-c*dddzic)/s + 1.0d0
c                    ddsdzic(1) = 0.0d0
c                    ddsdzic(2) = 0.0d0
c                    ddsdzic(3) = dz(3) / ds(3)
c                    ddssdzic(1) = ddxxdzic*dxxss(1) +
c    &                    ddyydzic*dyyss(1) + ddzzdzic(1)*dzzss(1)
c                    ddssdzic(2) = ddxxdzic*dxxss(2) +
c    &                    ddyydzic*dyyss(2) + ddzzdzic(2)*dzzss(2)
c                    ddssdzic(3) = ddxxdzic*dxxss(3) +
c    &                    ddyydzic*dyyss(3) + ddzzdzic(3)*dzzss(3)
                  end if
c
                  dcosodxi  = dx(k)*ddxxdxi - dxx(k) +
     &                        dy(k)*ddyydxi + dz(k)*ddzzdxi -
     &                  coso(k)*(ds(k)*ddssdxi(k)+dss(k)*ddsdxi(k))
                  dcosodyi  = dx(k)*ddxxdyi + dy(k)*ddyydyi -
     &                        dyy(k) + dz(k)*ddzzdyi  -
     &                  coso(k)*(ds(k)*ddssdyi(k)+dss(k)*ddsdyi(k))
                  dcosodzi  = dx(k)*ddxxdzi + dy(k)*ddyydzi +
     &                        dz(k)*ddzzdzi - dzz(k) -
     &                  coso(k)*(ds(k)*ddssdzi(k)+dss(k)*ddsdzi(k))
                  dcosodxia = dx(k)*ddxxdxia(k) + dxx(k)*ddxdxia(k) +
     &                        dy(k)*ddyydxia + dz(k)*ddzzdxia -
     &                  coso(k)*(ds(k)*ddssdxia(k)+dss(k)*ddsdxia(k))
                  dcosodyia = dx(k)*ddxxdyia + dy(k)*ddyydyia(k) +
     &                        dyy(k)*ddydyia(k) + dz(k)*ddzzdyia -
     &                  coso(k)*(ds(k)*ddssdyia(k)+dss(k)*ddsdyia(k))
                  dcosodzia = dx(k)*ddxxdzia + dy(k)*ddyydzia +
     &                        dz(k)*ddzzdzia(k) + dzz(k)*ddzdzia(k) -
     &                  coso(k)*(ds(k)*ddssdzia(k)+dss(k)*ddsdzia(k))
                  dcosodxib = dx(k)*ddxxdxib(k) + dxx(k)*ddxdxib(k) +
     &                        dy(k)*ddyydxib + dz(k)*ddzzdxib -
     &                  coso(k)*(ds(k)*ddssdxib(k)+dss(k)*ddsdxib(k))
                  dcosodyib = dx(k)*ddxxdyib + dy(k)*ddyydyib(k) +
     &                        dyy(k)*ddydyib(k) + dz(k)*ddzzdyib -
     &                  coso(k)*(ds(k)*ddssdyib(k)+dss(k)*ddsdyib(k))
                  dcosodzib = dx(k)*ddxxdzib + dy(k)*ddyydzib +
     &                        dz(k)*ddzzdzib(k) + dzz(k)*ddzdzib(k) -
     &                  coso(k)*(ds(k)*ddssdzib(k)+dss(k)*ddsdzib(k))
c                 dcosodxic = dx(k)*ddxxdxic(k) + dxx(k)*ddxdxic(k) +
c    &                        dy(k)*ddyydxic + dz(k)*ddzzdxic -
c    &                  coso(k)*(ds(k)*ddssdxic(k)+dss(k)*ddsdxic(k))
c                 dcosodyic = dx(k)*ddxxdyic + dy(k)*ddyydyic(k) +
c    &                        dyy(k)*ddydyic(k) + dz(k)*ddzzdyic -
c    &                  coso(k)*(ds(k)*ddssdyic(k)+dss(k)*ddsdyic(k))
c                 dcosodzic = dx(k)*ddxxdzic + dy(k)*ddyydzic +
c    &                        dz(k)*ddzzdzic(k) + dzz(k)*ddzdzic(k) -
c    &                  coso(k)*(ds(k)*ddssdzic(k)+dss(k)*ddsdzic(k))
                  dcosodxic = -dcosodxi - dcosodxia - dcosodxib
                  dcosodyic = -dcosodyi - dcosodyia - dcosodyib
                  dcosodzic = -dcosodzi - dcosodzia - dcosodzib
c
c     out-of-plane bend calculation; out-of-plane portion
c             hessian diagonal block elements
c
                  hea(1,i,1,i) = hea(1,i,1,i)  !+ deodcoso*doxixi
     &                              + d2eodcoso2*dcosodxi*dcosodxi
                  hea(1,i,2,i) = hea(1,i,2,i)  !+ deodcoso*doxiyi
     &                              + d2eodcoso2*dcosodxi*dcosodyi
                  hea(1,i,3,i) = hea(1,i,3,i)  !+ deodcoso*doxizi
     &                              + d2eodcoso2*dcosodxi*dcosodzi
                  hea(2,i,2,i) = hea(2,i,2,i)  !+ deodcoso*doyiyi
     &                              + d2eodcoso2*dcosodyi*dcosodyi
                  hea(2,i,3,i) = hea(2,i,3,i)  !+ deodcoso*doyizi
     &                              + d2eodcoso2*dcosodyi*dcosodzi
                  hea(3,i,3,i) = hea(3,i,3,i)  !+ deodcoso*dozizi
     &                              + d2eodcoso2*dcosodzi*dcosodzi
c
                  hea(1,ia,1,ia) = hea(1,ia,1,ia)  !+ deodcoso*doxiaxia
     &                                + d2eodcoso2*dcosodxia*dcosodxia
                  hea(1,ia,2,ia) = hea(1,ia,2,ia)  !+ deodcoso*doxiayia
     &                                + d2eodcoso2*dcosodxia*dcosodyia
                  hea(1,ia,3,ia) = hea(1,ia,3,ia)  !+ deodcoso*doxiazia
     &                                + d2eodcoso2*dcosodxia*dcosodzia
                  hea(2,ia,2,ia) = hea(2,ia,2,ia)  !+ deodcoso*doyiayia
     &                                + d2eodcoso2*dcosodyia*dcosodyia
                  hea(2,ia,3,ia) = hea(2,ia,3,ia)  !+ deodcoso*doyiazia
     &                                + d2eodcoso2*dcosodyia*dcosodzia
                  hea(3,ia,3,ia) = hea(3,ia,3,ia)  !+ deodcoso*doziazia
     &                                + d2eodcoso2*dcosodzia*dcosodzia
c
                  hea(1,ib,1,ib) = hea(1,ib,1,ib)  !+ deodcoso*doxibxib
     &                                + d2eodcoso2*dcosodxib*dcosodxib
                  hea(1,ib,2,ib) = hea(1,ib,2,ib)  !+ deodcoso*doxibyib
     &                                + d2eodcoso2*dcosodxib*dcosodyib
                  hea(1,ib,3,ib) = hea(1,ib,3,ib)  !+ deodcoso*doxibzib
     &                                + d2eodcoso2*dcosodxib*dcosodzib
                  hea(2,ib,2,ib) = hea(2,ib,2,ib)  !+ deodcoso*doyibyib
     &                                + d2eodcoso2*dcosodyib*dcosodyib
                  hea(2,ib,3,ib) = hea(2,ib,3,ib)  !+ deodcoso*doyibzib
     &                                + d2eodcoso2*dcosodyib*dcosodzib
                  hea(3,ib,3,ib) = hea(3,ib,3,ib)  !+ deodcoso*dozibzib
     &                                + d2eodcoso2*dcosodzib*dcosodzib
c
                  hea(1,ic,1,ic) = hea(1,ic,1,ic)  !+ deodcoso*doxicxic
     &                                + d2eodcoso2*dcosodxic*dcosodxic
                  hea(1,ic,2,ic) = hea(1,ic,2,ic)  !+ deodcoso*doxicyic
     &                                + d2eodcoso2*dcosodxic*dcosodyic
                  hea(1,ic,3,ic) = hea(1,ic,3,ic)  !+ deodcoso*doxiczic
     &                                + d2eodcoso2*dcosodxic*dcosodzic
                  hea(2,ic,2,ic) = hea(2,ic,2,ic)  !+ deodcoso*doyicyic
     &                                + d2eodcoso2*dcosodyic*dcosodyic
                  hea(2,ic,3,ic) = hea(2,ic,3,ic)  !+ deodcoso*doyiczic
     &                                + d2eodcoso2*dcosodyic*dcosodzic
                  hea(3,ic,3,ic) = hea(3,ic,3,ic)  !+ deodcoso*doziczic
     &                                + d2eodcoso2*dcosodzic*dcosodzic
c
c     out-of-plane bend calculation; in-the-plane portion
c
                  dt = angin(k) - anat(kang)
                  dt2 = dt * dt
                  deiddt = 0.02191418d0 * acon(kang) * dt
     &                       * (2.0d0 + 0.00006d0*sf*dt2*dt2)
                  d2eiddt2 = 0.02191418d0 * acon(kang)
     &                          * (2.0d0 + 0.0003d0*sf*dt2*dt2)
                  sine = dsqrt(1.0d0-cosi(k)*cosi(k))
                  ddtdcosi = -radian / (normi(k)*sine)
                  d2dtcosi2 = ddtdcosi * cosi(k)/(normi(k)*sine**2)
                  deidcosi = deiddt * ddtdcosi
                  d2eidcosi2 = d2eiddt2*ddtdcosi + deiddt*d2dtcosi2
c
                  dcosidxi  = (dxx(l)+dxx(m))*ddxxdxi +
     &                        (dyy(l)+dyy(m))*ddyydxi +
     &                        (dzz(l)+dzz(m))*ddzzdxi -
     &               cosi(k)*(dss(l)*ddssdxi(m)+dss(m)*ddssdxi(l))
                  dcosidyi  = (dxx(l)+dxx(m))*ddxxdyi +
     &                        (dyy(l)+dyy(m))*ddyydyi +
     &                        (dzz(l)+dzz(m))*ddzzdyi -
     &               cosi(k)*(dss(l)*ddssdyi(m)+dss(m)*ddssdyi(l))
                  dcosidzi  = (dxx(l)+dxx(m))*ddxxdzi +
     &                        (dyy(l)+dyy(m))*ddyydzi +
     &                        (dzz(l)+dzz(m))*ddzzdzi -
     &               cosi(k)*(dss(l)*ddssdzi(m)+dss(m)*ddssdzi(l))
                  dcosidxia = dxx(l)*ddxxdxia(m) + dxx(m)*ddxxdxia(l) +
     &                        (dyy(l)+dyy(m))*ddyydxia +
     &                        (dzz(l)+dzz(m))*ddzzdxia -
     &               cosi(k)*(dss(l)*ddssdxia(m)+dss(m)*ddssdxia(l))
                  dcosidyia = (dxx(l)+dxx(m))*ddxxdyia +
     &                        dyy(l)*ddyydyia(m) + dyy(m)*ddyydyia(l) +
     &                        (dzz(l)+dzz(m))*ddzzdyia -
     &               cosi(k)*(dss(l)*ddssdyia(m)+dss(m)*ddssdyia(l))
                  dcosidzia = (dxx(l)+dxx(m))*ddxxdzia +
     &                        (dyy(l)+dyy(m))*ddyydzia +
     &                        dzz(l)*ddzzdzia(m) + dzz(m)*ddzzdzia(l) -
     &               cosi(k)*(dss(l)*ddssdzia(m)+dss(m)*ddssdzia(l))
                  dcosidxib = dxx(l)*ddxxdxib(m) + dxx(m)*ddxxdxib(l) +
     &                        (dyy(l)+dyy(m))*ddyydxib +
     &                        (dzz(l)+dzz(m))*ddzzdxib -
     &               cosi(k)*(dss(l)*ddssdxib(m)+dss(m)*ddssdxib(l))
                  dcosidyib = (dxx(l)+dxx(m))*ddxxdyib +
     &                        dyy(l)*ddyydyib(m) + dyy(m)*ddyydyib(l) +
     &                        (dzz(l)+dzz(m))*ddzzdyib -
     &               cosi(k)*(dss(l)*ddssdyib(m)+dss(m)*ddssdyib(l))
                  dcosidzib = (dxx(l)+dxx(m))*ddxxdzib +
     &                        (dyy(l)+dyy(m))*ddyydzib +
     &                        dzz(l)*ddzzdzib(m) + dzz(m)*ddzzdzib(l) -
     &               cosi(k)*(dss(l)*ddssdzib(m)+dss(m)*ddssdzib(l))
c                 dcosidxic = dxx(l)*ddxxdxic(m) + dxx(m)*ddxxdxic(l) +
c    &                        (dyy(l)+dyy(m))*ddyydxic +
c    &                        (dzz(l)+dzz(m))*ddzzdxic -
c    &               cosi(k)*(dss(l)*ddssdxic(m)+dss(m)*ddssdxic(l))
c                 dcosidyic = (dxx(l)+dxx(m))*ddxxdyic +
c    &                        dyy(l)*ddyydyic(m) + dyy(m)*ddyydyic(l) +
c    &                        (dzz(l)+dzz(m))*ddzzdyic -
c    &               cosi(k)*(dss(l)*ddssdyic(m)+dss(m)*ddssdyic(l))
c                 dcosidzic = (dxx(l)+dxx(m))*ddxxdzic +
c    &                        (dyy(l)+dyy(m))*ddyydzic +
c    &                        dzz(l)*ddzzdzic(m) + dzz(m)*ddzzdzic(l) -
c    &               cosi(k)*(dss(l)*ddssdzic(m)+dss(m)*ddssdzic(l))
                  dcosidxic = -dcosidxi - dcosidxia - dcosidxib
                  dcosidyic = -dcosidyi - dcosidyia - dcosidyib
                  dcosidzic = -dcosidzi - dcosidzia - dcosidzib
c
c     out-of-plane bend calculation; in-the-plane portion
c             hessian diagonal block elements
c
                  hea(1,i,1,i) = hea(1,i,1,i)  !+ deidcosi*dixixi
     &                              + d2eidcosi2*dcosidxi*dcosidxi
                  hea(1,i,2,i) = hea(1,i,2,i)  !+ deidcosi*dixiyi
     &                              + d2eidcosi2*dcosidxi*dcosidyi
                  hea(1,i,3,i) = hea(1,i,3,i)  !+ deidcosi*dixizi
     &                              + d2eidcosi2*dcosidxi*dcosidzi
                  hea(2,i,2,i) = hea(2,i,2,i)  !+ deidcosi*diyiyi
     &                              + d2eidcosi2*dcosidyi*dcosidyi
                  hea(2,i,3,i) = hea(2,i,3,i)  !+ deidcosi*diyizi
     &                              + d2eidcosi2*dcosidyi*dcosidzi
                  hea(3,i,3,i) = hea(3,i,3,i)  !+ deidcosi*dizizi
     &                              + d2eidcosi2*dcosidzi*dcosidzi
c
                  hea(1,ia,1,ia) = hea(1,ia,1,ia)  !+ deidcosi*dixiaxia
     &                                + d2eidcosi2*dcosidxia*dcosidxia
                  hea(1,ia,2,ia) = hea(1,ia,2,ia)  !+ deidcosi*dixiayia
     &                                + d2eidcosi2*dcosidxia*dcosidyia
                  hea(1,ia,3,ia) = hea(1,ia,3,ia)  !+ deidcosi*dixiazia
     &                                + d2eidcosi2*dcosidxia*dcosidzia
                  hea(2,ia,2,ia) = hea(2,ia,2,ia)  !+ deidcosi*diyiayia
     &                                + d2eidcosi2*dcosidyia*dcosidyia
                  hea(2,ia,3,ia) = hea(2,ia,3,ia)  !+ deidcosi*diyiazia
     &                                + d2eidcosi2*dcosidyia*dcosidzia
                  hea(3,ia,3,ia) = hea(3,ia,3,ia)  !+ deidcosi*diziazia
     &                                + d2eidcosi2*dcosidzia*dcosidzia
c
                  hea(1,ib,1,ib) = hea(1,ib,1,ib)  !+ deidcosi*dixibxib
     &                                + d2eidcosi2*dcosidxib*dcosidxib
                  hea(1,ib,2,ib) = hea(1,ib,2,ib)  !+ deidcosi*dixibyib
     &                                + d2eidcosi2*dcosidxib*dcosidyib
                  hea(1,ib,3,ib) = hea(1,ib,3,ib)  !+ deidcosi*dixibzib
     &                                + d2eidcosi2*dcosidxib*dcosidzib
                  hea(2,ib,2,ib) = hea(2,ib,2,ib)  !+ deidcosi*diyibyib
     &                                + d2eidcosi2*dcosidyib*dcosidyib
                  hea(2,ib,3,ib) = hea(2,ib,3,ib)  !+ deidcosi*diyibzib
     &                                + d2eidcosi2*dcosidyib*dcosidzib
                  hea(3,ib,3,ib) = hea(3,ib,3,ib)  !+ deidcosi*dizibzib
     &                                + d2eidcosi2*dcosidzib*dcosidzib
c
                  hea(1,ic,1,ic) = hea(1,ic,1,ic)  !+ deidcosi*dixicxic
     &                                + d2eidcosi2*dcosidxic*dcosidxic
                  hea(1,ic,2,ic) = hea(1,ic,2,ic)  !+ deidcosi*dixicyic
     &                                + d2eidcosi2*dcosidxic*dcosidyic
                  hea(1,ic,3,ic) = hea(1,ic,3,ic)  !+ deidcosi*dixiczic
     &                                + d2eidcosi2*dcosidxic*dcosidzic
                  hea(2,ic,2,ic) = hea(2,ic,2,ic)  !+ deidcosi*diyicyic
     &                                + d2eidcosi2*dcosidyic*dcosidyic
                  hea(2,ic,3,ic) = hea(2,ic,3,ic)  !+ deidcosi*diyiczic
     &                                + d2eidcosi2*dcosidyic*dcosidzic
                  hea(3,ic,3,ic) = hea(3,ic,3,ic)  !+ deidcosi*diziczic
     &                                + d2eidcosi2*dcosidzic*dcosidzic
               end if
c
c     increment the angle definition counters
c
               m = m + 1
               if (m .gt. n12(i)) then
                  l = l + 1
                  m = l + 1
               end if
            end do
         end if
      end do
      return
      end

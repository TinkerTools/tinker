c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ebond4  --  bond stretch Hessian; full matrix  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ebond4" calculates the full second derivative matrix of the
c     bond stretching energy with respect to cartesian coordinates
c
c
      subroutine ebond4
      implicit none
      include 'sizes.i'
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
         hdiag(1,i) = hdiag(1,i) + d2e(1,1)
         hdiag(2,i) = hdiag(2,i) + d2e(2,2)
         hdiag(3,i) = hdiag(3,i) + d2e(3,3)
         hdiag(1,k) = hdiag(1,k) + d2e(1,1)
         hdiag(2,k) = hdiag(2,k) + d2e(2,2)
         hdiag(3,k) = hdiag(3,k) + d2e(3,3)

         j = n*(3*i-3) - (3*i-1)*(3*i-2)/2 + 3*i - 1
         h(j) = h(j) + d2e(1,2)
         h(j+1) = h(j+1) + d2e(1,3)
         h(j+n) = h(j+n) + d2e(2,3)

         j = n*(3*k-3) - (3*k-1)*(3*k-2)/2 + 3*k - 1
         h(j) = h(j) + d2e(1,2)
         h(j+1) = h(j+1) + d2e(1,3)
         h(j+n) = h(j+n) + d2e(2,3)

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

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program cmapedit  --  interconvert CMAP parameter formats  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "cmapedit" reads CMAP parameters in Amber format, as for the
c     ff19SB force field, and then outputs TORTORS parameters in the
c     format used by Tinker
c
c     note Amber to Tinker conversion assumes a 15 degree grid in both
c     torsions containing 24 x 24 = 576 Amber grid values, which then
c     converts to 25 x 25 = 625 grid values for use with Tinker
c
c
      program cmapedit
      implicit none
      integer i,j,k,m,n
      integer input,iout
      real*8 tamber(625)
      real*8 tx(625)
      real*8 ty(625)
      real*8 tf(625)
c
c
c     read the Amber-formatted CMAP torsion-torsion values
c
      input = 5
      iout = 6
      open (unit=input, file='cmap.txt', status='old')
      k = 0
      do i = 1, 72
         read (input,10)  (tamber(j),j=k+1,k+8)
   10    format (f9.5,7f10.5)
         k = k + 8
      end do
c
c     store the Tinker-formatted CMAP torsion-torsion values
c
      n = 0
      k = 0
      do i = -180, 165, 15
         k = k + 1
         m = 24 * (k-1)
         do j = -180, 165, 15
            m = m + 1
            n = n + 1
            tx(n) = dble(i)
            ty(n) = dble(j)
            tf(n) = tamber(m)
         end do
         m = m + 1
         n = n + 1
         tx(n) = dble(i)
         ty(n) = 180.0d0
         tf(n) = tamber(m-24)
      end do
      k = k + 1
      m = 24 * (k-1)
      do j = -180, 165, 15
         m = m + 1
         n = n + 1
         tx(n) = 180.0d0
         ty(n) = dble(j)
         tf(n) = tamber(m-576)
      end do
      m = m + 1
      n = n + 1
      tx(n) = 180.0d0
      ty(n) = 180.0d0
      tf(n) = tamber(m-576-24)
c
c     print the Tinker-formatted CMAP torsion-torsion values
c
      i = 0
      dowhile (i .lt. 624)
         write (iout,20)  tx(i+1),ty(i+1),tf(i+1),
     &                    tx(i+2),ty(i+2),tf(i+2),
     &                    tx(i+3),ty(i+3),tf(i+3)
   20    format (1x,2f7.1,f10.5,1x,2f7.1,f10.5,1x,2f7.1,f10.5)
         i = i + 3
      end do
      write (iout,30)  tx(625),ty(625),tf(625)
   30 format (1x,2f7.1,f10.5)
      end

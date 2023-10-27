c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine chkpole  --  check multipoles at chiral sites  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "chkpole" inverts multipole moments as necessary at atoms
c     with chiral local reference frame definitions
c
c
      subroutine chkpole
      use atoms
      use mpole
      use repel
      use xrepel
      implicit none
      integer i,k
      integer ia,ib,ic,id
      real*8 xad,yad,zad
      real*8 xbd,ybd,zbd
      real*8 xcd,ycd,zcd
      real*8 c1,c2,c3,vol
      logical dopol,dorep
      logical doxrep
      logical check
c
c
c     loop over multipoles and test for chirality inversion
c
      do i = 1, n
         dopol = .false.
         dorep = .false.
         if (allocated(pollist)) then
            if (pollist(i) .ne. 0)  dopol = .true.
         end if
         if (allocated(replist)) then
            if (replist(i) .ne. 0)  dorep = .true.
         else if (allocated(xreplist)) then
            if (xreplist(i) .ne. 0)  doxrep = .true.
         end if
         if (dopol .or. dorep .or. doxrep) then
            check = .true.
            if (polaxe(i) .ne. 'Z-then-X')  check = .false.
            if (yaxis(i) .eq. 0)  check = .false.
            if (check) then
               k = yaxis(i)
               ia = i
               ib = zaxis(i)
               ic = xaxis(i)
               id = abs(k)
c
c     compute the signed parallelpiped volume at chiral site
c
               xad = x(ia) - x(id)
               yad = y(ia) - y(id)
               zad = z(ia) - z(id)
               xbd = x(ib) - x(id)
               ybd = y(ib) - y(id)
               zbd = z(ib) - z(id)
               xcd = x(ic) - x(id)
               ycd = y(ic) - y(id)
               zcd = z(ic) - z(id)
               c1 = ybd*zcd - zbd*ycd
               c2 = ycd*zad - zcd*yad
               c3 = yad*zbd - zad*ybd
               vol = xad*c1 + xbd*c2 + xcd*c3
c
c     invert the multipole components involving the y-axis
c
               if ((k.lt.0.and.vol.gt.0.0d0) .or.
     &             (k.gt.0.and.vol.lt.0.0d0)) then
                  yaxis(i) = -k
                  if (dopol) then
                     pole(3,i) = -pole(3,i)
                     pole(6,i) = -pole(6,i)
                     pole(8,i) = -pole(8,i)
                     pole(10,i) = -pole(10,i)
                     pole(12,i) = -pole(12,i)
                  end if
                  if (dorep) then
                     repole(3,i) = -repole(3,i)
                     repole(6,i) = -repole(6,i)
                     repole(8,i) = -repole(8,i)
                     repole(10,i) = -repole(10,i)
                     repole(12,i) = -repole(12,i)
                  else if (doxrep) then
                     xrepole(3,i) = -xrepole(3,i)
                  end if
               end if
            end if
         end if
      end do
      return
      end

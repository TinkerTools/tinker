c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kstrbnd  --  assign stretch-bend parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kstrbnd" assigns the parameters for the stretch-bend
c     interactions and processes new or changed parameter values
c
c
      subroutine kstrbnd
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kstbnd.i'
      include 'potent.i'
      include 'strbnd.i'
      integer i,j,k,nsb
      integer ia,ib,ic
      integer ita,itb,itc
      integer nba,nbc
      integer size,next
      real*8 sb1,sb2,temp
      logical header
      character*4 pa,pb,pc
      character*12 blank,pt
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing stretch-bend parameters
c
      blank = '            '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'STRBND ') then
            ia = 0
            ib = 0
            ic = 0
            sb1 = 0.0d0
            sb2 = 0.0d0
            string = record(next:120)
            read (string,*,err=10,end=10)  ia,ib,ic,sb1,sb2
   10       continue
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Stretch-Bend Parameters :',
     &                 //,5x,'Atom Classes',6x,'K(SB)-1',5x,'K(SB)-2',/)
            end if
            write (iout,30)  ia,ib,ic,sb1,sb2
   30       format (4x,3i4,2x,2f12.3)
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
               temp = sb1
               sb1 = sb2
               sb2 = temp
            end if
            do j = 1, maxnsb
               if (ksb(j).eq.blank .or. ksb(j).eq.pt) then
                  ksb(j) = pt
                  stbn(1,j) = sb1
                  stbn(2,j) = sb2
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KSTRBND  --  Too many Stretch-Bend',
     &                 ' Interaction Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nsb = maxnsb
      do i = maxnsb, 1, -1
         if (ksb(i) .eq. blank)  nsb = i - 1
      end do
c
c     assign the stretch-bend parameters for each angle
c
      nstrbnd = 0
      if (nsb .ne. 0) then
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            ita = class(ia)
            itb = class(ib)
            itc = class(ic)
            size = 4
            call numeral (ita,pa,size)
            call numeral (itb,pb,size)
            call numeral (itc,pc,size)
            if (ita .le. itc) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            do j = 1, nsb
               if (ksb(j) .eq. pt) then
                  nstrbnd = nstrbnd + 1
                  do k = 1, n12(ib)
                     if (i12(k,ib) .eq. ia)  nba = bndlist(k,ib)
                     if (i12(k,ib) .eq. ic)  nbc = bndlist(k,ib)
                  end do
                  isb(1,nstrbnd) = i
                  isb(2,nstrbnd) = nba
                  isb(3,nstrbnd) = nbc
                  if (ita .le. itc) then
                     sbk(1,nstrbnd) = stbn(1,j)
                     sbk(2,nstrbnd) = stbn(2,j)
                  else
                     sbk(1,nstrbnd) = stbn(2,j)
                     sbk(2,nstrbnd) = stbn(1,j)
                  end if
                  goto 60
               end if
            end do
   60       continue
         end do
      end if
c
c     turn off the stretch-bend potential if it is not used
c
      if (nstrbnd .eq. 0)  use_strbnd = .false.
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kcharge  --  assign partial charge parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kcharge" assigns partial charges to the atoms within
c     the structure and processes any new or changed values
c
c
      subroutine kcharge
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kchrge.i'
      include 'keys.i'
      include 'potent.i'
      integer i,j,k,m
      integer ia,next
      integer nc12(maxatm)
      real*8 cg
      logical header
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing partial charge parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            cg = 0.0d0
            string = record(next:120)
            read (string,*,err=40,end=40)  ia,cg
            if (ia .gt. 0) then
               if (header) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Additional Atomic Partial Charge',
     &                       ' Parameters :',
     &                    //,5x,'Atom Type',10x,'Charge',/)
               end if
               if (ia .le. maxtyp) then
                  chg(ia) = cg
                  write (iout,20)  ia,cg
   20             format (4x,i6,8x,f12.4)
               else
                  write (iout,30)
   30             format (/,' KCHARGE  --  Too many Partial Charge',
     &                       ' Parameters')
                  abort = .true.
               end if
            end if
   40       continue
         end if
      end do
c
c     find and store all the atomic partial charges
c
      do i = 1, n
         pchg(i) = chg(type(i))
      end do
c
c     process keywords containing atom specific partial charges
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHARGE ') then
            ia = 0
            cg = 0.0d0
            string = record(next:120)
            read (string,*,err=70,end=70)  ia,cg
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Partial Charges for',
     &                       ' Specific Atoms :',
     &                    //,6x,'Atom',14x,'Charge',/)
               end if
               write (iout,60)  ia,cg
   60          format (4x,i6,8x,f12.4)
               pchg(ia) = cg
            end if
   70       continue
         end if
      end do
c
c     remove zero partial charges from the list of charges
c
      nion = 0
      do i = 1, n
         chglist(i) = 0
         if (pchg(i) .ne. 0.0d0) then
            nion = nion + 1
            iion(nion) = i
            jion(nion) = i
            kion(nion) = i
            pchg(nion) = pchg(i)
            chglist(i) = nion
         end if
      end do
c
c     optionally use neutral groups for neighbors and cutoffs
c
      if (neutnbr .or. neutcut) then
         do i = 1, n
            nc12(i) = 0
            do j = 1, n12(i)
               k = chglist(i12(j,i))
               if (k .ne. 0)  nc12(i) = nc12(i) + 1
            end do
         end do
         do i = 1, nion
            k = iion(i)
            if (n12(k) .eq. 1) then
               do j = 1, n12(k)
                  m = i12(j,k)
                  if (nc12(m) .gt. 1) then
                     if (neutnbr)  jion(i) = m
                     if (neutcut)  kion(i) = m
                  end if
               end do
            end if
         end do
      end if
c
c     turn off charge-charge and charge-dipole terms if not used
c
      if (nion .eq. 0) then
         use_charge = .false.
         use_chgdpl = .false.
      end if
      return
      end

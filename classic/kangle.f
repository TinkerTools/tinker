c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kangle  --  angle bend parameter assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kangle" assigns the force constants and ideal angles for
c     the bond angles; also processes new or changed parameters
c
c
      subroutine kangle
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kangs.i'
      include 'keys.i'
      include 'potent.i'
      include 'usage.i'
      integer i,j
      integer ia,ib,ic
      integer ita,itb,itc
      integer na,na5,na4
      integer na3,naf
      integer jen,ih,nh
      integer next,size
      integer minat,iring
      real*8 fc,an,pr
      real*8 an1,an2,an3
      logical header,done
      logical use_ring
      character*4 pa,pb,pc
      character*6 label
      character*12 blank,pt
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing angle bending parameters
c
      blank = '         '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:6) .eq. 'ANGLE ')  iring = 0
         if (keyword(1:7) .eq. 'ANGLE5 ')  iring = 5
         if (keyword(1:7) .eq. 'ANGLE4 ')  iring = 4
         if (keyword(1:7) .eq. 'ANGLE3 ')  iring = 3
         if (iring .ge. 0) then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an1 = 0.0d0
            an2 = 0.0d0
            an3 = 0.0d0
            jen = 0
            string = record(next:120)
            read (string,*,err=10,end=10)  ia,ib,ic,fc,an1,an2,an3
   10       continue
            if (an2.ne.0.0d0 .or. an3.ne.0.0d0)  jen = 1
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Angle Bending Parameters :',
     &                 //,5x,'Atom Classes',9x,'K(B)',7x,'Angle',/)
            end if
            if (iring .eq. 0) then
               if (jen .eq. 0) then
                  write (iout,30)  ia,ib,ic,fc,an1
   30             format (4x,3i4,2x,2f12.3)
               else if (an1 .ne. 0.0d0) then
                  write (iout,40)  ia,ib,ic,fc,an1
   40             format (4x,3i4,2x,2f12.3,3x,'0-H''s')
               end if
               if (an2 .ne. 0.0d0) then
                  write (iout,50)  ia,ib,ic,fc,an2
   50             format (4x,3i4,2x,2f12.3,3x,'1-H''s')
               end if
               if (an3 .ne. 0.0d0) then
                  write (iout,60)  ia,ib,ic,fc,an3
   60             format (4x,3i4,2x,2f12.3,3x,'2-H''s')
               end if
            else
               if (iring .eq. 5)  label = '5-Ring'
               if (iring .eq. 4)  label = '4-Ring'
               if (iring .eq. 3)  label = '3-Ring'
               if (jen .eq. 0) then
                  write (iout,70)  ia,ib,ic,fc,an1,label
   70             format (4x,3i4,2x,2f12.3,3x,a6)
               else if (an1 .ne. 0.0d0) then
                  write (iout,80)  ia,ib,ic,fc,an1,label
   80             format (4x,3i4,2x,2f12.3,3x,a6,3x,'0-H''s')
               end if
               if (an2 .ne. 0.0d0) then
                  write (iout,90)  ia,ib,ic,fc,an2,label
   90             format (4x,3i4,2x,2f12.3,3x,a6,3x,'1-H''s')
               end if
               if (an3 .ne. 0.0d0) then
                  write (iout,100)  ia,ib,ic,fc,an3,label
  100             format (4x,3i4,2x,2f12.3,3x,a6,3x,'2-H''s')
               end if
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            if (iring .eq. 0) then
               do j = 1, maxna
                  if (ka(j).eq.blank .or. ka(j).eq.pt) then
                     ka(j) = pt
                     acon(j) = fc
                     ang(1,j) = an1
                     ang(2,j) = an2
                     ang(3,j) = an3
                     goto 150
                  end if
               end do
               write (iout,110)
  110          format (/,' KANGLE  --  Too many Bond Angle',
     &                       ' Bending Parameters')
               abort = .true.
            else if (iring .eq. 5) then
               do j = 1, maxna5
                  if (ka5(j).eq.blank .or. ka5(j).eq.pt) then
                     ka5(j) = pt
                     acon5(j) = fc
                     ang5(1,j) = an1
                     ang5(2,j) = an2
                     ang5(3,j) = an3
                     goto 150
                  end if
               end do
               write (iout,120)
  120          format (/,' KANGLE  --  Too many 5-Ring Angle',
     &                       ' Bending Parameters')
               abort = .true.
            else if (iring .eq. 4) then
               do j = 1, maxna4
                  if (ka4(j).eq.blank .or. ka4(j).eq.pt) then
                     ka4(j) = pt
                     acon4(j) = fc
                     ang4(1,j) = an1
                     ang4(2,j) = an2
                     ang4(3,j) = an3
                     goto 150
                  end if
               end do
               write (iout,130)
  130          format (/,' KANGLE  --  Too many 4-Ring Angle',
     &                       ' Bending Parameters')
               abort = .true.
            else if (iring .eq. 3) then
               do j = 1, maxna3
                  if (ka3(j).eq.blank .or. ka3(j).eq.pt) then
                     ka3(j) = pt
                     acon3(j) = fc
                     ang3(1,j) = an1
                     ang3(2,j) = an2
                     ang3(3,j) = an3
                     goto 150
                  end if
               end do
               write (iout,140)
  140          format (/,' KANGLE  --  Too many 3-Ring Angle',
     &                       ' Bending Parameters')
               abort = .true.
            end if
  150       continue
         end if
      end do
c
c     process keywords containing Fourier angle bending parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:7) .eq. 'ANGLEF ') then
            ia = 0
            ib = 0
            ic = 0
            fc = 0.0d0
            an = 0.0d0
            pr = 0.0d0
            string = record(next:120)
            read (string,*,err=160,end=160)  ia,ib,ic,fc,an,pr
  160       continue
            if (header) then
               header = .false.
               write (iout,170)
  170          format (/,' Additional Fourier Angle Bending',
     &                    ' Parameters :',
     &                 //,5x,'Atom Classes',9x,'K(B)',7x,'Shift',
     &                    6x,'Period',/)
            end if
            write (iout,180)  ia,ib,ic,fc,an,pr
  180       format (4x,3i4,2x,3f12.3)
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt = pa//pb//pc
            else
               pt = pc//pb//pa
            end if
            do j = 1, maxnaf
               if (kaf(j).eq.blank .or. kaf(j).eq.pt) then
                  kaf(j) = pt
                  aconf(j) = fc
                  angf(1,j) = an
                  angf(2,j) = pr
                  goto 200
               end if
            end do
            write (iout,190)
  190       format (/,' KANGLE  --  Too many Fourier Angle',
     &                    ' Bending Parameters')
            abort = .true.
  200       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      na = maxna
      na5 = maxna5
      na4 = maxna4
      na3 = maxna3
      naf = maxnaf
      do i = maxna, 1, -1
         if (ka(i) .eq. blank)  na = i - 1
      end do
      do i = maxna5, 1, -1
         if (ka5(i) .eq. blank)  na5 = i - 1
      end do
      do i = maxna4, 1, -1
         if (ka4(i) .eq. blank)  na4 = i - 1
      end do
      do i = maxna3, 1, -1
         if (ka3(i) .eq. blank)  na3 = i - 1
      end do
      do i = maxnaf, 1, -1
         if (kaf(i) .eq. blank)  naf = i - 1
      end do
      use_ring = .false.
      if (min(na5,na4,na3) .ne. 0)  use_ring = .true.
c
c     set generic parameters for use with any number of hydrogens
c
      do i = 1, na
         if (ang(2,i).eq.0.0d0 .and. ang(3,i).eq.0.0d0) then
            ang(2,i) = ang(1,i)
            ang(3,i) = ang(1,i)
         end if
      end do
      do i = 1, na5
         if (ang5(2,i).eq.0.0d0 .and. ang5(3,i).eq.0.0d0) then
            ang5(2,i) = ang5(1,i)
            ang5(3,i) = ang5(1,i)
         end if
      end do
      do i = 1, na4
         if (ang4(2,i).eq.0.0d0 .and. ang4(3,i).eq.0.0d0) then
            ang4(2,i) = ang4(1,i)
            ang4(3,i) = ang4(1,i)
         end if
      end do
      do i = 1, na3
         if (ang3(2,i).eq.0.0d0 .and. ang3(3,i).eq.0.0d0) then
            ang3(2,i) = ang3(1,i)
            ang3(3,i) = ang3(1,i)
         end if
      end do
c
c     assign ideal bond angle and force constant for each angle
c
      header = .true.
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
         ak(i) = 0.0d0
         anat(i) = 0.0d0
         afld(i) = 0.0d0
         angtyp(i) = 'HARMONIC'
         done = .false.
c
c     count number of non-angle hydrogens on the central atom
c
         nh = 1
         do j = 1, n12(ib)
            ih = i12(j,ib)
            if (ih.ne.ia .and. ih.ne.ic .and. atomic(ih).eq.1)
     &         nh = nh + 1
         end do
c
c     make a check for bond angles contained inside small rings
c
         iring = 0
         if (use_ring) then
            call chkring (iring,ia,ib,ic,0)
            if (iring .eq. 6)  iring = 0
            if (iring.eq.5 .and. na5.eq.0)  iring = 0
            if (iring.eq.4 .and. na4.eq.0)  iring = 0
            if (iring.eq.3 .and. na3.eq.0)  iring = 0
         end if
c
c     assign angle bending parameters for bond angles
c
         if (iring .eq. 0) then
            do j = 1, na
               if (ka(j).eq.pt .and. ang(nh,j).ne.0.0d0) then
                  ak(i) = acon(j)
                  anat(i) = ang(nh,j)
                  done = .true.
                  goto 210
               end if
            end do
c
c     assign bending parameters for 5-membered ring angles
c
         else if (iring .eq. 5) then
            do j = 1, na5
               if (ka5(j).eq.pt .and. ang5(nh,j).ne.0.0d0) then
                  ak(i) = acon5(j)
                  anat(i) = ang5(nh,j)
                  done = .true.
                  goto 210
               end if
            end do
c
c     assign bending parameters for 4-membered ring angles
c
         else if (iring .eq. 4) then
            do j = 1, na4
               if (ka4(j).eq.pt .and. ang4(nh,j).ne.0.0d0) then
                  ak(i) = acon4(j)
                  anat(i) = ang4(nh,j)
                  done = .true.
                  goto 210
               end if
            end do
c
c     assign bending parameters for 3-membered ring angles
c
         else if (iring .eq. 3) then
            do j = 1, na3
               if (ka3(j).eq.pt .and. ang3(nh,j).ne.0.0d0) then
                  ak(i) = acon3(j)
                  anat(i) = ang3(nh,j)
                  done = .true.
                  goto 210
               end if
            end do
         end if
c
c     assign Fourier angle bending parameters for bond angles
c
         if (.not. done) then
            do j = 1, naf
               if (kaf(j) .eq. pt) then
                  ak(i) = aconf(j)
                  anat(i) = angf(1,j)
                  afld(i) = angf(2,j)
                  angtyp(i) = 'FOURIER'
                  done = .true.
                  goto 210
               end if
            end do
         end if
c
c     warning if suitable angle bending parameter not found
c
  210    continue
         minat = min(atomic(ia),atomic(ib),atomic(ic))
         if (minat .eq. 0)  done = .true.
         if (use_angle .and. .not.done) then
            if (use(ia) .or. use(ib) .or. use(ic))  abort = .true.
            if (header) then
               header = .false.
               write (iout,220)
  220          format (/,' Undefined Angle Bending Parameters :',
     &                 //,' Type',18x,'Atom Names',19x,
     &                    'Atom Classes',/)
            end if
            label = 'Angle '
            if (iring .eq. 5)  label = '5-Ring'
            if (iring .eq. 4)  label = '4-Ring'
            if (iring .eq. 3)  label = '3-Ring'
            write (iout,230)  label,ia,name(ia),ib,name(ib),
     &                        ic,name(ic),ita,itb,itc
  230       format (1x,a6,5x,3(i6,'-',a3),7x,3i5)
         end if
      end do
c
c     turn off the angle bending potential if it is not used
c
      if (nangle .eq. 0)  use_angle = .false.
      return
      end

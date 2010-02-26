c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine ktors  --  torsional parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ktors" assigns torsional parameters to each torsion in
c     the structure and processes any new or changed values
c
c
      subroutine ktors
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'ktorsn.i'
      include 'math.i'
      include 'potent.i'
      include 'tors.i'
      include 'usage.i'
      integer i,j
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer nt,nt5,nt4
      integer size,next
      integer iring,minat
      integer nlist,ilist
      integer kindex(maxnt)
      integer ft(6)
      real*8 angle,vt(6),st(6)
      logical header,done
      logical use_ring
      character*4 pa,pb,pc,pd
      character*4 zeros
      character*7 label
      character*16 blank
      character*16 pt,pt0
      character*16 pt1,pt2
      character*16 klist(maxnt)
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing torsional angle parameters
c
      blank = '                '
      zeros = '0000'
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         iring = -1
         if (keyword(1:8) .eq. 'TORSION ')  iring = 0
         if (keyword(1:9) .eq. 'TORSION5 ')  iring = 5
         if (keyword(1:9) .eq. 'TORSION4 ')  iring = 4
         if (iring .ge. 0) then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            do j = 1, 6
               vt(j) = 0.0d0
               st(j) = 0.0d0
               ft(j) = 0
            end do
            string = record(next:120)
            read (string,*,err=10,end=10)  ia,ib,ic,id,
     &                                     (vt(j),st(j),ft(j),j=1,6)
   10       continue
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            if (ib .lt. ic) then
               pt = pa//pb//pc//pd
            else if (ic .lt. ib) then
               pt = pd//pc//pb//pa
            else if (ia .le. id) then
               pt = pa//pb//pc//pd
            else if (id .lt. ia) then
               pt = pd//pc//pb//pa
            end if
            call torphase (ft,vt,st)
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Torsional Parameters :',
     &                 //,5x,'Atom Classes',4x,'1-Fold',4x,'2-Fold',
     &                    4x,'3-Fold',4x,'4-Fold',4x,'5-Fold',
     &                    4x,'6-Fold',/)
            end if
            if (iring .eq. 0) then
               write (iout,30)  ia,ib,ic,id,
     &                          (vt(j),nint(st(j)),j=1,6)
   30          format (1x,4i4,1x,6(f6.2,i4))
            else
               if (iring .eq. 5)  label = '5-Ring '
               if (iring .eq. 4)  label = '4-Ring '
               write (iout,40)  ia,ib,ic,id,
     &                          (vt(j),nint(st(j)),j=1,6),label(1:6)
   40          format (1x,4i4,1x,6(f6.2,i4),3x,a6)
            end if
            if (iring .eq. 0) then
               do j = 1, maxnt
                  if (kt(j).eq.blank .or. kt(j).eq.pt) then
                     kt(j) = pt
                     t1(1,j) = vt(1)
                     t1(2,j) = st(1)
                     t2(1,j) = vt(2)
                     t2(2,j) = st(2)
                     t3(1,j) = vt(3)
                     t3(2,j) = st(3)
                     t4(1,j) = vt(4)
                     t4(2,j) = st(4)
                     t5(1,j) = vt(5)
                     t5(2,j) = st(5)
                     t6(1,j) = vt(6)
                     t6(2,j) = st(6)
                     goto 60
                  end if
               end do
               write (iout,50)
   50          format (/,' KTORS  --  Too many Torsional Angle',
     &                    ' Parameters')
               abort = .true.
   60          continue
            else if (iring .eq. 5) then
               do j = 1, maxnt5
                  if (kt5(j).eq.blank .or. kt5(j).eq.pt) then
                     kt5(j) = pt
                     t15(1,j) = vt(1)
                     t15(2,j) = st(1)
                     t25(1,j) = vt(2)
                     t25(2,j) = st(2)
                     t35(1,j) = vt(3)
                     t35(2,j) = st(3)
                     t45(1,j) = vt(4)
                     t45(2,j) = st(4)
                     t55(1,j) = vt(5)
                     t55(2,j) = st(5)
                     t65(1,j) = vt(6)
                     t65(2,j) = st(6)
                     goto 80
                  end if
               end do
               write (iout,70)
   70          format (/,' KTORS  --  Too many 5-Ring Torsional',
     &                    ' Parameters')
               abort = .true.
   80          continue
            else if (iring .eq. 4) then
               do j = 1, maxnt4
                  if (kt4(j).eq.blank .or. kt4(j).eq.pt) then
                     kt4(j) = pt
                     t14(1,j) = vt(1)
                     t14(2,j) = st(1)
                     t24(1,j) = vt(2)
                     t24(2,j) = st(2)
                     t34(1,j) = vt(3)
                     t34(2,j) = st(3)
                     t44(1,j) = vt(4)
                     t44(2,j) = st(4)
                     t54(1,j) = vt(5)
                     t54(2,j) = st(5)
                     t64(1,j) = vt(6)
                     t64(2,j) = st(6)
                     goto 100
                  end if
               end do
               write (iout,90)
   90          format (/,' KTORS  --  Too many 4-Ring Torsional',
     &                    ' Parameters')
               abort = .true.
  100          continue
            end if
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nt = maxnt
      nt5 = maxnt5
      nt4 = maxnt4
      do i = maxnt, 1, -1
         if (kt(i) .eq. blank)  nt = i - 1
      end do
      do i = maxnt5, 1, -1
         if (kt5(i) .eq. blank)  nt5 = i - 1
      end do
      do i = maxnt4, 1, -1
         if (kt4(i) .eq. blank)  nt4 = i - 1
      end do
      use_ring = .false.
      if (min(nt5,nt4) .ne. 0)  use_ring = .true.
c
c     assign torsional parameters for each torsional angle
c     by putting the parameter values into the "tors" arrays
c
      header = .true.
      nlist = 0
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         call numeral (itd,pd,size)
         if (itb .lt. itc) then
            pt = pa//pb//pc//pd
         else if (itc .lt. itb) then
            pt = pd//pc//pb//pa
         else if (ita .le. itd) then
            pt = pa//pb//pc//pd
         else if (itd .lt. ita) then
            pt = pd//pc//pb//pa
         end if
         pt2 = zeros//pt(5:16)
         pt1 = pt(1:12)//zeros
         pt0 = zeros//pt(5:12)//zeros
         tors1(1,i) = 0.0d0
         tors1(2,i) = 0.0d0
         tors2(1,i) = 0.0d0
         tors2(2,i) = 0.0d0
         tors3(1,i) = 0.0d0
         tors3(2,i) = 0.0d0
         tors4(1,i) = 0.0d0
         tors4(2,i) = 0.0d0
         tors5(1,i) = 0.0d0
         tors5(2,i) = 0.0d0
         tors6(1,i) = 0.0d0
         tors6(2,i) = 0.0d0
         done = .false.
c
c     make a check for torsions inside small rings
c
         iring = 0
         if (use_ring) then
            call chkring (iring,ia,ib,ic,id)
            if (iring .eq. 6)  iring = 0
            if (iring.eq.5 .and. nt5.eq.0)  iring = 0
            if (iring.eq.4 .and. nt4.eq.0)  iring = 0
         end if
c
c     find parameters for this torsion; first check "klist"
c     to save time for angle types already located
c
         if (iring .eq. 0) then
            do j = 1, nlist
               if (klist(j) .eq. pt) then
                  ilist = kindex(j)
                  tors1(1,i) = tors1(1,ilist)
                  tors1(2,i) = tors1(2,ilist)
                  tors2(1,i) = tors2(1,ilist)
                  tors2(2,i) = tors2(2,ilist)
                  tors3(1,i) = tors3(1,ilist)
                  tors3(2,i) = tors3(2,ilist)
                  tors4(1,i) = tors4(1,ilist)
                  tors4(2,i) = tors4(2,ilist)
                  tors5(1,i) = tors5(1,ilist)
                  tors5(2,i) = tors5(2,ilist)
                  tors6(1,i) = tors6(1,ilist)
                  tors6(2,i) = tors6(2,ilist)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt
               if (kt(j) .eq. pt) then
                  nlist = nlist + 1
                  klist(nlist) = pt
                  kindex(nlist) = i
                  tors1(1,i) = t1(1,j)
                  tors1(2,i) = t1(2,j)
                  tors2(1,i) = t2(1,j)
                  tors2(2,i) = t2(2,j)
                  tors3(1,i) = t3(1,j)
                  tors3(2,i) = t3(2,j)
                  tors4(1,i) = t4(1,j)
                  tors4(2,i) = t4(2,j)
                  tors5(1,i) = t5(1,j)
                  tors5(2,i) = t5(2,j)
                  tors6(1,i) = t6(1,j)
                  tors6(2,i) = t6(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt
               if (kt(j).eq.pt1 .or. kt(j).eq.pt2) then
                  tors1(1,i) = t1(1,j)
                  tors1(2,i) = t1(2,j)
                  tors2(1,i) = t2(1,j)
                  tors2(2,i) = t2(2,j)
                  tors3(1,i) = t3(1,j)
                  tors3(2,i) = t3(2,j)
                  tors4(1,i) = t4(1,j)
                  tors4(2,i) = t4(2,j)
                  tors5(1,i) = t5(1,j)
                  tors5(2,i) = t5(2,j)
                  tors6(1,i) = t6(1,j)
                  tors6(2,i) = t6(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt
               if (kt(j) .eq. pt0) then
                  tors1(1,i) = t1(1,j)
                  tors1(2,i) = t1(2,j)
                  tors2(1,i) = t2(1,j)
                  tors2(2,i) = t2(2,j)
                  tors3(1,i) = t3(1,j)
                  tors3(2,i) = t3(2,j)
                  tors4(1,i) = t4(1,j)
                  tors4(2,i) = t4(2,j)
                  tors5(1,i) = t5(1,j)
                  tors5(2,i) = t5(2,j)
                  tors6(1,i) = t6(1,j)
                  tors6(2,i) = t6(2,j)
                  done = .true.
                  goto 110
               end if
            end do
c
c     find the parameters for a 5-ring torsion
c
         else if (iring .eq. 5) then
            do j = 1, nt5
               if (kt5(j) .eq. pt) then
                  tors1(1,i) = t15(1,j)
                  tors1(2,i) = t15(2,j)
                  tors2(1,i) = t25(1,j)
                  tors2(2,i) = t25(2,j)
                  tors3(1,i) = t35(1,j)
                  tors3(2,i) = t35(2,j)
                  tors4(1,i) = t45(1,j)
                  tors4(2,i) = t45(2,j)
                  tors5(1,i) = t55(1,j)
                  tors5(2,i) = t55(2,j)
                  tors6(1,i) = t65(1,j)
                  tors6(2,i) = t65(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt5
               if (kt5(j).eq.pt1 .or. kt5(j).eq.pt2) then
                  tors1(1,i) = t15(1,j)
                  tors1(2,i) = t15(2,j)
                  tors2(1,i) = t25(1,j)
                  tors2(2,i) = t25(2,j)
                  tors3(1,i) = t35(1,j)
                  tors3(2,i) = t35(2,j)
                  tors4(1,i) = t45(1,j)
                  tors4(2,i) = t45(2,j)
                  tors5(1,i) = t55(1,j)
                  tors5(2,i) = t55(2,j)
                  tors6(1,i) = t65(1,j)
                  tors6(2,i) = t65(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt5
               if (kt5(j) .eq. pt0) then
                  tors1(1,i) = t15(1,j)
                  tors1(2,i) = t15(2,j)
                  tors2(1,i) = t25(1,j)
                  tors2(2,i) = t25(2,j)
                  tors3(1,i) = t35(1,j)
                  tors3(2,i) = t35(2,j)
                  tors4(1,i) = t45(1,j)
                  tors4(2,i) = t45(2,j)
                  tors5(1,i) = t55(1,j)
                  tors5(2,i) = t55(2,j)
                  tors6(1,i) = t65(1,j)
                  tors6(2,i) = t65(2,j)
                  done = .true.
                  goto 110
               end if
            end do
c
c     find the parameters for a 4-ring torsion
c
         else if (iring .eq. 4) then
            do j = 1, nt4
               if (kt4(j) .eq. pt) then
                  tors1(1,i) = t14(1,j)
                  tors1(2,i) = t14(2,j)
                  tors2(1,i) = t24(1,j)
                  tors2(2,i) = t24(2,j)
                  tors3(1,i) = t34(1,j)
                  tors3(2,i) = t34(2,j)
                  tors4(1,i) = t44(1,j)
                  tors4(2,i) = t44(2,j)
                  tors5(1,i) = t54(1,j)
                  tors5(2,i) = t54(2,j)
                  tors6(1,i) = t64(1,j)
                  tors6(2,i) = t64(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt4
               if (kt4(j).eq.pt1 .or. kt4(j).eq.pt2) then
                  tors1(1,i) = t14(1,j)
                  tors1(2,i) = t14(2,j)
                  tors2(1,i) = t24(1,j)
                  tors2(2,i) = t24(2,j)
                  tors3(1,i) = t34(1,j)
                  tors3(2,i) = t34(2,j)
                  tors4(1,i) = t44(1,j)
                  tors4(2,i) = t44(2,j)
                  tors5(1,i) = t54(1,j)
                  tors5(2,i) = t54(2,j)
                  tors6(1,i) = t64(1,j)
                  tors6(2,i) = t64(2,j)
                  done = .true.
                  goto 110
               end if
            end do
            do j = 1, nt4
               if (kt4(j) .eq. pt0) then
                  tors1(1,i) = t14(1,j)
                  tors1(2,i) = t14(2,j)
                  tors2(1,i) = t24(1,j)
                  tors2(2,i) = t24(2,j)
                  tors3(1,i) = t34(1,j)
                  tors3(2,i) = t34(2,j)
                  tors4(1,i) = t44(1,j)
                  tors4(2,i) = t44(2,j)
                  tors5(1,i) = t54(1,j)
                  tors5(2,i) = t54(2,j)
                  tors6(1,i) = t64(1,j)
                  tors6(2,i) = t64(2,j)
                  done = .true.
                  goto 110
               end if
            end do
         end if
c
c     warning if suitable torsional parameter not found
c
  110    continue
         minat = min(atomic(ia),atomic(ib),atomic(ic),atomic(id))
         if (minat .eq. 0)  done = .true.
         if (.not. done) then
            if (use_tors) then
               if (use(ia) .or. use(ib) .or. use(ic) .or. use(id))
     &            abort = .true.
            end if
            if (header) then
               header = .false.
               write (iout,120)
  120          format (/,' Undefined Torsional Parameters :',
     &                 //,' Type',24x,'Atom Names',24x,
     &                    'Atom Classes',/)
            end if
            label = 'Torsion'
            if (iring .eq. 5)  label = '5-Ring '
            if (iring .eq. 4)  label = '4-Ring '
            write (iout,130)  label,ia,name(ia),ib,name(ib),ic,
     &                        name(ic),id,name(id),ita,itb,itc,itd
  130       format (1x,a7,4x,4(i6,'-',a3),5x,4i5)
         end if
      end do
c
c     find the cosine and sine of phase angle for each torsion
c
      do i = 1, ntors
         angle = tors1(2,i) / radian
         tors1(3,i) = cos(angle)
         tors1(4,i) = sin(angle)
         angle = tors2(2,i) / radian
         tors2(3,i) = cos(angle)
         tors2(4,i) = sin(angle)
         angle = tors3(2,i) / radian
         tors3(3,i) = cos(angle)
         tors3(4,i) = sin(angle)
         angle = tors4(2,i) / radian
         tors4(3,i) = cos(angle)
         tors4(4,i) = sin(angle)
         angle = tors5(2,i) / radian
         tors5(3,i) = cos(angle)
         tors5(4,i) = sin(angle)
         angle = tors6(2,i) / radian
         tors6(3,i) = cos(angle)
         tors6(4,i) = sin(angle)
      end do
c
c     turn off the torsional potential if it is not used
c
      if (ntors .eq. 0)  use_tors = .false.
      return
      end

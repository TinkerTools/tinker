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
      include 'ring.i'
      include 'usage.i'
cn
      include 'bond.i'
      include 'fields.i'
      include 'merck.i'
cn
      integer i,j,k
      integer ia,ib,ic,id
      integer ita,itb,itc
      integer na,na5,na4
      integer na3,naf
      integer jen,ih,nh
      integer next,size
      integer minat,iring
cn
      integer AT,collclass
      integer itta,ittb,ittc
      integer bnd_ab,bnd_bc
      integer l,m
      real*8 fc,an,pr
      real*8 an1,an2,an3
      real*8 Z2(53),C(53),D,beta
cn      
      logical header,done
      logical use_ring3
      logical use_ring4
      logical use_ring5
cn
      logical ring3,ring4
cn      
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
c     use small rings if present and parameters are available
c
      use_ring5 = .true.
      use_ring4 = .true.
      use_ring3 = .true.
      if (nring5.eq.0 .or. na5.eq.0)  use_ring5 = .false.
      if (nring4.eq.0 .or. na4.eq.0)  use_ring4 = .false.
      if (nring3.eq.0 .or. na3.eq.0)  use_ring3 = .false.
c
c     assign ideal bond angle and force constant for each angle
c
      header = .true.

cn

      if (forcefield .eq. 'MMFF94') then
cn
cn fills in the tables with the empirical rule parameters
cn (cf T. A. Halgren, "Merck Molecular Force Field. V.
cn Extension of MMFF94 Using Experimental Data,
cn Additional Computational Data, and Empirical Rules",
cn J. Comput. Chem. 17, pp. 616-641)
cn

       do i = 1, 53
cn if parameter doesn't exist, it will give "1000",
cn we thus can see that it doesn't exist.
        Z2(i) = 1000.0
        C(i) = 1000.0
       enddo

       Z2(1) = 1.395
       Z2(5) = 0.0
       Z2(6) = 2.494
       Z2(7) = 2.711
       Z2(8) = 3.045
       Z2(9) = 2.847
       Z2(14) = 2.350
       Z2(15) = 2.350
       Z2(16) = 2.980
       Z2(17) = 2.909
       Z2(35) = 3.017
       Z2(33) = 0.0
       Z2(53) = 3.086
       C(1) = 0.0
       C(5) = 0.704
       C(6) = 1.016
       C(7) = 1.113
       C(8) = 1.337
       C(9) = 0.0
       C(14) = 0.811
       C(15) = 1.068
       C(16) = 1.249
       C(17) = 1.078
       C(35) = 0.0
       C(33) = 0.825
       C(53) = 0.0


      do i = 1, nangle       
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itta = type(ia)
         ittb = type(ib)
         ittc = type(ic) 

cn
cn calculation of the the AT
cn

       AT = 0
       
cn looks at the BT

       do j=1,nlignes
        if (((ia .eq. BT_1(j,1)).AND.(ib.eq.BT_1(j,2)))
     &     .OR.((ib.eq.BT_1(j,1)).AND.(ia.eq.BT_1(j,2)))) then
        AT = AT+1
        endif
        if (((ic .eq. BT_1(j,1)).AND.(ib.eq.BT_1(j,2)))
     &     .OR.((ib.eq.BT_1(j,1)).AND.(ic.eq.BT_1(j,2)))) then
        AT = AT+1
        endif
       enddo
cn looks wether the atoms belong to a 3- or 4-membered ring

       ring3 = .false.
       ring4 = .false.
       do j = 1,nring3       
        do k = 1,3
         if (ia.eq.iring3(k,j)) then
          do l = 1,3
           if (ib.eq.iring3(l,j)) then
            do m = 1,3           
             if (ic.eq.iring3(m,j)) ring3 = .true.
            enddo
           endif
          enddo
         endif
        enddo
       enddo
       if (ring3.neqv..true.) then
       do j = 1,nring4
        do k = 1,4
         if (ia.eq.iring4(k,j)) then
          do l = 1,4
           if (ib.eq.iring4(l,j)) then
            do m = 1,4
             if (ic.eq.iring4(m,j)) ring4 = .true. 
            enddo
           endif
          enddo
         endif
        enddo
       enddo
       endif
       if (AT.eq.0.AND.ring4.eqv..true.) then
cn sum of BT = 0 and 4-membered ring     
        AT = 4
       else if (AT.eq.1.AND.ring4.eqv..true.) then
cn sum of BT = 1 and 4-membered ring     
        AT = 7
       else if (AT.eq.2.AND.ring4.eqv..true.) then
cn sum of BT = 2 and 4-membered ring
        AT = 8
       else if (AT.eq.0.AND.ring3.eqv..true.) then
cn sum of BT = 0 and 3-membered ring
        AT = 3
       else if (AT.eq.1.AND.ring3.eqv..true.) then
cn sum of BT = 1 and 3-membered ring
        AT = 5
       else if (AT.eq.2.AND.ring3.eqv..true.) then
cn sum of BT = 2 and 3-membered ring     
        AT = 6
       endif

cn
cn atom class equivalencies assignment
cn (cf T. A. Halgren, "Merck Molecular Force Field. I. 
cn Basis, Form, Scope, Parametrization, and Performance
cn of MMFF94", J. Comput. Chem. 17, pp. 490-519)
cn


       collclass = 0
  201  continue
       collclass = collclass+1
       if (collclass.eq.1) then
        ita = eqclass(itta,1)
        itb = eqclass(ittb,1)
        itc = eqclass(ittc,1)
       else if (collclass.eq.2) then
        ita = eqclass(itta,2)
        itb = eqclass(ittb,2)
        itc = eqclass(ittc,2)
       else if (collclass.eq.3) then
        ita = eqclass(itta,3)
        itb = eqclass(ittb,2)
        itc = eqclass(ittc,3)
       else if (collclass.eq.4) then
        ita = eqclass(itta,4)
        itb = eqclass(ittb,2)
        itc = eqclass(ittc,4)
       else if (collclass.eq.5) then
        ita = eqclass(itta,5)
        itb = eqclass(ittb,2)
        itc = eqclass(ittc,5)
       endif

       if (collclass .gt. 5) then
         goto 221
       
       else if (AT .eq. 0) then
        ak(i) = MMFF_ka(ita,itb,itc)
        anat(i) = MMFF_teta0(ita,itb,itc)
        if (collclass.eq.5) then
  202   continue
cn        
cn if "wild card" parameters : still need to calculate
cn the force constant by an empirical rule.
cn if (full) empirical rule : the reference angle
cn has already been calculated at mark 221.
cn
         if (Z2(atomic(ia)).eq.1000) goto 221
         if (Z2(atomic(ib)).eq.1000) goto 221
         if (Z2(atomic(ic)).eq.1000) goto 221
         if (C(atomic(ia)).eq.1000) goto 221
         if (C(atomic(ib)).eq.1000) goto 221
         if (C(atomic(ic)).eq.1000) goto 221
         do k = 1, maxbnd
cn take the number of the two bonds         
          if ((min(ia,ib).eq.ibnd(1,k))
     &        .AND.(max(ia,ib).eq.ibnd(2,k))) then
              bnd_ab = k
          endif
          if ((min(ic,ib).eq.ibnd(1,k))
     &       .AND.(max(ic,ib).eq.ibnd(2,k))) then
             bnd_bc = k
          endif
         enddo
         D = (bl(bnd_ab)-bl(bnd_bc))**2/
     &       (bl(bnd_ab)+bl(bnd_bc))**2
     
         beta = 1.0
         if (ring4.eqv..true.) beta = 0.85
         if (ring3.eqv..true.) beta = 0.05
         ak(i) = beta*1.75*Z2(atomic(ia))*Z2(atomic(ic))
     &           *C(atomic(ib))
     &           /((anat(i)*0.01745329252)**2
     &           *(bl(bnd_ab)+bl(bnd_bc))*exp(2*D))
        endif
                  done = .true.

               if (ak(i) .eq. 1000.0) done = .false.
               if (anat(i) .eq. 1000.0) done = .false.
         if (done.eqv..false.) then
          goto 201
         endif
                  goto 221
       else if (AT .eq. 1) then
       ak(i) = MMFF_ka1(ita,itb,itc)
       anat(i) = MMFF_teta01(ita,itb,itc)
                  done = .true.

               if (ak(i) .eq. 1000.0) done = .false.
               if (anat(i) .eq. 1000.0) done = .false.
         if (done.eqv..false..AND.collclass.lt.5) then
          goto 201
         else if (collclass.eq.5) then
          goto 202
         endif

                  goto 221
                  
       else if (AT .eq. 2) then
       ak(i) = MMFF_ka2(ita,itb,itc)
       anat(i) = MMFF_teta02(ita,itb,itc)
                  done = .true.
                  
               if (ak(i) .eq. 1000.0) done = .false.
               if (anat(i) .eq. 1000.0) done = .false.
         if (done.eqv..false..AND.collclass.lt.5) then
          goto 201
         else if (collclass.eq.5) then
          goto 202
         endif    

                  goto 221

       else if (AT .eq. 3) then
       ak(i) = MMFF_ka3(ita,itb,itc)
       anat(i) = MMFF_teta03(ita,itb,itc)
                  done = .true.

               if (ak(i) .eq. 1000.0) done = .false.
               if (anat(i) .eq. 1000.0) done = .false.
         if (done.eqv..false..AND.collclass.lt.5) then
          goto 201
         else if (collclass.eq.5) then
          goto 202
         endif

                  goto 221

       else if (AT .eq. 4) then
       ak(i) = MMFF_ka4(ita,itb,itc)
       anat(i) = MMFF_teta04(ita,itb,itc)

                  done = .true.

               if (ak(i) .eq. 1000.0) done = .false.
               if (anat(i) .eq. 1000.0) done = .false.
         if (done.eqv..false..AND.collclass.lt.5) then
          goto 201
         else if (collclass.eq.5) then
          goto 202
         endif

                  goto 221
       else if (AT .eq. 5) then
       ak(i) = MMFF_ka5(ita,itb,itc)
       anat(i) = MMFF_teta05(ita,itb,itc)

                  done = .true.

               if (ak(i) .eq. 1000.0) done = .false.
               if (anat(i) .eq. 1000.0) done = .false.
         if (done.eqv..false..AND.collclass.lt.5) then
          goto 201
         else if (collclass.eq.5) then
          goto 202
         endif

                  goto 221
       else if (AT .eq. 6) then
       ak(i) = MMFF_ka6(ita,itb,itc)
       anat(i) = MMFF_teta06(ita,itb,itc)

                  done = .true.

               if (ak(i) .eq. 1000.0) done = .false.
               if (anat(i) .eq. 1000.0) done = .false.
         if (done.eqv..false..AND.collclass.lt.5) then
          goto 201
         else if (collclass.eq.5) then
          goto 202
         endif

                  goto 221
       else if (AT .eq. 7) then
       ak(i) = MMFF_ka7(ita,itb,itc)
       anat(i) = MMFF_teta07(ita,itb,itc)

                  done = .true.

               if (ak(i) .eq. 1000.0) done = .false.
               if (anat(i) .eq. 1000.0) done = .false.
         if (done.eqv..false..AND.collclass.lt.5) then
          goto 201
         else if (collclass.eq.5) then         
          goto 202
         endif

                  goto 221
       else if (AT .eq. 8) then
       ak(i) = MMFF_ka8(ita,itb,itc)
       anat(i) = MMFF_teta08(ita,itb,itc)

                  done = .true.

               if (ak(i) .eq. 1000.0) done = .false.
               if (anat(i) .eq. 1000.0) done = .false.
         if (done.eqv..false..AND.collclass.lt.5) then
          goto 201
         else if (collclass.eq.5) then
          goto 202
         endif

                  goto 221
       endif
c
c     warning if suitable angle bending parameter not found
c
  221    continue

         minat = min(atomic(ia),atomic(ib),atomic(ic))
         if (minat .eq. 0)  done = .true.
         if (.not. done) then
            if (use_angle) then
cn when angle bending parameters are missing,
cn the empirical rules applies
             anat(i) = 120
             if (crd(itb).eq.4) anat(i) = 109.45
             if (crd(itb).eq.2) then
              if (atomic(ib).eq.8) then
               anat(i) = 105.0
              else if (atomic(ib).gt.10) then
               anat(i) = 95.0
              else if (lin(itb).eq.1) then
               anat(i) = 180.0
              endif
             endif
             if (crd(itb).eq.3.AND.val(itb).eq.3
     &                .AND.mltb(itb).eq.0) then
              if (atomic(ib).eq.7) then
               anat(i) = 107.0
              else
               anat(i) = 92.0
              endif
             endif
             if (ring3.eqv..true.) anat(i) = 60.0
             if (ring4.eqv..true.) anat(i) = 90.0
             goto 202
            end if
         end if

         angtyp(i) = 'HARMONIC'
         if (anat(i) .eq. 180.0d0)  angtyp(i) = 'LINEAR'

      end do

      else

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
c     make a check for bond angles inside small rings
c
         iring = 0
         if (n12(ia).eq.1 .or. n12(ic).eq.1)  goto 210
         if (use_ring3) then
            do j = 1, n12(ia)
               if (ic .eq. i12(j,ia)) then
                  iring = 3
                  goto 210
               end if
            end do
         end if
         if (use_ring4) then
            do j = 1, n12(ia)
               id = i12(j,ia)
               if (ib .ne. id) then
                  do k = 1, n12(ic)
                     if (id .eq. i12(k,ic)) then
                        iring = 4
                        goto 210
                     end if
                  end do
               end if
            end do
         end if
         if (use_ring5) then
            do j = 1, n12(ia)
               id = i12(j,ia)
               if (ib .ne. id) then
                  do k = 1, n13(ic)
                     if (id .eq. i13(k,ic)) then
                        iring = 5
                        goto 210
                     end if
                  end do
               end if
            end do
         end if
  210    continue
c
c     assign angle bending parameters for bond angles
c
         if (iring .eq. 0) then
            do j = 1, na
               if (ka(j).eq.pt .and. ang(nh,j).ne.0.0d0) then
                  ak(i) = acon(j)
                  anat(i) = ang(nh,j)
                  done = .true.
                  goto 220
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
                  goto 220
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
                  goto 220
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
                  goto 220
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
                  goto 220
               end if
            end do
         end if
c
c     warning if suitable angle bending parameter not found
c
  220    continue
         minat = min(atomic(ia),atomic(ib),atomic(ic))
         if (minat .eq. 0)  done = .true.
         if (.not. done) then
            if (use_angle) then
               if (use(ia) .or. use(ib) .or. use(ic))  abort = .true.
            end if
            if (header) then
               header = .false.
               write (iout,230)
  230          format (/,' Undefined Angle Bending Parameters :',
     &                 //,' Type',18x,'Atom Names',19x,
     &                    'Atom Classes',/)
            end if
            label = 'Angle '
            if (iring .eq. 5)  label = '5-Ring'
            if (iring .eq. 4)  label = '4-Ring'
            if (iring .eq. 3)  label = '3-Ring'
            write (iout,240)  label,ia,name(ia),ib,name(ib),
     &                        ic,name(ic),ita,itb,itc
  240       format (1x,a6,5x,3(i6,'-',a3),7x,3i5)
         end if
      end do

       endif

cn
c
c     turn off the angle bending potential if it is not used
c
      if (nangle .eq. 0)  use_angle = .false.

      return
      end

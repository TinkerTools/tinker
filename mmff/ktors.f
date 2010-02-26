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
      include 'fields.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'ktorsn.i'
      include 'math.i'
      include 'merck.i'
      include 'potent.i'
      include 'ring.i'
      include 'tors.i'
      include 'usage.i'
      integer i,j,k
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer nt,nt5,nt4
      integer size,next
      integer iring,minat
      integer nlist,ilist
      integer kindex(maxnt)
      integer ft(6)
cn
      integer AB,BC,CD
      integer itta,ittb,ittc,ittd
      integer collclass
      integer TT,l,m,o
      real*8 beta,pi_bc,Ub,Uc
      real*8 Vb,Vc,N_bc
      real*8 W_b,W_c
      integer rowtab_b,rowtab_c 
      logical ring4,ring5
      logical skipring
cn
      real*8 angle,vt(6),st(6)
      logical header,done
      logical use_ring5
      logical use_ring4
      character*4 pa,pb,pc,pd
      character*7 label
      character*16 blank,pt
      character*16 klist(maxnt)
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing torsional angle parameters
c
      blank = '                '
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
c
c     use small rings if present and parameters are available
c
      use_ring5 = .true.
      use_ring4 = .true.
      if (nring5.eq.0 .or. nt5.eq.0)  use_ring5 = .false.
      if (nring4.eq.0 .or. nt4.eq.0)  use_ring4 = .false.
c
c     assign torsional parameters for each torsional angle
c     by putting the parameter values into the "tors" arrays
c
      header = .true.
      nlist = 0

cn

      if (forcefield .eq. 'MMFF94') then      


      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         itta = type(ia)
         ittb = type(ib)
         ittc = type(ic)
         ittd = type(id)
         done = .false.
         collclass = 0
         skipring = .false.
  101    continue
         
cn
cn atom class equivalencies assignment
cn (cf T. A. Halgren, "Merck Molecular Force Field. I.
cn Basis, Form, Scope, Parametrization, and Performance
cn of MMFF94", J. Comput. Chem. 17, pp. 490-519)
cn
         
         collclass = collclass+1 
          if (collclass.eq.1) then
           ita = eqclass(itta,collclass)
           itb = eqclass(ittb,collclass)
           itc = eqclass(ittc,collclass)
           itd = eqclass(ittd,collclass)
          else if (collclass.eq.2) then
           ita = eqclass(itta,collclass)
           itb = eqclass(ittb,collclass)
           itc = eqclass(ittc,collclass)
           itd = eqclass(ittd,collclass)
          else if (collclass.eq.3) then
           ita = eqclass(itta,3)
           itb = eqclass(ittb,2)
           itc = eqclass(ittc,2)
           itd = eqclass(ittd,5)
          else if (collclass.eq.4) then
           ita = eqclass(itta,5)
           itb = eqclass(ittb,2)
           itc = eqclass(ittc,2)
           itd = eqclass(ittd,3)
          else if (collclass.eq.5) then
           ita = eqclass(itta,5)
           itb = eqclass(ittb,2)
           itc = eqclass(ittc,2)
           itd = eqclass(ittd,5)
          endif

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

cn
cn Torsion Type attribution
cn depending on the BT of each bond and the belonging
cn to a 4- or 5-membered ring
cn

cn looks at the BT
         AB=0
         if (ia .le. ib) then
          do j = 1 ,nlignes
           if ((ia .eq. BT_1(j,1)) .AND. (ib .eq. BT_1(j,2))) then
            AB=1
           endif
          enddo
         else if (ib .le. ia) then
          do j = 1 ,nlignes
           if ((ib .eq. BT_1(j,1)) .AND. (ia .eq. BT_1(j,2))) then
            AB=1
           endif
          enddo
         endif

         BC=0
         if (ib .le. ic) then
          do j = 1,nlignes
           if ((ib .eq. BT_1(j,1)) .AND. (ic .eq. BT_1(j,2))) then
            BC=1
           endif
          enddo
         else if (ic .le. ib) then
          do j = 1,nlignes
           if ((ic .eq. BT_1(j,1)) .AND. (ib .eq. BT_1(j,2))) then
            BC=1
           endif
          enddo
         endif

         CD=0
         if (ic .le. id) then
          do j = 1,nlignes
           if ((ic .eq. BT_1(j,1)) .AND. (id .eq. BT_1(j,2))) then
            CD=1
           endif
          enddo
         else if (id .le. ic) then
          do j = 1,nlignes
           if ((id .eq. BT_1(j,1)) .AND. (ic .eq. BT_1(j,2))) then
            CD=1
           endif
          enddo
         endif

c
c     make a check for torsions inside small rings
c
       ring4 = .false.
       ring5 = .false.
       
       do j = 1,nring4
        do k = 1,4
         if (ia.eq.iring4(k,j)) then
          do l = 1,4
           if (ib.eq.iring4(l,j)) then
            do m = 1,4
             if (ic.eq.iring4(m,j)) then
              do o = 1,4
               if(id.eq.iring4(o,j)) ring4 = .true.
              enddo
             endif
            enddo
           endif
          enddo
         endif
        enddo
       enddo
       do j = 1,nring5
        do k = 1,5
         if (ia.eq.iring5(k,j)) then
          do l = 1,5
           if (ib.eq.iring5(l,j)) then
            do m = 1,5
             if (ic.eq.iring5(m,j)) then
              do o = 1,5
               if (id.eq.iring5(o,j)) ring5 = .true.
              enddo
             endif
            enddo
           endif
          enddo
         endif
        enddo
       enddo

       if (skipring.eq..true.) then
        ring4 = .false.
        ring5 = .false.
       endif

cn looks at the small cycles

        if (ring4.eqv..true.) then
cn 4-memebred cycle
          TT = 4
            do j = 1, nt4
               if (kt4(j) .eq. pt) then
                  tors1(1,i) = t14(1,j)
                  tors1(2,i) = t14(2,j)
                  tors2(1,i) = t24(1,j)
                  tors2(2,i) = t24(2,j)
                  tors3(1,i) = t34(1,j)
                  tors3(2,i) = t34(2,j)
                  done = .true.
                  goto 123
               end if
            end do

         if ((done.eqv..false.).and.collclass.lt.5) then
           goto 101
         endif
        endif

        if ((ring5.eqv..true.).AND.
     &       (class(ia).eq.1.OR.class(ib).eq.1
     &   .OR.class(ic).eq.1.OR.class(id).eq.1)) then
cn 5-memebered cycle with at least one of the atoms of class 1
cn (sp3-hybridized atom)
cn
          TT = 5
          do j = 1, nt5
           if (kt5(j) .eq. pt) then
            tors1(1,i) = t15(1,j)
            tors1(2,i) = t15(2,j)
            tors2(1,i) = t25(1,j)
            tors2(2,i) = t25(2,j)
            tors3(1,i) = t35(1,j)
            tors3(2,i) = t35(2,j)
            done = .true.
           end if
          end do

         if ((done.eqv..false.).and.collclass.lt.5) then
           goto 101
         else if ((done.eqv..false.).and.collclass.eq.5) then
          collclass = 0
          skipring = .true.
          goto 101
         endif
        endif

cn if AB and/or CD =1, with BC=0, then TT = 2


         if ((AB .eq. 1.AND.(mltb(class(ic)).eq.0.OR.
     &      sbmb(class(ic)).eq.0))
     &      .OR.(CD .eq. 1.AND.(mltb(class(ib)).eq.0.OR.
     &      sbmb(class(ib)).eq.0))) then
         
cn this condition has not been written that way in MMFF94's
cn litterature, but this has been deduced when validating
cn the implementation of the force field, using
cn the validation suite :
cn http://server.ccl.net/cca/data/MMFF94/index.shtml

            TT = 2
            do j = 1, maxnt
               if (kt_2(j) .eq. pt) then
                  tors1(1,i) = t1_2(1,j)
                  tors1(2,i) = t1_2(2,j)
                  tors2(1,i) = t2_2(1,j)
                  tors2(2,i) = t2_2(2,j)
                  tors3(1,i) = t3_2(1,j)
                  tors3(2,i) = t3_2(2,j)
                  done = .true.
                  goto 123
               end if
            end do

         if ((done.eqv..false.).and.collclass.lt.5) then
           goto 101
         endif
         if ((done.eqv..false.).and.collclass.eq.5) then
           goto 111
cn => "fully wildcarded" parameters with TT=0 instead of TT=2
         endif


                if (tors1(1,i) .eq. 1000.0) done = .false.
                if (tors1(2,i) .eq. 1000.0) done = .false.
                if (tors2(1,i) .eq. 1000.0) done = .false.
                if (tors2(2,i) .eq. 1000.0) done = .false.
                if (tors3(1,i) .eq. 1000.0) done = .false.
                if (tors3(2,i) .eq. 1000.0) done = .false.

                  goto 123
cn if BC=1, then TT =1
         else if (BC .eq. 1) then
            TT = 1
            do j = 1, maxnt          
               if (kt_1(j) .eq. pt) then
                  tors1(1,i) = t1_1(1,j)
                  tors1(2,i) = t1_1(2,j)
                  tors2(1,i) = t2_1(1,j)
                  tors2(2,i) = t2_1(2,j)
                  tors3(1,i) = t3_1(1,j)
                  tors3(2,i) = t3_1(2,j)
                  done = .true.
                  goto 123

               end if
            end do


         if ((done.eqv..false.).and.collclass.lt.5) then
           goto 101
         endif

                if (tors1(1,i) .eq. 1000.0) done = .false.
                if (tors1(2,i) .eq. 1000.0) done = .false.
                if (tors2(1,i) .eq. 1000.0) done = .false.
                if (tors2(2,i) .eq. 1000.0) done = .false.
                if (tors3(1,i) .eq. 1000.0) done = .false.
                if (tors3(2,i) .eq. 1000.0) done = .false.

                  goto 123
cn in all other cases, TT = 0 
         else if (done.neqv..true.) then
cn
cn !!! pour TT4/5 !!!
cn
  111       continue 
            TT = 0
            do j = 1, maxnt
            
               if (kt(j) .eq. pt) then
                  tors1(1,i) = t1(1,j)
                  tors1(2,i) = t1(2,j)
                  tors2(1,i) = t2(1,j)
                  tors2(2,i) = t2(2,j)
                  tors3(1,i) = t3(1,j)
                  tors3(2,i) = t3(2,j)
                  done = .true.
                  goto 123
               end if
            end do

         if ((done.eqv..false.).and.collclass.lt.5) then
           goto 101
         endif

                if (tors1(1,i) .eq. 1000.0) done = .false.
                if (tors1(2,i) .eq. 1000.0) done = .false.
                if (tors2(1,i) .eq. 1000.0) done = .false.
                if (tors2(2,i) .eq. 1000.0) done = .false.
                if (tors3(1,i) .eq. 1000.0) done = .false.
                if (tors3(2,i) .eq. 1000.0) done = .false.

                  goto 123
         endif


c
c     warning if suitable torsional parameter not found
c
  123    continue


cn
cn if the parameters aren't available,
cn we use the "empirical rule".
cn (cf T. A. Halgren, "Merck Molecular Force Field. V.
cn Extension of MMFF94 Using Experimental Data,
cn Additional Computational Data, and Empirical Rules",
cn J. Comput. Chem. 17, pp. 616-641)
cn

         if (.not. done) then

          ita=class(ia)
          itb=class(ib)
          itc=class(ic)
          itd=class(id)
            if (atomic(ib).eq.6) then
             Ub=2.0
             Vb=2.12
            else if (atomic(ib).eq.7) then
             Ub=2.0
             Vb=1.5
            else if (atomic(ib).eq.8) then
             Ub=2.0
             Vb=0.2
            else if (atomic(ib).eq.14) then
             Ub=1.25
             Vb=1.22
            else if (atomic(ib).eq.15) then
             Ub=1.25
             Vb=2.4
            else if (atomic(ib).eq.16) then
             Ub=1.25
             Vb=0.49
            endif
            if (atomic(ic).eq.6) then
             Uc=2.0
             Vc=2.12
            else if (atomic(ic).eq.7) then
             Uc=2.0
             Vc=1.5
            else if (atomic(ic).eq.8) then
             Uc=2.0
             Vc=0.2
            else if (atomic(ic).eq.14) then
             Uc=1.25
             Vc=1.22
            else if (atomic(ic).eq.15) then
             Uc=1.25
             Vc=2.4
            else if (atomic(ic).eq.16) then
             Uc=1.25
             Vc=0.49
            endif

              N_bc=(crd(itb)-1)*(crd(itc)-1)

          if(atomic(ib).eq.1) rowtab_b=0
          if((atomic(ib).ge.3).and.(atomic(ib).le.10)) rowtab_b=1
          if((atomic(ib).ge.11).and.(atomic(ib).le.18)) rowtab_b=2
          if((atomic(ib).ge.19).and.(atomic(ib).le.36)) rowtab_b=3
          if((atomic(ib).ge.37).and.(atomic(ib).le.54)) rowtab_b=4
          if(atomic(ic).eq.1) rowtab_c=0
          if((atomic(ic).ge.3).and.(atomic(ic).le.10)) rowtab_c=1
          if((atomic(ic).ge.11).and.(atomic(ic).le.18)) rowtab_c=2
          if((atomic(ic).ge.19).and.(atomic(ic).le.36)) rowtab_c=3
          if((atomic(ic).ge.37).and.(atomic(ic).le.54)) rowtab_c=4

          if (lin(itb).eq.1.OR.lin(itc).eq.1) then
cn page 631, in "Merck Molecular Force Field V", a)
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) = 0.0
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0
                  done = .true.
                  goto 123
          else if (arom(itb).eq.1.and.arom(itc).eq.1) then
cn page 631, in "Merck Molecular Force Field V", b) 
            if(pilp(itb).eq.0.and.pilp(itc).eq.0) then
              pi_bc = 0.5
            else
              pi_bc = 0.3
            endif
            if((val(itb).eq.3.and.val(itc).eq.4).or.
     &         (val(itb).eq.4.and.val(itc).eq.3).or.     
     &         (val(itb).eq.4.and.val(itc).eq.34).or.
     &         (val(itb).eq.34.and.val(itc).eq.4).or.
     &         (val(itb).eq.34.and.val(itc).eq.3).or.
     &         (val(itb).eq.3.and.val(itc).eq.34).or.
     &         (val(itb).eq.34.and.val(itc).eq.34)) then
              beta=3
            else
              beta=6
            endif
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) =beta*pi_bc*sqrt(Ub*Uc) 
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if ((mltb(itb).eq.2.and.mltb(itc).eq.2).or.
     &            (mltb(itc).eq.2.and.mltb(itb).eq.2)) then
cn page 631, in "Merck Molecular Force Field V", c)
                   beta=6
                   pi_bc=1
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) =beta*pi_bc*sqrt(Ub*Uc)
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if (mltb(itb).eq.2.or.mltb(itc).eq.2) then
              beta=6
              pi_bc=0.4
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) =beta*pi_bc*sqrt(Ub*Uc)
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if (crd(itb).eq.4.and.crd(itc).eq.4) then
cn page 632, in "Merck Molecular Force Field V", d) 
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) = 0.0
                  tors2(2,i) = 180.0
                  tors3(1,i) = sqrt(Vb*Vc)/N_bc
                  tors3(2,i) = 0.0
                  done = .true.
                  goto 123
         else if ((crd(itb).eq.4.and.crd(itc).eq.3.and.
     &            ((val(itc).eq.4.or.val(itc).eq.34).or.
     &              mltb(itc).ne.0)).or.
     &            (crd(itc).eq.4.and.crd(itb).eq.3.and.
     &            ((val(itb).eq.4.or.val(itb).eq.34).or.
     &              mltb(itb).ne.0))) then
cn page 632, in "Merck Molecular Force Field V", e) & f)
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) = 0.0
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if ((crd(itb).eq.4.and.crd(itc).eq.2.and.
     &             (val(itc).eq.3.or.mltb(itc).ne.0)).or.
     &            (crd(itb).eq.4.and.crd(itc).eq.2.and.
     &             (val(itc).eq.3.or.mltb(itc).ne.0)))then
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) = 0.0
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if (crd(itb).eq.4.or.crd(itc).eq.4) then
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) = 0.0
                  tors2(2,i) = 180.0
                  tors3(1,i) = sqrt(Vb*Vc)/N_bc
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if (pilp(itb).eq.1.and.pilp(itc).eq.1)then
cn page 632, in "Merck Molecular Force Field V", g) 
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) = 0.0
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0
                  done = .true.
                  goto 123
         else if (pilp(itb).ne.0.and.mltb(itc).ne.0)then
              beta=6 
              if (mltb(itb).eq.1) then
               pi_bc=0.5
              else if (rowtab_b.eq.1.and.rowtab_c.eq.1) then
               pi_bc=0.3
              else if (rowtab_b.ne.1.or.rowtab_c.ne.1) then
               pi_bc=0.15
              endif
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) =beta*pi_bc*sqrt(Ub*Uc)
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if (pilp(itc).ne.0.and.mltb(itb).ne.0)then
              beta=6
              if (mltb(itc).eq.1) then
               pi_bc=0.5
              else if (rowtab_b.eq.1.and.rowtab_c.eq.1) then
               pi_bc=0.3
              else if (rowtab_b.ne.1.or.rowtab_c.ne.1) then
               pi_bc=0.15
              endif
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) =beta*pi_bc*sqrt(Ub*Uc)
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if ((mltb(itb).eq.1.or.mltb(itc).eq.1).and.
     &            (atomic(ib).ne.6.or.atomic(ic).ne.6)) then
              beta=6
              pi_bc=0.4
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) =beta*pi_bc*sqrt(Ub*Uc)
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if (mltb(itb).ne.0.and.mltb(itc).ne.0) then
c              write(*,*)'"an example being the central ',
c     &                  'C-C bond in butadiene or biphenyl"'
              beta=6
              pi_bc=0.15
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) =beta*pi_bc*sqrt(Ub*Uc)
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if (atomic(ib).eq.8.and.atomic(ic).eq.8) then
cn page 632, in "Merck Molecular Force Field V", h)
              W_b=2.0
              W_c=2.0
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) = -sqrt(W_b*W_c)
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if ((atomic(ib).eq.8.and.atomic(ic).eq.16).or.
     &            (atomic(ib).eq.16.and.atomic(ic).eq.8))then       
              W_b=2.0
              W_c=8.0
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) = -sqrt(W_b*W_c)
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else if (atomic(ib).eq.16.and.atomic(ic).eq.16) then
              W_b=8.0
              W_c=8.0
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) = -sqrt(W_b*W_c)
                  tors2(2,i) = 180.0
                  tors3(1,i) = 0.0
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         else
                  tors1(1,i) = 0.0
                  tors1(2,i) = 0.0
                  tors2(1,i) = 0.0
                  tors2(2,i) = 180.0
                  tors3(1,i) = sqrt(Vb*Vc)/N_bc
                  tors3(2,i) = 0.0

                  done = .true.
                  goto 123
         endif  
             

            if (use_tors) then
               if (use(ia) .or. use(ib) .or. use(ic) .or. use(id))
     &            abort = .true.
            end if
            if (header) then
               header = .false.
               write (iout,131)
  131          format (/,' Undefined Torsional Parameters :',
     &                 //,' Type',24x,'Atom Names',24x,
     &                    'Atom Classes',/)
            end if
            label = 'Torsion'
            if (iring .eq. 5)  label = '5-Ring '
            if (iring .eq. 4)  label = '4-Ring '
            write (iout,141)  label,ia,name(ia),ib,name(ib),ic,
     &                        name(ic),id,name(id),ita,itb,itc,itd
  141       format (1x,a7,4x,4(i6,'-',a3),5x,4i5)
         end if
      end do

      else 

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
         if (n12(ia).eq.1 .or. n12(id).eq.1)  goto 110
         if (use_ring4) then
            do j = 1, n12(id)
               if (ia .eq. i12(j,id)) then
                  iring = 4
                  goto 110
               end if
            end do
         end if
         if (use_ring5) then
            do j = 1, n12(id)
               if (ic .ne. i12(j,id)) then
                  do k = 1, n12(ia)
                     if (i12(j,id) .eq. i12(k,ia)) then
                        iring = 5
                        goto 110
                     end if
                  end do
               end if
            end do
         end if
  110    continue
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
                  goto 120
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
                  goto 120
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
                  goto 120
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
                  goto 120
               end if
            end do
         end if



c
c     warning if suitable torsional parameter not found
c
  120    continue
         minat = min(atomic(ia),atomic(ib),atomic(ic),atomic(id))
         if (minat .eq. 0)  done = .true.
         if (.not. done) then
            if (use_tors) then
               if (use(ia) .or. use(ib) .or. use(ic) .or. use(id))
     &            abort = .true.
            end if
            if (header) then
               header = .false.
               write (iout,130)
  130          format (/,' Undefined Torsional Parameters :',
     &                 //,' Type',24x,'Atom Names',24x,
     &                    'Atom Classes',/)
            end if
            label = 'Torsion'
            if (iring .eq. 5)  label = '5-Ring '
            if (iring .eq. 4)  label = '4-Ring '
            write (iout,140)  label,ia,name(ia),ib,name(ib),ic,
     &                        name(ic),id,name(id),ita,itb,itc,itd
  140       format (1x,a7,4x,4(i6,'-',a3),5x,4i5)
         end if
      end do

        endif
cn
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

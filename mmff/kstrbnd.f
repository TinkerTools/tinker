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
      include 'fields.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kstbnd.i'
      include 'merck.i'
      include 'potent.i'
      include 'ring.i'
      include 'strbnd.i'
      integer i,j,k,nsb
      integer size,next
      integer ia,ib,ic
      integer it,ita,itb,itc
      integer ina,inc
      integer nba,nbc
      integer nb1,nb2
      integer STBNT
      integer rowtab_a
      integer rowtab_b
      integer rowtab_c
      integer l,m
      integer AB,BC
      real*8 sb1,sb2,temp
      logical header
      logical ring3,ring4
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
c     assign stretch-bend parameters for each angle
c
      nstrbnd = 0

      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
cn
         if (forcefield .eq. 'MMFF94'
     &       .AND. lin(class(ib)).eq.0) then
cn stretch-bend interactions are omitted when
cn the bond angle is linear

          ita = class (ia)
          itb = class (ib)
          itc = class (ic)
          ina = atomic(ia)
          inc = atomic(ic)
          sbk(1,nstrbnd+1) = 0.0d0
          sbk(2,nstrbnd+1) = 0.0d0
          do k = 1, n12(ib)
            if (i12(k,ib) .eq. ia)  nb1 = bndlist(k,ib)
            if (i12(k,ib) .eq. ic)  nb2 = bndlist(k,ib)
          end do

          STBNT = 0
          AB = 0
          BC = 0
cn
cn STBNT : Stretch-Bend Type
cn in a "a-b-c" angle (with b being the central
cn atom and a having a class lower than c),
cn if the BT of a-b = 1, then STBNT =1
cn if the BT of b-c = 1, then STBNT =2
cn if both = 1, then STBNT = 3
cn 
cn if 3-membered ring, then STBNT =5
cn if 4-membered ring, then STBNT =4
cn if 3-membered ring with BT of a-b = 1, then STBNT =6
cn if 3-membered ring with BT of b-c = 1, then STBNT =7
cn if 3-membered ring with BT of both = 1, then STBNT =8
cn if 4-membered ring with BT of a-b = 1, then STBNT =9
cn if 4-membered ring with BT of b-c = 1, then STBNT =10
cn if 4-membered ring with BT of both = 1, then STBNT =11
cn
cn else, if all BT = 0 and no small ring, then STBNT =0
cn
cn looks if the 3 atoms belong to a same 3- or 4-membered cycle
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
       if (ring3.eqv..false.) then
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

          if (ita .lt. itc) then
cn looks at the BT          
           do j=1,nlignes
            if (((ia.eq.BT_1(j,1).AND.ib .eq. BT_1(j,2)).OR.
     &         (ib.eq.BT_1(j,1).AND.ia .eq. BT_1(j,2)))) then
               AB = 1
            endif
            if (((ic.eq.BT_1(j,1).AND.ib .eq. BT_1(j,2)).OR.
     &         (ib.eq.BT_1(j,1).AND.ic .eq. BT_1(j,2)))) then
               BC = 1
            endif
           enddo
           if (AB.eq.1.AND.BC.eq.0) STBNT = 1
           if (AB.eq.0.AND.BC.eq.1) STBNT = 2
           if (AB.eq.1.AND.BC.eq.1) STBNT = 3
cn looks at the small rings
           if (STBNT.eq.0.AND. 
     &          ring3.eqv..true.) then           
cn sum of BT = 0 and 3-membered ring
              STBNT = 5
           else if (STBNT.eq.1.AND.
     &          ring3.eqv..true.) then
cn BT a-b = 1 and 3-membered ring     
              STBNT = 6
           else if (STBNT.eq.2.AND.
     &          ring3.eqv..true.) then
cn BT b-c = 1 and 3-membered ring     
              STBNT = 7
           else if (STBNT.eq.3.AND.
     &          ring3.eqv..true.) then
cn sum of BT = 2 and 3-membered ring     
              STBNT = 8
           else if (STBNT.eq.0.AND.
     &      ring4.eqv..true.) then
cn sum of BT = 0 and 4-membered ring     
              STBNT = 4
           else if (STBNT.eq.1.AND.
     &      ring4.eqv..true.) then
              STBNT = 9
cn BT a-b = 1 and 4-membered ring 
           else if (STBNT.eq.2.AND.
     &      ring4.eqv..true.) then
cn BT b-c = 1 and 4-membered ring     
              STBNT = 10
           else if (STBNT.eq.3.AND.
     &      ring4.eqv..true.) then
cn sum of BT = 2 and 4-membered ring     
              STBNT = 11
           endif

          else if (ita .gt. itc) then
cn looks at the BT         
           do j=1,nlignes
            if (((ia.eq.BT_1(j,1).AND.ib .eq. BT_1(j,2)).OR.
     &         (ib.eq.BT_1(j,1).AND.ia .eq. BT_1(j,2)))) then
               AB = 1
            endif
            if (((ic.eq.BT_1(j,1).AND.ib .eq. BT_1(j,2)).OR.
     &         (ib.eq.BT_1(j,1).AND.ic .eq. BT_1(j,2)))) then
               BC = 1
            endif
           enddo
           if (AB.eq.1.AND.BC.eq.0) STBNT = 2
           if (AB.eq.0.AND.BC.eq.1) STBNT = 1
           if (AB.eq.1.AND.BC.eq.1) STBNT = 3

cn looks at the small rings
           if (STBNT.eq.0.AND.
     &          ring3.eqv..true.) then
cn sum of BT = 0 and 3-membered ring
              STBNT = 5
           else if (STBNT.eq.1.AND.
     &          ring3.eqv..true.) then
cn BT b-c = 1 and 3-membered ring
              STBNT = 6
           else if (STBNT.eq.2.AND.
     &          ring3.eqv..true.) then
cn BT a-b = 1 and 3-membered ring
              STBNT = 7
           else if (STBNT.eq.3.AND.
     &          ring3.eqv..true.) then
cn sum of BT = 2 and 3-membered ring
              STBNT = 8
           else if (STBNT.eq.0.AND.
     &            ring4.eqv..true.) then
cn sum of BT = 0 and 4-membered ring
              STBNT = 4
           else if (STBNT.eq.1.AND.
     &             ring4.eqv..true.) then
              STBNT = 9
cn BT b-c = 1 and 4-membered ring
           else if (STBNT.eq.2.AND.
     &              ring4.eqv..true.) then
cn BT a-b = 1 and 4-membered ring
              STBNT = 10
           else if (STBNT.eq.3.AND.
     &              ring4.eqv..true.) then
cn sum of BT = 2 and 4-membered ring
              STBNT = 11
           endif
           
          else if (ita .eq. itc) then
           do j=1,nlignes
            if (((ic.eq.BT_1(j,1).AND.ib .eq. BT_1(j,2)).OR.
     &         (ib.eq.BT_1(j,1).AND.ic .eq. BT_1(j,2)))) then
               BC = 1
            endif
            if (((ia.eq.BT_1(j,1).AND.ib .eq. BT_1(j,2)).OR.
     &         (ib.eq.BT_1(j,1).AND.ia .eq. BT_1(j,2)))) then
               AB = 1
            endif
           enddo
           if (AB.eq.1.AND.BC.eq.0) STBNT = 1 
           if (AB.eq.0.AND.BC.eq.1) STBNT = 2 
           if (AB.eq.1.AND.BC.eq.1) STBNT = 3

cn looks at the small rings
           if (STBNT.eq.0.AND.
     &          ring3.eqv..true.) then
              STBNT = 5
           else if (STBNT.eq.1.AND.
     &          ring3.eqv..true.) then
              STBNT = 6
           else if (STBNT.eq.2.AND.
     &          ring3.eqv..true.) then
              STBNT = 7
           else if (STBNT.eq.3.AND.
     &          ring3.eqv..true.) then
              STBNT = 8
           else if (STBNT.eq.0.AND.
     &              ring4.eqv..true.) then
              STBNT = 4
            else if (STBNT.eq.1.AND.
     &             ring4.eqv..true.) then
              STBNT = 9
           else if (STBNT.eq.2.AND.
     &              ring4.eqv..true.) then
              STBNT = 10
           else if (STBNT.eq.3.AND.
     &              ring4.eqv..true.) then
              STBNT = 11
           endif

          endif
cn
cn if the parameters aren't available,
cn we need to know the line of the elements in
cn the periodic table of the elements
cn in order to use the
cn "empirical rule"
cn (cf T. A. Halgren, "Merck Molecular Force Field. V.
cn Extension of MMFF94 Using Experimental Data,
cn Additional Computational Data, and Empirical Rules",
cn J. Comput. Chem. 17, pp. 616-641)
cn
          if(atomic(ia).eq.1) rowtab_a=0
          if((atomic(ia).ge.3).and.(atomic(ia).le.10)) rowtab_a=1
          if((atomic(ia).ge.11).and.(atomic(ia).le.18)) rowtab_a=2
          if((atomic(ia).ge.19).and.(atomic(ia).le.36)) rowtab_a=3
          if((atomic(ia).ge.37).and.(atomic(ia).le.54)) rowtab_a=4
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
          if (STBNT .eq. 11) then
               if ((stbn_abc11(ita,itb,itc).ne.1000.0d0)
     &          .AND. (stbn_cba11(ita,itb,itc).ne.1000.0d0)) then
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc11(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba11(ita,itb,itc)
               else
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 10) then
               if ((stbn_abc10(ita,itb,itc).ne.1000.0d0)
     &          .AND. (stbn_cba10(ita,itb,itc).ne.1000.0d0)) then
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc10(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba10(ita,itb,itc)
               else
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 9) then
               if ((stbn_abc9(ita,itb,itc).ne.1000.0d0)
     &          .AND. (stbn_cba9(ita,itb,itc).ne.1000.0d0)) then
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc9(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba9(ita,itb,itc)
               else
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2

                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 8) then
               if ((stbn_abc8(ita,itb,itc).ne.1000.0d0)
     &          .AND. (stbn_cba3(ita,itb,itc).ne.1000.0d0)) then
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc8(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba8(ita,itb,itc)
               else
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 7) then
               if ((stbn_abc7(ita,itb,itc).ne.1000.0d0)
     &          .AND. (stbn_cba7(ita,itb,itc).ne.1000.0d0)) then
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc7(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba7(ita,itb,itc)
               else
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 6) then
               if ((stbn_abc6(ita,itb,itc).ne.1000.0d0)
     &          .AND. (stbn_cba3(ita,itb,itc).ne.1000.0d0)) then
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc6(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba6(ita,itb,itc)
               else
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 5) then
               if (((stbn_abc5(ita,itb,itc).ne.1000.0d0)
     &          .AND. (stbn_cba3(ita,itb,itc).ne.1000.0d0))
     &          .OR.(ita.eq.22.AND.itb.eq.22.AND.itc.eq.22)) then
cn angle 22-22-22 with STBNT 5 has ksb_abc =  ksb_cba = 0     
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc5(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba5(ita,itb,itc)
               else
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 4) then
               if ((stbn_abc4(ita,itb,itc).ne.1000.0d0)
     &          .AND. (stbn_cba4(ita,itb,itc).ne.1000.0d0))  then
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc4(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba4(ita,itb,itc)
               else
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 3) then
               if ((stbn_abc3(ita,itb,itc).ne.1000.0d0)
     &          .AND. (stbn_cba3(ita,itb,itc).ne.1000.0d0))  then
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc3(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba3(ita,itb,itc)
               else
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                
                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 2) then
              if ((stbn_abc2(ita,itb,itc).ne.1000.0d0)
     &         .AND. (stbn_cba2(ita,itb,itc).ne.1000.0d0)) then
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc2(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba2(ita,itb,itc)
               else
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2

                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 1) then
               if ((stbn_abc1(ita,itb,itc).ne.1000.0d0)
     &          .AND. (stbn_cba1(ita,itb,itc).ne.1000.0d0)) then
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc1(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba1(ita,itb,itc)
               else 
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2

                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          else if (STBNT .eq. 0) then
               if (((stbn_abc(ita,itb,itc) .ne. 1000.0d0)
     &          .AND. (stbn_cba(ita,itb,itc).ne. 1000.0d0))
     &          .OR.(ita.eq.12.AND.itb.eq.20.AND.itc.eq.20)
     &          .OR.(ita.eq.20.AND.itb.eq.20.AND.itc.eq.12)) then
cn parameters for stretch-bend for types 12-20-20 (CBA) = 0 
cn (parameter not missing)
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = stbn_abc(ita,itb,itc)
                sbk(2,nstrbnd) = stbn_cba(ita,itb,itc)
               else
cn empirical rule for stretch-bend applies               
                nstrbnd = nstrbnd + 1
                isb(1,nstrbnd) = i
                isb(2,nstrbnd) = nb1
                isb(3,nstrbnd) = nb2
                sbk(1,nstrbnd) = defstbnd_abc(rowtab_a,
     &                             rowtab_b,rowtab_c)
                sbk(2,nstrbnd) = defstbnd_cba(rowtab_a,
     &                             rowtab_b,rowtab_c)
               end if
          endif
         else
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
         endif
      end do
c
c     turn off the stretch-bend potential if it is not used
c
      if (nstrbnd .eq. 0)  use_strbnd = .false.
      return
      end

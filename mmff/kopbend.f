c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kopbend  --  out-of-plane bending parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kopbend" assigns the force constants for out-of-plane
c     bending at trigonal centers via Wilson-Decius-Cross angle
c     bends; also processes any new or changed parameter values
c
c
      subroutine kopbend
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'fields.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'kopbnd.i'
      include 'merck.i'
      include 'opbend.i'
      include 'potent.i'
      integer i,j,nopb,it
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer itta,ittb
      integer ittc,ittd
      integer number
      integer next,size
      integer collclass
      real*8 fopb
      logical header
      logical jopb(maxclass)
      character*4 pa,pb,pc,pd
      character*4 zero4
      character*8 zero8
      character*16 blank,pt
      character*16 pt0,pt00
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     process keywords containing out-of-plane bend parameters
c
      blank = '                '
      zero4 = '0000'
      zero8 = '00000000'
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'OPBEND ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            fopb = 0.0d0
            string = record(next:120)
            read (string,*,err=10,end=10)  ia,ib,ic,id,fopb
   10       continue
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            if (ic .le. id) then
               pt = pa//pb//pc//pd
            else
               pt = pa//pb//pd//pc
            end if
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' Additional Out-of-Plane Bend Parameters :',
     &                 //,5x,'Atom Classes',19x,'K(OPB)',/)
            end if
            write (iout,30)  ia,ib,ic,id,fopb
   30       format (4x,4i4,10x,f12.3)
            size = 4
            do j = 1, maxnopb
               if (kopb(j).eq.blank .or. kopb(j).eq.pt) then
                  kopb(j) = pt
                  opbn(j) = fopb
                  goto 50
               end if
            end do
            write (iout,40)
   40       format (/,' KOPBEND --  Too many Out-of-Plane',
     &                 ' Angle Bending Parameters')
            abort = .true.
   50       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nopb = maxnopb
      do i = maxnopb, 1, -1
         if (kopb(i) .eq. blank)  nopb = i - 1
      end do
c
c     jump to the parameter assignment phase for MMFF
c
      if (forcefield .eq. 'MMFF94')  goto 100
c
c     make list of atom classes using out-of-plane bending
c
      do i = 1, maxclass
         jopb(i) = .false.
      end do
      do i = 1, maxnopb
         if (kopb(i) .eq. blank)  goto 60
         it = number (kopb(i)(5:8))
         jopb(it) = .true.
      end do
   60 continue
c
c     assign out-of-plane bending parameters for each angle
c
      nopbend = 0
      if (nopb .ne. 0) then
         header = .true.
         do i = 1, nangle
            ib = iang(2,i)
            itb = class(ib)
            if (jopb(itb) .and. n12(ib).eq.3) then
               ia = iang(1,i)
               ita = class(ia)
               ic = iang(3,i)
               itc = class(ic)
               id = iang(4,i)
               itd = class(id)
               size = 4
               call numeral (ita,pa,size)
               call numeral (itb,pb,size)
               call numeral (itc,pc,size)
               call numeral (itd,pd,size)
               if (ita .le. itc) then
                  pt = pd//pb//pa//pc
               else
                  pt = pd//pb//pc//pa
               end if
               pt0 = pd//pb//zero8
               pt00 = zero4//pb//zero8
               do j = 1, nopb
                  if (kopb(j) .eq. pt) then
                     nopbend = nopbend + 1
                     iopb(nopbend) = i
                     opbk(nopbend) = opbn(j)
                     goto 90
                  end if
               end do
               do j = 1, nopb
                  if (kopb(j) .eq. pt0) then
                     nopbend = nopbend + 1
                     iopb(nopbend) = i
                     opbk(nopbend) = opbn(j)
                     goto 90
                  end if
               end do
               do j = 1, nopb
                  if (kopb(j) .eq. pt00) then
                     nopbend = nopbend + 1
                     iopb(nopbend) = i
                     opbk(nopbend) = opbn(j)
                     goto 90
                  end if
               end do
               abort = .true.
               if (header) then
                  header = .false.
                  write (iout,70)
   70             format (/,' Undefined Out-of-Plane Bend Parameters :',
     &                    //,' Type',24x,'Atom Names',24x,
     &                       'Atom Classes',/)
               end if
               write (iout,80)  id,name(id),ib,name(ib),ia,name(ia),
     &                          ic,name(ic),itd,itb,ita,itc
   80          format (' Angle-OP',3x,4(i6,'-',a3),5x,4i5)
   90          continue
            else
               iang(4,i) = ib
            end if
         end do
      end if
c
c     mark angles at trigonal sites to use projected in-plane values
c
      do i = 1, nopbend
         j = iopb(i)
      if (angtyp(j) .eq. 'HARMONIC')  angtyp(j) = 'IN-PLANE'
      end do
c
c     assign out-of-plane bending parameters for MMFF
c
  100 continue
      if (forcefield .eq. 'MMFF94') then
         nopbend = 0
         header = .true.
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            id = iang(4,i)
            itta = type(ia)
            ittb = type(ib)
            ittc = type(ic)
            ittd = type(id)
            collclass = 0
 110        continue
            collclass = collclass + 1
            if (collclass .eq. 1) then
               ita = eqclass(itta,1)
               itb = eqclass(ittb,1)
               itc = eqclass(ittc,1)
               itd = eqclass(ittd,1)
            else if (collclass .eq. 2) then
               ita = eqclass(itta,2)
               itb = eqclass(ittb,2)
               itc = eqclass(ittc,2)
               itd = eqclass(ittd,2)
            else if (collclass .eq. 3) then
               ita = eqclass(itta,3)
               itb = eqclass(ittb,2)
               itc = eqclass(ittc,3)
               itd = eqclass(ittd,3)
            else if (collclass .eq. 4) then
               ita = eqclass(itta,4)
               itb = eqclass(ittb,2)
               itc = eqclass(ittc,4)
               itd = eqclass(ittd,4)
            else if (collclass .eq. 5) then
               ita = eqclass(itta,5)
               itb = eqclass(ittb,2)
               itc = eqclass(ittc,5)
               itd = eqclass(ittd,5)
            end if
            if (ia.ne.0 .and. ib.ne.0 .and.
     &          ic.ne.0 .and. id.ne.0) then
                if (collclass.gt.5) then
                   nopbend = nopbend + 1
                   iopb(nopbend) = i
                   opbk(nopbend) = 0.0d0
                else
                   size = 4
                   call numeral (ita,pa,size)
                   call numeral (itb,pb,size)
                   call numeral (itc,pc,size)
                   call numeral (itd,pd,size)
                   if (itd.le.ita .and. itd.le.itc) then
                      if (ita .le. itc) then
                         pt = pd//pb//pa//pc
                      else
                         pt = pd//pb//pc//pa
                      end if
                   else if (ita.le.itc .and. ita.le.itd) then
                      if (itd .le. itc) then
                         pt = pa//pb//pd//pc
                      else
                         pt = pa//pb//pc//pd
                      end if
                   else if (itc.le.ita .and. itc.le.itd) then
                      if (ita .le. itd) then
                         pt = pc//pb//pa//pd
                      else
                         pt = pc//pb//pd//pa
                      end if
                   end if
                   do j = 1, nopb
                      if (kopb(j) .eq. pt) then
                         nopbend = nopbend + 1
                         iopb(nopbend) = i
                         opbk(nopbend) = opbn(j)
                         goto 120
                      end if
                   end do
                   if (class(ib).eq.8 .or. class(ib).eq.17 .or.
     &                 class(ib).eq.26 .or. class(ib).eq.43 .or.
     &                 class(ib).eq.49 .or. class(ib).eq.73 .or.
     &                 class(ib).eq.82) then
                      nopbend = nopbend + 1
                      iopb(nopbend) = i
                      opbk(nopbend) = 0.0d0
                      goto 120
                   end if
                   goto 110
  120              continue
                end if
             end if
         end do
      end if
c
c     turn off the out-of-plane bending term if it is not used
c
      if (nopbend .eq. 0)  use_opbend = .false.
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtprm  --  output of force field parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtprm" writes out a formatted listing of the default
c     set of potential energy parameters for a force field
c
c
      subroutine prtprm (itxt)
      use angpot
      use bndpot
      use chgpot
      use fields
      use kanang
      use kangs
      use kantor
      use katoms
      use kbonds
      use kcflux
      use kchrge
      use kcpen
      use kctrn
      use kdipol
      use kdsp
      use kexpl
      use khbond
      use kiprop
      use kitors
      use kmulti
      use kopbnd
      use kopdst
      use korbs
      use kpitor
      use kpolr
      use krepl
      use ksolut
      use kstbnd
      use ksttor
      use ktorsn
      use ktrtor
      use kurybr
      use kvdws
      use kvdwpr
      use mplpot
      use polpot
      use sizes
      use urypot
      use vdwpot
      implicit none
      integer i,j,k,itxt
      integer number,npg
      integer k1,k2,k3
      integer k4,k5
      integer fold(6)
      real*8 ampli(6)
      real*8 phase(6)
      logical exist
      character*1 formfeed
      character*3 blank3
      character*8 blank8
      character*12 blank12
      character*16 blank16
      character*20 blank20
c
c
c     define blank character strings of various lengths
c
      blank3 = '   '
      blank8 = '        '
      blank12 = '            '
      blank16 = '                '
      blank20 = '                    '
c
c     set the string value of the formfeed character (Ctrl-L)
c
      formfeed = char(12)
c
c     force field atom type definitions
c
      exist = .false.
      do i = 1, maxtyp
         if (symbol(i) .ne. blank3)  exist = .true.
      end do
      if (exist) then
         write (itxt,10)  forcefield
   10    format (//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,20)
   20    format (//,15x,'Force Field Atom Definitions',
     &           //,54x,'Atomic',4x,'Atomic',
     &           /,5x,'Type',3x,'Class',3x,'Symbol',3x,'Description',
     &              14x,'Number',4x,'Weight',3x,'Valence',/)
         do i = 1, maxtyp
            if (symbol(i) .ne. blank3) then
               write (itxt,30)  i,atmcls(i),symbol(i),describe(i),
     &                          atmnum(i),weight(i),ligand(i)
   30          format (3x,i5,3x,i5,5x,a3,5x,a24,i5,f12.3,i7)
            end if
         end do
      end if
c
c     bond stretching parameters
c
      if (kb(1) .ne. blank8) then
         write (itxt,40)  formfeed,forcefield
   40    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,50)
   50    format (//,15x,'Bond Stretching Parameters',
     &           ///,22x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb
            if (kb(i) .eq. blank8)  goto 70
            k1 = number(kb(i)(1:4))
            k2 = number(kb(i)(5:8))
            write (itxt,60)  i,k1,k2,bcon(i),blen(i)
   60       format (8x,i7,5x,i4,'-',i4,6x,f12.3,f12.4)
         end do
   70    continue
      end if
c
c     bond stretching parameters for 5-membered rings
c
      if (kb5(1) .ne. blank8) then
         write (itxt,80)
   80    format (//,15x,'5-Membered Ring Stretch Parameters',
     &           ///,22x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb5
            if (kb5(i) .eq. blank8)  goto 100
            k1 = number(kb5(i)(1:4))
            k2 = number(kb5(i)(5:8))
            write (itxt,90)  i,k1,k2,bcon5(i),blen5(i)
   90       format (8x,i7,5x,i4,'-',i4,6x,f12.3,f12.4)
         end do
  100    continue
      end if
c
c     bond stretching parameters for 4-membered rings
c
      if (kb4(1) .ne. blank8) then
         write (itxt,110)
  110    format (//,15x,'4-Membered Ring Stretch Parameters',
     &           ///,22x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb4
            if (kb4(i) .eq. blank8)  goto 130
            k1 = number(kb4(i)(1:4))
            k2 = number(kb4(i)(5:8))
            write (itxt,120)  i,k1,k2,bcon4(i),blen4(i)
  120       format (8x,i7,5x,i4,'-',i4,6x,f12.3,f12.4)
         end do
  130    continue
      end if
c
c     bond stretching parameters for 3-membered rings
c
      if (kb3(1) .ne. blank8) then
         write (itxt,140)
  140    format (//,15x,'3-Membered Ring Stretch Parameters',
     &           ///,22x,'Classes',15x,'KS',7x,'Length',/)
         do i = 1, maxnb3
            if (kb3(i) .eq. blank8)  goto 160
            k1 = number(kb3(i)(1:4))
            k2 = number(kb3(i)(5:8))
            write (itxt,150)  i,k1,k2,bcon3(i),blen3(i)
  150       format (8x,i7,5x,i4,'-',i4,6x,f12.3,f12.4)
         end do
  160    continue
      end if
c
c     cubic and quartic bond stretching parameters
c
      if (cbnd.ne.0.0d0 .or. qbnd.ne.0.0d0) then
         write (itxt,170)  cbnd,qbnd
  170    format (//,15x,'Higher Order Stretching Constants',
     &           ///,20x,'Cubic',f17.3,/,20x,'Quartic',f15.3)
      end if
c
c     electronegativity bond length correction parameters
c
      if (kel(1) .ne. blank12) then
         write (itxt,180)
  180    format (//,15x,'Electronegativity Bond Length Parameters',
     &           ///,25x,'Classes',21x,'dLength',/)
         do i = 1, maxnel
            if (kel(i) .eq. blank12)  goto 200
            k1 = number(kel(i)(1:4))
            k2 = number(kel(i)(5:8))
            k3 = number(kel(i)(9:12))
            write (itxt,190)  i,k1,k2,k3,dlen(i)
  190       format (8x,i7,5x,i4,'-',i4,'-',i4,14x,f12.4)
         end do
  200    continue
      end if
c
c     bond angle bending parameters
c
      if (ka(1) .ne. blank12) then
         write (itxt,210)  formfeed,forcefield
  210    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,220)
  220    format (//,15x,'Angle Bending Parameters',
     &           ///,18x,'Classes',11x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,44x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do i = 1, maxna
            if (ka(i) .eq. blank12)  goto 250
            k1 = number(ka(i)(1:4))
            k2 = number(ka(i)(5:8))
            k3 = number(ka(i)(9:12))
            if (ang(2,i).eq.0.0d0 .and. ang(3,i).eq.0.0d0) then
               write (itxt,230)  i,k1,k2,k3,acon(i),ang(1,i)
  230          format (3x,i5,5x,i4,'-',i4,'-',i4,2f12.3)
            else
               write (itxt,240)  i,k1,k2,k3,acon(i),(ang(j,i),j=1,3)
  240          format (3x,i5,5x,i4,'-',i4,'-',i4,4f12.3)
            end if
         end do
  250    continue
      end if
c
c     angle bending parameters for 5-membered rings
c
      if (ka5(1) .ne. blank12) then
         write (itxt,260)
  260    format (//,17x,'5-Membered Ring Bend Parameters',
     &           ///,18x,'Classes',11x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,44x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do i = 1, maxna5
            if (ka5(i) .eq. blank12)  goto 290
            k1 = number(ka5(i)(1:4))
            k2 = number(ka5(i)(5:8))
            k3 = number(ka5(i)(9:12))
            if (ang5(2,i).eq.0.0d0 .and. ang5(3,i).eq.0.0d0) then
               write (itxt,270)  i,k1,k2,k3,acon5(i),ang5(1,i)
  270          format (3x,i5,5x,i4,'-',i4,'-',i4,2f12.3)
            else
               write (itxt,280)  i,k1,k2,k3,acon5(i),(ang5(j,i),j=1,3)
  280          format (3x,i5,5x,i4,'-',i4,'-',i4,4f12.3)
            end if
         end do
  290    continue
      end if
c
c     angle bending parameters for 4-membered rings
c
      if (ka4(1) .ne. blank12) then
         write (itxt,300)
  300    format (//,15x,'4-Membered Ring Bend Parameters',
     &           ///,18x,'Classes',11x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,44x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do i = 1, maxna4
            if (ka4(i) .eq. blank12)  goto 330
            k1 = number(ka4(i)(1:4))
            k2 = number(ka4(i)(5:8))
            k3 = number(ka4(i)(9:12))
            if (ang4(2,i).eq.0.0d0 .and. ang4(3,i).eq.0.0d0) then
               write (itxt,310)  i,k1,k2,k3,acon4(i),ang4(1,i)
  310          format (3x,i5,5x,i4,'-',i4,'-',i4,2f12.3)
            else
               write (itxt,320)  i,k1,k2,k3,acon4(i),(ang4(j,i),j=1,3)
  320          format (3x,i5,5x,i4,'-',i4,'-',i4,4f12.3)
            end if
         end do
  330    continue
      end if
c
c     angle bending parameters for 3-membered rings
c
      if (ka3(1) .ne. blank12) then
         write (itxt,340)
  340    format (//,15x,'3-Membered Ring Bend Parameters',
     &           ///,18x,'Classes',11x,'KB',6x,'Value 1',
     &              5x,'Value 2',5x,'Value 3',
     &           /,44x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         do  i = 1, maxna3
            if (ka3(i) .eq. blank12)  goto 370
            k1 = number(ka3(i)(1:4))
            k2 = number(ka3(i)(5:8))
            k3 = number(ka3(i)(9:12))
            if (ang3(2,i).eq.0.0d0 .and. ang3(3,i).eq.0.0d0) then
               write (itxt,350)  i,k1,k2,k3,acon3(i),ang3(1,i)
  350          format (3x,i5,5x,i4,'-',i4,'-',i4,2f12.3)
            else
               write (itxt,360)  i,k1,k2,k3,acon3(i),(ang3(j,i),j=1,3)
  360          format (3x,i5,5x,i4,'-',i4,'-',i4,4f12.3)
            end if
         end do
  370    continue
      end if
c
c     in-plane projected angle bending parameters
c
      if (kap(1) .ne. blank12) then
         write (itxt,380)
  380    format (//,15x,'In-Plane Angle Bending Parameters',
     &           ///,18x,'Classes',11x,'KB',6x,'Value 1',5x,'Value 2',
     &           /,45x,'(X-R)',7x,'(X-H)'/)
         do  i = 1, maxnap
            if (kap(i) .eq. blank12)  goto 400
            k1 = number(kap(i)(1:4))
            k2 = number(kap(i)(5:8))
            k3 = number(kap(i)(9:12))
            write (itxt,390)  i,k1,k2,k3,aconp(i),(angp(j,i),j=1,2)
  390       format (3x,i5,5x,i4,'-',i4,'-',i4,3f12.3)
         end do
  400    continue
      end if
c
c     Fourier bond angle bending parameters
c
      if (kaf(1) .ne. blank12) then
         write (itxt,410)
  410    format (//,15x,'Fourier Angle Bending Parameters',
     &           ///,18x,'Classes',11x,'KB',8x,'Shift',6x,'Period',/)
         do  i = 1, maxnaf
            if (kaf(i) .eq. blank12)  goto 430
            k1 = number(kaf(i)(1:4))
            k2 = number(kaf(i)(5:8))
            k3 = number(kaf(i)(9:12))
            write (itxt,420)  i,k1,k2,k3,aconf(i),(angf(j,i),j=1,2)
  420       format (3x,i5,5x,i4,'-',i4,'-',i4,3f12.3)
         end do
  430    continue
      end if
c
c     cubic through sextic bond angle bending parameters
c
      if (cang.ne.0.0d0 .or. qang.ne.0.0d0 .or.
     &    pang.ne.0.0d0 .or. sang.ne.0.0d0) then
         write (itxt,440)  cang,qang,pang,sang
  440    format (//,15x,'Higher Order Bending Constants',
     &           ///,20x,'Cubic',d17.3,/,20x,'Quartic',d15.3,
     &           /,20x,'Pentic',d16.3,/,20x,'Sextic',d16.3)
      end if
c
c     stretch-bend parameters
c
      if (ksb(1) .ne. blank12) then
         write (itxt,450)  formfeed,forcefield
  450    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,460)
  460    format (//,15x,'Stretch-Bend Parameters',
     &           ///,18x,'Classes',18x,'KSB1',8x,'KSB2',/)
         do i = 1, maxnsb
            if (ksb(i) .eq. blank12)  goto 480
            k1 = number(ksb(i)(1:4))
            k2 = number(ksb(i)(5:8))
            k3 = number(ksb(i)(9:12))
            write (itxt,470)  i,k1,k2,k3,stbn(1,i),stbn(2,i)
  470       format (3x,i5,5x,i4,'-',i4,'-',i4,8x,2f12.3)
         end do
  480    continue
      end if
c
c     Urey-Bradley parameters
c
      if (ku(1) .ne. blank12) then
         write (itxt,490)  formfeed,forcefield
  490    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,500)
  500    format (//,15x,'Urey-Bradley Parameters',
     &           ///,18x,'Classes',19x,'KB',6x,'Distance',/)
         do i = 1, maxnu
            if (ku(i) .eq. blank12)  goto 520
            k1 = number(ku(i)(1:4))
            k2 = number(ku(i)(5:8))
            k3 = number(ku(i)(9:12))
            write (itxt,510)  i,k1,k2,k3,ucon(i),dst13(i)
  510       format (3x,i5,5x,i4,'-',i4,'-',i4,8x,f12.3,f12.4)
         end do
  520    continue
      end if
c
c     cubic and quartic Urey-Bradley parameters
c
      if (cury.ne.0.0d0 .or. qury.ne.0.0d0) then
         write (itxt,530)  cury,qury
  530    format (//,15x,'Higher Order Urey-Bradley Constants',
     &           ///,20x,'Cubic',f17.3,/,20x,'Quartic',f15.3)
      end if
c
c     angle-angle parameters
c
      exist = .false.
      do i = 1, maxclass
         do k = 1, 3
            if (anan(k,i) .ne. 0.0d0)  exist = .true.
         end do
      end do
      if (exist) then
         write (itxt,540)  formfeed,forcefield
  540    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,550)
  550    format (//,15x,'Angle-Angle Parameters',
     &           ///,20x,'Class',9x,'KAA 1',7x,'KAA 2',7x,'KAA 3',
     &           /,33x,'(R-X-R)',5x,'(R-X-H)',5x,'(H-X-H)',/)
         k = 0
         do i = 1, maxclass
            if (anan(1,i).ne.0.0d0 .or. anan(2,i).ne.0.0d0
     &               .or. anan(3,i).ne.0.0d0) then
               k = k + 1
               write (itxt,560)  k,i,(anan(j,i),j=1,3)
  560          format (6x,i7,4x,i7,3x,3f12.3)
            end if
         end do
      end if
c
c     out-of-plane bending parameters
c
      if (kopb(1) .ne. blank16) then
         write (itxt,570)  formfeed,forcefield
  570    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,580)
  580    format (//,15x,'Out-of-Plane Bend Parameters',
     &           ///,26x,'Classes',11x,'KOPB',/)
         do i = 1, maxnopb
            if (kopb(i) .eq. blank16)  goto 600
            k1 = number(kopb(i)(1:4))
            k2 = number(kopb(i)(5:8))
            k3 = number(kopb(i)(9:12))
            k4 = number(kopb(i)(13:16))
            write (itxt,590)  i,k1,k2,k3,k4,opbn(i)
  590       format (6x,i7,5x,i4,'-',i4,'-',i4,'-',i4,f12.3)
         end do
  600    continue
      end if
c
c     out-of-plane distance parameters
c
      if (kopd(1) .ne. blank16) then
         write (itxt,610)  formfeed,forcefield
  610    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,620)
  620    format (//,15x,'Out-of-Plane Distance Parameters',
     &           ///,26x,'Classes',11x,'KOPD',/)
         do i = 1, maxnopd
            if (kopd(i) .eq. blank16)  goto 640
            k1 = number(kopd(i)(1:4))
            k2 = number(kopd(i)(5:8))
            k3 = number(kopd(i)(9:12))
            k4 = number(kopd(i)(13:16))
            write (itxt,630)  i,k1,k2,k3,k4,opds(i)
  630       format (6x,i7,5x,i4,'-',i4,'-',i4,'-',i4,f12.3)
         end do
  640    continue
      end if
c
c     improper dihedral parameters
c
      if (kdi(1) .ne. blank16) then
         write (itxt,650)  formfeed,forcefield
  650    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,660)
  660    format (//,15x,'Improper Dihedral Parameters',
     &           ///,20x,'Classes',12x,'KID',7x,'Target',/)
         do i = 1, maxndi
            if (kdi(i) .eq. blank16)  goto 680
            k1 = number(kdi(i)(1:4))
            k2 = number(kdi(i)(5:8))
            k3 = number(kdi(i)(9:12))
            k4 = number(kdi(i)(13:16))
            write (itxt,670)  i,k1,k2,k3,k4,dcon(i),tdi(i)
  670       format (2x,i5,5x,i4,'-',i4,'-',i4,'-',i4,f12.3,f12.4)
         end do
  680    continue
      end if
c
c     improper torsional parameters
c
      if (kti(1) .ne. blank16) then
         write (itxt,690)  formfeed,forcefield
  690    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,700)
  700    format (//,15x,'Improper Torsion Parameters',
     &           ///,17x,'Classes',15x,'KTI Values',/)
         do i = 1, maxnti
            if (kti(i) .eq. blank16)  goto 720
            k1 = number(kti(i)(1:4))
            k2 = number(kti(i)(5:8))
            k3 = number(kti(i)(9:12))
            k4 = number(kti(i)(13:16))
            j = 0
            if (ti1(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = ti1(1,i)
               phase(j) = ti1(2,i)
            end if
            if (ti2(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = ti2(1,i)
               phase(j) = ti2(2,i)
            end if
            if (ti3(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = ti3(1,i)
               phase(j) = ti3(2,i)
            end if
            write (itxt,710)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  710       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,2x,3(f8.3,f6.1,i2))
         end do
  720    continue
      end if
c
c     torsional angle parameters
c
      if (kt(1) .ne. blank16) then
         write (itxt,730)  formfeed,forcefield
  730    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,740)
  740    format (//,15x,'Torsional Parameters',
     &           ///,17x,'Classes',15x,'KT Values',/)
         do i = 1, maxnt
            if (kt(i) .eq. blank16)  goto 760
            k1 = number(kt(i)(1:4))
            k2 = number(kt(i)(5:8))
            k3 = number(kt(i)(9:12))
            k4 = number(kt(i)(13:16))
            j = 0
            if (t1(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = t1(1,i)
               phase(j) = t1(2,i)
            end if
            if (t2(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = t2(1,i)
               phase(j) = t2(2,i)
            end if
            if (t3(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = t3(1,i)
               phase(j) = t3(2,i)
            end if
            if (t4(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 4
               ampli(j) = t4(1,i)
               phase(j) = t4(2,i)
            end if
            if (t5(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 5
               ampli(j) = t5(1,i)
               phase(j) = t5(2,i)
            end if
            if (t6(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 6
               ampli(j) = t6(1,i)
               phase(j) = t6(2,i)
            end if
            write (itxt,750)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  750       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,2x,6(f8.3,f6.1,i2))
         end do
  760    continue
      end if
c
c     torsional angle parameters for 5-membered rings
c
      if (kt5(1) .ne. blank16) then
         write (itxt,770)
  770    format (//,15x,'5-Membered Ring Torsion Parameters',
     &           ///,17x,'Classes',15x,'KT Values',/)
         do i = 1, maxnt5
            if (kt5(i) .eq. blank16)  goto 790
            k1 = number(kt5(i)(1:4))
            k2 = number(kt5(i)(5:8))
            k3 = number(kt5(i)(9:12))
            k4 = number(kt5(i)(13:16))
            j = 0
            if (t15(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = t15(1,i)
               phase(j) = t15(2,i)
            end if
            if (t25(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = t25(1,i)
               phase(j) = t25(2,i)
            end if
            if (t35(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = t35(1,i)
               phase(j) = t35(2,i)
            end if
            if (t45(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 4
               ampli(j) = t45(1,i)
               phase(j) = t45(2,i)
            end if
            if (t55(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 5
               ampli(j) = t55(1,i)
               phase(j) = t55(2,i)
            end if
            if (t65(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 6
               ampli(j) = t65(1,i)
               phase(j) = t65(2,i)
            end if
            write (itxt,780)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  780       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,2x,6(f8.3,f6.1,i2))
         end do
  790    continue
      end if
c
c     torsional angle parameters for 4-membered rings
c
      if (kt4(1) .ne. blank16) then
         write (itxt,800)
  800    format (//,15x,'4-Membered Ring Torsion Parameters',
     &           ///,17x,'Classes',15x,'KT Values',/)
         do i = 1, maxnt4
            if (kt4(i) .eq. blank16)  goto 820
            k1 = number(kt4(i)(1:4))
            k2 = number(kt4(i)(5:8))
            k3 = number(kt4(i)(9:12))
            k4 = number(kt4(i)(13:16))
            j = 0
            if (t14(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 1
               ampli(j) = t14(1,i)
               phase(j) = t14(2,i)
            end if
            if (t24(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 2
               ampli(j) = t24(1,i)
               phase(j) = t24(2,i)
            end if
            if (t34(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 3
               ampli(j) = t34(1,i)
               phase(j) = t34(2,i)
            end if
            if (t44(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 4
               ampli(j) = t44(1,i)
               phase(j) = t44(2,i)
            end if
            if (t54(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 5
               ampli(j) = t54(1,i)
               phase(j) = t54(2,i)
            end if
            if (t64(1,i) .ne. 0.0d0) then
               j = j + 1
               fold(j) = 6
               ampli(j) = t64(1,i)
               phase(j) = t64(2,i)
            end if
            write (itxt,810)  i,k1,k2,k3,k4,(ampli(k),
     &                        phase(k),fold(k),k=1,j)
  810       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,2x,6(f8.3,f6.1,i2))
         end do
  820    continue
      end if
c
c     pi-system torsion parameters
c
      if (kpt(1) .ne. blank8) then
         write (itxt,830)  formfeed,forcefield
  830    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,840)
  840    format (//,15x,'Pi-Orbital Torsion Parameters',
     &           ///,18x,'Classes',15x,'KPT',/)
         do i = 1, maxnpt
            if (kpt(i) .eq. blank8)  goto 860
            k1 = number(kpt(i)(1:4))
            k2 = number(kpt(i)(5:8))
            write (itxt,850)  i,k1,k2,ptcon(i)
  850       format (4x,i7,5x,i4,'-',i4,6x,f12.3)
         end do
  860    continue
      end if
c
c     stretch-torsion parameters
c
      if (kbt(1) .ne. blank16) then
         write (itxt,870)  formfeed,forcefield
  870    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,880)
  880    format (//,15x,'Stretch-Torsion Parameters',
     &           ///,17x,'Classes',12x,'Bond',8x,'KST1',
     &              8x,'KST2',8x,'KST3',/)
         do i = 1, maxnbt
            if (kbt(i) .eq. blank16)  goto 900
            k1 = number(kbt(i)(1:4))
            k2 = number(kbt(i)(5:8))
            k3 = number(kbt(i)(9:12))
            k4 = number(kbt(i)(13:16))
            write (itxt,890)  i,k1,k2,k3,k4,(btcon(j,i),j=1,9)
  890       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,9x,'1st',3f12.3,
     &              /,37x,'2nd',3f12.3,/,37x,'3rd',3f12.3)
         end do
  900    continue
      end if
c
c     angle-torsion parameters
c
      if (kat(1) .ne. blank16) then
         write (itxt,910)  formfeed,forcefield
  910    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,920)
  920    format (//,15x,'Angle-Torsion Parameters',
     &           ///,17x,'Classes',12x,'Angle',7x,'KAT1',
     &              8x,'KAT2',8x,'KAT3',/)
         do i = 1, maxnat
            if (kat(i) .eq. blank16)  goto 940
            k1 = number(kat(i)(1:4))
            k2 = number(kat(i)(5:8))
            k3 = number(kat(i)(9:12))
            k4 = number(kat(i)(13:16))
            write (itxt,930)  i,k1,k2,k3,k4,(atcon(j,i),j=1,6)
  930       format (2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,9x,'1st',3f12.3
     &              /,37x,'2nd',3f12.3)
         end do
  940    continue
      end if
c
c     torsion-torsion parameters
c
      if (ktt(1) .ne. blank20) then
         write (itxt,950)  formfeed,forcefield
  950    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,960)
  960    format (//,15x,'Torsion-Torsion Parameters',
     &           ///,19x,'Classes',18x,'KNX',9x,'KNY')
         do i = 1, maxntt
            if (ktt(i) .eq. blank20)  goto 990
            k1 = number(ktt(i)(1:4))
            k2 = number(ktt(i)(5:8))
            k3 = number(ktt(i)(9:12))
            k4 = number(ktt(i)(13:16))
            k5 = number(ktt(i)(17:20))
            write (itxt,970)  i,k1,k2,k3,k4,k5,tnx(i),tny(i)
  970       format (/,2x,i5,2x,i4,'-',i4,'-',i4,'-',i4,'-',i4,2x,2i12,/)
            k = tnx(i) * tny(i)
            write (itxt,980)  (tbf(j,i),j=1,k)
  980       format (3x,6f12.4)
         end do
  990    continue
      end if
c
c     van der Waals parameters
c
      exist = .false.
      do i = 1, maxtyp
         if (rad(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1000)  formfeed,forcefield
 1000    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         if (vdwindex .eq. 'CLASS') then
            write (itxt,1010)
 1010       format (//,15x,'Van der Waals Parameters',
     &              ///,21x,'Class',6x,'Radius',6x,'Epsilon',
     &                    4x,'Reduction',/)
         else
            write (itxt,1020)
 1020       format (//,15x,'Van der Waals Parameters',
     &              ///,22x,'Type',6x,'Radius',6x,'Epsilon',
     &                    4x,'Reduction',/)
         end if
         k = 0
         do i = 1, maxtyp
            if (rad(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1030)  k,i,rad(i),eps(i),reduct(i)
 1030          format (8x,i7,4x,i7,3f12.3)
            end if
         end do
c
c     van der Waals scaling parameters
c
         write (itxt,1040)  v2scale,v3scale,v4scale,v5scale
 1040    format (//,15x,'Van der Waals Scaling Factors',
     &           ///,20x,'1-2 Atoms',f17.3,/,20x,'1-3 Atoms',f17.3,
     &           /,20x,'1-4 Atoms',f17.3,/,20x,'1-5 Atoms',f17.3)
      end if
c
c     van der Waals 1-4 parameters for atom types
c
      exist = .false.
      do i = 1, maxtyp
         if (rad4(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         if (vdwindex .eq. 'CLASS') then
            write (itxt,1050)
 1050       format (//,15x,'Van der Waals Parameters for 1-4',
     &                 ' Interactions',
     &              ///,20x,'Class',7x,'Radius',6x,'Epsilon',/)
         else
            write (itxt,1060)
 1060       format (//,15x,'Van der Waals Parameters for 1-4',
     &                 ' Interactions',
     &              ///,20x,'Type',8x,'Radius',6x,'Epsilon',/)
         end if
         k = 0
         do i = 1, maxtyp
            if (rad4(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1070)  k,i,rad4(i),eps4(i)
 1070          format (8x,i7,2x,i7,2x,2f12.3)
            end if
         end do
      end if
c
c     van der Waals parameters for specific atom pairs
c
      if (kvpr(1) .ne. blank8) then
         if (vdwindex .eq. 'CLASS') then
            write (itxt,1080)
 1080       format (//,15x,'Van der Waals Parameters for Atom Pairs',
     &              ///,22x,'Classes',7x,'Radii Sum',4x,'Epsilon',/)
         else
            write (itxt,1090)
 1090       format (//,15x,'Van der Waals Parameters for Atom Pairs',
     &              ///,23x,'Types',8x,'Radii Sum',4x,'Epsilon',/)
         end if
         do i = 1, maxnvp
            if (kvpr(i) .eq. blank8)  goto 1110
            k1 = number(kvpr(i)(1:4))
            k2 = number(kvpr(i)(5:8))
            write (itxt,1100)  i,k1,k2,radpr(i),epspr(i)
 1100       format (8x,i7,5x,i4,'-',i4,2x,2f12.3)
         end do
 1110    continue
      end if
c
c     hydrogen bonding parameters for specific atom pairs
c
      if (khb(1) .ne. blank8) then
         if (vdwindex .eq. 'CLASS') then
            write (itxt,1120)
 1120       format (//,15x,'Hydrogen Bonding Parameters for Atom Pairs',
     &              ///,22x,'Classes',7x,'Radii Sum',4x,'Epsilon',/)
         else
            write (itxt,1130)
 1130       format (//,15x,'Hydrogen Bonding Parameters for Atom Pairs',
     &              ///,23x,'Types',8x,'Radii Sum',4x,'Epsilon',/)
         end if
         do i = 1, maxnhb
            if (khb(i) .eq. blank8)  goto 1150
            k1 = number(khb(i)(1:4))
            k2 = number(khb(i)(5:8))
            write (itxt,1140)  i,k1,k2,radhb(i),epshb(i)
 1140       format (8x,i7,5x,i4,'-',i4,2x,2f12.3)
         end do
 1150    continue
      end if
c
c     Pauli repulsion parameters
c
      exist = .false.
      do i = 1, maxclass
         if (prsiz(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1160)  formfeed,forcefield
 1160    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,1170)
 1170    format (//,15x,'Pauli Repulsion Parameters',
     &           ///,24x,'Class',14x,'Size',8x,'Damp',5x,'Valence'/)
         k = 0
         do i = 1, maxclass
            if (prsiz(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1180)  k,i,prsiz(i),prdmp(i),prele(i)
 1180          format (10x,i7,3x,i7,8x,2f12.4,f12.3)
            end if
         end do
      end if
c
c     damped dispersion parameters
c
      exist = .false.
      do i = 1, maxclass
         if (dspsix(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1190)  formfeed,forcefield
 1190    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,1200)
 1200    format (//,15x,'Damped Dispersion Parameters',
     &           ///,24x,'Class',15x,'C6',9x,'Damp',/)
         k = 0
         do i = 1, maxclass
            if (dspsix(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1210)  k,i,dspsix(i),dspdmp(i)
 1210          format (10x,i7,3x,i7,8x,2f12.4)
            end if
         end do
      end if
c
c     atomic partial charge parameters
c
      exist = .false.
      do i = 1, maxtyp
         if (chg(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1220)  formfeed,forcefield
 1220    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,1230)
 1230    format (//,15x,'Atomic Partial Charge Parameters',
     &           ///,24x,'Type',9x,'Partial Chg',/)
         k = 0
         do i = 1, maxtyp
            if (chg(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1240)  k,i,chg(i)
 1240          format (10x,i7,3x,i7,6x,f12.3)
            end if
         end do
c
c     atomic partial charge scaling parameters
c
         write (itxt,1250)  c1scale,c2scale,c3scale,c4scale,c5scale
 1250    format (//,15x,'Atomic Partial Charge Scaling Factors',
     &           ///,20x,'1-1 Atoms',f17.3,/,20x,'1-2 Atoms',f17.3,
     &           /,20x,'1-3 Atoms',f17.3,/,20x,'1-4 Atoms',f17.3,
     &           /,20x,'1-5 Atoms',f17.3)
      end if
c
c     bond dipole moment parameters
c
      if (kd(1) .ne. blank8) then
         write (itxt,1260)  formfeed,forcefield
 1260    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,1270)
 1270    format (//,15x,'Bond Dipole Moment Parameters',
     &           ///,25x,'Types',10x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd
            if (kd(i) .eq. blank8)  goto 1290
            k1 = number(kd(i)(1:4))
            k2 = number(kd(i)(5:8))
            write (itxt,1280)  i,k1,k2,dpl(i),pos(i)
 1280       format (10x,i7,5x,i4,'-',i4,6x,2f12.3)
         end do
 1290    continue
      end if
c
c     bond dipole moment parameters for 5-membered rings
c
      if (kd5(1) .ne. blank8) then
         write (itxt,1300)
 1300    format (//,15x,'5-Membered Ring Bond Dipole Parameters',
     &           ///,25x,'Types',10x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd5
            if (kd5(i) .eq. blank8)  goto 1320
            k1 = number(kd5(i)(1:4))
            k2 = number(kd5(i)(5:8))
            write (itxt,1310)  i,k1,k2,dpl5(i),pos5(i)
 1310       format (10x,i7,5x,i4,'-',i4,6x,2f12.3)
         end do
 1320    continue
      end if
c
c     bond dipole moment parameters for 4-membered rings
c
      if (kd4(1) .ne. blank8) then
         write (itxt,1330)
 1330    format (//,15x,'4-Membered Ring Bond Dipole Parameters',
     &           ///,25x,'Types',10x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd4
            if (kd4(i) .eq. blank8)  goto 1350
            k1 = number(kd4(i)(1:4))
            k2 = number(kd4(i)(5:8))
            write (itxt,1340)  i,k1,k2,dpl4(i),pos4(i)
 1340       format (10x,i7,5x,i4,'-',i4,6x,2f12.3)
         end do
 1350    continue
      end if
c
c     bond dipole moment parameters for 3-membered rings
c
      if (kd3(1) .ne. blank8) then
         write (itxt,1360)
 1360    format (//,15x,'3-Membered Ring Bond Dipole Parameters',
     &           ///,25x,'Types',10x,'Bond Dipole',4x,'Position',/)
         do i = 1, maxnd3
            if (kd3(i) .eq. blank8)  goto 1380
            k1 = number(kd3(i)(1:4))
            k2 = number(kd3(i)(5:8))
            write (itxt,1370)  i,k1,k2,dpl3(i),pos3(i)
 1370       format (10x,i7,5x,i4,'-',i4,6x,2f12.3)
         end do
 1380    continue
      end if
c
c     atomic multipole electrostatic parameters
c
      if (kmp(1) .ne. blank16) then
         write (itxt,1390)  formfeed,forcefield
 1390    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,1400)
 1400    format (//,17x,'Atomic Multipole Parameters',
     &           ///,11x,'Type',7x,'Axis Types',8x,'Frame',
     &              9x,'Multipoles (M-D-Q)',/)
         do i = 1, maxnmp
            if (kmp(i) .eq. blank16)  goto 1420
            k1 = number(kmp(i)(1:4))
            k2 = number(kmp(i)(5:8))
            k3 = number(kmp(i)(9:12))
            k4 = number(kmp(i)(13:16))
            write (itxt,1410)  i,k1,k2,k3,k4,mpaxis(i),multip(1,i),
     &                         multip(2,i),multip(3,i),multip(4,i),
     &                         multip(5,i),multip(8,i),multip(9,i),
     &                         multip(11,i),multip(12,i),multip(13,i)
 1410       format (2x,i5,3x,i4,3x,i4,2x,i4,2x,i4,5x,a8,2x,f10.5,
     &                 /,48x,3f10.5,/,48x,f10.5,
     &                 /,48x,2f10.5,/,48x,3f10.5)
         end do
 1420    continue
c
c     atomic multipole scaling parameters
c
         write (itxt,1430)  m2scale,m3scale,m4scale,m5scale
 1430    format (//,15x,'Atomic Multipole Scale Factors',
     &           ///,20x,'1-2 Atoms',f17.3,/,20x,'1-3 Atoms',f17.3,
     &           /,20x,'1-4 Atoms',f17.3,/,20x,'1-5 Atoms',f17.3)
      end if
c
c     charge penetration parameters
c
      exist = .false.
      do i = 1, maxclass
         if (cpele(i).ne.0.0d0 .or. cpalp(i).ne.0.0d0)  exist = .true
     &.
      end do
      if (exist) then
         write (itxt,1440)  formfeed,forcefield
 1440    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,1450)
 1450    format (//,15x,'Charge Penetration Parameters',
     &           ///,24x,'Class',10x,'Core Chg',8x,'Damp',/)
         k = 0
         do i = 1, maxclass
            if (cpele(i).ne.0.0d0 .or. cpalp(i).ne.0.0d0) then
               k = k + 1
               write (itxt,1460)  k,i,cpele(i),cpalp(i)
 1460          format (10x,i7,3x,i7,8x,2f12.4)
            end if
         end do
      end if
c
c     atomic dipole polarizability parameters
c
      exist = .false.
      use_thole = .false.
      use_tholed = .false.
      do i = 1, maxclass
         if (polr(i) .ne. 0.0d0)  exist = .true.
         if (athl(i) .ne. 0.0d0)  use_thole = .true.
         if (dthl(i) .ne. 0.0d0)  use_tholed = .true.
      end do
      if (exist) then
         write (itxt,1470)  formfeed,forcefield
 1470    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         if (use_tholed) then
            write (itxt,1480)
 1480       format (//,15x,'Dipole Polarizability Parameters',
     &              ///,22x,'Class',7x,'Alpha',5x,'Thole',
     &                 4x,'TholeD',6x,'Group Types',/)
         else if (use_thole) then
            write (itxt,1490)
 1490       format (//,15x,'Dipole Polarizability Parameters',
     &              ///,22x,'Class',7x,'Alpha',5x,'Thole',
     &                 6x,'Group Atom Types',/)
         else
            write (itxt,1500)
 1500       format (//,15x,'Dipole Polarizability Parameters',
     &              ///,22x,'Class',7x,'Alpha',6x,'Group Atom Types',/)
         end if
         k = 0
         do i = 1, maxclass
            if (polr(i) .ne. 0.0d0) then
               k = k + 1
               npg = 0
               do j = 1, maxval
                  if (pgrp(j,i) .ne. 0)  npg = npg + 1
               end do
               if (use_tholed) then
                  if (npg .eq. 0) then
                     write (itxt,1510)  k,i,polr(i),athl(i),dthl(i)
 1510                format (8x,i7,4x,i7,3x,3f10.3)
                  else
                     write (itxt,1520)  k,i,polr(i),athl(i),dthl(i),
     &                                  (pgrp(j,i),j=1,npg)
 1520                format (8x,i7,4x,i7,3x,3f10.3,4x,6i5)
                  end if
               else if (use_thole) then
                  if (npg .eq. 0) then
                     write (itxt,1530)  k,i,polr(i),athl(i)
 1530                format (8x,i7,4x,i7,3x,2f10.3)
                  else
                     write (itxt,1540)  k,i,polr(i),athl(i),
     &                                  (pgrp(j,i),j=1,npg)
 1540                format (8x,i7,4x,i7,3x,2f10.3,4x,6i5)
                  end if
               else
                  if (npg .eq. 0) then
                     write (itxt,1550)  k,i,polr(i)
 1550                format (8x,i7,4x,i7,3x,f10.3)
                  else
                     write (itxt,1560)  k,i,polr(i),(pgrp(j,i),j=1,npg)
 1560                format (8x,i7,4x,i7,3x,f10.3,4x,6i4)
                  end if
               end if
            end if
         end do
c
c     dipole polarizability scaling parameters
c
         write (itxt,1570)  d1scale,d2scale,d3scale,d4scale
 1570    format (//,15x,'Direct Induction Scale Factors',
     &           ///,20x,'1-1 Groups',f15.3,/,20x,'1-2 Groups',f15.3,
     &           /,20x,'1-3 Groups',f15.3,/,20x,'1-4 Groups',f15.3)
         write (itxt,1580)  u1scale,u2scale,u3scale,u4scale
 1580    format (//,15x,'Mutual Induction Scale Factors',
     &           ///,20x,'1-1 Groups',f15.3,/,20x,'1-2 Groups',f15.3,
     &           /,20x,'1-3 Groups',f15.3,/,20x,'1-4 Groups',f15.3)
         write (itxt,1590)  p2scale,p3scale,p4scale,p5scale
 1590    format (//,15x,'Inter-Group Polarizability Scale Factors',
     &           ///,20x,'1-2 Atoms',f16.3,/,20x,'1-3 Atoms',f16.3,
     &           /,20x,'1-4 Atoms',f16.3,/,20x,'1-5 Atoms',f16.3)
         write (itxt,1600)  p2iscale,p3iscale,p4iscale,p5iscale
 1600    format (//,15x,'Intra-Group Polarizability Scale Factors',
     &           ///,20x,'1-2 Atoms',f16.3,/,20x,'1-3 Atoms',f16.3,
     &           /,20x,'1-4 Atoms',f16.3,/,20x,'1-5 Atoms',f16.3)
         write (itxt,1610)  w2scale,w3scale,w4scale,w5scale
 1610    format (//,15x,'Induced Dipole Interaction Scale Factors',
     &           ///,20x,'1-2 Atoms',f16.3,/,20x,'1-3 Atoms',f16.3,
     &           /,20x,'1-4 Atoms',f16.3,/,20x,'1-5 Atoms',f16.3)
      end if
c
c     exchange polarization parameters
c
      exist = .false.
      do i = 1, maxclass
         if (pepdmp(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1620)  formfeed,forcefield
 1620    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,1630)
 1630    format (//,15x,'Exchange Polarization Parameters',
     &           ///,22x,'Class',8x,'Spring',8x,'Size',8x,'Damp',
     &              8x,'Use'/)
         k = 0
         do i = 1, maxclass
            if (pepdmp(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1640)  k,i,pepk(i),peppre(i),
     &                            pepdmp(i),pepl(k)
 1640          format (10x,i7,1x,i7,4x,2f12.4,f12.3,9x,l1)
            end if
         end do
      end if
c
c     charge transfer parameters
c
      exist = .false.
      do i = 1, maxclass
         if (ctchg(i).ne.0.0d0 .or. ctdmp(i).ne.0.0d0)  exist = .true
     &.
      end do
      if (exist) then
         write (itxt,1650)  formfeed,forcefield
 1650    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,1660)
 1660    format (//,15x,'Charge Transfer Parameters',
     &           ///,24x,'Class',12x,'Charge',7x,'Alpha',/)
         k = 0
         do i = 1, maxclass
            if (ctchg(i).ne.0.0d0 .or. ctdmp(i).ne.0.0d0) then
               k = k + 1
               write (itxt,1670)  k,i,ctchg(i),ctdmp(i)
 1670          format (10x,i7,3x,i7,8x,2f12.4)
            end if
         end do
      end if
c
c     bond charge flux parameters
c
      if (kcfb(1) .ne. blank8) then
         write (itxt,1680)  formfeed,forcefield
 1680    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,1690)
 1690    format (//,15x,'Bond Charge Flux Parameters',
     &           ///,22x,'Classes',14x,'KCFB',/)
         do i = 1, maxncfb
            if (kcfb(i) .eq. blank8)  goto 1710
            k1 = number(kcfb(i)(1:4))
            k2 = number(kcfb(i)(5:8))
            write (itxt,1700)  i,k1,k2,cflb(i)
 1700       format (8x,i7,5x,i4,'-',i4,6x,f12.4)
         end do
 1710    continue
      end if
c
c     angle charge flux parameters
c
      if (kcfa(1) .ne. blank12) then
         write (itxt,1720)  formfeed,forcefield
 1720    format (a1,//,15x,'Tinker Force Field Parameters for ',a20)
         write (itxt,1730)
 1730    format (//,15x,'Angle Charge Flux Parameters',
     &           ///,18x,'Classes',10x,'KCFA1',7x,'KCFA2',
     &              7x,'KCFB1',7x,'KCFB2',/)
         do i = 1, maxncfa
            if (kcfa(i) .eq. blank12)  goto 1750
            k1 = number(kcfa(i)(1:4))
            k2 = number(kcfa(i)(5:8))
            k3 = number(kcfa(i)(9:12))
            write (itxt,1740)  i,k1,k2,k3,cfla(1,i),cfla(2,i),
     &                        cflab(1,i),cflab(2,i)
 1740       format (1x,i7,5x,i4,'-',i4,'-',i4,1x,4f12.4)
         end do
 1750    continue
      end if
c
c     implicit solvation parameters
c
      exist = .false.
      do i = 1, maxtyp
         if (pbr(i).ne.0.0d0 .or. csr(i).ne.0.0d0
     &          .or. gkr(i).ne.0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1760)  formfeed,forcefield
 1760    format (a1,//,15x,'Tinker Force Field Parameters for ',a2
     &0)
         write (itxt,1770)
 1770    format (//,15x,'Implicit Solvation Parameters',
     &           ///,22x,'Type',6x,'PB Size',
     &              5x,'ddCOSMO',5x,'GK Size',/)
         k = 0
         do i = 1, maxtyp
            if (pbr(i).ne.0.0d0 .or. csr(i).ne.0.0d0
     &             .or. gkr(i).ne.0.0d0) then
               k = k + 1
               write (itxt,1780)  k,i,pbr(i),csr(i),gkr(i)
 1780          format (8x,i7,4x,i7,1x,3f12.4)
            end if
         end do
      end if
c
c     conjugated pisystem atom parameters
c
      exist = .false.
      do i = 1, maxclass
         if (ionize(i) .ne. 0.0d0)  exist = .true.
      end do
      if (exist) then
         write (itxt,1790)  formfeed,forcefield
 1790    format (a1,//,15x,'Tinker Force Field Parameters for ',a2
     &0)
         write (itxt,1800)
 1800    format (//,15x,'Conjugated Pisystem Atom Parameters',
     &           ///,20x,'Class',3x,'Electron',
     &              3x,'Ionization',3x,'Repulsion',/)
         k = 0
         do i = 1, maxclass
            if (ionize(i) .ne. 0.0d0) then
               k = k + 1
               write (itxt,1810)  k,i,electron(i),ionize(i),repulse(i)
 1810          format (6x,i7,4x,i7,f10.1,2x,2f12.3)
            end if
         end do
      end if
c
c     conjugated pisystem bond parameters
c
      if (kpi(1) .ne. blank8) then
         write (itxt,1820)
 1820    format (//,15x,'Conjugated Pisystem Bond Parameters',
     &           ///,20x,'Classes',8x,'d Force',4x,'d Length',/)
         do i = 1, maxnpi
            if (kpi(i) .eq. blank8)  goto 1840
            k1 = number(kpi(i)(1:4))
            k2 = number(kpi(i)(5:8))
            write (itxt,1830)  i,k1,k2,sslope(i),tslope(i)
 1830       format (6x,i7,5x,i4,'-',i4,3x,2f12.3)
         end do
 1840    continue
      end if
c
c     conjugated pisystem bond parameters for 5-membered rings
c
      if (kpi5(1) .ne. blank8) then
         write (itxt,1850)
 1850    format (//,15x,'5-Membered Ring Pisystem Bond Parameters'
     &,           ///,20x,'Classes',8x,'d Force',4x,'d Length',/)
         do i = 1, maxnpi5
            if (kpi5(i) .eq. blank8)  goto 1870
            k1 = number(kpi5(i)(1:4))
            k2 = number(kpi5(i)(5:8))
            write (itxt,1860)  i,k1,k2,sslope5(i),tslope5(i)
 1860       format (6x,i7,5x,i4,'-',i4,3x,2f12.3)
         end do
 1870    continue
      end if
c
c     conjugated pisystem bond parameters for 4-membered rings
c
      if (kpi4(1) .ne. blank8) then
         write (itxt,1880)
 1880    format (//,15x,'4-Membered Ring Pisystem Bond Parameters'
     &,           ///,20x,'Classes',8x,'d Force',4x,'d Length',/)
         do i = 1, maxnpi4
            if (kpi4(i) .eq. blank8)  goto 1900
            k1 = number(kpi4(i)(1:4))
            k2 = number(kpi4(i)(5:8))
            write (itxt,1890)  i,k1,k2,sslope4(i),tslope4(i)
 1890       format (6x,i7,5x,i4,'-',i4,3x,2f12.3)
         end do
 1900    continue
      end if
      return
      end

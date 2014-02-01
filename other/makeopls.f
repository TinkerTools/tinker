c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2009  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program makeopls  --  parse BOSS-format OPLS parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
      program makeopls
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'kchrge.i'
      integer i,j,k
      integer j1,j2,j3,j4
      integer k1,k2,k3,k4
      integer m1,m2,m3,m4
      integer it,ie,temp
      integer iboss,next
      integer ntyp,ncls
      integer nprm,digits
      integer attach
      integer freeunit
      integer typnum(maxtyp)
      integer clsnum(maxtyp)
      integer nlst(maxclass)
      integer clslst(20,maxclass)
      real*8 big,dummy
      real*8 mass,delta
      real*8 cg,rd,ep
      real*8 fc,bd,an
      real*8 vt1,vt2
      real*8 vt3,vt4
      real*8 radtyp(maxtyp)
      real*8 epstyp(maxtyp)
      real*8 prm1(maxprm)
      real*8 prm2(maxprm)
      real*8 prm3(maxprm)
      real*8 prm4(maxprm)
      logical exist,header,done
      character*3 id
      character*3 sm1,sm2
      character*3 sm3,sm4
      character*3 cls1(maxprm)
      character*3 cls2(maxprm)
      character*3 cls3(maxprm)
      character*3 cls4(maxprm)
      character*3 clsnam(maxclass)
      character*120 bossfile
      character*120 record
      character*120 string
c
c
c     setup for calculation and initialize types and classes
c
      call initial
c
c     get the name of the BOSS atom type definition file
c
      call nextarg (bossfile,exist)
      if (exist) then
         call basefile (bossfile)
         call suffix (bossfile,'atom')
         call version (bossfile,'old')
         inquire (file=bossfile,exist=exist)
      end if
      dowhile (.not. exist)
         write (iout,10)
   10    format (/,' Enter the BOSS Atom Definition File Name :  ',$)
         read (input,20)  bossfile
   20    format (a120)
         call basefile (bossfile)
         call suffix (bossfile,'atom')
         call version (bossfile,'old')
         inquire (file=bossfile,exist=exist)
      end do
c
c     open and then read the BOSS atom definition file;
c     the alternative "if" tests below choose between use
c     of atom types or classes for indexing vdw parameters
c
      ntyp = 0
      ncls = 0
      header = .true.
      iboss = freeunit ()
      open (unit=iboss,file=bossfile,status='old')
      rewind (unit=iboss)
      dowhile (.true.)
         read (iboss,30,err=80,end=80)  record
   30    format (a120)
         next = 1
         call getnumb (record,it,next)
         call getnumb (record,ie,next)
         call getword (record,id,next)
         string = record(next:120)
         read (string,*)  cg,rd,ep
         ntyp = ntyp + 1
         typnum(ntyp) = it
         atmnum(ntyp) = ie
         symbol(ntyp) = id
         chg(ntyp) = cg
         radtyp(ntyp) = rd
         epstyp(ntyp) = ep
         do i = 1, ntyp-1
            if (symbol(i) .eq. id) then
c           if (symbol(i).eq.id .and. radtyp(i).eq.rd
c    &                .and. epstyp(i).eq.ep) then
               clsnum(ntyp) = clsnum(i)
               goto 40
            end if
         end do
         ncls = ncls + 1
         clsnum(ntyp) = ncls
   40    continue
         if (header) then
            header = .false.
            write (iout,50)
   50       format (/,' Summary of OPLS Atom Definitions :')
            write (iout,60)
   60       format (/,6x,'OPLS',3x,'Type',2x,'Class',3x,'Atomic',
     &                 2x,'Symbol',2x,'Charge',4x,'Radius',6x,'Eps',/)
         end if
         write (iout,70)  it,ntyp,clsnum(ntyp),ie,id,cg,rd,ep
   70    format (3x,4i7,5x,a3,2f10.4,f12.6)
      end do
   80 continue
      close (unit=iboss)
c
c     find TINKER classes associated with each OPLS class
c
      ncls = 0
      header = .true.
      do i = 1, ntyp
         done = .false.
         do j = 1, i-1
            if (symbol(j) .eq. symbol(i))  done = .true.
         end do
         if (.not. done) then
            ncls = ncls + 1
            clsnam(ncls) = symbol(i)
            nlst(ncls) = 1
            clslst(nlst(ncls),ncls) = clsnum(i)
            do j = i+1, ntyp
               if (symbol(j) .eq. symbol(i)) then
                  done = .false.
                  do k = 1, nlst(ncls)
                     if (clsnum(j) .eq. clslst(k,ncls))  done = .true.
                  end do
                  if (.not. done) then
                     nlst(ncls) = nlst(ncls) + 1
                     if (nlst(ncls) .gt. 20) then
                        write (iout,90)
   90                   format (/,' MAKEOPLS  --  Too many TINKER',
     &                             ' Atom Classes in List')
                        call fatal
                     end if
                     clslst(nlst(ncls),ncls) = clsnum(j)
                  end if
               end if
            end do
            if (header) then
               header = .false.
               write (iout,100)
  100          format (/,' TINKER Classes for Each OPLS Class :',/)
            end if
            write (iout,110)  ncls,clsnam(ncls),
     &                        (clslst(j,ncls),j=1,nlst(ncls))
  110       format (6x,i5,3x,a3,3x,50i6)
         end if
      end do
c
c     print the atom type definitions in TINKER format
c
      write (iout,120)
  120 format (/,' Atom Type Definitions :',/)
      do i = 1, ntyp
         attach = 0
         mass = 0.0d0
         if (atmnum(i) .eq. 1)  mass = 1.008d0
         if (atmnum(i) .eq. 2)  mass = 4.003d0
         if (atmnum(i) .eq. 3)  mass = 6.941d0
         if (atmnum(i) .eq. 6)  mass = 12.011d0
         if (atmnum(i) .eq. 7)  mass = 14.007d0
         if (atmnum(i) .eq. 8)  mass = 15.999d0
         if (atmnum(i) .eq. 9)  mass = 18.998d0
         if (atmnum(i) .eq. 10)  mass = 20.179d0
         if (atmnum(i) .eq. 11)  mass = 22.990d0
         if (atmnum(i) .eq. 12)  mass = 24.305d0
         if (atmnum(i) .eq. 14)  mass = 28.086d0
         if (atmnum(i) .eq. 15)  mass = 30.974d0
         if (atmnum(i) .eq. 16)  mass = 32.06d0
         if (atmnum(i) .eq. 17)  mass = 35.453d0
         if (atmnum(i) .eq. 18)  mass = 39.948d0
         if (atmnum(i) .eq. 19)  mass = 39.098d0
         if (atmnum(i) .eq. 20)  mass = 40.08d0
         if (atmnum(i) .eq. 35)  mass = 79.904d0
         if (atmnum(i) .eq. 36)  mass = 83.80d0
         if (atmnum(i) .eq. 37)  mass = 85.468d0
         if (atmnum(i) .eq. 38)  mass = 87.62d0
         if (atmnum(i) .eq. 53)  mass = 126.905d0
         if (atmnum(i) .eq. 54)  mass = 131.30d0
         if (atmnum(i) .eq. 55)  mass = 132.905d0
         if (atmnum(i) .eq. 56)  mass = 137.33d0
         if (atmnum(i) .eq. 57)  mass = 138.906d0
         if (atmnum(i) .eq. 60)  mass = 144.24d0
         if (atmnum(i) .eq. 63)  mass = 151.96d0
         if (atmnum(i) .eq. 64)  mass = 157.25d0
         if (atmnum(i) .eq. 70)  mass = 173.04d0
         if (atmnum(i) .eq. 89)  mass = 227.0d0
         if (atmnum(i) .eq. 90)  mass = 232.038d0
         if (atmnum(i) .eq. 92)  mass = 238.029d0
         if (atmnum(i) .eq. 95)  mass = 243.0d0
         write (iout,130)  i,clsnum(i),symbol(i),atmnum(i),mass,attach
  130    format ('atom',6x,2i5,4x,a3,4x,'"',16x,'"',6x,i5,f11.3,i5)
      end do
c
c     print the van der Waals parameters in TINKER format;
c     the "if" test on "clsnum" is needed if atom classes,
c     instead of types, are used in index the vdw parameters
c
      write (iout,140)
  140 format (/,' Van der Waals Parameters :',/)
      delta = 0.00000001d0
      ncls = 0
      do i = 1, ntyp
c        if (clsnum(i) .gt. ncls) then
         ncls = ncls + 1
         rd = radtyp(i)
         ep = epstyp(i)
         digits = 6
         big = 100000.0d0
         if (abs(rd-dble(nint(big*rd))/big).lt.delta .and.
     &             abs(ep-dble(nint(big*ep))/big).lt.delta)  digits = 5
         big = 10000.0d0
         if (abs(rd-dble(nint(big*rd))/big).lt.delta .and.
     &             abs(ep-dble(nint(big*ep))/big).lt.delta)  digits = 4
         if (digits .eq. 4) then
            write (iout,150)  ncls,rd,ep
  150       format ('vdw',7x,i5,9x,2f11.4)
         else if (digits .eq. 5) then
            write (iout,160)  ncls,rd,ep
  160       format ('vdw',7x,i5,8x,2f12.5)
         else if (digits .eq. 6) then
            write (iout,170)  ncls,rd,ep
  170       format ('vdw',7x,i5,7x,2f13.6)
         end if
c        end if
      end do
c
c     print the partial charge parameters in TINKER format
c
      write (iout,180)
  180 format (/,' Partial Charge Parameters :',/)
      do i = 1, ntyp
         write (iout,190)  i,chg(i)
  190    format ('charge',4x,i5,3x,f11.4)
      end do
c
c     get the name of the BOSS bond parameter file
c
      call nextarg (bossfile,exist)
      if (exist) then
         call basefile (bossfile)
         call suffix (bossfile,'bond')
         call version (bossfile,'old')
         inquire (file=bossfile,exist=exist)
      end if
      dowhile (.not. exist)
         write (iout,200)
  200    format (/,' Enter the BOSS Bond Parameter File Name :  ',$)
         read (input,210)  bossfile
  210    format (a120)
         call basefile (bossfile)
         call suffix (bossfile,'bond')
         call version (bossfile,'old')
         inquire (file=bossfile,exist=exist)
      end do
c
c     open and then read the BOSS bond parameter file
c
      nprm = 0
      header = .true.
      iboss = freeunit ()
      open (unit=iboss,file=bossfile,status='old')
      rewind (unit=iboss)
      dowhile (.true.)
         sm1 = '   '
         sm2 = '   '
         fc = 0.0d0
         bd = 0.0d0
         read (iboss,220,err=270,end=270)  record
  220    format (a120)
         sm1 = record(1:2)
         sm2 = record(4:5)
         string = record(6:120)
         read (string,*,err=230,end=230)  fc,bd
  230    continue
         nprm = nprm + 1
         cls1(nprm) = sm1
         cls2(nprm) = sm2
         prm1(nprm) = fc
         prm2(nprm) = bd
         if (header) then
            header = .false.
            write (iout,240)
  240       format (/,' Summary of OPLS Bond Parameters :')
            write (iout,250)
  250       format (/,12x,'OPLS Classes',9x,'Force',6x,'Length',/)
         end if
         write (iout,260)  nprm,sm1,sm2,fc,bd
  260    format (4x,i5,5x,a3,3x,a3,3x,f12.1,f12.4)
      end do
  270 continue
      close (unit=iboss)
c
c     print the bond stretch parameters in TINKER format
c
      write (iout,280)
  280 format (/,' Bond Stretching Parameters :',/)
      do i = 1, nprm
         do j1 = 1, ncls
            if (clsnam(j1) .ne. cls1(i))  goto 310
            do j2 = 1, ncls
               if (clsnam(j2) .ne. cls2(i))  goto 300
               do k1 = 1, nlst(j1)
                  do k2 = 1, nlst(j2)
                     m1 = clslst(k1,j1)
                     m2 = clslst(k2,j2)
                     write (iout,290)  min(m1,m2),max(m1,m2),
     &                                 prm1(i),prm2(i)
  290                format ('bond',6x,2i5,4x,f11.1,f11.4)
                  end do
               end do
  300          continue
            end do
  310       continue
         end do
      end do
c
c     get the name of the BOSS angle parameter file
c
      call nextarg (bossfile,exist)
      if (exist) then
         call basefile (bossfile)
         call suffix (bossfile,'angle')
         call version (bossfile,'old')
         inquire (file=bossfile,exist=exist)
      end if
      dowhile (.not. exist)
         write (iout,320)
  320    format (/,' Enter the BOSS Angle Parameter File Name :  ',$)
         read (input,330)  bossfile
  330    format (a120)
         call basefile (bossfile)
         call suffix (bossfile,'angle')
         call version (bossfile,'old')
         inquire (file=bossfile,exist=exist)
      end do
c
c     open and then read the BOSS angle parameter file
c
      nprm = 0
      header = .true.
      iboss = freeunit ()
      open (unit=iboss,file=bossfile,status='old')
      rewind (unit=iboss)
      dowhile (.true.)
         sm1 = '   '
         sm2 = '   '
         sm3 = '   '
         fc = 0.0d0
         an = 0.0d0
         read (iboss,340,err=390,end=390)  record
  340    format (a120)
         sm1 = record(1:2)
         sm2 = record(4:5)
         sm3 = record(7:8)
         string = record(9:120)
         read (string,*,err=350,end=350)  fc,an
  350    continue
         nprm = nprm + 1
         cls1(nprm) = sm1
         cls2(nprm) = sm2
         cls3(nprm) = sm3
         prm1(nprm) = fc
         prm2(nprm) = an
         if (header) then
            header = .false.
            write (iout,360)
  360       format (/,' Summary of OPLS Angle Parameters :')
            write (iout,370)
  370       format (/,15x,'OPLS Classes',12x,'Force',7x,'Angle',/)
         end if
         write (iout,380)  nprm,sm1,sm2,sm3,fc,an
  380    format (4x,i5,5x,a3,3x,a3,3x,a3,3x,f12.2,f12.2)
      end do
  390 continue
      close (unit=iboss)
c
c     print the angle bend parameters in TINKER format
c
      write (iout,400)
  400 format (/,' Angle Bending Parameters :',/)
      do i = 1, nprm
         do j1 = 1, ncls
            if (clsnam(j1) .ne. cls1(i))  goto 440
            do j2 = 1, ncls
               if (clsnam(j2) .ne. cls2(i))  goto 430
               do j3 = 1, ncls
                  if (clsnam(j3) .ne. cls3(i))  goto 420
                  do k1 = 1, nlst(j1)
                     do k2 = 1, nlst(j2)
                        do k3 = 1, nlst(j3)
                           m1 = clslst(k1,j1)
                           m2 = clslst(k2,j2)
                           m3 = clslst(k3,j3)
                           write (iout,410)  min(m1,m3),m2,max(m1,m3),
     &                                       prm1(i),prm2(i)
  410                      format ('angle',5x,3i5,4x,f11.2,f11.2)
                        end do
                     end do
                  end do
  420             continue
               end do
  430          continue
            end do
  440       continue
         end do
      end do
c
c     get the name of the BOSS torsion parameter file
c
      call nextarg (bossfile,exist)
      if (exist) then
         call basefile (bossfile)
         call suffix (bossfile,'tors')
         call version (bossfile,'old')
         inquire (file=bossfile,exist=exist)
      end if
      dowhile (.not. exist)
         write (iout,450)
  450    format (/,' Enter the BOSS Torsion Parameter File Name :  ',$)
         read (input,460)  bossfile
  460    format (a120)
         call basefile (bossfile)
         call suffix (bossfile,'tors')
         call version (bossfile,'old')
         inquire (file=bossfile,exist=exist)
      end do
c
c     open and then read the BOSS torsion parameter file
c
      nprm = 0
      header = .true.
      iboss = freeunit ()
      open (unit=iboss,file=bossfile,status='old')
      rewind (unit=iboss)
      dowhile (.true.)
         sm1 = '   '
         sm2 = '   '
         sm3 = '   '
         sm4 = '   '
         vt1 = 0.0d0
         vt2 = 0.0d0
         vt3 = 0.0d0
         vt4 = 0.0d0
         read (iboss,470,err=520,end=520)  record
  470    format (a120)
         sm1 = record(48:49)
         sm2 = record(51:52)
         sm3 = record(54:55)
         sm4 = record(57:58)
         read (record,*,err=480,end=480)  dummy,vt1,vt2,vt3,vt4
  480    continue
         nprm = nprm + 1
         cls1(nprm) = sm1
         cls2(nprm) = sm2
         cls3(nprm) = sm3
         cls4(nprm) = sm4
         prm1(nprm) = vt1
         prm2(nprm) = vt2
         prm3(nprm) = vt3
         prm4(nprm) = vt4
         if (header) then
            header = .false.
            write (iout,490)
  490       format (/,' Summary of OPLS Torsion Parameters :')
            write (iout,500)
  500       format (/,10x,'OPLS Classes',15x,'V1',8x,'V2',
     &                 8x,'V3',8x,'V4',/)
         end if
         write (iout,510)  nprm,sm1,sm2,sm3,sm4,vt1,vt2,vt3,vt4
  510    format (4x,i5,6x,a3,3x,a3,3x,a3,3x,a3,3x,4f10.3)
      end do
  520 continue
      close (unit=iboss)
c
c     print the torsion parameters in TINKER format
c
      write (iout,530)
  530 format (/,' Torsional Parameters :',/)
      do i = 1, nprm
         do j1 = 1, ncls
            if (clsnam(j1) .ne. cls1(i))  goto 580
            do j2 = 1, ncls
               if (clsnam(j2) .ne. cls2(i))  goto 570
               do j3 = 1, ncls
                  if (clsnam(j3) .ne. cls3(i))  goto 560
                  do j4 = 1, ncls
                     if (clsnam(j4) .ne. cls4(i))  goto 550
                     do k1 = 1, nlst(j1)
                        do k2 = 1, nlst(j2)
                           do k3 = 1, nlst(j3)
                              do k4 = 1, nlst(j4)
                                 m1 = clslst(k1,j1)
                                 m2 = clslst(k2,j2)
                                 m3 = clslst(k3,j3)
                                 m4 = clslst(k4,j4)
                                 if ((m2.eq.m3 .and. m1.gt.m4) .or.
     &                                  m2.gt.m3) then
                                    temp = m1
                                    m1 = m4
                                    m4 = temp
                                    temp = m2
                                    m2 = m3
                                    m3 = temp
                                 end if
                                 if (prm4(i) .eq. 0.0d0) then
                                    write (iout,540)  m1,m2,m3,m4,
     &                                                prm1(i),prm2(i),
     &                                                prm3(i)
  540                               format ('torsion',3x,4i5,2x,f9.3,
     &                                      ' 0.0 1',f9.3,' 180.0 2',
     &                                      f9.3,' 0.0 3')
                                 else
                                    write (iout,545)  m1,m2,m3,m4,
     &                                                prm1(i),prm2(i),
     &                                                prm3(i),prm4(i)
  545                               format ('torsion',3x,4i5,2x,f9.3,
     &                                      ' 0.0 1',f9.3,' 180.0 2',
     &                                      f9.3,' 0.0 3',f9.3,
     &                                      ' 180.0 4')
                                 end if
                              end do
                           end do
                        end do
                     end do
  550                continue
                  end do
  560             continue
               end do
  570          continue
            end do
  580       continue
         end do
      end do
c
c     handle wildcard torsion parameters as a special case
c
      do i = 1, nprm
         if (cls1(i)(2:2).eq.'?' .and. cls4(i)(2:2).eq.'?') then
            do j2 = 1, ncls
               if (clsnam(j2) .ne. cls2(i))  goto 610
               do j3 = 1, ncls
                  if (clsnam(j3) .ne. cls3(i))  goto 600
                  do k2 = 1, nlst(j2)
                     do k3 = 1, nlst(j3)
                        m1 = 0
                        m2 = clslst(k2,j2)
                        m3 = clslst(k3,j3)
                        m4 = 0
                        if (m2 .gt. m3) then
                           temp = m1
                           m1 = m4
                           m4 = temp
                           temp = m2
                           m2 = m3
                           m3 = temp
                        end if
                        if (prm4(i) .eq. 0.0d0) then
                           write (iout,541)  m1,m2,m3,m4,prm1(i),
     &                                       prm2(i),prm3(i)
  541                      format ('torsion',3x,4i5,2x,f9.3,' 0.0 1',
     &                             f9.3,' 180.0 2',f9.3,' 0.0 3')
                        else
                           write (iout,546)  m1,m2,m3,m4,prm1(i),
     &                                       prm2(i),prm3(i),prm4(i)
  546                      format ('torsion',3x,4i5,2x,f9.3,' 0.0 1',
     &                             f9.3,' 180.0 2',f9.3,' 0.0 3',f9.3,
     &                             ' 180.0 4')
                        end if
                     end do
                  end do
  600             continue
               end do
  610          continue
            end do
         else if (cls1(i)(2:2) .eq. '?') then
            do j2 = 1, ncls
               if (clsnam(j2) .ne. cls2(i))  goto 650
               do j3 = 1, ncls
                  if (clsnam(j3) .ne. cls3(i))  goto 640
                  do j4 = 1, ncls
                     if (clsnam(j4) .ne. cls4(i))  goto 630
                     do k2 = 1, nlst(j2)
                        do k3 = 1, nlst(j3)
                           do k4 = 1, nlst(j4)
                              m1 = 0
                              m2 = clslst(k2,j2)
                              m3 = clslst(k3,j3)
                              m4 = clslst(k4,j4)
                              if ((m2.eq.m3 .and. m1.gt.m4) .or.
     &                               m2.gt.m3) then
                                 temp = m1
                                 m1 = m4
                                 m4 = temp
                                 temp = m2
                                 m2 = m3
                                 m3 = temp
                              end if
                              if (prm4(i) .eq. 0.0d0) then
                                 write (iout,542)  m1,m2,m3,m4,prm1(i),
     &                                             prm2(i),prm3(i)
  542                            format ('torsion',3x,4i5,2x,f9.3,
     &                                   ' 0.0 1',f9.3,' 180.0 2',
     &                                   f9.3,' 0.0 3')
                              else
                                 write (iout,547)  m1,m2,m3,m4,prm1(i),
     &                                             prm2(i),prm3(i),
     &                                             prm4(i)
  547                            format ('torsion',3x,4i5,2x,f9.3,
     &                                   ' 0.0 1',f9.3,' 180.0 2',
     &                                   f9.3,' 0.0 3',f9.3,' 180.0 4')
                              end if
                           end do
                        end do
                     end do
  630                continue
                  end do
  640             continue
               end do
  650          continue
            end do
         else if (cls4(i)(2:2) .eq. '?') then
            do j1 = 1, ncls
               if (clsnam(j1) .ne. cls1(i))  goto 690
               do j2 = 1, ncls
                  if (clsnam(j2) .ne. cls2(i))  goto 680
                  do j3 = 1, ncls
                     if (clsnam(j3) .ne. cls3(i))  goto 670
                     do k1 = 1, nlst(j1)
                        do k2 = 1, nlst(j2)
                           do k3 = 1, nlst(j3)
                              m1 = clslst(k1,j1)
                              m2 = clslst(k2,j2)
                              m3 = clslst(k3,j3)
                              m4 = 0
                              if ((m2.eq.m3 .and. m1.gt.m4) .or.
     &                               m2.gt.m3) then
                                 temp = m1
                                 m1 = m4
                                 m4 = temp
                                 temp = m2
                                 m2 = m3
                                 m3 = temp
                              end if
                              if (prm4(i) .eq. 0.0d0) then
                                 write (iout,543)  m1,m2,m3,m4,prm1(i),
     &                                             prm2(i),prm3(i)
  543                            format ('torsion',3x,4i5,2x,f9.3,
     &                                   ' 0.0 1',f9.3,' 180.0 2',
     &                                   f9.3,' 0.0 3')
                              else
                                 write (iout,548)  m1,m2,m3,m4,prm1(i),
     &                                             prm2(i),prm3(i),
     &                                             prm4(i)
  548                            format ('torsion',3x,4i5,2x,f9.3,
     &                                   ' 0.0 1',f9.3,' 180.0 2',
     &                                   f9.3,' 0.0 3',f9.3,' 180.0 4')
                              end if
                           end do
                        end do
                     end do
  670                continue
                  end do
  680             continue
               end do
  690          continue
            end do
         end if
      end do
      end

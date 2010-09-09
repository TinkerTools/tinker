c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine readgau  --  read data from G09 output file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "readgau" reads an ab initio optimized structure, forces,
c     Hessian and frequencies from a Gaussian 09 output file
c
c
      subroutine readgau
      implicit none
      include 'sizes.i'
      include 'ascii.i'
      include 'iounit.i'
      include 'qmstuf.i'
      include 'units.i'
      integer i,j
      integer igau,code
      integer nfreq,nghess
      integer itmp,jtmp,ktmp
      integer length,next
      integer freeunit
      integer trimtext
      logical hasinputxyz
      logical hasMP2
      logical waiter
      logical exist
      real*8 frcunit,hessunit
      character*4 arcstart
      character*6 gname
      character*120 keyword
      character*120 gaufile
      character*120 record
      character*120 string
      character*120 word
c
c
c     initialize some values prior to opening the log file
c
      exist = .false.
      hasinputxyz = .false.
      ngatom = 0
      nfreq = 0
      arcstart = '1'//char(backslash)//'1'//char(backslash)
      waiter = .false.
c
c     specify and open the Gaussian 09 output log file
c
      call nextarg (gaufile,exist)
      if (exist) then
         inquire (file=gaufile,exist=exist)
         igau = freeunit()
         call basefile (gaufile)
         call suffix (gaufile,'log')
         call version (gaufile,'old')
         inquire (file=gaufile,exist=exist)
         if (.not. exist) then
            call basefile (gaufile)
            call suffix (gaufile,'out')
            call version (gaufile,'old')
            inquire (file=gaufile,exist=exist)
         end if
      end if
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter the Name of the Gaussian Output File :  ',$)
         read (input,20)  gaufile
   20    format (a120)
         igau = freeunit ()
         call basefile (gaufile)
         call suffix (gaufile,'log')
         call version (gaufile,'old')
         inquire (file=gaufile,exist=exist)
         if (.not. exist) then
            call basefile (gaufile)
            call suffix (gaufile,'out')
            call version (gaufile,'old')
            inquire (file=gaufile,exist=exist)
         end if
      end do
c
c     read structure, forces and frequencies from Gaussian output
c
      open (unit=igau,file=gaufile,status='old')
      rewind (unit=igau)
c     do while (.true. .and. .not.eof(igau))
      do while (.true.)
         read (igau,30,err=160,end=160)  record
   30    format (a120)
         next = 1
         keyword = record
         call trimhead (keyword)
         length = trimtext (keyword)
         call upcase (keyword)
         if (waiter.and.keyword (1:12) .ne. 'JOB CPU TIME') then
            goto 160
         else
            waiter = .false.
         end if
         if (keyword(1:20) .eq. 'STANDARD ORIENTATION') then
            do i = 1, 4
               read (igau,40,err=160,end=160)  record
   40          format (a120)
            end do
            i = 1
            do while (.true.)
               read (igau,50,err=160,end=160)  record
   50          format (a120)
               read (record,*,err=60,end=60)  itmp,jtmp,ktmp,
     &                                        gx(i),gy(i),gz(i)
               if (jtmp .le. 0)  goto 60
               i = i + 1
            end do
   60       continue
            ngatom = i - 1
         else if (keyword(37:58) .eq. 'FORCES (HARTREES/BOHR)') then
            read (igau,70,err=160,end=160)  record
   70       format (a120)
            read (igau,80,err=160,end=160)  record
   80       format (a120)
            frcunit = hartree / bohr
            do i = 1, ngatom
               read (igau,90,err=160,end=160)  record
   90          format (a120)
               read (record,*,err=100,end=100)  itmp,jtmp,gforce(1,i),
     &                                          gforce(2,i),gforce(3,i)
               do j = 1, 3
                  gforce(j,i) = frcunit * gforce(j,i)
               end do
  100          continue
            end do
         else if (keyword(1:14) .eq. 'FREQUENCIES --') then
            string = keyword(15:120)
            read (string,*,err=110,end=110)  gfreq(nfreq+1),
     &                                       gfreq(nfreq+2),
     &                                       gfreq(nfreq+3)
  110       continue
            nfreq = nfreq + 3
c
c     read the Hessian from archive section at bottom of output
c
         else if (keyword(1:4) .eq. arcstart) then
            itmp = 0
c           do while (.true. .and. .not.eof(igau))
            do while (.true.)
               if (next .gt. 73) then
                  read (igau,120,err=160,end=160)  record
  120             format (a120)
                  next = 1
               end if
               call readarcword (igau,record,word,length,next)
               if (word(1:1) .eq. char(backslash))  itmp = itmp + 1
               if (itmp.eq.16 .and. hasinputxyz) then
                  do i = 1, ngatom
                     do j = 1, 5
                        if (next .gt. 73) then
                           read (igau,130,err=160,end=160)  record
  130                      format (a120)
                           next = 1
                        end if
                        call readarcword (igau,record,word,length,next)
                        if (j .eq. 1)  read(word(1:length),*)  gname
                        if (j .eq. 2)  read(word(1:length),*)  gx(i)
                        if (j .eq. 3)  read(word(1:length),*)  gy(i)
                        if (j .eq. 4)  read(word(1:length),*)  gz(i)
                     end do
                  end do
               end if
               if (itmp.gt.16 .and. word(1:2).eq.'HF') then
                  do i = 1, 2
                     if (next .gt. 73) then
                        read (igau,140,err=160,end=160)  record
  140                   format (a120)
                        next = 1
                     end if
                     call readarcword (igau,record,word,length,next)
                  end do
                  read (word(1:length),*)  egau
                  egau = hartree * egau
               else if (itmp.gt.16 .and. word(1:3).eq.'MP2') then
                  hasmp2 = .true.
                  do i = 1, 2
                     if (next .gt. 73) then
                        read (igau,150,err=160,end=160)  record
  150                   format (a120)
                        next = 1
                     end if
                     call readarcword (igau,record,word,length,next)
                  end do
                  read (word(1:length),*)  egau
                  egau = hartree * egau
               else if (word(1:5) .eq. 'NImag') then
                  do i = 1, 4
                     call readarcword (igau,record,word,length,next)
                  end do
                  hessunit = hartree / bohr**2
                  nghess = (3*ngatom*(3*ngatom+1)) / 2
                  do i = 1, nghess
                     call readarcword (igau,record,word,length,next)
                     read (word(1:length),*)  gh(i)
                     gh(i) = hessunit * gh(i)
                  end do
                  goto 160
               end if
               code = ichar(word(1:1))
               if (code .eq. atsign)  then
                  waiter = .true.
                  goto 160
               end if
            end do
         end if
      end do
  160 continue
      close (unit=igau)
c
c     zero out the frequencies if none were in Gaussian output
c
      if (nfreq .eq. 0) then
         do i = 1, 3*ngatom
            gfreq(i) = 0.0d0
         end do
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine readarcword  --  read Gaussian archive section  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "readarcword" reads data from Gaussian archive section; each
c     entry is terminated with a backslash symbol
c
c     igau     file unit of the Gaussian output file
c     word     information to be read
c     length   length of the word
c
c
      subroutine readarcword (igau,string,word,length,next)
      implicit none
      include 'ascii.i'
      integer i,igau,code
      integer next,length
      character*1 letter
      character*120 word
      character*120 string
c
c
c     initialize some values prior to parsing the test string
c
      length = 1
      letter = ' '
      do i = 1, 120
         word(i:i) = ' '
      end do
c
c     attempt to read a text word entry from the input string
c
      letter = string (next:next)
      code = ichar(letter)
      if (code.eq.backslash .or. code.eq.equal
     &       .or. code.eq.space) then
         word(1:1) = letter
         next = next + 1
         length = 1
         return
      end if
   10 continue
      do i = next, 75
         if (code.eq.backslash .or. code.eq.equal
     &          .or. code.eq.space)  return
         if (next .gt. 70) then
            read (igau,20,err=30,end=30)  string
   20       format (a120)
            next = 1
            goto 10
         end if
         if (code .eq. comma) then
            next = next + 1
            return
         end if
         if (code.eq.backslash .or. code.eq.equal
     &          .or. code.eq.space)  return
         word(length:length) = letter
         next = next + 1
         letter = string(next:next)
         code = ichar(letter)
         length = length + 1
      end do
      if (code .eq. atsign) then
         word(1:1) = letter
         length = 1
      end if
   30 continue
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine trimhead  --  remove spaces before first text  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "trimhead" removes the spaces before the first non-blank
c     character in a text string
c
c
      subroutine trimhead (string)
      implicit none
      integer i,j,k
      character*120 string
      character*120 temp
c
c
c     loop over characters, removing blank beginning spaces
c
      do i = 1, 120
         temp(i:i) = ' '
      end do
      j = 0
      k = 0
      do i = 1, 120
         if (string(i:i) .ne. ' ')  j = 1
         if (j .eq. 1) then
            k = k + 1
            temp(k:k) = string(i:i)
         end if
      end do
      do i = 1, 120
         string(i:i) = temp(i:i)
      end do
      return
      end

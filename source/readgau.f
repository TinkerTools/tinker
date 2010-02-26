c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine readgau  --  read data from G03 output file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "readgau" reads an ab initio optimized structure, forces,
c     Hessian and frequencies from a Gaussian 03 output file
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
      logical exist
      real*8 frcunit,hessunit
      character*4 arcstart
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
      ngatom = 0
      nfreq = 0
      arcstart = '1'//char(backslash)//'1'//char(backslash)
c
c     specify and open the Gaussian 03 output log file
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
      dowhile (.not. exist)
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
      open (unit=igau,file=gaufile,status='old')
      rewind (unit=igau)
c
c     read each line of Gaussian output and find nonblank string
c
      dowhile (.true.)
         read (igau,30,err=120,end=120)  keyword
   30    format (a120)
         j = 0
         do i = 120, 1, -1
            if (keyword(i:i) .ne. ' ')  j = i
         end do
         j = j - 1
         do i = 1, 120-j
            keyword(i:i) = keyword(i+j:i+j)
         end do
         do i = 121-j, 120
            keyword(i:i) = ' '
         end do
         length = trimtext (keyword)
         call upcase (keyword)
c
c     get structure, forces and frequencies from Gaussian output
c
         if (keyword(1:20). eq. 'STANDARD ORIENTATION') then
            do i = 1, 4
               read (igau,30,err=120,end=120)  record
            end do
            i = 1
            dowhile (.true.)
               read (igau,40,err=120,end=120)  record
   40          format (a120)
               read (record,*,err=50,end=50)  itmp,jtmp,ktmp,
     &                                        gx(i),gy(i),gz(i)
               if (jtmp .le. 0)  goto 50
               i = i + 1
            end do
   50       continue
            ngatom = i - 1
         else if (keyword(37:58) .eq. 'FORCES (HARTREES/BOHR)') then
            frcunit = hartree / bohr
            read (igau,60,err=120,end=120)  record
   60       format (a120)
            read (igau,70,err=120,end=120)  record
   70       format (a120)
            do i = 1, ngatom
               read (igau,80,err=120,end=120)  record
   80          format (a120)
               read (record,*,err=90,end=90)  itmp,jtmp,gforce(1,i),
     &                                          gforce(2,i),gforce(3,i)
               do j = 1, 3
                  gforce(j,i) = gforce(j,i) * frcunit
               end do
   90          continue
            end do
         else if (keyword(1:14) .eq. 'FREQUENCIES --') then
            string = keyword(15:120)
            read (string,*,err=100,end=100)  gfreq(nfreq+1),
     &                                       gfreq(nfreq+2),
     &                                       gfreq(nfreq+3)
  100       continue
            nfreq = nfreq + 3
c
c     get the Hessian from archive section at bottom of output
c
         else if (keyword(1:4) .eq. arcstart) then
            hessunit = hartree / bohr**2
            next = 1
            dowhile (.true.)
               call readarcword (igau,record,word,length,next)
               if (word(1:5) .eq. 'NImag') then
                  do i = 1, 4
                     call readarcword (igau,record,word,length,next)
                  end do
                  nghess = (ngatom*3*(ngatom*3+1)) / 2
                  do i = 1, nghess
                     call readarcword (igau,record,word,length,next)
                     read (word(1:length),*)  gh(i)
                     gh(i) = gh(i) * hessunit
                  end do
                  goto 110
               end if
               code = ichar(word(1:1))
               if (code .eq. atsign)  goto 110
            end do
  110       continue
         end if
      end do
  120 continue
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
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getarcword  --  read Gaussian archive section  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getarcword" reads data from Gaussian archive section; each
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
      if (code.eq.backslash .or. code.eq.equal) then
         word(1:1) = letter
         next = next + 1
         length = 1
         return
      end if
   10 continue
      do i = next, 75
         if (code.eq.backslash .or. code.eq.equal)  return
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
         if (code.eq.backslash .or. code.eq.equal)  return
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

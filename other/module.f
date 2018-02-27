c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2014  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program module  --  convert include common to use module  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "module" takes a Fortran 77 source file with include statements
c     accessing common blocks and outputs a Fortran 90 source file
c     with use statements accessing modules
c
c
      program module
      implicit none
      integer i,j,isrc
      integer nline
      integer maxline
      integer length
      integer trimtext
      integer freeunit
      integer start,stop
      logical exist
      character*240 srcfile
      character*240 record
      character*240, allocatable :: line(:)
c
c
c     setup the use of the Tinker subroutine library
c
      call initial
c
c     try to get a filename from the command line arguments
c
      call nextarg (srcfile,exist)
      if (exist) then
         call basefile (srcfile)
         call suffix (srcfile,'f','old')
         inquire (file=srcfile,exist=exist)
      end if
c
c     ask for the user specified input source file name
c
      do while (.not. exist)
         write (*,10)
   10    format (/,' Enter Fortran Source File Name :  ',$)
         read (*,20)  srcfile
   20    format (a240)
         call basefile (srcfile)
         call suffix (srcfile,'f','old')
         inquire (file=srcfile,exist=exist)
      end do
c
c     perform dynamic allocation of some local arrays
c
      maxline = 50000
      allocate (line(maxline))
c
c     first open and then read the Fortran source file
c
      isrc = freeunit ()
      open (unit=isrc,file=srcfile,status='old')
      rewind (unit=isrc)
      nline = 0
      dowhile (.true.)
         read (isrc,30,err=40,end=40)  record
   30    format (a240)
         nline = nline + 1
         line(nline) = record
      end do
   40 continue
      close (unit=isrc,status='delete')
c
c     replace included common blocks with used modules
c
      do i = 1, nline
         record = line(i)
         if (record(1:19) .eq. '      implicit none') then
            start = i
            do j = i+1, nline
               record = line(j)
               if (record(1:13) .ne. '      include') then
                  stop = j - 1
                  goto 50
               end if
            end do
   50       continue
            do j = start, stop-1
               record = line(j+1)
               length = trimtext (record)
               length = length - 3
               line(j) = '      use '//record(16:length)
            end do
            line(stop) = '      implicit none'
         end if
      end do
c
c     example to convert one module name into another
c
      do i = 1, nline
         record = line(i)
         if (record(1:16) .eq. '      use atmtyp') then
            length = trimtext (record)
            line(i) = '      use atomid'//record(17:length)
         end if
         if (record(1:15) .eq. '      use angle') then
            length = trimtext (record)
            line(i) = '      use angbnd'//record(17:length)
         end if
         if (record(1:14) .eq. '      use bond') then
            length = trimtext (record)
            line(i) = '      use bndstr'//record(17:length)
         end if
      end do
c
c     open the output and copy modified source to a file
c
      isrc = freeunit ()
      call suffix (srcfile,'f','new')
      open (unit=isrc,file=srcfile,status='new')
      rewind (unit=isrc)
      do i = 1, nline
         record = line(i)
         length = trimtext (record)
         write (isrc,60)  record(1:length)
   60    format (a)
      end do
      close (unit=isrc)
c
c     perform deallocation of some local arrays
c
      deallocate (line)
      end

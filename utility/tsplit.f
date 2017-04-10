c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program tsplit  --  split source code into separate files  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "tsplit" takes a concatenated source listing of the TINKER code
c     and splits it into individual source files; a TINKER-specific
c     version of the standard fsplit routine
c
c
      program tsplit
      implicit none
      integer maxline
      parameter (maxline=500000)
      integer i,j,k,m
      integer input,iout
      integer itxt,leng
      integer start,stop
      integer nline,module
      character*240 filename
      character*240 record
      character*240 code(maxline)
c
c
c     get the name of the source file to split into modules
c
      input = 5
      iout = 6
      write (iout,10)
   10 format (/,' Enter TINKER Source Listing File Name :  ',$)
      read (input,20)  filename
   20 format (a240)
c
c     read and store a concatenated listing of TINKER source
c
      itxt = 1
c     itxt = freeunit ()
      open (unit=itxt,file=filename,status='old')
      do i = 1, maxline
         read (itxt,30,err=40,end=40)  code(i)
   30    format (a240)
         nline = i
      end do
   40 continue
      close (unit=itxt)
c
c     search for and write out the individual source modules
c
      i = 0
      module = 0
      dowhile (i .le. nline)
         i = i + 1
         record = code(i)
         if (index(record,'COPYRIGHT') .ne. 0) then
            stop = i - 4
            if (module .gt. 0) then
               open (unit=itxt,name=filename,status='new')
               do j = start, stop
                  record = code(j)
                  do m = 80, 1, -1
                     if (record(m:m) .ne. ' ') then
                        leng = m
                        goto 50
                     end if
                  end do
   50             continue
                  write (itxt,60)  record(1:leng)
   60             format (a)
               end do
               close (unit=itxt)
            end if
            module = module + 1
            start = i - 3
            do k = 5, 50
               record = code(i+k)
               if (index(record,'##  program') .ne. 0) then
                  do m = 19, 80
                     if (record(m:m) .eq. ' ') then
                        filename = record(19:m-1)//'.f'
                        goto 70
                     end if
                  end do
               else if (index(record,'##  subroutine') .ne. 0) then
                  do m = 22, 80
                     if (record(m:m) .eq. ' ') then
                        filename = record(22:m-1)//'.f'
                        goto 70
                     end if
                  end do
               else if (index(record,'##  function') .ne. 0) then
                  do m = 20, 80
                     if (record(m:m) .eq. ' ') then
                        filename = record(20:m-1)//'.f'
                        goto 70
                     end if
                  end do
               else if (index(record,'.i') .ne. 0) then
                  m = index(record,'.i') - 1
                  filename = record(11:m)//'.i'
                  goto 70
               end if
            end do
   70       continue
            i = i + k
         end if
      end do
      stop = nline
      open (unit=itxt,name=filename,status='new')
      do j = start, stop
         record = code(j)
         do m = 80, 1, -1
            if (record(m:m) .ne. ' ') then
               leng = m
               goto 80
            end if
         end do
   80    continue
         write (itxt,90)  record(1:leng)
   90    format (a)
      end do
      close (unit=itxt)
      end

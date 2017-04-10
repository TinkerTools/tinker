c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program header  --  change symbol used in file headers  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "header" alters a source file to use "#" instead of "*"
c     as the decorative symbol in the file headers
c
c
      program header
      implicit none
      include 'argue.i'
      include 'files.i'
      include 'iounit.i'
      integer maxline
      parameter (maxline=100000)
      integer i,k,isrc,itop,nline
      integer length,trimtext,freeunit
      character*240 string
      character*240 record(maxline)
c
c
c     get a filename supplied as a command line argument
c
      call initial
      if (narg .ge. 1) then
         filename = arg(1)
c
c     ask for the user specified input source filename
c
      else
   10    continue
         write (iout,20)
   20    format (' Enter the Name of the Source File :  ',$)
         read (input,30)  filename
   30    format (a240)
      end if
      isrc = freeunit ()
      open (unit=isrc,file=filename,status='old',err=10)
      rewind (unit=isrc)
c
c     change the header symbol in the first few lines
c
      nline = 0
      do i = 1, maxline
         read (isrc,40,err=50,end=50)  record(i)
   40    format (a240)
         nline = nline + 1
      end do
c
c     close the file with deletion, then reopen it for writing
c
   50 continue
      close (unit=isrc,status='delete')
      open (unit=isrc,file=filename,status='new')
c
c     write out the altered source file with a new header
c
      itop = 0
      do i = 1, nline
         itop = itop + 1
         string = record(i)
         length = trimtext (string)
         if (itop .le. 11) then
            do k = 1, length
               if (string(k:k) .eq. '*')  string(k:k) = '#'
            end do
         end if
         if (length .eq. 9) then
            if (string(7:9) .eq. 'end')  itop = 0
         end if
         write (isrc,60)  (string(k:k),k=1,length)
   60    format (<length>a1)
      end do
      close (unit=isrc)
c
c     perform any final tasks before program exit
c
      call final
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  program reorder  --  put scan minima in energy order  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "reorder" takes the coordinate file from a potential energy
c     surface scan and writes a new set of files numbered in the
c     order of energy rank
c
c     the list file required as input contains an original file
c     number and its energy on each line; the file should be
c     sorted by energy value since new files are written in the
c     order they appear in the list file
c
c     note the original scan coordinate files will be deleted!
c
c
      program reorder
      implicit none
      include 'iounit.i'
      include 'titles.i'
      integer maxmin
      parameter (maxmin=100000)
      integer i,nmin,lext,next
      integer ilist,iold,inew
      integer trimtext,freeunit
      integer imin(maxmin)
      character*3 oldname,newname
      character*7 ext
      character*9 emin(maxmin)
      character*240 listfile
      character*240 oldfile
      character*240 newfile
      character*240 record
c
c
c     get the filename for the energy ordered list of minima
c
      call initial
      write (iout,10)
   10 format (/,' Enter Name of File Listing Sorted Minima :  ',$)
      read (input,20)  listfile
   20 format (a240)
      call basefile (listfile)
      call suffix (listfile,'list')
      call version (listfile,'old')
c
c     read the list storing the energy ranks and values
c
      ilist = freeunit ()
      open (unit=ilist,file=listfile,status='old')
      rewind (unit=ilist)
      nmin = 0
      dowhile (.true.)
         read (ilist,30,err=40,end=40)  record
   30    format (a240)
         nmin = nmin + 1
         next = 1
         call getnumb (record,imin(nmin),next)
         call gettext (record,emin(nmin),next)
      end do
   40 continue
      close (unit=ilist)
c
c     set the old and new base filenames; change as needed
c
      oldname = 'map'
      newname = 'chd'
c
c     write new coordinate files numbered in rank order
c
      write (iout,50)
   50 format (/,' Minimum Reordering following Potential',
     &           ' Surface Scan :',/)
      do i = 1, nmin
         lext = 3
         call numeral (imin(i),ext,lext)
         oldfile = oldname//'.'//ext(1:lext)
         iold = freeunit ()
         open (unit=iold,file=oldfile,status='old')
         rewind (unit=iold)
         call readxyz (iold)
         close (iold,status='delete')
         lext = 3
         call numeral (i,ext,lext)
         newfile = newname//'.'//ext(1:lext)
         title = 'Cycloheptadecane  (Minimum '//ext(1:lext)//
     &              ', Energy '//emin(i)//')'
         ltitle = trimtext (title)
         inew = freeunit ()
         open (unit=inew,file=newfile,status='new')
         call prtxyz (inew)
         close (inew)
         if (mod(i,100) .eq. 0) then
            write (iout,60)  i
   60       format (' Processing Complete for File',i8)
         end if
      end do
      end

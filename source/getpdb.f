c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine getpdb  --  get a Protein Data Bank file  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "getpdb" asks for a Protein Data Bank file name, gets the
c     format type, and then reads in the coordinates file
c
c
      subroutine getpdb
      use files
      use inform
      use iounit
      use pdb
      implicit none
      integer ipdb,nask
      integer next,last,temp
      integer freeunit
      logical exist,done
      real*8 crd
      character*1 letter1
      character*3 letter3
      character*4 letter4
      character*6 remark
      character*240 pdbfile
      character*240 record
      character*240 string
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (pdbfile,exist)
      if (exist) then
         call basefile (pdbfile)
         call suffix (pdbfile,'pdb','old')
         inquire (file=pdbfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      nask = 0
      do while (.not.exist .and. nask.lt.maxask)
         nask = nask + 1
         pdbfile = ' '
         write (iout,10)
   10    format (/,' Enter Protein Data Bank File Name :  ',$)
         read (input,20)  pdbfile
   20    format (a240)
         call basefile (pdbfile)
         call suffix (pdbfile,'pdb','old')
         inquire (file=pdbfile,exist=exist)
      end do
      if (.not. exist)  call fatal
c
c     open the coordinates file with PDB format as default
c
      pdbtyp = 'PDB'
      filename = pdbfile
      ipdb = freeunit ()
      open (unit=ipdb,file=pdbfile,status='old')
      rewind (unit=ipdb)
c
c     check format by trying to read coordinates in PDB format
c
      done = .false.
      do while (.not. done)
         read (ipdb,30,err=130,end=130)  record
   30    format (a240)
         remark = record(1:6)
         call upcase (remark)
         if (remark(1:5).eq.'ATOM ' .or. remark.eq.'HETATM' ) then
            done = .true.
            pdbtyp = 'CIF'
            next = 6
            if (remark .eq. 'HETATM')  next = 7
            call getnumb (record,temp,next)
            string = record(next+1:next+4)
            read (string,40,err=120,end=120)  letter4
   40       format (a4)
            string = record(next+5:next+5)
            read (string,50,err=120,end=120)  letter1
   50       format (a1)
            string = record(next+6:next+8)
            read (string,60,err=120,end=120)  letter3
   60       format (a3)
            string = record(next+10:next+10)
            read (string,70,err=120,end=120)  letter1
   70       format (a1)
            next = next + 11
            last = next
            call getnumb (record,temp,next)
            if (next .eq. last) then
               string = record(next:next+3)
               read (string,80,err=120,end=120)  temp
   80          format (i4)
               next = next + 4
            end if
            string = record(next:next)
            read (string,90,err=120,end=120)  letter1
   90       format (a1)
            string = record(next+1:240)
            read (string,*,err=100,end=100)  crd,crd,crd
            goto 110
  100       continue
            string = record(31:38)
            read (string,*,err=120,end=120)  crd
            string = record(39:46)
            read (string,*,err=120,end=120)  crd
            string = record(47:54)
            read (string,*,err=120,end=120)  crd
  110       continue
            pdbtyp = 'PDB'
  120       continue
         end if
      end do
  130 continue
      rewind (unit=ipdb)
c
c     read the coordinates file in either PDB or CIF format
c
      if (pdbtyp .eq. 'PDB')  call readpdb (ipdb)
      if (pdbtyp .eq. 'CIF')  call readcif (ipdb)
      close (unit=ipdb)
      return
      end

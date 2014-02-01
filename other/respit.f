c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2005  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program respit  --  make RESP input files from G03 output  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "respit" gets atomic coordinates, electrostatic potential grid
c     and other values from a Gaussian 03 log file and creates input
c     files for use with the Amber RESP program
c
c
      program respit
      implicit none
      include 'files.i'
      include 'iounit.i'
      integer maxatm,maxgrid
      parameter (maxatm=100)
      parameter (maxgrid=10000)
      integer i,size
      integer charge,vary
      integer ilog,iesp,idat
      integer natm,ngrd,npot
      integer freeunit,trimtext
      integer atnum(maxatm)
      real*8 xa(maxatm),xg(maxgrid)
      real*8 ya(maxatm),yg(maxgrid)
      real*8 za(maxatm),zg(maxgrid)
      real*8 pot(maxgrid)
      logical exist
      character*80 record
      character*80 string
      character*80 title
      character*120 logfile
      character*120 espfile
      character*120 datfile
c
c
c     try to get a filename from the command line arguments
c
      call initial
      call nextarg (logfile,exist)
      if (exist) then
         call basefile (logfile)
         call suffix (logfile,'log')
         call version (logfile,'old')
         inquire (file=logfile,exist=exist)
      end if
c
c     ask for the user specified Gaussian 03 log file
c
      dowhile (.not. exist)
         write (iout,10)
   10    format (/,' Enter Gaussian 03 Log File Name :  ',$)
         read (input,20)  logfile
   20    format (a120)
         call basefile (logfile)
         call suffix (logfile,'log')
         call version (logfile,'old')
         inquire (file=logfile,exist=exist)
      end do
c
c     initialize counters for atoms and potential grid
c
      natm = 0
      ngrd = 0
      npot = 0
c
c     get coordinates and potential from G03 log file
c
      ilog = freeunit ()
      open (unit=ilog,file=logfile,status='old')
      rewind (unit=ilog)
      dowhile (.true.)
         read (ilog,30,err=40,end=40)  record
   30    format (a80)
         if (record(8:20) .eq. 'Atomic Center') then
            natm = natm + 1
            string = record(32:80)
            read (string,*)  xa(natm),ya(natm),za(natm)
         else if (record(7:20) .eq. 'ESP Fit Center') then
            ngrd = ngrd + 1
            string = record(32:80)
            read (string,*)  xg(ngrd),yg(ngrd),zg(ngrd)
         else if (record(7:11) .eq. 'Fit  ') then
            npot = npot + 1
            string = record(12:80)
            read (string,*)  pot(npot)
         end if
      end do
   40 continue
c
c     create RESP input file with coordinates and potential
c
      iesp = freeunit ()
      espfile = filename(1:leng)//'.esp'
      call version (espfile,'new')
      open (unit=iesp,file=espfile,status='new')
      write (iesp,50)  natm,ngrd
   50 format (2i5)
      do i = 1, natm
         write (iesp,60)  xa(i),ya(i),za(i)
   60    format (16x,3e16.6)
      end do
      do i = 1, ngrd
         write (iesp,70)  pot(i),xg(i),yg(i),zg(i)
   70    format (4e16.6)
      end do
      close (unit=iesp)
c
c     give name of RESP coordinate and potential file
c
      write (iout,80)  natm,ngrd
   80 format (/,' Number of Atoms :'i8,6x,
     &           'Number of Grid Points :',i8)
      size = trimtext (espfile)
      write (iout,90)  espfile(1:size)
   90 format (/,' RESP Input Data File Written To :  ',a)
c
c     get title, charge and atomic numbers from G03 log file
c
      rewind (unit=ilog)
      dowhile (.true.)
         read (ilog,100,err=130,end=130)  record
  100    format (a80)
         if (record(2:9) .eq. 'Charge =') then
            charge = 0
            string = record(10:12)
            read (string,*)  charge
         else if (record(27:44) .eq. 'Input orientation:') then
            read (ilog,110,err=130,end=130)
  110       format (///)
            do i = 1, natm
               read (ilog,120,err=130,end=130)  atnum(i)
  120          format (10x,i10)
            end do
         end if
      end do
  130 continue
      close (unit=ilog)
c
c     create RESP input file with control options
c
      idat = freeunit ()
      datfile = filename(1:leng)//'.dat'
      call version (datfile,'new')
      open (unit=idat,file=datfile,status='new')
      title = 'RESP Charge Fitting'
      leng = trimtext (title)
      write (idat,140)  title(1:leng)
  140 format (a)
      write (idat,150)  charge,natm
  150 format (2i5)
      vary = 0
      do i = 1, natm
         write (idat,160)  atnum(i),vary
  160    format (2i5)
      end do
      write (idat,170)
  170 format ()
      close (unit=idat)
c
c     give name of RESP program control options file
c
      size = trimtext (datfile)
      write (iout,180)  datfile(1:size)
  180 format (/,' RESP Control File Written To :  ',a)
      end

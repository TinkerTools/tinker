c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program download  --  parse http log for TINKER downloads  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "download" parses an Apache http logfile and extracts the number
c     and name of machines downloading the complete TINKER molecular
c     modeling package; duplicates are removed so each machine/site is
c     listed only once
c
c
      program download
      implicit none
      integer maxhit
      parameter (maxhit=20000)
      integer i,k,ihit,ilist
      integer ilog,imatch
      integer iload,ifull
      integer iload1,iload2
      integer iload3,iload4
      integer iload5
      integer leng,next
      integer freeunit,trimtext
      integer itag(maxhit)
      character*4 year(maxhit)
      character*5 time(maxhit)
      character*6 date(maxhit)
      character*7 tag(15)
      character*7 cargo(maxhit)
      character*120 logfile
      character*120 site(maxhit)
      character*256 record
      logical exist,match
      data tag  / 'Source ', 'Linux  ', 'SrcLin ', 'MacOSX ',
     &            'SrcMac ', 'LinMac ', 'SrcLM  ', 'Windows',
     &            'SrcWin ', 'LinWin ', 'SrcLW  ', 'MacWin ',
     &            'SrcMW  ', 'LMW    ', 'SrcLMW ' /
c
c
c     initialize the number of http-based TINKER downloads
c
      ihit = 0
c
c     try to get http logfile from command line arguments
c
      call command
      call nextarg (logfile,exist)
      if (exist)  inquire (file=logfile,exist=exist)
c
c     ask the user for the name of the http logfile
c
      dowhile (.not. exist)
         write (*,10)
   10    format (/,' Enter the http Log File Name :  ',$)
         read (*,20)  logfile
   20    format (a120)
         inquire (file=logfile,exist=exist)
      end do
c
c     open the http logfile for parsing
c
      ilog = freeunit ()
      open (unit=ilog,file=logfile,status='old')
      rewind (unit=ilog)
c
c     write a header prior to analyzing the logfile
c
      write (*,30)
   30 format (/,' Downloads of the TINKER Molecular',
     &           ' Modeling Package',/)
c
c     read each line looking for source or executable downloads
c
      dowhile (.true.)
         read (ilog,40,end=50)  record
   40    format (a256)
         iload1 = index(record,'tinker.tar.gz')
         iload1 = max(iload1,index(record,'tinker-5.1'))
         if (iload1 .ne. 0)  iload1 = index(record,'.tar.gz')
         iload2 = index(record,'tinker.zip')
         iload2 = max(iload2,index(record,'tinker-5.1'))
         if (iload2 .ne. 0)  iload2 = index(record,'.zip')
         iload3 = index(record,'linux.tar.gz')
         iload4 = index(record,'macosx.tar.gz')
         iload5 = index(record,'windows.zip')
         iload = max(iload1,iload2,iload3,iload4,iload5)
         ifull = index(record,' 200 ')
         match = (ifull.ne.0 .and. iload.ne.0)
         if (match) then
            ihit = ihit + 1
            next = 1
            call gettext (record,site(ihit),next)
            date(ihit) = '      '
            time(ihit) = '     '
            year(ihit) = '    '
            next = index(record,'- - [') + 4
            time(ihit) = record(next+13:next+17)
            date(ihit) = record(next+1:next+6)
            date(ihit)(3:3) = ' '
            year(ihit) = record(next+8:next+11)
            cargo(ihit) = '       '
            if (iload1.ne.0 .or. iload2.ne.0)  itag(ihit) = 1
            if (iload3 .ne. 0)  itag(ihit) = 2
            if (iload4 .ne. 0)  itag(ihit) = 4
            if (iload5 .ne. 0)  itag(ihit) = 8
         end if
      end do
   50 continue
c
c     account for multiple downloads from the same site
c
      do i = 1, ihit
         do k = i+1, ihit
            if (site(i) .eq. site(k)) then
               if (itag(k).eq.1 .and. mod(itag(i),2).eq.0)
     &            itag(i) = itag(i) + 1
               if (itag(k).eq.2 .and. mod(itag(i),4).le.1)
     &            itag(i) = itag(i) + 2
               if (itag(k).eq.4 .and. mod(itag(i),8).le.3)
     &            itag(i) = itag(i) + 4
               if (itag(k).eq.8 .and. itag(i).le.7)
     &            itag(i) = itag(i) + 8
            end if
         end do
         cargo(i) = tag(itag(i))
      end do
c
c     print the list of download sites, skipping duplicates
c
      ilist = 0
      do i = 1, ihit
         do k = 1, i-1
            if (site(i) .eq. site(k))  goto 70
         end do
         ilist = ilist + 1
         leng = trimtext (site(i))
         write (*,60)  ilist,cargo(i),time(i),date(i),
     &                 year(i),site(i)(1:leng)
   60    format (i5,2x,a7,3x,a5,3x,a6,1x,a4,3x,a)
   70    continue
      end do
      end

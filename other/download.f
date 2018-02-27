c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program download  --  parse ftplog for Tinker downloads  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "download" parses an ftp logfile produced under Digital Unix
c     and extracts the number and name of machines downloading the
c     complete Tinker molecular modeling package; duplicates are
c     removed so each machine/site is listed only once
c
c     if downloads of Tinker are not matched with a site, then the
c     logfile contains an unexpected format or the maxtag parameter
c     needs to be set to a larger value 
c
c
      program download
      implicit none
      integer maxtag,maxhit
      parameter (maxtag=10000)
      parameter (maxhit=10000)
      integer i,k,itag,ihit,ilist
      integer iconn,imatch,ilog
      integer start,stop
      integer size,sizemin
      integer leng,next
      integer freeunit,trimtext
      integer process(maxhit)
      character*3 day(maxhit)
      character*4 year(maxhit)
      character*5 time(maxhit)
      character*6 numeral,tag(0:maxtag)
      character*6 date(maxhit)
      character*240 logfile
      character*256 record
      character*256 string(0:maxtag)
      character*256 site(maxhit)
      logical exist,connect,match
c
c
c     initialize the number of ftp connections and downloads
c
      itag = 0
      ihit = 0
      sizemin = 4000000
c
c     try to get ftp logfile from command line arguments
c
      call command
      call nextarg (logfile,exist)
      if (exist)  inquire (file=logfile,exist=exist)
c
c     ask the user for the name of the ftp logfile
c
      dowhile (.not. exist)
         write (*,10)
   10    format (/,' Enter the ftp Log File Name :  ',$)
         read (*,20)  logfile
   20    format (a240)
         inquire (file=logfile,exist=exist)
      end do
c
c     open the ftp logfile for parsing
c
      ilog = freeunit ()
      open (unit=ilog,file=logfile,status='old')
      rewind (unit=ilog)
c
c     write a header prior to analyzing the logfile
c
      write (*,30)
   30 format (/,' Downloads of the Tinker Molecular',
     &           ' Modeling Package',/)
c
c     read each line looking for connects and downloads
c
      dowhile (.true.)
         read (ilog,40,end=60)  record
   40    format (a256)
         iconn = index(record,'connection from')
         imatch = max(index(record,'tinker.tar.gz succeeded'),
     &                index(record,'tinker.zip succeeded'),
     &                index(record,'tinker-linux.sh.gz succeeded'),
     &                index(record,'tinker-macosx.sit succeeded'),
     &                index(record,'tinker-windows.exe succeeded'))
         if (imatch .ne. 0)  imatch = index(record,'retrieve')
         connect = (iconn .ne. 0)
         match = (imatch .ne. 0)
         if (match) then
            start = index(record,'succeeded') + 10
            call getnumb (record,size,start)
            if (size .lt. sizemin)  match = .false.
         end if
         if (connect .or. match) then
            start = index(record,'[') + 1
            stop = index(record,']') - 1
            numeral = record(start:stop)
         end if
         if (connect) then
            itag = itag + 1
            k = mod(itag,maxtag)
            tag(k) = numeral
            start = iconn + 16
            string(k) = record(start:)
         end if
         if (match) then
            ihit = ihit + 1
            next = 1
            call getnumb (numeral,process(ihit),next)
            site(ihit) = '*** Tinker Download Not Matched ***'
            day(ihit) = '   '
            date(ihit) = '      '
            time(ihit) = '     '
            year(ihit) = '    '
            start = itag
            stop = max(1,itag-maxtag+1)
            do i = start, stop, -1
               k = mod(i,maxtag)
               if (numeral .eq. tag(k)) then
                  next = 1
                  call gettext (string(k),site(ihit),next)
                  call lowcase (site(ihit))
                  day(ihit) = string(k)(next+4:next+6)
                  date(ihit) = string(k)(next+8:next+13)
                  time(ihit) = string(k)(next+15:next+19)
                  year(ihit) = string(k)(next+24:next+27)
                  goto 50
               end if
            end do
   50       continue
         end if
      end do
   60 continue
c
c     scan the list of hits, skipping over duplicate sites
c
      ilist = 0
      do i = 1, ihit
         if (site(i) .ne. '*** Tinker Download Not Matched ***') then
            do k = 1, i-1
               if (site(i) .eq. site(k))  goto 80
            end do
         end if
c
c     print the final set of unique download sites
c
         ilist = ilist + 1
         leng = trimtext (site(i))
         write (*,70)  ilist,process(i),day(i),time(i),date(i),
     &                 year(i),site(i)(1:leng)
   70    format (i5,2x,i6,3x,a3,1x,a5,3x,a6,1x,a4,3x,a)
   80    continue
      end do
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine getprm  --  get force field parameter files  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "getprm" finds any potential energy parameter files and
c     then opens and reads the parameters
c
c
      subroutine getprm
      use argue
      use files
      use inform
      use iounit
      use keys
      use params
      implicit none
      integer maxfile
      parameter (maxfile=12)
      integer i,j,iprm
      integer nfile,nask
      integer trimtext,next
      integer freeunit
      logical exist,useprm
      character*4 none
      character*20 keyword
      character*240 prmfile
      character*240 prefix
      character*240 record
      character*240 string
      character*240 ifile(maxfile)
c
c
c     set default usage and number of parameter files
c
      useprm = .true.
      nfile = 0
c
c     check for parameter file with base name of current system
c
      prmfile = filename(1:leng)//'.prm'
      call version (prmfile,'old')
      inquire (file=prmfile,exist=exist)
      if (exist) then
         nfile = nfile + 1
         ifile(nfile) = prmfile
      end if
c
c     check for the existence of a generic parameter file
c
      if (ldir .eq. 0) then
         prmfile = 'tinker.prm'
      else
         prmfile = filename(1:ldir)//'tinker.prm'
      end if
      call version (prmfile,'old')
      inquire (file=prmfile,exist=exist)
      if (exist) then
         nfile = nfile + 1
         ifile(nfile) = prmfile
      end if
c
c     try to get a parameter filename from the command line
c
      do i = 1, narg-1
         string = arg(i)
         call upcase (string)
         if (string(1:2) .eq. '-P') then
            prmfile = arg(i+1)
            if (prmfile(1:2) .eq. '~/') then
               call getenv ('HOME',prefix)
               prmfile = prefix(1:trimtext(prefix))//
     &                      prmfile(2:trimtext(prmfile))
            end if
            call suffix (prmfile,'prm','old')
            inquire (file=prmfile,exist=exist)
            if (exist) then
               nfile = nfile + 1
               ifile(nfile) = prmfile
            else
               write (iout,10)
   10          format (/,' GETPRM  --  Parameter File Named on',
     &                    ' Command Line not Found')
               call fatal
            end if
         end if
      end do
c
c     search the keyword list for the parameter filename
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11).eq.'PARAMETERS '
     &          .or. keyword(1:10).eq.'PARAMETER ') then
            string = record(next:240)
            next = 1
            call getstring (string,prmfile,next)
            if (next .eq. 1)  call gettext (string,prmfile,next)
            if (prmfile(1:2) .eq. '~/') then
               call getenv ('HOME',prefix)
               prmfile = prefix(1:trimtext(prefix))//
     &                      prmfile(2:trimtext(prmfile))
            end if
            call suffix (prmfile,'prm','old')
            inquire (file=prmfile,exist=exist)
            if (exist) then
               do j = 1, nfile
                  if (prmfile .eq. ifile(j))  goto 20
               end do
               nfile = nfile + 1
               ifile(nfile) = prmfile
   20          continue
            else
               none = prmfile(1:4)
               call upcase (none)
               if (none .eq. 'NONE')  useprm = .false.
            end if
         end if
      end do
      if (.not. useprm)  nfile = 0
c
c     if necessary, ask for the parameter filename
c
      if (useprm .and. nfile.eq.0) then
         nask = 0
         exist = .false.
         do while (.not.exist .and. nask.lt.maxask)
            nask = nask + 1
            write (iout,30)
   30       format (/,' Enter Parameter File Name [<Enter>=NONE] :  ',$)
            read (input,40)  prmfile
   40       format (a240)
            next = 1
            call getword (prmfile,none,next)
            call upcase (none)
            if (next.eq.1 .or. none.eq.'NONE') then
               exist = .true.
               useprm = .false.
            else
               if (prmfile(1:2) .eq. '~/') then
                  call getenv ('HOME',prefix)
                  prmfile = prefix(1:trimtext(prefix))//
     &                         prmfile(2:trimtext(prmfile))
               end if
               call suffix (prmfile,'prm','old')
               inquire (file=prmfile,exist=exist)
               if (exist) then
                  nfile = nfile + 1
                  ifile(nfile) = prmfile
               end if
            end if
         end do
      end if
c
c     check to make sure a parameter file is available
c
      if (useprm .and. nfile.eq.0) then
         write (iout,50)
   50    format (/,' GETPRM  --  A Valid Parameter File',
     &              ' was not Provided')
         call fatal
      end if
c
c     read the parameter files and count the total lines
c
      nprm = 0
      if (useprm) then
         do i = 1, nfile
            iprm = freeunit ()
            prmfile = ifile(i)
            open (unit=iprm,file=prmfile,status='old')
            rewind (unit=iprm)
            do while (.true.)
               read (iprm,60,err=70,end=70)
   60          format ()
               nprm = nprm + 1
            end do
   70       continue
            close (unit=iprm)
         end do
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(prmline))  deallocate (prmline)
      allocate (prmline(nprm))
c
c     reread the parameter files and store for latter use
c
      nprm = 0
      do i = 1, nfile
         iprm = freeunit ()
         prmfile = ifile(i)
         open (unit=iprm,file=prmfile,status='old')
         rewind (unit=iprm)
         dowhile (.true.)
            read (iprm,80,err=90,end=90)  record
   80       format (a240)
            nprm = nprm + 1
            prmline(nprm) = record
         end do
   90    continue
         close (unit=iprm)
      end do
c
c     convert underbar characters to dashes in all keywords
c
      do i = 1, nprm
         next = 1
         record = prmline(i)
         call gettext (record,keyword,next)
         do j = 1, next-1
            if (record(j:j) .eq. '_')  record(j:j) = '-'
         end do
         prmline(i) = record
      end do
c
c     count and allocate memory for the parameter values
c
      call setprm
c
c     initialize force field control and parameter values
c
      call initprm
c
c     get control and parameter values from the parameter file
c
      if (useprm)  call readprm
      return
      end

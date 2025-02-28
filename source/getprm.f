c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine getprm  --  get force field parameter file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "getprm" finds the potential energy parameter file
c     and then opens and reads the parameters
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
      integer i,j,iprm
      integer nask,next
      integer freeunit
      integer trimtext
      logical exist,useprm
      character*4 none
      character*20 keyword
      character*240 prmfile
      character*240 prefix
      character*240 record
      character*240 string
c
c
c     set the default name for the parameter file
c
      useprm = .true.
      prmfile = filename(1:leng)//'.prm'
c
c     try to get a parameter filename from the command line
c
      exist = .false.
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
            if (.not. exist) then
               write (iout,10)
   10          format (/,' GETPRM  --  Parameter File Specified',
     &                    ' on Command Line not Found')
               call fatal
            end if
         end if
      end do
c
c     search the keyword list for the parameter filename
c
      if (.not. exist) then
         do i = 1, nkey
            next = 1
            record = keyline(i)
            call gettext (record,keyword,next)
            call upcase (keyword)
            if (keyword(1:11).eq.'PARAMETERS '
     &             .or. keyword(1:10).eq.'PARAMETER ') then
               string = record(next:240)
               next = 1
               call getstring (string,prmfile,next)
               if (next .eq. 1)  call gettext (string,prmfile,next)
               if (prmfile(1:2) .eq. '~/') then
                  call getenv ('HOME',prefix)
                  prmfile = prefix(1:trimtext(prefix))//
     &                         prmfile(2:trimtext(prmfile))
               end if
               call suffix (prmfile,'prm','old')
               inquire (file=prmfile,exist=exist)
            end if
         end do
      end if
c
c     test for user specified absence of a parameter file
c
      if (.not. exist) then
         none = prmfile(1:4)
         call upcase (none)
         if (none .eq. 'NONE') then
            exist = .true.
            useprm = .false.
         end if
      end if
c
c     if necessary, ask for the parameter filename
c
      nask = 0
      do while (.not.exist .and. nask.lt.maxask)
         nask = nask + 1
         write (iout,20)
   20    format (/,' Enter Parameter File Name [<Enter>=NONE] :  ',$)
         read (input,30)  prmfile
   30    format (a240)
         next = 1
         call getword (prmfile,none,next)
         call upcase (none)
         if (next .eq. 1) then
            exist = .true.
            useprm = .false.
         else if (none.eq.'NONE' .and. next.eq.5) then
            exist = .true.
            useprm = .false.
         else
            if (prmfile(1:2) .eq. '~/') then
               call getenv ('HOME',prefix)
               prmfile = prefix(1:trimtext(prefix))//
     &                      prmfile(2:trimtext(prmfile))
            end if
            call suffix (prmfile,'prm','old')
            inquire (file=prmfile,exist=exist)
         end if
      end do
      if (.not. exist)  call fatal
c
c     read the parameter file to get number of lines
c
      nprm = 0
      if (useprm) then
         iprm = freeunit ()
         open (unit=iprm,file=prmfile,status='old')
         rewind (unit=iprm)
         do while (.true.)
            read (iprm,40,err=50,end=50)
   40       format ()
            nprm = nprm + 1
         end do
   50    continue
         rewind (unit=iprm)
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(prmline))  deallocate (prmline)
      allocate (prmline(nprm))
c
c     reread the parameter file and store for latter use
c
      do i = 1, nprm
         read (iprm,60,err=70,end=70)  record
   60    format (a240)
         prmline(i) = record
      end do
   70 continue
      close (unit=iprm)
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

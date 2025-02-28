c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine getkey  --  find and store contents of keyfile  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "getkey" finds a valid keyfile and stores its contents as
c     line images for subsequent keyword parameter searching
c
c
      subroutine getkey
      use argue
      use files
      use iounit
      use keys
      use openmp
      implicit none
      integer i,j,ikey
      integer next,length
      integer freeunit
      integer trimtext
      logical exist,header
      character*20 keyword
      character*240 keyfile
      character*240 comment
      character*240 record
      character*240 string
c
c
c     check for a keyfile specified on command line
c
      exist = .false.
      do i = 1, narg-1
         string = arg(i)
         call upcase (string)
         if (string(1:2) .eq. '-K') then
            keyfile = arg(i+1)
            call suffix (keyfile,'key','old')
            inquire (file=keyfile,exist=exist)
            if (.not. exist) then
               write (iout,10)
   10          format (/,' GETKEY  --  Keyfile Specified',
     &                    ' on Command Line was not Found')
               call fatal
            end if
         end if
      end do
c
c     try to get keyfile from base name of current system
c
      if (.not. exist) then
         keyfile = filename(1:leng)//'.key'
         call version (keyfile,'old')
         inquire (file=keyfile,exist=exist)
      end if
c
c     check for the existence of a generic keyfile
c
      if (.not. exist) then
         if (ldir .eq. 0) then
            keyfile = 'tinker.key'
         else
            keyfile = filename(1:ldir)//'tinker.key'
         end if
         call version (keyfile,'old')
         inquire (file=keyfile,exist=exist)
      end if
c
c     read the keyfile to get number of lines
c
      nkey = 0
      if (exist) then
         ikey = freeunit ()
         open (unit=ikey,file=keyfile,status='old')
         rewind (unit=ikey)
         do while (.true.)
            read (ikey,20,err=30,end=30)
   20       format ()
            nkey = nkey + 1
         end do
   30    continue
         rewind (unit=ikey)
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(keyline))  deallocate (keyline)
      allocate (keyline(nkey))
c
c     reread the keyfile and store for latter use
c
      do i = 1, nkey
         read (ikey,40,err=50,end=50)  record
   40    format (a240)
         keyline(i) = record
      end do
   50 continue
      close (unit=ikey)
c
c     convert underbar characters to dashes in all keywords
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         do j = 1, next-1
            if (record(j:j) .eq. '_')  record(j:j) = '-'
         end do
         keyline(i) = record
      end do
c
c     check for comment lines to be echoed to the output
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ECHO ') then
            comment = record(next:240)
            length = trimtext (comment)
            if (header) then
               header = .false.
               write (iout,60)
   60          format ()
            end if
            if (length .eq. 0) then
               write (iout,70)
   70          format ()
            else
               write (iout,80)  comment(1:length)
   80          format (a)
            end if
         end if
      end do
c
c     set number of OpenMP threads for parallelization
c
!$    do i = 1, nkey
!$       next = 1
!$       record = keyline(i)
!$       call upcase (record)
!$       call gettext (record,keyword,next)
!$       string = record(next:240)
!$       if (keyword(1:15) .eq. 'OPENMP-THREADS ') then
!$          read (string,*,err=90,end=90)  nthread
!$          call omp_set_num_threads (nthread)
!$       end if
!$ 90    continue
!$    end do
c
c     check for number of OpenMP threads on command line
c
!$    do i = 1, narg-1
!$       string = arg(i)
!$       call upcase (string)
!$       if (string(1:2) .eq. '-T') then
!$          next = 1
!$          string = arg(i+1)
!$          call getnumb (string,nthread,next)
!$          if (nthread .eq. 0)  nthread = 1
!$          call omp_set_num_threads (nthread)
!$       end if
!$    end do
      return
      end

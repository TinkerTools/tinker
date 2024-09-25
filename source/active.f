c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine active  --  set the list of active atoms  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "active" sets the list of atoms that are used during
c     coordinate manipulation or potential energy calculations
c
c
      subroutine active
      use atoms
      use inform
      use iounit
      use keys
      use usage
      implicit none
      integer i,j,next
      integer nmobile,nfixed
      integer center,nsphere
      integer, allocatable :: mobile(:)
      integer, allocatable :: fixed(:)
      real*8 xcenter,ycenter,zcenter
      real*8 radius,radius2,dist2
      character*20 keyword
      character*240 record
      character*240 string
      logical header
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(iuse))  deallocate (iuse)
      if (allocated(use))  deallocate (use)
      allocate (iuse(n))
      allocate (use(0:n))
c
c     perform dynamic allocation of some local arrays
c
      allocate (mobile(n))
      allocate (fixed(n))
c
c     set defaults for the numbers and lists of active atoms
c
      nuse = n
      use(0) = .false.
      do i = 1, n
         use(i) = .true.
      end do
      nmobile = 0
      nfixed = 0
      do i = 1, n
         mobile(i) = 0
         fixed(i) = 0
      end do
      nsphere = 0
c
c     get any keywords containing active atom parameters
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
c
c     get any lists of atoms whose coordinates are active
c
         if (keyword(1:7) .eq. 'ACTIVE ') then
            read (string,*,err=10,end=10)  (mobile(i),i=nmobile+1,n)
   10       continue
            do while (mobile(nmobile+1) .ne. 0)
               nmobile = nmobile + 1
            end do
c
c     get any lists of atoms whose coordinates are inactive
c
         else if (keyword(1:9) .eq. 'INACTIVE ') then
            read (string,*,err=20,end=20)  (fixed(i),i=nfixed+1,n)
   20       continue
            do while (fixed(nfixed+1) .ne. 0)
               nfixed = nfixed + 1
            end do
c
c     get the center and radius of the sphere of active atoms
c
         else if (keyword(1:14) .eq. 'ACTIVE-SPHERE ') then
            center = 0
            xcenter = 0.0d0
            ycenter = 0.0d0
            zcenter = 0.0d0
            radius = 0.0d0
            read (string,*,err=30,end=30)  xcenter,ycenter,
     &                                     zcenter,radius
   30       continue
            if (radius .eq. 0.0d0) then
               read (string,*,err=60,end=60)  center,radius
               xcenter = x(center)
               ycenter = y(center)
               zcenter = z(center)
            end if
            nsphere = nsphere + 1
            if (nsphere .eq. 1) then
               nuse = 0
               do i = 1, n
                  use(i) = .false.
               end do
               if (verbose) then
                  write (iout,40)
   40             format (/,' Spheres used to Select Active Atoms :',
     &                    //,3x,'Atom Center',11x,'Coordinates',
     &                       12x,'Radius',6x,'# Active Atoms')
               end if
            end if
            radius2 = radius * radius
            do i = 1, n
               if (.not. use(i)) then
                  dist2 = (x(i)-xcenter)**2 + (y(i)-ycenter)**2
     &                            + (z(i)-zcenter)**2
                  if (dist2 .le. radius2) then
                     nuse = nuse + 1
                     use(i) = .true.
                  end if
               end if
            end do
            if (verbose) then
               write (iout,50)  center,xcenter,ycenter,
     &                          zcenter,radius,nuse
   50          format (2x,i8,6x,3f9.2,2x,f9.2,7x,i8)
            end if
   60       continue
         end if
      end do
c
c     remove active or inactive atoms not in the system
c
      header = .true.
      do i = 1, n
         if (abs(mobile(i)) .gt. n) then
            mobile(i) = 0
            if (header) then
               header = .false.
               write (iout,70)
   70          format (/,' ACTIVE  --  Warning, Illegal Atom Number',
     &                    ' in ACTIVE Atom List')
            end if
         end if
      end do
      header = .true.
      do i = 1, n
         if (abs(fixed(i)) .gt. n) then
            fixed(i) = 0
            if (header) then
               header = .false.
               write (iout,80)
   80          format (/,' ACTIVE  --  Warning, Illegal Atom Number',
     &                    ' in INACTIVE Atom List')
            end if
         end if
      end do
c
c     set active atoms to those marked as not inactive
c
      i = 1
      do while (fixed(i) .ne. 0)
         if (fixed(i) .gt. 0) then
            j = fixed(i)
            if (use(j)) then
               use(fixed(i)) = .false.
               nuse = nuse - 1
            end if
            i = i + 1
         else
            do j = abs(fixed(i)), abs(fixed(i+1))
               if (use(j)) then
                  use(j) = .false.
                  nuse = nuse - 1
               end if
            end do
            i = i + 2
         end if
      end do
c
c     set active atoms to only those marked as active
c
      i = 1
      do while (mobile(i) .ne. 0)
         if (i .eq. 1) then
            nuse = 0
            do j = 1, n
               use(j) = .false.
            end do
         end if
         if (mobile(i) .gt. 0) then
            j = mobile(i)
            if (.not. use(j)) then
               use(j) = .true.
               nuse = nuse + 1
            end if
            i = i + 1
         else
            do j = abs(mobile(i)), abs(mobile(i+1))
               if (.not. use(j)) then
                  use(j) = .true.
                  nuse = nuse + 1
               end if
            end do
            i = i + 2
         end if
      end do
c
c     use logical array to set the list of active atoms
c
      j = 0
      do i = 1, n
         if (use(i)) then
            j = j + 1
            iuse(j) = i
         end if
      end do
c
c     output the final list of the active atoms
c
      if (debug .and. nuse.gt.0 .and. nuse.lt.n) then
         write (iout,90)
   90    format (/,' List of Active Atoms for Energy',
     &              ' Calculations :',/)
         write (iout,100)  (iuse(i),i=1,nuse)
  100    format (3x,10i7)
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mobile)
      deallocate (fixed)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine saveonly  --  set the list of save coord atoms  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "saveonly" sets the list of atoms that are used during
c     coordinate saving routines
c
c
      subroutine saveonly
      use atoms
      use iounit
      use keys
      use output
      implicit none
      integer i,j,next
      integer nfixed
      integer, allocatable :: fixed(:)
      character*20 keyword
      character*240 record
      character*240 string
      logical header
      logical, allocatable :: saved(:)
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ionly))  deallocate (ionly)
      if (allocated(ionlyinv))  deallocate (ionlyinv)
      allocate (ionly(n))
      allocate (ionlyinv(n))
c
c     perform dynamic allocation of some local arrays
c
      allocate (fixed(n))
      allocate (saved(n))
c
c     set defaults for the numbers and lists of saved atoms
c
      onlysave = .false.
      nonly = 0
      do i = 1, n
         ionly(i) = 0
         ionlyinv(i) = 0
      end do
      nfixed = 0
      do i = 1, n
         fixed(i) = 0
         saved(i) = .false.
      end do
c
c     get any keywords containing save-only atom parameters
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
c
c     get any lists of atoms whose coordinates should be saved
c
         if (keyword(1:10) .eq. 'SAVE-ONLY ') then
            read (string,*,err=10,end=10)  (fixed(i),i=nfixed+1,n)
   10       continue
            do while (fixed(nfixed+1) .ne. 0)
               nfixed = nfixed + 1
            end do
         end if
      end do
c
c     remove saved atoms not in the system
c
      header = .true.
      do i = 1, n
         if (abs(fixed(i)) .gt. n) then
            fixed(i) = 0
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' SAVEONLY  --  Warning, Illegal Atom Number',
     &                    ' in SAVE-ONLY Atom List')
            end if
         end if
      end do
c
c     set saved atoms to only those marked as save
c
      i = 1
      do while (fixed(i) .ne. 0)
         if (fixed(i) .gt. 0) then
            j = fixed(i)
            saved(j) = .true.
            i = i + 1
         else
            do j = abs(fixed(i)), abs(fixed(i+1))
               saved(j) = .true.
            end do
            i = i + 2
         end if
      end do
      do i = 1, n
         if (saved(i)) then
            nonly = nonly + 1
            ionly(nonly) = i
            ionlyinv(i) = nonly
         end if
      end do
      if (nonly > 0)  onlysave = .true.
c
c     perform deallocation of some local arrays
c
      deallocate (fixed)
      deallocate (saved)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine msystem  --  set exclusion for moment of system  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "msystem" sets the list of atoms that are excluded while
c     computing moment of system
c
c
      subroutine msystem
      use atoms
      use iounit
      use keys
      use moment
      implicit none
      integer i,j,next
      integer nfixed
      integer, allocatable :: fixed(:)
      character*20 keyword
      character*240 record
      character*240 string
      logical header
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(momuse))  deallocate (momuse)
      allocate (momuse(n))
c
c     perform dynamic allocation of some local arrays
c
      allocate (fixed(n))
c
c     set defaults for the numbers and lists of atoms to be used
c
      do i = 1, n
         momuse(i) = .true.
      end do
      nfixed = 0
      do i = 1, n
         fixed(i) = 0
      end do
c
c     get any keywords containing exc-moment atom parameters
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
c
c     get any lists of atoms whose coordinates should be used
c
         if (keyword(1:13) .eq. 'EXC-MOMENT ') then
            read (string,*,err=10,end=10)  (fixed(i),i=nfixed+1,n)
   10       continue
            do while (fixed(nfixed+1) .ne. 0)
               nfixed = nfixed + 1
            end do
         end if
      end do
c
c     remove used atoms not in the system
c
      header = .true.
      do i = 1, n
         if (abs(fixed(i)) .gt. n) then
            fixed(i) = 0
            if (header) then
               header = .false.
               write (iout,20)
   20          format (/,' MSYSTEM  --  Warning, Illegal Atom Number',
     &                    ' in EXC-MOMENT Atom List')
            end if
         end if
      end do
c
c     set inactive atoms to false
c
      i = 1
      do while (fixed(i) .ne. 0)
         if (fixed(i) .gt. 0) then
            j = fixed(i)
            momuse(j) = .false.
            i = i + 1
         else
            do j = abs(fixed(i)), abs(fixed(i+1))
               momuse(j) = .false.
            end do
            i = i + 2
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fixed)
      return
      end

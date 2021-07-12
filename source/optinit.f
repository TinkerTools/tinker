c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2018  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine optinit  --  initialize structure optimization  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "optinit" initializes values and keywords used by multiple
c     structure optimization methods
c
c
      subroutine optinit
      use inform
      use keys
      use output
      use potent
      implicit none
      integer i,next
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set default values for optimization parameters
c
      iprint = -1
      iwrite = -1
      frcsave = .false.
      uindsave = .false.
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:9) .eq. 'PRINTOUT ') then
            read (string,*,err=10,end=10)  iprint
         else if (keyword(1:9) .eq. 'WRITEOUT ') then
            read (string,*,err=10,end=10)  iwrite
         else if (keyword(1:11) .eq. 'SAVE-FORCE ') then
            frcsave = .true.
         else if (keyword(1:13) .eq. 'SAVE-INDUCED ') then
            uindsave = .true.
         end if
   10    continue
      end do
c
c     check for use of induced dipole prediction methods
c
      if (use_polar)  call predict
      return
      end

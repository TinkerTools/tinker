c
c
c     ###################################################
c     ##                                               ##
c     ##                                               ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine stresinit  --  set stress tenor output         ##
c     ##                                                            ##
c     ################################################################
c
c
c     "stresinit" checks if user requests outputting the stress tensor  
c     and if so, how often.
c
c     NOTE: This is the minimum amount of code required for functionality.
c     A future "stresinit" update could add logic to verify that correct 
c     ensemble (NVT) and thermostat (Velocity-Scaling) are selected.
c
c
      subroutine stresinit (dt)
      use keys
      use strvar
      implicit none
      integer i,next
      real*8 dt
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set default values for stress tensor output variables
c
      stresav = .false.
      stresfrq = 10.0d0
c
c     search keywords for stress output parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:12) .eq. 'SAVE-STRESS ') then
            stresav = .true.
         else if (keyword(1:12) .eq. 'STRESS-FREQ ') then
            read (string,*,err=10,end=10) stresfrq
         end if
   10    continue
      end do
      if (stresav) then
         stresfrq = 0.001d0 * stresfrq
         istress = nint(stresfrq/dt)
      end if
      return
      end

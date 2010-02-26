c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine calendar  --  find the current date and time  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "calendar" returns the current time as a set of integer values
c     representing the year, month, day, hour, minute and second
c
c     note only one of the various implementations below should
c     be activated by removing comment characters
c
c
      subroutine calendar (year,month,day,hour,minute,second)
      implicit none
      integer year,month
      integer day,hour
      integer minute,second
c
c
c     use the standard "date_and_time" intrinsic function
c
c     integer values(8)
c     character*5 zone
c     character*8 date
c     character*10 time
c     call date_and_time (date,time,zone,values)
c     year = values(1)
c     month = values(2)
c     day = values(3)
c     hour = values(5)
c     minute = values(6)
c     second = values(7)
c
c     use the obsolete "itime" and "idate" intrinsic functions
c
      integer hms(3)
      call itime (hms)
      hour = hms(1)
      minute = hms(2)
      second = hms(3)
      call idate (month,day,year)
      return
      end

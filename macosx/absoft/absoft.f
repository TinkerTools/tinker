c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  function etime  --  system-level elapsed time routine  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "etime" returns actual total CPU time since the start of
c     executation for the calling process as well an an array with
c     the separate user and system times
c
c
      function etime (tarray)
      implicit none
      integer times,result
      integer array(4)
      real etime,tarray(2)
c
c
c     set the user and system PCU times via a system call
c
      result = times (array)
      tarray(1) = array(1) / 100.0
      tarray(2) = array(2) / 100.0
c
c     set returned CPU time to sum of system and user times
c
      etime = tarray(1) + tarray(2)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine itime  --  system-level wall-clock time routine  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "itime" returns the current wall-clock time in an array with
c     the order: hour, minute, second
c
c
      subroutine itime (iarray)
      implicit none
      integer*4 iarray(3)
      integer localtime,now
      structure /tm/
         integer*4 tm_sec
         integer*4 tm_min
         integer*4 tm_hour
         integer*4 tm_mday
         integer*4 tm_mon
         integer*4 tm_year
         integer*4 tm_wday
         integer*4 tm_yday
         integer*4 tm_isdst
         integer*4 tm_gmtoff
         integer*4 tm_zone
      end structure
      record /tm/ clock
      pointer (p_clock,clock)
c
c
c     set the current time via a system call
c
      call time (now)
      p_clock = localtime (now)
      iarray(1) = clock.tm_hour
      iarray(2) = clock.tm_min
      iarray(3) = clock.tm_sec
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine idate  --  system-level calendar date routine  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "idate" returns the current calendar date in an array with
c     the order: day, month, year
c
c
      subroutine idate (tarray)
      implicit none
      integer*4 tarray(3)
      integer localtime,now
      structure /tm/
         integer*4 tm_sec
         integer*4 tm_min
         integer*4 tm_hour
         integer*4 tm_mday
         integer*4 tm_mon
         integer*4 tm_year
         integer*4 tm_wday
         integer*4 tm_yday
         integer*4 tm_isdst
         integer*4 tm_gmtoff
         integer*4 tm_zone
      end structure
      record /tm/ clock
      pointer (p_clock,clock)
c
c
c     set the current date via a system call
c
      call time (now)
      p_clock = localtime(now)
      tarray(1) = clock.tm_mday
      tarray(2) = clock.tm_mon + 1
      tarray(3) = clock.tm_year + 1900
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine getarg  --  system-level argument list routine  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "getarg" returns the k-th command line argument present at
c     program startup as a character string
c
c     Note: This function will not work if the main program was
c     not built under Absoft Fortran
c
c
      subroutine getarg (k,argument)
      implicit none
      integer*4 i,k
      character*(*) argument
      integer*4 argc
      integer*4 C_string_pointers(0:z'7fffffff')
      pointer (argv,C_string_pointers)
      integer*4 environ
      common /ABSOFT__ARGS/ argc,argv,environ
      character the_arg(z'7fffffff')
      pointer (arg_pointer,the_arg)
c
c
c     initialize the argument and get a pointer to its value
c
      argument = ' '
      if (k .ge. argc)  return
      arg_pointer = C_string_pointers (k)
c
c     copy the argument into the return character string
c
      do i = 1, len(argument)
         if (the_arg(i) .eq. char(0))  goto 10
         argument(i:i) = the_arg(i)
      end do
   10 continue
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function iargc  --  system-level argument count routine  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "iargc" returns the number of arguments which were entered
c     on the command line at the start of program execution
c
c     Note: This function will not work if the main program was
c     not built under Absoft Fortran
c
c
      function iargc ()
      implicit none
      integer*4 iargc
      integer*4 argc,argv,environ
      common /ABSOFT__ARGS/ argc,argv,environ
c
c
c     get the number of command line arguments; decrement
c     by one so the initial value is the zeroth argument
c
      iargc = argc - 1
      return
      end

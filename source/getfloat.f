c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine getfloat  --  extract float from a string  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "getfloat" searches an input string for the first floating
c     point number and puts the numeric value in "number"; returns
c     zero with "next" unchanged if no floating point value is found
c
c     variables and parameters:
c
c     string    input character string to be searched
c     number    output with the first number in the variable
c     next      input with first position of search string;
c                 output with the position following the number
c
c
      subroutine getfloat (string,number,next)
      implicit none
      integer next
      integer initial
      real*8 number
      logical numeral
      character*40 text
      character*(*) string
c
c
c     initialize number and flag for presence of a number
c
      number = 0.0d0
      numeral = .false.
c
c     search the string for the first floating point number
c
      initial = next
      call gettext (string,text,next)
      read (text,*,err=10,end=10)  number
      numeral = .true.
   10 continue
      if (.not. numeral)  next = initial
      return
      end

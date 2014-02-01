c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2005  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine getreal  --  extract a real from a string  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "getreal" searchs an input string from left to right for a
c     real number and puts the numeric value in "number"; returns
c     zero with "next" unchanged if no real number value is found
c
c     variables and parameters:
c
c     string    input character string to be searched
c     number    output with the first real number found
c     next      input with first position of search string;
c                 output with the position following text
c
c
      subroutine getreal (string,number,next)
      implicit none
      include 'ascii.i'
      integer i,j
      integer len,length
      integer next,code
      integer first,last
      integer initial,final
      real*8 number
      logical numeral
      character*(*) string
c
c
c     initialize number and get the input text string length
c
      number = 0.0d0
      numeral = .false.
      length = len(string(next:))
c
c     move through the string one character at a time,
c     searching for the first non-blank character
c
      first = next
      last = 0
      initial = next
      final = next + length - 1
      do i = initial, final
         code = ichar(string(i:i))
         if (code.ne.space .and. code.ne.tab) then
            first = i
            do j = i+1, final
               code = ichar(string(j:j))
               if (code.eq.space .or. code.eq.tab .or.
     &             code.eq.comma .or. code.eq.semicolon .or.
     &             code.eq.colon .or. code.eq.underbar) then
                  last = j - 1
                  next = j
                  read (string(first:last),*,err=10)  number
                  numeral = .true.
                  goto 20
   10             continue
               end if
            end do
            last = final
            next = last + 1
         end if
      end do
   20 continue
c
c     reset to the start of the string if no value was found
c
      if (.not. numeral)  next = initial
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine suffix  --  test for default file extension  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "suffix" checks a filename for the presence of an extension,
c     and appends an extension and version if none is found
c
c
      subroutine suffix (string,extension,status)
      use ascii
      implicit none
      integer i,k
      integer leng,lext
      integer trimtext
      logical exist
      character*1 letter
      character*3 status
      character*(*) string
      character*(*) extension
c
c
c     get the full length of the current filename
c
      leng = trimtext (string)
      lext = trimtext (extension)
c
c     check for an extension on the current filename
c
      k = leng
      do i = 1, leng
         letter = string(i:i)
         if (letter .eq. '/')  k = leng
c        if (letter .eq. '\')  k = leng
         if (ichar(letter) .eq. backslash)  k = leng
         if (letter .eq. ']')  k = leng
         if (letter .eq. ':')  k = leng
         if (letter .eq. '~')  k = leng
         if (letter .eq. '.')  k = i - 1
      end do
c
c     append an extension or version as appropriate
c
      if (k .eq. leng) then
         exist = .false.
         if (leng .ne. 0) then
            inquire (file=string(1:leng),exist=exist)
         end if
         if (.not. exist) then
            string = string(1:leng)//'.'//extension(1:lext)
            call version (string,status)
         end if
      else if (status .eq. 'new') then
         call version (string,status)
      end if
      return
      end

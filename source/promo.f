c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine promo  --  copywrite notice and version info  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "promo" writes a short message containing information
c     about the TINKER version number and the copyright notice
c
c
      subroutine promo
      implicit none
      include 'iounit.i'
c
c
c     print out the informational header message
c
      write (iout,10)
   10 format (/,' ',78('#'),
     &        /,' ',78('#'),
     &        /,' ##',74x,'##',
     &        /,' ##',13x,'TINKER  ---  Software Tools for',
     &           ' Molecular Design',13x,'##',
     &        /,' ##',74x,'##',
     &        /,' ##',24x,'Version 6.0   October 2011',24x,'##',
     &        /,' ##',74x,'##',
     &        /,' ##',15x,'Copyright (c)  Jay William Ponder',
     &           '  1990-2011',15x,'##',
     &        /,' ##',28x,'All Rights Reserved',27x,'##',
     &        /,' ##',74x,'##',
     &        /,' ',78('#'),
     &        /,' ',78('#'),/)
      return
      end

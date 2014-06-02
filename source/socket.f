c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module socket  --  socket communication control parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     runtyp      calculation type for passing socket information
c     cstep       current optimization or dynamics step number
c     cdt         current dynamics cumulative simulation time
c     cenergy     current potential energy from simulation
c     cdx         current gradient components along the x-axis
c     cdy         current gradient components along the y-axis
c     cdz         current gradient components along the z-axis
c     use_socket  logical flag governing use of external sockets
c     skt_init    logical flag to indicate socket initialization
c     skt_close   logical flag to indicate socket shutdown
c
c
      module socket
      implicit none
      integer runtyp
      integer cstep
      real*8 cdt
      real*8 cenergy
      real*8, allocatable :: cdx(:)
      real*8, allocatable :: cdy(:)
      real*8, allocatable :: cdz(:)
      logical use_socket
      logical skt_init
      logical skt_close
      save
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  socket.i  --  control parameters for socket communication  ##
c     ##                                                             ##
c     #################################################################
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
      integer runtyp,cstep
      real*8  cdt,cenergy
      real*8  cdx,cdy,cdz
      logical use_socket
      logical skt_init
      logical skt_close
      common /socket/ runtyp,cstep,cdt,cenergy,cdx(maxatm),cdy(maxatm),
     &                cdz(maxatm),use_socket,skt_init,skt_close

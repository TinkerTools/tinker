c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  iounit.i  --  Fortran input/output (I/O) unit numbers  ##
c     ##                                                         ##
c     #############################################################
c
c
c     iout    Fortran I/O unit for main output (default=6)
c     input   Fortran I/O unit for main input (default=5)
c
c
      integer iout,input
      common /iounit/ iout,input

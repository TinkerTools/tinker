c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  mdstuf.i  --  control of molecular dynamics trajectory  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     nfree       total number of degrees of freedom for a system
c     velsave     flag to save velocity vector components to a file
c     frcsave     flag to save force vector components to a file
c     uindsave    flag to save induced atomic dipoles to a file
c     integrate   type of molecular dynamics integration algorithm
c
c
      integer nfree
      logical velsave
      logical frcsave
      logical uindsave
      character*10 integrate
      common /mdstuf/ nfree,velsave,frcsave,uindsave,integrate

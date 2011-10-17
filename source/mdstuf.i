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
c     irest       steps between removal of COM inertia (0=no removal)
c     bmnsplit    splitting coefficient for Beeman integrator
c     dorest      logical flag to remove center of mass inertia
c     velsave     logical flag to save velocity vector components
c     frcsave     logical flag to save force vector components
c     uindsave    logical flag to save induced atomic dipoles
c     integrate   type of molecular dynamics integration algorithm
c
c
      integer nfree,irest
      integer bmnsplit
      logical dorest
      logical velsave
      logical frcsave
      logical uindsave
      character*10 integrate
      common /mdstuf/ nfree,irest,bmnsplit,dorest,velsave,frcsave,
     &                uindsave,integrate

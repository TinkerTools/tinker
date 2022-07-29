c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2022  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readcart  --  input of Cartesian coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readcart" gets a set of Cartesian coordinates from either
c     a formatted or binary disk file
c
c
      subroutine readcart (ixyz,first)
      use output
      implicit none
      integer ixyz
      logical first
c
c
c     get next coordinates set from formatted or binary file
c
      if (archive) then
         call readxyz (ixyz)
      else if (binary) then
         call readdcd (ixyz,first)
      end if
      return
      end

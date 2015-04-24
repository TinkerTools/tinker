c
c
c     ############################################################
c     ##                  COPYRIGHT (C) 2015                    ##
c     ##     by Mark Friedrichs, Lee-Ping Wang & Jay Ponder     ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  module openmm  --  opaque handle for OpenMM objects  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     ommhandle   opaque handle pointing to OpenMM data structure
c
c
      module openmm
      implicit none
      integer*8 ommhandle
      save
      end

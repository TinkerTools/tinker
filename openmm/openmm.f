c
c
c     ############################################################
c     ##                  COPYRIGHT (C) 2015                    ##
c     ##     by Mark Friedrichs, Lee-Ping Wang & Jay Ponder     ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module openmm  --  OpenMM-related objects and variables  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     ommHandle      opaque handle pointing to OpenMM data structure
c     cudaPrecision  string with CUDA precision (SINGLE, MIXED, DOUBLE)
c     ommPlatform    string with OpenMM platform type (REFERENCE, CUDA)
c     cudaDevice     string with names/numbers of the CUDA GPU cards
c
c
      module openmm
      implicit none
      integer*8 ommHandle
      character*6 cudaPrecision
      character*9 ommPlatform
      character*16 cudaDevice
      save
      end

c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module thrmint  --  thermodynamic integration variables  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     tinbcount      total number of blocks recorded so far
c     tinbsave       total number of blocks already written out
c     tinbtot        total number of blocks the schedule can record
c     tinstepavg     steps averaged into one dU/dlambda sample
c     tidedllist     dU/dlambda values saved within the current block
c     tilmdadedl     block averaged dU/dlambda in the order recorded
c     tilmdadedlstd  standard deviation within each averaged block
c     tilmdahist     main lambda in effect when each block was recorded
c
c
      module thrmint
      implicit none
      integer tinbcount
      integer tinbsave
      integer tinbtot
      integer tinstepavg
      real*8, allocatable :: tidedllist(:)
      real*8, allocatable :: tilmdadedl(:)
      real*8, allocatable :: tilmdadedlstd(:)
      real*8, allocatable :: tilmdahist(:)
      save
      end

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
c     tibin          index of the current lambda window, 1 to tinbin
c     tinbcount      total number of blocks recorded so far
c     tinbin         number of lambda windows spanning the [0,1] range
c     tinblock       number of averaged blocks in the current window
c     tinbsave       total number of blocks already written out
c     tinbtot        total number of blocks the schedule can record
c     tinequil       equilibration steps in the current window
c     tinstepavg     steps averaged into one dU/dlambda sample
c     tiwindow       total dynamics steps in the current window
c     tiwinend       last dynamics step belonging to each window
c     tieqratio      fraction of each window discarded as equilibration
c     tilmda         main lambda value of the current window
c     tidedllist     dU/dlambda values saved within the current block
c     tifraclist     fraction of the run spent in each lambda window
c     tilmdadedl     block averaged dU/dlambda in the order recorded
c     tilmdadedlstd  standard deviation within each averaged block
c     tilmdahist     main lambda in effect when each block was recorded
c     tilmdalist     main lambda value of each window in schedule order
c     tifile         name of the file holding the block averages
c
c
      module thrmint
      implicit none
      integer tibin
      integer tinbcount
      integer tinbin
      integer tinblock
      integer tinbsave
      integer tinbtot
      integer tinequil
      integer tinstepavg
      integer tiwindow
      integer, allocatable :: tiwinend(:)
      real*8 tieqratio
      real*8 tilmda
      real*8, allocatable :: tidedllist(:)
      real*8, allocatable :: tifraclist(:)
      real*8, allocatable :: tilmdadedl(:)
      real*8, allocatable :: tilmdadedlstd(:)
      real*8, allocatable :: tilmdahist(:)
      real*8, allocatable :: tilmdalist(:)
      character*240 tifile
      save
      end

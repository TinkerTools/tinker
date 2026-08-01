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
c     tinbin         number of lambda windows spanning the [0,1] range
c     tinblock       number of averaged blocks per lambda window
c     tinequil       equilibration steps per window, tiwindow*tieqratio
c     tinstepavg     steps averaged into one dU/dlambda sample
c     tiwindow       total dynamics steps per lambda window
c     tinbcount      number of blocks recorded so far in each window
c     tinbsave       number of blocks already written out per window
c     tieqratio      fraction of each window discarded as equilibration
c     tilmda         main lambda value of the current window
c     tidedllist     dU/dlambda values saved within the current block
c     tilmdalist     main lambda value of each window in schedule order
c     tilmdadedl     block averaged dU/dlambda for each lambda window
c     tilmdadedlstd  standard deviation within each averaged block
c     tifile         name of the file holding the block averages
c
c
      module thrmint
      implicit none
      integer tibin
      integer tinbin
      integer tinblock
      integer tinequil
      integer tinstepavg
      integer tiwindow
      integer, allocatable :: tinbcount(:)
      integer, allocatable :: tinbsave(:)
      real*8 tieqratio
      real*8 tilmda
      real*8, allocatable :: tidedllist(:)
      real*8, allocatable :: tilmdalist(:)
      real*8, allocatable :: tilmdadedl(:,:)
      real*8, allocatable :: tilmdadedlstd(:,:)
      character*240 tifile
      save
      end

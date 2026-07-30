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
c     tibin         index of the current lambda window, 0 to tinbin-1
c     tinbin        number of lambda windows spanning the [0,1] range
c     tinequil      equilibration steps per window, tiwindow*tieqratio
c     tinstepavg    steps averaged into one dU/dlambda sample
c     tiwindow      total dynamics steps per lambda window
c     tieqratio     fraction of each window discarded as equilibration
c     tilmda        main lambda value of the current window
c     tidedllist    dU/dlambda values saved within the current block
c     tilmdadedl    block averaged dU/dlambda for each lambda window
c     tilmdadedlstd standard deviation within each averaged block
c
c
      module thrmint
      implicit none
      integer tibin
      integer tinbin
      integer tinequil
      integer tinstepavg
      integer tiwindow
      real*8 tieqratio
      real*8 tilmda
      real*8, allocatable :: tidedllist(:)
      real*8, allocatable :: tilmdadedl(:,:)
      real*8, allocatable :: tilmdadedlstd(:,:)
      save
      end

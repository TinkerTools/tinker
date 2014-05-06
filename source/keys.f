c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module keys  --  contents of the keyword control file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     nkey      number of nonblank lines in the keyword file
c     keyline   contents of each individual keyword file line
c
c
      module keys
      use sizes
      implicit none
      integer nkey
      character*120 keyline(maxkey)
      save
      end

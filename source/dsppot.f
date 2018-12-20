c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module dsppot  --  dispersion interaction scale factors  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     dsp2scale   scale factor for 1-2 dispersion energy interactions
c     dsp3scale   scale factor for 1-3 dispersion energy interactions
c     dsp4scale   scale factor for 1-4 dispersion energy interactions
c     dsp5scale   scale factor for 1-5 dispersion energy interactions
c     use_dcorr   flag to use long range dispersion correction
c
c
      module dsppot
      implicit none
      real*8 dsp2scale
      real*8 dsp3scale
      real*8 dsp4scale
      real*8 dsp5scale
      logical use_dcorr
      save
      end

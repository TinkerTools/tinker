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
c     dhcut       start of vdw switch handing hal to Ewald dispersion
c     dhoff       end of vdw switch handing hal to Ewald dispersion
c     use_dcorr   flag to use long range dispersion correction
c     use_dhal    flag for an Ewald dispersion tail mapped from hal
c     dhdamp      flag for damping on the hal-mapped dispersion sites
c
c
      module dsppot
      implicit none
      real*8 dsp2scale
      real*8 dsp3scale
      real*8 dsp4scale
      real*8 dsp5scale
      real*8 dhcut,dhoff
      logical use_dcorr
      logical use_dhal
      logical dhdamp
      save
      end

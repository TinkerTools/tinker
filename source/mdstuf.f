c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module mdstuf  --  molecular dynamics trajectory controls  ##
c     ##                                                             ##
c     #################################################################
c
c
c     nfree       total number of degrees of freedom for a system
c     irest       steps between removal of COM inertia (0=no removal)
c     bmnmix      mixing coefficient for use with Beeman integrator
c     dorest      logical flag to remove center of mass inertia
c     velsave     logical flag to save velocity vector components
c     frcsave     logical flag to save force vector components
c     uindsave    logical flag to save induced atomic dipoles
c     integrate   type of molecular dynamics integration algorithm
c     nrespa      number of small time steps per large time step in 2-level RESPA
c     nres_bond
c     nres_tors
c     nres_short
c     mult_bond
c     mult_tors
c     mult_short

      module mdstuf
      implicit none
      real*8 mult_bond
      real*8 mult_tors
      real*8 mult_short      
      integer nfree
      integer irest
      integer bmnmix
      integer nrespa
      integer nres_bond
      integer nres_tors
      integer nres_short
      logical dorest
      logical velsave
      logical frcsave
      logical uindsave
      logical XO_RESPA
      logical XI_RESPA
      logical use_mpole_switch
      character*11 integrate
      save
      end

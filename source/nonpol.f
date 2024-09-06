c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module nonpol  --  nonpolar cavity & dispersion parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     epso      water oxygen eps for implicit dispersion term
c     epsh      water hydrogen eps for implicit dispersion term
c     rmino     water oxygen Rmin for implicit dispersion term
c     rminh     water hydrogen Rmin for implicit dispersion term
c     awater    water number density at standard temp & pressure
c     slevy     enthalpy-to-free energy scale factor for dispersion
c     shctd     HCT overlap scale factor for the dispersion integral
c     dspoff    radius offset for the start of dispersion integral
c
c     cavprb    probe radius for use in computing cavitation energy
c     solvprs   limiting microscopic solvent pressure value
c     surften   limiting macroscopic surface tension value
c     spcut     starting radius for solvent pressure tapering
c     spoff     cutoff radius for solvent pressure tapering
c     stcut     starting radius for surface tension tapering
c     stoff     cutoff radius for surface tension tapering
c     radcav    atomic radius of each atom for cavitation energy
c     raddsp    atomic radius of each atom for dispersion energy
c     epsdsp    vdw well depth of each atom for dispersion energy
c     cdsp      maximum dispersion energy for each atom
c
c
      module nonpol
      implicit none
      real*8 epso,epsh
      real*8 rmino,rminh
      real*8 awater,slevy
      real*8 shctd,dspoff
      parameter (epso=0.1100d0)
      parameter (epsh=0.0135d0)
      parameter (rmino=1.7025d0)
      parameter (rminh=1.3275d0)
      parameter (awater=0.033428d0)
      parameter (slevy=1.0d0)
      parameter (shctd=0.75d0)
      parameter (dspoff=1.056d0)
      real*8 cavprb
      real*8 solvprs
      real*8 surften
      real*8 spcut,spoff
      real*8 stcut,stoff
      real*8, allocatable :: radcav(:)
      real*8, allocatable :: raddsp(:)
      real*8, allocatable :: epsdsp(:)
      real*8, allocatable :: cdsp(:)
      save
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module limits  --  interaction taper & cutoff distances  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     vdwcut      cutoff distance for van der Waals interactions
c     repcut      cutoff distance for Pauli repulsion interactions
c     dispcut     cutoff distance for dispersion interactions
c     chgcut      cutoff distance for charge-charge interactions
c     dplcut      cutoff distance for dipole-dipole interactions
c     mpolecut    cutoff distance for atomic multipole interactions
c     ctrncut     cutoff distance for charge transfer interactions
c     vdwtaper    distance at which van der Waals switching begins
c     reptaper    distance at which Pauli repulsion switching begins
c     disptaper   distance at which dispersion switching begins
c     chgtaper    distance at which charge-charge switching begins
c     dpltaper    distance at which dipole-dipole switching begins
c     mpoletaper  distance at which atomic multipole switching begins
c     ctrntaper   distance at which charge transfer switching begins
c     ewaldcut    cutoff distance for real space Ewald electrostatics
c     dewaldcut   cutoff distance for real space Ewald dispersion
c     usolvcut    cutoff distance for dipole solver preconditioner
c     use_ewald   logical flag governing use of electrostatic Ewald
c     use_dewald  logical flag governing use of dispersion Ewald
c     use_lights  logical flag governing use of method of lights
c     use_list    logical flag governing use of any neighbor lists
c     use_vlist   logical flag governing use of van der Waals list
c     use_dlist   logical flag governing use of dispersion list
c     use_clist   logical flag governing use of charge list
c     use_mlist   logical flag governing use of multipole list
c     use_ulist   logical flag governing use of preconditioner list
c
c
      module limits
      implicit none
      real*8 vdwcut,repcut
      real*8 dispcut,chgcut
      real*8 dplcut,mpolecut
      real*8 ctrncut
      real*8 vdwtaper,reptaper
      real*8 disptaper,chgtaper
      real*8 dpltaper,mpoletaper
      real*8 ctrntaper
      real*8 ewaldcut,dewaldcut
      real*8 usolvcut
      logical use_ewald,use_dewald
      logical use_lights,use_list
      logical use_vlist,use_dlist
      logical use_clist,use_mlist
      logical use_ulist
      save
      end

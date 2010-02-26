c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  fields.i  --  molecular mechanics force field description  ##
c     ##                                                             ##
c     #################################################################
c
c
c     biotyp       force field atom type of each biopolymer type
c     forcefield   string used to describe the current forcefield
c
c
      integer biotyp
      character*20 forcefield
      common /fields/ biotyp(maxbio),forcefield

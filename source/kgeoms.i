c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  kgeoms.i  --  parameters for the geometrical restraints  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     xpfix      x-coordinate target for each restrained position
c     ypfix      y-coordinate target for each restrained position
c     zpfix      z-coordinate target for each restrained position
c     pfix       force constant and flat-well range for each position
c     dfix       force constant and target range for each distance
c     afix       force constant and target range for each angle
c     tfix       force constant and target range for each torsion
c     gfix       force constant and target range for each group distance
c     chir       force constant and target range for chiral centers
c     depth      depth of shallow Gaussian basin restraint
c     width      exponential width coefficient of Gaussian basin
c     rwall      radius of spherical droplet boundary restraint
c     npfix      number of position restraints to be applied
c     ipfix      atom number involved in each position restraint
c     kpfix      flags to use x-, y-, z-coordinate position restraints
c     ndfix      number of distance restraints to be applied
c     idfix      atom numbers defining each distance restraint
c     nafix      number of angle restraints to be applied
c     iafix      atom numbers defining each angle restraint
c     ntfix      number of torsional restraints to be applied
c     itfix      atom numbers defining each torsional restraint
c     ngfix      number of group distance restraints to be applied
c     igfix      group numbers defining each group distance restraint
c     nchir      number of chirality restraints to be applied
c     ichir      atom numbers defining each chirality restraint
c     use_basin  logical flag governing use of Gaussian basin
c     use_wall   logical flag governing use of droplet boundary
c
c
      integer npfix,ipfix
      integer kpfix
      integer ndfix,idfix
      integer nafix,iafix
      integer ntfix,itfix
      integer ngfix,igfix
      integer nchir,ichir
      real*8 xpfix,ypfix,zpfix
      real*8 pfix,dfix,afix
      real*8 tfix,gfix,chir
      real*8 depth,width,rwall
      logical use_basin,use_wall
      common /kgeoms/ xpfix(maxfix),ypfix(maxfix),zpfix(maxfix),
     &                pfix(2,maxfix),dfix(3,maxfix),afix(3,maxfix),
     &                tfix(3,maxfix),gfix(3,maxfix),chir(3,maxfix),
     &                depth,width,rwall,npfix,ipfix(maxfix),
     &                kpfix(3,maxfix),ndfix,idfix(2,maxfix),nafix,
     &                iafix(3,maxfix),ntfix,itfix(4,maxfix),ngfix,
     &                igfix(2,maxfix),nchir,ichir(4,maxfix),use_basin,
     &                use_wall

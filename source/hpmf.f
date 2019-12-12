c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module hpmf  --  hydrophobic potential of mean force term  ##
c     ##                                                             ##
c     #################################################################
c
c
c     rcarbon      radius of a carbon atom for use with HPMF
c     rwater       radius of a water molecule for use with HPMF
c     acsurf       surface area of a hydrophobic carbon atom
c     safact       constant for calculation of atomic surface area
c     tgrad        tanh slope (set very steep, default=100)
c     toffset      shift the tanh plot along the x-axis (default=6)
c     hpmfcut      cutoff distance for pairwise HPMF interactions
c     hd1,hd2,hd3  hydrophobic PMF well depth parameter
c     hc1,hc2,hc3  hydrophobic PMF well center point
c     hw1,hw2,hw3  reciprocal of the hydrophobic PMF well width
c
c     npmf         number of hydrophobic carbon atoms in the system
c     ipmf         number of the atom for each HPMF carbon atom site
c     rpmf         radius of each atom for use with hydrophobic PMF
c     acsa         SASA value for each hydrophobic PMF carbon atom
c
c
      module hpmf
      implicit none
      real*8 rcarbon,rwater
      real*8 acsurf,safact
      real*8 tgrad,toffset
      real*8 hpmfcut
      real*8 hd1,hd2,hd3
      real*8 hc1,hc2,hc3
      real*8 hw1,hw2,hw3
      parameter (rcarbon=1.7d0)
      parameter (rwater=1.4d0)
      parameter (acsurf=120.7628d0)
      parameter (safact=0.3516d0)
      parameter (tgrad=100.0d0)
      parameter (toffset=6.0d0)
      parameter (hpmfcut=11.0d0)
      parameter (hd1=-0.7308004860404441194d0)
      parameter (hd2=0.2001645051578760659d0)
      parameter (hd3=-0.0905499953418473502d0)
      parameter (hc1=3.8167879266271396155d0)
      parameter (hc2=5.4669162286016419472d0)
      parameter (hc3=7.1167694861385353278d0)
      parameter (hw1=1.6858993102248638341d0)
      parameter (hw2=1.3906405621629980285d0)
      parameter (hw3=1.5741657341338335385d0)
      integer npmf
      integer, allocatable :: ipmf(:)
      real*8, allocatable :: rpmf(:)
      real*8, allocatable :: acsa(:)
      save
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2011  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program apbstest  --  compute the APBS solvation energy  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "apbstest" computes the APBS Poisson-Boltzmann solvation energy
c     for an input structure
c
c
      program apbstest
      include 'energi.i'
      include 'inform.i'
      call initial
      call getxyz
      call mechanic
      call empole3
      verbose = .false.
      call esolv3
      write (*,10)  es
 10   format (/,' Total Solvation Energy',14x,f12.4)
      end

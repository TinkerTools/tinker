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
c     for a structure; requires an APBS-enabled TINKER version
c
c
      program apbstest
      use energi
      use inform
      use iounit
      implicit none
c
c
c     compute and report the APBS-based solvation energy
c
      call initial
      call getxyz
      call mechanic
      verbose = .false.
      call empole3
      call esolv3
      write (iout,10)  es
 10   format (/,' Total Solvation Energy',14x,f12.4)
      end

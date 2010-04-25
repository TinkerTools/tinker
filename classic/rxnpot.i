c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  rxnpot.i  --  specifics of reaction field functional form  ##
c     ##                                                             ##
c     #################################################################
c
c
c     rfsize    radius of reaction field sphere centered at origin
c     rfbulkd   bulk dielectric constant of reaction field continuum
c     rfterms   number of terms to use in reaction field summation
c
c
      integer rfterms
      real*8 rfsize,rfbulkd
      common /rxnpot/ rfsize,rfbulkd,rfterms

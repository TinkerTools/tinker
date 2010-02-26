c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  params.i  --  contents of force field parameter file  ##
c     ##                                                        ##
c     ############################################################
c
c
c     nprm      number of nonblank lines in the parameter file
c     prmline   contents of each individual parameter file line
c
c
      integer nprm
      character*120 prmline
      common /params/ nprm,prmline(maxprm)

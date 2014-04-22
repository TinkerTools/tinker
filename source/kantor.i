c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2014 by Chao Lv & Jay William Ponder  ##
c     ##                 All Rights Reserved                  ##
c     ##########################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  kantor.i  --  forcefield parameters for angle-torsions  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxnat   maximum number of angle-torsion parameter entries
c
c     atcon    torsional amplitude parameters for angle-torsion
c     kat      string of atom classes for angle-torsion terms
c
c
      integer maxnat
      parameter (maxnat=500)
      real*8 atcon
      character*16 kat
      common /kbntor/ atcon(6,maxnat),kat(maxnat)

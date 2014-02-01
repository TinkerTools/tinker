c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  angang.i  --  angle-angle terms in current structure  ##
c     ##                                                        ##
c     ############################################################
c
c
c     kaa       force constant for angle-angle cross terms
c     nangang   total number of angle-angle interactions
c     iaa       angle numbers used in each angle-angle term
c
c
      integer nangang,iaa
      real*8 kaa
      common /angang/ kaa(maxtors),nangang,iaa(2,maxtors)

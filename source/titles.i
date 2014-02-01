c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  titles.i  --  title for the current molecular system  ##
c     ##                                                        ##
c     ############################################################
c
c
c     ltitle   length in characters of the nonblank title string
c     title    title used to describe the current structure
c
c
      integer ltitle
      character*120 title
      common /titles/ ltitle,title

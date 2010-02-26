c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  ring.i  --  number and location of small ring structures  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nring3   total number of 3-membered rings in the system
c     iring3   numbers of the atoms involved in each 3-ring
c     nring4   total number of 4-membered rings in the system
c     iring4   numbers of the atoms involved in each 4-ring
c     nring5   total number of 5-membered rings in the system
c     iring5   numbers of the atoms involved in each 5-ring
c     nring6   total number of 6-membered rings in the system
c     iring6   numbers of the atoms involved in each 6-ring
c
c
      integer nring3,iring3
      integer nring4,iring4
      integer nring5,iring5
      integer nring6,iring6
      common /ring/ nring3,iring3(3,maxring),nring4,iring4(4,maxring),
     &              nring5,iring5(5,maxring),nring6,iring6(6,maxring)

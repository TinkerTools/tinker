c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  atoms.i  --  number, position and type of current atoms  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     x       current x-coordinate for each atom in the system
c     y       current y-coordinate for each atom in the system
c     z       current z-coordinate for each atom in the system
c     n       total number of atoms in the current system
c     type    atom type number for each atom in the system
c
c
      integer n,type
      real*8 x,y,z
      common /atoms/ x(maxatm),y(maxatm),z(maxatm),n,type(maxatm)

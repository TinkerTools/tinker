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
c     LPW added some stuff for storing atomic coordinates here
c     xstor   stored x-coordinate for each atom in the system
c     ystor   stored y-coordinate for each atom in the system
c     zstor   stored z-coordinate for each atom in the system
c     n       total number of atoms in the current system
c     type    atom type number for each atom in the system
c
c
      integer n,type
      real*8 x,y,z,xstor,ystor,zstor
      common /atoms/ x(maxatm),y(maxatm),z(maxatm),
     &               xstor(maxatm),ystor(maxatm),zstor(maxatm),
     &               n,type(maxatm)

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  improp.i  --  improper dihedrals in the current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     kprop    force constant values for improper dihedral angles
c     vprop    ideal improper dihedral angle value in degrees
c     niprop   total number of improper dihedral angles in the system
c     iiprop   numbers of the atoms in each improper dihedral angle
c
c
      integer niprop,iiprop
      real*8 kprop,vprop
      common /improp/ kprop(maxtors),vprop(maxtors),niprop,
     &                iiprop(4,maxtors)

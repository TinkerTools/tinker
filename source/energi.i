c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  energi.i  --  individual potential energy components  ##
c     ##                                                        ##
c     ############################################################
c
c
c     esum   total potential energy of the system
c     eb     bond stretch potential energy of the system
c     ea     angle bend potential energy of the system
c     eba    stretch-bend potential energy of the system
c     eub    Urey-Bradley potential energy of the system
c     eaa    angle-angle potential energy of the system
c     eopb   out-of-plane bend potential energy of the system
c     eopd   out-of-plane distance potential energy of the system
c     eid    improper dihedral potential energy of the system
c     eit    improper torsion potential energy of the system
c     et     torsional potential energy of the system
c     ept    pi-orbital torsion potential energy of the system
c     ebt    stretch-torsion potential energy of the system
c     ett    torsion-torsion potential energy of the system
c     ev     van der Waals potential energy of the system
c     ec     charge-charge potential energy of the system
c     ecd    charge-dipole potential energy of the system
c     ed     dipole-dipole potential energy of the system
c     em     atomic multipole potential energy of the system
c     ep     polarization potential energy of the system
c     er     reaction field potential energy of the system
c     es     solvation potential energy of the system
c     elf    metal ligand field potential energy of the system
c     eg     geometric restraint potential energy of the system
c     ex     extra term potential energy of the system
c
c
      real*8 esum,eb,ea,eba
      real*8 eub,eaa,eopb,eopd
      real*8 eid,eit,et,ept
      real*8 ebt,ett,ev,ec
      real*8 ecd,ed,em,ep,er
      real*8 es,elf,eg,ex
      common /energi/ esum,eb,ea,eba,eub,eaa,eopb,eopd,eid,eit,et,ept,
     &                ebt,ett,ev,ec,ecd,ed,em,ep,er,es,elf,eg,ex

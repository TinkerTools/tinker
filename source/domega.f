c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module domega  --  derivative components over torsions  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     tesum   total energy derivatives over torsions
c     teb     bond stretch derivatives over torsions
c     tea     angle bend derivatives over torsions
c     teba    stretch-bend derivatives over torsions
c     teub    Urey-Bradley derivatives over torsions
c     teaa    angle-angle derivatives over torsions
c     teopb   out-of-plane bend derivatives over torsions
c     teopd   out-of-plane distance derivatives over torsions
c     teid    improper dihedral derivatives over torsions
c     teit    improper torsion derivatives over torsions
c     tet     torsional derivatives over torsions
c     tept    pi-orbital torsion derivatives over torsions
c     tebt    stretch-torsion derivatives over torsions
c     teat    angle-torsion derivatives over torsions
c     tett    torsion-torsion derivatives over torsions
c     tev     van der Waals derivatives over torsions
c     tec     charge-charge derivatives over torsions
c     tecd    charge-dipole derivatives over torsions
c     ted     dipole-dipole derivatives over torsions
c     tem     atomic multipole derivatives over torsions
c     tep     polarization derivatives over torsions
c     ter     reaction field derivatives over torsions
c     tes     solvation derivatives over torsions
c     telf    metal ligand field derivatives over torsions
c     teg     geometric restraint derivatives over torsions
c     tex     extra energy term derivatives over torsions
c
c
      module domega
      use sizes
      implicit none
      real*8 tesum(maxrot)
      real*8 teb(maxrot)
      real*8 tea(maxrot)
      real*8 teba(maxrot)
      real*8 teub(maxrot)
      real*8 teaa(maxrot)
      real*8 teopb(maxrot)
      real*8 teopd(maxrot)
      real*8 teid(maxrot)
      real*8 teit(maxrot)
      real*8 tet(maxrot)
      real*8 tept(maxrot)
      real*8 tebt(maxrot)
      real*8 teat(maxrot)
      real*8 tett(maxrot)
      real*8 tev(maxrot)
      real*8 tec(maxrot)
      real*8 tecd(maxrot)
      real*8 ted(maxrot)
      real*8 tem(maxrot)
      real*8 tep(maxrot)
      real*8 ter(maxrot)
      real*8 tes(maxrot)
      real*8 telf(maxrot)
      real*8 teg(maxrot)
      real*8 tex(maxrot)
      save
      end

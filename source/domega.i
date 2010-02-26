c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #########################################################
c     ##                                                     ##
c     ##  domega.i  --  derivative components over torsions  ##
c     ##                                                     ##
c     #########################################################
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
      real*8 tesum,teb,tea,teba
      real*8 teub,teaa,teopb,teopd
      real*8 teid,teit,tet,tept
      real*8 tebt,tett,tev,tec
      real*8 tecd,ted,tem,tep,ter
      real*8 tes,telf,teg,tex
      common /domega/ tesum(maxrot),teb(maxrot),tea(maxrot),
     &                teba(maxrot),teub(maxrot),teaa(maxrot),
     &                teopb(maxrot),teopd(maxrot),teid(maxrot),
     &                teit(maxrot),tet(maxrot),tept(maxrot),
     &                tebt(maxrot),tett(maxrot),tev(maxrot),
     &                tec(maxrot),tecd(maxrot),ted(maxrot),
     &                tem(maxrot),tep(maxrot),ter(maxrot),
     &                tes(maxrot),telf(maxrot),teg(maxrot),
     &                tex(maxrot)

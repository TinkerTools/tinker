c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  deriv.i  --  Cartesian coordinate derivative components  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     desum   total energy Cartesian coordinate derivatives
c     deb     bond stretch Cartesian coordinate derivatives
c     dea     angle bend Cartesian coordinate derivatives
c     deba    stretch-bend Cartesian coordinate derivatives
c     deub    Urey-Bradley Cartesian coordinate derivatives
c     deaa    angle-angle Cartesian coordinate derivatives
c     deopb   out-of-plane bend Cartesian coordinate derivatives
c     deopd   out-of-plane distance Cartesian coordinate derivatives
c     deid    improper dihedral Cartesian coordinate derivatives
c     deit    improper torsion Cartesian coordinate derivatives
c     det     torsional Cartesian coordinate derivatives
c     dept    pi-orbital torsion Cartesian coordinate derivatives
c     debt    stretch-torsion Cartesian coordinate derivatives
c     dett    torsion-torsion Cartesian coordinate derivatives
c     dev     van der Waals Cartesian coordinate derivatives
c     dec     charge-charge Cartesian coordinate derivatives
c     decd    charge-dipole Cartesian coordinate derivatives
c     ded     dipole-dipole Cartesian coordinate derivatives
c     dem     multipole Cartesian coordinate derivatives
c     dep     polarization Cartesian coordinate derivatives
c     der     reaction field Cartesian coordinate derivatives
c     des     solvation Cartesian coordinate derivatives
c     delf    metal ligand field Cartesian coordinate derivatives
c     deg     geometric restraint Cartesian coordinate derivatives
c     dex     extra energy term Cartesian coordinate derivatives
c
c
      real*8 desum,deb,dea,deba
      real*8 deub,deaa,deopb,deopd
      real*8 deid,det,dept,deit
      real*8 debt,dett,dev,dec
      real*8 decd,ded,dem,dep,der
      real*8 des,delf,deg,dex
      common /deriv/ desum(3,maxatm),deb(3,maxatm),dea(3,maxatm),
     &               deba(3,maxatm),deub(3,maxatm),deaa(3,maxatm),
     &               deopb(3,maxatm),deopd(3,maxatm),deid(3,maxatm),
     &               deit(3,maxatm),det(3,maxatm),dept(3,maxatm),
     &               debt(3,maxatm),dett(3,maxatm),dev(3,maxatm),
     &               dec(3,maxatm),decd(3,maxatm),ded(3,maxatm),
     &               dem(3,maxatm),dep(3,maxatm),der(3,maxatm),
     &               des(3,maxatm),delf(3,maxatm),deg(3,maxatm),
     &               dex(3,maxatm)

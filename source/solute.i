c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  solute.i  --  parameters for continuum solvation models  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     rsolv     atomic radius of each atom for continuum solvation
c     asolv     atomic surface area solvation parameters
c     rborn     Born radius of each atom for GB/SA solvation
c     drb       solvation derivatives with respect to Born radii
c     drbp      GK polarization derivatives with respect to Born radii
c     drobc     chain rule term for Onufriev-Bashford-Case radii
c     doffset   dielectric offset to continuum solvation atomic radii
c     p1        single-atom scale factor for analytical Still radii
c     p2        1-2 interaction scale factor for analytical Still radii
c     p3        1-3 interaction scale factor for analytical Still radii
c     p4        nonbonded scale factor for analytical Still radii
c     p5        soft cutoff parameter for analytical Still radii
c     gpol      polarization self-energy values for each atom
c     shct      overlap scale factors for Hawkins-Cramer-Truhlar radii
c     aobc      alpha values for Onufriev-Bashford-Case radii
c     bobc      beta values for Onufriev-Bashford-Case radii
c     gobc      gamma values for Onufriev-Bashford-Case radii
c     vsolv     atomic volume of each atom for use with ACE
c     wace      "omega" values for atom class pairs for use with ACE
c     s2ace     "sigma^2" values for atom class pairs for use with ACE
c     uace      "mu" values for atom class pairs for use with ACE
c     solvtyp   type of continuum solvation energy model in use
c     borntyp   method to be used for the Born radius computation
c
c
      real*8 rsolv,asolv
      real*8 rborn,drb,drbp
      real*8 drobc,doffset
      real*8 p1,p2,p3,p4,p5
      real*8 gpol,shct
      real*8 aobc,bobc,gobc
      real*8 vsolv,wace
      real*8 s2ace,uace
      character*8 solvtyp
      character*8 borntyp
      common /solute/ rsolv(maxatm),asolv(maxatm),rborn(maxatm),
     &                drb(maxatm),drbp(maxatm),drobc(maxatm),
     &                doffset,p1,p2,p3,p4,p5,gpol(maxatm),shct(maxatm),
     &                aobc(maxatm),bobc(maxatm),gobc(maxatm),
     &                vsolv(maxatm),wace(maxclass,maxclass),
     &                s2ace(maxclass,maxclass),uace(maxclass,maxclass),
     &                solvtyp,borntyp

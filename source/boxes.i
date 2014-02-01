c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  boxes.i  --  parameters for periodic boundary conditions  ##
c     ##                                                            ##
c     ################################################################
c
c
c     xbox        length of a-axis of periodic box in Angstroms
c     ybox        length of b-axis of periodic box in Angstroms
c     zbox        length of c-axis of periodic box in Angstroms
c     alpha       angle between b- and c-axes of box in degrees
c     beta        angle between a- and c-axes of box in degrees
c     gamma       angle between a- and b-axes of box in degrees
c     xbox2       half of the a-axis length of periodic box
c     ybox2       half of the b-axis length of periodic box
c     zbox2       half of the c-axis length of periodic box
c     box34       three-fourths axis length of truncated octahedron
c     lvec        real space lattice vectors as matrix rows
c     recip       reciprocal lattice vectors as matrix columns
c     volbox      volume in Ang**3 of the periodic box
c     beta_sin    sine of the beta periodic box angle
c     beta_cos    cosine of the beta periodic box angle
c     gamma_sin   sine of the gamma periodic box angle
c     gamma_cos   cosine of the gamma periodic box angle
c     beta_term   term used in generating triclinic box
c     gamma_term  term used in generating triclinic box
c     orthogonal  flag to mark periodic box as orthogonal
c     monoclinic  flag to mark periodic box as monoclinic
c     triclinic   flag to mark periodic box as triclinic
c     octahedron  flag to mark box as truncated octahedron
c     spacegrp    space group symbol for the unitcell type
c
c
      real*8 xbox,ybox,zbox
      real*8 alpha,beta,gamma
      real*8 xbox2,ybox2,zbox2
      real*8 box34,volbox
      real*8 lvec,recip
      real*8 beta_sin,beta_cos
      real*8 gamma_sin,gamma_cos
      real*8 beta_term,gamma_term
      logical orthogonal,monoclinic
      logical triclinic,octahedron
      character*10 spacegrp
      common /boxes/ xbox,ybox,zbox,alpha,beta,gamma,xbox2,ybox2,zbox2,
     &               box34,lvec(3,3),recip(3,3),volbox,beta_sin,
     &               beta_cos,gamma_sin,gamma_cos,beta_term,gamma_term,
     &               orthogonal,monoclinic,triclinic,octahedron,spacegrp

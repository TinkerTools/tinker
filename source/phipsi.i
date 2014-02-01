c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  phipsi.i  --  phi-psi-omega-chi angles for a protein  ##
c     ##                                                        ##
c     ############################################################
c
c
c     phi      value of the phi angle for each amino acid residue
c     psi      value of the psi angle for each amino acid residue
c     omega    value of the omega angle for each amino acid residue
c     chi      values of the chi angles for each amino acid residue
c     chiral   chirality of each amino acid residue (1=L, -1=D)
c     disulf   residue joined to each residue via a disulfide link
c
c
      integer chiral,disulf
      real*8 phi,psi,omega,chi
      common /phipsi/ phi(maxres),psi(maxres),omega(maxres),
     &                chi(4,maxres),chiral(maxres),disulf(maxres)

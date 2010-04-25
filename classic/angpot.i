c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  angpot.i  --  specifics of angle bend functional forms  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     angunit    convert angle bending energy to kcal/mole
c     stbnunit   convert stretch-bend energy to kcal/mole
c     aaunit     convert angle-angle energy to kcal/mole
c     opbunit    convert out-of-plane bend energy to kcal/mole
c     opdunit    convert out-of-plane distance energy to kcal/mole
c     cang       cubic coefficient in angle bending potential
c     qang       quartic coefficient in angle bending potential
c     pang       quintic coefficient in angle bending potential
c     sang       sextic coefficient in angle bending potential
c     copb       cubic coefficient in out-of-plane bend potential
c     qopb       quartic coefficient in out-of-plane bend potential
c     popb       quintic coefficient in out-of-plane bend potential
c     sopb       sextic coefficient in out-of-plane bend potential
c     copd       cubic coefficient in out-of-plane distance potential
c     qopd       quartic coefficient in out-of-plane distance potential
c     popd       quintic coefficient in out-of-plane distance potential
c     sopd       sextic coefficient in out-of-plane distance potential
c     angtyp     type of angle bending function for each bond angle
c     opbtyp     type of out-of-plane bend potential energy function
c
c
      real*8 angunit,stbnunit,aaunit
      real*8 opbunit,opdunit
      real*8 cang,qang,pang,sang
      real*8 copb,qopb,popb,sopb
      real*8 copd,qopd,popd,sopd
      character*8 angtyp,opbtyp
      common /angpot/ angunit,stbnunit,aaunit,opbunit,opdunit,cang,
     &                qang,pang,sang,copb,qopb,popb,sopb,copd,qopd,
     &                popd,sopd,angtyp(maxang),opbtyp

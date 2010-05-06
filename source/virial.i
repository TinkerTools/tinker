c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  virial.i  --  components of internal virial tensor  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     vir    total internal virial Cartesian tensor components
c
c
      real*8 vir, viri
      common /virial/ vir(3,3),viri(3,3)
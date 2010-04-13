c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  vdw.i  --  van der Waals parameters for current structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     radmin     minimum energy distance for each atom class pair
c     epsilon    well depth parameter for each atom class pair
c     radmin4    minimum energy distance for 1-4 interaction pairs
c     epsilon4   well depth parameter for 1-4 interaction pairs
c     radhbnd    minimum energy distance for hydrogen bonding pairs
c     epshbnd    well depth parameter for hydrogen bonding pairs
c     kred       value of reduction factor parameter for each atom
c     ired       attached atom from which reduction factor is applied
c     nvdw       total number van der Waals active sites in the system
c     ivdw       number of the atom for each van der Waals active site
c     jvdw       type or class index into vdw parameters for each atom
c     nvt        number of distinct van der Waals types in the system
c     ivt        number of each distinct vdw type/class in the system
c     jvt        frequency of each vdw type or class in the system
c
c
      integer ired,nvdw
      integer ivdw,jvdw
      integer nvt,ivt,jvt
      real*8 radmin,epsilon
      real*8 radmin4,epsilon4
      real*8 radhbnd,epshbnd,kred
      common /vdw/ radmin(maxclass,maxclass),
     &             epsilon(maxclass,maxclass),
     &             radmin4(maxclass,maxclass),
     &             epsilon4(maxclass,maxclass),
     &             radhbnd(maxclass,maxclass),
     &             epshbnd(maxclass,maxclass),
     &             kred(maxatm),ired(maxatm),
     &             nvdw,ivdw(maxatm),jvdw(maxatm),
     &             nvt,ivt(maxatm),jvt(maxatm)

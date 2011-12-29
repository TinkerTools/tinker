c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2011  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  uprior.i  --  previous values of induced dipole moments  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxualt   maximum number of sets of induced dipoles to save
c
c     baspc     coefficients for always stable predictor-corrector
c     udalt     prior values for induced dipoles at each site
c     upalt     prior values for induced dipoles in energy field
c     nualt     number of prior sets of induced dipoles in storage
c     use_aspc  flag to control use of ASPC for induced dipoles
c
c
      integer maxualt
      parameter (maxualt=6)
      integer nualt
      real*8 baspc
      real*8, pointer :: udalt(:,:,:)
      real*8, pointer :: upalt(:,:,:)
      logical use_aspc
      common /uprior/ baspc(maxualt),udalt,upalt,nualt,use_aspc

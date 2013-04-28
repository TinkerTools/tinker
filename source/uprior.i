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
c     bpred     coefficients for induced dipole predictor polynomial
c     udalt     prior values for induced dipoles at each site
c     upalt     prior values for induced dipoles in energy field
c     usalt     prior values for induced dipoles for PB/GK solvation
c     upsalt    prior values for induced dipoles in PB/GK energy field
c     nualt     number of prior sets of induced dipoles in storage
c     use_pred  flag to control use of induced dipole prediction
c
c
      integer maxualt
      parameter (maxualt=6)
      integer nualt
      real*8 bpred
      real*8, pointer :: udalt(:,:,:)
      real*8, pointer :: upalt(:,:,:)
      real*8, pointer :: usalt(:,:,:)
      real*8, pointer :: upsalt(:,:,:)
      logical use_pred
      common /uprior/ bpred(maxualt),udalt,upalt,usalt,upsalt,nualt,
     &                use_pred

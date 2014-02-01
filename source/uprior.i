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
c     gear      coefficients for Gear predictor binomial method
c     aspc      coefficients for always stable predictor-corrector
c     bpred     coefficients for induced dipole predictor polynomial
c     bpredp    coefficients for predictor polynomial in energy field
c     bpreds    coefficients for predictor for PB/GK solvation
c     bpredps   coefficients for predictor in PB/GK energy field
c     udalt     prior values for induced dipoles at each site
c     upalt     prior values for induced dipoles in energy field
c     usalt     prior values for induced dipoles for PB/GK solvation
c     upsalt    prior values for induced dipoles in PB/GK energy field
c     nualt     number of prior sets of induced dipoles in storage
c     use_pred  flag to control use of induced dipole prediction
c     polpred   type of predictor polynomial (Gear, ASPC or LSQR)
c
c
      integer maxualt
      parameter (maxualt=7)
      integer nualt
      real*8 gear,aspc
      real*8 bpred,bpredp
      real*8 bpreds,bpredps
      real*8, pointer :: udalt(:,:,:)
      real*8, pointer :: upalt(:,:,:)
      real*8, pointer :: usalt(:,:,:)
      real*8, pointer :: upsalt(:,:,:)
      logical use_pred
      character*4 polpred
      common /uprior/ gear(maxualt),aspc(maxualt),bpred(maxualt),
     &                bpredp(maxualt),bpreds(maxualt),bpredps(maxualt),
     &                udalt,upalt,usalt,upsalt,nualt,use_pred,polpred

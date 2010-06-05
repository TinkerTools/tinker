c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2010  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  vibs.i  --  components of iterative vibrational analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxroot    maximum number of roots converged simultaneously
c     maxbasis   maximum number of basis vectors for diagonalization
c
c     phi        trial vectors for iterative vibrational analysis
c     phik       alternate vectors for iterative vibrational analysis
c     pwork      temporary work array for eigenvector transformation
c
c
      integer maxroot
      integer maxbasis
      parameter (maxroot=50)
      parameter (maxbasis=3*maxroot)
      real*8 phi,phik,pwork
      common /vibs/ phi(maxvar,maxbasis),phik(maxvar,maxbasis),
     &              pwork(maxvar,maxbasis)

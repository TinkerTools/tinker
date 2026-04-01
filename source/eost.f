c
c
c     ##################################################################
c     ##     COPYRIGHT (C) 2026 by  Moses Chung, Lianqing Zheng,      ##
c     ##        Jayvee Abella, Michael Schnieders, Pengyu Ren,        ##
c     ##                     Wei Yang, Jay Ponder                     ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine eost - orthogonal space tempering algorithm  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "eost" calculates free energies using the orthogonal space
c     tempering algorithm using biasing gaussians
c
c
      subroutine eost
      use dlmda
      use mutant
      use ost
      implicit none
c
c
c     xxx
c
      return 
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine lmdachain - chain rule for main lambda deriv  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "lmdachain" applies chain rule to the main lambda derivative
c     to compute the global lambda derivative of energy, energy^2,
c     force, and virial
c
c
      subroutine lmdachain
      use atoms
      use dlmda
      use mutant
      use ost
      implicit none
      integer i,j
c
c
c     apply chain rule
c
      deplmda2 = deplmda2 * dplambda*dplambda
     &           + deplmda * d2plambda
      deplmda = deplmda * dplambda
      devlmda2 = devlmda2 * dvlambda*dvlambda
     &           + devlmda * d2vlambda
      devlmda = devlmda * dvlambda
      demlmda2 = demlmda2 * delambda*delambda
     &           + demlmda * d2elambda
      demlmda = demlmda * delambda
      do i = 1, n
         do j = 1, 3
            dldep(j,i) = dldep(j,i) * dplambda
            dldem(j,i) = dldem(j,i) * delambda
            dldev(j,i) = dldev(j,i) * dvlambda
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            dldepvir(j,i) = dldepvir(j,i) * dplambda
            dldemvir(j,i) = dldemvir(j,i) * delambda
            dldevvir(j,i) = dldevvir(j,i) * dvlambda
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine mapsublmda - map from lambda to sublambda  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "mapsublmda" maps from main lambda to sublambdas
c
c
      subroutine mapsublmda
      use dlmda
      use mutant
      use ost
      implicit none
      real*8 x
      real*8 taper
      real*8 dtaper
      real*8 d2taper
      character*6 mode
c
c
c     map from lambda to sublambdas
c
      mode = 'OSTPOL'
      call sublmdataper (mode,ostlambda,taper,dtaper,d2taper)
      plambda = 1.0d0 - taper
      dplambda = -dtaper
      d2plambda = -d2taper
      mode = 'OSTELE'
      call sublmdataper (mode,ostlambda,taper,dtaper,d2taper)
      elambda = 1.0d0 - taper
      delambda = -dtaper
      d2elambda = -d2taper
      mode = 'OSTVDW'
      call sublmdataper (mode,ostlambda,taper,dtaper,d2taper)
      vlambda = 1.0d0 - taper
      dvlambda = -dtaper
      d2vlambda = -d2taper
c
c     set flags to compute polarization lambda derivative
c
      use_pol4i = (plambda .le. ostplmda1)
      use_pol4f = (plambda .ge. ostplmda0)
      return
      end
c
c
c     ######################################################
c     ##                                                  ##
c     ##  subroutine sublmdataper - tapers the sublambda  ##
c     ##                                                  ##
c     ######################################################
c
c
c     "sublmdataper" tapers the mapping from main lambda to
c     sublambda at the endpoints
c
c
      subroutine sublmdataper (mode,x,taper,dtaper,d2taper)
      use shunt
      implicit none
      real*8 taper
      real*8 dtaper
      real*8 d2taper
      real*8 x,x2,x3
      real*8 x4,x5
      character*6 mode
c
c
c     get taper coefficients from existing Tinker switch routine
c
      call switch (mode)
c
c     return if outside switching window
c
      if (x .le. cut) then
         taper = 1.0d0
         dtaper = 0.0d0
         d2taper = 0.0d0
         return
      else if (x .ge. off) then
         taper = 0.0d0
         dtaper = 0.0d0
         d2taper = 0.0d0
         return
      end if
c
c     compute the quintic taper and derivative
c
      x2 = x*x
      x3 = x2*x
      x4 = x2*x2
      x5 = x2*x3
      taper = c5*x5 + c4*x4 + c3*x3 + c2*x2 + c1*x + c0
      dtaper = 5.0d0*c5*x4 + 4.0d0*c4*x3 + 3.0d0*c3*x2
     &            + 2.0d0*c2*x + c1
      d2taper = 20.0d0*c5*x3 + 12.0d0*c4*x2 + 6.0d0*c3*x
     &             + 2.0d0*c2
      return
      end

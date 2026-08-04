c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine mapsublmda -- map from lambda to sublambda  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "mapsublmda" maps from the main lambda value "lmda" to the
c     electrostatic, polarization and van der Waals sublambdas
c
c
      subroutine mapsublmda (lmda)
      use dlmda
      use mutant
      implicit none
      real*8 lmda
      real*8 taper
      real*8 dtaper
      real*8 d2taper
c
c
c     the staged relative schedule maps the sublambdas its own way
c
      if (use_relstage) then
         call maprelstage (lmda)
         return
      end if
c
c     map from lambda to sublambdas
c
      if (use_plmdamap) then
         if (plmdamap .eq. 'EXP') then
            call sublmdaexp (lmda,plmdaexp,plambda,
     &                       dpldlmda,d2pldlmda2)
         else if (plmdamap .eq. 'INV') then
            call sublmdainvpower (lmda,plmdainvn,plmdainveps,plambda,
     &                            dpldlmda,d2pldlmda2)
         else
            call quintaper (lmda,qntplmda0,qntplmda1,
     &                      taper,dtaper,d2taper)
            plambda = 1.0d0 - taper
            dpldlmda = -dtaper
            d2pldlmda2 = -d2taper
         end if
      end if
      if (use_elmdamap) then
         if (elmdamap .eq. 'EXP') then
            call sublmdaexp (lmda,elmdaexp,elambda,
     &                       deldlmda,d2eldlmda2)
         else if (elmdamap .eq. 'INV') then
            call sublmdainvpower (lmda,elmdainvn,elmdainveps,elambda,
     &                            deldlmda,d2eldlmda2)
         else
            call quintaper (lmda,qntelmda0,qntelmda1,
     &                      taper,dtaper,d2taper)
            elambda = 1.0d0 - taper
            deldlmda = -dtaper
            d2eldlmda2 = -d2taper
         end if
      end if
      if (use_vlmdamap) then
         if (vlmdamap .eq. 'EXP') then
            call sublmdaexp (lmda,vlmdaexp,vlambda,
     &                       dvldlmda,d2vldlmda2)
         else if (vlmdamap .eq. 'INV') then
            call sublmdainvpower (lmda,vlmdainvn,vlmdainveps,vlambda,
     &                            dvldlmda,d2vldlmda2)
         else
            call quintaper (lmda,qntvlmda0,qntvlmda1,
     &                      taper,dtaper,d2taper)
            vlambda = 1.0d0 - taper
            dvldlmda = -dtaper
            d2vldlmda2 = -d2taper
         end if
      end if
c
c     select the dual topology endpoints needed by QNT maps; both
c     endpoints remain active on and within each switching window
c
      if (use_elmdamap .and. elmdamap.eq.'QNT') then
         use_ele4i = (lmda .le. qntelmda1)
         use_ele4f = (lmda .ge. qntelmda0)
      end if
      if (use_plmdamap .and. plmdamap.eq.'QNT') then
         use_pol4i = (lmda .le. qntplmda1)
         use_pol4f = (lmda .ge. qntplmda0)
      end if
      if (use_vlmdamap .and. vlmdamap.eq.'QNT') then
         use_vdw4i = (lmda .le. qntvlmda1)
         use_vdw4f = (lmda .ge. qntvlmda0)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine maprelstage -- staged relative lambda mapping  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "maprelstage" maps the main lambda "lmda" onto the sublambdas for
c     the staged relative free energy schedule, in which the two ligands
c     are discharged and recharged one at a time while the van der Waals
c     terms morph between them in the middle window
c
c        lmda > relstg2lmda0    ligand 1 electrostatics, weight 0 -> 1
c        lmda < relstg1lmda1    ligand 0 electrostatics, weight 1 -> 0
c        otherwise              both ligands electrostatically decoupled
c
c     the electrostatic weight is the quintic taper, so the mixing
c     exponent used by the energy routines is one and the whole main
c     lambda chain rule is carried by "deldlmda"; polarization stages
c     with the multipoles while van der Waals morphs over its own window
c
c
      subroutine maprelstage (lmda)
      use dlmda
      use mutant
      implicit none
      real*8 lmda
      real*8 w
      real*8 taper
      real*8 dtaper
      real*8 d2taper
c
c
c     find the leg the main lambda sits in and its weight
c
      if (lmda .gt. relstg2lmda0) then
         call quintaper (lmda,relstg2lmda0,relstg2lmda1,
     &                   taper,dtaper,d2taper)
         relstage = 'LIG1'
         w = 1.0d0 - taper
         deldlmda = -dtaper
         d2eldlmda2 = -d2taper
      else if (lmda .lt. relstg1lmda1) then
         call quintaper (lmda,relstg1lmda0,relstg1lmda1,
     &                   taper,dtaper,d2taper)
         relstage = 'LIG0'
         w = taper
         deldlmda = dtaper
         d2eldlmda2 = d2taper
      else
         relstage = 'VDWM'
         w = 0.0d0
         deldlmda = 0.0d0
         d2eldlmda2 = 0.0d0
      end if
c
c     numerical guard for w = 1.0d0 - taper
c
      elambda = min(1.0d0,max(0.0d0,w))
      if (elambda .eq. 0.0d0)  relstage = 'VDWM'
      relstagemix = (elambda .gt. 0.0d0 .and. elambda .lt. 1.0d0)
c
c     polarization stages with the multipoles, same states same weight
c
      plambda = elambda
      dpldlmda = deldlmda
      d2pldlmda2 = d2eldlmda2
c
c     van der Waals morphs between the two ligands
c
      call quintaper (lmda,qntvlmda0,qntvlmda1,
     &                taper,dtaper,d2taper)
      vlambda = 1.0d0 - taper
      dvldlmda = -dtaper
      d2vldlmda2 = -d2taper
c
c     the staged routines branch on "relstagemix" instead of the
c     quantized endpoint flags, so leave the flags fully open
c
      use_ele4i = .true.
      use_ele4f = .true.
      use_pol4i = .true.
      use_pol4f = .true.
      use_vdw4i = .true.
      use_vdw4f = .true.
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine quintaper -- quintic taper over a window  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "quintaper" evaluates the quintic switching polynomial and its
c     first two derivatives over an arbitrary window, the taper falling
c     smoothly from one at "cut" to zero at "off"; the polynomial is
c     evaluated in the reduced coordinate (x-cut)/(off-cut) rather than
c     from the coefficients of "switch", which carry 1/(off-cut)**5 and
c     lose most of their precision once the window is narrow
c
c
      subroutine quintaper (x,cut,off,taper,dtaper,d2taper)
      implicit none
      real*8 x,cut,off
      real*8 taper
      real*8 dtaper
      real*8 d2taper
      real*8 rinv,u,u2,v
c
c
c     return if outside the switching window
c
      dtaper = 0.0d0
      d2taper = 0.0d0
      if (x .le. cut) then
         taper = 1.0d0
         return
      else if (x .ge. off) then
         taper = 0.0d0
         return
      end if
c
c     compute the quintic taper and its derivatives
c
      rinv = 1.0d0 / (off-cut)
      u = (x-cut) * rinv
      u2 = u * u
      v = 1.0d0 - u
      taper = 1.0d0 - u2*u*(10.0d0-15.0d0*u+6.0d0*u2)
      dtaper = -30.0d0 * u2 * v * v * rinv
      d2taper = -60.0d0 * u * v * (1.0d0-2.0d0*u) * rinv * rinv
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sublmdaexp -- exponential sublambda mapping  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sublmdaexp" maps from main lambda to a sublambda using a
c     power law and returns the first and second lambda derivatives
c
c
      subroutine sublmdaexp (x,exponent,lmda,dlmda,d2lmda)
      implicit none
      integer exponent
      real*8 x
      real*8 lmda
      real*8 dlmda
      real*8 d2lmda
      real*8 expnt
c
c
c     initialize values and handle endpoints explicitly
c
      expnt = dble(exponent)
      if (x .le. 0.0d0) then
         lmda = 0.0d0
         if (exponent .eq. 1) then
            dlmda = 1.0d0
            d2lmda = 0.0d0
         else if (exponent .eq. 2) then
            dlmda = 0.0d0
            d2lmda = 2.0d0
         else
            dlmda = 0.0d0
            d2lmda = 0.0d0
         end if
         return
      else if (x .ge. 1.0d0) then
         lmda = 1.0d0
         dlmda = expnt
         d2lmda = expnt * (expnt-1.0d0)
         return
      end if
c
c     compute lambda^exponent and its derivatives
c
      lmda = x**exponent
      dlmda = expnt * x**(exponent-1)
      if (exponent .eq. 1) then
         d2lmda = 0.0d0
      else
         d2lmda = expnt * (expnt-1.0d0) * x**(exponent-2)
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine sublmdainvpower -- inverse-power mapping  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "sublmdainvpower" maps from main lambda to sublambda using a
c     shifted inverse-power law and returns first and second lambda
c     derivatives
c
c
      subroutine sublmdainvpower (x,n,eps,lmda,dlmda,d2lmda)
      implicit none
      integer n
      real*8 x
      real*8 eps
      real*8 lmda
      real*8 dlmda
      real*8 d2lmda
      real*8 xval
      real*8 shift
      real*8 power
      real*8 root0
      real*8 denom
      real*8 base
c
c
c     set a bounded coordinate and handle the identity map directly
c
      xval = x
      if (xval .lt. 0.0d0)  xval = 0.0d0
      if (xval .gt. 1.0d0)  xval = 1.0d0
      if (n .le. 1) then
         lmda = xval
         dlmda = 1.0d0
         d2lmda = 0.0d0
         return
      end if
c
c     compute normalized shifted inverse-power map
c
      shift = eps
      if (shift .le. 0.0d0)  shift = 0.1d0
      power = 1.0d0 / dble(n)
      root0 = shift**power
      denom = (1.0d0+shift)**power - root0
      base = xval + shift
      lmda = (base**power-root0) / denom
      dlmda = power * base**(power-1.0d0) / denom
      d2lmda = power * (power-1.0d0)
     &           * base**(power-2.0d0) / denom
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine lmdachain -- chain rule for main lambda deriv  ##
c     ##                                                            ##
c     ################################################################
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
      implicit none
      integer i,j
c
c
c     apply chain rule for derivative of energy wrt global lambda
c
      d2epdl2 = d2epdl2 * dpldlmda*dpldlmda + depdl * d2pldlmda2
      depdl = depdl * dpldlmda
      d2evdl2 = d2evdl2 * dvldlmda*dvldlmda + devdl * d2vldlmda2
      devdl = devdl * dvldlmda
      d2emdl2 = d2emdl2 * deldlmda*deldlmda + demdl * d2eldlmda2
      demdl = demdl * deldlmda
      do i = 1, n
         do j = 1, 3
            dfpdl(j,i) = dfpdl(j,i) * dpldlmda
            dfmdl(j,i) = dfmdl(j,i) * deldlmda
            dfvdl(j,i) = dfvdl(j,i) * dvldlmda
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            depvirdl(j,i) = depvirdl(j,i) * dpldlmda
            demvirdl(j,i) = demvirdl(j,i) * deldlmda
            devvirdl(j,i) = devvirdl(j,i) * dvldlmda
         end do
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine avgstd -- average and std deviation kernel  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "avgstd" computes the average and population standard deviation
c     of the "count" samples of "list" starting at index "begin"
c
c
      subroutine avgstd (list,begin,count,avg,std)
      implicit none
      integer i
      integer begin
      integer count
      real*8 avg,std
      real*8 delta
      real*8 list(*)
c
c
c     an empty sample range has no average or deviation
c
      if (count .lt. 1) then
         avg = 0.0d0
         std = 0.0d0
         return
      end if
c
c     compute the average from the collected samples
c
      avg = 0.0d0
      do i = begin, begin+count-1
         avg = avg + list(i)
      end do
      avg = avg / dble(count)
c
c     compute the population standard deviation of the samples
c
      std = 0.0d0
      do i = begin, begin+count-1
         delta = list(i) - avg
         std = std + delta*delta
      end do
      std = sqrt(std/dble(count))
      return
      end

c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine refreshsublmda  --  refresh active sublambdas  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "refreshsublmda" maps the main lambda onto the component
c     sublambdas and installs the resulting absolute-topology
c     electrostatic state
c
c
      subroutine refreshsublmda
      use dlmda
      use mutant
      implicit none
c
c
c     return when no main lambda drives the sublambda maps
c
      if (.not. use_mainlmda)  return
c
c     update sublambdas, mapping derivatives and endpoint flags
c
      call mapsublmda (lambda)
c
c     ordinary absolute-topology energy routines consume installed
c     parameter arrays instead of applying elambda themselves; relative
c     dual topology routines install their own subsystem endpoint states
c
      if (.not. use_rel)  call altelec
      return
      end
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
c
c
c     the staged relative schedule maps the sublambdas its own way
c
      if (use_relstage) then
         call maprelstage (lmda)
         return
      end if
c
c     map the main lambda onto each sublambda a map drives
c
      if (use_plmdamap) then
         call sublmdamap (lmda,plmdamap,plmdaexp,plmdainvn,plmdainveps,
     &                    qntplmda0,qntplmda1,plambda,dpldlmda,
     &                    d2pldlmda2)
      end if
      if (use_elmdamap) then
         call sublmdamap (lmda,elmdamap,elmdaexp,elmdainvn,elmdainveps,
     &                    qntelmda0,qntelmda1,elambda,deldlmda,
     &                    d2eldlmda2)
      end if
      if (use_vlmdamap) then
         call sublmdamap (lmda,vlmdamap,vlmdaexp,vlmdainvn,vlmdainveps,
     &                    qntvlmda0,qntvlmda1,vlambda,dvldlmda,
     &                    d2vldlmda2)
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sublmdamap -- map one sublambda from lambda  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sublmdamap" maps the main lambda "lmda" onto the sublambda of a
c     single term, taking the exponential, inverse power or quintic
c     taper form named by "map", and returning that sublambda with its
c     first two derivatives with respect to the main lambda
c
c
      subroutine sublmdamap (lmda,map,nexp,invn,inveps,qnt0,qnt1,
     &                       sub,dsub,d2sub)
      implicit none
      integer nexp,invn
      real*8 lmda,inveps
      real*8 qnt0,qnt1
      real*8 sub,dsub,d2sub
      real*8 taper,dtaper,d2taper
      character*3 map
c
c
      if (map .eq. 'EXP') then
         call sublmdaexp (lmda,nexp,sub,dsub,d2sub)
      else if (map .eq. 'INV') then
         call sublmdainvpower (lmda,invn,inveps,sub,dsub,d2sub)
      else
         call quintaper (lmda,qnt0,qnt1,taper,dtaper,d2taper)
         sub = 1.0d0 - taper
         dsub = -dtaper
         d2sub = -d2taper
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
c     "maprelstage" maps the main lambda "lmda" onto the sublambdas of
c     the one staged relative leg with the following configuration:
c
c        LIG2   charge ligand 2 against the decoupled reference, its
c                 weight rising with the main lambda
c        VDWM   both ligands electrostatically decoupled while van der
c                 Waals morphs from ligand 2 onto ligand 1
c        LIG1   charge ligand 1 against the decoupled reference, its
c                 weight rising with the main lambda
c
c
      subroutine maprelstage (lmda)
      use dlmda
      use mutant
      implicit none
      real*8 lmda
c
c
c     van der Waals interpolates between the two coupled states on every
c     leg, morphing over its own map in the middle and held at one
c     end or the other while a ligand is being charged
c
      vrelst0 = rellig2
      vrelst1 = rellig1
c
c     the middle leg holds both ligands decoupled, so electrostatics
c     and polarization sit at the reference state and leave the chain
c     rule while van der Waals morphs across its map
c
      if (relstage .eq. 'VDWM') then
         erelst0 = relnone
         erelst1 = relnone
         elambda = 0.0d0
         deldlmda = 0.0d0
         d2eldlmda2 = 0.0d0
         call sublmdamap (lmda,vlmdamap,vlmdaexp,vlmdainvn,vlmdainveps,
     &                    qntvlmda0,qntvlmda1,vlambda,dvldlmda,
     &                    d2vldlmda2)
c
c     the ligand 1 leg charges ligand 1 against the decoupled reference
c     with van der Waals already morphed onto it
c
      else if (relstage .eq. 'LIG1') then
         erelst0 = relnone
         erelst1 = rellig1
         call sublmdamap (lmda,elmdamap,elmdaexp,elmdainvn,elmdainveps,
     &                    qntelmda0,qntelmda1,elambda,deldlmda,
     &                    d2eldlmda2)
         vlambda = 1.0d0
         dvldlmda = 0.0d0
         d2vldlmda2 = 0.0d0
c
c     the ligand 2 leg discharges ligand 2 as the main lambda rises, so
c     its weight is the complement of the map, with van der Waals still
c     on it
c
      else
         erelst0 = relnone
         erelst1 = rellig2
         call sublmdamap (lmda,elmdamap,elmdaexp,elmdainvn,elmdainveps,
     &                    qntelmda0,qntelmda1,elambda,deldlmda,
     &                    d2eldlmda2)
         elambda = 1.0d0 - elambda
         deldlmda = -deldlmda
         d2eldlmda2 = -d2eldlmda2
         vlambda = 0.0d0
         dvldlmda = 0.0d0
         d2vldlmda2 = 0.0d0
      end if
c
c     numerical guard on the map complement
c
      elambda = min(1.0d0,max(0.0d0,elambda))
c
c     polarization stages with the multipoles, same states same weight
c
      prelst0 = erelst0
      prelst1 = erelst1
      plambda = elambda
      dpldlmda = deldlmda
      d2pldlmda2 = d2eldlmda2
      return
      end
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine relpowerwt -- power law weight and derivs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "relpowerwt" evaluates the power law interpolation weight "x" to
c     the "nexp" and its first two derivatives, the linear case being
c     taken separately so that a zero sublambda never reaches a zero
c     power
c
c
      subroutine relpowerwt (x,nexp,w,dw,d2w)
      implicit none
      integer nexp
      real*8 x,w,dw,d2w
c
c
      w = x**nexp
      dw = 0.0d0
      d2w = 0.0d0
      if (nexp .eq. 1) then
         dw = 1.0d0
      else if (nexp .ge. 2) then
         dw = dble(nexp) * x**(nexp-1)
         d2w = dble(nexp) * dble(nexp-1) * x**(nexp-2)
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine relneed -- live dual topology endpoint test  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "relneed" decides which of the two dual topology endpoint states
c     a term has to build, given its interpolation weight "w", the first
c     and second derivatives "dw" and "d2w" of that weight with respect
c     to its sublambda, and the chain rule factors "chain" and "d2chain"
c     carrying the sublambda back to the main lambda
c
c
      subroutine relneed (w,dw,d2w,chain,d2chain,need0,need1)
      implicit none
      real*8 w,dw,d2w
      real*8 chain,d2chain
      real*8 c1,c2
      logical need0,need1
c
c
      c1 = dw * chain
      c2 = d2w*chain*chain + dw*d2chain
      need1 = (w .ne. 0.0d0) .or. (c1 .ne. 0.0d0) .or. (c2 .ne. 0.0d0)
      need0 = (w .ne. 1.0d0) .or. (c1 .ne. 0.0d0) .or. (c2 .ne. 0.0d0)
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
c     compute map
c
      lmda = x**exponent
      if (exponent .eq. 1) then
         dlmda = 1.0d0
         d2lmda = 0.0d0
      else if (exponent .eq. 2) then
         dlmda = 2.0d0 * x
         d2lmda = 2.0d0
      else
         expnt = dble(exponent)
         dlmda = expnt * x**(exponent-1)
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
c     compute map
c
      xval = x
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
      use mutant
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
c
c     sum up to get the total lambda derivative
c
      dedl = devdl + demdl + depdl
      d2edl2 = d2evdl2 + d2emdl2 + d2epdl2
      do i = 1, n
         do j = 1, 3
            dfsumdl(j,i) = dfvdl(j,i) + dfmdl(j,i) + dfpdl(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            dvirdl(j,i) = devvirdl(j,i) + demvirdl(j,i) + depvirdl(j,i)
         end do
      end do
c
c     zero d2edl2 and dfdl if ast polarization is used
c
      if (use_pdlmda .and. use_past) then
         d2evdl2 = 0.0d0
         d2emdl2 = 0.0d0
         d2epdl2 = 0.0d0
         d2edl2 = 0.0d0
         do i = 1, n
            do j = 1, 3
               dfvdl(j,i) = 0.0d0
               dfmdl(j,i) = 0.0d0
               dfpdl(j,i) = 0.0d0
               dfsumdl(j,i) = 0.0d0
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               devvirdl(j,i) = 0.0d0
               demvirdl(j,i) = 0.0d0
               depvirdl(j,i) = 0.0d0
               dvirdl(j,i) = 0.0d0
            end do
         end do
      end if
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

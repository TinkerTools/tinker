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
     &                    plmdaapmn,plmdaapmrho,qntplmda0,qntplmda1,
     &                    plambda,dpldlmda,d2pldlmda2)
      end if
      if (use_elmdamap) then
         call sublmdamap (lmda,elmdamap,elmdaexp,elmdainvn,elmdainveps,
     &                    elmdaapmn,elmdaapmrho,qntelmda0,qntelmda1,
     &                    elambda,deldlmda,d2eldlmda2)
      end if
      if (use_vlmdamap) then
         call sublmdamap (lmda,vlmdamap,vlmdaexp,vlmdainvn,vlmdainveps,
     &                    vlmdaapmn,vlmdaapmrho,qntvlmda0,qntvlmda1,
     &                    vlambda,dvldlmda,d2vldlmda2)
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
c     single term, taking the exponential, inverse power, asymmetric
c     power or quintic taper form named by "map", and returning that
c     sublambda with its first two derivatives with respect to the
c     main lambda
c
c
      subroutine sublmdamap (lmda,map,nexp,invn,inveps,apmn,apmrho,
     &                       qnt0,qnt1,sub,dsub,d2sub)
      implicit none
      integer nexp,invn
      integer apmn
      real*8 lmda,inveps
      real*8 apmrho
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
      else if (map .eq. 'APM') then
         call sublmdaapm (lmda,apmn,apmrho,sub,dsub,d2sub)
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
     &                    vlmdaapmn,vlmdaapmrho,qntvlmda0,qntvlmda1,
     &                    vlambda,dvldlmda,d2vldlmda2)
c
c     the ligand 1 leg charges ligand 1 against the decoupled reference
c     with van der Waals already morphed onto it
c
      else if (relstage .eq. 'LIG1') then
         erelst0 = relnone
         erelst1 = rellig1
         call sublmdamap (lmda,elmdamap,elmdaexp,elmdainvn,elmdainveps,
     &                    elmdaapmn,elmdaapmrho,qntelmda0,qntelmda1,
     &                    elambda,deldlmda,d2eldlmda2)
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
     &                    elmdaapmn,elmdaapmrho,qntelmda0,qntelmda1,
     &                    elambda,deldlmda,d2eldlmda2)
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
c     #############################################################
c     ##                                                         ##
c     ##  subroutine sublmdaapm -- asymmetric power law mapping  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "sublmdaapm" maps from main lambda to sublambda using a
c     normalized asymmetric power law and returns first and second
c     lambda derivatives; the map runs from zero to one over the unit
c     interval with a slope of "rho" at the decoupled end and a slope
c     of one at the coupled end, so that "rho" sets how much faster
c     the sublambda leaves the decoupled end than it arrives at the
c     coupled end, while "n" sets how quickly that head start decays
c
c
      subroutine sublmdaapm (x,n,rho,lmda,dlmda,d2lmda)
      implicit none
      integer n
      real*8 x
      real*8 rho
      real*8 lmda
      real*8 dlmda
      real*8 d2lmda
      real*8 xval
      real*8 omx
      real*8 rn
      real*8 left
      real*8 right
      real*8 denom
c
c
c     a flat slope ratio or a degenerate power gives the identity map,
c     the upper bound holding the normalization off of its pole
c
      xval = x
      if (n.lt.2 .or. rho.le.1.0d0 .or. rho.ge.dble(n)) then
         lmda = xval
         dlmda = 1.0d0
         d2lmda = 0.0d0
         return
      end if
c
c     compute the normalized asymmetric power map
c
      rn = dble(n)
      left = rn * (rho-1.0d0) / (rn-rho)
      right = left / rn
      denom = 1.0d0 + right
      omx = 1.0d0 - xval
      lmda = (xval + left*(1.0d0-omx**(n+1))/(rn+1.0d0)
     &           + right*xval**(n+1)/(rn+1.0d0)) / denom
      dlmda = (1.0d0 + left*omx**n + right*xval**n) / denom
      d2lmda = rn * (right*xval**(n-1)-left*omx**(n-1)) / denom
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
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine rdbiashead -- input lambda bias header  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "rdbiashead" reads the fixed-size header of a lambda bias history
c     file from an open unit into the scalar histogram and lambda
c     particle state, setting the flambda grid only for ost; the title
c     record is returned so the caller can tell which method wrote the
c     file, and the caller reads the rows
c
c
      subroutine rdbiashead (ihis,histfile,title)
      use bath
      use dlmda
      use iounit
      use mutant
      use ost
      implicit none
      integer ihis
      integer trimtext
      integer lmdastep0,lmdaintv0
      integer nlmda0,nflmda0,fli00
      integer nlmdahist0,sizelmdahist0
      real*8 wlmda0,wflmda0
      real*8 lambda0
      real*8 lmdatheta0,lmdavtheta0
      real*8 lmdamass0,lmdafric0,lmdadt0
      real*8 lmdadeltag0,oststdev0
      real*8 kelvin0
      character*240 record
      character*(*) histfile
      character*(*) title
c
c
c     read the title, scalar state and history label records
c
      read (ihis,10,err=90,end=90)  record
      title = record
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  lmdastep0,lmdaintv0,
     &   nlmda0,nflmda0,fli00,nlmdahist0,sizelmdahist0
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  wlmda0,wflmda0,oststdev0,
     &   kelvin0
      read (ihis,10,err=90,end=90)  record
      read (ihis,10,err=90,end=90)  record
      read (record,*,err=90,end=90)  lambda0,lmdatheta0,lmdavtheta0,
     &   lmdamass0,lmdafric0,lmdadt0,lmdadeltag0
      read (ihis,10,err=90,end=90)  record
   10 format (a240)
c
c     validate the stored dimensions
c
      if (nlmda0 .lt. 2)  goto 90
      if (lmdaintv0 .lt. 1)  goto 90
      if (kelvin0 .le. 0.0d0)  goto 90
      if (nlmdahist0 .lt. 0)  goto 90
      if (sizelmdahist0 .lt. nlmdahist0)  sizelmdahist0 = nlmdahist0
      if (sizelmdahist0 .lt. 1)  sizelmdahist0 = 1
      if (use_ost) then
         if (nflmda0 .lt. 1)  goto 90
         if (fli00 .lt. 1 .or. fli00 .gt. nflmda0)  goto 90
      end if
c
c     set scalar state from the history file
c
      lmdastep = lmdastep0
      lmdaintv = lmdaintv0
      call setlmdaphase
      nlmda = nlmda0
      nlmdahist = nlmdahist0
      sizelmdahist = sizelmdahist0
      wlmda = wlmda0
      wlmda2 = 0.5d0 * wlmda
      lambda = lambda0
      lmdaavg = 0.0d0
      lmdastd = 0.0d0
      dedlavg = 0.0d0
      dedlstd = 0.0d0
      lmdatheta = lmdatheta0
      lmdavtheta = lmdavtheta0
      lmdamass = lmdamass0
      lmdafric = lmdafric0
      lmdadt = lmdadt0
      lmdadeltag = lmdadeltag0
      kelvin = kelvin0
c
c     set the flambda grid and gaussian cutoff used only by ost
c
      if (use_ost) then
         nflmda = nflmda0
         fli0 = fli00
         wflmda = wflmda0
         wflmda2 = 0.5d0 * wflmda
         oststdev = oststdev0
      end if
      return
c
c     malformed history header
c
   90 continue
      close (unit=ihis)
      write (iout,20)  histfile(1:trimtext(histfile))
   20 format (/,' RDBIASHEAD  --  Error while Reading Lambda Bias',
     &           ' History Header',
     &        /,'                File Name :  ',a)
      call fatal
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine prtbiashead  --  output lambda bias header  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "prtbiashead" writes the fixed-size lambda bias history header
c     under the given title and history label from the current ost
c     histogram state
c
c
      subroutine prtbiashead (ihis,title,label)
      use bath
      use dlmda
      use mutant
      use ost
      implicit none
      integer ihis
      integer trimtext
      character*(*) title
      character*(*) label
c
c
c     write the title and the scalar histogram state
c
      write (ihis,10)  title(1:trimtext(title))
      write (ihis,20)
      write (ihis,30)  lmdastep,lmdaintv,nlmda,nflmda,
     &                 fli0,nlmdahist,sizelmdahist
      write (ihis,40)
      write (ihis,50)  wlmda,wflmda,oststdev,kelvin
      write (ihis,60)
      write (ihis,70)  lambda,lmdatheta,lmdavtheta,
     &                 lmdamass,lmdafric,lmdadt,lmdadeltag
      write (ihis,10)  label(1:trimtext(label))
   10 format (a)
   20 format (' Integer State :')
   30 format (7i12)
   40 format (' Grid State :')
   50 format (4d26.16)
   60 format (' Lambda State :')
   70 format (7d26.16)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine updbiashead  --  update lambda bias header  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "updbiashead" overwrites the fixed-size lambda bias history
c     header in place under the given title and history label;
c     unformatted stream output avoids truncating the appended history
c
c
      subroutine updbiashead (ihis,title,label)
      use bath
      use dlmda
      use mutant
      use ost
      implicit none
      integer ihis
      integer ieol
      integer leol
      character*240 record
      character*2 newline
      character*(*) title
      character*(*) label
c
c     preserve the file's existing line-ending convention
c
      read (ihis,pos=1)  record
      ieol = index(record,achar(10))
      newline = achar(10)//' '
      leol = 1
      if (ieol .gt. 1) then
         if (record(ieol-1:ieol-1) .eq. achar(13)) then
            newline = achar(13)//achar(10)
            leol = 2
         end if
      end if
c
c     format each header record internally and write its raw bytes
c
      write (record,10)  title
      write (ihis,pos=1)  record(1:len_trim(record)),newline(1:leol)
      write (record,20)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,30)  lmdastep,lmdaintv,nlmda,nflmda,
     &                   fli0,nlmdahist,sizelmdahist
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,40)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,50)  wlmda,wflmda,oststdev,kelvin
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,60)
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,70)  lambda,lmdatheta,lmdavtheta,
     &                   lmdamass,lmdafric,lmdadt,lmdadeltag
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
      write (record,10)  label
      write (ihis)  record(1:len_trim(record)),newline(1:leol)
   10 format (a)
   20 format (' Integer State :')
   30 format (7i12)
   40 format (' Grid State :')
   50 format (4d26.16)
   60 format (' Lambda State :')
   70 format (7d26.16)
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  function lmdabin -- get bin index for lambda  ##
c     ##                                                ##
c     ####################################################
c
c
c     "lmdabin" computes the lambda bin index for a lambda value
c
c
      function lmdabin (lambda)
      use dlmda
      implicit none
      integer lmdabin
      real*8 lambda
c
c
c     set lmdabin value
c
      lmdabin = nint(lambda / wlmda) + 1
      if (lmdabin .lt. 1)  lmdabin = 1
      if (lmdabin .gt. nlmda)  lmdabin = nlmda
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine setlmdaphase -- split the sample interval  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "setlmdaphase" divides the lambda sample interval into the
c     phase that propagates the lambda particle, the phase that
c     equilibrates at the frozen lambda and the phase that averages
c     dU/dlambda at that same fixed lambda; the phase counts are the
c     authoritative split, while lmdapcratio is only the nominal
c     fraction left over before truncation to whole samples
c
c
      subroutine setlmdaphase
      use dlmda
      implicit none
c
c
c     divide the interval, keeping at least one propagation step and
c     at least two samples to average
c
      if (lmdaintv .lt. 1)  lmdaintv = 1
      lmdapcratio = 1.0d0 - (lmdaparatio+lmdapbratio)
      lmdanpa = int(lmdaparatio*dble(lmdaintv))
      lmdanpb = int(lmdapbratio*dble(lmdaintv))
      lmdanpa = max(1,min(lmdanpa,lmdaintv-1))
      lmdanpb = max(0,min(lmdanpb,lmdaintv-lmdanpa))
      lmdanpc = lmdaintv - lmdanpa - lmdanpb
      do while (lmdanpc.lt.2 .and. lmdanpb.gt.0)
         lmdanpb = lmdanpb - 1
         lmdanpc = lmdanpc + 1
      end do
      do while (lmdanpc.lt.2 .and. lmdanpa.gt.1)
         lmdanpa = lmdanpa - 1
         lmdanpc = lmdanpc + 1
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine lmdalangevin -- propagate lambda particle  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "lmdalangevin" propagates the auxiliary lambda particle in
c     theta space, where lambda = sin(theta)**2
c
c
      subroutine lmdalangevin
      use bath
      use dlmda
      use math
      use mutant
      use units
      implicit none
      real*8 c
      real*8 force
      real*8 gamma
      real*8 normal
      real*8 sigma
      real*8 sinth
      real*8 ktm
      external normal
c
c
c     return if lambda dynamics parameters are invalid
c
      if (lmdadt .le. 0.0d0)  return
      if (lmdamass .le. 0.0d0)  return
c
c     force on theta from dU/dlambda and lambda = sin(theta)**2
c
      force = -deffdl * sin(2.0d0*lmdatheta)
c
c     propagate theta velocity with Langevin friction and noise
c
      gamma = max(0.0d0,lmdafric)
      if (gamma .gt. 0.0d0) then
         c = exp(-gamma*lmdadt)
         ktm = boltzmann * kelvin / lmdamass
         sigma = sqrt(ktm*(1.0d0-c*c))
         lmdavtheta = c*lmdavtheta
     &                 + (1.0d0-c)*force/(gamma*lmdamass)
     &                 + sigma*normal()
      else
         lmdavtheta = lmdavtheta + lmdadt*force/lmdamass
      end if
c
c     update theta and wrap it into the principal periodic interval
c
      lmdatheta = lmdatheta + lmdadt*lmdavtheta
      do while (lmdatheta .gt. pi)
         lmdatheta = lmdatheta - 2.0d0*pi
      end do
      do while (lmdatheta .le. -pi)
         lmdatheta = lmdatheta + 2.0d0*pi
      end do
c
c     map theta back to the main lambda
c
      sinth = sin(lmdatheta)
      lambda = sinth * sinth
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine efreelmda -- free energy at current lambda  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "efreelmda" computes the free energy at the current lambda by
c     integrating the mean force of the lambda bins with linear
c     interpolation, and also returns its lambda derivative
c
c
      subroutine efreelmda (eflmda,dfdl)
      use dlmda
      use mutant
      implicit none
      integer ilmda0,ilmda1
      real*8 eflmda,dfdl
      real*8 fl0,fl1
      real*8 lmda0,lmda1
      real*8 slope
      real*8 x
c
c
c     initialize free energy and derivative
c
      eflmda = 0.0d0
      dfdl = 0.0d0
c
c     handle endpoint at lambda = 0
c
      if (lambda .le. 0.0d0) then
         dfdl = lmdafmean(1)
         return
      end if
c
c     integrate over lambda intervals
c
      do ilmda0 = 1, nlmda-1
         ilmda1 = ilmda0 + 1
         lmda0 = dble(ilmda0-1) * wlmda
         lmda1 = dble(ilmda1-1) * wlmda
         fl0 = lmdafmean(ilmda0)
         fl1 = lmdafmean(ilmda1)
         slope = (fl1-fl0) / wlmda
c
c     integrate only to lambda if it lies in this interval
c
         if (lambda .le. lmda1) then
            x = lambda - lmda0
            eflmda = eflmda + fl0*x + 0.5d0*slope*x*x
            dfdl = fl0 + slope*x
            return
         end if
c
c     otherwise integrate the full interval
c
         eflmda = eflmda
     &            + 0.5d0*(fl0+fl1)*wlmda
      end do
c
c     handle endpoint at lambda = 1
c
      dfdl = lmdafmean(nlmda)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function efreetot -- total free energy from mean force  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "efreetot" computes the total free energy change by integrating
c     the mean force of the lambda bins over lambda using the
c     trapezoid rule
c
c
      function efreetot ()
      use dlmda
      implicit none
      integer ilmda
      real*8 efreetot
c
c
c     initialize free energy
c
      efreetot = 0.0d0
c
c     integrate over lambda bins
c
      do ilmda = 1, nlmda-1
         efreetot = efreetot + 0.5d0
     &                *(lmdafmean(ilmda)+lmdafmean(ilmda+1))*wlmda
      end do
      return
      end

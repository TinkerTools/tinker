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
c     sublambdas and installs the resulting electrostatic state
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
c     the ordinary energy routines consume the installed parameter
c     arrays instead of applying elambda themselves
c
      call altelec
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine grpscale  --  coupling scale per group  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "grpscale" returns the scale applied to the electrostatic and
c     polarization parameters of the environment and of each ligand
c     group at the coupling value "lmda", along with its derivative
c     with respect to "lmda"; the environment is unscaled, the charging
c     ligand carries "lmda", and in a staged relative leg the other
c     ligand is annihilated
c
c
      subroutine grpscale (lmda,sc,dsc)
      use dlmda
      use mutant
      implicit none
      real*8 lmda
      real*8 sc(0:2)
      real*8 dsc(0:2)
c
c
c     the first ligand group charges unless the leg charges the second
c
      sc(0) = 1.0d0
      dsc(0) = 0.0d0
      sc(1) = lmda
      dsc(1) = 1.0d0
      sc(2) = 0.0d0
      dsc(2) = 0.0d0
      if (use_rel .and. relstage.eq.'LIG2') then
         sc(1) = 0.0d0
         dsc(1) = 0.0d0
         sc(2) = lmda
         dsc(2) = 1.0d0
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine grpvlmda  --  van der Waals lambda of a group  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "grpvlmda" returns the van der Waals lambda "vlmd" seen by the
c     pairs of ligand group "ig" and the sign "vsgn" of its derivative
c     with respect to vlambda, where the second ligand of a relative
c     dual topology couples as the complement of vlambda
c
c
      subroutine grpvlmda (ig,vlmd,vsgn)
      use mutant
      implicit none
      integer ig
      real*8 vlmd,vsgn
c
c
c     the first ligand group couples as vlambda itself
c
      vlmd = vlambda
      vsgn = 1.0d0
      if (ig .eq. 2) then
         vlmd = 1.0d0 - vlambda
         vsgn = -1.0d0
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine grpsoft  --  soft core lambda terms by group  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "grpsoft" sets the soft core scale "lsc" and shift "shf" of each
c     ligand group for the form shared by the buffered 14-7 potential
c     and the Pauli repulsion, with scale lambda^n and shift equal to
c     alpha*lambda^m*(1-lambda)^p, where "pw" is n, "alpha" is alpha,
c     "pm" is m and "pp" is p
c
c
      subroutine grpsoft (pw,alpha,pm,pp,lsc,shf)
      implicit none
      integer ig
      real*8 pw,alpha,pm,pp
      real*8 vlmd,vsgn
      real*8 shf1,shf2
      real*8 lsc(2)
      real*8 shf(2)
c
c
c     set the soft core lambda terms for each ligand group
c
      do ig = 1, 2
         call grpvlmda (ig,vlmd,vsgn)
         lsc(ig) = vlmd**pw
         call lmdapow (vlmd,alpha,pm,pp,.false.,shf(ig),shf1,shf2)
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine grpsoft4  --  soft core terms and lambda derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "grpsoft4" sets the soft core terms of each ligand group as in
c     "grpsoft" along with the factors needed for their first and
c     second derivatives with respect to the group lambda, where
c     "lsc1" and "lsc2" hold the scale over lambda and lambda squared,
c     and "shf1" and "shf2" hold lambda and lambda squared times the
c     shift derivatives, so that no power of lambda goes negative;
c     "lsc1" also carries the sign of the derivative of the group
c     lambda with respect to vlambda, which is negative for a second
c     ligand group, and the second derivative factors are only set
c     when requested
c
c
      subroutine grpsoft4 (pw,alpha,pm,pp,lsc,shf,lsc1,lsc2,shf1,shf2)
      use dlmda
      implicit none
      integer ig
      real*8 pw,alpha,pm,pp
      real*8 vlmd,vsgn
      real*8 lsc(2),shf(2)
      real*8 lsc1(2),lsc2(2)
      real*8 shf1(2),shf2(2)
c
c
c     set the soft core terms and their lambda derivative factors
c
      do ig = 1, 2
         call grpvlmda (ig,vlmd,vsgn)
         lsc(ig) = vlmd**pw
         lsc1(ig) = vsgn * vlmd**(pw-1.0d0)
         lsc2(ig) = 0.0d0
         if (use_d2lmda)  lsc2(ig) = vlmd**(pw-2.0d0)
         call lmdapow (vlmd,alpha,pm,pp,use_d2lmda,shf(ig),shf1(ig),
     &                 shf2(ig))
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine grppcls  --  van der Waals pair class by group  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "grppcls" sets the class "mutcls" of the van der Waals type
c     pairs between the environment (group 0) and the ligand groups,
c     where -1 marks a pair that is skipped, 0 a pair at full weight,
c     and 1 or 2 a pair scaled by the soft core of that ligand group
c
c     a ligand to environment pair takes the soft core of its ligand,
c     the pairs within a ligand keep full weight under decoupling but
c     take its soft core under annihilation, and the two ligands of a
c     relative dual topology never interact
c
c
      subroutine grppcls
      use mutant
      implicit none
      integer g
c
c
c     set the class of each pair of atom groups
c
      mutcls(0,0) = 0
      do g = 1, 2
         mutcls(0,g) = g
         mutcls(g,0) = g
         mutcls(g,g) = 0
         if (vcouple .eq. 1)  mutcls(g,g) = g
      end do
      mutcls(1,2) = -1
      mutcls(2,1) = -1
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine grpvwt  --  van der Waals pair weights by group  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "grpvwt" finds the weights "wg" of the van der Waals type pairs
c     between the environment (group 0) and the ligand groups for the
c     soft core lambda power "pw", along with their first and second
c     derivatives "dwg" and "d2wg" with respect to vlambda
c
c     a pair carries full weight, the lambda power of a ligand group
c     or no weight as its class from "grppcls" sets, which gives the
c     pair loops their soft core group as well
c
c
      subroutine grpvwt (pw,wg,dwg,d2wg)
      use dlmda
      use mutant
      implicit none
      integer g,h,c
      real*8 pw,vlmd,vsgn
      real*8 v(2),dv(2),d2v(2)
      real*8 wg(0:2,0:2)
      real*8 dwg(0:2,0:2)
      real*8 d2wg(0:2,0:2)
c
c
c     get the lambda power of each ligand group and its derivatives
c
      do g = 1, 2
         call grpvlmda (g,vlmd,vsgn)
         call lmdapwr (vlmd,pw,use_d2lmda,v(g),dv(g),d2v(g))
         dv(g) = vsgn * dv(g)
      end do
c
c     weight each pair of groups by its pair class
c
      do g = 0, 2
         do h = 0, 2
            c = mutcls(h,g)
            wg(h,g) = 0.0d0
            dwg(h,g) = 0.0d0
            d2wg(h,g) = 0.0d0
            if (c .eq. 0) then
               wg(h,g) = 1.0d0
            else if (c .gt. 0) then
               wg(h,g) = v(c)
               dwg(h,g) = dv(c)
               d2wg(h,g) = d2v(c)
            end if
         end do
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function lmdadrv  --  lambda derivative of a term used  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "lmdadrv" tests whether a term whose lambda derivatives are
c     enabled by "use" must compute them, given the first and second
c     derivatives "dl" and "d2l" of its sublambda map; they are needed
c     while the map has slope, or curvature when second derivatives
c     are requested, since the chain rule zeroes them otherwise
c
c
      function lmdadrv (use,dl,d2l)
      use dlmda
      implicit none
      real*8 dl,d2l
      logical use,lmdadrv
c
c
      lmdadrv = use .and. (dl.ne.0.0d0
     &             .or. (use_d2lmda .and. d2l.ne.0.0d0))
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine vdlflags  --  van der Waals lambda deriv use  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "vdlflags" sets whether the terms that follow the van der Waals
c     lambda compute their first lambda derivatives, as "lmdadrv"
c     decides, and whether they also compute the second, force and
c     virial ones
c
c
      subroutine vdlflags (dodl,dod2)
      use dlmda
      implicit none
      logical dodl,dod2
      logical lmdadrv
c
c
c     the second derivatives are only needed alongside the first
c
      dodl = lmdadrv (use_vdlmda,dvldlmda,d2vldlmda2)
      dod2 = (dodl .and. use_d2lmda)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine dlcpair  --  pair lambda force and virial  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "dlcpair" adds the lambda derivative "dldx", "dldy", "dldz" of
c     the gradient of a central pair of atoms "i" and "k" separated
c     by "xr", "yr", "zr" to the force array "dfl" and its virial to
c     "dvir"
c
c
      subroutine dlcpair (i,k,xr,yr,zr,dldx,dldy,dldz,dfl,dvir)
      implicit none
      integer i,k
      real*8 xr,yr,zr
      real*8 dldx,dldy,dldz
      real*8 dfl(3,*)
      real*8 dvir(3,3)
c
c
      dfl(1,i) = dfl(1,i) + dldx
      dfl(2,i) = dfl(2,i) + dldy
      dfl(3,i) = dfl(3,i) + dldz
      dfl(1,k) = dfl(1,k) - dldx
      dfl(2,k) = dfl(2,k) - dldy
      dfl(3,k) = dfl(3,k) - dldz
      call dlcvir (xr,yr,zr,dldx,dldy,dldz,dvir)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine dlcvir  --  pair lambda derivative virial  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "dlcvir" adds to "dvir" the virial of the lambda derivative
c     "dldx", "dldy", "dldz" of the gradient of a central pair whose
c     separation is "xr", "yr", "zr"
c
c
      subroutine dlcvir (xr,yr,zr,dldx,dldy,dldz,dvir)
      implicit none
      real*8 xr,yr,zr
      real*8 dldx,dldy,dldz
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
      real*8 dvir(3,3)
c
c
      vxx = xr * dldx
      vyx = yr * dldx
      vzx = zr * dldx
      vyy = yr * dldy
      vzy = zr * dldy
      vzz = zr * dldz
      dvir(1,1) = dvir(1,1) + vxx
      dvir(2,1) = dvir(2,1) + vyx
      dvir(3,1) = dvir(3,1) + vzx
      dvir(1,2) = dvir(1,2) + vyx
      dvir(2,2) = dvir(2,2) + vyy
      dvir(3,2) = dvir(3,2) + vzy
      dvir(1,3) = dvir(1,3) + vzx
      dvir(2,3) = dvir(2,3) + vzy
      dvir(3,3) = dvir(3,3) + vzz
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ctlmdasc  --  charge transfer lambda scale  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ctlmdasc" returns the scale applied to charge transfer pairs
c     for each pair of atom groups, along with its first and second
c     derivatives with respect to the polarization lambda, and flags
c     the group pairs that carry energy or, when "dodl" is set, a
c     lambda derivative; single topology scales a pair by the product
c     of the two group scales from "grpscale", as for the multipoles,
c     and dual topology interpolates between the pair scales of the
c     two endpoint states with the power law weight of polarization,
c     so in a staged relative leg the other ligand is annihilated and
c     the two ligands never interact
c
c
      subroutine ctlmdasc (dodl,sc,dsc,d2sc,live)
      use dlmda
      use mutant
      implicit none
      integer i,k
      real*8 w,dw,d2w
      real*8 s0,s1
      real*8 g(0:2),dg(0:2)
      real*8 g0(0:2),dg0(0:2)
      real*8 g1(0:2),dg1(0:2)
      real*8 sc(0:2,0:2)
      real*8 dsc(0:2,0:2)
      real*8 d2sc(0:2,0:2)
      logical dodl
      logical live(0:2,0:2)
c
c
c     dual topology interpolates the pair scales of the endpoints
c
      if (use_epdt) then
         call lmdapwr (plambda,dble(epdtexp),.true.,w,dw,d2w)
         call grpscale (0.0d0,g0,dg0)
         call grpscale (1.0d0,g1,dg1)
         do i = 0, 2
            do k = 0, 2
               s0 = g0(k) * g0(i)
               s1 = g1(k) * g1(i)
               sc(k,i) = s0 + w*(s1-s0)
               dsc(k,i) = dw * (s1-s0)
               d2sc(k,i) = d2w * (s1-s0)
            end do
         end do
c
c     single topology takes the product of the two group scales
c
      else
         call grpscale (plambda,g,dg)
         do i = 0, 2
            do k = 0, 2
               sc(k,i) = g(k) * g(i)
               dsc(k,i) = dg(k)*g(i) + g(k)*dg(i)
               d2sc(k,i) = 2.0d0 * dg(k) * dg(i)
            end do
         end do
      end if
c
c     a group pair is live while it carries energy or a needed
c     lambda derivative
c
      do i = 0, 2
         do k = 0, 2
            live(k,i) = (sc(k,i) .ne. 0.0d0)
            if (dodl) then
               live(k,i) = (live(k,i) .or. dsc(k,i).ne.0.0d0
     &                         .or. d2sc(k,i).ne.0.0d0)
            end if
         end do
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine eleclive  --  live electrostatic groups  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "eleclive" flags the atom groups that carry multipoles,
c     polarizability or a lambda slope in the current evaluation;
c     the environment is always live, the charging ligand is live
c     unless it is held uncharged with no lambda derivative left to
c     take, and in a staged relative leg the other ligand never is;
c     a pair with a site outside the live groups contributes nothing,
c     so the pair loops skip it instead of forming zero times an
c     inverse distance for two sites that may overlap; "doskip" tells
c     the pair loops whether any atom lies outside the live groups
c
c
      subroutine eleclive (live,doskip)
      use atoms
      use dlmda
      use mutant
      implicit none
      integer i
      real*8 g(0:2),dg(0:2)
      logical live(0:2)
      logical doskip
      logical edrv,pdrv,dead
      logical lmdadrv
c
c
c     a program that bypasses the mutation setup has no ligand groups,
c     so every atom is put in the environment for the pair loops
c
      if (allocated(mutg)) then
         if (size(mutg) .lt. n)  deallocate (mutg)
      end if
      if (.not. allocated(mutg)) then
         allocate (mutg(n))
         do i = 1, n
            mutg(i) = 0
         end do
      end if
c
c     a term keeps its lambda derivative as "lmdadrv" decides, as in
c     the "gradient" routine
c
      edrv = lmdadrv (use_edlmda,deldlmda,d2eldlmda2)
      pdrv = lmdadrv (use_pdlmda,dpldlmda,d2pldlmda2)
      dead = (elambda.eq.0.0d0 .and. plambda.eq.0.0d0
     &           .and. .not.edrv .and. .not.pdrv)
c
c     the environment is always live, and a ligand group only if it
c     is the one that charges, as "grpscale" sets it
c
      call grpscale (1.0d0,g,dg)
      live(0) = .true.
      do i = 1, 2
         live(i) = (g(i).ne.0.0d0 .and. .not.dead)
      end do
c
c     the pair loops test group membership only when some atom
c     belongs to a group that is not live
c
      doskip = .false.
      do i = 1, n
         if (.not. live(mutg(i)))  doskip = .true.
      end do
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
      if (use_rel) then
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
         call lmdapwr (lmda,dble(nexp),.true.,sub,dsub,d2sub)
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
c     the staged relative schedule, where each declared leg spans its
c     own range of the main lambda with the following configuration:
c
c        LIG2   charge ligand 2 with ligand 1 annihilated, its
c                 weight falling as the main lambda rises
c        VDWM   both ligands electrostatically annihilated while van
c                 der Waals morphs from ligand 2 onto ligand 1
c        LIG1   charge ligand 1 with ligand 2 annihilated, its
c                 weight rising with the main lambda
c
c     the leg holding "lmda" is walked on its local coordinate, the
c     term it drives following its own map across the whole leg, and
c     polarization always carries the electrostatic weight
c
c
      subroutine maprelstage (lmda)
      use dlmda
      use mutant
      implicit none
      real*8 lmda
      real*8 s,ds
c
c
c     find the active leg and the local coordinate along it
c
      call relstgloc (lmda,s,ds)
c
c     the middle leg holds both ligands uncharged, so electrostatics
c     sits at the reference state and leaves the chain rule while
c     van der Waals morphs across its map
c
      if (relstage .eq. 'VDWM') then
         elambda = 0.0d0
         deldlmda = 0.0d0
         d2eldlmda2 = 0.0d0
         call sublmdamap (s,vlmdamap,vlmdaexp,vlmdainvn,vlmdainveps,
     &                    vlmdaapmn,vlmdaapmrho,0.0d0,1.0d0,
     &                    vlambda,dvldlmda,d2vldlmda2)
c
c     the ligand legs follow the electrostatic map across the leg
c
      else
         call sublmdamap (s,elmdamap,elmdaexp,elmdainvn,elmdainveps,
     &                    elmdaapmn,elmdaapmrho,0.0d0,1.0d0,
     &                    elambda,deldlmda,d2eldlmda2)
c
c     the ligand 1 leg charges ligand 1 with ligand 2 annihilated and
c     van der Waals already morphed onto ligand 1
c
         if (relstage .eq. 'LIG1') then
            vlambda = 1.0d0
            dvldlmda = 0.0d0
            d2vldlmda2 = 0.0d0
c
c     the ligand 2 leg discharges ligand 2 as the main lambda rises, so
c     its weight is the complement of the map, with van der Waals still
c     on it
c
         else
            elambda = 1.0d0 - elambda
            deldlmda = -deldlmda
            d2eldlmda2 = -d2eldlmda2
            vlambda = 0.0d0
            dvldlmda = 0.0d0
            d2vldlmda2 = 0.0d0
         end if
      end if
c
c     numerical guard on the map complement
c
      elambda = min(1.0d0,max(0.0d0,elambda))
c
c     carry the chain rule from the local coordinate of the leg back
c     to the main lambda
c
      deldlmda = deldlmda * ds
      d2eldlmda2 = d2eldlmda2 * ds * ds
      dvldlmda = dvldlmda * ds
      d2vldlmda2 = d2vldlmda2 * ds * ds
c
c     polarization stages with the multipoles, same state same weight
c
      plambda = elambda
      dpldlmda = deldlmda
      d2pldlmda2 = d2eldlmda2
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine relstgloc -- active leg and local coordinate  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "relstgloc" finds the declared leg of the staged relative
c     schedule that holds the main lambda "lmda" and makes it the
c     active leg, returning the local coordinate "s" running from zero
c     to one across that leg and its derivative "ds" with respect to
c     the main lambda; a main lambda outside every declared range holds
c     the end state of the adjacent leg, where the derivative is zero
c
c
      subroutine relstgloc (lmda,s,ds)
      use dlmda
      implicit none
      integer i,k
      real*8 lmda,s,ds
      real*8 lo,hi
      character*4 leg(3)
      data leg  / 'LIG2', 'VDWM', 'LIG1' /
c
c
c     take the last declared leg starting at or below the main lambda,
c     or the first declared leg when the main lambda lies below all
c
      k = 0
      do i = 1, 3
         if (use_relstg(i)) then
            if (k.eq.0 .or. lmda.ge.relstglo(i))  k = i
         end if
      end do
c
c     with no leg declared the active leg is left as it was set, and
c     spans the whole main lambda
c
      s = lmda
      ds = 1.0d0
      if (k .eq. 0)  return
c
c     set the active leg and the local coordinate along it
c
      relstage = leg(k)
      lo = relstglo(k)
      hi = relstghi(k)
      ds = 1.0d0 / (hi-lo)
      s = min(1.0d0,max(0.0d0,(lmda-lo)*ds))
      if (lmda .ge. hi)  s = 1.0d0
      if (lmda.lt.lo .or. lmda.gt.hi)  ds = 0.0d0
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine lmdapwr  --  lambda power and its derivs  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "lmdapwr" evaluates "x" to the power "a" along with its first
c     derivative with respect to "x", and its second derivative only
c     if "do2" is set; a derivative whose coefficient vanishes is set
c     to zero, so that a zero "x" never reaches a negative power in
c     the constant or linear case
c
c
      subroutine lmdapwr (x,a,do2,f,df,d2f)
      implicit none
      real*8 x,a
      real*8 f,df,d2f
      logical do2
c
c
      f = x**a
      df = 0.0d0
      d2f = 0.0d0
      if (a .ne. 0.0d0)  df = a * x**(a-1.0d0)
      if (do2 .and. a.ne.0.0d0 .and. a.ne.1.0d0) then
         d2f = a * (a-1.0d0) * x**(a-2.0d0)
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine lmdapow  --  soft core shift and its derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "lmdapow" evaluates the soft core shift "c" times "x" to the "a"
c     times one minus "x" to the "b", along with "x" times its first
c     derivative and, only if "do2" is set, "x" squared times its
c     second derivative; scaling by powers of "x" keeps both finite at
c     zero for any nonnegative "a", and a term whose coefficient
c     vanishes is skipped, so that neither factor reaches a negative
c     power at an endpoint
c
c
      subroutine lmdapow (x,c,a,b,do2,f,f1,f2)
      implicit none
      real*8 x,c,a,b
      real*8 f,f1,f2
      real*8 omx
      real*8 xa,xa1,xa2
      logical do2
c
c
c     find the shift along with its scaled lambda derivatives
c
      omx = 1.0d0 - x
      xa = c * x**a
      xa1 = xa * x
      xa2 = xa1 * x
      f = xa * omx**b
      f1 = 0.0d0
      if (a .ne. 0.0d0)  f1 = a * xa * omx**b
      if (b .ne. 0.0d0)  f1 = f1 - b*xa1*omx**(b-1.0d0)
      f2 = 0.0d0
      if (do2) then
         if (a .ne. 0.0d0) then
            f2 = a*(a-1.0d0) * xa * omx**b
            if (b .ne. 0.0d0) then
               f2 = f2 - 2.0d0*a*b * xa1 * omx**(b-1.0d0)
            end if
         end if
         if (b.ne.0.0d0 .and. b.ne.1.0d0) then
            f2 = f2 + b*(b-1.0d0) * xa2 * omx**(b-2.0d0)
         end if
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
c     carrying the sublambda back to the main lambda; the curvature
c     only matters when second derivatives are requested
c
c
      subroutine relneed (w,dw,d2w,chain,d2chain,need0,need1)
      implicit none
      real*8 w,dw,d2w
      real*8 chain,d2chain
      real*8 c1,c2
      logical need0,need1
      logical slope
      logical lmdadrv
c
c
      c1 = dw * chain
      c2 = d2w*chain*chain + dw*d2chain
      slope = lmdadrv (.true.,c1,c2)
      need1 = (w.ne.0.0d0 .or. slope)
      need0 = (w.ne.1.0d0 .or. slope)
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
      real*8 shift
      real*8 power
      real*8 root0
      real*8 denom
      real*8 base
c
c
c     compute map
c
      if (n .le. 1) then
         lmda = x
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
      base = x + shift
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
      if (n.lt.2 .or. rho.le.1.0d0 .or. rho.ge.dble(n)) then
         lmda = x
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
      omx = 1.0d0 - x
      lmda = (x + left*(1.0d0-omx**(n+1))/(rn+1.0d0)
     &           + right*x**(n+1)/(rn+1.0d0)) / denom
      dlmda = (1.0d0 + left*omx**n + right*x**n) / denom
      d2lmda = rn * (right*x**(n-1)-left*omx**(n-1)) / denom
c
c     the end states map exactly onto themselves, since rounding can
c     put the sublambda an ulp past one, where a soft core power of
c     one minus the sublambda is undefined
c
      if (x .eq. 0.0d0)  lmda = 0.0d0
      if (x .eq. 1.0d0)  lmda = 1.0d0
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
c     apply chain rule for second derivatives, force and virial wrt
c     global lambda, using the first derivatives before their rescale
c
      if (use_d2lmda) then
         d2epdl2 = d2epdl2 * dpldlmda*dpldlmda + depdl * d2pldlmda2
         d2evdl2 = d2evdl2 * dvldlmda*dvldlmda + devdl * d2vldlmda2
         d2emdl2 = d2emdl2 * deldlmda*deldlmda + demdl * d2eldlmda2
         d2edl2 = d2evdl2 + d2emdl2 + d2epdl2
         do i = 1, n
            do j = 1, 3
               dfpdl(j,i) = dfpdl(j,i) * dpldlmda
               dfmdl(j,i) = dfmdl(j,i) * deldlmda
               dfvdl(j,i) = dfvdl(j,i) * dvldlmda
               dfsumdl(j,i) = dfvdl(j,i) + dfmdl(j,i) + dfpdl(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               depvirdl(j,i) = depvirdl(j,i) * dpldlmda
               demvirdl(j,i) = demvirdl(j,i) * deldlmda
               devvirdl(j,i) = devvirdl(j,i) * dvldlmda
               dvirdl(j,i) = devvirdl(j,i) + demvirdl(j,i)
     &                          + depvirdl(j,i)
            end do
         end do
      end if
c
c     apply chain rule for derivative of energy wrt global lambda
c     and sum up to get the total lambda derivative
c
      depdl = depdl * dpldlmda
      devdl = devdl * dvldlmda
      demdl = demdl * deldlmda
      dedl = devdl + demdl + depdl
c
c     zero d2edl2 and dfdl unless second lambda derivs are requested
c
      if (.not. use_d2lmda) then
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
c     ##########################################################
c     ##                                                      ##
c     ##  function depcriteria -- interval sample acceptance  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "depcriteria" decides whether the samples collected over the
c     last interval are converged enough to be added to the lambda
c     bias, by comparing their deviation against a tolerance with
c     both absolute and relative parts; every interval is accepted
c     unless the convergence gate is turned on
c
c
      function depcriteria (avg,std)
      use dlmda
      implicit none
      logical depcriteria
      real*8 avg,std
      real*8 tolerance
c
c
c     accept every interval when the gate is off
c
      depcriteria = .true.
      if (.not. use_lmdacv)  return
c
c     accept only a deviation strictly inside the tolerance
c
      depcriteria = .false.
      tolerance = lmdacvstd + lmdacvrat*abs(avg)
      if (tolerance .gt. 0.0d0) then
         if (std/tolerance .lt. 1.0d0)  depcriteria = .true.
      end if
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
c     particle state; an abf header holds only the sample and lambda
c     particle state, a metadynamics header holds the gaussian counts,
c     the bin width, the gaussian cutoff and the temperature, while an
c     ost header also sets the flambda grid; the title record is
c     returned so the caller can tell which method wrote the file, and
c     the caller reads the rows
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
      logical abfhead
      logical metahdr
      character*240 record
      character*(*) histfile
      character*(*) title
c
c
c     read the title and the integer state label records; an abf
c     history keeps only the sample and lambda particle state
c
      read (ihis,10,err=90,end=90)  record
   10 format (a240)
      title = record
      abfhead = (index(title,abftitle(1:trimtext(abftitle))) .gt. 0)
      metahdr = (index(title,metatitle(1:trimtext(metatitle))) .gt. 0)
      read (ihis,10,err=90,end=90)  record
c
c     read the integer state, and the grid state of a metadynamics
c     or an ost history; a metadynamics history stores only the bin
c     width, so recover the bin count that goes with it
c
      read (ihis,10,err=90,end=90)  record
      if (abfhead) then
         read (record,*,err=90,end=90)  lmdastep0,lmdaintv0,
     &      nlmda0,nlmdahist0
         sizelmdahist0 = nlmdahist0
      else if (metahdr) then
         read (record,*,err=90,end=90)  lmdastep0,lmdaintv0,
     &      nlmdahist0,sizelmdahist0
         read (ihis,10,err=90,end=90)  record
         read (ihis,10,err=90,end=90)  record
         read (record,*,err=90,end=90)  wlmda0,oststdev0,kelvin0
         if (.not. (wlmda0 .gt. 0.0d0))  goto 90
         if (.not. (oststdev0 .gt. 0.0d0))  goto 90
         nlmda0 = max(2,nint(1.0d0/wlmda0)+1)
      else
         read (record,*,err=90,end=90)  lmdastep0,lmdaintv0,
     &      nlmda0,nflmda0,fli00,nlmdahist0,sizelmdahist0
         read (ihis,10,err=90,end=90)  record
         read (ihis,10,err=90,end=90)  record
         read (record,*,err=90,end=90)  wlmda0,wflmda0,oststdev0,
     &      kelvin0
      end if
c
c     read the lambda particle state and the history label
c
      read (ihis,10,err=90,end=90)  record
      if (index(record,'Lambda State') .eq. 0)  goto 90
      read (ihis,10,err=90,end=90)  record
      if (abfhead) then
         read (record,*,err=90,end=90)  lambda0,lmdatheta0,lmdavtheta0
      else
         read (record,*,err=90,end=90)  lambda0,lmdatheta0,lmdavtheta0,
     &      lmdamass0,lmdafric0,lmdadt0,lmdadeltag0
      end if
      read (ihis,10,err=90,end=90)  record
c
c     validate the stored dimensions
c
      if (nlmda0 .lt. 2)  goto 90
      if (lmdaintv0 .lt. 1)  goto 90
      if (nlmdahist0 .lt. 0)  goto 90
      if (sizelmdahist0 .lt. nlmdahist0)  sizelmdahist0 = nlmdahist0
      if (sizelmdahist0 .lt. 1)  sizelmdahist0 = 1
      if (.not. abfhead) then
         if (kelvin0 .le. 0.0d0)  goto 90
         if (use_ost .and. .not.metahdr) then
            if (nflmda0 .lt. 1)  goto 90
            if (fli00 .lt. 1 .or. fli00 .gt. nflmda0)  goto 90
         end if
      end if
c
c     an abf history derives its bin width from the bin count
c
      if (abfhead)  wlmda0 = 1.0d0 / dble(nlmda0-1)
c
c     set scalar state from the history file
c
      lmdastep = lmdastep0
      lmdaintv = lmdaintv0
      call setlmdaphase
      nlmda = nlmda0
      if (metahdr) then
         nmetahist = nlmdahist0
         sizemetahist = sizelmdahist0
      else
         nlmdahist = nlmdahist0
         sizelmdahist = sizelmdahist0
      end if
      wlmda = wlmda0
      wlmda2 = 0.5d0 * wlmda
      lambda = lambda0
      lmdaavg = 0.0d0
      lmdastd = 0.0d0
      dedlavg = 0.0d0
      dedlstd = 0.0d0
      lmdatheta = lmdatheta0
      lmdavtheta = lmdavtheta0
      if (abfhead)  return
c
c     set the lambda particle dynamics state kept by ost and by
c     metadynamics
c
      lmdamass = lmdamass0
      lmdafric = lmdafric0
      lmdadt = lmdadt0
      lmdadeltag = lmdadeltag0
      kelvin = kelvin0
      if (metahdr) then
         oststdev = oststdev0
         return
      end if
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
c     ##  subroutine fmtbiashead  --  format lambda bias header  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "fmtbiashead" formats the "nrec" records of the fixed-size
c     lambda bias history header under the given title and history
c     label; an abf title keeps only the sample and lambda particle
c     state, a metadynamics title keeps the gaussian counts, the bin
c     width, the gaussian cutoff and the temperature, otherwise the
c     full ost histogram state is formatted
c
c
      subroutine fmtbiashead (title,label,nrec,record)
      use bath
      use dlmda
      use mutant
      use ost
      implicit none
      integer nrec
      integer trimtext
      logical abfhead
      logical metahdr
      character*240 record(*)
      character*(*) title
      character*(*) label
c
c
c     format the title and the scalar histogram state
c
      abfhead = (index(title,abftitle(1:trimtext(abftitle))) .gt. 0)
      metahdr = (index(title,metatitle(1:trimtext(metatitle))) .gt. 0)
      nrec = 1
      write (record(nrec),10)  title
      nrec = nrec + 1
      write (record(nrec),20)
      nrec = nrec + 1
      if (abfhead) then
         write (record(nrec),80)  lmdastep,lmdaintv,nlmda,nlmdahist
      else if (metahdr) then
         write (record(nrec),80)  lmdastep,lmdaintv,nmetahist,
     &                            sizemetahist
         nrec = nrec + 1
         write (record(nrec),40)
         nrec = nrec + 1
         write (record(nrec),100)  wlmda,oststdev,kelvin
      else
         write (record(nrec),30)  lmdastep,lmdaintv,nlmda,nflmda,
     &                            fli0,nlmdahist,sizelmdahist
         nrec = nrec + 1
         write (record(nrec),40)
         nrec = nrec + 1
         write (record(nrec),50)  wlmda,wflmda,oststdev,kelvin
      end if
c
c     format the lambda particle state and the history label
c
      nrec = nrec + 1
      write (record(nrec),60)
      nrec = nrec + 1
      if (abfhead) then
         write (record(nrec),90)  lambda,lmdatheta,lmdavtheta
      else
         write (record(nrec),70)  lambda,lmdatheta,lmdavtheta,
     &                            lmdamass,lmdafric,lmdadt,lmdadeltag
      end if
      nrec = nrec + 1
      write (record(nrec),10)  label
   10 format (a)
   20 format (' Integer State :')
   30 format (7i12)
   40 format (' Grid State :')
   50 format (4d26.16)
   60 format (' Lambda State :')
   70 format (7d26.16)
   80 format (4i12)
   90 format (3d26.16)
  100 format (3d26.16)
      return
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
c     under the given title and history label, as "fmtbiashead"
c     formats it
c
c
      subroutine prtbiashead (ihis,title,label)
      implicit none
      integer i,ihis
      integer nrec
      integer trimtext
      character*240 record(8)
      character*(*) title
      character*(*) label
c
c
c     write the formatted header records
c
      call fmtbiashead (title,label,nrec,record)
      do i = 1, nrec
         write (ihis,10)  record(i)(1:trimtext(record(i)))
   10    format (a)
      end do
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
c     header in place under the given title and history label, as
c     "fmtbiashead" formats it; unformatted stream output avoids
c     truncating the appended history
c
c
      subroutine updbiashead (ihis,title,label)
      implicit none
      integer i,ihis
      integer ieol
      integer leol
      integer nrec
      integer nsize
      integer trimtext
      character*240 record(8)
      character*2 newline
      character*(*) title
      character*(*) label
c
c
c     preserve the file's existing line-ending convention, reading
c     no further than the end of a short header-only file
c
      inquire (unit=ihis,size=nsize)
      nsize = min(nsize,len(record(1)))
      record(1) = ' '
      if (nsize .gt. 0)  read (ihis,pos=1)  record(1)(1:nsize)
      ieol = index(record(1),achar(10))
      newline = achar(10)//' '
      leol = 1
      if (ieol .gt. 1) then
         if (record(1)(ieol-1:ieol-1) .eq. achar(13)) then
            newline = achar(13)//achar(10)
            leol = 2
         end if
      end if
c
c     write the raw bytes of each formatted header record
c
      call fmtbiashead (title,label,nrec,record)
      write (ihis,pos=1)  record(1)(1:trimtext(record(1))),
     &                    newline(1:leol)
      do i = 2, nrec
         write (ihis)  record(i)(1:trimtext(record(i))),newline(1:leol)
      end do
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
c     dU/dlambda at that same fixed lambda, the averaging phase taking
c     whatever the other two leave
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
c     #######################################################
c     ##                                                   ##
c     ##  subroutine elmdadyn -- adaptive lambda sampling  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "elmdadyn" advances the adaptive lambda bias by one dynamics
c     step for orthogonal space tempering, metadynamics or adaptive
c     biasing force; it builds the effective lambda derivative, saves
c     the interval samples, deposits the interval average at the end
c     of each interval and propagates the lambda particle; when the
c     lambda is moved by Monte Carlo trials the particle is left alone
c     and "lmdamonte" changes the lambda once the step is complete
c
c
      subroutine elmdadyn
      use dlmda
      use mutant
      use ost
      implicit none
      integer isamp
      integer nskip
      logical depcriteria
c
c
c     increment lmdastep step counter
c
      lmdastep = lmdastep + 1
c
c     build the effective lambda derivative from the unbiased value
c     summed by lmdachain and the bias saved this step by eostbias
c     or eabfbias
c
      if (use_ost) then
         ostdgdl = ostbdgdl + ostbdgdfl*d2edl2
         lmdaddgdl = lmdadfdl
         deffdl = dedl + ostdgdl - lmdaddgdl
      else if (use_meta) then
         deffdl = dedl + ostbdgdl
      else if (use_abf) then
         lmdaddgdl = lmdadfdl
         deffdl = dedl - lmdaddgdl
      end if
c
c     save all values in the interval, but average only after the
c     propagation and equilibration phases; only ost and abf keep
c     dU/dlambda, which they average at a fixed lambda
c
      isamp = mod(lmdastep-1,lmdaintv) + 1
      lmdallist(isamp) = lambda
      if (use_ost .or. use_abf)  lmdaflist(isamp) = dedl
c
c     deposit the interval average every lmdaintv steps; an ost or
c     abf interval is kept only when its samples are converged
c
      if (isamp .eq. lmdaintv) then
         nskip = lmdanpa + lmdanpb
         call avgstd (lmdallist,nskip+1,lmdanpc,lmdaavg,lmdastd)
         if (use_ost .or. use_abf) then
            call avgstd (lmdaflist,nskip+1,lmdanpc,dedlavg,dedlstd)
            if (depcriteria(dedlavg,dedlstd)) then
               if (use_ost)  call ostdeposit
               if (use_abf)  call abfdeposit
            end if
         else if (use_meta) then
            call metadeposit
         end if
      end if
c
c     propagate the lambda particle; ost and abf hold lambda fixed
c     after the propagation phase
c
      if (.not. use_lmdamc) then
         if (use_ost .or. use_abf) then
            if (isamp .le. lmdanpa)  call lmdalangevin
         else
            call lmdalangevin
         end if
      end if
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
c     theta space, where lambda follows the theta map "lmdathmap";
c     it draws the random kick and has "lmdalangstep" take the step
c
c
      subroutine lmdalangevin
      use dlmda
      implicit none
      real*8 gauss
      real*8 normal
      external normal
c
c
c     return if lambda dynamics parameters are invalid
c
      if (lmdadt .le. 0.0d0)  return
      if (lmdamass .le. 0.0d0)  return
c
c     draw the random kick only when the friction puts it to use
c
      gauss = 0.0d0
      if (lmdafric .gt. 0.0d0)  gauss = normal ()
      call lmdalangstep (gauss)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine lmdalangstep -- lambda particle time step  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "lmdalangstep" takes one Langevin step of the lambda particle
c     in theta space, with the standard normal deviate "gauss" as
c     its random kick; the force on theta is in kcal/mole, so the
c     mass is in kcal/mole*ps**2 and the size of the kick follows
c     from the thermal energy in kcal/mole
c
c
      subroutine lmdalangstep (gauss)
      use bath
      use dlmda
      use math
      use mutant
      use units
      implicit none
      real*8 c
      real*8 force
      real*8 gamma
      real*8 gauss
      real*8 sigma
      real*8 dldth
      real*8 ktm
c
c
c     force on theta from dU/dlambda and the theta map derivative
c
      call lmdathetamap (lmdatheta,lambda,dldth)
      force = -deffdl * dldth
c
c     propagate theta velocity with Langevin friction and noise
c
      gamma = max(0.0d0,lmdafric)
      if (gamma .gt. 0.0d0) then
         c = exp(-gamma*lmdadt)
         ktm = gasconst * kelvin / lmdamass
         sigma = sqrt(ktm*(1.0d0-c*c))
         lmdavtheta = c*lmdavtheta
     &                 + (1.0d0-c)*force/(gamma*lmdamass)
     &                 + sigma*gauss
      else
         lmdavtheta = lmdavtheta + lmdadt*force/lmdamass
      end if
c
c     update theta and wrap it into the periodic interval [0,pi)
c
      lmdatheta = lmdatheta + lmdadt*lmdavtheta
      lmdatheta = lmdatheta - pi*floor(lmdatheta/pi)
      if (lmdatheta .ge. pi)  lmdatheta = 0.0d0
      if (lmdatheta .lt. 0.0d0)  lmdatheta = 0.0d0
c
c     map theta back to the main lambda
c
      call lmdathetamap (lmdatheta,lambda,dldth)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine elmdabias  --  bias that follows lambda alone  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "elmdabias" finds the adaptive bias energy at the current
c     lambda for the methods whose bias depends on lambda alone,
c     adaptive biasing force and metadynamics, from the history as
c     it now stands; "slope" returns the value that each method saves
c     to drive the lambda particle, the mean force for abf and the
c     bias derivative for metadynamics; no saved value is changed
c
c
      subroutine elmdabias (vbias,slope)
      use dlmda
      use mutant
      implicit none
      real*8 vbias,slope
c
c
      if (use_abfdyn) then
         call efreelmda (lambda,vbias,slope)
         vbias = -vbias
      else
         call emetabias (lambda,vbias,slope)
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine lmdamonte  --  Monte Carlo move of lambda  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "lmdamonte" attempts a Monte Carlo change of the main lambda
c     at fixed coordinates, in place of propagating the lambda
c     particle; a trial is made at random, on average once every
c     "lmdamctrial" steps of the propagation phase of the sample
c     interval, by a step drawn uniformly up to "lmdamcmove" in
c     either direction and reflected at the ends of the lambda range
c
c     the trial is accepted by the Metropolis test on the change in
c     the potential energy plus the adaptive bias, both taken at the
c     coordinates, box and bias history in place when the move is
c     tried; "epot" is the energy of the completed step on entry, and
c     is replaced by the energy at the new lambda only when the move
c     is accepted
c
c     the induced dipoles, the lambda dependent parameters and the
c     forces and multipoles read by the saved output are returned
c     to what the step left after a rejected move, while
c     an accepted move keeps those of the trial; prediction of the
c     induced dipoles is assumed to be off, as is enforced for the
c     adaptive bias methods by "predict", so no history is changed
c
c
      subroutine lmdamonte (epot)
      use atoms
      use bath
      use bound
      use deriv
      use dlmda
      use mdstuf
      use mpole
      use mutant
      use ost
      use output
      use polar
      use potent
      use units
      use virial
      implicit none
      integer isamp
      real*8 epot
      real*8 eold,enew
      real*8 de,kt
      real*8 lold,ltrial
      real*8 vbias,slope
      real*8 energy,random
      real*8 dedlold,d2edl2old
      real*8 vbiasold,dfdlold
      real*8 bdgdlold,bdgdflold
      real*8 virsave(3,3)
      real*8, allocatable :: dtrial(:,:)
      real*8, allocatable :: desumold(:,:)
      real*8, allocatable :: rpoleold(:,:)
      real*8, allocatable :: polold(:,:,:)
      real*8, allocatable :: solvold(:,:,:)
      logical accept,stale
      logical savepol,savesolv
      logical savepole
      external random
c
c
c     ost and abf move lambda only in the propagation phase of the
c     sample interval, while metadynamics moves it at any step
c
      isamp = mod(lmdastep-1,lmdaintv) + 1
      if (use_ostdyn .or. use_abfdyn) then
         if (isamp .gt. lmdanpa)  return
      end if
c
c     decide whether to attempt a lambda change at this step
c
      if (random() .ge. 1.0d0/dble(lmdamctrial))  return
c
c     draw the trial lambda and reflect it into the lambda range
c
      kt = gasconst * kelvin
      ltrial = lambda + lmdamcmove*(2.0d0*random()-1.0d0)
      if (ltrial .lt. 0.0d0)  ltrial = -ltrial
      if (ltrial .gt. 1.0d0)  ltrial = 2.0d0 - ltrial
      nlmdamctry = nlmdamctry + 1
c
c     the energy of the step no longer holds once the coordinates or
c     the box have changed since its forces were found, as they have
c     after a scaling barostat or after a rigid body step, which finds
c     its forces before it moves; nor does the ost energy hold once
c     this step has added to the history
c
      stale = (integrate .eq. 'RIGIDBODY')
      if (isobaric .and. use_bounds) then
         if (barostat.eq.'BERENDSEN' .or. barostat.eq.'BUSSI')
     &      stale = .true.
      end if
      if (use_ostdyn .and. isamp.eq.lmdaintv)  stale = .true.
c
c     perform dynamic allocation of some local arrays
c
      savepol = use_polar
      savesolv = (use_polar .and. use_solv .and. allocated(uinds))
      savepole = (allocated(rpole) .and.
     &               (ustcsave .or. uchgsave .or. usyssave))
      if (use_ostdyn) then
         allocate (dtrial(3,n))
         allocate (desumold(3,n))
      end if
      if (savepole)  allocate (rpoleold(maxpole,n))
      if (savepol)  allocate (polold(3,n,4))
      if (savesolv)  allocate (solvold(3,n,4))
c
c     save the native lambda state before any energy is found, so
c     that a rejected move restores just what the step has left
c
      lold = lambda
      dedlold = dedl
      d2edl2old = d2edl2
      vbiasold = lmdavbias
      dfdlold = lmdadfdl
      bdgdlold = ostbdgdl
      bdgdflold = ostbdgdfl
      virsave = vir
      if (use_ostdyn)  desumold = desum
      if (savepole)  rpoleold = rpole
      if (savepol) then
         polold(:,:,1) = uind
         polold(:,:,2) = uinp
         polold(:,:,3) = udir
         polold(:,:,4) = udirp
      end if
      if (savesolv) then
         solvold(:,:,1) = uinds
         solvold(:,:,2) = uinps
         solvold(:,:,3) = udirs
         solvold(:,:,4) = udirps
      end if
c
c     get the old energy at the native lambda, using the energy of
c     the step unless it is stale; a bias that follows lambda alone
c     is taken from the history as it now stands, which may hold a
c     deposit made after the forces of the step were found
c
      if (use_ostdyn) then
         eold = epot
         if (stale)  call gradient (eold,dtrial)
      else
         if (stale) then
            eold = energy ()
         else
            eold = epot - lmdavbias
         end if
         call elmdabias (vbias,slope)
         eold = eold + vbias
      end if
c
c     get the new energy at the trial lambda; the ost bias follows
c     dU/dlambda, so it is found along with the lambda derivative
c
      lambda = ltrial
      if (use_ostdyn) then
         call gradient (enew,dtrial)
         vbias = lmdavbias
      else
         call elmdabias (vbias,slope)
         enew = energy () + vbias
      end if
c
c     acceptance from the change in potential energy plus bias
c
      de = enew - eold
      accept = (de .le. 0.0d0)
      if (.not.accept .and. kt.gt.0.0d0) then
         if (random() .lt. exp(-de/kt))  accept = .true.
      end if
c
c     accept the step, and keep the lambda particle at the new lambda
c
      if (accept) then
         nlmdamcacc = nlmdamcacc + 1
         epot = enew
         lmdavbias = vbias
         if (use_ostdyn) then
            ostdgdl = ostbdgdl + ostbdgdfl*d2edl2
            lmdaddgdl = lmdadfdl
         else if (use_abfdyn) then
            lmdadfdl = slope
            lmdaddgdl = slope
         else
            ostbdgdl = slope
         end if
         call lmdathetainv (lambda,lmdatheta)
         lmdavtheta = 0.0d0
c
c     reject the step, and restore values prior to the trial change
c
      else
         lambda = lold
         dedl = dedlold
         d2edl2 = d2edl2old
         lmdavbias = vbiasold
         lmdadfdl = dfdlold
         ostbdgdl = bdgdlold
         ostbdgdfl = bdgdflold
         if (savepol) then
            uind = polold(:,:,1)
            uinp = polold(:,:,2)
            udir = polold(:,:,3)
            udirp = polold(:,:,4)
         end if
         if (savesolv) then
            uinds = solvold(:,:,1)
            uinps = solvold(:,:,2)
            udirs = solvold(:,:,3)
            udirps = solvold(:,:,4)
         end if
c
c     reinstall the native lambda dependent parameters, along with
c     the charge flux monopoles and the global frame multipoles
c
         call refreshsublmda
         if (use_mpole .or. use_polar)  call alteprst
c
c     put back the forces of the step, and the global frame
c     multipoles as its last routine left them, which the saved
c     output reads; the multipoles are kept only if they are saved
c
         if (use_ostdyn)  desum = desumold
         if (savepole)  rpole = rpoleold
      end if
c
c     the virial of the completed step is left in place
c
      vir = virsave
c
c     perform deallocation of some local arrays
c
      if (use_ostdyn) then
         deallocate (dtrial)
         deallocate (desumold)
      end if
      if (savepole)  deallocate (rpoleold)
      if (savepol)  deallocate (polold)
      if (savesolv)  deallocate (solvold)
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  function lmdanbin  --  odd count of lambda bins  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "lmdanbin" returns the requested number of lambda bins kept to
c     an odd count of at least three, so bins are always centered at
c     lambda values of 0, 0.5 and 1
c
c
      function lmdanbin (n)
      implicit none
      integer n,lmdanbin
c
c
      lmdanbin = max(3,n)
      if (mod(lmdanbin,2) .eq. 0)  lmdanbin = lmdanbin + 1
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine lmdathkey  --  get theta map keywords  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "lmdathkey" sets the theta map of the lambda particle and its
c     sharpness from the keyfile, replacing an unknown map or a
c     sharpness outside (0,1) by the defaults
c
c
      subroutine lmdathkey
      use dlmda
      use keys
      implicit none
      integer i,next
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set the default smoothed triangle map and its sharpness
c
      lmdathmap = 'TRI'
      lmdathalpha = 0.999999999d0
c
c     search keywords for the theta map options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:17) .eq. 'LAMBDA-THETA-MAP ') then
            call getword (record,lmdathmap,next)
            call upcase (lmdathmap)
         else if (keyword(1:19) .eq. 'LAMBDA-THETA-ALPHA ') then
            read (string,*,err=10,end=10)  lmdathalpha
         end if
   10    continue
      end do
c
c     validate the theta map of the lambda particle
c
      if (lmdathmap.ne.'SIN' .and. lmdathmap.ne.'TRI') then
         lmdathmap = 'TRI'
      end if
      if (lmdathalpha.le.0.0d0 .or. lmdathalpha.ge.1.0d0) then
         lmdathalpha = 0.999999999d0
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine lmdathetamap  --  map theta onto main lambda  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "lmdathetamap" maps the lambda particle coordinate "theta" onto
c     the main lambda, taking the sine squared or smoothed triangle
c     form named by "lmdathmap", and returns the derivative of lambda
c     with respect to theta
c
c     the smoothed triangle map is lambda = 1/2 - asin(a*cos(2*theta))
c     / (2*asin(a)), which becomes sin(theta)**2 as the sharpness "a"
c     goes to zero and a triangle wave, flat in lambda, as "a" goes to
c     one; both maps have period pi in theta
c
c
      subroutine lmdathetamap (theta,lmda,dldth)
      use dlmda
      implicit none
      real*8 theta,lmda,dldth
      real*8 alpha,asina
      real*8 sinth,costh
c
c
      if (lmdathmap .eq. 'TRI') then
         alpha = lmdathalpha
         asina = asin(alpha)
         costh = cos(2.0d0*theta)
         lmda = 0.5d0 - asin(alpha*costh)/(2.0d0*asina)
         dldth = alpha * sin(2.0d0*theta)
     &              / (asina*sqrt(1.0d0-alpha*alpha*costh*costh))
      else
         sinth = sin(theta)
         lmda = sinth * sinth
         dldth = sin(2.0d0*theta)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine lmdathetainv  --  map main lambda onto theta  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "lmdathetainv" inverts the theta map named by "lmdathmap",
c     returning the theta in [0,pi/2] that gives the main lambda "lmda"
c
c
      subroutine lmdathetainv (lmda,theta)
      use dlmda
      implicit none
      real*8 lmda,theta
      real*8 lclip,alpha,asina
      real*8 arg
c
c
      lclip = min(1.0d0,max(0.0d0,lmda))
      if (lmdathmap .eq. 'TRI') then
         alpha = lmdathalpha
         asina = asin(alpha)
         arg = sin((1.0d0-2.0d0*lclip)*asina) / alpha
         arg = min(1.0d0,max(-1.0d0,arg))
         theta = 0.5d0 * acos(arg)
      else
         theta = asin(sqrt(lclip))
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function lmdathjac  --  smoothed theta map lambda density  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "lmdathjac" returns the factor dtheta/dlambda of the theta map,
c     by which the lambda particle visits each lambda, smoothed by a
c     normalized gaussian of width "sig" centered at lambda "lmda" and
c     cut off at the distance "cut", along with its two images
c     reflected at the lambda end points, as a metadynamics gaussian
c     deposited at "lmda" sees it
c
c     note dtheta/dlambda is infinite at the ends, but written as an
c     integral over theta the smoothed factor is finite, and it is
c     found by Simpson's rule over the theta range of each image
c
c
      function lmdathjac (lmda,sig,cut)
      use math
      implicit none
      integer i,img
      integer nstep
      parameter (nstep=1000)
      real*8 lmdathjac
      real*8 lmda,sig,cut
      real*8 lo,hi
      real*8 th0,th1,h
      real*8 theta,lth,dldth
      real*8 src,delta
      real*8 pref,wt,total
c
c
c     a zero width gaussian gives no smoothed factor
c
      lmdathjac = 0.0d0
      if (sig .le. 0.0d0)  return
      pref = 1.0d0 / (sig*sqrt(2.0d0*pi))
c
c     find the lambda window whose image lies within the cutoff
c
      do img = 1, 3
         if (img .eq. 1) then
            lo = lmda - cut
            hi = lmda + cut
         else if (img .eq. 2) then
            lo = -lmda - cut
            hi = -lmda + cut
         else
            lo = 2.0d0 - lmda - cut
            hi = 2.0d0 - lmda + cut
         end if
         lo = max(0.0d0,lo)
         hi = min(1.0d0,hi)
c
c     integrate the gaussian over the theta range of the window
c
         if (lo .lt. hi) then
            call lmdathetainv (lo,th0)
            call lmdathetainv (hi,th1)
            h = (th1-th0) / dble(nstep)
            total = 0.0d0
            do i = 0, nstep
               theta = th0 + dble(i)*h
               call lmdathetamap (theta,lth,dldth)
               if (img .eq. 1) then
                  src = lth
               else if (img .eq. 2) then
                  src = -lth
               else
                  src = 2.0d0 - lth
               end if
               delta = lmda - src
               wt = 2.0d0
               if (mod(i,2) .eq. 1)  wt = 4.0d0
               if (i.eq.0 .or. i.eq.nstep)  wt = 1.0d0
               total = total + wt*exp(-0.5d0*delta*delta/(sig*sig))
            end do
            lmdathjac = lmdathjac + pref*total*h/3.0d0
         end if
      end do
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
c     "efreelmda" computes the free energy at lambda "lmda" by
c     integrating the mean force of the lambda bins with linear
c     interpolation, and also returns its lambda derivative
c
c
      subroutine efreelmda (lmda,eflmda,dfdl)
      use dlmda
      implicit none
      integer ilmda0,ilmda1
      real*8 lmda
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
      if (lmda .le. 0.0d0) then
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
c     integrate only to lambda if it lies inside this interval
c
         if (lmda .lt. lmda1) then
            x = lmda - lmda0
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
c     the mean force of the lambda bins over the whole lambda range,
c     which "efreelmda" does by the trapezoid rule
c
c
      function efreetot ()
      implicit none
      real*8 efreetot
      real*8 dfdl
c
c
      call efreelmda (1.0d0,efreetot,dfdl)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine setlmdasched  --  build the lambda schedule  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "setlmdasched" fills the list of lambda values visited by the
c     sequential lambda windows, either from the explicit
c     "LAMBDA-WINDOW" keyword lines already parsed into "lmdawinlist"
c     or, by default, as the "nlmdawin" evenly spaced values that
c     descend from one to zero; the fraction of the run spent in
c     each window is resolved into "lmdawinfrac" summing to one
c
c
      subroutine setlmdasched (nwin,nwinset)
      use dlmda
      use iounit
      use mutant
      implicit none
      integer i
      integer nwin
      integer nunspec
      real*8 fsum,frem
      real*8, allocatable :: tlist(:)
      logical nwinset
      logical up,down
c
c
c     with no explicit windows, space them evenly over [0,1]
c
      if (nwin .eq. 0) then
         if (nlmdawin .lt. 2) then
            write (iout,10)
   10       format (/,' SETLMDASCHED  --  LAMBDA-NWINDOW must be',
     &                 ' at least 2')
            call fatal
         end if
         if (allocated(lmdawinlist))  deallocate (lmdawinlist)
         if (allocated(lmdawinfrac))  deallocate (lmdawinfrac)
         allocate (lmdawinlist(nlmdawin))
         allocate (lmdawinfrac(nlmdawin))
         do i = 1, nlmdawin
            lmdawinlist(i) = 1.0d0 - dble(i-1)/dble(nlmdawin-1)
            lmdawinfrac(i) = 1.0d0 / dble(nlmdawin)
         end do
         lmdawinlist(nlmdawin) = 0.0d0
      else
c
c     an explicit schedule sets the number of windows itself
c
         if (nwinset) then
            write (iout,20)
   20       format (/,' SETLMDASCHED  --  LAMBDA-NWINDOW and',
     &                 ' LAMBDA-WINDOW cannot both be used')
            call fatal
         end if
         nlmdawin = nwin
c
c     compact the parsed values down to the number of windows
c
         allocate (tlist(nlmdawin))
         tlist = lmdawinlist(1:nlmdawin)
         call move_alloc (tlist,lmdawinlist)
         allocate (tlist(nlmdawin))
         tlist = lmdawinfrac(1:nlmdawin)
         call move_alloc (tlist,lmdawinfrac)
c
c     a time share that was given must be a positive fraction
c
         fsum = 0.0d0
         nunspec = 0
         do i = 1, nlmdawin
            if (lmdawinfrac(i) .lt. 0.0d0) then
               nunspec = nunspec + 1
            else if (lmdawinfrac(i) .gt. 0.0d0) then
               fsum = fsum + lmdawinfrac(i)
            else
               write (iout,30)  i
   30          format (/,' SETLMDASCHED  --  LAMBDA-WINDOW',i5,
     &                    ' was given a time fraction that is',
     &                    ' not positive')
               call fatal
            end if
         end do
c
c     share whatever time is left among the windows that did not
c     ask for a fraction of their own
c
         if (nunspec .gt. 0) then
            frem = 1.0d0 - fsum
            if (frem .le. 0.0d0) then
               write (iout,40)  fsum
   40          format (/,' SETLMDASCHED  --  LAMBDA-WINDOW fractions',
     &                    ' total',f12.6,' leaving no time for the',
     &                 /,'                   windows without an',
     &                    ' explicit value')
               call fatal
            end if
            do i = 1, nlmdawin
               if (lmdawinfrac(i) .lt. 0.0d0) then
                  lmdawinfrac(i) = frem / dble(nunspec)
               end if
            end do
         end if
c
c     rescale the time shares so that they span the whole run
c
         fsum = 0.0d0
         do i = 1, nlmdawin
            fsum = fsum + lmdawinfrac(i)
         end do
         do i = 1, nlmdawin
            lmdawinfrac(i) = lmdawinfrac(i) / fsum
         end do
c
c     each window must sit within the physical lambda range
c
         do i = 1, nlmdawin
            if (.not. (lmdawinlist(i).ge.0.0d0 .and.
     &                 lmdawinlist(i).le.1.0d0)) then
               write (iout,50)  i,lmdawinlist(i)
   50          format (/,' SETLMDASCHED  --  LAMBDA-WINDOW',i5,
     &                    ' value',f12.6,' is outside [0,1]')
               call fatal
            end if
         end do
c
c     the schedule must walk in one direction without repeating
c
         if (nlmdawin .ge. 2) then
            up = .true.
            down = .true.
            do i = 2, nlmdawin
               if (lmdawinlist(i) .le. lmdawinlist(i-1))  up = .false.
               if (lmdawinlist(i) .ge. lmdawinlist(i-1))  down = .false.
            end do
            if (.not.up .and. .not.down) then
               write (iout,60)
   60          format (/,' SETLMDASCHED  --  LAMBDA-WINDOW values must',
     &                    ' increase or decrease monotonically')
               call fatal
            end if
         end if
      end if
c
c     start the schedule at its first window
c
      lmdawin = 1
      lambda = lmdawinlist(1)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine initlmdawin  --  set up lambda window layout  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "initlmdawin" divides a dynamics run of "nstep" steps among the
c     sequential lambda windows according to their requested time
c     fractions, and puts the main lambda at the start of the schedule
c
c
      subroutine initlmdawin (nstep)
      use dlmda
      use iounit
      use mutant
      implicit none
      integer i
      integer nstep
      integer istart
      real*8 cum
c
c
c     give each window the share of the run that it asked for, and
c     pin the last boundary so the whole trajectory is covered
c
      if (allocated(lmdawinend))  deallocate (lmdawinend)
      allocate (lmdawinend(nlmdawin))
      cum = 0.0d0
      do i = 1, nlmdawin
         cum = cum + lmdawinfrac(i)
         lmdawinend(i) = nint(cum*dble(nstep))
      end do
      lmdawinend(nlmdawin) = nstep
c
c     every window must receive at least one dynamics step
c
      istart = 0
      do i = 1, nlmdawin
         if (lmdawinend(i) .le. istart) then
            write (iout,10)  i,lmdawinlist(i)
   10       format (/,' INITLMDAWIN  --  LAMBDA-WINDOW',i5,
     &                 ' at lambda',f12.6,' was given no dynamics',
     &                 ' steps')
            call fatal
         end if
         istart = lmdawinend(i)
      end do
c
c     start the schedule at its first window
c
      lmdawin = 1
      lambda = lmdawinlist(1)
      call lmdawinsize (lmdawin,lmdawinlen,lmdawineq)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine lmdawinsize  --  size of one lambda window  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "lmdawinsize" finds the number of dynamics steps in one of
c     the lambda windows, and the number of those steps at the start
c     of the window that are discarded as equilibration
c
c
      subroutine lmdawinsize (iwin,nlen,neq)
      use dlmda
      implicit none
      integer iwin
      integer nlen,neq
c
c
c     measure the window against the preceding boundary
c
      nlen = lmdawinend(iwin)
      if (iwin .gt. 1)  nlen = lmdawinend(iwin) - lmdawinend(iwin-1)
      neq = int(dble(nlen) * lmdawinratio)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine lmdawinphase  --  find step in lambda window  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "lmdawinphase" finds the position of a dynamics step within the
c     production portion of the current lambda window; the value is
c     zero while the window is equilibrating, once the schedule has
c     run out, or when the window layout was never set up
c
c
      subroutine lmdawinphase (istep,iprod)
      use dlmda
      implicit none
      integer istep
      integer iprod
      integer istart
c
c
c     no step is in production outside of the window schedule
c
      iprod = 0
      if (.not. allocated(lmdawinend))  return
      if (lmdawin.lt.1 .or. lmdawin.gt.nlmdawin)  return
c
c     count the steps taken since the window finished equilibrating
c
      istart = 0
      if (lmdawin .gt. 1)  istart = lmdawinend(lmdawin-1)
      iprod = max(0,istep-istart-lmdawineq)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine lmdawindyn  --  end of step lambda work  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "lmdawindyn" advances the adaptive lambda bias and the lambda
c     particle by one dynamics step, takes the sample that a
c     sequential lambda method collects at each step, and attempts
c     the Monte Carlo move of an adaptively biased lambda; it is
c     called by the integrators once the step is complete and before
c     the trajectory is saved, so anything finished on a step is
c     written out on that step; "epot" is the potential energy of
c     the step, which an accepted lambda move replaces
c
c
      subroutine lmdawindyn (istep,epot)
      use dlmda
      implicit none
      integer istep
      real*8 epot
c
c
c     advance the adaptive lambda bias and the lambda particle, now
c     that the step and any Monte Carlo box change are complete
c
      if (use_ostdyn .or. use_metadyn .or. use_abfdyn)  call elmdadyn
c
c     collect dU/dlambda for thermodynamic integration
c
      if (use_ti)  call etidyn (istep)
c
c     find trial lambda energies for free energy perturbation
c
      if (use_fep)  call efepdyn (istep)
c
c     try a Monte Carlo move of an adaptively biased lambda
c
      if (use_lmdamc) then
         if (use_ostdyn .or. use_metadyn .or. use_abfdyn) then
            call lmdamonte (epot)
         end if
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine nextlmdawin  --  advance the lambda window  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "nextlmdawin" moves the main lambda to the value of the next
c     window in the schedule; past the final window the lambda is
c     left where it is, and further calls have no effect
c
c     note the accelerations saved by the integrator are not updated
c     for the new lambda, so the first half step of the velocities
c     in a new window still uses forces from the previous window
c
c
      subroutine nextlmdawin
      use dlmda
      use mutant
      implicit none
c
c
c     take the next lambda value
c
      if (lmdawin .gt. nlmdawin)  return
      lmdawin = lmdawin + 1
      if (lmdawin .le. nlmdawin) then
         lambda = lmdawinlist(lmdawin)
         call lmdawinsize (lmdawin,lmdawinlen,lmdawineq)
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine lmdawinstep  --  test for lambda window end  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "lmdawinstep" moves on to the next lambda window once the final
c     dynamics step of the current window has been taken
c
c
      subroutine lmdawinstep (istep)
      use dlmda
      implicit none
      integer istep
c
c
c     nothing is left to advance outside of the window schedule
c
      if (.not. allocated(lmdawinend))  return
      if (lmdawin.lt.1 .or. lmdawin.gt.nlmdawin)  return
c
c     move on to the next lambda window at the window boundary
c
      if (istep .eq. lmdawinend(lmdawin))  call nextlmdawin
      return
      end

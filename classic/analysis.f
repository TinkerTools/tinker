c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine analysis  --  energy components and analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "analysis" calls the series of routines needed to calculate
c     the potential energy and perform energy partitioning analysis
c     in terms of type of interaction or atom number
c
c
      subroutine analysis (energy)
      implicit none
      include 'sizes.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'potent.i'
      include 'vdwpot.i'
      integer i
      real*8 energy
      real*8 cutoff
c
c
c     zero out each of the potential energy components
c
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      er = 0.0d0
      es = 0.0d0
      elf = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     zero out energy partitioning components for each atom
c
      do i = 1, n
         aeb(i) = 0.0d0
         aea(i) = 0.0d0
         aeba(i) = 0.0d0
         aeub(i) = 0.0d0
         aeaa(i) = 0.0d0
         aeopb(i) = 0.0d0
         aeopd(i) = 0.0d0
         aeid(i) = 0.0d0
         aeit(i) = 0.0d0
         aet(i) = 0.0d0
         aept(i) = 0.0d0
         aebt(i) = 0.0d0
         aett(i) = 0.0d0
         aev(i) = 0.0d0
         aec(i) = 0.0d0
         aecd(i) = 0.0d0
         aed(i) = 0.0d0
         aem(i) = 0.0d0
         aep(i) = 0.0d0
         aer(i) = 0.0d0
         aes(i) = 0.0d0
         aelf(i) = 0.0d0
         aeg(i) = 0.0d0
         aex(i) = 0.0d0
      end do
c
c     zero out the total intermolecular energy
c
      einter = 0.0d0
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_group)  call bounds
c
c     remove any previous use of the replicates method
c
      cutoff = 0.0d0
      call replica (cutoff)
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call piscf
c
c     call the local geometry energy component routines
c
      if (use_bond)  call ebond3
      if (use_angle)  call eangle3
      if (use_strbnd)  call estrbnd3
      if (use_urey)  call eurey3
      if (use_angang)  call eangang3
      if (use_opbend)  call eopbend3
      if (use_opdist)  call eopdist3
      if (use_improp)  call eimprop3
      if (use_imptor)  call eimptor3
      if (use_tors)  call etors3
      if (use_pitors)  call epitors3
      if (use_strtor)  call estrtor3
      if (use_tortor)  call etortor3
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj3
         if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck3
         if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb3
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal3
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss3
      end if
c
c     call the electrostatic energy component routines
c
      if (use_charge)  call echarge3
      if (use_chgdpl)  call echgdpl3
      if (use_dipole)  call edipole3
      if (use_mpole .or. use_polar)  call empole3
      if (use_rxnfld)  call erxnfld3
c
c     call any miscellaneous energy component routines
c
      if (use_solv)  call esolv3
      if (use_metal)  call emetal3
      if (use_geom)  call egeom3
      if (use_extra)  call extra3
c
c     sum up to give the total potential energy
c
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
     &          + ep + er + es + elf + eg + ex
      energy = esum
c
c     sum up to give the total potential energy per atom
c
      do i = 1, n
         aesum(i) = aeb(i) + aea(i) + aeba(i) + aeub(i) + aeaa(i)
     &                 + aeopb(i) + aeopd(i) + aeid(i) + aeit(i)
     &                 + aet(i) + aept(i) + aebt(i) + aett(i) + aev(i)
     &                 + aec(i) + aecd(i) + aed(i) + aem(i) + aep(i)
     &                 + aer(i) + aes(i) + aelf(i) + aeg(i) + aex(i)
      end do
      return
      end

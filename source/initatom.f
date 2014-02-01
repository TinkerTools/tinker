c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2012  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine initatom  --  setup atoms in periodic table  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "initatom" sets the atomic symbol for each element in the
c     periodic table
c
c
      subroutine initatom
      implicit none
      include 'sizes.i'
      include 'ptable.i'
      integer i
      character*3 atmsym(maxele)
c
c     supported atomic symbols for the elements
c
      data atmsym  / 'H  ', 'He ', 'Li ', 'Be ', 'B  ', 'C  ', 'N  ',
     &               'O  ', 'F  ', 'Ne ', 'Na ', 'Mg ', 'Al ', 'Si ',
     &               'P  ', 'S  ', 'Cl ', 'Ar ', 'K  ', 'Ca ', 'Sc ',
     &               'Ti ', 'V  ', 'Cr ', 'Mn ', 'Fe ', 'Co ', 'Ni ',
     &               'Cu ', 'Zn ', 'Ga ', 'Ge ', 'As ', 'Se ', 'Br ',
     &               'Kr ', 'Rb ', 'Sr ', 'Y  ', 'Zr ', 'Nb ', 'Mo ',
     &               'Tc ', 'Ru ', 'Rh ', 'Pd ', 'Ag ', 'Cd ', 'In ',
     &               'Sn ', 'Sb ', 'Te ', 'I  ', 'Xe ', 'Cs ', 'Ba ',
     &               'La ', 'Ce ', 'Pr ', 'Nd ', 'Pm ', 'Sm ', 'Eu ',
     &               'Gd ', 'Tb ', 'Dy ', 'Ho ', 'Er ', 'Tm ', 'Yb ',
     &               'Lu ', 'Hf ', 'Ta ', 'W  ', 'Re ', 'Os ', 'Ir ',
     &               'Pt ', 'Au ', 'Hg ', 'Tl ', 'Pb ', 'Bi ', 'Po ',
     &               'At ', 'Rn ', 'Fr ', 'Ra ', 'Ac ', 'Th ', 'Pa ',
     &               'U  ', 'Np ', 'Pu ', 'Am ', 'Cm ', 'Bk ', 'Cf ',
     &               'Es ', 'Fm ', 'Md ', 'No ', 'Lr ', 'Rf ', 'Db ',
     &               'Sg ', 'Bh ', 'Hs ', 'Mt ', 'Ds ', 'Rg ', 'Cn ' /
c
c
c     set atomic symbols for each element in periodic table
c
      do i = 1, maxele
         elemnt(i) = atmsym(i)
      end do
      return
      end

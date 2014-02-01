c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2007  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program nuclion  --  insert ions to neutralize phosphates  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "nuclion" neutralizes the phosphate backbone of nucleic acids
c     by adding a cation adjacent to each phosphate group
c
c
      program nuclion
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'katoms.i'
      integer i,m,ixyz
      integer ctype
      integer freeunit
      character*3 cname
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     find the atom type and name for a monovalent cation;
c     current version uses a sodium ion as the cation
c
      ctype = 0
      cname = '   '
      do i = 1, maxtyp
         if (atmnum(i) .eq. 11) then
            ctype = i
            cname = symbol(i)
         end if
      end do
c
c     place a cation at each phosphate and append to atom list;
c     assumes monovalent oxygens directly follow phosphorus
c
      m = n
      do i = 1, m
         if (atomic(i) .eq. 15) then
            n = n + 1
            type(n) = ctype
            name(n) = cname
            n12(n) = 0
            x(n) = x(i) + 2.0d0*(x(i+1)+x(i+2)-2.0d0*x(i))
            y(n) = y(i) + 2.0d0*(y(i+1)+y(i+2)-2.0d0*y(i))
            z(n) = z(i) + 2.0d0*(z(i+1)+z(i+2)-2.0d0*z(i))
         end if
      end do
c
c     write out a coordinates file containing the ions
c
      ixyz = freeunit ()
      call prtxyz (ixyz)
      call final
      end

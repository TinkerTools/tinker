c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine zatom  --  adds a single atom to Z-matrix  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "zatom" adds an atom to the end of the current Z-matrix
c     and then increments the atom counter; atom type, defining
c     atoms and internal coordinates are passed as arguments
c
c
      subroutine zatom (bionum,bond,angle,dihed,iz1,iz2,iz3,iz4)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'fields.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'zclose.i'
      include 'zcoord.i'
      integer bionum
      integer iz1,iz2
      integer iz3,iz4
      real*8 bond,angle
      real*8 dihed
c
c
c     fill various arrays with information for this atom
c
      if (bionum .gt. 0) then
         type(n) = biotyp(bionum)
         if (type(n) .ne. 0) then
            name(n) = symbol(type(n))
         else
            name(n) = '   '
         end if
         zbond(n) = bond
         zang(n) = angle
         ztors(n) = dihed
         if (ztors(n) .lt. -180.0d0) then
            ztors(n) = ztors(n) + 360.0d0
         else if (ztors(n) .gt. 180.0d0) then
            ztors(n) = ztors(n) - 360.0d0
         end if
         iz(1,n) = iz1
         iz(2,n) = iz2
         iz(3,n) = iz3
         iz(4,n) = iz4
c
c     increment atom counter and check for too many atoms
c
         n = n + 1
         if (n .gt. maxatm) then
            write (iout,10)  maxatm
   10       format (/,' ZATOM  --  The Maximum of',i8,' Atoms',
     &                 ' has been Exceeded')
            call fatal
         end if
c
c     add an extra bond to make a ring closure
c
      else if (bionum .eq. -1) then
         nadd = nadd + 1
         iadd(1,nadd) = iz1
         iadd(2,nadd) = iz2
c
c     delete an extra bond to make separate molecules
c
      else if (bionum .eq. -2) then
         ndel = ndel + 1
         idel(1,ndel) = iz1
         idel(2,ndel) = iz2
      end if
      return
      end

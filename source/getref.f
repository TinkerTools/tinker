c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getref  --  get structure from reference area  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getref" copies structure information from the reference area
c     into the standard variables for the current system structure
c
c
      subroutine getref (iref)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'refer.i'
      include 'titles.i'
      integer i,j,iref
c
c
c     retrieve the filename and title line for the structure
c
      filename = reffile(iref)
      leng = refleng(iref)
      title = reftitle(iref)
      ltitle = refltitle(iref)
c
c     retrieve the coordinates, type and connectivity of each atom
c
      n = nref(iref)
      do i = 1, n
         name(i) = refnam(i,iref)
         x(i) = xref(i,iref)
         y(i) = yref(i,iref)
         z(i) = zref(i,iref)
         type(i) = reftyp(i,iref)
         n12(i) = n12ref(i,iref)
         do j = 1, n12(i)
            i12(j,i) = i12ref(j,i,iref)
         end do
      end do
      return
      end

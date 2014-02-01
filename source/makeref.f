c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1996  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine makeref  --  copy structure to reference area  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "makeref" copies the information contained in the "xyz" file
c     of the current structure into corresponding reference areas
c
c
      subroutine makeref (iref)
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
c     copy the filename and title line for the structure
c
      reffile(iref) = filename
      refleng(iref) = leng
      reftitle(iref) = title
      refltitle(iref) = ltitle
c
c     copy the coordinates, type and connectivity of each atom
c
      nref(iref) = n
      do i = 1, n
         refnam(i,iref) = name(i)
         xref(i,iref) = x(i)
         yref(i,iref) = y(i)
         zref(i,iref) = z(i)
         reftyp(i,iref) = type(i)
         n12ref(i,iref) = n12(i)
         do j = 1, n12(i)
            i12ref(j,i,iref) = i12(j,i)
         end do
      end do
      return
      end

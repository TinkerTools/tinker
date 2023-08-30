c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2022 by Moses Chung, Zhi Wang & Jay Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kgaussvol  --  GaussVol parameter assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kgaussvol" allocates global arrays and assigns atomic radii and
c     volume
c
c
      subroutine kgaussvol
      use atomid
      use atoms
      use gssvol
      use kvdws
      use math
      implicit none
      integer i
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(gvradius))  deallocate (gvradius)
      if (allocated(gvradius2))  deallocate (gvradius2)
      if (allocated(gvvol))  deallocate (gvvol)
      if (allocated(gvvol2))  deallocate (gvvol2)
      if (allocated(gvdr))  deallocate (gvdr)
      if (allocated(gvdv))  deallocate (gvdv)
      if (allocated(gvfree_volume))  deallocate (gvfree_volume)
      if (allocated(gvself_volume))  deallocate (gvself_volume)
      allocate (gvradius(n))
      allocate (gvradius2(n))
      allocate (gvvol(n))
      allocate (gvvol2(n))
      allocate (gvdr(3,n))
      allocate (gvdv(n))
      allocate (gvfree_volume(n))
      allocate (gvself_volume(n))
c
c     set atomic radii
c
      do i = 1, n
         gvradius(i) = rad(class(i))
         gvradius2(i) = rad(class(i)) + rad_offset
      end do
c
c     set atomic volume
c
      do i = 1, n
         if (use(i) .and. atomic(i) .ne. 1) then
            gvvol(i) = 4.0d0/3.0d0*PI*gvradius(i)**3
            gvvol2(i) = 4.0d0/3.0d0*PI*gvradius2(i)**3
         else
            gvvol(i) = 0.0d0
            gvvol2(i) = 0.0d0
         end if
      end do
      return
      end

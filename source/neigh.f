c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module neigh  --  pairwise neighbor list indices & storage  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     maxvlst     maximum size of van der Waals pair neighbor lists
c     maxelst     maximum size of electrostatic pair neighbor lists
c     maxulst     maximum size of dipole preconditioner pair lists
c     nvlst       number of sites in list for each vdw site
c     vlst        site numbers in neighbor list of each vdw site
c     nelst       number of sites in list for each electrostatic site
c     elst        site numbers in list of each electrostatic site
c     nulst       number of sites in list for each preconditioner site
c     ulst        site numbers in list of each preconditioner site
c     lbuffer     width of the neighbor list buffer region
c     pbuffer     width of the preconditioner list buffer region
c     lbuf2       square of half the neighbor list buffer width
c     pbuf2       square of half the preconditioner list buffer width
c     vbuf2       square of van der Waals cutoff plus the list buffer
c     vbufx       square of van der Waals cutoff plus 2X list buffer
c     dbuf2       square of dispersion cutoff plus the list buffer
c     dbufx       square of dispersion cutoff plus 2X list buffer
c     cbuf2       square of charge cutoff plus the list buffer
c     cbufx       square of charge cutoff plus 2X list buffer
c     mbuf2       square of multipole cutoff plus the list buffer
c     mbufx       square of multipole cutoff plus 2X list buffer
c     ubuf2       square of preconditioner cutoff plus the list buffer
c     ubufx       square of preconditioner cutoff plus 2X list buffer
c     xvold       x-coordinate at last vdw/dispersion list update
c     yvold       y-coordinate at last vdw/dispersion list update
c     zvold       z-coordinate at last vdw/dispersion list update
c     xeold       x-coordinate at last electrostatic list update
c     yeold       y-coordinate at last electrostatic list update
c     zeold       z-coordinate at last electrostatic list update
c     xuold       x-coordinate at last preconditioner list update
c     yuold       y-coordinate at last preconditioner list update
c     zuold       z-coordinate at last preconditioner list update
c     dovlst      logical flag to rebuild vdw neighbor list
c     dodlst      logical flag to rebuild dispersion neighbor list
c     doclst      logical flag to rebuild charge neighbor list
c     domlst      logical flag to rebuild multipole neighbor list
c     doulst      logical flag to rebuild preconditioner neighbor list
c
c
      module neigh
      implicit none
      integer maxvlst
      integer maxelst
      integer maxulst
      integer, allocatable :: nvlst(:)
      integer, allocatable :: vlst(:,:)
      integer, allocatable :: nelst(:)
      integer, allocatable :: elst(:,:)
      integer, allocatable :: nulst(:)
      integer, allocatable :: ulst(:,:)
      real*8 lbuffer,pbuffer
      real*8 lbuf2,pbuf2
      real*8 vbuf2,vbufx
      real*8 dbuf2,dbufx
      real*8 cbuf2,cbufx
      real*8 mbuf2,mbufx
      real*8 ubuf2,ubufx
      real*8, allocatable :: xvold(:)
      real*8, allocatable :: yvold(:)
      real*8, allocatable :: zvold(:)
      real*8, allocatable :: xeold(:)
      real*8, allocatable :: yeold(:)
      real*8, allocatable :: zeold(:)
      real*8, allocatable :: xuold(:)
      real*8, allocatable :: yuold(:)
      real*8, allocatable :: zuold(:)
      logical dovlst,dodlst
      logical doclst,domlst
      logical doulst
      save
      end

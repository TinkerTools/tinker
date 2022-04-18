c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module kangs  --  bond angle bend forcefield parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxna    maximum number of harmonic angle bend parameter entries
c     maxna5   maximum number of 5-membered ring angle bend entries
c     maxna4   maximum number of 4-membered ring angle bend entries
c     maxna3   maximum number of 3-membered ring angle bend entries
c     maxnap   maximum number of in-plane angle bend parameter entries
c     maxnaf   maximum number of Fourier angle bend parameter entries
c
c     acon     force constant parameters for harmonic angle bends
c     acon5    force constant parameters for 5-ring angle bends
c     acon4    force constant parameters for 4-ring angle bends
c     acon3    force constant parameters for 3-ring angle bends
c     aconp    force constant parameters for in-plane angle bends
c     aconf    force constant parameters for Fourier angle bends
c     ang      bond angle parameters for harmonic angle bends
c     ang5     bond angle parameters for 5-ring angle bends
c     ang4     bond angle parameters for 4-ring angle bends
c     ang3     bond angle parameters for 3-ring angle bends
c     angp     bond angle parameters for in-plane angle bends
c     angf     phase shift angle and periodicity for Fourier bends
c     ka       string of atom classes for harmonic angle bends
c     ka5      string of atom classes for 5-ring angle bends
c     ka4      string of atom classes for 4-ring angle bends
c     ka3      string of atom classes for 3-ring angle bends
c     kap      string of atom classes for in-plane angle bends
c     kaf      string of atom classes for Fourier angle bends
c
c
      module kangs
      implicit none
      integer maxna
      integer maxna5
      integer maxna4
      integer maxna3
      integer maxnap
      integer maxnaf
      real*8, allocatable :: acon(:)
      real*8, allocatable :: acon5(:)
      real*8, allocatable :: acon4(:)
      real*8, allocatable :: acon3(:)
      real*8, allocatable :: aconp(:)
      real*8, allocatable :: aconf(:)
      real*8, allocatable :: ang(:,:)
      real*8, allocatable :: ang5(:,:)
      real*8, allocatable :: ang4(:,:)
      real*8, allocatable :: ang3(:,:)
      real*8, allocatable :: angp(:,:)
      real*8, allocatable :: angf(:,:)
      character*12, allocatable :: ka(:)
      character*12, allocatable :: ka5(:)
      character*12, allocatable :: ka4(:)
      character*12, allocatable :: ka3(:)
      character*12, allocatable :: kap(:)
      character*12, allocatable :: kaf(:)
      save
      end

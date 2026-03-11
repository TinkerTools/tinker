c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module output  --  output file format control parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     nonly       total number of save sites in the system
c     ionly       number of the atom for each save site
c     ionlyinv    inverse map of ionly
c     print3n     array for printing atomic 3D vectors
c     archive     logical flag for coordinates in Tinker XYZ format
c     binary      logical flag for coordinates in DCD binary format
c     noversion   logical flag governing use of filename versions
c     overwrite   logical flag to overwrite intermediate files inplace
c     coordsave   logical flag to save coordinates
c     dynsave     logical flag to save dynamics (.dyn) file
c     onlysave    logical flag to only save certain coordinates
c     arcsave     logical flag to save coordinates in Tinker XYZ format
c     dcdsave     logical flag to save coordinates in DCD binary format
c     cyclesave   logical flag to mark use of numbered cycle files
c     velsave     logical flag to save velocity vector components
c     frcsave     logical flag to save force vector components
c     uindsave    logical flag to save induced atomic dipoles
c     ustcsave    logical flag to save static atomic dipoles
c     uchgsave    logical flag to save charge atomic dipoles
c     usyssave    logical flag to save unique atom type dipole moment
c     vsyssave    logical flag to save unique atom type velocity
c     udirsave    logical flag to save the direct induced atomic dipoles
c     defsave     logical flag to save the direct electric field
c     tefsave     logical flag to save the total electric field
c     coordtype   selects Cartesian, internal, rigid body or none
c
c
      module output
      implicit none
      integer nonly
      integer, allocatable :: ionly(:)
      integer, allocatable :: ionlyinv(:)
      real*8, allocatable :: print3n(:,:)
      logical archive
      logical binary
      logical noversion
      logical overwrite
      logical coordsave
      logical dynsave
      logical cyclesave
      logical onlysave
      logical arcsave
      logical dcdsave
      logical velsave
      logical frcsave
      logical uindsave
      logical ustcsave
      logical uchgsave
      logical usyssave
      logical vsyssave
      logical udirsave
      logical defsave
      logical tefsave
      character*9 coordtype
      save
      end

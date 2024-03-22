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
c     archive     logical flag for coordinates in Tinker XYZ format
c     binary      logical flag for coordinates in DCD binary format
c     noversion   logical flag governing use of filename versions
c     overwrite   logical flag to overwrite intermediate files inplace
c     coordsave   logical flag to save coordinates
c     arcsave     logical flag to save coordinates in Tinker XYZ format
c     dcdsave     logical flag to save coordinates in DCD binary format
c     cyclesave   logical flag to mark use of numbered cycle files
c     velsave     logical flag to save velocity vector components
c     frcsave     logical flag to save force vector components
c     uindsave    logical flag to save induced atomic dipoles
c     ustcsave    logical flag to save static atomic dipoles
c     usyssave    logical flag to save total dipole of the system
c     coordtype   selects Cartesian, internal, rigid body or none
c
c
      module output
      implicit none
      logical archive
      logical binary
      logical noversion
      logical overwrite
      logical coordsave
      logical cyclesave
      logical arcsave
      logical dcdsave
      logical velsave
      logical frcsave
      logical uindsave
      logical ustcsave
      logical usyssave
      character*9 coordtype
      save
      end

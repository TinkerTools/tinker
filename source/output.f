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
c     archive     logical flag to save structures as formatted archive
c     binary      logical flag to save structures as binary archive
c     noversion   logical flag governing use of filename versions
c     overwrite   logical flag to overwrite intermediate files inplace
c     cyclesave   logical flag to mark use of numbered cycle files
c     velsave     logical flag to save velocity vector components
c     frcsave     logical flag to save force vector components
c     uindsave    logical flag to save induced atomic dipoles
c     coordtype   selects Cartesian, internal, rigid body or none
c
c
      module output
      implicit none
      logical archive
      logical binary
      logical noversion
      logical overwrite
      logical cyclesave
      logical velsave
      logical frcsave
      logical uindsave
      character*9 coordtype
      save
      end

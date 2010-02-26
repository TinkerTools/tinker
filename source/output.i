c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  output.i  --  control of coordinate output file format  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     archive    logical flag to save structures in an archive
c     noversion  logical flag governing use of filename versions
c     overwrite  logical flag to overwrite intermediate files inplace
c     cyclesave  logical flag to mark use of numbered cycle files
c     coordtype  selects Cartesian, internal, rigid body or none
c
c
      logical archive,noversion
      logical overwrite,cyclesave
      character*9 coordtype
      common /output/ archive,noversion,overwrite,cyclesave,coordtype

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2023  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  module shapes  --  UnionBall area and volume variables  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxedge         maximum number of edges between ball centers
c     maxtetra        maximum number of tetrahedra in the system
c     npoint          total number of balls (points) in the system
c     nvertex         total number of vertices in the system
c     ntetra          total number of tetrahedra in the system
c     nnew            total number of entries on new tetrahedra list
c     nfree           total number of spaces on free tetrahedra list
c     nkill           total number of tetrahedra on list to kill
c     nlinkfacet      total number of triangle facets in the system
c     newlist         list with index numbers of the new tetrahedra
c     freespace       list of the tetrahedra currently in free space
c     killspace       list of the existing tetrahedra to be killed
c     vinfo           information value for each of the vertices
c     tedge           number of an edge found in each tetrahedron
c     tinfo           orientation information for each tetrahedron
c     tnindex         index related to tetrahedron orientation
c     tetra           numbers of the four balls in each tetrahedron
c     tneighbor       store the four neighbors of each tetrahedron
c     linkfacet       numbers of two tetrahedra defining each facet
c     linkindex       vertex numbers opposite each facet triangle
c     epsln2          minimal value of determinant over two balls
c     epsln3          minimal value of determinant over three balls
c     epsln4          minimal value of determinant over four balls
c     epsln5          minimal value of determinant over five balls
c     crdball         coordinates in Angstroms of balls as 1-D array
c     radball         radius value for each ball in Angstroms
c     wghtball        weight value assigned for each ball
c
c
      module shapes
      integer maxedge
      integer maxtetra
      integer npoint,nvertex
      integer ntetra,nnew
      integer nfree,nkill
      integer nlinkfacet
      integer, allocatable :: newlist(:)
      integer, allocatable :: freespace(:)
      integer, allocatable :: killspace(:)
      integer, allocatable :: vinfo(:)
      integer, allocatable :: tedge(:)
      integer, allocatable :: tinfo(:)
      integer, allocatable :: tnindex(:)
      integer, allocatable :: tetra(:,:)
      integer, allocatable :: tneighbor(:,:)
      integer, allocatable :: linkfacet(:,:)
      integer, allocatable :: linkindex(:,:)
      real*8 epsln2,epsln3
      real*8 epsln4,epsln5
      real*8, allocatable :: crdball(:)
      real*8, allocatable :: radball(:)
      real*8, allocatable :: wghtball(:)
      save
      end

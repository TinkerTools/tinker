c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program h2otest  --  water dimer distances and angles  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "h2otest" computes the O-O distance, the O-H distance, the
c     nonlinearity angle (H-O...O) and acceptor flap angle for
c     the canonical water dimer
c
c
      program h2otest
      use sizes
      use atoms
      use iounit
      implicit none
      real*8 geometry
      real*8 distOO,distOH
      real*8 angle,flap
c
c
c     read the structure for the water dimer
c
      call initial
      call getxyz
c
c     create an atom at midpoint of acceptor H-H
c
      n = 7
      x(7) = 0.5d0 * (x(5) + x(6))
      y(7) = 0.5d0 * (y(5) + y(6))
      z(7) = 0.5d0 * (z(5) + z(6))
c
c     compute the distance between oxygen atoms
c
      distOO = geometry (1,4,0,0)
c
c     compute the hydrogen bond O-H distance
c
      distOH = min(geometry(2,4,0,0),geometry(3,4,0,0))
c
c     compute the flap angle of the acceptor molecule
c
      angle = min(geometry(2,1,4,0),geometry(3,1,4,0))
c
c     compute the flap angle of the acceptor molecule
c
      flap = 180.0d0 - geometry (1,4,7,0)
c
c     print out the water dimer geometry values
c
      write (iout,10)  distOO,distOH,angle,flap
   10 format (/,' Dimer O-O Distance :',10x,f10.4,
     &        /,' Dimer O-H Distance :',10x,f10.4,
     &        /,' Nonlinearity Angle :',10x,f10.2,
     &        /,' Acceptor Flap Angle :',10x,f9.2)
      end

c
c
c     ################################################################
c     ## COPYRIGHT (C) 2013 by Xiao Zhu, Pengyu Ren & Jay W. Ponder ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  tarray.i  --  storage of dipole-dipole matrix elements  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     tdipdip    stored dipole-dipole matrix element values
c     tindex     index into stored dipole-dipole matrix values
c     ntpair     number of stored dipole-dipole matrix elements
c
c
      integer ntpair
      integer, pointer :: tindex(:,:)
      real*8, pointer :: tdipdip(:,:)
      common /tarray/ tdipdip,tindex,ntpair

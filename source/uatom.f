c
c
c     ##############################################################
c     ##  COPYRIGHT (C)  2024  by  Moses KJ Chung & Jay W Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     #######################################################
c     ##                                                   ##
c     ##  module uatom  --  properties based on atom type  ##
c     ##                                                   ##
c     #######################################################
c
c
c     nunique    number of unique atom types in the system
c     utype      map from unique type to atom type
c     utypeinv   map from atom type to unique type
c     utv1       unique type vector 1
c     utv2       unique type vector 2
c
c
      module uatom
      use sizes
      implicit none
      integer nunique
      integer utype(maxtyp)
      integer utypeinv(maxtyp)
      real*8 utv1(3,maxtyp)
      real*8 utv2(3,maxtyp)
      save
      end

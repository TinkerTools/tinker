c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  neigh.i  --  pairwise neighbor list indices and storage  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     xvold       x-coordinate at last vdw neighbor list update
c     yvold       y-coordinate at last vdw neighbor list update
c     zvold       z-coordinate at last vdw neighbor list update
c     xcold       x-coordinate at last charge neighbor list update
c     ycold       y-coordinate at last charge neighbor list update
c     zcold       z-coordinate at last charge neighbor list update
c     xmold       x-coordinate at last multipole neighbor list update
c     ymold       y-coordinate at last multipole neighbor list update
c     zmold       z-coordinate at last multipole neighbor list update
c     lbuffer     width of the neighbor list buffer region
c     lbuf2       square of half the neighbor list buffer width
c     vbuf2       square of vdw cutoff plus neighbor list buffer
c     cbuf2       square of charge cutoff plus neighbor list buffer
c     mbuf2       square of multipole cutoff plus neighbor list buffer
c     vbufx       square of vdw cutoff plus twice the list buffer
c     cbufx       square of charge cutoff plus twice the list buffer
c     mbufx       square of multipole cutoff plus twice the list buffer
c     nvlst       number of sites in list for each vdw site
c     vlst        site numbers in neighbor list of each vdw site
c     nelst       number of sites in list for each electrostatic site
c     elst        site numbers in list of each electrostatic site
c     dovlst      logical flag to rebuild vdw neighbor list
c     doclst      logical flag to rebuild charge neighbor list
c     domlst      logical flag to rebuild multipole neighbor list
c
c
      integer nvlst,vlst
      integer nelst,elst
      real*8 xvold,yvold,zvold
      real*8 xcold,ycold,zcold
      real*8 xmold,ymold,zmold
      real*8 lbuffer,lbuf2
      real*8 vbuf2,cbuf2,mbuf2
      real*8 vbufx,cbufx,mbufx
      logical dovlst,doclst,domlst
      common /neigh/ xvold(maxatm),yvold(maxatm),zvold(maxatm),
     &               xcold(maxatm),ycold(maxatm),zcold(maxatm),
     &               xmold(maxatm),ymold(maxatm),zmold(maxatm),
     &               lbuffer,lbuf2,vbuf2,cbuf2,mbuf2,vbufx,cbufx,
     &               mbufx,nvlst(maxatm),vlst(maxvlst,maxatm),
     &               nelst(maxatm),elst(maxelst,maxatm),dovlst,
     &               doclst,domlst

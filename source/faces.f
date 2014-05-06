c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module faces  --  Connolly area and volume variables  ##
c     ##                                                        ##
c     ############################################################
c
c
c     maxcls   maximum number of neighboring atom pairs
c     maxtt    maximum number of temporary tori
c     maxt     maximum number of total tori
c     maxp     maximum number of probe positions
c     maxv     maximum number of vertices
c     maxen    maximum number of concave edges
c     maxfn    maximum number of concave faces
c     maxc     maximum number of circles
c     maxep    maximum number of convex edges
c     maxfs    maximum number of saddle faces
c     maxcy    maximum number of cycles
c     mxcyep   maximum number of cycle convex edges
c     maxfp    maximum number of convex faces
c     mxfpcy   maximum number of convex face cycles
c
c
      module faces
      use sizes
      implicit none
      integer maxcls,maxtt
      integer maxt,maxp,maxv
      integer maxen,maxfn,maxc
      integer maxep,maxfs
      integer maxcy,mxcyep
      integer maxfp,mxfpcy
      parameter (maxcls=50*maxatm)
      parameter (maxtt=25*maxatm)
      parameter (maxt=3*maxatm)
      parameter (maxp=2*maxatm)
      parameter (maxv=5*maxatm)
      parameter (maxen=5*maxatm)
      parameter (maxfn=2*maxatm)
      parameter (maxc=5*maxatm)
      parameter (maxep=5*maxatm)
      parameter (maxfs=3*maxatm)
      parameter (maxcy=maxatm)
      parameter (mxcyep=30)
      parameter (maxfp=maxatm)
      parameter (mxfpcy=10)
c
c
c     na       number of atoms
c     pr       probe radius
c     ar       atomic radii
c     a        atomic coordinates
c
c
      integer na
      real*8 pr
      real*8 ar(maxatm)
      real*8 a(3,maxatm)
c
c
c     skip     if true, atom is not used
c     nosurf   if true, atom has no free surface
c     afree    atom free of neighbors
c     abur     atom buried
c
c
      logical skip(maxatm)
      logical nosurf(maxatm)
      logical afree(maxatm)
      logical abur(maxatm)
c
c
c     cls      atom numbers of neighbors
c     clst     pointer from neighbor to torus
c     acls     begin and end pointers for atoms neighbors
c
c
      integer cls(maxcls)
      integer clst(maxcls)
      integer acls(2,maxatm)
c
c
c     ntt      number of temporary tori
c     ttfe     first edge of each temporary torus
c     ttle     last edge of each temporary torus
c     enext    pointer to next edge of temporary torus
c     tta      temporary torus atom numbers
c     ttbur    temporary torus buried
c     ttfree   temporary torus free
c
c
      integer ntt
      integer ttfe(maxtt)
      integer ttle(maxtt)
      integer enext(maxen)
      integer tta(2,maxtt)
      logical ttbur(maxtt)
      logical ttfree(maxtt)
c
c
c     nt       number of tori
c     tfe      torus first edge
c     ta       torus atom numbers
c     tr       torus radius
c     t        torus center
c     tax      torus axis
c     tfree    torus free of neighbors
c
c
      integer nt
      integer tfe(maxt)
      integer ta(2,maxt)
      real*8 tr(maxt)
      real*8 t(3,maxt)
      real*8 tax(3,maxt)
      logical tfree(maxt)
c
c
c     np       number of probe positions
c     pa       probe position atom numbers
c     p        probe position coordinates
c
c
      integer np
      integer pa(3,maxp)
      real*8 p(3,maxp)
c
c
c     v        vertex coordinates
c     nv       number of vertices
c     va       vertex atom number
c     vp       vertex probe number
c
c
      integer nv
      integer va(maxv)
      integer vp(maxv)
      real*8 v(3,maxv)
c
c
c     nen      number of concave edges
c     nfn      number of concave faces
c     env      vertex numbers for each concave edge
c     fnen     concave face concave edge numbers
c
c
      integer nen
      integer nfn
      integer env(2,maxen)
      integer fnen(3,maxfn)
c
c
c     nc       number of circles
c     ca       circle atom number
c     ct       circle torus number
c     cr       circle radius
c     c        circle center
c
c
      integer nc
      integer ca(maxc)
      integer ct(maxc)
      real*8 cr(maxc)
      real*8 c(3,maxc)
c
c
c     nep      number of convex edges
c     epc      convex edge circle number
c     epv      convex edge vertex numbers
c     afe      first convex edge of each atom
c     ale      last convex edge of each atom
c     epnext   pointer to next convex edge of atom
c
c
      integer nep
      integer epc(maxep)
      integer epv(2,maxep)
      integer afe(maxatm)
      integer ale(maxatm)
      integer epnext(maxep)
c
c
c     nfs      number of saddle faces
c     fsen     saddle face concave edge numbers
c     fsep     saddle face convex edge numbers
c
c
      integer nfs
      integer fsen(2,maxfs)
      integer fsep(2,maxfs)
c
c
c     ncy      number of cycles
c     cynep    number of convex edges in cycle
c     cyep     cycle convex edge numbers
c
c
      integer ncy
      integer cynep(maxcy)
      integer cyep(mxcyep,maxcy)
c
c
c     nfp      number of convex faces
c     fpa      atom number of convex face
c     fpncy    number of cycles bounding convex face
c     fpcy     convex face cycle numbers
c
c
      integer nfp
      integer fpa(maxfp)
      integer fpncy(maxfp)
      integer fpcy(mxfpcy,maxfp)
      save
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  faces.i  --  variables for Connolly area and volume  ##
c     ##                                                       ##
c     ###########################################################
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
      integer maxcls,maxtt,maxt,maxp,maxv
      integer maxen,maxfn,maxc,maxep,maxfs
      integer maxcy,mxcyep,maxfp,mxfpcy
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
c     a        atomic coordinates
c     ar       atomic radii
c     pr       probe radius
c     na       number of atoms
c
c
      integer na
      real*8 a,ar,pr
      common /face01/ a(3,maxatm),ar(maxatm),pr,na
c
c
c     skip     if true, atom is not used
c     nosurf   if true, atom has no free surface
c     afree    atom free of neighbors
c     abur     atom buried
c
c
      logical skip,nosurf,afree,abur
      common /face02/ skip(maxatm),nosurf(maxatm),afree(maxatm),
     &                abur(maxatm)
c
c
c     acls     begin and end pointers for atoms neighbors
c     cls      atom numbers of neighbors
c     clst     pointer from neighbor to torus
c
c
      integer acls,cls,clst
      common /face03/ acls(2,maxatm),cls(maxcls),clst(maxcls)
c
c
c     ntt      number of temporary tori
c     tta      temporary torus atom numbers
c     ttfe     first edge of each temporary torus
c     ttle     last edge of each temporary torus
c     enext    pointer to next edge of temporary torus
c     ttbur    temporary torus buried
c     ttfree   temporary torus free
c
c
      integer ntt,tta,ttfe,ttle,enext
      logical ttbur,ttfree
      common /face04/ ntt,tta(2,maxtt),ttfe(maxtt),ttle(maxtt),
     &                enext(maxen),ttbur(maxtt),ttfree(maxtt)
c
c
c     t        torus center
c     tr       torus radius
c     tax      torus axis
c     nt       number of tori
c     ta       torus atom numbers
c     tfe      torus first edge
c     tfree    torus free of neighbors
c
c
      integer nt,ta,tfe
      real*8 t,tr,tax
      logical tfree
      common /face05/ t(3,maxt),tr(maxt),tax(3,maxt),nt,
     &                ta(2,maxt),tfe(maxt),tfree(maxt)
c
c
c     p        probe position coordinates
c     np       number of probe positions
c     pa       probe position atom numbers
c
c
      integer np,pa
      real*8 p
      common /face06/ p(3,maxp),np,pa(3,maxp)
c
c
c     v        vertex coordinates
c     nv       number of vertices
c     va       vertex atom number
c     vp       vertex probe number
c
c
      integer nv,va,vp
      real*8 v
      common /face07/ v(3,maxv),nv,va(maxv),vp(maxv)
c
c
c     nen      number of concave edges
c     env      vertex numbers for each concave edge
c     nfn      number of concave faces
c     fnen     concave face concave edge numbers
c
c
      integer nen,env,nfn,fnen
      common /face08/ nen,env(2,maxen),nfn,fnen(3,maxfn)
c
c
c     c        circle center
c     cr       circle radius
c     nc       number of circles
c     ca       circle atom number
c     ct       circle torus number
c
c
      integer nc,ca,ct
      real*8 c,cr
      common /face09/ c(3,maxc),cr(maxc),nc,ca(maxc),ct(maxc)
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
      integer nep,epc,epv,afe,ale,epnext
      common /face10/ nep,epc(maxep),epv(2,maxep),afe(maxatm),
     &                ale(maxatm),epnext(maxep)
c
c
c     nfs      number of saddle faces
c     fsen     saddle face concave edge numbers
c     fsep     saddle face convex edge numbers
c
c
      integer nfs,fsen,fsep
      common /face11/ nfs,fsen(2,maxfs),fsep(2,maxfs)
c
c
c     ncy      number of cycles
c     cynep    number of convex edges in cycle
c     cyep     cycle convex edge numbers
c
c
      integer ncy,cynep,cyep
      common /face12/ ncy,cynep(maxcy),cyep(mxcyep,maxcy)
c
c
c     nfp      number of convex faces
c     fpa      atom number of convex face
c     fpcy     convex face cycle numbers
c     fpncy    number of cycles bounding convex face
c
c
      integer nfp,fpa,fpcy,fpncy
      common /face13/ nfp,fpa(maxfp),fpcy(mxfpcy,maxfp),fpncy(maxfp)

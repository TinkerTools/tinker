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
c     maxeq    maximum number of convex edges
c     maxfs    maximum number of saddle faces
c     maxfq    maximum number of convex faces
c     maxcy    maximum number of cycles
c     mxcyeq   maximum number of convex edge cycles
c     mxfqcy   maximum number of convex face cycles
c
c
      module faces
      implicit none
      integer maxcls,maxtt
      integer maxt,maxp
      integer maxv,maxen
      integer maxfn,maxc
      integer maxeq,maxfs
      integer maxfq,maxcy
      integer mxcyeq,mxfqcy
c
c
c     na       number of atoms
c     pr       probe radius
c     ar       atomic radii
c     axyz     atomic coordinates
c
c
      integer na
      real*8 pr
      real*8, allocatable :: ar(:)
      real*8, allocatable :: axyz(:,:)
c
c
c     skip     if true, atom is not used
c     nosurf   if true, atom has no free surface
c     afree    atom free of neighbors
c     abur     atom buried
c
c
      logical, allocatable :: skip(:)
      logical, allocatable :: nosurf(:)
      logical, allocatable :: afree(:)
      logical, allocatable :: abur(:)
c
c
c     cls      atom numbers of neighbors
c     clst     pointer from neighbor to torus
c     acls     begin and end pointers for atoms neighbors
c
c
      integer, allocatable :: cls(:)
      integer, allocatable :: clst(:)
      integer, allocatable :: acls(:,:)
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
      integer, allocatable :: ttfe(:)
      integer, allocatable :: ttle(:)
      integer, allocatable :: enext(:)
      integer, allocatable :: tta(:,:)
      logical, allocatable :: ttbur(:)
      logical, allocatable :: ttfree(:)
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
      integer, allocatable :: tfe(:)
      integer, allocatable :: ta(:,:)
      real*8, allocatable :: tr(:)
      real*8, allocatable :: t(:,:)
      real*8, allocatable :: tax(:,:)
      logical, allocatable :: tfree(:)
c
c
c     np       number of probe positions
c     pa       probe position atom numbers
c     p        probe position coordinates
c
c
      integer np
      integer, allocatable :: pa(:,:)
      real*8, allocatable :: p(:,:)
c
c
c     nv       number of vertices
c     va       vertex atom number
c     vp       vertex probe number
c     vxyz     vertex coordinates
c
c
      integer nv
      integer, allocatable :: va(:)
      integer, allocatable :: vp(:)
      real*8, allocatable :: vxyz(:,:)
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
      integer, allocatable :: env(:,:)
      integer, allocatable :: fnen(:,:)
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
      integer, allocatable :: ca(:)
      integer, allocatable :: ct(:)
      real*8, allocatable :: cr(:)
      real*8, allocatable :: c(:,:)
c
c
c     neq      number of convex edges
c     eqc      convex edge circle number
c     eqv      convex edge vertex numbers
c     afe      first convex edge of each atom
c     ale      last convex edge of each atom
c     eqnext   pointer to next convex edge of atom
c
c
      integer neq
      integer, allocatable :: eqc(:)
      integer, allocatable :: eqv(:,:)
      integer, allocatable :: afe(:)
      integer, allocatable :: ale(:)
      integer, allocatable :: eqnext(:)
c
c
c     nfs      number of saddle faces
c     fsen     saddle face concave edge numbers
c     fseq     saddle face convex edge numbers
c
c
      integer nfs
      integer, allocatable :: fsen(:,:)
      integer, allocatable :: fseq(:,:)
c
c
c     ncy      number of cycles
c     cyneq    number of convex edges in cycle
c     cyeq     cycle convex edge numbers
c
c
      integer ncy
      integer, allocatable :: cyneq(:)
      integer, allocatable :: cyeq(:,:)
c
c
c     nfq      number of convex faces
c     fqa      atom number of convex face
c     fqncy    number of cycles bounding convex face
c     fqcy     convex face cycle numbers
c
c
      integer nfq
      integer, allocatable :: fqa(:)
      integer, allocatable :: fqncy(:)
      integer, allocatable :: fqcy(:,:)
      save
      end

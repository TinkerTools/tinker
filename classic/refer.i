c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  refer.i  --  storage of reference atomic coordinate set  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     xref        reference x-coordinates for atoms in each system
c     yref        reference y-coordinates for atoms in each system
c     zref        reference z-coordinates for atoms in each system
c     nref        total number of atoms in each reference system
c     reftyp      atom types of the atoms in each reference system
c     n12ref      number of atoms bonded to each reference atom
c     i12ref      atom numbers of atoms 1-2 connected to each atom
c     refleng     length in characters of each reference filename
c     refltitle   length in characters of each reference title line
c     refnam      atom names of the atoms in each reference system
c     reffile     base filename for each reference system
c     reftitle    title used to describe each reference system
c
c
      integer nref,reftyp
      integer n12ref,i12ref
      integer refleng,refltitle
      real*8 xref,yref,zref
      character*3 refnam
      character*120 reffile,reftitle
      common /refer/ xref(maxatm,maxref),yref(maxatm,maxref),
     &               zref(maxatm,maxref),nref(maxref),
     &               reftyp(maxatm,maxref),n12ref(maxatm,maxref),
     &               i12ref(maxval,maxatm,maxref),refleng(maxref),
     &               refltitle(maxref),refnam(maxatm,maxref),
     &               reffile(maxref),reftitle(maxref)

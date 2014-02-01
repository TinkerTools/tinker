c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  group.i  --  partitioning of system into atom groups  ##
c     ##                                                        ##
c     ############################################################
c
c
c     grpmass     total mass of all the atoms in each group
c     wgrp        weight for each set of group-group interactions
c     ngrp        total number of atom groups in the system
c     kgrp        contiguous list of the atoms in each group
c     igrp        first and last atom of each group in the list
c     grplist     number of the group to which each atom belongs
c     use_group   flag to use partitioning of system into groups
c     use_intra   flag to include only intragroup interactions
c     use_inter   flag to include only intergroup interactions
c
c
      integer ngrp,kgrp
      integer igrp,grplist
      real*8 grpmass,wgrp
      logical use_group
      logical use_intra
      logical use_inter
      common /group/ grpmass(maxgrp),wgrp(0:maxgrp,0:maxgrp),ngrp,
     &               kgrp(maxatm),igrp(2,0:maxgrp),grplist(maxatm),
     &               use_group,use_intra,use_inter

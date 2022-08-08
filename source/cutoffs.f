c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine cutoffs  --  distance cutoffs & neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "cutoffs" initializes and stores spherical energy cutoff
c     distance windows, Hessian element and Ewald sum cutoffs,
c     and allocates pairwise neighbor lists
c
c
      subroutine cutoffs
      use atoms
      use bound
      use hescut
      use keys
      use limits
      use neigh
      use polpot
      use tarray
      implicit none
      integer i,next
      integer limit
      real*8 big,value
      logical truncate
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for spherical energy cutoff distances
c
      big = 1.0d12
      if (use_bounds) then
         vdwcut = 9.0d0
         dispcut = 9.0d0
         chgcut = 9.0d0
         dplcut = 9.0d0
         mpolecut = 9.0d0
      else
         vdwcut = big
         dispcut = big
         chgcut = big
         dplcut = big
         mpolecut = big
      end if
      repcut = 6.0d0
      ctrncut = 6.0d0
      ewaldcut = 7.0d0
      dewaldcut = 7.0d0
      usolvcut = 4.5d0
c
c     set defaults for tapering, Hessian cutoff and neighbor buffers
c
      vdwtaper = 0.90d0
      reptaper = 0.90d0
      disptaper = 0.90d0
      chgtaper = 0.65d0
      dpltaper = 0.75d0
      mpoletaper = 0.65d0
      ctrntaper = 0.90d0
      hesscut = 0.0d0
      lbuffer = 2.0d0
      pbuffer = 2.0d0
c
c     set defaults for Ewald sum, tapering style and neighbor method
c
      use_ewald = .false.
      use_dewald = .false.
      truncate = .false.
      use_lights = .false.
      use_list = .false.
      use_vlist = .false.
      use_dlist = .false.
      use_clist = .false.
      use_mlist = .false.
      use_ulist = .false.
      dovlst = .true.
      dodlst = .true.
      doclst = .true.
      domlst = .true.
      doulst = .true.
c
c     search the keywords for various cutoff parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
c
c     get values related to use of Ewald for electrostatics
c
         if (keyword(1:6) .eq. 'EWALD ') then
            use_ewald = .true.
         else if (keyword(1:13) .eq. 'EWALD-CUTOFF ') then
            read (string,*,err=10,end=10)  ewaldcut
c
c     get values related to use of Ewald for dispersion
c
         else if (keyword(1:7) .eq. 'DEWALD ') then
            use_dewald = .true.
         else if (keyword(1:14) .eq. 'DEWALD-CUTOFF ') then
            read (string,*,err=10,end=10)  dewaldcut
c
c     get values for the tapering style and neighbor method
c
         else if (keyword(1:9) .eq. 'TRUNCATE ') then
            truncate = .true.
         else if (keyword(1:7) .eq. 'LIGHTS ') then
            use_lights = .true.
         else if (keyword(1:14) .eq. 'NEIGHBOR-LIST ') then
            use_list = .true.
            use_vlist = .true.
            use_dlist = .true.
            use_clist = .true.
            use_mlist = .true.
            use_ulist = .true.
         else if (keyword(1:9) .eq. 'VDW-LIST ') then
            use_list = .true.
            use_vlist = .true.
         else if (keyword(1:10) .eq. 'DISP-LIST ') then
            use_list = .true.
            use_dlist = .true.
         else if (keyword(1:12) .eq. 'CHARGE-LIST ') then
            use_list = .true.
            use_clist = .true.
         else if (keyword(1:11) .eq. 'MPOLE-LIST ') then
            use_list = .true.
            use_mlist = .true.
            use_ulist = .true.
c
c     get values for the dipole solver preconditioner
c
         else if (keyword(1:12) .eq. 'USOLVE-LIST ') then
            use_list = .true.
            use_ulist = .true.
         else if (keyword(1:14) .eq. 'USOLVE-CUTOFF ') then
            if (usolvcut .ne. 0.0d0) then
               read (string,*,err=10,end=10)  usolvcut
            end if
         else if (keyword(1:16) .eq. 'USOLVE-DIAGONAL ') then
            usolvcut = 0.0d0
c
c     get cutoff for the magnitude of Hessian elements
c
         else if (keyword(1:15) .eq. 'HESSIAN-CUTOFF ') then
            read (string,*,err=10,end=10)  hesscut
c
c     get the cutoff radii for potential energy functions
c
         else if (keyword(1:7) .eq. 'CUTOFF ') then
            read (string,*,err=10,end=10)  value
            vdwcut = value
            repcut = value
            dispcut = value
            chgcut = value
            dplcut = value
            mpolecut = value
            ewaldcut = value
            dewaldcut = value
            ctrncut = value
         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
            read (string,*,err=10,end=10)  vdwcut
         else if (keyword(1:14) .eq. 'REPULS-CUTOFF ') then
            read (string,*,err=10,end=10)  repcut
         else if (keyword(1:12) .eq. 'DISP-CUTOFF ') then
            read (string,*,err=10,end=10)  dispcut
         else if (keyword(1:14) .eq. 'CHARGE-CUTOFF ') then
            read (string,*,err=10,end=10)  chgcut
         else if (keyword(1:14) .eq. 'DIPOLE-CUTOFF ') then
            read (string,*,err=10,end=10)  dplcut
         else if (keyword(1:13) .eq. 'MPOLE-CUTOFF ') then
            read (string,*,err=10,end=10)  mpolecut
         else if (keyword(1:14) .eq. 'CHGTRN-CUTOFF ') then
            read (string,*,err=10,end=10)  ctrncut
c
c     get distance for initialization of energy switching
c
         else if (keyword(1:6) .eq. 'TAPER ') then
            read (string,*,err=10,end=10)  value
            vdwtaper = value
            reptaper = value
            disptaper = value
            chgtaper = value
            dpltaper = value
            mpoletaper = value
            ctrntaper = value
         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
            read (string,*,err=10,end=10)  vdwtaper
         else if (keyword(1:13) .eq. 'REPULS-TAPER ') then
            read (string,*,err=10,end=10)  reptaper
         else if (keyword(1:11) .eq. 'DISP-TAPER ') then
            read (string,*,err=10,end=10)  disptaper
         else if (keyword(1:13) .eq. 'CHARGE-TAPER ') then
            read (string,*,err=10,end=10)  chgtaper
         else if (keyword(1:13) .eq. 'DIPOLE-TAPER ') then
            read (string,*,err=10,end=10)  dpltaper
         else if (keyword(1:12) .eq. 'MPOLE-TAPER ') then
            read (string,*,err=10,end=10)  mpoletaper
         else if (keyword(1:13) .eq. 'CHGTRN-TAPER ') then
            read (string,*,err=10,end=10)  ctrntaper
c
c     get buffer width for use with pairwise neighbor lists
c
         else if (keyword(1:12) .eq. 'LIST-BUFFER ') then
            read (string,*,err=10,end=10)  lbuffer
         else if (keyword(1:14) .eq. 'USOLVE-BUFFER ') then
            read (string,*,err=10,end=10)  pbuffer
         end if
   10    continue
      end do
c
c     check to see if preconditioner list should be disabled
c
      if (poltyp .eq. 'DIRECT')  use_ulist = .false.
      if (usolvcut .le. 0.0d0)  use_ulist = .false.
      if (use_list)  usolvcut = usolvcut - pbuffer
c
c     apply any Ewald cutoff to dispersion and electrostatics
c
      if (use_ewald) then
         chgcut = ewaldcut
         mpolecut = ewaldcut
      end if
      if (use_dewald) then
         dispcut = dewaldcut
      end if
c
c     convert any tapering percentages to absolute distances
c
      if (vdwtaper .lt. 1.0d0)  vdwtaper = vdwtaper * vdwcut
      if (reptaper .lt. 1.0d0)  reptaper = reptaper * repcut
      if (disptaper .lt. 1.0d0)  disptaper = disptaper * dispcut
      if (chgtaper .lt. 1.0d0)  chgtaper = chgtaper * chgcut
      if (dpltaper .lt. 1.0d0)  dpltaper = dpltaper * dplcut
      if (mpoletaper .lt. 1.0d0)  mpoletaper = mpoletaper * mpolecut
      if (ctrntaper .lt. 1.0d0)  ctrntaper = ctrntaper * ctrncut
c
c     apply truncation cutoffs if they were requested
c
      if (truncate) then
         vdwtaper = big
         reptaper = big
         disptaper = big
         chgtaper = big
         dpltaper = big
         mpoletaper = big
         ctrntaper = big
      end if
c
c     set buffer region limits for pairwise neighbor lists
c
      lbuf2 = (0.5d0*lbuffer)**2
      pbuf2 = (0.5d0*pbuffer)**2
      vbuf2 = (vdwcut+lbuffer)**2
      dbuf2 = (dispcut+lbuffer)**2
      cbuf2 = (chgcut+lbuffer)**2
      mbuf2 = (mpolecut+lbuffer)**2
      ubuf2 = (usolvcut+pbuffer)**2
      vbufx = (vdwcut+2.0d0*lbuffer)**2
      dbufx = (dispcut+2.0d0*lbuffer)**2
      cbufx = (chgcut+2.0d0*lbuffer)**2
      mbufx = (mpolecut+2.0d0*lbuffer)**2
      ubufx = (usolvcut+2.0d0*pbuffer)**2
c
c     specify maximum size for each of the neighbor lists
c
      maxvlst = 2500
      if (vdwcut.ne.big .and. dispcut.ne.big) then
         limit = int(sqrt(max(vbuf2,dbuf2))**3) + 100
         maxvlst = min(limit,maxvlst)
      else if (vdwcut .ne. big) then
         limit = int(sqrt(vbuf2)**3) + 100
         maxvlst = min(limit,maxvlst)
      else if (dispcut .ne. big) then
         limit = int(sqrt(dbuf2)**3) + 100
         maxvlst = min(limit,maxvlst)
      end if
      maxelst = 2500
      if (chgcut.ne.big .and. mpolecut.ne.big) then
         limit = int(sqrt(max(cbuf2,mbuf2))**3) + 100
         maxelst = min(limit,maxelst)
      else if (chgcut .ne. big) then
         limit = int(sqrt(cbuf2)**3) + 100
         maxelst = min(limit,maxelst)
      else if (mpolecut .ne. big) then
         limit = int(sqrt(mbuf2)**3) + 100
         maxelst = min(limit,maxelst)
      end if
      maxulst = 500
      limit = int(sqrt(ubuf2)**3) + 100
      maxulst = min(limit,maxulst)
c
c     perform dynamic allocation of some global arrays
c
      if (use_vlist .or. use_dlist) then
         if (allocated(nvlst))  deallocate (nvlst)
         if (allocated(vlst))  deallocate (vlst)
         if (allocated(xvold))  deallocate (xvold)
         if (allocated(yvold))  deallocate (yvold)
         if (allocated(zvold))  deallocate (zvold)
         allocate (nvlst(n))
         allocate (vlst(maxvlst,n))
         allocate (xvold(n))
         allocate (yvold(n))
         allocate (zvold(n))
      end if
      if (use_clist .or. use_mlist) then
         if (allocated(nelst))  deallocate (nelst)
         if (allocated(elst))  deallocate (elst)
         if (allocated(xeold))  deallocate (xeold)
         if (allocated(yeold))  deallocate (yeold)
         if (allocated(zeold))  deallocate (zeold)
         allocate (nelst(n))
         allocate (elst(maxelst,n))
         allocate (xeold(n))
         allocate (yeold(n))
         allocate (zeold(n))
         if (poltyp .ne. 'DIRECT') then
            if (allocated(tindex))  deallocate (tindex)
            if (allocated(tdipdip))  deallocate (tdipdip)
            allocate (tindex(2,n*maxelst))
            allocate (tdipdip(6,n*maxelst))
         end if
      end if
      if (use_ulist) then
         if (allocated(nulst))  deallocate (nulst)
         if (allocated(ulst))  deallocate (ulst)
         if (allocated(xuold))  deallocate (xuold)
         if (allocated(yuold))  deallocate (yuold)
         if (allocated(zuold))  deallocate (zuold)
         allocate (nulst(n))
         allocate (ulst(maxulst,n))
         allocate (xuold(n))
         allocate (yuold(n))
         allocate (zuold(n))
      end if
      return
      end

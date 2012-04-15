c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine hessian  --  atom-by-atom Hessian elements  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "hessian" calls subroutines to calculate the Hessian elements
c     for each atom in turn with respect to Cartesian coordinates
c
c
      subroutine hessian (h,hinit,hstop,hindex,hdiag)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'couple.i'
      include 'cutoff.i'
      include 'hescut.i'
      include 'hessn.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'potent.i'
      include 'rigid.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k
      integer ii,nhess
      integer hindex(*)
      integer hinit(3,*)
      integer hstop(3,*)
      real*8 percent,filled
      real*8 cutoff,rdn,hmax
      real*8 h(*)
      real*8 hdiag(3,*)
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      logical, allocatable :: keep(:)
c
c
c     zero out total number of indexed Hessian elements
c
      nhess = 0
      do i = 1, n
         do j = 1, 3
            hinit(j,i) = 1
            hstop(j,i) = 0
            hdiag(j,i) = 0.0d0
         end do
      end do
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_rigid)  call bounds
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc
c
c     compute the induced dipoles at polarizable atoms
c
      if (use_polar) then
         call chkpole
         call rotpole
         call induce
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (keep(n))
c
c     calculate the reduced atomic coordinates
c
      if (use_vdw) then
         do i = 1, n
            ii = ired(i)
            rdn = kred(i)
            xred(i) = rdn*(x(i)-x(ii)) + x(ii)
            yred(i) = rdn*(y(i)-y(ii)) + y(ii)
            zred(i) = rdn*(z(i)-z(ii)) + z(ii)
         end do
      end if
c
c     zero out the Hessian elements for the current atom
c
      do i = 1, n
         if (use(i)) then
            do k = 1, n
               do j = 1, 3
                  hessx(j,k) = 0.0d0
                  hessy(j,k) = 0.0d0
                  hessz(j,k) = 0.0d0
               end do
            end do
c
c     remove any previous use of the replicates method
c
            cutoff = 0.0d0
            call replica (cutoff)
c
c     call the local geometry Hessian component routines
c
            if (use_bond)  call ebond2 (i)
            if (use_angle)  call eangle2 (i)
            if (use_strbnd)  call estrbnd2 (i)
            if (use_urey)  call eurey2 (i)
            if (use_angang)  call eangang2 (i)
            if (use_opbend)  call eopbend2 (i)
            if (use_opdist)  call eopdist2 (i)
            if (use_improp)  call eimprop2 (i)
            if (use_imptor)  call eimptor2 (i)
            if (use_tors)  call etors2 (i)
            if (use_pitors)  call epitors2 (i)
            if (use_strtor)  call estrtor2 (i)
            if (use_tortor)  call etortor2 (i)
c
c     call the van der Waals Hessian component routines
c
            if (use_vdw) then
               if (vdwtyp .eq. 'LENNARD-JONES') then
                  call elj2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'BUCKINGHAM') then
                  call ebuck2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'MM3-HBOND') then
                  call emm3hb2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'BUFFERED-14-7') then
                  call ehal2 (i,xred,yred,zred)
               else if (vdwtyp .eq. 'GAUSSIAN') then
                  call egauss2 (i,xred,yred,zred)
               end if
            end if
c
c     call the electrostatic Hessian component routines
c
            if (use_charge) call echarge2 (i)
            if (use_chgdpl)  call echgdpl2 (i)
            if (use_dipole)  call edipole2 (i)
            if (use_mpole .or. use_polar)   call empole2 (i)
            if (use_rxnfld)   call erxnfld2 (i)
c
c     call any miscellaneous Hessian component routines
c
            if (use_solv)  call esolv2 (i)
            if (use_metal)  call emetal2 (i)
            if (use_geom)  call egeom2 (i)
            if (use_extra)  call extra2 (i)
c
c     set the diagonal Hessian matrix elements
c
            hdiag(1,i) = hdiag(1,i) + hessx(1,i)
            hdiag(2,i) = hdiag(2,i) + hessy(2,i)
            hdiag(3,i) = hdiag(3,i) + hessz(3,i)
c
c     search each 3x3 block to see which blocks will be kept
c
            do k = i+1, n
               keep(k) = .false.
               if (use(k)) then
                  hmax = max(abs(hessx(1,k)),abs(hessx(2,k)),
     &                       abs(hessx(3,k)),abs(hessy(1,k)),
     &                       abs(hessy(2,k)),abs(hessy(3,k)),
     &                       abs(hessz(1,k)),abs(hessz(2,k)),
     &                       abs(hessz(3,k)))
                  if (hmax .ge. hesscut)  keep(k) = .true.
               end if
            end do
c
c     copy selected off-diagonal Hessian elements for current
c     atom into an indexed master list of Hessian elements;
c     if any elements of 3x3 block are kept, keep them all
c
            hinit(1,i) = nhess + 1
            do j = 2, 3
               nhess = nhess + 1
               hindex(nhess) = 3*i + j - 3
               h(nhess) = hessx(j,i)
            end do
            do k = i+1, n
               if (keep(k)) then
                  do j = 1, 3
                     nhess = nhess + 1
                     hindex(nhess) = 3*k + j - 3
                     h(nhess) = hessx(j,k)
                  end do
               end if
            end do
            hstop(1,i) = nhess
            hinit(2,i) = nhess + 1
            nhess = nhess + 1
            hindex(nhess) = 3*i
            h(nhess) = hessy(3,i)
            do k = i+1, n
               if (keep(k)) then
                  do j = 1, 3
                     nhess = nhess + 1
                     hindex(nhess) = 3*k + j - 3
                     h(nhess) = hessy(j,k)
                  end do
               end if
            end do
            hstop(2,i) = nhess
            hinit(3,i) = nhess + 1
            do k = i+1, n
               if (keep(k)) then
                  do j = 1, 3
                     nhess = nhess + 1
                     hindex(nhess) = 3*k + j - 3
                     h(nhess) = hessz(j,k)
                  end do
               end if
            end do
            hstop(3,i) = nhess
c
c     check for storage of too many Hessian elements
c
            if (nhess .gt. maxhess) then
               write (iout,10)  i,n,nhess,maxhess,hesscut
   10          format (/,' Current Atom :',i7,6x,'Total Atoms :',i7,
     &                 /,' Current Required Hessian Storage : ',i12,
     &                 /,' Maximum Allowed Hessian Storage :  ',i12,
     &                 /,' Minimum Significant Hessian Value :',f12.5)
               write (iout,20)
   20          format (/,' HESSIAN  --  Increase MAXHESS',
     &                    ' and/or HESSCUT')
               call fatal
            end if
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (keep)
c
c     print message telling how much storage was finally used
c
      if (verbose) then
         percent = 100.0d0 * dble(nhess)/dble(3*n*(3*n-1)/2)
         filled = 100.0d0 * dble(nhess)/dble(maxhess)
         write (iout,30)  nhess,percent,filled
   30    format (' HESSIAN  --',i11,' Elements',f9.2,
     &              ' % Off-diag H',f8.2,' % Storage')
      end if
      return
      end

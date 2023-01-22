c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine echgtrn3  --  charge transfer energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echgtrn3" calculates the charge transfer energy; also partitions
c     the energy among the atoms
c
c
      subroutine echgtrn3
      use limits
      implicit none
c
c
c     choose method for summing over charge transfer interactions
c
      if (use_mlist) then
         call echgtrn3c
      else if (use_lights) then
         call echgtrn3b
      else
         call echgtrn3a
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine echgtrn3a  --  double loop chgtrn analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "echgtrn3a" calculates the charge transfer interaction energy
c     and also partitions the energy among the atoms using a pairwise
c     double loop
c
c
      subroutine echgtrn3a
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpot
      use chgtrn
      use cell
      use couple
      use ctrpot
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mplpot
      use mpole
      use mutant
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 e,f,fgrp
      real*8 r,r2,r3
      real*8 r4,r5
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 chgi,chgk
      real*8 chgik
      real*8 alphai,alphak
      real*8 expi,expk
      real*8 expik
      real*8 taper
      real*8, allocatable :: mscale(:)
      logical proceed,usei
      logical muti,mutk
      logical header,huge
      character*6 mode
c
c
c     zero out the charge transfer energy and partitioning terms
c
      nect = 0
      ect = 0.0d0
      do i = 1, n
         aect(i) = 0.0d0
      end do
      if (nct .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHGTRN'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nct.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge Transfer Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     calculate the charge transfer energy term
c
      do ii = 1, npole-1
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         chgi = chgct(ii)
         alphai = dmpct(ii)
         if (alphai .eq. 0.0d0)  alphai = 1000.0d0
         usei = use(i)
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = m2scale
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = m3scale
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = m4scale
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  chgk = chgct(kk)
                  alphak = dmpct(kk)
                  if (alphak .eq. 0.0d0)  alphak = 1000.0d0
                  if (ctrntyp .eq. 'SEPARATE') then
                     expi = exp(-alphai*r)
                     expk = exp(-alphak*r)
                     e = -chgi*expk - chgk*expi
                  else
                     chgik = sqrt(abs(chgi*chgk))
                     expik = exp(-0.5d0*(alphai+alphak)*r)
                     e = -chgik * expik
                  end if
                  e = f * e * mscale(k)
c
c     apply lambda scaling for interaction annihilation 
c
                  if (muti .or. mutk) then
                     e = e * elambda
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c     
c     increment the overall charge transfer energy components
c
                  if (e .ne. 0.0d0) then
                     nect = nect + 1
                     ect = ect + e
                     aect(i) = aect(i) + 0.5d0*e
                     aect(k) = aect(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Charge Transfer',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,e
   30                format (' ChgTrn',4x,2(i7,'-',a3),
     &                          9x,f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            chgi = chgct(ii)
            alphai = dmpct(ii)
            if (alphai .eq. 0.0d0)  alphai = 1000.0d0
            usei = use(i)
            muti = mut(i)
c
c     set exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               mscale(i12(j,i)) = m2scale
            end do
            do j = 1, n13(i)
               mscale(i13(j,i)) = m3scale
            end do
            do j = 1, n14(i)
               mscale(i14(j,i)) = m4scale
            end do
            do j = 1, n15(i)
               mscale(i15(j,i)) = m5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do kk = i, npole
               k = ipole(kk)
               mutk = mut(k)
               proceed = .true.
               if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
               if (.not. use_intra)  proceed = .true.
               if (proceed)  proceed = (usei .or. use(k))
               if (proceed) then
                  do jcell = 2, ncell
                     xr = x(k) - xi
                     yr = y(k) - yi
                     zr = z(k) - zi
                     call imager (xr,yr,zr,jcell)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        chgk = chgct(kk)
                        alphak = dmpct(kk)
                        if (alphak .eq. 0.0d0)  alphak = 1000.0d0
                        if (ctrntyp .eq. 'SEPARATE') then
                           expi = exp(-alphai*r)
                           expk = exp(-alphak*r)
                           e = -chgi*expk - chgk*expi
                        else
                           chgik = sqrt(abs(chgi*chgk))
                           expik = exp(-0.5d0*(alphai+alphak)*r)
                           e = -chgik * expik
                        end if
                        e = f * e * mscale(k)
c
c     apply lambda scaling for interaction annihilation 
c
                        if (muti .or. mutk) then
                           e = e * elambda
                        end if
c
c     use energy switching if near the cutoff distance
c
                        if (r2 .gt. cut2) then
                           r3 = r2 * r
                           r4 = r2 * r2
                           r5 = r2 * r3
                           taper = c5*r5 + c4*r4 + c3*r3
     &                                + c2*r2 + c1*r + c0
                           e = e * taper
                        end if
c
c     scale the interaction based on its group membership
c
                        if (use_group)  e = e * fgrp
c     
c     increment the overall charge transfer energy components
c
                        if (e .ne. 0.0d0) then
                           nect = nect + 1
                           if (i .eq. k) then
                              ect = ect + 0.5d0*e
                              aect(i) = aect(i) + 0.5d0*e
                           else
                              ect = ect + e
                              aect(i) = aect(i) + 0.5d0*e
                              aect(k) = aect(k) + 0.5d0*e
                           end if
                        end if
c
c     increment the total intermolecular energy
c
                        if (molcule(i) .ne. molcule(k)) then
                           einter = einter + e
                        end if
c
c     print a message if the energy of this interaction is large
c
                        huge = (e .gt. 10.0d0)
                        if ((debug.and.e.ne.0.0d0)
     &                        .or. (verbose.and.huge)) then
                           if (header) then
                              header = .false.
                              write (iout,40)
   40                         format (/,' Individual Charge Transfer',
     &                                   ' Interactions :',
     &                                //,' Type',14x,'Atom Names',
     &                                   15x,'Distance',8x,'Energy',/)
                           end if
                           write (iout,50)  i,name(i),k,name(k),r,e
   50                      format (' ChgTrn',4x,2(i7,'-',a3),2x,
     &                                '(XTAL)',1x,f10.4,2x,f12.4)
                        end if
                     end if
                  end do
               end if
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               mscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               mscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               mscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               mscale(i15(j,i)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echgtrn3b  --  method lights chgtrn analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echgtrn3b" calculates the charge transfer interaction energy
c     and also partitions the energy among the atoms using the method
c     of lights
c
c
      subroutine echgtrn3b
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use chgpot
      use chgtrn
      use cell
      use couple
      use ctrpot
      use energi
      use group
      use inform
      use inter
      use iounit
      use light
      use molcul
      use mplpot
      use mpole
      use mutant
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer kgy,kgz
      integer start,stop
      real*8 e,f,fgrp
      real*8 r,r2,r3
      real*8 r4,r5
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 chgi,chgk
      real*8 chgik
      real*8 alphai,alphak
      real*8 expi,expk
      real*8 expik
      real*8 taper
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei
      logical unique,repeat
      logical muti,mutk
      logical header,huge
      character*6 mode
c
c
c     zero out the charge transfer energy
c
      nect = 0
      ect = 0.0d0
      do i = 1, n
         aect(i) = 0.0d0
      end do
      if (nct .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (xsort(8*n))
      allocate (ysort(8*n))
      allocate (zsort(8*n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHGTRN'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nct.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge Transfer Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     transfer the interaction site coordinates to sorting arrays
c
      do ii = 1, npole
         i = ipole(ii)
         xsort(i) = x(i)
         ysort(i) = y(i)
         zsort(i) = z(i)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .true.
      call lights (off,nct,xsort,ysort,zsort,unique)
c
c     calculate the charge transfer energy term
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         chgi = chgct(ii)
         alphai = dmpct(ii)
         if (alphai .eq. 0.0d0)  alphai = 1000.0d0
         usei = use(i)
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = m2scale
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = m3scale
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = m4scale
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = m5scale
         end do
c
c     loop over method of lights neighbors of current atom
c
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   20    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 50
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 50
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 50
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 50
            end if
            k = ipole(kk-((kk-1)/npole)*npole)
            mutk = mut(k)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_bounds) then
                  if (abs(xr) .gt. xcell2)  xr = xr - sign(xcell,xr)
                  if (abs(yr) .gt. ycell2)  yr = yr - sign(ycell,yr)
                  if (abs(zr) .gt. zcell2)  zr = zr - sign(zcell,zr)
                  if (monoclinic) then
                     xr = xr + zr*beta_cos
                     zr = zr * beta_sin
                  else if (triclinic) then
                     xr = xr + yr*gamma_cos + zr*beta_cos
                     yr = yr*gamma_sin + zr*beta_term
                     zr = zr * gamma_term
                  end if
               end if
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  chgk = chgct(kk)
                  alphak = dmpct(kk)
                  if (alphak .eq. 0.0d0)  alphak = 1000.0d0
                  if (ctrntyp .eq. 'SEPARATE') then
                     expi = exp(-alphai*r)
                     expk = exp(-alphak*r)
                     e = -chgi*expk - chgk*expi
                  else
                     chgik = sqrt(abs(chgi*chgk))
                     expik = exp(-0.5d0*(alphai+alphak)*r)
                     e = -chgik * expik
                  end if
                  e = f * e * mscale(k)
c
c     apply lambda scaling for interaction annihilation 
c
                  if (muti .or. mutk) then
                     e = e * elambda
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c     
c     increment the overall charge transfer energy components
c
                  if (e .ne. 0.0d0) then
                     nect = nect + 1
                     ect = ect + e
                     aect(i) = aect(i) + 0.5d0*e
                     aect(k) = aect(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,30)
   30                   format (/,' Individual Charge Transfer',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,40)  i,name(i),k,name(k),r,e
   40                format (' ChgTrn',4x,2(i7,'-',a3),
     &                          9x,f10.4,2x,f12.4)
                  end if
               end if
            end if
   50       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 20
         end if
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echgtrn3c  --  neighbor list chgtrn analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echgtrn3c" calculates the charge transfer interaction energy
c     and also partitions the energy among the atoms using a pairwise
c     neighbor list
c
c
      subroutine echgtrn3c
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpot
      use chgtrn
      use cell
      use couple
      use ctrpot
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use mplpot
      use mpole
      use mutant
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,f,fgrp
      real*8 r,r2,r3
      real*8 r4,r5
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 chgi,chgk
      real*8 chgik
      real*8 alphai,alphak
      real*8 expi,expk
      real*8 expik
      real*8 taper
      real*8, allocatable :: mscale(:)
      logical proceed,usei
      logical muti,mutk
      logical header,huge
      character*6 mode
c
c     zero out the charge transfer energy and partitioning terms
c
      nect = 0
      ect = 0.0d0
      do i = 1, n
         aect(i) = 0.0d0
      end do
      if (nct .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHGTRN'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nct.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge Transfer Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,chgct,dmpct,n12,i12,n13,i13,
!$OMP& n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,nelst,
!$OMP& elst,use,use_group,use_intra,use_bounds,ctrntyp,f,off2,
!$OMP& elambda,mut,cut2,molcule,c0,c1,c2,c3,c4,c5,name,verbose,
!$OMP& debug,header,iout)
!$OMP& firstprivate(mscale) shared(ect,nect,aect,einter)
!$OMP DO reduction(+:ect,nect,aect,einter) schedule(guided)
c
c     compute the charge transfer energy
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         chgi = chgct(ii)
         alphai = dmpct(ii)
         if (alphai .eq. 0.0d0)  alphai = 1000.0d0
         usei = use(i)
         muti = mut(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = m2scale
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = m3scale
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = m4scale
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            mutk = mut(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  chgk = chgct(kk)
                  alphak = dmpct(kk)
                  if (alphak .eq. 0.0d0)  alphak = 1000.0d0
                  if (ctrntyp .eq. 'SEPARATE') then
                     expi = exp(-alphai*r)
                     expk = exp(-alphak*r)
                     e = -chgi*expk - chgk*expi
                  else
                     chgik = sqrt(abs(chgi*chgk))
                     expik = exp(-0.5d0*(alphai+alphak)*r)
                     e = -chgik * expik
                  end if
                  e = f * e * mscale(k)
c
c     apply lambda scaling for interaction annihilation 
c
                  if (muti .or. mutk) then
                     e = e * elambda
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group)  e = e * fgrp
c     
c     increment the overall charge transfer energy components
c
                  if (e .ne. 0.0d0) then
                     nect = nect + 1
                     ect = ect + e
                     aect(i) = aect(i) + 0.5d0*e
                     aect(k) = aect(k) + 0.5d0*e
                  end if
c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (e .gt. 10.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Charge Transfer',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,e
   30                format (' ChgTrn',4x,2(i7,'-',a3),
     &                          9x,f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      return
      end

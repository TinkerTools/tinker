c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge3  --  charge-charge energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge3" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms
c
c
      subroutine echarge3
      use limits
      use warp
      implicit none
c
c
c     choose the method for summing over pairwise interactions
c
      if (use_smooth) then
         call echarge3g
      else if (use_ewald) then
         if (use_clist) then
            call echarge3f
         else if (use_lights) then
            call echarge3e
         else
            call echarge3d
         end if
      else if (use_clist) then
         call echarge3c
      else if (use_lights) then
         call echarge3b
      else
         call echarge3a
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine echarge3a  --  charge analysis via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "echarge3a" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a pairwise
c     double loop
c
c
      subroutine echarge3a
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,in,ic,im
      integer kk,kn,kc,km
      real*8 e,fgrp
      real*8 r,r2,rb
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      logical header,huge
      character*6 mode
c
c
c     zero out the charge interaction energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHARGE'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nion.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge-Charge Interactions :',
     &           //,' Type',14x,'Atom Names',17x,'Charges',
     &              5x,'Distance',6x,'Energy',/)
      end if
c
c     compute and partition the charge interaction energy
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         im = molcule(i)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         usei = (use(i) .or. use(ic))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
            km = molcule(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
            if (proceed)  proceed = (cscale(kn) .ne. 0.0d0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               if (use_bounds)  call image (xc,yc,zc)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kk) * cscale(kn)
                  e = fik / rb
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e*taper + trans
                  end if
c
c     increment the overall charge-charge energy components
c
                  if (e .ne. 0.0d0) then
                     if (use_group)  e = e * fgrp
                     nec = nec + 1
                     ec = ec + e
                     aec(i) = aec(i) + 0.5d0*e
                     aec(k) = aec(k) + 0.5d0*e
                     if (im .ne. km)  einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Charge-Charge',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             17x,'Charges',5x,'Distance',
     &                             6x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),
     &                                pchg(ii),pchg(kk),r,e
   30                format (' Charge',4x,2(i7,'-',a3),8x,
     &                          2f7.2,f11.4,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         usei = (use(i) .or. use(ic))
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nion
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 2, ncell
                  xc = xic - x(kc)
                  yc = yic - y(kc)
                  zc = zic - z(kc)
                  call imager (xc,yc,zc,j)
                  rc2 = xc*xc + yc*yc + zc*zc
                  if (rc2 .le. off2) then
                     xr = xc + xi - x(k) + x(kc)
                     yr = yc + yi - y(k) + y(kc)
                     zr = zc + zi - z(k) + z(kc)
                     r2 = xr*xr + yr*yr + zr*zr
                     r = sqrt(r2)
                     rb = r + ebuffer
                     fik = fi * pchg(kk)
                     if (use_polymer) then
                        if (r2 .le. polycut2)  fik = fik * cscale(kn)
                     end if
                     e = fik / rb
c
c     use shifted energy switching if near the cutoff distance
c
                     shift = fik / (0.5d0*(off+cut))
                     e = e - shift
                     if (rc2 .gt. cut2) then
                        rc = sqrt(rc2)
                        rc3 = rc2 * rc
                        rc4 = rc2 * rc2
                        rc5 = rc2 * rc3
                        rc6 = rc3 * rc3
                        rc7 = rc3 * rc4
                        taper = c5*rc5 + c4*rc4 + c3*rc3
     &                             + c2*rc2 + c1*rc + c0
                        trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                                  + f3*rc3 + f2*rc2 + f1*rc + f0)
                        e = e*taper + trans
                     end if
c
c     increment the overall charge-charge energy components
c
                     if (e .ne. 0.0d0) then
                        if (use_group)  e = e * fgrp
                        if (i .eq. k)  e = 0.5d0 * e
                        nec = nec + 1
                        ec = ec + e
                        aec(i) = aec(i) + 0.5d0*e
                        aec(k) = aec(k) + 0.5d0*e
                        einter = einter + e
                     end if
c
c     print a message if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 100.0d0)
                     if ((debug.and.e.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Charge-Charge',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                17x,'Charges',5x,'Distance',
     &                                6x,'Energy',/)
                        end if
                        write (iout,50)  i,name(i),k,name(k),
     &                                   pchg(ii),pchg(kk),r,e
   50                   format (' Charge',4x,2(i7,'-',a3),1x,
     &                             '(XTAL)',1x,2f7.2,f11.4,f12.4)
                     end if
                  end if
               end do
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine echarge3b  --  method of lights charge analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echarge3b" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using the method
c     of lights
c
c
      subroutine echarge3b
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use light
      use molcul
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,in,ic,im
      integer kk,kn,kc,km
      integer kgy,kgz,kmap
      integer start,stop
      integer ikmin,ikmax
      real*8 e,fgrp
      real*8 r,r2,rb
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8, allocatable :: cscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical unique,repeat
      logical header,huge
      character*6 mode
c
c
c     zero out the charge interaction energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
      allocate (xsort(8*n))
      allocate (ysort(8*n))
      allocate (zsort(8*n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHARGE'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nion.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge-Charge Interactions :',
     &           //,' Type',14x,'Atom Names',17x,'Charges',
     &              5x,'Distance',6x,'Energy',/)
      end if
c
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nion
         k = kion(i)
         xsort(i) = x(k)
         ysort(i) = y(k)
         zsort(i) = z(k)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .true.
      call lights (off,nion,xsort,ysort,zsort,unique)
c
c     loop over all atoms computing the interactions
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         im = molcule(i)
         usei = (use(i) .or. use(ic))
         xic = xsort(rgx(ii))
         yic = ysort(rgy(ii))
         zic = zsort(rgz(ii))
         xi = x(i) - x(ic)
         yi = y(i) - y(ic)
         zi = z(i) - z(ic)
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
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
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 60
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 60
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 60
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 60
            end if
            kmap = kk - ((kk-1)/nion)*nion
            k = iion(kmap)
            kn = jion(kmap)
            kc = kion(kmap)
            km = molcule(k)
            prime = (kk .le. nion)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - xsort(j)
               yc = yic - ysort(kgy)
               zc = zic - zsort(kgz)
               if (use_bounds) then
                  if (abs(xc) .gt. xcell2)  xc = xc - sign(xcell,xc)
                  if (abs(yc) .gt. ycell2)  yc = yc - sign(ycell,yc)
                  if (abs(zc) .gt. zcell2)  zc = zc - sign(zcell,zc)
                  if (monoclinic) then
                     xc = xc + zc*beta_cos
                     zc = zc * beta_sin
                  else if (triclinic) then
                     xc = xc + yc*gamma_cos + zc*beta_cos
                     yc = yc*gamma_sin + zc*beta_term
                     zc = zc * gamma_term
                  end if
               end if
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kmap)
                  if (prime)  fik = fik * cscale(kn)
                  if (use_polymer) then
                     if (r2 .gt. polycut2)  fik = fi * pchg(kmap)
                  end if
                  e = fik / rb
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e*taper + trans
                  end if
c
c     increment the overall charge-charge energy components
c
                  if (e .ne. 0.0d0) then
                     if (use_group)  e = e * fgrp
                     nec = nec + 1
                     ec = ec + e
                     aec(i) = aec(i) + 0.5d0*e
                     aec(k) = aec(k) + 0.5d0*e
                     if (.not.prime .or. im.ne.km)  einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,30)
   30                   format (/,' Individual Charge-Charge',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             17x,'Charges',5x,'Distance',
     &                             6x,'Energy',/)
                     end if
                     ikmin = min(i,k)
                     ikmax = max(i,k)
                     if (prime) then
                        write (iout,40)  ikmin,name(ikmin),ikmax,
     &                                   name(ikmax),pchg(ii),
     &                                   pchg(kmap),r,e
   40                   format (' Charge',4x,2(i7,'-',a3),8x,
     &                             2f7.2,f11.4,f12.4)
                     else
                        write (iout,50)  ikmin,name(ikmin),ikmax,
     &                                   name(ikmax),pchg(ii),
     &                                   pchg(kmap),r,e
   50                   format (' Charge',4x,2(i7,'-',a3),1x,
     &                             '(XTAL)',1x,2f7.2,f11.4,f12.4)
                     end if
                  end if
               end if
            end if
   60       continue
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
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge3c  --  neighbor list charge analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge3c" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a pairwise
c     neighbor list
c
c
      subroutine echarge3c
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use charge
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use neigh
      use shunt
      use usage
      implicit none
      integer i,j,k,kkk
      integer ii,in,ic,im
      integer kk,kn,kc,km
      real*8 e,fgrp
      real*8 r,r2,rb
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 shift,taper,trans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      logical header,huge
      character*6 mode
c
c
c     zero out the charge interaction energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHARGE'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nion.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge-Charge Interactions :',
     &           //,' Type',14x,'Atom Names',17x,'Charges',
     &              5x,'Distance',6x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nion,iion,jion,kion,use,
!$OMP& x,y,z,f,pchg,nelst,elst,n12,n13,n14,n15,i12,i13,i14,i15,
!$OMP& c1scale,c2scale,c3scale,c4scale,c5scale,use_group,use_bounds,
!$OMP& off,off2,cut,cut2,c0,c1,c2,c3,c4,c5,f0,f1,f2,f3,f4,f5,f6,f7,
!$OMP% molcule,ebuffer,name,verbose,debug,header,iout)
!$OMP& firstprivate(cscale) shared (ec,nec,aec,einter)
!$OMP DO reduction(+:ec,nec,aec,einter) schedule(guided)
c
c     compute and partition the charge interaction energy
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         ic = kion(ii)
         im = molcule(i)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(ii)
         usei = (use(i) .or. use(ic))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
            km = molcule(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
            if (proceed)  proceed = (cscale(kn) .ne. 0.0d0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               if (use_bounds)  call image (xc,yc,zc)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kk) * cscale(kn)
                  e = fik / rb
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     e = e*taper + trans
                  end if
c
c     increment the overall charge-charge energy components
c
                  if (e .ne. 0.0d0) then
                     if (use_group)  e = e * fgrp
                     nec = nec + 1
                     ec = ec + e
                     aec(i) = aec(i) + 0.5d0*e
                     aec(k) = aec(k) + 0.5d0*e
                     if (im .ne. km)  einter = einter + e
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Charge-Charge',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             17x,'Charges',5x,'Distance',
     &                             6x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),
     &                                pchg(ii),pchg(kk),r,e
   30                format (' Charge',4x,2(i7,'-',a3),8x,
     &                          2f7.2,f11.4,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
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
      deallocate (cscale)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge3d  --  Ewald charge analysis via loop  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge3d" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a particle
c     mesh Ewald summation
c
c
      subroutine echarge3d
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use energi
      use ewald
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use pme
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,in,im
      integer kk,kn,km
      real*8 e,efull
      real*8 f,fi,fik
      real*8 fs,fgrp
      real*8 r,r2,rb,rew
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      logical header,huge
      character*6 mode
      external erfc
c
c
c     zero out the Ewald summation energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bseorder
      aewald = aeewald
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nion.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge-Charge Interactions :',
     &           //,' Type',14x,'Atom Names',17x,'Charges',
     &              5x,'Distance',6x,'Energy',/)
      end if
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / rootpi
      do ii = 1, nion
         i = iion(ii)
         e = fs * pchg(ii)**2
         ec = ec + e
         nec = nec + 1
         aec(i) = aec(i) + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, nion
            i = iion(ii)
            xd = xd + pchg(ii)*x(i)
            yd = yd + pchg(ii)*y(i)
            zd = zd + pchg(ii)*z(i)
         end do
         e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
         nec = nec + 1
         do ii = 1, nion
            i = iion(ii)
            aec(i) = aec(i) + e/dble(nion)
         end do
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      call ecrecip3
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         im = molcule(i)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            km = molcule(k)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
c
c     find energy for interactions within real space cutoff
c
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kk)
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  e = (fik/rb) * (erfterm+scale-1.0d0)
c
c     increment the overall charge-charge energy components
c
                  ec = ec + e
                  efull = (fik/rb) * scale
                  if (efull .ne. 0.0d0) then
                     nec = nec + 1
                     aec(i) = aec(i) + 0.5d0*efull
                     aec(k) = aec(k) + 0.5d0*efull
                     if (im .ne. km)  einter = einter + efull
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(efull) .gt. 100.0d0)
                  if ((debug.and.efull.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Real Space Ewald',
     &                             ' Charge-Charge Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             17x,'Charges',5x,'Distance',
     &                             6x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),
     &                                pchg(ii),pchg(kk),r,efull
   30                format (' Charge',4x,2(i7,'-',a3),8x,
     &                          2f7.2,f11.4,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate real space portion involving other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nion
            k = iion(kk)
            kn = jion(kk)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 2, ncell
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  call imager (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     rb = r + ebuffer
                     fik = fi * pchg(kk)
                     rew = aewald * r
                     erfterm = erfc (rew)
                     scale = 1.0d0
                     if (use_group)  scale = scale * fgrp
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           scale = scale * cscale(kn)
                        end if
                     end if
                     e = (fik/rb) * (erfterm+scale-1.0d0)
                     if (i .eq. k)  e = 0.5d0 * e
c
c     increment the overall charge-charge energy components
c
                     ec = ec + e
                     aec(i) = aec(i) + 0.5d0*e
                     aec(k) = aec(k) + 0.5d0*e
                     efull = (fik/rb) * scale
                     if (i .eq. k)  efull = 0.5d0 * efull
                     if (efull .ne. 0.0d0) then
                        nec = nec + 1
                        einter = einter + efull
                     end if
c
c     print a message if the energy of this interaction is large
c
                     huge = (abs(efull) .gt. 100.0d0)
                     if ((debug.and.efull.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Real Space Ewald',
     &                                ' Charge-Charge Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                17x,'Charges',5x,'Distance',
     &                                6x,'Energy',/)
                        end if
                        write (iout,50)  i,name(i),k,name(k),
     &                                   pchg(ii),pchg(kk),r,efull
   50                   format (' Charge',4x,2(i7,'-',a3),1x,
     &                             '(XTAL)',1x,2f7.2,f11.4,f12.4)
                     end if
                  end if
               end do
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine echarge3e  --  Ewald charge analysis via lights  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echarge3e" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a particle
c     mesh Ewald summation and the method of lights
c
c
      subroutine echarge3e
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use energi
      use ewald
      use group
      use inform
      use inter
      use iounit
      use light
      use math
      use molcul
      use pme
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,in,im
      integer kk,kn,km
      integer kgy,kgz,kmap
      integer start,stop
      real*8 e,efull
      real*8 f,fi,fik
      real*8 fs,fgrp
      real*8 r,r2,rb,rew
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale
      real*8, allocatable :: cscale(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical unique,repeat
      logical header,huge
      character*6 mode
      external erfc
c
c
c     zero out the Ewald summation energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bseorder
      aewald = aeewald
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
      allocate (xsort(8*n))
      allocate (ysort(8*n))
      allocate (zsort(8*n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nion.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge-Charge Interactions :',
     &           //,' Type',14x,'Atom Names',17x,'Charges',
     &              5x,'Distance',6x,'Energy',/)
      end if
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / rootpi
      do ii = 1, nion
         i = iion(ii)
         e = fs * pchg(ii)**2
         ec = ec + e
         nec = nec + 1
         aec(i) = aec(i) + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, nion
            i = iion(ii)
            xd = xd + pchg(ii)*x(i)
            yd = yd + pchg(ii)*y(i)
            zd = zd + pchg(ii)*z(i)
         end do
         e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
         nec = nec + 1
         do ii = 1, nion
            i = iion(ii)
            aec(i) = aec(i) + e/dble(nion)
         end do
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      call ecrecip3
c
c     compute the real space portion of the Ewald summation;
c     transfer the interaction site coordinates to sorting arrays
c
      do i = 1, nion
         k = iion(i)
         xsort(i) = x(k)
         ysort(i) = y(k)
         zsort(i) = z(k)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .true.
      call lights (off,nion,xsort,ysort,zsort,unique)
c
c     loop over all atoms computing the interactions
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         im = molcule(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         fi = f * pchg(ii)
         usei = use(i)
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
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
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 60
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 60
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 60
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 60
            end if
            kmap = kk - ((kk-1)/nion)*nion
            k = iion(kmap)
            kn = jion(kmap)
            km = molcule(k)
            prime = (kk .le. nion)
c
c     decide whether to compute the current interaction
c
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
c
c     find energy for interactions within real space cutoff
c
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
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = 1.0d0
                  if (prime)  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  fik = fi * pchg(kmap)
                  if (use_polymer) then
                     if (r2 .gt. polycut2)  fik = fi * pchg(kmap)
                  end if
                  e = (fik/rb) * (erfterm+scale-1.0d0)
c
c     increment the overall charge-charge energy components
c
                  ec = ec + e
                  aec(i) = aec(i) + 0.5d0*e
                  aec(k) = aec(k) + 0.5d0*e
                  efull = (fik/rb) * scale
                  if (efull .ne. 0.0d0) then
                     nec = nec + 1
                     if (.not.prime .or. im.ne.km)
     &                  einter = einter + efull
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(efull) .gt. 100.0d0)
                  if ((debug.and.efull.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,30)
   30                   format (/,' Individual Real Space Ewald',
     &                             ' Charge-Charge Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             17x,'Charges',5x,'Distance',
     &                             6x,'Energy',/)
                     end if
                     if (prime) then
                        write (iout,40)  i,name(i),k,name(k),pchg(ii),
     &                                   pchg(kmap),r,efull
   40                   format (' Charge',4x,2(i7,'-',a3),8x,
     &                             2f7.2,f11.4,f12.4)
                     else
                        write (iout,50)  i,name(i),k,name(k),pchg(ii),
     &                                   pchg(kmap),r,efull
   50                   format (' Charge',4x,2(i7,'-',a3),1x,
     &                             '(XTAL)',1x,2f7.2,f11.4,f12.4)
                     end if
                  end if
               end if
            end if
   60       continue
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
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge3f  --  Ewald charge analysis via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge3f" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms using a particle
c     mesh Ewald summation and a pairwise neighbor list
c
c
      subroutine echarge3f
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use energi
      use ewald
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use neigh
      use pme
      use shunt
      use usage
      implicit none
      integer i,j,k,kkk
      integer ii,in,im
      integer kk,kn,km
      real*8 e,efull
      real*8 f,fi,fik
      real*8 fs,fgrp
      real*8 r,r2,rb,rew
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xd,yd,zd
      real*8 erfc,erfterm
      real*8 scale
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      logical header,huge
      character*6 mode
      external erfc
c
c
c     zero out the Ewald summation energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bseorder
      aewald = aeewald
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nion.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge-Charge Interactions :',
     &           //,' Type',14x,'Atom Names',17x,'Charges',
     &              5x,'Distance',6x,'Energy',/)
      end if
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / rootpi
      do ii = 1, nion
         i = iion(ii)
         e = fs * pchg(ii)**2
         ec = ec + e
         nec = nec + 1
         aec(i) = aec(i) + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, nion
            i = iion(ii)
            xd = xd + pchg(ii)*x(i)
            yd = yd + pchg(ii)*y(i)
            zd = zd + pchg(ii)*z(i)
         end do
         e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
         nec = nec + 1
         do ii = 1, nion
            i = iion(ii)
            aec(i) = aec(i) + e/dble(nion)
         end do
      end if
c
c     compute the reciprocal space part of the Ewald summation
c
      call ecrecip3
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nion,iion,jion,use,
!$OMP& x,y,z,f,pchg,nelst,elst,n12,n13,n14,n15,i12,i13,i14,i15,
!$OMP& c1scale,c2scale,c3scale,c4scale,c5scale,use_group,off2,
!$OMP& aewald,molcule,ebuffer,name,verbose,debug,header,iout)
!$OMP& firstprivate(cscale) shared (ec,nec,aec,einter)
!$OMP DO reduction(+:ec,nec,aec,einter) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(ii)
         im = molcule(i)
         usei = use(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = iion(kk)
            kn = jion(kk)
            km = molcule(k)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
c
c     find energy for interactions within real space cutoff
c
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  fik = fi * pchg(kk)
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  e = (fik/rb) * (erfterm+scale-1.0d0)
c
c     increment the overall charge-charge energy components
c
                  ec = ec + e
                  aec(i) = aec(i) + 0.5d0*e
                  aec(k) = aec(k) + 0.5d0*e
                  efull = (fik/rb) * scale
                  if (efull .ne. 0.0d0) then
                     nec = nec + 1
                     if (im .ne. km)  einter = einter + efull
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (abs(efull) .gt. 100.0d0)
                  if ((debug.and.efull.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Real Space Ewald',
     &                             ' Charge-Charge Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             17x,'Charges',5x,'Distance',
     &                             6x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),
     &                                pchg(ii),pchg(kk),r,efull
   30                format (' Charge',4x,2(i7,'-',a3),8x,
     &                          2f7.2,f11.4,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
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
      deallocate (cscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge3g  --  charge analysis for smoothing  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge3g" calculates the charge-charge interaction energy
c     and partitions the energy among the atoms for use with
c     potential smoothing methods
c
c
      subroutine echarge3g
      use action
      use analyz
      use atomid
      use atoms
      use charge
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use molcul
      use usage
      use warp
      implicit none
      integer i,j,k
      integer ii,in,im
      integer kk,kn,km
      real*8 e,fgrp
      real*8 r,r2,rb,rb2
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 erf,wterm,width
      real*8 width2,width3
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      logical header,huge
      external erf
c
c
c     zero out the charge interaction energy and partitioning
c
      nec = 0
      ec = 0.0d0
      do i = 1, n
         aec(i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. nion.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Charge-Charge Interactions :',
     &           //,' Type',14x,'Atom Names',17x,'Charges',
     &              5x,'Distance',6x,'Energy',/)
      end if
c
c     set the energy units conversion factor
c
      f = electric / dielec
c
c     set the extent of smoothing to be performed
c
      width = deform * diffc
      if (use_dem) then
         if (width .gt. 0.0d0)  width = 0.5d0 / sqrt(width)
      else if (use_gda) then
         wterm = sqrt(3.0d0/(2.0d0*diffc))
      end if
      width2 = width * width
      width3 = width * width2
c
c     compute and partition the charge interaction energy
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(ii)
         im = molcule(i)
         usei = (use(i))
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(ii)
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(kk)
            km = molcule(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed)  proceed = (cscale(kn) .ne. 0.0d0)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
               rb = r + ebuffer
               fik = fi * pchg(kk) * cscale(kn)
               e = fik / rb
c
c     transform the potential function via smoothing
c
               if (use_dem) then
                  if (width .gt. 0.0d0) then
                     e = e * erf(width*rb)
                  end if
               else if (use_gda) then
                  width = m2(i) + m2(k)
                  if (width .gt. 0.0d0) then
                     width = wterm / sqrt(width)
                     e = e * erf(width*rb)
                  end if
               else if (use_tophat) then
                  if (width .gt. rb) then
                     rb2 = rb * rb
                     e = fik * (3.0d0*width2-rb2) / (2.0d0*width3)
                  end if
               else if (use_stophat) then
                  e = fik / (rb+width)
               end if
c
c     increment the overall charge-charge energy components
c
               if (e .ne. 0.0d0) then
                  if (use_group)  e = e * fgrp
                  nec = nec + 1
                  ec = ec + e
                  aec(i) = aec(i) + 0.5d0*e
                  aec(k) = aec(k) + 0.5d0*e
                  if (im .ne. km)  einter = einter + e
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 100.0d0)
               if ((debug.and.e.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Charge-Charge',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          17x,'Charges',5x,'Distance',
     &                          6x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),k,name(k),pchg(ii),
     &                             pchg(kk),r,e
   30             format (' Charge',4x,2(i7,'-',a3),8x,
     &                       2f7.2,f11.4,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ecrecip3  --  PME reciprocal charge analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ecrecip3" evaluates the reciprocal space portion of the
c     particle mesh Ewald energy due to partial charges, and
c     partitions the energy among the atoms
c
c     literature reference:
c
c     U. Essmann, L. Perera, M. L Berkowitz, T. Darden, H. Lee and
c     L. G. Pedersen, "A Smooth Particle Mesh Ewald Method", Journal
c     of Chemical Physics, 103, 8577-8593 (1995)
c
c     W. Smith and D. Fincham, "The Ewald Sum in Truncated Octahedral
c     and Rhombic Dodecahedral Boundary Conditions", Molecular
c     Simulation, 10, 67-71 (1993)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine ecrecip3
      use analyz
      use bound
      use boxes
      use charge
      use chgpot
      use energi
      use ewald
      use math
      use pme
      implicit none
      integer i,j,k
      integer k1,k2,k3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff,ntot
      real*8 e,f,hsq,denom
      real*8 term,expterm
      real*8 pterm,volterm
      real*8 h1,h2,h3
      real*8 r1,r2,r3
      real*8, allocatable :: fphi(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = 0.5d0 * electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      ntot = nfft1 * nfft2 * nfft3
      if (allocated(qgrid)) then
         if (size(qgrid) .ne. 2*ntot)  call fftclose
      end if
      if (.not. allocated(qgrid))  call fftsetup
c
c     setup spatial decomposition and B-spline coefficients
c
      call getchunk
      call moduli
      call bspline_fill
      call table_fill
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_pchg
      call fftfront
c
c     use scalar sum to get the reciprocal space energy
c
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      nff = nfft1 * nfft2
      ntot = nff * nfft3
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (nonprism) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
         end if
         qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
         qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
      end do
c
c     account for zeroth grid point for nonperiodic system
c
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qgrid(1,1,1,1) = expterm * qgrid(1,1,1,1)
         qgrid(2,1,1,1) = expterm * qgrid(2,1,1,1)
      end if
c
c     perform the 3-D FFT backward transformation
c
      call fftback
c
c     perform dynamic allocation of some local arrays
c
      allocate (fphi(4,nion))
c
c     extract the partial charge electrostatic potential
c
      call fphi_pchg (fphi)
c
c     sum over charges and increment total charge energy
c
      e = 0.0d0
      do i = 1, nion
         j = iion(i)
         term = f * pchg(i) * fphi(1,i)
         e = e + term
         aec(j) = aec(j) + term
      end do
      ec = ec + e
c
c     perform deallocation of some local arrays
c
      deallocate (fphi)
      return
      end

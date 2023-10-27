c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2006 by David Gohara & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine nblist  --  maintain pairwise neighbor lists  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nblist" builds and maintains nonbonded pair neighbor lists
c     for vdw, dispersion, electrostatic and polarization terms
c
c
      subroutine nblist
      use limits
      use neigh
      use potent
      implicit none
c
c
c     rebuild list if using both vdw and dispersion terms
c
      if (use_vdw .and. use_disp) then
         dovlst = .true.
         dodlst = .true.
      end if
c
c     rebuild list if using both charge and multipole terms
c
      if (use_charge .and. use_mpole) then
         doclst = .true.
         domlst = .true.
      end if
c
c     update the appropriate nonbonded neighbor lists
c
      if (use_vdw .and. use_vlist)  call vlist
      if (use_disp .and. use_dlist)  call dlist
      if ((use_charge.or.use_solv) .and. use_clist)  call clist
      if ((use_repel.or.use_xrepel.or.use_mpole.or.use_polar
     &       .or.use_chgtrn.or.use_solv) .and. use_mlist)  call mlist
      if (use_polar .and. use_ulist)  call ulist
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine vlist  --  get van der Waals neighbor lists  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "vlist" performs an update or a complete rebuild of the
c     nonbonded neighbor lists for vdw sites
c
c
      subroutine vlist
      use atoms
      use bound
      use boxes
      use iounit
      use neigh
      use vdw
      implicit none
      integer i,j,k
      integer ii,kk,iv
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius
      real*8 rdn,r2
      logical, allocatable :: update(:)
c
c
c     apply reduction factors to find coordinates for each site
c
      do ii = 1, nvdw
         i = ivdw(ii)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(vbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' VLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (dovlst) then
         dovlst = .false.
         if (nonprism) then
            call vbuild
         else
            call vlight
         end if
         return
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(n))
c
c     test sites for displacement exceeding half the buffer
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do ii = 1, nvdw
         i = ivdw(ii)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         xr = xi - xvold(i)
         yr = yi - yvold(i)
         zr = zi - zvold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            update(i) = .true.
            xvold(i) = xi
            yvold(i) = yi
            zvold(i) = zi
         end if
      end do
!$OMP END DO
c
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do ii = 1, nvdw
         i = ivdw(ii)
         if (update(i)) then
            xi = xvold(i)
            yi = yvold(i)
            zi = zvold(i)
            nvlst(i) = 0
            do kk = ii+1, nvdw
               k = ivdw(kk)
               xr = xi - xvold(k)
               yr = yi - yvold(k)
               zr = zi - zvold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. vbuf2) then
                  nvlst(i) = nvlst(i) + 1
                  vlst(nvlst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists for lower numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do ii = 1, nvdw
         i = ivdw(ii)
         if (update(i)) then
            xi = xvold(i)
            yi = yvold(i)
            zi = zvold(i)
            do kk = 1, ii-1
               k = ivdw(kk)
               if (.not. update(k)) then
                  xr = xi - xvold(k)
                  yr = yi - yvold(k)
                  zr = zi - zvold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. vbuf2) then
!$OMP CRITICAL
                     do j = 1, nvlst(k)
                        if (vlst(j,k) .eq. i)  goto 20
                     end do
                     nvlst(k) = nvlst(k) + 1
                     vlst(nvlst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. vbufx) then
!$OMP CRITICAL
                     do j = 1, nvlst(k)
                        if (vlst(j,k) .eq. i) then
                           vlst(j,k) = vlst(nvlst(k),k)
                           nvlst(k) = nvlst(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do ii = 1, nvdw
         i = ivdw(ii)
         if (nvlst(i) .ge. maxvlst) then
            write (iout,40)
   40       format (/,' VLIST  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (update)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine vbuild  --  build vdw list for all sites  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "vbuild" performs a complete rebuild of the van der Waals
c     pair neighbor list for all sites
c
c
      subroutine vbuild
      use bound
      use iounit
      use neigh
      use vdw
      implicit none
      integer i,k,ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store coordinates to reflect update of the site
c
      do ii = 1, nvdw
         i = ivdw(ii)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         xvold(i) = xi
         yvold(i) = yi
         zvold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nvlst(i) = 0
         do kk = ii+1, nvdw
            k = ivdw(kk)
            xr = xi - xred(k)
            yr = yi - yred(k)
            zr = zi - zred(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. vbuf2) then
               nvlst(i) = nvlst(i) + 1
               vlst(nvlst(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nvlst(i) .ge. maxvlst) then
            write (iout,10)
   10       format (/,' VBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine vlight  --  build vdw pair list via lights  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "vlight" performs a complete rebuild of the van der Waals
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine vlight
      use atoms
      use bound
      use cell
      use iounit
      use light
      use neigh
      use vdw
      implicit none
      integer i,j,k
      integer ii,kk
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical unique,repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(nvdw))
      allocate (ysort(nvdw))
      allocate (zsort(nvdw))
c
c     transfer interaction site coordinates to sorting arrays
c
      do ii = 1, nvdw
         i = ivdw(ii)
         nvlst(i) = 0
         xvold(i) = xred(i)
         yvold(i) = yred(i)
         zvold(i) = zred(i)
         xsort(ii) = xred(i)
         ysort(ii) = yred(i)
         zsort(ii) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .false.
      off = sqrt(vbuf2)
      call lights (off,nvdw,xsort,ysort,zsort,unique)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,
!$OMP& xi,yi,zi,xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do ii = 1, nvdw
         i = ivdw(ii)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii)
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            if (kk .le. ii)  goto 20
            k = ivdw(kk)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            xr = xi - xred(k)
            yr = yi - yred(k)
            zr = zi - zred(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. vbuf2) then
               nvlst(i) = nvlst(i) + 1
               vlst(nvlst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii)
            stop = nvdw
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nvlst(i) .ge. maxvlst) then
            write (iout,30)
   30       format (/,' VLIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine dlist  --  get damped dispersion neighbor lists  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "dlist" performs an update or a complete rebuild of the
c     nonbonded neighbor lists for damped dispersion sites
c
c
      subroutine dlist
      use atoms
      use bound
      use boxes
      use disp
      use iounit
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(n))
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(dbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' DLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (dodlst) then
         dodlst = .false.
         if (nonprism) then
            call dbuild
         else
            call dlight
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do ii = 1, ndisp
         i = idisp(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         xr = xi - xvold(i)
         yr = yi - yvold(i)
         zr = zi - zvold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            update(i) = .true.
            xvold(i) = xi
            yvold(i) = yi
            zvold(i) = zi
         end if
      end do
!$OMP END DO
c
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do ii = 1, ndisp
         i = idisp(ii)
         if (update(i)) then
            xi = xvold(i)
            yi = yvold(i)
            zi = zvold(i)
            nvlst(i) = 0
            do kk = ii+1, ndisp
               k = idisp(kk)
               xr = xi - xvold(k)
               yr = yi - yvold(k)
               zr = zi - zvold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. dbuf2) then
                  nvlst(i) = nvlst(i) + 1
                  vlst(nvlst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists for lower numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do ii = 1, ndisp
         i = idisp(ii)
         if (update(i)) then
            xi = xvold(i)
            yi = yvold(i)
            zi = zvold(i)
            do k = 1, i-1
               if (.not. update(k)) then
                  xr = xi - xvold(k)
                  yr = yi - yvold(k)
                  zr = zi - zvold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. dbuf2) then
!$OMP CRITICAL
                     do j = 1, nvlst(k)
                        if (vlst(j,k) .eq. i)  goto 20
                     end do
                     nvlst(k) = nvlst(k) + 1
                     vlst(nvlst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. dbufx) then
!$OMP CRITICAL
                     do j = 1, nvlst(k)
                        if (vlst(j,k) .eq. i) then
                           vlst(j,k) = vlst(nvlst(k),k)
                           nvlst(k) = nvlst(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do ii = 1, ndisp
         i = idisp(ii)
         if (nvlst(i) .ge. maxvlst) then
            write (iout,40)
   40       format (/,' DLIST  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (update)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine dbuild  --  build dispersion list for all sites  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "dbuild" performs a complete rebuild of the damped dispersion
c     neighbor list for all sites
c
c
      subroutine dbuild
      use atoms
      use bound
      use disp
      use iounit
      use neigh
      implicit none
      integer i,k,ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store new coordinates to reflect update of the site
c
      do ii = 1, ndisp
         i = idisp(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         xvold(i) = xi
         yvold(i) = yi
         zvold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nvlst(i) = 0
         do kk = ii+1, ndisp
            k = idisp(kk)
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. dbuf2) then
               nvlst(i) = nvlst(i) + 1
               vlst(nvlst(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nvlst(i) .ge. maxvlst) then
            write (iout,10)
   10       format (/,' DBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine dlight  --  get damp dispersion list via lights  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "dlight" performs a complete rebuild of the damped dispersion
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine dlight
      use atoms
      use bound
      use cell
      use disp
      use iounit
      use light
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical unique,repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(ndisp))
      allocate (ysort(ndisp))
      allocate (zsort(ndisp))
c
c     transfer interaction site coordinates to sorting arrays
c
      do ii = 1, ndisp
         i = idisp(ii)
         nvlst(i) = 0
         xvold(i) = x(i)
         yvold(i) = y(i)
         zvold(i) = z(i)
         xsort(ii) = x(i)
         ysort(ii) = y(i)
         zsort(ii) = z(i)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .false.
      off = sqrt(dbuf2)
      call lights (off,ndisp,xsort,ysort,zsort,unique)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,
!$OMP& xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do ii = 1, ndisp
         i = idisp(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii)
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            if (kk .le. ii)  goto 20
            k = idisp(kk)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. dbuf2) then
               nvlst(i) = nvlst(i) + 1
               vlst(nvlst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii)
            stop = ndisp
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nvlst(i) .ge. maxvlst) then
            write (iout,30)
   30       format (/,' DLIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine clist  --  get partial charge neighbor lists  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "clist" performs an update or a complete rebuild of the
c     nonbonded neighbor lists for partial charges
c
c
      subroutine clist
      use atoms
      use bound
      use boxes
      use charge
      use iounit
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk,ic
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(n))
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(cbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' CLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (doclst) then
         doclst = .false.
         if (nonprism) then
            call cbuild
         else
            call clight
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,ic,
!$OMP& xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do ii = 1, nion
         i = iion(ii)
         ic = kion(i)
         xi = x(ic)
         yi = y(ic)
         zi = z(ic)
         xr = xi - xeold(i)
         yr = yi - yeold(i)
         zr = zi - zeold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            update(i) = .true.
            xeold(i) = xi
            yeold(i) = yi
            zeold(i) = zi
         end if
      end do
!$OMP END DO
c
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do ii = 1, nion
         i = iion(ii)
         if (update(i)) then
            xi = xeold(i)
            yi = yeold(i)
            zi = zeold(i)
            nelst(i) = 0
            do kk = ii+1, nion
               k = iion(kk)
               xr = xi - xeold(k)
               yr = yi - yeold(k)
               zr = zi - zeold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. cbuf2) then
                  nelst(i) = nelst(i) + 1
                  elst(nelst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists for lower numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do ii = 1, nion
         i = iion(ii)
         if (update(i)) then
            xi = xeold(i)
            yi = yeold(i)
            zi = zeold(i)
            do kk = 1, ii-1
               k = iion(kk)
               if (.not. update(k)) then
                  xr = xi - xeold(k)
                  yr = yi - yeold(k)
                  zr = zi - zeold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. cbuf2) then
!$OMP CRITICAL
                     do j = 1, nelst(k)
                        if (elst(j,k) .eq. i)  goto 20
                     end do
                     nelst(k) = nelst(k) + 1
                     elst(nelst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. cbufx) then
!$OMP CRITICAL
                     do j = 1, nelst(k)
                        if (elst(j,k) .eq. i) then
                           elst(j,k) = elst(nelst(k),k)
                           nelst(k) = nelst(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do ii = 1, nion
         i = iion(ii)
         if (nelst(i) .ge. maxelst) then
            write (iout,40)
   40       format (/,' CLIST  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (update)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine cbuild  --  build charge list for all sites  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "cbuild" performs a complete rebuild of the partial charge
c     electrostatic neighbor list for all sites
c
c
      subroutine cbuild
      use atoms
      use bound
      use charge
      use iounit
      use neigh
      implicit none
      integer i,k
      integer ii,kk
      integer ic,kc
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,ii,kk,ic,kc,
!$OMP& xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store new coordinates to reflect update of the site
c
      do ii = 1, nion
         i = iion(ii)
         ic = kion(i)
         xi = x(ic)
         yi = y(ic)
         zi = z(ic)
         xeold(i) = xi
         yeold(i) = yi
         zeold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nelst(i) = 0
         do kk = ii+1, nion
            k = iion(kk)
            kc = kion(k)
            xr = xi - x(kc)
            yr = yi - y(kc)
            zr = zi - z(kc)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. cbuf2) then
               nelst(i) = nelst(i) + 1
               elst(nelst(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nelst(i) .ge. maxelst) then
            write (iout,10)
   10       format (/,' CBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine clight  --  get partial charge list via lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "clight" performs a complete rebuild of the partial charge
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine clight
      use atoms
      use bound
      use cell
      use charge
      use iounit
      use light
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      integer ic,kc
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical unique,repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(nion))
      allocate (ysort(nion))
      allocate (zsort(nion))
c
c     transfer interaction site coordinates to sorting arrays
c
      do ii = 1, nion
         i = iion(ii)
         ic = kion(i)
         nelst(i) = 0
         xeold(i) = x(ic)
         yeold(i) = y(ic)
         zeold(i) = z(ic)
         xsort(ii) = x(ic)
         ysort(ii) = y(ic)
         zsort(ii) = z(ic)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .false.
      off = sqrt(cbuf2)
      call lights (off,nion,xsort,ysort,zsort,unique)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,ic,kc,
!$OMP& xi,yi,zi,xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do ii = 1, nion
         i = iion(ii)
         ic = kion(i)
         xi = x(ic)
         yi = y(ic)
         zi = z(ic)
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii)
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            if (kk .le. ii)  goto 20
            k = iion(kk)
            kc = kion(k)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            xr = xi - x(kc)
            yr = yi - y(kc)
            zr = zi - z(kc)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. cbuf2) then
               nelst(i) = nelst(i) + 1
               elst(nelst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii)
            stop = nion
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nelst(i) .ge. maxelst) then
            write (iout,30)
   30       format (/,' CLIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mlist  --  get atomic multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mlist" performs an update or a complete rebuild of the
c     nonbonded neighbor lists for atomic multipoles
c
c
      subroutine mlist
      use atoms
      use bound
      use boxes
      use iounit
      use mpole
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(n))
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(mbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' MLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (domlst) then
         domlst = .false.
         if (nonprism) then
            call mbuild
         else
            call mlight
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         xr = xi - xeold(i)
         yr = yi - yeold(i)
         zr = zi - zeold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. lbuf2) then
            update(i) = .true.
            xeold(i) = xi
            yeold(i) = yi
            zeold(i) = zi
         end if
      end do
!$OMP END DO
c
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do ii = 1, npole
         i = ipole(ii)
         if (update(i)) then
            xi = xeold(i)
            yi = yeold(i)
            zi = zeold(i)
            nelst(i) = 0
            do kk = ii+1, npole
               k = ipole(kk)
               xr = xi - xeold(k)
               yr = yi - yeold(k)
               zr = zi - zeold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. mbuf2) then
                  nelst(i) = nelst(i) + 1
                  elst(nelst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists for lower numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do ii = 1, npole
         i = ipole(ii)
         if (update(i)) then
            xi = xeold(i)
            yi = yeold(i)
            zi = zeold(i)
            do kk = 1, ii-1
               k = ipole(kk)
               if (.not. update(k)) then
                  xr = xi - xeold(k)
                  yr = yi - yeold(k)
                  zr = zi - zeold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. mbuf2) then
!$OMP CRITICAL
                     do j = 1, nelst(k)
                        if (elst(j,k) .eq. i)  goto 20
                     end do
                     nelst(k) = nelst(k) + 1
                     elst(nelst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. mbufx) then
!$OMP CRITICAL
                     do j = 1, nelst(k)
                        if (elst(j,k) .eq. i) then
                           elst(j,k) = elst(nelst(k),k)
                           nelst(k) = nelst(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do ii = 1, npole
         i = ipole(ii)
         if (nelst(i) .ge. maxelst) then
            write (iout,40)
   40       format (/,' MLIST  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (update)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine mbuild  --  build mpole list for all sites  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "mbuild" performs a complete rebuild of the atomic multipole
c     electrostatic neighbor list for all sites
c
c
      subroutine mbuild
      use atoms
      use bound
      use iounit
      use mpole
      use neigh
      implicit none
      integer i,k,ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store new coordinates to reflect update of the site
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         xeold(i) = xi
         yeold(i) = yi
         zeold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nelst(i) = 0
         do kk = ii+1, npole
            k = ipole(kk)
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. mbuf2) then
               nelst(i) = nelst(i) + 1
               elst(nelst(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nelst(i) .ge. maxelst) then
            write (iout,10)
   10       format (/,' MBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mlight  --  get multipole pair list via lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mlight" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine mlight
      use atoms
      use bound
      use cell
      use iounit
      use light
      use mpole
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical unique,repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(npole))
      allocate (ysort(npole))
      allocate (zsort(npole))
c
c     transfer interaction site coordinates to sorting arrays
c
      do ii = 1, npole
         i = ipole(ii)
         nelst(i) = 0
         xeold(i) = x(i)
         yeold(i) = y(i)
         zeold(i) = z(i)
         xsort(ii) = x(i)
         ysort(ii) = y(i)
         zsort(ii) = z(i)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .false.
      off = sqrt(mbuf2)
      call lights (off,npole,xsort,ysort,zsort,unique)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,
!$OMP& xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii)
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            if (kk .le. ii)  goto 20
            k = ipole(kk)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. mbuf2) then
               nelst(i) = nelst(i) + 1
               elst(nelst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii)
            stop = npole
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nelst(i) .ge. maxelst) then
            write (iout,30)
   30       format (/,' MLIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ulist  --  get preconditioner neighbor lists  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ulist" performs an update or a complete rebuild of the
c     neighbor lists for the polarization preconditioner
c
c
      subroutine ulist
      use atoms
      use bound
      use boxes
      use iounit
      use mpole
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
      logical, allocatable :: update(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (update(n))
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(ubuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' ULIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (doulst) then
         doulst = .false.
         if (nonprism) then
            call ubuild
         else
            call ulight
         end if
         return
      end if
c
c     test sites for displacement exceeding half the buffer
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         xr = xi - xuold(i)
         yr = yi - yuold(i)
         zr = zi - zuold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         update(i) = .false.
         if (r2 .ge. pbuf2) then
            update(i) = .true.
            xuold(i) = xi
            yuold(i) = yi
            zuold(i) = zi
         end if
      end do
!$OMP END DO
c
c     rebuild the higher numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do ii = 1, npole
         i = ipole(ii)
         if (update(i)) then
            xi = xuold(i)
            yi = yuold(i)
            zi = zuold(i)
            nulst(i) = 0
            do kk = ii+1, npole
               k = ipole(kk)
               xr = xi - xuold(k)
               yr = yi - yuold(k)
               zr = zi - zuold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. ubuf2) then
                  nulst(i) = nulst(i) + 1
                  ulst(nulst(i),i) = k
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     adjust lists for lower numbered neighbors of updated sites
c
!$OMP DO schedule(guided)
      do ii = 1, npole
         i = ipole(ii)
         if (update(i)) then
            xi = xuold(i)
            yi = yuold(i)
            zi = zuold(i)
            do kk = 1, ii-1
               k = ipole(kk)
               if (.not. update(k)) then
                  xr = xi - xuold(k)
                  yr = yi - yuold(k)
                  zr = zi - zuold(k)
                  call imagen (xr,yr,zr)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. ubuf2) then
!$OMP CRITICAL
                     do j = 1, nulst(k)
                        if (ulst(j,k) .eq. i)  goto 20
                     end do
                     nulst(k) = nulst(k) + 1
                     ulst(nulst(k),k) = i
   20                continue
!$OMP END CRITICAL
                  else if (r2 .le. ubufx) then
!$OMP CRITICAL
                     do j = 1, nulst(k)
                        if (ulst(j,k) .eq. i) then
                           ulst(j,k) = ulst(nulst(k),k)
                           nulst(k) = nulst(k) - 1
                           goto 30
                        end if
                     end do
   30                continue
!$OMP END CRITICAL
                  end if
               end if
            end do
         end if
      end do
!$OMP END DO
c
c     check to see if any neighbor lists are too long
c
!$OMP DO schedule(guided)
      do ii = 1, npole
         i = ipole(ii)
         if (nulst(i) .ge. maxulst) then
            write (iout,40)
   40       format (/,' ULIST  --  Too many Neighbors;',
     &                 ' Increase MAXULST')
            call fatal
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (update)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ubuild  --  preconditioner list for all sites  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ubuild" performs a complete rebuild of the polarization
c     preconditioner neighbor list for all sites
c
c
      subroutine ubuild
      use atoms
      use bound
      use iounit
      use mpole
      use neigh
      implicit none
      integer i,k,ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,k,ii,kk,xi,yi,zi,xr,yr,zr,r2)
!$OMP DO schedule(guided)
c
c     store new coordinates to reflect update of the site
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         xuold(i) = xi
         yuold(i) = yi
         zuold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
         nulst(i) = 0
         do kk = ii+1, npole
            k = ipole(kk)
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. ubuf2) then
               nulst(i) = nulst(i) + 1
               ulst(nulst(i),i) = k
            end if
         end do
c
c     check to see if the neighbor list is too long
c
         if (nulst(i) .ge. maxulst) then
            write (iout,10)
   10       format (/,' UBUILD  --  Too many Neighbors;',
     &                 ' Increase MAXULST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ulight  --  get preconditioner list via lights  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ulight" performs a complete rebuild of the polarization
c     preconditioner pair neighbor list for all sites using the
c     method of lights
c
c
      subroutine ulight
      use atoms
      use bound
      use cell
      use iounit
      use light
      use mpole
      use neigh
      implicit none
      integer i,j,k
      integer ii,kk
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical unique,repeat
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xsort(npole))
      allocate (ysort(npole))
      allocate (zsort(npole))
c
c     transfer interaction site coordinates to sorting arrays
c
      do ii = 1, npole
         i = ipole(ii)
         nulst(i) = 0
         xuold(i) = x(i)
         yuold(i) = y(i)
         zuold(i) = z(i)
         xsort(ii) = x(i)
         ysort(ii) = y(i)
         zsort(ii) = z(i)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .false.
      off = sqrt(ubuf2)
      call lights (off,npole,xsort,ysort,zsort,unique)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,k,ii,kk,xi,yi,zi,
!$OMP& xr,yr,zr,r2,kgy,kgz,start,stop,repeat)
!$OMP DO schedule(guided)
c
c     loop over all atoms computing the neighbor lists
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii)
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            if (kk .le. ii)  goto 20
            k = ipole(kk)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            xr = xi - x(k)
            yr = yi - y(k)
            zr = zi - z(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. ubuf2) then
               nulst(i) = nulst(i) + 1
               ulst(nulst(i),i) = k
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii)
            stop = npole
            goto 10
         end if
c
c     check to see if the neighbor list is too long
c
         if (nulst(i) .ge. maxulst) then
            write (iout,30)
   30       format (/,' ULIGHT  --  Too many Neighbors;',
     &                 ' Increase MAXULST')
            call fatal
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end

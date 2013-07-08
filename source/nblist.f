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
c     "nblist" constructs and maintains nonbonded pair neighbor lists
c     for vdw and electrostatic interactions
c
c
      subroutine nblist
      implicit none
      include 'cutoff.i'
      include 'potent.i'
c
c
c     update the vdw and electrostatic neighbor lists
c
      if (use_vdw .and. use_vlist)  call vlist
      if (use_charge .and. use_clist)  call clist
      if ((use_mpole.or.use_polar) .and. use_mlist)  call mlist
      if (use_polar .and. use_ulist)  call ulist
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine vlist  --  build van der Waals neighbor lists  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "vlist" performs an update or a complete rebuild of the
c     van der Waals neighbor list
c
c
      subroutine vlist
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
c      include 'mpi.i'
      include 'neigh.i'
      include 'vdw.i'
      integer i,j,k
      integer ii,iv
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius
      real*8 rdn,r2
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
c      real*8 lbuffer_2
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
c
c     apply reduction factors to find coordinates for each site
c
      do i = 1, nvdw
         ii = ivdw(i)
         iv = ired(ii)
         rdn = kred(ii)
         xred(i) = rdn*(x(ii)-x(iv)) + x(iv)
         yred(i) = rdn*(y(ii)-y(iv)) + y(iv)
         zred(i) = rdn*(z(ii)-z(iv)) + z(iv)
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
c      lbuffer_2 = lbuffer/2
c
c     perform a complete list build instead of an update
c
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c      if (dovlst) then
c         dovlst = .false. 
c         call setup_nlist_builder(%val(nvdw), %val(nnode), 
c     &                            %val(radius), 
c     &                            %val(lbuffer_2), %val(xbox),
c     &                            %val(ybox), %val(zbox), %val(maxvlst))
c         call initial_nlist_build (xred,yred,zred,vlst,nvlst)
c         return
c      else
c         call subsequent_nlist_build (xred, yred, zred, vlst, nvlst)
c        return
c      end if
c
c
cccccccccccccccccccccccccccccccccccccccccccccccc
      if (dovlst) then
         dovlst = .false.
         if (octahedron) then
            do i = 1, nvdw
               call vbuild (i,xred,yred,zred)
            end do
         else
            call vlight (xred,yred,zred)
         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, nvdw
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         xr = xi - xvold(i)
         yr = yi - yvold(i)
         zr = zi - zvold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
c
c     store coordinates to reflect update of this site
c
         if (r2 .ge. lbuf2) then
            xvold(i) = xi
            yvold(i) = yi
            zvold(i) = zi
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, nvdw
               xr = xi - xvold(k)
               yr = yi - yvold(k)
               zr = zi - zvold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. vbuf2) then
                  j = j + 1
                  vlst(j,i) = k
               end if
            end do
            nvlst(i) = j
c
c     check to see if the neighbor list is too long
c
            if (nvlst(i) .ge. maxvlst) then
               write (iout,20)
   20          format (/,' VLIST  --  Too many Neighbors;',
     &                    ' Increase MAXVLST')
               call fatal
            end if
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               xr = xi - xvold(k)
               yr = yi - yvold(k)
               zr = zi - zvold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. vbuf2) then
                  do j = 1, nvlst(k)
                     if (vlst(j,k) .eq. i)  goto 30
                  end do
                  nvlst(k) = nvlst(k) + 1
                  vlst(nvlst(k),k) = i
   30             continue
               else if (r2 .le. vbufx) then
                  do j = 1, nvlst(k)
                     if (vlst(j,k) .eq. i) then
                        vlst(j,k) = vlst(nvlst(k),k)
                        nvlst(k) = nvlst(k) - 1
                        goto 40
                     end if
                  end do
   40             continue
               end if
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine vbuild  --  make vdw pair list for one site  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "vbuild" performs a complete rebuild of the van der Waals
c     pair neighbor list for a single site
c
c
      subroutine vbuild (i,xred,yred,zred)
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'iounit.i'
      include 'neigh.i'
      include 'vdw.i'
      integer i,j,k
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
      real*8 xred(*)
      real*8 yred(*)
      real*8 zred(*)
c
c
c     store coordinates to reflect update of the site
c
      xi = xred(i)
      yi = yred(i)
      zi = zred(i)
      xvold(i) = xi
      yvold(i) = yi
      zvold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
      j = 0
      do k = i+1, nvdw
         xr = xi - xred(k)
         yr = yi - yred(k)
         zr = zi - zred(k)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .le. vbuf2) then
            j = j + 1
            vlst(j,i) = k
         end if
      end do
      nvlst(i) = j
c
c     check to see if the neighbor list is too long
c
      if (nvlst(i) .ge. maxvlst) then
         write (iout,10)
   10    format (/,' VBUILD  --  Too many Neighbors;',
     &              ' Increase MAXVLST')
         call fatal
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine vlight  --  make vdw pair list for all sites  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "vlight" performs a complete rebuild of the van der Waals
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine vlight (xred,yred,zred)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'neigh.i'
      include 'vdw.i'
      integer i,j,k
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8 xred(*)
      real*8 yred(*)
      real*8 zred(*)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
      nlight = (ncell+1) * nvdw
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         nvlst(i) = 0
         xvold(i) = xred(i)
         yvold(i) = yred(i)
         zvold(i) = zred(i)
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      off = sqrt(vbuf2)
      call lights (off,nvdw,xsort,ysort,zsort)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     loop over all atoms computing the neighbor lists
c
      do i = 1, nvdw
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i) + 1
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - xred(k)
            yr = yi - yred(k)
            zr = zi - zred(k)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. vbuf2) then
               if (i .lt. k) then
                  nvlst(i) = nvlst(i) + 1
                  vlst(nvlst(i),i) = k
               else
                  nvlst(k) = nvlst(k) + 1
                  vlst(nvlst(k),k) = i
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i) + 1
            stop = nlight
            goto 10
         end if
      end do
c
c     check to see if the neighbor lists are too long
c
      do i = 1, nvdw
         if (nvlst(i) .ge. maxvlst) then
            write (iout,30)
   30       format (/,' VFULL  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine clist  --  build partial charge neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "clist" performs an update or a complete rebuild of the
c     electrostatic neighbor lists for partial charges
c
c
      subroutine clist
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'charge.i'
      include 'iounit.i'
      include 'neigh.i'
      integer i,j,k,ii
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
c
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
         if (octahedron) then
            do i = 1, nion
               call cbuild (i)
            end do
         else
            call clight
         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, nion
         ii = kion(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xcold(i)
         yr = yi - ycold(i)
         zr = zi - zcold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
c
c     store coordinates to reflect update of this site
c
         if (r2 .ge. lbuf2) then
            xcold(i) = xi
            ycold(i) = yi
            zcold(i) = zi
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, nion
               xr = xi - xcold(k)
               yr = yi - ycold(k)
               zr = zi - zcold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. cbuf2) then
                  j = j + 1
                  elst(j,i) = k
               end if
            end do
            nelst(i) = j
c
c     check to see if the neighbor list is too long
c
            if (nelst(i) .ge. maxelst) then
               write (iout,20)
   20          format (/,' CLIST  --  Too many Neighbors;',
     &                    ' Increase MAXELST')
               call fatal
            end if
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               xr = xi - xcold(k)
               yr = yi - ycold(k)
               zr = zi - zcold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. cbuf2) then
                  do j = 1, nelst(k)
                     if (elst(j,k) .eq. i)  goto 30
                  end do
                  nelst(k) = nelst(k) + 1
                  elst(nelst(k),k) = i
   30             continue
               else if (r2 .le. cbufx) then
                  do j = 1, nelst(k)
                     if (elst(j,k) .eq. i) then
                        elst(j,k) = elst(nelst(k),k)
                        nelst(k) = nelst(k) - 1
                        goto 40
                     end if
                  end do
   40             continue
               end if
            end do
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine cbuild  --  make charge pair list for one site  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "cbuild" performs a complete rebuild of the partial charge
c     electrostatic neighbor list for a single site
c
c
      subroutine cbuild (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'charge.i'
      include 'iounit.i'
      include 'neigh.i'
      integer i,j,k,ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     store new coordinates to reflect update of the site
c
      ii = kion(i)
      xi = x(ii)
      yi = y(ii)
      zi = z(ii)
      xcold(i) = xi
      ycold(i) = yi
      zcold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
      j = 0
      do k = i+1, nion
         kk = kion(k)
         xr = xi - x(kk)
         yr = yi - y(kk)
         zr = zi - z(kk)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .le. cbuf2) then
            j = j + 1
            elst(j,i) = k
         end if
      end do
      nelst(i) = j
c
c     check to see if the neighbor list is too long
c
      if (nelst(i) .ge. maxelst) then
         write (iout,10)
   10    format (/,' CBUILD  --  Too many Neighbors;',
     &              ' Increase MAXELST')
         call fatal
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine clight  --  make charge pair list for all sites  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "clight" performs a complete rebuild of the partial charge
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine clight
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'iounit.i'
      include 'light.i'
      include 'neigh.i'
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
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
      nlight = (ncell+1) * nion
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nion
         nelst(i) = 0
         ii = kion(i)
         xcold(i) = x(ii)
         ycold(i) = y(ii)
         zcold(i) = z(ii)
         xsort(i) = x(ii)
         ysort(i) = y(ii)
         zsort(i) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      off = sqrt(cbuf2)
      call lights (off,nion,xsort,ysort,zsort)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     loop over all atoms computing the neighbor lists
c
      do i = 1, nion
         ii = kion(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i) + 1
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            kk = kion(k)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. cbuf2) then
               if (i .lt. k) then
                  nelst(i) = nelst(i) + 1
                  elst(nelst(i),i) = k
               else
                  nelst(k) = nelst(k) + 1
                  elst(nelst(k),k) = i
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i) + 1
            stop = nlight
            goto 10
         end if
      end do
c
c     check to see if the neighbor lists are too long
c
      do i = 1, nion
         if (nelst(i) .ge. maxelst) then
            write (iout,30)
   30       format (/,' CFULL  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mlist  --  build atom multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mlist" performs an update or a complete rebuild of the
c     electrostatic neighbor lists for atomic multipoles
c
c
      subroutine mlist
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      integer i,j,k,ii
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
c
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
         if (octahedron) then
            do i = 1, npole
               call mbuild (i)
            end do
         else
            call mlight
         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xmold(i)
         yr = yi - ymold(i)
         zr = zi - zmold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
c
c     store coordinates to reflect update of this site
c
         if (r2 .ge. lbuf2) then
            xmold(i) = xi
            ymold(i) = yi
            zmold(i) = zi
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, npole
               xr = xi - xmold(k)
               yr = yi - ymold(k)
               zr = zi - zmold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. mbuf2) then
                  j = j + 1
                  elst(j,i) = k
               end if
            end do
            nelst(i) = j
c
c     check to see if the neighbor list is too long
c
            if (nelst(i) .ge. maxelst) then
               write (iout,20)
   20          format (/,' MLIST  --  Too many Neighbors;',
     &                    ' Increase MAXELST')
               call fatal
            end if
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               xr = xi - xmold(k)
               yr = yi - ymold(k)
               zr = zi - zmold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. mbuf2) then
                  do j = 1, nelst(k)
                     if (elst(j,k) .eq. i)  goto 30
                  end do
                  nelst(k) = nelst(k) + 1
                  elst(nelst(k),k) = i
   30             continue
               else if (r2 .le. mbufx) then
                  do j = 1, nelst(k)
                     if (elst(j,k) .eq. i) then
                        elst(j,k) = elst(nelst(k),k)
                        nelst(k) = nelst(k) - 1
                        goto 40
                     end if
                  end do
   40             continue
               end if
            end do
         end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mbuild  --  make mpole pair list for one site  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mbuild" performs a complete rebuild of the atomic multipole
c     electrostatic neighbor list for a single site
c
c
      subroutine mbuild (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      integer i,j,k,ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     store new coordinates to reflect update of the site
c
      ii = ipole(i)
      xi = x(ii)
      yi = y(ii)
      zi = z(ii)
      xmold(i) = xi
      ymold(i) = yi
      zmold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
      j = 0
      do k = i+1, npole
         kk = ipole(k)
         xr = xi - x(kk)
         yr = yi - y(kk)
         zr = zi - z(kk)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .le. mbuf2) then
            j = j + 1
            elst(j,i) = k
         end if
      end do
      nelst(i) = j
c
c     check to see if the neighbor list is too long
c
      if (nelst(i) .ge. maxelst) then
         write (iout,10)
   10    format (/,' MBUILD  --  Too many Neighbors;',
     &              ' Increase MAXELST')
         call fatal
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mlight  --  make mpole pair list for all sites  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mlight" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine mlight
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
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
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
      nlight = (ncell+1) * npole
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, npole
         nelst(i) = 0
         ii = ipole(i)
         xmold(i) = x(ii)
         ymold(i) = y(ii)
         zmold(i) = z(ii)
         xsort(i) = x(ii)
         ysort(i) = y(ii)
         zsort(i) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      off = sqrt(mbuf2)
      call lights (off,npole,xsort,ysort,zsort)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     loop over all atoms computing the neighbor lists
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i) + 1
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            kk = ipole(k)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. mbuf2) then
               if (i .lt. k) then
                  nelst(i) = nelst(i) + 1
                  elst(nelst(i),i) = k
               else
                  nelst(k) = nelst(k) + 1
                  elst(nelst(k),k) = i
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i) + 1
            stop = nlight
            goto 10
         end if
      end do
c
c     check to see if the neighbor lists are too long
c
      do i = 1, npole
         if (nelst(i) .ge. maxelst) then
            write (iout,30)
   30       format (/,' MFULL  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ulist  --  build preconditioner neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ulist" performs an update or a complete rebuild of the
c     neighbor lists for the polarization preconditioner
c
c
      subroutine ulist
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      integer i,j,k,ii
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
c
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
         if (octahedron) then
            do i = 1, npole
               call ubuild (i)
            end do
         else
            call ulight
         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xuold(i)
         yr = yi - yuold(i)
         zr = zi - zuold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
c
c     store coordinates to reflect update of this site
c
         if (r2 .ge. pbuf2) then
            xuold(i) = xi
            yuold(i) = yi
            zuold(i) = zi
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, npole
               xr = xi - xuold(k)
               yr = yi - yuold(k)
               zr = zi - zuold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. ubuf2) then
                  j = j + 1
                  ulst(j,i) = k
               end if
            end do
            nulst(i) = j
c
c     check to see if the neighbor list is too long
c
            if (nulst(i) .ge. maxulst) then
               write (iout,20)
   20          format (/,' ULIST  --  Too many Neighbors;',
     &                    ' Increase MAXULST')
               call fatal
            end if
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               xr = xi - xuold(k)
               yr = yi - yuold(k)
               zr = zi - zuold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. ubuf2) then
                  do j = 1, nulst(k)
                     if (ulst(j,k) .eq. i)  goto 30
                  end do
                  nulst(k) = nulst(k) + 1
                  ulst(nulst(k),k) = i
   30             continue
               else if (r2 .le. ubufx) then
                  do j = 1, nulst(k)
                     if (ulst(j,k) .eq. i) then
                        ulst(j,k) = ulst(nulst(k),k)
                        nulst(k) = nulst(k) - 1
                        goto 40
                     end if
                  end do
   40             continue
               end if
            end do
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ubuild  --  preconditioner list for one site  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ubuild" performs a complete rebuild of the polarization
c     preconditioner neighbor list for a single site
c
c
      subroutine ubuild (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      integer i,j,k,ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     store new coordinates to reflect update of the site
c
      ii = ipole(i)
      xi = x(ii)
      yi = y(ii)
      zi = z(ii)
      xuold(i) = xi
      yuold(i) = yi
      zuold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
      j = 0
      do k = i+1, npole
         kk = ipole(k)
         xr = xi - x(kk)
         yr = yi - y(kk)
         zr = zi - z(kk)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .le. ubuf2) then
            j = j + 1
            ulst(j,i) = k
         end if
      end do
      nulst(i) = j
c
c     check to see if the neighbor list is too long
c
      if (nulst(i) .ge. maxulst) then
         write (iout,10)
   10    format (/,' UBUILD  --  Too many Neighbors;',
     &              ' Increase MAXULST')
         call fatal
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ulight  --  preconditioner list for all sites  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ulight" performs a complete rebuild of the polarization
c     preconditioner pair neighbor list for all sites using the
c     method of lights
c
c
      subroutine ulight
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
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
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
      nlight = (ncell+1) * npole
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, npole
         nulst(i) = 0
         ii = ipole(i)
         xuold(i) = x(ii)
         yuold(i) = y(ii)
         zuold(i) = z(ii)
         xsort(i) = x(ii)
         ysort(i) = y(ii)
         zsort(i) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      off = sqrt(ubuf2)
      call lights (off,npole,xsort,ysort,zsort)
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     loop over all atoms computing the neighbor lists
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i) + 1
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            kk = ipole(k)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - x(kk)
            yr = yi - y(kk)
            zr = zi - z(kk)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. ubuf2) then
               if (i .lt. k) then
                  nulst(i) = nulst(i) + 1
                  ulst(nulst(i),i) = k
               else
                  nulst(k) = nulst(k) + 1
                  ulst(nulst(k),k) = i
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i) + 1
            stop = nlight
            goto 10
         end if
      end do
c
c     check to see if the neighbor lists are too long
c
      do i = 1, npole
         if (nulst(i) .ge. maxulst) then
            write (iout,30)
   30       format (/,' UFULL  --  Too many Neighbors;',
     &                 ' Increase MAXULST')
            call fatal
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine imagen  --  neighbor minimum image distance  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "imagen" takes the components of pairwise distance between
c     two points and converts to the components of the minimum
c     image distance; fast version for neighbor list generation
c     which only returns the correct component magnitudes
c
c
      subroutine imagen (xr,yr,zr)
      implicit none
      include 'boxes.i'
      real*8 xr,yr,zr
c
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         xr = abs(xr)
         yr = abs(yr)
         zr = abs(zr)
         if (xr .gt. xbox2)  xr = xr - xbox
         if (yr .gt. ybox2)  yr = yr - ybox
         if (zr .gt. zbox2)  zr = zr - zbox
c
c     for monoclinic lattice, convert "xr" and "zr" specially
c
      else if (monoclinic) then
         zr = zr / beta_sin
         yr = abs(yr)
         xr = xr - zr*beta_cos
         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
         if (yr .gt. ybox2)  yr = yr - ybox
         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
         xr = xr + zr*beta_cos
         zr = zr * beta_sin
c
c     for triclinic lattice, use general conversion equations
c
      else if (triclinic) then
         zr = zr / gamma_term
         yr = (yr - zr*beta_term) / gamma_sin
         xr = xr - yr*gamma_cos - zr*beta_cos
         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
         if (abs(yr) .gt. ybox2)  yr = yr - sign(ybox,yr)
         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
         xr = xr + yr*gamma_cos + zr*beta_cos
         yr = yr*gamma_sin + zr*beta_term
         zr = zr * gamma_term
c
c     for truncated octahedron, remove the corner pieces
c
      else if (octahedron) then
         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
         if (abs(yr) .gt. ybox2)  yr = yr - sign(ybox,yr)
         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
      return
      end

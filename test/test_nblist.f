c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_nblist  --  neighbor lists for box shapes  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_nblist" checks that every pairwise neighbor list holds all
c     the pairs it has to, for orthogonal, monoclinic, triclinic,
c     truncated octahedron and rhombic dodecahedron boxes: after a
c     complete build every pair within the cutoff plus the buffer, and
c     after an update that follows small moves, two longer moves or a
c     scaling of the box, every pair within the cutoff
c
c     the reference is a brute-force search over all pairs that builds
c     its own lattice vectors and applies the image rules of "image"
c
c     the lights method builds the lists of the prism cells and a
c     double loop those of the octahedron and dodecahedron
c
c
      subroutine test_nblist
      implicit none
c
c
      call test_nblist_key ('test_nblist_amoeba','octamoeba.key',3)
      call test_nblist_key ('test_nblist_hippo','octhippo.key',3)
      call test_nblist_key ('test_nblist_hippoc6','octhippoc6.key',3)
      call test_nblist_key ('test_nblist_tip3p','octtip3p.key',2)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_nblist_key  --  lists of one force field  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_nblist_key" loads the octahedral water box with the given
c     key, then for every box builds, updates and checks the lists in
c     use; "nlist" is the number of lists the key must use
c
c
      subroutine test_nblist_key (tname,key,nlist)
      use atoms
      use boxes
      use math
      use neigh
      implicit none
      integer i,l,nlist
      integer nlat,nused
      integer seed
      parameter (nlat=10)
      integer ishape(nlat)
      real*8 s,u
      real*8 abc(3,nlat)
      real*8 ang(3,nlat)
      real*8 test_nblist_rand
      logical skiptest
      character*(*) tname,key
      character*3 tag
c
c     boxes as shape (1 prism, 2 octahedron, 3 dodecahedron),
c     edge lengths and angles
c
      data ishape / 1, 1, 1, 1, 1, 1, 1, 2, 2, 3 /
      data abc / 26.0d0, 26.0d0, 26.0d0,
     &           25.0d0, 27.0d0, 28.0d0,
     &           26.0d0, 26.0d0, 26.0d0,
     &           26.0d0, 27.0d0, 25.0d0,
     &           26.0d0, 26.0d0, 26.0d0,
     &           26.0d0, 26.0d0, 26.0d0,
     &           26.0d0, 26.0d0, 27.0d0,
     &           26.0d0, 26.0d0, 26.0d0,
     &           30.0d0, 30.0d0, 30.0d0,
     &           26.0d0, 26.0d0, 26.0d0 /
      data ang / 90.0d0, 90.0d0, 90.0d0,
     &           90.0d0, 90.0d0, 90.0d0,
     &           90.0d0, 100.0d0, 90.0d0,
     &           90.0d0, 115.0d0, 90.0d0,
     &           85.0d0, 95.0d0, 100.0d0,
     &           70.0d0, 80.0d0, 65.0d0,
     &           90.0d0, 90.0d0, 120.0d0,
     &           90.0d0, 90.0d0, 90.0d0,
     &           90.0d0, 90.0d0, 90.0d0,
     &           90.0d0, 90.0d0, 90.0d0 /
c
c
      if (skiptest(tname,'amoeba'))  return
      call pushdir ('file/truncated8')
      call loadfix ('octwater',key)
      seed = 2026
      do l = 1, nlat
         write (tag,'(i3)')  l
c
c     set the box
c
         xbox = abc(1,l)
         ybox = abc(2,l)
         zbox = abc(3,l)
         alpha = ang(1,l)
         beta = ang(2,l)
         gamma = ang(3,l)
         octahedron = (ishape(l) .eq. 2)
         dodecadron = (ishape(l) .eq. 3)
         if (dodecadron)  zbox = xbox * root2
         call lattice
c
c     complete build: pairs within the cutoff plus the buffer
c
         dovlst = .true.
         dodlst = .true.
         doclst = .true.
         domlst = .true.
         doulst = .true.
         call nblist
         call test_nblist_all (tname//' box'//tag//' build',.true.,
     &                         nused)
         call assert_int (nused,nlist,tname//' lists in use')
         call test_nblist_skip (tname//' box'//tag,.true.)
c
c     update after moves within half the buffer and longer moves
c     of the first and the last oxygen
c
         do i = 1, n
            u = test_nblist_rand (seed)
            x(i) = x(i) + 0.2d0*u - 0.1d0
            u = test_nblist_rand (seed)
            y(i) = y(i) + 0.2d0*u - 0.1d0
            u = test_nblist_rand (seed)
            z(i) = z(i) + 0.2d0*u - 0.1d0
         end do
         x(1) = x(1) + 2.0d0
         x(n-2) = x(n-2) - 2.0d0
         call nblist
         call test_nblist_all (tname//' box'//tag//' moved',.false.,
     &                         nused)
         call test_nblist_skip (tname//' box'//tag//' moved',.false.)
c
c     update after scaling the box and the coordinates
c
         s = 1.003d0
         xbox = s * xbox
         ybox = s * ybox
         zbox = s * zbox
         call lattice
         do i = 1, n
            x(i) = s * x(i)
            y(i) = s * y(i)
            z(i) = s * z(i)
         end do
         call nblist
         call test_nblist_all (tname//' box'//tag//' scaled',.false.,
     &                         nused)
         call test_nblist_skip (tname//' box'//tag//' scaled',.false.)
      end do
      call popdir
      call final
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_nblist_all  --  check every list in use  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_nblist_all" checks each neighbor list that "nblist" keeps,
c     to the cutoff plus the buffer right after a complete build, or
c     else to the cutoff, and returns the number of lists checked
c
c
      subroutine test_nblist_all (label,built,nused)
      use atoms
      use charge
      use disp
      use limits
      use mpole
      use neigh
      use potent
      use vdw
      implicit none
      integer i,ii,nused
      real*8 cut
      real*8, allocatable :: xs(:)
      real*8, allocatable :: ys(:)
      real*8, allocatable :: zs(:)
      logical built
      character*(*) label
c
c
      allocate (xs(n))
      allocate (ys(n))
      allocate (zs(n))
      nused = 0
c
c     van der Waals sites at their reduced positions
c
      if (use_vdw .and. use_vlist) then
         cut = sqrt(vbuf2)
         if (.not. built)  cut = cut - lbuffer
         do ii = 1, nvdw
            i = ivdw(ii)
            xs(i) = xred(i)
            ys(i) = yred(i)
            zs(i) = zred(i)
         end do
         call test_nblist_one (label//' vdw',nvdw,ivdw,xs,ys,zs,
     &                         nvlst,vlst,maxvlst,cut)
         nused = nused + 1
      end if
c
c     dispersion sites
c
      if (use_disp .and. use_dlist) then
         cut = sqrt(dbuf2)
         if (.not. built)  cut = cut - lbuffer
         call test_nblist_one (label//' disp',ndisp,idisp,x,y,z,
     &                         ndlst,dlst,maxdlst,cut)
         nused = nused + 1
      end if
c
c     partial charges at their cutoff switching sites
c
      if ((use_charge.or.use_solv) .and. use_clist) then
         cut = sqrt(cbuf2)
         if (.not. built)  cut = cut - lbuffer
         do ii = 1, nion
            i = iion(ii)
            xs(i) = x(kion(i))
            ys(i) = y(kion(i))
            zs(i) = z(kion(i))
         end do
         call test_nblist_one (label//' charge',nion,iion,xs,ys,zs,
     &                         nelst,elst,maxelst,cut)
         nused = nused + 1
      end if
c
c     atomic multipoles
c
      if ((use_repel.or.use_mpole.or.use_polar
     &       .or.use_chgtrn.or.use_solv) .and. use_mlist) then
         cut = sqrt(mbuf2)
         if (.not. built)  cut = cut - lbuffer
         call test_nblist_one (label//' mpole',npole,ipole,x,y,z,
     &                         nelst,elst,maxelst,cut)
         nused = nused + 1
      end if
c
c     polarization preconditioner
c
      if (use_polar .and. use_ulist) then
         cut = sqrt(ubuf2)
         if (.not. built)  cut = cut - pbuffer
         call test_nblist_one (label//' precond',npole,ipole,x,y,z,
     &                         nulst,ulst,maxulst,cut)
         nused = nused + 1
      end if
      deallocate (xs)
      deallocate (ys)
      deallocate (zs)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_nblist_one  --  check one list by search  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_nblist_one" counts the pairs of sites within "cut" of each
c     other that the list misses; a list holds, for the atom of each
c     site, the atoms of the sites that follow it within reach
c
c
      subroutine test_nblist_one (label,nsite,isite,xs,ys,zs,
     &                            nlst,lst,maxlst,cut)
      use atoms
      use boxes
      use math
      implicit none
      integer i,j,k,ii,kk
      integer nsite,maxlst,nmiss
      integer isite(*)
      integer nlst(*)
      integer lst(maxlst,*)
      integer, allocatable :: mark(:)
      real*8 cut,c2,r2
      real*8 a,b,c,cb,cg,sg
      real*8 cy,cz,f1,f2,f3
      real*8 xr,yr,zr
      real*8 xs(*)
      real*8 ys(*)
      real*8 zs(*)
      character*(*) label
c
c
c     lattice vectors as columns (a,0,0), (b cg,b sg,0), (c cb,cy,cz)
c     of the prism, or the orthorhombic parent cell of a non-prism box
c
      a = xbox
      b = ybox
      c = zbox
      cb = 0.0d0
      cg = 0.0d0
      sg = 1.0d0
      cy = 0.0d0
      cz = c
      if (.not. nonprism) then
         cb = cos(beta/radian)
         cg = cos(gamma/radian)
         sg = sin(gamma/radian)
         cy = c * (cos(alpha/radian)-cb*cg) / sg
         cz = sqrt(c*c - (c*cb)**2 - cy*cy)
      end if
c
c     mark the neighbors of each site, then search all later sites
c
      allocate (mark(n))
      do i = 1, n
         mark(i) = 0
      end do
      c2 = cut * cut * (1.0d0-1.0d-10)
      nmiss = 0
      do ii = 1, nsite
         i = isite(ii)
         do j = 1, nlst(i)
            mark(lst(j,i)) = ii
         end do
         do kk = ii+1, nsite
            k = isite(kk)
            xr = xs(i) - xs(k)
            yr = ys(i) - ys(k)
            zr = zs(i) - zs(k)
c
c     round the fractional coordinates to the nearest cell
c
            f3 = zr / cz
            f2 = (yr - cy*f3) / (b*sg)
            f1 = (xr - b*cg*f2 - c*cb*f3) / a
            f1 = f1 - dnint(f1)
            f2 = f2 - dnint(f2)
            f3 = f3 - dnint(f3)
c
c     move the corners of the octahedron and the dodecahedron
c
            if (octahedron) then
               if (abs(f1)+abs(f2)+abs(f3) .gt. 0.75d0) then
                  f1 = f1 - sign(0.5d0,f1)
                  f2 = f2 - sign(0.5d0,f2)
                  f3 = f3 - sign(0.5d0,f3)
               end if
            else if (dodecadron) then
               if (abs(f1)+abs(f2)+2.0d0*abs(f3) .gt. 1.0d0) then
                  f1 = f1 - sign(0.5d0,f1)
                  f2 = f2 - sign(0.5d0,f2)
                  f3 = f3 - sign(0.5d0,f3)
               end if
            end if
            xr = a*f1 + b*cg*f2 + c*cb*f3
            yr = b*sg*f2 + cy*f3
            zr = cz*f3
            r2 = xr*xr + yr*yr + zr*zr
            if (r2.le.c2 .and. mark(k).ne.ii)  nmiss = nmiss + 1
         end do
      end do
      deallocate (mark)
      call assert_int (nmiss,0,label//' missing pairs')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_nblist_skip  --  atoms without a C6 site  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_nblist_skip" checks that a dispersion list update leaves
c     alone the atoms that have no dispersion site; with "init" set
c     it empties their lists and stores their current coordinates,
c     else it counts those whose list is no longer empty
c
c
      subroutine test_nblist_skip (label,init)
      use atoms
      use disp
      use limits
      use neigh
      use potent
      implicit none
      integer i,ii,nbad
      logical init
      logical, allocatable :: skip(:)
      character*(*) label
c
c
      if (.not. (use_disp.and.use_dlist))  return
      if (ndisp .eq. n)  return
      allocate (skip(n))
      do i = 1, n
         skip(i) = .true.
      end do
      do ii = 1, ndisp
         skip(idisp(ii)) = .false.
      end do
      nbad = 0
      do i = 1, n
         if (skip(i)) then
            if (init) then
               ndlst(i) = 0
               xdold(i) = x(i)
               ydold(i) = y(i)
               zdold(i) = z(i)
            else if (ndlst(i) .ne. 0) then
               nbad = nbad + 1
            end if
         end if
      end do
      deallocate (skip)
      if (.not. init) then
         call assert_int (nbad,0,label//' skipped atoms')
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function test_nblist_rand  --  repeatable random number  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_nblist_rand" returns a number in [0,1) from a linear
c     congruential sequence, so that the moves do not depend on the
c     random number seed of the run
c
c
      function test_nblist_rand (seed)
      implicit none
      integer seed
      integer*8 next
      real*8 test_nblist_rand
c
c
      next = mod(16807_8*int(seed,8),2147483647_8)
      seed = int(next)
      test_nblist_rand = dble(seed) / 2147483647.0d0
      return
      end

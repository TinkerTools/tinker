c
c
c     ##################################################################
c     ##  COPYRIGHT (C)  1990  by  Shawn Huston & Jay William Ponder  ##
c     ##              All Rights Reserved                             ##
c     ##################################################################
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine nblist  --  create pair neighbor list  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "nblist" creates a neighbor table of switching-sites via
c     the Method of Lights due to Sullivan, Mountain and O'Connell
c     (J. Comput. Phys., Vol. 61, pp 138-153, 1985)
c
c     variables:
c
c     i,j                        index for atoms
c     {ib,jb,kb,ie,je,ke,ic}     index for switching-sites
c     ik                         compressed index for switching
c     ksw                        number of switching sites
c
c     usage of routines COMPRS and GATHER:
c
c     COMPRS(v,b,u,n,ju)  loads u with the XYZ locations from
c                         array v of sites indicated by array b
c
c     GATHER(v,i,u,n)     loads u with the SX-SY-SZ locations from
c                         array v of sites indicated by array i
c
c
      subroutine nblist
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'cutoff.i'
      include 'neigh.i'
      include 'shunt.i'
      integer i,ii,j,k,ik,kk
      integer na,nnx,nnxy,nnxyz
      integer ic,icprev,prev,ksw,kc
      integer ib,jb,kb,ie,je,ke
      integer ibu,jbu,kbu,ieu,jeu,keu
      integer nfill
      integer ibs(maxatm),jbs(maxatm),kbs(maxatm)
      integer ies(maxatm),jes(maxatm),kes(maxatm)
      integer locx(maxatm),locy(maxatm),locz(maxatm)
      integer dnbrx(maxatm),dnbrxy(maxatm),dnbrxyz(maxatm)
      integer dlnbr(maxatm)
      integer rgx(maxatm),rgy(maxatm),rgz(maxatm)
      integer bt(maxatm)
      real*8 sx(maxatm),sy(maxatm),sz(maxatm)
      real*8 rx,ry,rz,add,dist
      real*8 elapsed
c
c
c     temporary timing stuff
c
      call setime
c
c     get rectangle size based on cutoff and bwidth
c
      rx = chgcut + bwidth
      ry = chgcut + bwidth
      rz = chgcut + bwidth
      if ( (rx .gt. xbox2) .or. (ry .gt. ybox2)
     &                    .or. (rz .gt. zbox2) ) then
         write (6,*)  'STOP : NBLIST - box size conflict'
         stop
      end if
c
c     create compressed x-coordinate array for switching sites
c
      ik = 0
      icprev = 0
      do i = 1, n
         ic = swsite(i)
         if (ic .ne. icprev) then
            ik = ik + 1
            sx(ik) = x(ic)
            icprev = ic
            locx(ik) = ic
         end if
      end do
c
c     sort the record {sx,locx} on the key SX
c
      call sisort(sx,locx,ik)
c
c     create compressed y-coordinate array for switching sites
c
      ik = 0
      icprev = 0
      do i = 1, n
         ic = swsite(i)
         if (ic .ne. icprev) then
            ik = ik + 1
            sy(ik) = y(ic)
            icprev = ic
            locy(ik) = ic
         end if
      end do
c
c     sort the record {sy,locy} on the key SY
c
      call sisort(sy,locy,ik)
c
c     create compressed  z-coordinate array for switching sites
c
      ik = 0
      icprev = 0
      do i = 1, n
         ic = swsite(i)
         if (ic .ne. icprev) then
            ik = ik + 1
            sz(ik) = z(ic)
            icprev = ic
            locz(ik) = ic
         end if
      end do
c
c     sort the record {sz,locz} on the key SY
c
      call sisort(sz,locz,ik)
c
c     create the arrays RGY and RGZ which point from
c     array Y into SY, and Z into SZ
c
      ksw = ik
      do i = 1, ksw
         rgx(locx(i)) = i
         rgy(locy(i)) = i
         rgz(locz(i)) = i
      end do
c
c     initialize IB pointers into SX array
c
      ik = 1
      add = -xbox
      ib = 1
   10 continue
      ib = ib + 1
      if (ib .eq. ik) goto 10
   20 continue
      dist = sx(ib) + add - sx(ik)
      if ( dist .gt. -rx ) then
         ibs(locx(ik)) = ib
         goto 40
      end if
   30 continue
      ib = ib + 1
      if (ib .gt. ksw) then
         write (6,*)  '2 OOPS!'
         stop
      end if
      goto 20
   40 continue
c
c     continue with remaining switching-sites
c
      do ik = 2, ksw
   50    continue
         dist = sx(ib) + add - sx(ik)
         if ( dist .gt. -rx) then
            ibs(locx(ik)) = ib
            goto 70
         end if
   60    continue
         ib = ib + 1
         if (ib .gt. ksw) then
            ib = 1
            add = 0
         end if
c
c     continue looping over ib for the current ik
c
         goto 50
   70    continue
      end do
c
c     initialize IE pointers into SX array
c
      ik = ksw
      add = xbox
      ie = ksw
   80 continue
      ie = ie - 1
      if (ie .eq. ik) goto 80
   90 continue
      if ( (sx(ie) + add - sx(ik)) .lt. rx ) then
         ies(locx(ik)) = ie
         goto 110
      end if
  100 continue
      ie = ie - 1
      if (ie .lt. 1) then
         write (6,*)  '6 OOPS!'
         stop
      end if
      goto 90
  110 continue
c
c     continue with remaining switching-sites
c
      do ik = ksw-1, 1, -1
  120    continue
         dist = sx(ie) + add - sx(ik)
         if ( (sx(ie) + add - sx(ik)) .lt. rx) then
            ies(locx(ik)) = ie
            goto 140
         end if
  130    continue
         ie = ie - 1
         if (ie .lt. 1) then
            ie = ksw
            add = 0
         end if
c
c     continue looping over ie for the current ik
c
         goto 120
  140    continue
      end do
c
c     initialize JB pointers into SY array
c
      ik = 1
      add = -ybox
      jb = 1
  150 continue
      jb = jb + 1
      if (jb .eq. ik) goto 150
  160 continue
      if ( (sy(jb) + add - sy(ik)) .gt. -ry ) then
         jbs(locy(ik)) = jb
         goto 180
      end if
  170 continue
      jb = jb + 1
      if (jb .gt. ksw) then
         write (6,*)  '12 OOPS!'
         stop
      end if
      goto 160
  180 continue
c
c     continue with remaining switching-sites
c
      do ik = 2, ksw
  190    continue
         if ( (sy(jb) + add - sy(ik)) .gt. -ry) then
            jbs(locy(ik)) = jb
            goto 210
         end if
  200    continue
         jb = jb + 1
         if (jb .gt. ksw) then
            jb = 1
            add = 0
         end if
c
c     continue looping over jb for the current ik
c
         goto 190
  210    continue
      end do
c
c     initialize JE pointers into SY array
c
      ik = ksw
      add = ybox
      je = ksw
  220 continue
      je = je - 1
      if (je .eq. ik) goto 220
  230 continue
      if ( (sy(je) + add - sy(ik)) .lt. ry ) then
         jes(locy(ik)) = je
         goto 250
      end if
  240 continue
      je = je - 1
      if (je .lt. 1) then
         write (6,*)  '16 OOPS!'
         stop
      end if
      goto 230
  250 continue
c
c     continue with remaining switching-sites
c
      do ik = ksw-1, 1, -1
  260    continue
         if ( (sy(je) + add - sy(ik)) .lt. ry) then
            jes(locy(ik)) = je
            goto 280
         end if
  270    continue
         je = je - 1
         if (je .lt. 1) then
            je = ksw
            add = 0
         end if
c
c     continue looping over je for the current ik
c
         goto 260
  280    continue
      end do
c
c     initialize KB pointers into SZ array
c
      ik = 1
      add = -zbox
      kb = 1
  290 continue
      kb = kb + 1
      if (kb .eq. ik) goto 290
  300 continue
      if ( (sz(kb) + add - sz(ik)) .gt. -rz ) then
         kbs(locz(ik)) = kb
         goto 320
      end if
  310 continue
      kb = kb + 1
      if (kb .gt. ksw) then
         write (6,*)  '22 OOPS!'
         stop
      end if
      goto 300
  320 continue
c
c     continue with remaining switching-sites
c
      do ik = 2, ksw
  330    continue
         if ( (sz(kb) + add - sz(ik)) .gt. -rz) then
            kbs(locz(ik)) = kb
            goto 350
         end if
  340    continue
         kb = kb + 1
         if (kb .gt. ksw) then
            kb= 1
            add = 0
         end if
c
c     continue looping over kb for the current ik
c
         goto 330
  350    continue
      end do
c
c     initialize KE pointers into SZ array
c
      ik = ksw
      add = zbox
      ke = ksw
  360 continue
      ke = ke - 1
      if (ke .eq. ik) goto 360
  370 continue
      if ( (sz(ke) + add - sz(ik)) .lt. rz ) then
         kes(locz(ik)) = ke
         goto 390
      end if
  380 continue
      ke = ke - 1
      if (ke .lt. 1) then
         write (6,*)  '116 OOPS!'
         stop
      end if
      goto 370
  390 continue
c
c     continue with remaining switching-sites
c
      do ik = ksw-1, 1, -1
  400    continue
         if ( (sz(ke) + add - sz(ik)) .lt. rz) then
            kes(locz(ik)) = ke
            goto 410
         end if
         ke = ke - 1
         if (ke .lt. 1) then
            ke = ksw
            add = 0
         end if
c
c     continue looping over ke for the current ik
c
         goto 400
  410    continue
      end do
c
c     get intersection of X,Y, and Z atom-atom
c     neighbors based on switching-site neighbors
c     N.B. an atom will be a neighbor of itself here
c
      nfill = 0
      icprev = 0
      ik = 0
      do i = 1, n-1
         do j = 1, ksw
            bt(j) = 0
         end do
         ic = swsite(i)
c
c     get compressed switching site index for particle i
c
         if (ic .ne. icprev) ik = rgx(ic)
         ibu = ibs(locx(ik))
         ieu = ies(locx(ik))
         if (ieu .gt. ibu) then
            do ii = ibu, ieu
               bt(ii) = 1
            end do
         else
            do ii = 1, ieu
               bt(ii) = 1
            end do
            do ii = ibu, ksw
               bt(ii) = 1
            end do
         end if
c
c     consolidate x-neighbors of switching-site ic
c
         call comprs (locx,bt,dnbrx,ksw,nnx)
c
c     get the SY-index of the x-neighbors
c
         call gather (rgy,dnbrx,dlnbr,ksw)
c
c     select the x-neighbor sites that are also y-neighbors
c
         jbu = jbs(locx(ik))
         jeu = jes(locx(ik))
         if (jeu .gt. jbu) then
            do j = 1, nnx
               bt(j) = ( (jbu.le.dlnbr(j)) .and. (jeu.ge.dlnbr(j)) )
            end do
         else
            do j = 1, nnx
               bt(j) = ( (jbu.le.dlnbr(j)) .or. (jeu.ge.dlnbr(j)) )
            end do
         end if
         bt(nnx+1) = 0
c
c     consolidate xy-neighbors of switching-site ic
c
         call comprs (dnbrx,bt,dnbrxy,nnx,nnxy)
c
c     get the SZ index of the xy-neighbors
c
         call gather (rgz,dnbrxy,dlnbr,ksw)
c
c     select the xy-neighbor sites that are also z-neighbors
c
         kbu = kbs(locx(ik))
         keu = kes(locx(ik))
         if (keu .gt. kbu) then
            do k = 1, nnxy
               bt(k) = ( (kbu.le.dlnbr(k)) .and. (keu.ge.dlnbr(k)) )
            end do
         else
            do k = 1, nnxy
               bt(k) = ( (kbu.le.dlnbr(k)) .or. (keu.ge.dlnbr(k)) )
            end do
         end if
c
c     consolidate xyz-neighbors of switching-site ic
c
         call comprs (dnbrxy,bt,dnbrxyz,nnxy,nnxyz)
c
c     add the neighbor atoms to particle i to the output
c     list (assumes each site has at least 1 neighbor)
c
         do j = 1, maxatm
            bt(j) = 0
         end do
         do j = 1, nnxyz
            bt(dnbrxyz(j)) = 1
         end do
         npoint(1,i) = nfill +1
c
c     fill upper triangle neighbor-list omitting diagonal
c
         do k = i+1, n
            kc = swsite(k)
            if (bt(kc) .eq. 1) then
               nfill = nfill + 1
               nlist(nfill) = k
            end if
         end do
         npoint(2,i) = nfill
         icprev = ic
      end do
      call getime(elapsed)
      write (6,*)  'NBLIST: Elapsed CPU time:',elapsed
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine gather  --  create an array u from v  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "gather" loads elements from array v into array u
c     according to the index array i
c
c
      subroutine gather (v,i,u,n)
      implicit none
      integer j,n,i(n),u(n),v(n)
c
c
      do j = 1, n
         u(j) = v(i(j))
      end do
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine comprs  --  create an array u from v  ##
c     ##                                                   ##
c     #######################################################
c
c     "comprs" loads elements from array v into array u
c     according to the bit array b
c
c     variable ju : number of elements in u (returned)
c
c
      subroutine comprs (v,b,u,n,ju)
      implicit none
      integer j,n,b(n),u(n),v(n)
      integer ib,jb,ju
c
c
      ju = 0
      do j = 1, n
         ib = b(j)
         do jb = 1, ib
            ju = ju + 1
            u(ju) = v(j)
         end do
      end do
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine switchsite  --  assign switching site  ##
c     ##                                                    ##
c     ########################################################
c
c
      subroutine switchsite
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'charge.i'
      include 'couple.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'shunt.i'
      integer i,j,k,ikey,next,freeunit
      integer mcut_typ(10),mcut_site(10)
      integer nmctype,nmcsite,soltyp
      integer molsite
      logical exist,header,use_molcut
      character*20 keyword
      character*240 keyfile
      character*240 record,string
c
c
c     open the keyword file for the current computation
c
      keyfile = filename(1:leng)//'.key'
      inquire (file=keyfile,exist=exist)
      if (exist) then
         ikey = freeunit ()
         open (unit=ikey,file=keyfile,status='old')
      end if
c
c     store default switching sites for each atom
c
      do i = 1, n
         swsite(i) = i
      end do
c
c     molecule based cutoffs may be optionally used
c     Store in kion(nion), and swsite(n) arrays
c
      use_molcut = .false.
      nmctype = 0
      nmcsite = 0
      do i = 1, 10
         mcut_typ(i) = 0
         mcut_site(i) = 0
      end do
      if (exist) then
         rewind (unit=ikey)
         dowhile (.true.)
            read (ikey,10,err=20,end=20)  record
   10       format (a240)
            next = 1
            call gettext (record,keyword,next)
            call upcase (keyword)
c
c     use MOLCUT to specify molecular cutoff (keep old keyword WATER)
c
            if ( (keyword(1:6) .eq. 'WATER ') .or.
     &           (keyword(1:7) .eq. 'MOLCUT ' ) ) then
               use_molcut = .true.
            end if
            if (keyword(1:10) .eq. 'MCUT_TYPE ') then
               nmctype = nmctype + 1
               string = record(next:240)
               read (string,*,err=20,end=20)  soltyp
               mcut_typ(nmctype) = soltyp
            end if
            if (keyword(1:10) .eq. 'MCUT_SITE ') then
               nmcsite = nmcsite + 1
               string = record(next:240)
               read (string,*,err=20,end=20)  soltyp
               mcut_site(nmcsite) = soltyp
            end if
         end do
   20    continue
      end if
      if (exist)  close (unit=ikey)
c
c     store molecular switching sites for based on
c     switching site-type and switching site atom
c
      if (use_molcut) then
         do j = 1, nmctype
            do i = 1, nion
               k = i12(1,iion(i))
               if (itype(k) .eq. mcut_typ(j)) then
                  kion(i) = k
                  swsite(iion(i)) = k
               end if
            end do
         end do
         do j = 1, nmcsite
            molsite = mcut_site(j)
            k = molcule(molsite)
            do i = imol(1,k), imol(2,k)
               swsite(i) = molsite
            end do
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine enlj1  --  Lennard-Jones energy & derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "enlj1" calculates the van der Waals energy and its first
c     derivatives with respect to Cartesian coordinates using
c     the Lennard-Jones 6-12 formalism. Interactions are counted
c     using neighbor a list. PBC, and tapering function
c     are functions of the set of switching-sites
c
c     NOTES:
c
c       1)  the replicates code at the end is not outfitted with
c     residue based switching. To implement this need to do something
c     with CELLS code.
c       2)  the reduce stuff may or may not work correctly with
c     residue based switching.
c
c
      subroutine enlj1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'couple.i'
      include 'deriv1.i'
      include 'einter.i'
      include 'energi.i'
      include 'latice.i'
      include 'molcul.i'
      include 'neigh.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
      include 'virial.i'
      integer i,j,k,ii,kk,ij,it,kt,skip(maxatm)
      integer ic,kc,kindx,k0,k1
      real*8 xi,yi,zi,xr,yr,zr
      real*8 xic,yic,zic,xc,yc,zc
      real*8 e,p6,p12,rv,eps,rdn
      real*8 redi,redii,redk,redkk
      real*8 dedx,dedy,dedz,de,dc,dedxc,dedyc,dedzc
      real*8 x_red(maxatm),y_red(maxatm),z_red(maxatm)
      real*8 rc,rc2,rc3,rc4,rc5
      real*8 rik,rik2,rik3,rik4,rik5
      real*8 taper,dtaper
c
c
c     zero out the van der Waals energy and first derivatives
c
      ev = 0.0d0
      e14 = 0.0d0
      do i = 1, n
         dev(1,i) = 0.0d0
         dev(2,i) = 0.0d0
         dev(3,i) = 0.0d0
         de14(1,i) = 0.0d0
         de14(2,i) = 0.0d0
         de14(3,i) = 0.0d0
      end do
c
c     set the coefficients for the tapering function
c
      call smooth ('VDW')
c
c     calculate "reduced" atomic coordinates
c
      do i = 1, n
         skip(i) = 0
         rdn = reduce(i)
         if (rdn.eq.0.0d0 .or. n12(i).eq.0) then
            x_red(i) = x(i)
            y_red(i) = y(i)
            z_red(i) = z(i)
         else
            ii = i12(1,i)
            x_red(i) = rdn*(x(i)-x(ii)) + x(ii)
            y_red(i) = rdn*(y(i)-y(ii)) + y(ii)
            z_red(i) = rdn*(z(i)-z(ii)) + z(ii)
         end if
      end do
c
c     find van der Waals energy and derivatives via double loop
c
      do ij = 1, navdw-1
         i = ivdw(ij)
         ic = swsite(i)
         do k = 1, n12(i)
            skip(i12(k,i)) = i
         end do
         do k = 1, n13(i)
            skip(i13(k,i)) = i
         end do
         do k = 1, n14(i)
            skip(i14(k,i)) = -i
         end do
         redi = reduce(i)
         if (redi .eq. 0.0d0) then
            ii = i
         else
            ii = i12(1,i)
            redii = 1.0d0 - redi
         end if
         it = itype(i)
         xi = x_red(i)
         yi = y_red(i)
         zi = z_red(i)
         xic = x_red(ic)
         yic = y_red(ic)
         zic = z_red(ic)
c
c     can still improve this - will have lots of zeros
c
         k0 = npoint(1,i)
         k1 = npoint(2,i)
         do kindx = k0,k1
            k = nlist(kindx)
            kc = swsite(k)
            redk = reduce(k)
            if (redk .eq. 0.0d0) then
               kk = k
            else
               kk = i12(1,k)
               redkk = 1.0d0 - redk
            end if
            if (use(i) .or. use(ii) .or. use(k) .or. use(kk)) then
               if (skip(k) .ne. i) then
                  kt = itype(k)
                  xr = xi - x_red(k)
                  yr = yi - y_red(k)
                  zr = zi - z_red(k)
                  xc = xic - x_red(kc)
                  yc = yic - y_red(kc)
                  zc = zic - z_red(kc)
                  if (use_images)  call images2 (xc,yc,zc,xr,yr,zr)
                  if (rc2 .le. off2) then
                     rik2 = xr*xr + yr*yr + zr*zr
                     rik = sqrt(rik2)
                     rv = radii(kt,it)
                     eps = epsilon(kt,it)
                     if (skip(k) .eq. -i)  eps = eps / v14scale
c
c     compute energy and derivatives for this interaction
c
                     p6 = rv**6 / rik2**3
                     p12 = p6 * p6
                     e = eps * (p12 - 2.0d0 * p6)
                     de = eps * (p12 - p6) * (-12.0d0/rik)
                     dc = 0.0d0
c
c     fifth order multiplicative tapering if near cutoff
c
                     if (rc2 .gt. cut2) then
                        rc = sqrt(rc2)
                        rc3 = rc2 * rc
                        rc4 = rc2 * rc2
                        rc5 = rc2 * rc3
                        taper = c5*rc5 + c4*rc4 + c3*rc3
     &                             + c2*rc2 + c1*rc + c0
                        dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                              + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                        dc = e * dtaper / rc
                        de = de * taper / rik
                        e = e * taper
                     end if
c
c     find the chain rule terms for derivative components
c
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     dedxc = dc * xc
                     dedyc = dc * yc
                     dedzc = dc * zc
c
c     increment the total van der Waals energy and derivatives
c
                     if (skip(k) .eq. -i) then
                        e14 = e14 + e
                        if (redi .eq. 0.0d0) then
                           de14(1,i) = de14(1,i) + dedx
                           de14(2,i) = de14(2,i) + dedy
                           de14(3,i) = de14(3,i) + dedz
                           de14(1,ic) = de14(1,ic) + dedxc
                           de14(2,ic) = de14(2,ic) + dedyc
                           de14(3,ic) = de14(3,ic) + dedzc
                        else
                           de14(1,i) = de14(1,i) + dedx*redi
                           de14(2,i) = de14(2,i) + dedy*redi
                           de14(3,i) = de14(3,i) + dedz*redi
                           de14(1,ic) = de14(1,ic) + dedxc*redi
                           de14(2,ic) = de14(2,ic) + dedyc*redi
                           de14(3,ic) = de14(3,ic) + dedzc*redi
                           de14(1,ii) = de14(1,ii) + dedx*redii
                           de14(2,ii) = de14(2,ii) + dedy*redii
                           de14(3,ii) = de14(3,ii) + dedz*redii
                        end if
                        if (redk .eq. 0.0d0) then
                           de14(1,k) = de14(1,k) - dedx
                           de14(2,k) = de14(2,k) - dedy
                           de14(3,k) = de14(3,k) - dedz
                           de14(1,kc) = de14(1,kc) - dedxc
                           de14(2,kc) = de14(2,kc) - dedyc
                           de14(3,kc) = de14(3,kc) - dedzc
                        else
                           de14(1,k) = de14(1,k) - dedx*redk
                           de14(2,k) = de14(2,k) - dedy*redk
                           de14(3,k) = de14(3,k) - dedz*redk
                           de14(1,kc) = de14(1,kc) - dedxc*redk
                           de14(2,kc) = de14(2,kc) - dedyc*redk
                           de14(3,kc) = de14(3,kc) - dedzc*redk
                           de14(1,kk) = de14(1,kk) - dedx*redkk
                           de14(2,kk) = de14(2,kk) - dedy*redkk
                           de14(3,kk) = de14(3,kk) - dedz*redkk
                        end if
                     else
                        ev = ev + e
                        if (redi .eq. 0.0d0) then
                           dev(1,i) = dev(1,i) + dedx
                           dev(2,i) = dev(2,i) + dedy
                           dev(3,i) = dev(3,i) + dedz
                           dev(1,ic) = dev(1,ic) + dedxc
                           dev(2,ic) = dev(2,ic) + dedyc
                           dev(3,ic) = dev(3,ic) + dedzc
                        else
                           dev(1,i) = dev(1,i) + dedx*redi
                           dev(2,i) = dev(2,i) + dedy*redi
                           dev(3,i) = dev(3,i) + dedz*redi
                           dev(1,ic) = dev(1,ic) + dedxc*redi
                           dev(2,ic) = dev(2,ic) + dedyc*redi
                           dev(3,ic) = dev(3,ic) + dedzc*redi
                           dev(1,ii) = dev(1,ii) + dedx*redii
                           dev(2,ii) = dev(2,ii) + dedy*redii
                           dev(3,ii) = dev(3,ii) + dedz*redii
                        end if
                        if (redk .eq. 0.0d0) then
                           dev(1,k) = dev(1,k) - dedx
                           dev(2,k) = dev(2,k) - dedy
                           dev(3,k) = dev(3,k) - dedz
                           dev(1,kc) = dev(1,kc) - dedxc
                           dev(2,kc) = dev(2,kc) - dedyc
                           dev(3,kc) = dev(3,kc) - dedzc
                        else
                           dev(1,k) = dev(1,k) - dedx*redk
                           dev(2,k) = dev(2,k) - dedy*redk
                           dev(3,k) = dev(3,k) - dedz*redk
                           dev(1,kc) = dev(1,kc) - dedxc*redk
                           dev(2,kc) = dev(2,kc) - dedyc*redk
                           dev(3,kc) = dev(3,kc) - dedzc*redk
                           dev(1,kk) = dev(1,kk) - dedx*redkk
                           dev(2,kk) = dev(2,kk) - dedy*redkk
                           dev(3,kk) = dev(3,kk) - dedz*redkk
                        end if
                     end if
c
c     increment the total intermolecular energy
c
                     if (molcule(i) .ne. molcule(k)) then
                        einter = einter + e
                     end if
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        virx = virx + xr*dedx
                        viry = viry + yr*dedy
                        virz = virz + zr*dedz
                     end if
                  end if
               end if
            end if
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replicate)  return
c
c     calculate interaction energy with other unit cells
c
      do i = 1, n
         ic = swsite(i)
         redi = reduce(i)
         if (redi .eq. 0.0d0) then
            ii = i
         else
            ii = i12(1,i)
            redii = 1.0d0 - redi
         end if
         it = itype(i)
         xi = x_red(i)
         yi = y_red(i)
         zi = z_red(i)
         k0 = npoint(1,i)
         k1 = npoint(2,i)
         do kindx = k0,k1
            k = nlist(kindx)
            kc = swsite(k)
            redk = reduce(k)
            if (redk .eq. 0.0d0) then
               kk = k
            else
               kk = i12(1,k)
               redkk = 1.0d0 - redk
            end if
            if (use(i) .or. use(ii) .or. use(k) .or. use(kk)) then
               do j = 2, ncell
                  kt = itype(k)
                  xr = xi - x_red(k)
                  yr = yi - y_red(k)
                  zr = zi - z_red(k)
                  call cells (xr,yr,zr,j)
                  rik2 = xr*xr + yr*yr + zr*zr
                  if (rik2 .le. off2) then
                     rv = radii(kt,it)
                     eps = epsilon(kt,it)
                     rik = sqrt(rik2)
c
c     compute energy and derivatives for this interaction
c
                     p6 = rv**6 / rik2**3
                     p12 = p6 * p6
                     e = eps * (p12 - 2.0d0 * p6)
                     de = eps * (p12 - p6) * (-12.0d0/rik)
c
c     fifth order multiplicative tapering if near cutoff
c
                     if (rik2 .gt. cut2) then
                        rik3 = rik2 * rik
                        rik4 = rik2 * rik2
                        rik5 = rik2 * rik3
                        taper = c5*rik5 + c4*rik4 + c3*rik3
     &                             + c2*rik2 + c1*rik + c0
                        dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                              + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                        de = e*dtaper + de*taper
                        e = e * taper
                     end if
c
c     find the chain rule terms for derivative components
c
                     de = de / rik
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ev = ev + e
                     if (redi .eq. 0.0d0) then
                        dev(1,i) = dev(1,i) + dedx
                        dev(2,i) = dev(2,i) + dedy
                        dev(3,i) = dev(3,i) + dedz
                     else
                        dev(1,i) = dev(1,i) + dedx*redi
                        dev(2,i) = dev(2,i) + dedy*redi
                        dev(3,i) = dev(3,i) + dedz*redi
                        dev(1,ii) = dev(1,ii) + dedx*redii
                        dev(2,ii) = dev(2,ii) + dedy*redii
                        dev(3,ii) = dev(3,ii) + dedz*redii
                     end if
                     if (i .ne. k) then
                        if (redk .eq. 0.0d0) then
                           dev(1,k) = dev(1,k) - dedx
                           dev(2,k) = dev(2,k) - dedy
                           dev(3,k) = dev(3,k) - dedz
                        else
                           dev(1,k) = dev(1,k) - dedx*redk
                           dev(2,k) = dev(2,k) - dedy*redk
                           dev(3,k) = dev(3,k) - dedz*redk
                           dev(1,kk) = dev(1,kk) - dedx*redkk
                           dev(2,kk) = dev(2,kk) - dedy*redkk
                           dev(3,kk) = dev(3,kk) - dedz*redkk
                        end if
                     end if
c
c     increment the total intermolecular energy
c
                     einter = einter + e
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        virx = virx + xr*dedx
                        viry = viry + yr*dedy
                        virz = virz + zr*dedz
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine encharge1  --  charge-charge energy & derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "encharge1" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     This routine loops over ATOMS (n) with use of NEIGH list
c
c
      subroutine encharge1
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'charge.i'
      include 'couple.i'
      include 'deriv1.i'
      include 'electr.i'
      include 'einter.i'
      include 'energi.i'
      include 'latice.i'
      include 'molcul.i'
      include 'neigh.i'
      include 'shunt.i'
      include 'usage.i'
      include 'virial.i'
      integer i,ii,j,k,k0,k1,kk,ic,kc,skip(maxatm)
      real*8 e,de,dc,f,fi,fik,rik,rik2
      real*8 xi,yi,zi,xr,yr,zr
      real*8 xic,yic,zic,xc,yc,zc
      real*8 dedx,dedy,dedz,dedxc,dedyc,dedzc
      real*8 rc,rc2,rc3,rc4,rc5,taper,dtaper
c     real*8 elapsed
c     call setime
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
         skip(i) = 0
      end do
      if (nion .eq. 0)  return
c
c     set conversion factor and tapering function coefficients
c
      f = 332.05382d0 / dielec
      call smooth ('CHARGE')
c
c     compute charge interaction energy and first derivatives
c
      do ii = 1, nion-1
         i = iion(ii)
         ic = swsite(i)
         do k = 1, n12(i)
            skip(i12(k,i)) = i
         end do
         do k = 1, n13(i)
            skip(i13(k,i)) = i
         end do
         do k = 1, n14(i)
            skip(i14(k,i)) = -i
         end do
         xi = x(i)
         yi = y(i)
         zi = z(i)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         fi = f * pchg(ii)
         k0 = npoint(1,i)
         k1 = npoint(2,i)
         do kk = k0,k1
            k = nlist(kk)
            kc = swsite(k)
            if (use(i) .or. use(k)) then
               if (skip(k) .ne. i) then
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
                  xc = xic - x(kc)
                  yc = yic - y(kc)
                  zc = zic - z(kc)
                  if (use_images)  call images2 (xc,yc,zc,xr,yr,zr)
                  rc2 = xc*xc + yc*yc + zc*zc
                  if (rc2 .le. off2) then
                     rik2 = xr*xr + yr*yr + zr*zr
                     rik = sqrt(rik2)
                     fik = fi * pchg(nqloc(k))
                     if (skip(k) .eq. -i)  fik = fik / chgscale
                     e = fik / rik
                     de = -fik / rik2
                     dc = 0.0d0
c
c     fifth order multiplicative tapering if near cutoff
c
                     if (rc2 .gt. cut2) then
                        rc = sqrt(rc2)
                        rc3 = rc2 * rc
                        rc4 = rc2 * rc2
                        rc5 = rc2 * rc3
                        taper = c5*rc5 + c4*rc4 + c3*rc3
     &                             + c2*rc2 + c1*rc + c0
                        dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                              + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                        dc = e * dtaper / rc
                        de = de * taper
                        e = e * taper
                     end if
c
c     form the chain rule terms for derivative expressions
c
                     de = de / rik
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     dedxc = dc * xc
                     dedyc = dc * yc
                     dedzc = dc * zc
c
c     increment the overall energy and derivative expressions
c
                     ec = ec + e
                     dec(1,i) = dec(1,i) + dedx
                     dec(2,i) = dec(2,i) + dedy
                     dec(3,i) = dec(3,i) + dedz
                     dec(1,ic) = dec(1,ic) + dedxc
                     dec(2,ic) = dec(2,ic) + dedyc
                     dec(3,ic) = dec(3,ic) + dedzc
                     dec(1,k) = dec(1,k) - dedx
                     dec(2,k) = dec(2,k) - dedy
                     dec(3,k) = dec(3,k) - dedz
                     dec(1,kc) = dec(1,kc) - dedxc
                     dec(2,kc) = dec(2,kc) - dedyc
                     dec(3,kc) = dec(3,kc) - dedzc
c
c     increment the total intermolecular energy
c
                     if (molcule(i) .ne. molcule(k)) then
                        einter = einter + e
                     end if
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        virx = virx + xr*dedx + xc*dedxc
                        viry = viry + yr*dedy + yc*dedyc
                        virz = virz + zr*dedz + zc*dedzc
                     end if
                  end if
               end if
            end if
         end do
      end do
c     call getime(elapsed)
c     write(6,*) 'ENCHARGE1: Elapsed CPU time:',elapsed
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replicate)  return
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nion-1
         i = iion(ii)
         ic = swsite(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         fi = f * pchg(ii)
         k0 = npoint(1,i)
         k1 = npoint(2,i)
         do kk = k0,k1
            k = nlist(kk)
            kc = swsite(k)
            if (use(i) .or. use(k)) then
               do j = 2, ncell
                  xc = xic - x(kc)
                  yc = yic - y(kc)
                  zc = zic - z(kc)
                  call cells (xc,yc,zc,j)
                  rc2 = xc*xc + yc*yc + zc*zc
                  if (rc2 .le. off2) then
                     xr = xi - x(k)
                     yr = yi - y(k)
                     zr = zi - z(k)
                     call cells (xr,yr,zr,j)
                     rik2 = xr*xr + yr*yr + zr*zr
                     rik = sqrt(rik2)
                     fik = fi * pchg(kk)
                     e = fik / rik
                     de = -fik / rik2
                     dc = 0.0d0
c
c     fifth order multiplicative tapering if near cutoff
c
                     if (rc2 .gt. cut2) then
                        rc = sqrt(rc2)
                        rc3 = rc2 * rc
                        rc4 = rc2 * rc2
                        rc5 = rc2 * rc3
                        taper = c5*rc5 + c4*rc4 + c3*rc3
     &                             + c2*rc2 + c1*rc + c0
                        dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                              + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                        dc = e * dtaper / rc
                        de = de * taper
                        e = e * taper
                     end if
c
c     form the chain rule terms for derivative expressions
c
                     de = de / rik
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     dedxc = dc * xc
                     dedyc = dc * yc
                     dedzc = dc * zc
c
c     increment the energy and gradient values
c
                     if (i .eq. k)  e = 0.5d0 * e
                     ec = ec + e
                     dec(1,i) = dec(1,i) + dedx
                     dec(2,i) = dec(2,i) + dedy
                     dec(3,i) = dec(3,i) + dedz
                     dec(1,ic) = dec(1,ic) + dedxc
                     dec(2,ic) = dec(2,ic) + dedyc
                     dec(3,ic) = dec(3,ic) + dedzc
                     if (i .ne. k) then
                        dec(1,k) = dec(1,k) - dedx
                        dec(2,k) = dec(2,k) - dedy
                        dec(3,k) = dec(3,k) - dedz
                        dec(1,kc) = dec(1,kc) - dedxc
                        dec(2,kc) = dec(2,kc) - dedyc
                        dec(3,kc) = dec(3,kc) - dedzc
                     end if
c
c     increment the total intermolecular energy
c
                     einter = einter + e
c
c     increment the virial for use in pressure computation
c
                     if (isobaric) then
                        virx = virx + xr*dedx + xc*dedxc
                        viry = viry + yr*dedy + yc*dedyc
                        virz = virz + zr*dedz + zc*dedzc
                     end if
                  end if
               end do
            end if
         end do
      end do
      return
      end

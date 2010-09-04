c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program distgeom  --  produce distance geometry structures  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "distgeom" uses a metric matrix distance geometry procedure to
c     generate structures with interpoint distances that lie within
c     specified bounds, with chiral centers that maintain chirality,
c     and with torsional angles restrained to desired values; the
c     user also has the ability to interactively inspect and alter
c     the triangle smoothed bounds matrix prior to embedding
c
c
      program distgeom
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'couple.i'
      include 'disgeo.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'kgeoms.i'
      include 'math.i'
      include 'refer.i'
      include 'tors.i'
      integer i,j,k,ja,kb
      integer ia,ib,ic,id
      integer swap,lext
      integer next,freeunit
      integer igeo,ngeo,nhydro
      integer r1,r2,r3,r4,b1,b2
      real*8 angle,weigh
      real*8 hbond1,hbond2
      real*8 bndfac,angfac
      real*8 radi,rmsvalue
      real*8 wall,cpu
      real*8 rab,rbc,rac,t1,t2
      real*8 qab,qbc,qcd,qac,qbd
      real*8 bndmin,bndmax
      real*8 tormin,tormax
      real*8 uppermax,big1,big2
      real*8 cosmin,cosmax
      real*8 cosabc,sinabc
      real*8 cosbcd,sinbcd
      logical exist,header,done
      logical query,info,quit
      character*1 answer
      character*1 letter
      character*7 ext
      character*120 geofile
      character*120 title
      character*120 record
      character*120 string
c
c
c     get the input structure file for the embedding
c
      call initial
      call getxyz
c
c     quit if there are too many atoms for distance geometry
c
      if (n .gt. maxgeo) then
         write (iout,10)
   10    format (/,' DISTGEOM  --  Too many Distance Geometry Atoms;',
     &              ' Increase MAXGEO')
         call fatal
      end if
c
c     set the lists of attached atoms and local interactions
c
      call attach
      call active
      call bonds
      call angles
      call torsions
c
c     get distance bound and torsional angle restraints
c
      call kgeom
c
c     store the input structure for later comparison
c
      call makeref (1)
c
c     assign approximate radii to each of the atoms
c
      do i = 1, n
         letter = name(i)(1:1)
         if (name(i) .eq. 'CH ') then
            vdwrad(i) = 1.5d0
         else if (name(i) .eq. 'CH2') then
            vdwrad(i) = 1.6d0
         else if (name(i) .eq. 'CH3') then
            vdwrad(i) = 1.7d0
         else if (letter .eq. 'H') then
            vdwrad(i) = 0.95d0
         else if (letter .eq. 'C') then
            vdwrad(i) = 1.45d0
         else if (letter .eq. 'N') then
            vdwrad(i) = 1.35d0
         else if (letter .eq. 'O') then
            vdwrad(i) = 1.35d0
         else if (letter .eq. 'P') then
            vdwrad(i) = 1.8d0
         else if (letter .eq. 'S') then
            vdwrad(i) = 1.8d0
         else
            vdwrad(i) = 0.5d0
         end if
      end do
c
c     find maximum value of vdw radii sum for an atom pair
c
      big1 = 0.0d0
      big2 = 0.0d0
      do i = 1, n
         radi = vdwrad(i)
         if (radi .gt. big1) then
            big2 = big1
            big1 = radi
         else if (radi .gt. big2) then
            big2 = radi
         end if
      end do
      vdwmax = big1 + big2
c
c     set number of distance geometry structures to generate
c
      ngeo = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=20,end=20)  ngeo
   20 continue
      if (ngeo .le. 0) then
         write (iout,30)
   30    format (/,' Number of Distance Geometry Structures',
     &              ' Desired [1] :  ',$)
         read (input,40)  ngeo
   40    format (i10)
         if (ngeo .le. 0)  ngeo = 1
      end if
c
c     enforce the original chirality of tetravalent atoms
c
      nchir = 0
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,50)
   50    format (/,' Impose Chirality Constraints on Tetrahedral',
     &              ' Atoms [Y] :  ',$)
         read (input,60)  record
   60    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .ne. 'N') then
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,70)
   70       format (/,' Use "Floating" Chirality for -XH2- and -XH3',
     &                 ' Groups [N] :  ',$)
            read (input,80)  record
   80       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         do i = 1, n
            if (n12(i) .eq. 4) then
               nhydro = 0
               if (answer .eq. 'Y') then
                  do j = 1, 4
                     letter = name(i12(j,i))(1:1)
                     if (letter .eq. 'H')  nhydro = nhydro + 1
                  end do
               end if
               if (nhydro .lt. 2) then
                  nchir = nchir + 1
                  do j = 1, 4
                     ichir(j,nchir) = i12(j,i)
                  end do
               end if
            end if
         end do
      end if
c
c     enforce the planarity or chirality of trigonal centers
c
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,90)
   90    format (/,' Impose Planarity and/or Chirality of Trigonal',
     &              ' Atoms [Y] :  ',$)
         read (input,100)  record
  100    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .ne. 'N') then
         do i = 1, n
            if (n12(i) .eq. 3) then
               nchir = nchir + 1
               do j = 1, 3
                  ichir(j,nchir) = i12(j,i)
               end do
               ichir(4,nchir) = i
            end if
         end do
      end if
c
c     enforce torsional planarity on adjacent trigonal sites
c
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,110)
  110    format (/,' Impose Torsional Planarity on Adjacent Trigonal',
     &              ' Atoms [Y] :  ',$)
         read (input,120)  record
  120    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .ne. 'N') then
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (n12(ia).eq.3 .and. n12(ib).eq.3) then
               do j = 1, n12(ia)
                  ja = i12(j,ia)
                  do k = 1, n12(ib)
                     kb = i12(k,ib)
                     if (ja.ne.ib .and. kb.ne.ia) then
                        nchir = nchir + 1
                        ichir(1,nchir) = ja
                        ichir(2,nchir) = kb
                        ichir(3,nchir) = ia
                        ichir(4,nchir) = ib
                     end if
                  end do
               end do
            end if
         end do
      end if
c
c     optionally inspect and alter the smoothed bounds matrix
c
      query = .false.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,130)
  130    format (/,' Do You Wish to Examine or Alter the Bounds',
     &              ' Matrix [N] :  ',$)
         read (input,140)  record
  140    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'Y')  query = .true.
c
c     choose the global enantiomer nearest to the original
c
      use_invert = .false.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,150)
  150    format (/,' Select the Enantiomer Closest to the Input',
     &              ' Structure [N] :  ',$)
         read (input,160)  record
  160    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'Y')  use_invert = .true.
c
c     set the type of refinement to be used after embedding
c
      use_anneal = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,170)
  170    format (/,' Refinement via Minimization or Annealing',
     &              ' [M or A, <CR>=A] :  ',$)
         read (input,180)  record
  180    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'M')  use_anneal = .false.
c
c     initialize chirality and planarity restraint values
c
      call kchiral
c
c     change the default distance restraint force constant
c
      do i = 1, ndfix
         if (dfix(1,i) .eq. 100.0d0)  dfix(1,i) = 1.0d0
      end do
c
c     print lists of the interatomic distance restraints
c
      if (verbose) then
         header = .true.
         do i = 1, ndfix
            ia = idfix(1,i)
            ib = idfix(2,i)
            weigh = dfix(1,i)
            bndmin = dfix(2,i)
            bndmax = dfix(3,i)
            if (header) then
               header = .false.
               write (iout,190)
  190          format (/,' Interatomic Distance Bound Restraints :',
     &                 //,12x,'Atom Numbers',7x,'LowerBound',
     &                    4x,'UpperBound',7x,'Weight',/)
            end if
            if (weigh .eq. 1.0d0) then
               write (iout,200)  i,ia,ib,bndmin,bndmax
  200          format (i6,5x,2i6,3x,2f14.4)
            else
               write (iout,210)  i,ia,ib,bndmin,bndmax,weigh
  210          format (i6,5x,2i6,3x,3f14.4)
            end if
         end do
c
c     print lists of the torsional angle restraints
c
         header = .true.
         do i = 1, ntfix
            ia = itfix(1,i)
            ib = itfix(2,i)
            ic = itfix(3,i)
            id = itfix(4,i)
            weigh = tfix(1,i)
            tormin = tfix(2,i)
            tormax = tfix(3,i)
            if (header) then
               header = .false.
               write (iout,220)
  220          format (/,' Intramolecular Torsional Angle Restraints :',
     &                 //,18x,'Atom Numbers',16x,'Torsion Range',
     &                    9x,'Weight',/)
            end if
            write (iout,230)  i,ia,ib,ic,id,tormin,tormax,weigh
  230       format (i6,5x,4i6,3x,3f12.4)
         end do
      end if
c
c     initialize the upper and lower bounds matrix
c
      do i = 1, n
         bnd(i,i) = 0.0d0
      end do
      uppermax = 1000000.0d0
      do i = 1, n
         do j = 1, i-1
            bnd(j,i) = uppermax
         end do
      end do
      do i = 1, n-1
         radi = vdwrad(i)
         do j = i+1, n
            bnd(j,i) = radi + vdwrad(j)
         end do
      end do
c
c     set the upper and lower bounds for 1-2 distances
c
      bndfac = 0.0d0
c     bndfac = 0.01d0
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         rab = sqrt((x(ia)-x(ib))**2 + (y(ia)-y(ib))**2
     &                     + (z(ia)-z(ib))**2)
         bndmin = (1.0d0 - bndfac) * rab
         bndmax = (1.0d0 + bndfac) * rab
         if (ia .gt. ib) then
            bnd(ia,ib) = bndmin
            bnd(ib,ia) = bndmax
         else
            bnd(ia,ib) = bndmax
            bnd(ib,ia) = bndmin
         end if
      end do
c
c     set the upper and lower bounds for 1-3 distances
c
      angfac = 0.0d0
c     angfac = 0.01d0
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         rab = sqrt((x(ia)-x(ib))**2 + (y(ia)-y(ib))**2
     &                    + (z(ia)-z(ib))**2)
         rbc = sqrt((x(ib)-x(ic))**2 + (y(ib)-y(ic))**2
     &                    + (z(ib)-z(ic))**2)
         rac = sqrt((x(ia)-x(ic))**2 + (y(ia)-y(ic))**2
     &                    + (z(ia)-z(ic))**2)
         angle = acos((rab**2+rbc**2-rac**2) / (2.0d0*rab*rbc))
         cosmin = cos(angle*(1.0d0-angfac))
         cosmax = cos(min(pi,angle*(1.0d0+angfac)))
         qab = min(bnd(ia,ib),bnd(ib,ia))
         qbc = min(bnd(ic,ib),bnd(ib,ic))
         bndmin = qab**2 + qbc**2 - 2.0d0*qab*qbc*cosmin
         bndmin = sqrt(max(0.0d0,bndmin))
         qab = max(bnd(ia,ib),bnd(ib,ia))
         qbc = max(bnd(ic,ib),bnd(ib,ic))
         bndmax = qab**2 + qbc**2 - 2.0d0*qab*qbc*cosmax
         bndmax = sqrt(max(0.0d0,bndmax))
         if (ia .gt. ic) then
            bnd(ia,ic) = bndmin
            bnd(ic,ia) = bndmax
         else
            bnd(ia,ic) = bndmax
            bnd(ic,ia) = bndmin
         end if
      end do
c
c     set the upper and lower bounds for 1-4 distances
c
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         cosmin = 1.0d0
         cosmax = -1.0d0
         do j = 1, ntfix
            r1 = itfix(1,j)
            r2 = itfix(2,j)
            r3 = itfix(3,j)
            r4 = itfix(4,j)
            if ((ia.eq.r1 .and. ib.eq.r2 .and.
     &           ic.eq.r3 .and. id.eq.r4) .or.
     &          (ia.eq.r4 .and. ib.eq.r3 .and.
     &           ic.eq.r2 .and. id.eq.r1)) then
               t1 = tfix(2,j) / radian
               t2 = tfix(3,j) / radian
               if (t2.ge.0.0d0 .and. t1.le.0.0d0) then
                  cosmin = 1.0d0
                  cosmax = min(cos(t1),cos(t2))
               else if (t1.ge.0.0d0 .and. t2.le.0.0d0) then
                  cosmin = max(cos(t1),cos(t2))
                  cosmax = -1.0d0
               else if (t1.ge.0.0d0 .and. t2.ge.t1) then
                  cosmin = cos(t1)
                  cosmax = cos(t2)
               else if (t2.le.0.0d0 .and. t1.le.t2) then
                  cosmin = cos(t2)
                  cosmax = cos(t1)
               end if
               goto 240
            end if
         end do
  240    continue
         qab = min(bnd(ia,ib),bnd(ib,ia))
         qbc = min(bnd(ib,ic),bnd(ic,ib))
         qcd = min(bnd(ic,id),bnd(id,ic))
         qac = min(bnd(ia,ic),bnd(ic,ia))
         qbd = min(bnd(ib,id),bnd(id,ib))
         cosabc = (qab**2+qbc**2-qac**2)/(2.0d0*qab*qbc)
         sinabc = sqrt(max(0.0d0,1.0d0-cosabc**2))
         cosbcd = (qbc**2+qcd**2-qbd**2)/(2.0d0*qbc*qcd)
         sinbcd = sqrt(max(0.0d0,1.0d0-cosbcd**2))
         bndmin = qab**2 + qbc**2 + qcd**2
     &               + 2.0d0*qab*qcd*cosabc*cosbcd
     &               - 2.0d0*qab*qcd*sinabc*sinbcd*cosmin
     &               - 2.0d0*qbc*(qab*cosabc+qcd*cosbcd)
         bndmin = sqrt(max(0.0d0,bndmin))
         qab = max(bnd(ia,ib),bnd(ib,ia))
         qbc = max(bnd(ib,ic),bnd(ic,ib))
         qcd = max(bnd(ic,id),bnd(id,ic))
         qac = max(bnd(ia,ic),bnd(ic,ia))
         qbd = max(bnd(ib,id),bnd(id,ib))
         cosabc = (qab**2+qbc**2-qac**2)/(2.0d0*qab*qbc)
         sinabc = sqrt(max(0.0d0,1.0d0-cosabc**2))
         cosbcd = (qbc**2+qcd**2-qbd**2)/(2.0d0*qbc*qcd)
         sinbcd = sqrt(max(0.0d0,1.0d0-cosbcd**2))
         bndmax = qab**2 + qbc**2 + qcd**2
     &               + 2.0d0*qab*qcd*cosabc*cosbcd
     &               - 2.0d0*qab*qcd*sinabc*sinbcd*cosmax
     &               - 2.0d0*qbc*(qab*cosabc+qcd*cosbcd)
         bndmax = sqrt(max(0.0d0,bndmax))
         if (ia .gt. id) then
            bnd(ia,id) = bndmin
            bnd(id,ia) = bndmax
         else
            bnd(ia,id) = bndmax
            bnd(id,ia) = bndmin
         end if
      end do
c
c     convert distance restraints into bounds matrix elements
c
      do i = 1, ndfix
         ia = idfix(1,i)
         ib = idfix(2,i)
         bndmin = dfix(2,i)
         bndmax = dfix(3,i)
         if (ia .gt. ib) then
            bnd(ia,ib) = bndmin
            bnd(ib,ia) = bndmax
         else
            bnd(ia,ib) = bndmax
            bnd(ib,ia) = bndmin
         end if
      end do
c
c     modify lower bounds to allow hydrogen bond formation
c
      hbond1 = 1.7d0
      hbond2 = 2.55d0
      do i = 1, n
         letter = name(i)(1:1)
         if (letter.eq.'N' .or. letter.eq.'O') then
            do j = 1, n
               letter = name(j)(1:1)
               if (letter .eq. 'H') then
                  k = i12(1,j)
                  letter = name(k)(1:1)
                  if (letter.eq.'N' .or. letter.eq.'O') then
                     if (j .gt. i) then
                        bnd(j,i) = min(hbond1,bnd(j,i))
                     else
                        bnd(i,j) = min(hbond1,bnd(i,j))
                     end if
                     if (k .gt. i) then
                        bnd(k,i) = min(hbond2,bnd(k,i))
                     else
                        bnd(i,k) = min(hbond2,bnd(i,k))
                     end if
                  end if
               end if
            end do
         end if
      end do
c
c     use the triangle inequalities to smooth the bounds
c
      if (verbose .and. n.le.130) then
         title = 'Input Distance Bounds :'
         call grafic (n,maxgeo,bnd,title)
      end if
      write (iout,250)
  250 format (/,' Bounds Smoothing via Triangle and Inverse',
     &           ' Triangle Inequality :')
      if (verbose)  call settime
      call geodesic
c     call triangle
      if (verbose) then
         call gettime (wall,cpu)
         write (iout,260)  wall
  260    format (/,' Time Required for Bounds Smoothing :',4x,
     &              f12.2,' seconds')
      end if
c
c     allow interactive alteration of the bounds matrix
c
      done = .false.
      do while (query .and. .not.done)
         done = .true.
         write (iout,270)
  270    format (/,' Enter an Atom Pair to Display Bounds',
     &              ' [<CR> When Done] :  ',$)
         read (input,280)  record
  280    format (a120)
         read (record,*,err=330,end=330)  b1,b2
         done = .false.
         if (b1.lt.1 .or. b2.gt.n .or. b1.eq.b2)  goto 330
         if (b1 .gt. b2) then
            swap = b1
            b1 = b2
            b2 = swap
         end if
         write (iout,290)  bnd(b2,b1),bnd(b1,b2)
  290    format (/,' Lower Bound :',f8.3,8x,'Upper Bound :',f8.3)
  300    continue
         write (iout,310)
  310    format (/,' Enter New Bounds or <CR> to Leave Unchanged :  ',$)
         read (input,320)  record
  320    format (a120)
         read (record,*,err=330,end=330)  bndmin,bndmax
         if (bndmin .gt. bndmax)  goto 300
         bnd(b2,b1) = bndmin
         bnd(b1,b2) = bndmax
         call trifix (b1,b2)
  330    continue
      end do
c
c     display the smoothed upper and lower bounds matrix
c
      if (verbose .and. n.le.130) then
         title = 'Triangle Smoothed Bounds :'
         call grafic (n,maxgeo,bnd,title)
      end if
c
c     find the largest value of an upper bound between atoms
c
      pathmax = 0.0d0
      do i = 1, n
         do j = 1, i-1
            if (pathmax .lt. bnd(j,i))  pathmax = bnd(j,i)
         end do
      end do
      write (iout,340)  pathmax
  340 format (/,' Largest Upper Bound Distance :',15x,f15.4)
c
c     check for any atoms that have no distance restraints
c
      quit = .false.
      do i = 1, n
         do j = 1, i-1
            if (bnd(j,i) .ne. uppermax)  goto 360
         end do
         do j = i+1, n
            if (bnd(i,j) .ne. uppermax)  goto 360
         end do
         quit = .true.
         write (iout,350)  i
  350    format (/,' DISTGEOM  --  Atom',i6,' has no Distance',
     &              ' Constraints')
  360    continue
      end do
      if (quit)  call fatal
c
c     generate the desired number of distance geometry structures
c
      do j = 1, ngeo
         write (iout,370)  j
  370    format (/,' Generation via Distance Geometry of Structure',i5)
         call embed
c
c     superpose the distance geometry structure on input structure
c
         info = verbose
         verbose = .false.
         call impose (n,xref,yref,zref,n,x,y,z,rmsvalue)
         verbose = info
c
c     write out the final optimized distance geometry structure
c
         lext = 3
         call numeral (j,ext,lext)
         geofile = filename(1:leng)//'.'//ext(1:lext)
         call version (geofile,'new')
         igeo = freeunit ()
         open (unit=igeo,file=geofile,status='new')
         call prtxyz (igeo)
         close (unit=igeo)
      end do
c
c     perform any final tasks before program exit
c
      call final
      end

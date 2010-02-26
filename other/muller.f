c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program muller  --  location of transition states and MERP  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "muller" is a program to implement an algorithm for the location
c     of saddle points and minimum energy paths devised Muller and Brown
c
c     literature references:
c
c     Klaus Muller and Leo D Brown, Theoretica Chimica Acta,
c     53, 75-93, (1979)
c
c     Klaus Muller, Angewandte Chemie Int. Edition, 19, 1-13, (1980)
c
c
      program muller
      include 'zcommon.for'
      integer maxpnt
      parameter (maxpnt=25)
      integer numatm,numatm2
      integer mmtype(maxatm),mmtype2(maxatm)
      integer attach(4,maxatm),attach2(4,maxatm),nattach(maxatm)
      integer q,p1,p2,blink(maxpnt),flink(maxpnt),pathpnts
      integer idim,itermax,iprnt,fc,gc,input,iout
      integer fcalls,gcalls,ncycle,maxcycle
      real*8 scale
      real*8 evalue(maxpnt),coord(3,maxatm,maxpnt)
      real*8 jdiag,eps,grdmax,fctdec,stpmax,f0min,gnorm
      logical diffcon,maxcoin,valdump
      character*1  reply
      character*2  symbol(maxatm),symbol2(maxatm)
      character*4  calctype
      character*9 coordtype
      character*80 ptitle,fname1,fname2
      character*80 conout,title1,title2
      external euclid,funct,gradt,optsave
      common /calls / fcalls,gcalls
      common /davprm/ idim,jdiag,eps,grdmax,fctdec,
     &                stpmax,fomin,itermax,iprnt
      common /iounit/ input,iout
      common /muller/ numatm,q,p1,p2,radius,itestq
      common /output/ coordtype
      common /scale / scale
      common /titles/ ptitle,calctype,fname1,fname2
      common /tstate/ evalue,coord,blink,flink
      data input  / 5 /
      data iout   / 6 /
c
c
c     *********************************************
c     *  read in the necessary input information  *
c     *********************************************
c
c     read in the title and type (ts or merp) of the calculation
c
      coordtype = 'none'
      write (iout,10)
   10 format (/'$Enter title :  ')
      read (input,20)  ptitle
   20 format (a80)
   30 continue
      write (iout,40)
   40 format (' Type of calculation to be performed :'
     &        /'     either location of a transition state (type "TS")'
     &        /'     or a minimum energy reaction pathway (type "MERP")'
     &        /'$Enter calculation type :  ')
      read (input,50)  calctype
c
c     bring in the first input structure
c
      if (calctype.ne.'ts' .and. calctype.ne.'merp')  goto 30
   50 format (a4)
      if (calctype .eq. 'ts')  then
         write (iout,60)
   60    format ('$Enter file name for structure "A" :  ')
      else
         write (iout,70)
   70    format ('$Enter file name for transition state :  ')
      end if
      read (input,80)  fname1
   80 format (a20)
      open (unit=7,name=fname1,status='old',readonly,
     &      defaultfile='.con',carriagecontrol='list')
      read (7,90)  numatm,title1
   90 format (i3,a80)
      do i = 1, numatm
         read (7,100)  symbol(i),coord(1,i,1),coord(2,i,1),
     &                 coord(3,i,1),mmtype(i),attach(1,i),
     &                 attach(2,i),attach(3,i),attach(4,i)
  100    format (1x,a2,5x,3f12.6,5i5)
         nattach(i) = 0
         do j = 1, 4
            if (attach(j,i) .ne. 0)  nattach(i) = nattach(i) + 1
         end do
      end do
      close (unit=7)
c
c     now, get the second input structure
c
      if (calctype .eq. 'ts')  then
         write (iout,110)
  110    format ('$Enter file name for structure "B" :  ')
      else
         write (iout,120)
  120    format ('$Enter file name for adjacent minimum :  ')
      end if
      read (input,130)  fname2
  130 format (a20)
      open (unit=7,name=fname2,status='old',readonly,
     &      defaultfile='.con',carriagecontrol='list')
      read (7,140)  numatm2,title2
  140 format (i3,a80)
      if (numatm .ne. numatm2) then
         stop 'The molecules contain different numbers of atoms....'
      end if
      diffcon = .false.
      do i = 1, numatm2
         read (7,150)  symbol2(i),coord(1,i,2),coord(2,i,2),
     &                 coord(3,i,2),mmtype2(i),attach2(1,i),
     &                 attach2(2,i),attach2(3,i),attach2(4,i)
  150    format (1x,a2,5x,3f12.6,5i5)
         if (symbol(i) .ne. symbol2(i))  diffcon = .true.
         if (mmtype(i) .ne. mmtype2(i))  diffcon = .true.
         do j = 1, 4
            if (attach(j,i) .ne. attach2(j,i))  diffcon = .true.
         end do
      end do
      if (diffcon) then
         stop 'The two ".CON" files have different formats....'
      end if
      close (unit=7)
      write (iout,160)
  160 format ('$Enter "Y" if molecules are to be superimposed :  ')
      read (input,170)  reply
  170 format (a1)
      maxcoin = .true.
      if (reply.eq.'n' .or. reply.eq.'N')  maxcoin = .false.
      if (calctype .eq. 'merp')  then
         write (iout,180)
  180    format ('$Enter minimum number of MERP points desired :  ')
         read (input,190)  pathpnts
  190    format (i3)
      else
         write (iout,200)
  200    format ('$Enter file name for transition state :  ')
         read (input,210) conout
  210    format (a80)
         write (iout,220)
  220    format ('$Should valley points be written out [n] :  ')
         read (input,230)  reply
  230    format (a1)
         valdump = .false.
         if (reply.eq.'y' .or. reply.eq.'Y')  valdump = .true.
      end if
c
c     now superimpose the two input structures, if user so desires
c
c     ....do rigidfit here....
c
c     ***************************************************************
c     *  initialize muller algorithms and the minimization routine  *
c     ***************************************************************
c
      natom = numatm
      do i = 1, natom
         itype(i) = mmtype(i)
         n12(i) = nattach(i)
         do j = 1, nattach(i)
            i12(j,i) = attach(j,i)
         end do
      end do
      call initmm
c
c     now, set up parameters needed by the minimizer "davidon"
c
      scale = 5.0d0
      idim = 3 * natom
      jdiag = 0.4d0
      eps = 1.0d-5
      grdmax = min(1.0d0,dfloat(natom)/50.0d0)
      fctdec = 0.01d0
      stpmax = max(1.0d0,dfloat(natom)/20.0d0)
      f0min = 0.0d0
      itermax = 25
      iprnt = 1
      fc = 0
      gc = 0
      fcalls = 0
      gcalls = 0
      maxcycle = 5
c
c     finally, initialize the doubly-linked list of path points
c
      q = 2
      p1 = 1
      p2 = 2
      blink(1) = 1
      flink(1) = 2
      blink(2) = 1
      flink(2) = 2
      call setcoord (coord(1,1,1))
      evalue(1) = energy()
      call setcoord (coord(1,1,2))
      evalue(2) = energy()
      if (calctype .eq. 'merp')  then
         write (iout,240)  evalue(1),(j,symbol(j),
     &                     (coord(i,j,1),i=1,3),mmtype(j),j=1,natom)
  240    format (/' Transition State Structure :  Energy and Location',
     &           /'     energy :  ',f10.4,
     &           /'     coords :  ',i5,3x,a2,3f12.6,i5,
     &            <natom-1>(/'               ',i5,3x,a2,3f12.6,i5))
         write (iout,250)  evalue(2),(j,symbol2(j),
     &                     (coord(i,j,2),i=1,3),mmtype2(j),j=1,natom)
  250    format (/' Adjacent Local Minimum :  Energy and Location',
     &           /'     energy :  ',f10.4,
     &           /'     coords :  ',i5,3x,a2,3f12.6,i5,
     &            <natom-1>(/'               ',i5,3x,a2,3f12.6,i5))
         goto 380
      else
         if (evalue(1) .lt. evalue(2))  then
            p1 = 2
            p2 = 1
         end if
         write (iout,260)  1,evalue(1),(j,symbol(j),
     &                     (coord(i,j,1),i=1,3),mmtype(j),j=1,natom)
         write (iout,260)  2,evalue(2),(j,symbol2(j),
     &                     (coord(i,j,2),i=1,3),mmtype2(j),j=1,natom)
  260    format (/' Starting Point ',i3,' :  Energy and Location',
     &           /'     Energy :  ',f10.4,
     &           /'     Coords :  ',i5,3x,a2,3f12.6,i5,
     &            <natom-1>(/'               ',i5,3x,a2,3f12.6,i5))
      end if
c
c     *********************************************
c     *  transition state localization algorithm  *
c     *********************************************
c
      radfac = 2.0d0 / 3.0d0
      quit = 0.1d0
c
c     if the distance between the 'p1' and 'p2' points is less than
c     'quit', then we are done and can print the results
c
  270 radius = euclid(coord(1,1,p1),coord(1,1,p2)) * radfac
      if (radius .lt. quit*radfac)  then
         call results
         if (conout .ne. ' ') then
            open (unit=7,file=conout,status='new',
     &            defaultfile='.con',carriagecontrol='list')
            write (7,280)  natom,ptitle(1:60)
  280       format (i3,a60)
            do i = 1, natom
               write (7,290)  symbol(i),i,coord(1,i,p1),coord(2,i,p1),
     &                        coord(3,i,p1),mmtype(i),
     &                        (attach(j,i),j=1,nattach(i))
  290          format (1x,a2,i5,3f12.6,5i5)
            end do
            close (unit=7)
         end if
         goto 440
      end if
c
c     set up the next starting point for a hypersphere search
c
      q = q + 1
      itestq = 0
      do j = 1, natom
         do i = 1, 3
            coord(i,j,q) = (1.0d0-radfac) * coord(i,j,p1) +
     &                      radfac * coord(i,j,p2)
         end do
      end do
      call setcoord (coord(1,1,q))
      evalue(q) = energy()
      distq1 = euclid(coord(1,1,q),coord(1,1,1))
      distq2 = euclid(coord(1,1,q),coord(1,1,2))
      write (iout,300)  q,p1,p2,evalue(q),distq1,distq2
  300 format (/' Searching for Valley Point ',i3,
     &         ' between ',i3,' and ',i3,
     &        //' Starting Point for Hypersphere Search :',
     &        /'     Energy :  ',f11.4,
     &        /'     Distance from 1 :  ',f8.4,5x,
     &         '     Distance from 2 :  ',f8.4)
      do j = 1, natom
         do i = 1, 3
            coord(i,j,q) = coord(i,j,q) * scale
         end do
      end do
c
c     perform the search and locate the next 'valley point'
c
      ncycle = 1
      call davidon (idim,coord(1,1,q),evalue(q),
     &              grdmax,funct,gradt,optsave)
      do j = 1, natom
         do i = 1, 3
            gnorm = gnorm + g(i,j)**2
         end do
      end do
      dowhile (gnorm.gt.grdmax .and. ncycle.lt.maxcycle)
         ncycle = ncycle + 1
         call davidon (idim,coord(1,1,q),evalue(q),
     &                 grdmin,funct,gradt,optsave)
         evalue(q) = davidon (idim,coord(1,1,q),jdiag,eps,
     &                        grdmax,fctdec,stpmax,f0min,
     &                        itermax,iprnt,fc,gc,gnorm)
      end do
      do j = 1, natom
         do i = 1, 3
            coord(i,j,q) = coord(i,j,q) / scale
         end do
      end do
      call project1 (coord(1,1,q))
      call testq
      call insert (q,p1,p2)
      distq1 = euclid(coord(1,1,q),coord(1,1,1))
      distq2 = euclid(coord(1,1,q),coord(1,1,2))
      write (iout,310)  q,evalue(q),distq1,distq2
  310 format (/' Valley Point ',i3,' :  Energy and Location',
     &        /'     Energy :  ',f11.4,
     &        /'     Distance from 1 :  ',f8.4,5x,
     &         '     Distance from 2 :  ',f8.4)
      if (valdump) then
         open (unit=7,file='valpnt.con',status='new',
     &         carriagecontrol='list')
         write (7,320)  natom,q
  320    format (i3,'Valley Point ',i3)
         do i = 1, natom
            write (7,330)  symbol(i),i,coord(1,i,q),coord(2,i,q),
     &                     coord(3,i,q),mmtype(i),
     &                     (attach(j,i),j=1,nattach(i))
  330       format (1x,a2,i5,3f12.6,5i5)
         end do
         close (unit=7)
      end if
c
c     now, find out which type of 'valley point' we have just
c     located; but first check "itestq" to make sure that
c     we didn't use the remedy procedure to get this point
c
      if (itestq .eq. 1)  then
         p1 = q
c
c     case a:  evalue(q) less than either evalue(p1) or evalue(p2)
c
      else if (evalue(q) .lt. evalue(p1) .and.
     &         evalue(q) .lt. evalue(p2))  then
         write (iout,340)  q
  340    format (/' For your Information:  "MULLER" has determined'
     &           /'     that a Local Minimum exists in the vicinity'
     &           /'     of Valley Point ',i3)
         write (iout,350)  ((coord(i,j,q),i=1,3),j=1,natom)
  350    format (/'     coords :  ',3f12.6,
     &           <natom-1>(/'               ',3f12.6))
         if (conout .ne. ' ') then
            open (unit=7,file='locmin.con',
     &            carriagecontrol='list',status='new')
            write (7,360)  natom,q
  360       format (i3,'Valley Point ',i3,'  --  near a Local Minimum')
            do i = 1, natom
               write (7,370)  symbol(i),i,coord(1,i,q),coord(2,i,q),
     &                        coord(3,i,q),mmtype(i),
     &                        (attach(j,i),j=1,nattach(i))
  370          format (1x,a2,i5,3f12.6,5i5)
            end do
            close (unit=7)
         end if
c
c     case b:  evalue(q) greater than either evalue(p1) and evalue(p2)
c
      else if (evalue(q) .ge. evalue(p1) .and.
     &         evalue(q) .ge. evalue(p2))  then
         p2 = p1
         p1 = q
         goto 270
      end if
c
c     case c:  evalue(q) between evalue(p1) and evalue(p2)
c               (also case a:  evalue(q) less than both)
c
      distblink = euclid(coord(1,1,p1),coord(1,1,blink(p1)))
      distflink = euclid(coord(1,1,p1),coord(1,1,flink(p1)))
      if (distblink .gt. distflink)  then
         p2 = blink(p1)
      else
         p2 = flink(p1)
      end if
      goto 270
c
c     ***************************************************
c     *  minimum energy reaction path (merp) algorithm  *
c     ***************************************************
c
  380 radfac = 1.0d0 / (pathpnts + 1.0d0)
      radius = euclid(coord(1,1,p1),coord(1,1,p2)) * radfac
      shutoff = (6.0d0 * radius) / 5.0d0
  390 q = q + 1
      itestq = 0
      if (q .eq. 3) then
         do j = 1, natom
            do i = 1, 3
               coord(i,j,3) = (1.0d0-radfac) * coord(i,j,p1) +
     &                         radfac * coord(i,j,p2)
            end do
         end do
      else
         radfac = radius / euclid(coord(1,1,p1),coord(1,1,p2))
         do j = 1, natom
            do i = 1, 3
               coord(i,j,q) = (1.0d0-radfac) * coord(i,j,p1) +
     &                         radfac * coord(i,j,p2)
            end do
         end do
      end if
      call setcoord (coord(1,1,q))
      evalue(q) = energy()
      distq1 = euclid(coord(1,1,q),coord(1,1,1))
      distq2 = euclid(coord(1,1,q),coord(1,1,2))
      write (iout,400)  q,evalue(q),distq1,distq2
  400 format (/' Searching for Minimum Energy Path Point ',i3,' :',
     &        //' Starting Point for Hypersphere Search :',
     &        /'     Energy :  ',f11.4,
     &        /'     Distance from TS :  ',f8.4,5x,
     &         '     Distance from Min ;  ',f8.4)
      do j = 1, natom
         do i = 1, 3
            coord(i,j,q) = coord(i,j,q) * scale
         end do
      end do
      evalue(q) = davidon (idim,coord(1,1,q),jdiag,eps,grdmax,fctdec,
     &                     stpmax,f0min,itermax,iprnt,fc,gc,gnorm)
      do j = 1, natom
         do i = 1, 3
            coord(i,j,q) = coord(i,j,q) / scale
         end do
      end do
      call project1 (coord(1,1,q))
      call testq
      call insert (q,p1,p2)
      distq1 = euclid(coord(1,1,q),coord(1,1,1))
      distq2 = euclid(coord(1,1,q),coord(1,1,2))
      write (iout,410)  q,evalue(q),distq1,distq2
  410 format (/' MERP Point ',i3,' :  Energy and Location',
     &        /'     Energy :  ',f11.4,
     &        /'     Distance from TS :  ',f8.4,5x,
     &         '     Distance from Min :  ',f8.4)
      if (euclid(coord(1,1,q),coord(1,1,p2)) .lt. shutoff)  then
         call results
         goto 440
      end if
      if (itestq .eq. 1)  then
         write (iout,420)
  420    format (/' For Your Information :  "MULLER" has determined'
     &           /'     that the Reaction Path is moving away from the'
     &           /'     initially supplied Local Minimum, indicating'
     &           /'     that this Minimum is not "adjacent" to the'
     &           /'     Transition State -- run will be aborted now !!')
         call results
         goto 440
      end if
      if (evalue(q) .gt. evalue(p1))  then
         write (iout,430)
  430    format (/' For Your Information :  the Minimum Energy Path'
     &           /'     point "MULLER" just found is higher in energy'
     &           /'     than the immediately preceeding MERP Point,'
     &           /'     indicating that the initially supplied Local'
     &           /'     Minimum is not "adjacent" to the Transition'
     &           /'     State -- the run will be aborted now ....')
         call results
         goto 440
      end if
c
c     now set merp point just found to be next p1 point, and return to
c     the top of the algorithm to start on a new merp point search
c
      p1 = q
      goto 390
c
c     perform any final tasks before program exit
c
  440 continue
      call final
      end
c
c
c     ######################
c     ##                  ##
c     ##  function funct  ##
c     ##                  ##
c     ######################
c
c
c     function funct --- a function, called by the minimization
c     routine "davidon", which evaluates the energy at a
c     potential surface point specified by "davidon"
c
c
      function funct (x)
      include 'sizes.for'
      integer numatm,q,p1,p2,itestq
      integer fcalls,gcalls
      real*8 xx(3,maxatm),x(3*maxatm)
      real*8 funct
      real*8 radius
      real*8 scale
      common /calls / fcalls,gcalls
      common /muller/ numatm,q,p1,p2,radius,itestq
      common /scale / scale
c
c
c     increment the number of function calls
c
      fcalls = fcalls + 1
c
c     copy the linear input vector into a (3 x numatm) array
c
      l = 0
      do j = 1, numatm
         do i = 1, 3
            l = l + 1
            xx(i,j) = x(l) / scale
         end do
      end do
c
c     project the point onto the appropriate hypersphere
c     (or hyperplane if we are doing "remedy")
c
      if (itestq .eq. 0)  then
         call project1 (xx)
      else
         call project2 (xx)
      end if
c
c     copy coordinates and get the energy
c
      call setcoord (xx)
      funct = energy()
c
c     project the point passed to/from "davidon"
c     onto the current hypershpere
c
      l = 0
      do j = 1, numatm
         do i = 1, 3
            l = l + 1
            x(l) = xx(i,j) * scale
         end do
      end do
      return
      end
c
c
c     ########################
c     ##                    ##
c     ##  subroutine gradt  ##
c     ##                    ##
c     ########################
c
c
c     subroutine gradt --- evaluates the gradient components of
c     the energy at a potential surface point specified by
c     the minimization routine "davidon"
c
c
      subroutine gradt (x,g)
      include 'sizes.for'
      integer maxpnt
      parameter (maxpnt=25)
      integer numatm,q,p1,p2,itestq
      integer fcalls,gcalls,input,iout
      real*8 xx(3,maxatm),x(*)
      real*8 gg(3,maxatm),g(*),derivs(3,maxatm)
      real*8 radius,dot,length
      real*8 evalue(maxpnt),coord(3,maxatm,maxpnt)
      real*8 midpoint(3,maxatm),norm(3,maxatm)
      real*8 actual,projec
      real*8 scale
      common /calls / fcalls,gcalls
      common /iounit/ input,iout
      common /muller/ numatm,q,p1,p2,radius,itestq
      common /remedy/ midpoint,norm
      common /scale / scale
      common /tstate/ evalue,coord
c
c
c     increment the number of gradient calls
c
      gcalls = gcalls + 1
c
c     copy the linear input vector into a 3xNumatm array
c
      l = 0
      do j = 1, numatm
         do i = 1, 3
            l = l + 1
            xx(i,j) = x(l) / scale
         end do
      end do
c
c     project the point onto the appropriate hypersphere
c     (or hyperplane if we are doing "remedy")
c
      if (itestq .eq. 0)  then
         call project1 (xx)
      else
         call project2 (xx)
      end if
c
c     get the gradients and copy them into the 'gg' array
c
      call setcoord (xx)
      call gradient (derivs)
      do j = 1, numatm
         do i = 1, 3
            gg(i,j) = derivs(i,j)
         end do
      end do
c
c     project the gradient onto the hyperplane tangent to the
c     current hypersphere (or onto the current hyperplane if
c     we are doing "remedy")
c
      l = 0
      if (itestq .eq. 0) then
         dot = 0.0d0
         do j = 1, numatm
            do i = 1, 3
               dot = dot + gg(i,j) * (xx(i,j)-coord(i,j,p1))
            end do
         end do
         length = dot / radius**2
         do j = 1, numatm
            do i = 1, 3
               l = l + 1
               g(l) = gg(i,j) - length * (xx(i,j)-coord(i,j,p1))
               g(l) = g(l) / scale
            end do
         end do
      else
         dot = 0.0d0
         do j = 1, numatm
            do i = 1, 3
               dot = dot + gg(i,j) * norm(i,j)
            end do
         end do
         length = dot
         do j = 1, numatm
            do i = 1, 3
               l = l + 1
               g(l) = gg(i,j) - length * norm(i,j)
               g(l) = g(l) / scale
            end do
         end do
      end if
c
c     section to compute the actual and projected gradient
c     norms as a check on the gradient projection method
c
c     l = 0
c     actual = 0.0d0
c     projec = 0.0d0
c     do j = 1, numatm
c        do i = 1, 3
c           actual = actual + gg(i,j) * gg(i,j)
c           l = l + 1
c           projec = projec + g(l) * g(l)
c        end do
c     end do
c     actual = sqrt(actual)
c     projec = sqrt(projec)
c     write (iout,10)  actual,projec
c  10 format (' Actual Gradient Norm :',f11.4,8x,
c    1        ' Projected Gradient Norm :',f11.4)
      return
      end
c
c
c     ###########################
c     ##                       ##
c     ##  subroutine setcoord  ##
c     ##                       ##
c     ###########################
c
c
      subroutine setcoord (xx)
      include 'zcommon.for'
      real*8 xx(3,maxatm)
c
c
      do j = 1, natom
         x(j) = xx(1,j)
         y(j) = xx(2,j)
         z(j) = xx(3,j)
      end do
      return
      end
c
c
c     ###########################
c     ##                       ##
c     ##  subroutine project1  ##
c     ##                       ##
c     ###########################
c
c
c     subroutine project1 --- moves a point "x" from open hyperspace
c     back onto the hypersphere defined by the current "p1" point
c     and the value of radius; in the case where "x" is the current
c     "q" point, then "x" is just "angles(1,q)"
c
c
      subroutine project1 (x)
      include 'sizes.for'
      integer maxpnt
      parameter (maxpnt=25)
      integer numatm,q,p1,p2
      real*8 x(3,maxatm),length,normal,radius
      real*8 evalue(maxpnt),coord(3,maxatm,maxpnt)
      common /muller/ numatm,q,p1,p2,radius
      common /tstate/ evalue,coord
c
c
      length = euclid(x,coord(1,1,p1))
      normal = radius / length
      do j = 1, numatm
         do i = 1, 3
            x(i,j) = (1.0d0-normal)*coord(i,j,p1) + normal*x(i,j)
         end do
      end do
      return
      end
c
c
c     #######################
c     ##                   ##
c     ##  function euclid  ##
c     ##                   ##
c     #######################
c
c
c     function euclid --- finds the euclidian distance between two
c     points, "x" and "y", whose angles are stored in the arrays
c     x(i) and y(i) where i = 1, idim
c
c
      function euclid (x,y)
      integer numatm,idim
      real*8 x(*),y(*)
      common /muller/ numatm
c
c
      idim = 3 * numatm
      sumsqu = 0.0d0
      do i = 1, idim
         a = abs(x(i) - y(i))
         sumsqu = sumsqu + a**2
      end do
      euclid = sqrt(sumsqu)
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine insert  ##
c     ##                     ##
c     #########################
c
c
c     subroutine insert --- places the current "q" point in its proper
c     place in the doubly linked list of valley points found to date;
c     ie, "q" is inserted between the current "p1" and "p2" points
c
c
      subroutine insert (q,p1,p2)
      include 'sizes.for'
      integer maxpnt
      parameter (maxpnt=25)
      integer q,p1,p2,blink(maxpnt),flink(maxpnt),t
      real*8 evalue(maxpnt),coord(3,maxatm,maxpnt)
      common /tstate/ evalue,coord,blink,flink
c
c
      if (flink(p1) .eq. p2)  then
         t = p1
      else
         t = p2
      end if
      flink(q) = flink(t)
      blink(q) = t
      flink(t) = q
      blink(flink(q)) = q
      return
      end
c
c
c     ########################
c     ##                    ##
c     ##  subroutine testq  ##
c     ##                    ##
c     ########################
c
c
c     subroutine testq --- tests the current "q" point to see if it is
c     between the current "p1" and "p2" points; ie, the test is passed
c     if "q" is closer to "p2" than it is to the pole of "p2" reflected
c     through "p1" - otherwise we set itestq = 1 and call "remedy" if
c     "calctype" is "ts", return to "muller" and abort if "calctype"
c     is "merp"
c
c
      subroutine testq
      include 'sizes.for'
      integer maxpnt
      parameter (maxpnt=25)
      integer numatm,q,p1,p2,itestq
      real*8 hypotnus,radius
      real*8 evalue(maxpnt),coord(3,maxatm,maxpnt)
      character*4  calctype
      character*80 ptitle
      common /muller/ numatm,q,p1,p2,radius,itestq
      common /titles/ ptitle,calctype
      common /tstate/ evalue,coord
c
c
      hypotnus = radius**2 + euclid(coord(1,1,p1),coord(1,1,p2))**2
      hypotnus = sqrt(hypotnus)
      if (euclid(coord(1,1,q),coord(1,1,p2)) .ge. hypotnus)  then
         itestq = 1
         if (calctype .eq. 'ts')  call remedy
      end if
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine remedy  ##
c     ##                     ##
c     #########################
c
c
c     subrouitne remedy --- if "testq" has determined that the current
c     "q" point is not between the "p1" and "p2" points, then "remedy"
c     fixes the situation by minimizing in a plane orthogonal to (and
c     containing the midpoint of the segment p1p2 ; both p1 and p2 are
c     subsequently deleted from the linked list of valley points and
c     the minimum we just found is then used is the next "p1" point
c
c
      subroutine remedy
      include 'sizes.for'
      integer maxpnt
      real*8 scale
      parameter (maxpnt=25)
      parameter (scale=5.0)
      integer numatm,q,p1,p2,input,iout,idim
      integer itermax,iprnt,fc,gc,ncycle,maxcycle
      integer blink(maxpnt),flink(maxpnt)
      real*8 midpoint(3,maxatm),norm(3,maxatm)
      real*8 evalue(maxpnt),coord(3,maxatm,maxpnt)
      real*8 jdiag,eps,grdmax,fctdec,stpmax,fomin,gnorm
      common /davprm/ idim,jdiag,eps,grdmax,fctdec,
     &                stpmax,f0min,itermax,iprnt
      common /iounit/ input,iout
      common /muller/ numatm,q,p1,p2
      common /remedy/ midpoint,norm
      common /tstate/ evalue,coord,blink,flink
c
c
      write (iout,10)
   10 format (' ***** Remedy Procedure Invoked *****')
      dist = euclid(coord(1,1,p1),coord(1,1,p2))
      do j = 1, numatm
         do i = 1, 3
            midpoint(i,j) = (coord(i,j,p1) + coord(i,j,p2)) / 2.0d0
            coord(i,j,q) = midpoint(i,j)
            norm(i,j) = (coord(i,j,p1) - coord(i,j,p2)) / dist
         end do
      end do
      do j = 1, numatm
         do i = 1, 3
            coord(i,j,q) = coord(i,j,q) * scale
         end do
      end do
      maxcycle = 5
      ncycle = 1
      evalue(q) = davidon (idim,coord(1,1,q),jdiag,eps,
     &                     grdmax,fctdec,stpmax,f0min,
     &                     itermax,iprnt,fc,gc,gnorm)
      dowhile (gnorm.gt.grdmax .and. ncycle.lt.maxcycle)
         ncycle = ncycle + 1
         evalue(q) = davidon (idim,coord(1,1,q),jdiag,eps,
     &                        grdmax,fctdec,stpmax,f0min,
     &                        itermax,iprnt,fc,gc,gnorm)
      end do
      do j = 1, numatm
         do i = 1, 3
            coord(i,j,q) = coord(i,j,q) / scale
         end do
      end do
      call project2 (coord(1,1,q))
      if (flink(p1) .eq. p2)  then
         p1 = blink(p1)
         p2 = flink(p2)
         flink(p1) = p2
         blink(p2) = p1
      else
         p1 = flink(p1)
         p2 = blink(p2)
         flink(p2) = p1
         blink(p1) = p2
      end if
      return
      end
c
c
c     ##########################
c     ##                      ##
c     ##  subroutine results  ##
c     ##                      ##
c     ##########################
c
c
      subroutine results
      include 'sizes.for'
      integer maxpnt
      parameter (maxpnt=25)
      integer numatm,q,p1
      integer t,flink(maxpnt),blink(maxpnt)
      integer fcalls,gcalls,input,iout
      real*8 evalue(maxpnt),coord(3,maxatm,maxpnt)
      character*4  calctype
      character*80 ptitle,fname1,fname2
      common /calls / fcalls,gcalls
      common /iounit/ input,iout
      common /muller/ numatm,q,p1
      common /titles/ ptitle,calctype,fname1,fname2
      common /tstate/ evalue,coord,blink,flink
c
c
      if (calctype .eq. 'ts')  then
         write (iout,10)  ptitle(1:55),fname1(1:20),fname2(1:20)
   10    format (//'1',75('*'),
     &            /1x,5('*'),65x,5('*'),
     &            /1x,5('*'),5x,a55,5x,5('*'),
     &            /1x,5('*'),65x,5('*'),
     &            /1x,5('*'),5x,'Location of a Transition State ',
     &                'between the Structures :',5x,5('*'),
     &            /1x,5('*'),65x,5('*'),
     &            /1x,5('*'),10x,a20,10x,a20,5x,5('*'),
     &            /1x,5('*'),65x,5('*'),
     &            /1x,75('*'))
      else
         write (iout,20)  ptitle(1:55),fname1(1:20),fname2(1:20)
   20    format (//'1',78('*'),
     &            /1x,5('*'),68x,5('*'),
     &            /1x,5('*'),5x,a58,5x,5('*'),
     &            /1x,5('*'),68x,5('*'),
     &            /1x,5('*'),5x,'Determination of the Minimum Energy ',
     &                'Reaction Path (MERP) :',5x,5('*'),
     &            /1x,5('*'),68x,5('*'),
     &            /1x,5('*'),5x,'Transition State = ',a17,
     &                'Minimum = ',a17,5('*'),
     &            /1x,5('*'),68x,5('*'),
     &            /1x,78('*'))
      end if
      t = 1
      if (calctype .eq. 'ts') then
         write (iout,30)
   30    format (//' MULLER''s Final Results :  "TS" Mode',
     &           //' Point        Energy       Dist from TS',
     &             '    Dist from 1     Dist from 2')
      else
         write (iout,40)
   40    format (//' MULLER''s Final Results :  "MERP" Mode',
     &           //' Point        Energy       Dist from TS',
     &             '    Dist from Min')
      end if
      dowhile (t .ne. 0)
         distt1 = euclid(coord(1,1,t),coord(1,1,1))
         distt2 = euclid(coord(1,1,t),coord(1,1,2))
         if (calctype .eq. 'ts') then
            disttp1 = euclid(coord(1,1,t),coord(1,1,p1))
            write (iout,50)  t,evalue(t),disttp1,distt1,distt2
   50       format (i4,5x,f11.4,5x,f11.4,5x,f11.4,5x,f11.4)
         else
            write (iout,60)  t,evalue(t),distt1,distt2
   60       format (i4,5x,f11.4,5x,f11.4,5x,f11.4)
         end if
         if (t .eq. flink(t))  then
            t = 0
         else
            t = flink(t)
         end if
      end do
      if (calctype .eq. 'ts')  then
         write (iout,70)
   70    format (//' Transition State Structure :')
         write (iout,80)  (j,(coord(i,j,p1),i=1,3),j=1,numatm)
   80    format (/'     Coords :  ',i5,3f12.6,
     &           <numatm-1>(/'               ',i5,3f12.6))
      end if
      write (iout,90)  fcalls,gcalls
   90 format (/' Total Function Calls :',i6,5x,
     &         ' Total Gradient Calls :',i6)
      if (calctype .eq. 'ts')  then
         write (iout,100)
  100    format (//2(1h,75('*')/)//)
      else
         write (iout,110)
  110    format (//2(1h,78('*')/)/)
      end if
      return
      end
c
c
c     ###########################
c     ##                       ##
c     ##  subroutine project2  ##
c     ##                       ##
c     ###########################
c
c
c     subroutine project2 --- used when we are doing "remedy", this
c     routine projects the current "q" point found in open hyperspace
c     onto a plane orthogonal to (and containing the midpoint of) the
c     segment p1p2 ; upon return the current "q" point is set to the
c     projection point we just found
c
c
      subroutine project2 (x)
      include 'sizes.for'
      integer numatm,idim
      real*8 x(3,maxatm)
      real*8 midpoint(3,maxatm),norm(3,maxatm),shadow
      common /muller/ numatm
      common /remedy/ midpoint,norm
c
c
      shadow = 0.0d0
      do j = 1, numatm
         do i = 1, 3
            shadow = shadow + (x(i,j)-midpoint(i,j)) * norm(i,j)
         end do
      end do
      do j = 1, numatm
         do i = 1, 3
            x(i,j) = x(i,j) - (shadow * norm(i,j))
         end do
      end do
      return
      end

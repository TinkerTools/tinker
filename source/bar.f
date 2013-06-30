c
c
c     ###########################################################
c     ##  COPYRIGHT (C) 2009 by C. Wu, P. Ren & Jay W. Ponder  ##
c     ##                  All Rights Reserved                  ##
c     ###########################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program bar  --  free energy perturbation via BAR method  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "bar" computes the free energy difference between two states
c     via the Bennett acceptance ratio (BAR) method
c
c     literature reference:
c
c     C. H. Bennett, "Efficient Estimation of Free Energy Differences
c     from Monte Carlo Data", Journal of Computational Physics, 22,
c     245-268 (1976)
c
c
      program bar
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'files.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'mutant.i'
      include 'polar.i'
      include 'potent.i'
      include 'units.i'
      integer maxstep,maxlam
      parameter (maxstep=20000)
      parameter (maxlam=25)
      integer ibar,iarc
      integer i,j,k
      integer start,stop,step
      integer now,nsnap
      integer next,size
      integer nlamda
      integer bailout
      integer freeunit,trimtext
      integer nframe(maxlam)
      real*8 temp,rt
      real*8 nom,denom,diff
      real*8 eold,enew,fe
      real*8 energy
      real*8 vlamda(maxlam)
      real*8 elamda(maxlam)
      real*8 ewin(3,maxstep,maxlam)
      real*8 delta0,delta(maxlam)
      real*8 oldpole(maxpole,maxatm)
      real*8 oldpolar(maxatm)
      real*8 oldpchg(maxatm)
      character*120 record
      character*120 string
      character*120 barinfile
      character*120 arcfile(maxlam)
      logical exist,query
      logical quit,abort
      logical getcoord,arcok
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
      if (vlambda.lt.1.0d0 .or. elambda.lt.1.0d0) then
         write (iout,10)
   10    format (/,1x,'BAR -- Error! Please set VLAMBDA and ELAMBDA',/
     &             1x,'       to 1 or REMOVE them from the keyfile.',/
     &             1x,' The input parameters need to be original.'/
     &             1x,' ...Quit!')
         call fatal
      end if
c
c     store the original electrostatic parameters
c
      if (use_mpole) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               do j = 1, 13
                  oldpole(j,i) = pole(j,i)
               end do
            end if
         end do
      end if
      if (use_polar) then
         do i = 1, npolar
            if (mut(i))  oldpolar(i) = polarity(i)
         end do
      end if
      if (use_charge) then
         do i = 1, nion
            if (mut(i))  oldpchg(i) = pchg(i)
         end do
      end if
c
c     obtain the target temperature value
c
      temp = 0.0d0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=20,end=20)  temp
         query = .false.
      end if
   20 continue
      if (query) then
         write (iout,30)
   30    format (/,' Enter the System Temperature [298 K] :  ',$)
         read (input,40)  temp
   40    format (f20.0)
      end if
      if (temp .eq. 0.0d0)  temp = 298.0d0
      rt = gasconst * temp
c
c     read the file containing in the following format
c     *****
c     nlamda
c     elamda vlamda archivefilename
c     elamda vlamda archivefilename
c     ...
c     ...
c     elamda vlamda archivefilenane (totally nlamda such lines)
c     *****
c
      call nextarg (barinfile,exist)
      if (.not. exist) then
         write (iout,50)
         read (input,60,err=70,end=70)  barinfile
      end if
   50 format (/,' Enter Name of the BAR Input File :',$)
   60 format (a120)
   70 continue
c
c     frames to be used for free energy estimate
c
      start = 0
      stop = 0
      step = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=80,end=80)  start
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  stop
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  step
   80 continue
      if (query) then
         write (iout,90)
   90    format (/,' First & Last Snapshot and Step',
     &              ' [<CR>=Exit] :  ',$)
         read (input,100)  record
  100    format (a120)
         read (record,*,err=110,end=110)  start,stop,step
  110    continue
      end if
      if (stop .eq. 0)  stop = start
      if (step .eq. 0)  step = 1
c
c     initialize values in BAR input file
c
      nlamda = 0
      do i = 1, maxlam
         vlamda(i) = 1.0d0
         elamda(i) = 1.0d0
      end do
c
c     read data from BAR input file
c
      ibar = freeunit ()
      open (unit=ibar,file=barinfile,status='old')
      rewind (unit=ibar)
      size = 0
      dowhile (size .eq. 0)
         read (ibar,120,err=160,end=160)  record
  120    format (a120)
         size = trimtext (record)
      end do
      abort = .false.
      quit = .true.
c
c     parse the title line to get the number of atoms
c
      i = 0
      next = 1
      call gettext (record,string,next)
      read (string,*,err=160,end=160)  nlamda
c
c     check for too many total windows/lamdas in the file
c
      if (nlamda .gt. maxlam) then
         write (iout,130)  maxlam
  130    format (/,' BAR -- The Maximum of',i8,' LAMBDAs',
     &              ' has been Exceeded')
         call fatal
      end if
c
c     read the elamda, vlamda and archive file names
c
      do i = 1, nlamda
         next = 1
         size = 0
         dowhile (size .eq. 0)
            read (ibar,140,err=160,end=160)  record
  140       format (a120)
            size = trimtext (record)
         end do
         call gettext (record,string,next)
         read (string,*,err=160,end=160) elamda(i)
         call gettext (record,string,next)
         read (string,*,err=160,end=160) vlamda(i)
         call gettext (record,arcfile(i),next)
  150    continue
      end do
      quit = .false.
  160 continue
      close (ibar)
c
c     initialize energies
c
      do i = 1, 3
         do j = 1, maxstep
            do k = 1, maxlam
               ewin (i,j,k) = 0.0d0
            end do
         end do
         nframe(i) = 0
      end do
c
c     calculate energies ewin(1,k): u(0,-1):ewin(2,k):u(0-0),ewin(3,k):(
c
      do i = 1, nlamda
         write (iout,170)  i
  170    format (' Processing the No. ',i3,' Window ...')
         iarc = freeunit()
         open (unit=iarc,file=arcfile(i),status='old')
         rewind (iarc)
         arcok = .true.
         if (start .ne. 0) then
            do k = 1, start-1
               arcok = getcoord (iarc)
            end do
            j = start
            nsnap = 0
            do while (j.ge.start .and. j.le.stop)
               arcok = getcoord (iarc)
               if (n .eq. 0)  then
                  if (j .eq. 1) then
                     write (iout,*)  ' No Snapshot Found'
                     nsnap = 0
                  else
                     write(iout,190) j-step, stop
                  end if
                  goto 200
               end if
               if (.not. arcok) then
                  write (iout,180) nsnap+1
  180             format (' BAR -- Stopped Reading Archive',
     &                       ' File at Snapshot ',i5)
                  goto 200
               end if
               nsnap = nsnap + 1
c
c     e(lamda) at previous lamda
c
               if(i .gt. 1) then
                  elambda = elamda(i-1)
                  vlambda = vlamda(i-1)
                  call resume (oldpole,oldpolar,oldpchg)
                  call altelec
                  ewin(1,j,i) = energy()
               end if
c
c     e(lamda) at curret lamda
c
               elambda = elamda(i)
               vlambda = vlamda(i)
               call resume (oldpole,oldpolar,oldpchg)
               call altelec
               ewin(2,j,i) = energy()
c
c     e(lamda) at next lamda
c
               if (i .lt. nlamda) then
                  elambda = elamda(i+1)
                  vlambda = vlamda(i+1)
                  call resume (oldpole,oldpolar,oldpchg)
                  call altelec
                  ewin(3,j,i) = energy ()
               end if
               j = j + step
               do k = 1, step-1
                  arcok = getcoord(iarc)
               end do
            end do
         end if
  190    format (//1x,i8," ( <",I8,
     &    ") is the last snapshot in the archive file",
     &    /,10x,"Reading archive file stopped!")
  200    continue
         now = stop
         close (unit=iarc)
         nframe(i) = nsnap
         write (iout,210)  i,nframe(i)
  210    format (1x,'Window ',i3,' has ',i5,' Snapshots')
      end do
c
c     calculate forwarding free energies for each window
c
      fe = 0.0d0
      do i = 1, nlamda - 1
         nom = 0.0d0
         denom = 0.0d0
         delta0 = 1.0d0 / dble(nframe(i))
         eold = delta0
         enew = delta0
         diff = 1.0d0 / dble(nframe(i))
         bailout = 0
         dowhile (diff .ge. 0.001/dble(nframe(i))
     &               .and. bailout .lt. 100)
            eold = enew
            do j = 1, nframe(i)
               denom = denom
     &         + 1.0d0/(1+exp((ewin(3,j,i)-ewin(2,j,i)-eold)/rt))
            end do
            denom = denom / dble(nframe(i))
            do j = 1, nframe(i+1)
               nom = nom
     &         + 1.0d0/(1+exp((ewin(1,j,i+1)-ewin(2,j,i+1)+eold)/rt))
            end do
            nom = nom/dble(nframe(i+1))
            delta(i) = rt*log(nom/denom) + eold
            enew = delta(i)
            diff = abs(enew-eold)
            bailout = bailout + 1
            write (iout,220)  i,bailout,enew,eold,diff
  220       format (' Window ',i3,' Run ',i3,' FE_new ',
     &                 f10.4,' FE_old ',f10.4,' Diff ',f10.4)
            if (bailout .eq. 100) then
               write (iout,230)
  230          format (/,' The simulation completed 100'
     &                    ' iterations without converging',/)
            end if
         end do
         write (iout,*)
         fe = fe + delta(i)
      end do
      write (iout,240)
  240 format(/,' The free energies for different windows:')
      do i = 1, nlamda-1
         write (iout,250)  i,i+1,delta(i)
  250    format (' Window ',i3,'  --',i3,f10.4,' kcal/mole')
      end do
      write (iout,260)  fe
  260 format (/,' The final free energy is : ',f15.8,' kcal/mole')
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine resume --  resume electrostatic parameters  ##
c     ##                                                         ##
c     #############################################################
c
c
      subroutine resume (oldpole,oldpolar,oldpchg)
      implicit none
      include 'sizes.i'
      include 'charge.i'
      include 'mpole.i'
      include 'mutant.i'
      include 'polar.i'
      include 'potent.i'
      integer i,j,k
      real*8 oldpole(maxpole,maxatm)
      real*8 oldpolar(maxatm)
      real*8 oldpchg(maxatm)
c
c
      if (use_mpole) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               do j = 1, 13
                  pole(j,i) = oldpole(j,i)
               end do
            end if
         end do
      end if
      if (use_polar) then
         do i = 1, npolar
            if (mut(i))  polarity(i) = oldpolar(i)
         end do
      end if
      if (use_charge) then
         do i = 1, nion
            if (mut(i))  pchg(i) = oldpchg(i)
         end do
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function getcoord   --  input of Cartesian coordinates  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "getcoord" gets a set of Cartesian coordinates from
c     an external disk file, and return a bool value
c
c
      function getcoord (ixyz)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'titles.i'
      integer i,j,k,m,ixyz
      integer next,size
      integer first,last
      integer nexttext
      integer trimtext
      integer list(maxatm)
      logical exist,opened
      logical quit,reorder
      logical clash,getcoord
      character*120 xyzfile
      character*120 record
      character*120 string
c
c
c     initialize the total number of atoms in the system
c
      n = 0
c
c     open the input file if it has not already been done
c
      inquire (unit=ixyz,opened=opened)
      if (.not. opened) then
         xyzfile = filename(1:leng)//'.xyz'
         call version (xyzfile,'old')
         inquire (file=xyzfile,exist=exist)
         if (exist) then
            open (unit=ixyz,file=xyzfile,status='old')
            rewind (unit=ixyz)
         else
            write (iout,10)
   10       format (/,' GETCOORD -- Unable to Find the Cartesian',
     &                 ' Coordinates File')
            getcoord = .false.
         end if
      end if
c
c     read first line and return if already at end of file
c
      quit = .false.
      abort = .true.
      size = 0
      dowhile (size .eq. 0)
         read (ixyz,20,err=60,end=60)  record
   20    format (a120)
         size = trimtext (record)
      end do
      abort = .false.
      quit = .true.
      getcoord = .true.
c
c     parse the title line to get the number of atoms
c
      i = 0
      next = 1
      call gettext (record,string,next)
      read (string,*,err=60,end=60)  n
c
c     extract the title and determine its length
c
      string = record(next:120)
      first = nexttext (string)
      last = trimtext (string)
      if (last .eq. 0) then
         title = ' '
         ltitle = 0
      else
         title = string(first:last)
         ltitle = trimtext (title)
      end if
c
c     check for too many total atoms in the file
c
      if (n .gt. maxatm) then
         write (iout,30)  maxatm
   30    format (/,' GETCOORD -- The Maximum of',i8,' Atoms',
     &              ' has been Exceeded')
         getcoord = .false.
      end if
c
c     initialize coordinates and connectivities for each atom
c
      do i = 1, n
         x(i) = 0.0d0
         y(i) = 0.0d0
         z(i) = 0.0d0
      end do
c
c     read the coordinates and connectivities for each atom
c
      do i = 1, n
         next = 1
         size = 0
         dowhile (size .eq. 0)
            read (ixyz,40,err=60,end=60)  record
   40       format (a120)
            size = trimtext (record)
         end do
         read (record,*,err=60,end=60)  tag(i)
         call getword (record,name(i),next)
         string = record(next:120)
         read (string,*,err=50,end=50)  x(i),y(i),z(i)
   50    continue
      end do
      quit = .false.
   60 continue
      if (.not. opened)  close (unit=ixyz)
c
c     an error occurred in reading the coordinate file
c
      if (quit) then
         write (iout,70)  i
   70    format (/,' GETCOORD -- Error in Coordinate File at Atom',i6)
         getcoord = .false.
      end if
      return
      end

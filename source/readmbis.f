c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2024  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine readmbis  --  input of MBIS multipole analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "readmbis" takes the Minimal Basis Iterative Stockholder (MBIS)
c     output as Cartesian multipoles from the Multiwfn program and
c     converts to Tinker format
c
c     this version assumes Multiwfn was invoked using Frank Jensen's
c     MBIS method to compute atomic multipoles
c
c
      subroutine readmbis (ichg,imbis)
      use atomid
      use atoms
      use files
      use iounit
      use mpole
      use units
      implicit none
      integer i,j,k,next
      integer ichg,imbis
      integer freeunit
      logical exist,done
      logical openchg
      logical openmbis
      character*3 atmnam
      character*240 record
      character*240 string
      character*240 chgfile
      character*240 mbisfile
c
c
c     open the charge and coordinates file if not already done
c
      inquire (unit=ichg,opened=openchg)
      if (.not. openchg) then
         chgfile = filename(1:leng)//'.chg'
         call version (chgfile,'old')
         inquire (file=chgfile,exist=exist)
         if (exist) then
            open (unit=ichg,file=chgfile,status='old')
            rewind (unit=ichg)
         else
            call nextarg (chgfile,exist)
            if (exist) then
               call basefile (chgfile)
               call suffix (chgfile,'chg','old')
               inquire (file=chgfile,exist=exist)
            end if
            do while (.not. exist)
               write (iout,10)
   10          format (/,' Enter CHG Output File Name :  ',$)
               read (input,20)  chgfile
   20          format (a240)
               call basefile (chgfile)
               call suffix (chgfile,'chg','old')
               inquire (file=chgfile,exist=exist)
            end do
         end if
      end if
c
c     open the MBIS atomic multipole file if not already done
c
      inquire (unit=imbis,opened=openmbis)
      if (.not. openmbis) then
         mbisfile = filename(1:leng)//'.mbis'
         call version (mbisfile,'old')
         inquire (file=mbisfile,exist=exist)
         if (exist) then
            open (unit=imbis,file=mbisfile,status='old')
            rewind (unit=imbis)
         else
            call nextarg (mbisfile,exist)
            if (exist) then
               call basefile (mbisfile)
               call suffix (mbisfile,'mbis_mpl','old')
               inquire (file=mbisfile,exist=exist)
            end if
            do while (.not. exist)
               write (iout,30)
   30          format (/,' Enter MBIS Output File Name :  ',$)
               read (input,40)  mbisfile
   40          format (a240)
               call basefile (mbisfile)
               call suffix (mbisfile,'mbis_mpl','old')
               inquire (file=mbisfile,exist=exist)
            end do
         end if
      end if
c
c     first open and then read the charge output file
c
      ichg = freeunit ()
      open (unit=ichg,file=chgfile,status='old')
      rewind (unit=ichg)
c
c     get the number of atoms and the atomic coordinates
c
      i = 0
      do while (.true.)
         read (ichg,50,err=60,end=60)  record
   50    format (a240)
         i = i + 1
         next = 1
         call gettext (record,name(i),next)
         string = record(next:240)
         read (string,*,err=60,end=60)  x(i),y(i),z(i)
      end do
   60 continue
      n = i
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(rpole))  allocate (rpole(maxpole,n))
c
c     now open and then read the MBIS output file
c
      imbis = freeunit ()
      open (unit=imbis,file=mbisfile,status='old')
      rewind (unit=imbis)
c
c     get the atomic multipole values from MBIS output file
c
      do while (.true.)
         read (imbis,70,err=110,end=110)  record
   70    format (a240)
         if (record(3:16) .eq. 'Atomic charges') then
            do i = 1, n
               read (imbis,80,err=110,end=110)  record
   80          format (a240)
               next = 1
               call gettext (record,name(i),next)
               string = record(next:240)
               read (string,*,err=110,end=110)  rpole(1,i)
            end do
         end if
         if (record(3:16) .eq. 'Atomic dipoles') then
            do i = 1, n
               read (imbis,90,err=110,end=110)  record
   90          format (a240)
               next = 1
               call gettext (record,name(i),next)
               string = record(next:240)
               read (string,*,err=110,end=110)  (rpole(j,i),j=2,4)
            end do
         end if
         if (record(3:31) .eq. 'Atomic quadrupoles, Traceless') then
            do i = 1, n
               read (imbis,100,err=110,end=110)  record
  100          format (a240)
               next = 1
               call gettext (record,name(i),next)
               string = record(next:240)
               read (string,*,err=110,end=110)  rpole(5,i),rpole(6,i),
     &                                          rpole(7,i),rpole(9,i),
     &                                          rpole(10,i),rpole(13,i)
               rpole(8,i) = rpole(6,i)
               rpole(11,i) = rpole(7,i)
               rpole(12,i) = rpole(10,i)
            end do
         end if
      end do
  110 continue
c
c     attempt to get atomic numbers from Multiwfn atom names
c
      do i = 1, n
         atomic(i) = 0
         atmnam = name(i)
         call upcase (atmnam)
         if (atmnam(1:2) .eq. 'SI') then
            atomic(i) = 14
         else if (atmnam(1:2) .eq. 'CL') then
            atomic(i) = 17
         else if (atmnam(1:2) .eq. 'BR') then
            atomic(i) = 35
         else if (atmnam(1:1) .eq. 'H') then
            atomic(i) = 1
         else if (atmnam(1:1) .eq. 'B') then
            atomic(i) = 5
         else if (atmnam(1:1) .eq. 'C') then
            atomic(i) = 6
         else if (atmnam(1:1) .eq. 'N') then
            atomic(i) = 7
         else if (atmnam(1:1) .eq. 'O') then
            atomic(i) = 8
         else if (atmnam(1:1) .eq. 'F') then
            atomic(i) = 9
         else if (atmnam(1:1) .eq. 'P') then
            atomic(i) = 15
         else if (atmnam(1:1) .eq. 'S') then
            atomic(i) = 16
         else if (atmnam(1:1) .eq. 'I') then
            atomic(i) = 53
         else
            read (atmnam,*,err=120,end=120)  atomic(i)
  120       continue
         end if
      end do
c
c     print the global frame Cartesian atomic multipoles
c
      write (iout,130)
  130 format (/,' Global Frame Cartesian Multipole Moments :')
      do i = 1, n
         write (iout,140)  i,name(i),atomic(i)
  140    format (/,' Atom:',i8,9x,'Name:',3x,a3,7x,'Atomic Number:',i8)
         write (iout,150)  x(i),y(i),z(i)
  150    format (/,' Coordinates:',5x,3f15.6)
         write (iout,160)  rpole(1,i)
  160    format (/,' Charge:',10x,f15.5)
         write (iout,170)  rpole(2,i),rpole(3,i),rpole(4,i)
  170    format (' Dipole:',10x,3f15.5)
         write (iout,180)  rpole(5,i)
  180    format (' Quadrupole:',6x,f15.5)
         write (iout,190)  rpole(8,i),rpole(9,i)
  190    format (18x,2f15.5)
         write (iout,200)  rpole(11,i),rpole(12,i),rpole(13,i)
  200    format (18x,3f15.5)
      end do
c
c     convert the dipole and quadrupole moments to Angstroms,
c     quadrupole divided by 3 for use as traceless values
c
      do i = 1, n
         do k = 2, 4
            rpole(k,i) = rpole(k,i) * bohr
         end do
         do k = 5, 13
            rpole(k,i) = rpole(k,i) * bohr**2 / 3.0d0
         end do
      end do
c
c     close the MBIS multipole analysis output file
c
      if (.not. openchg)  close (unit=imbis)
      if (.not. openmbis)  close (unit=imbis)
      return
      end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2008  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine readgdma  --  input of GDMA multipole analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "readgdma" takes the DMA output in spherical harmonics from
c     the GDMA program and converts to Cartesian multipoles in
c     the global coordinate frame
c
c     this version is compatible with the formatted output from
c     GDMA 2.2.04 released by Anthony Stone in Fall 2008; it also
c     reads the DMA output from Psi4 as it calls the GDMA code
c
c
      subroutine readgdma (idma)
      use atomid
      use atoms
      use dma
      use files
      use iounit
      use mpole
      use units
      implicit none
      integer i,j,k
      integer idma,next
      integer freeunit
      real*8 term
      logical exist,opened
      logical done,use_bohr
      character*3 atmnam
      character*240 record
      character*240 dmafile
c
c
c     open the input file if it has not already been done
c
      inquire (unit=idma,opened=opened)
      if (.not. opened) then
         dmafile = filename(1:leng)//'.dma'
         call version (dmafile,'old')
         inquire (file=dmafile,exist=exist)
         if (exist) then
            open (unit=idma,file=dmafile,status='old')
            rewind (unit=idma)
         else
            call nextarg (dmafile,exist)
            if (exist) then
               call basefile (dmafile)
               call suffix (dmafile,'dma','old')
               inquire (file=dmafile,exist=exist)
            end if
            do while (.not. exist)
               write (iout,10)
   10          format (/,' Enter GDMA Output File Name :  ',$)
               read (input,20)  dmafile
   20          format (a240)
               call basefile (dmafile)
               call suffix (dmafile,'dma','old')
               inquire (file=dmafile,exist=exist)
            end do
         end if
      end if
c
c     first open and then read the GDMA output file
c
      idma = freeunit ()
      open (unit=idma,file=dmafile,status='old')
c
c     count the number of atoms in the GDMA output file
c
      i = 0
      rewind (unit=idma)
      do while (.true.)
         read (idma,30,err=40,end=40)  record
   30    format (a240)
         if (record(12:14) .eq. 'x =') then
            i = i + 1
         else if (record(1:16) .eq. 'Total multipoles') then
            goto 40
         end if
      end do
   40 continue
      n = i
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(mp))  allocate (mp(n))
      if (.not. allocated(dpx))  allocate (dpx(n))
      if (.not. allocated(dpy))  allocate (dpy(n))
      if (.not. allocated(dpz))  allocate (dpz(n))
      if (.not. allocated(q20))  allocate (q20(n))
      if (.not. allocated(q21c))  allocate (q21c(n))
      if (.not. allocated(q21s))  allocate (q21s(n))
      if (.not. allocated(q22c))  allocate (q22c(n))
      if (.not. allocated(q22s))  allocate (q22s(n))
c
c     zero out the atomic coordinates and DMA values
c
      do i = 1, n
         x(i) = 0.0d0
         y(i) = 0.0d0
         z(i) = 0.0d0
         mp(i) = 0.0d0
         dpx(i) = 0.0d0
         dpy(i) = 0.0d0
         dpz(i) = 0.0d0
         q20(i) = 0.0d0
         q21c(i) = 0.0d0
         q21s(i) = 0.0d0
         q22c(i) = 0.0d0
         q22s(i) = 0.0d0
      end do
c
c     get coordinates and multipoles from GDMA output file
c
      i = 0
      rewind (unit=idma)
      do while (.true.)
         read (idma,50,err=70,end=70)  record
   50    format (a240)
         if (i .ne. 0)  call match1 (i,record)
         if (record(12:14) .eq. 'x =') then
            i = i + 1
            next = 1
            call gettext (record,name(i),next)
            read (record(15:24),*)  x(i)
            read (record(30:39),*)  y(i)
            read (record(45:54),*)  z(i)
            read (idma,60,err=70,end=70)
   60       format ()
         else if (record(1:16) .eq. 'Total multipoles') then
            goto 70
         end if
      end do
   70 continue
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(rpole))  allocate (rpole(maxpole,n))
c
c     convert quadrupole from spherical harmonic to Cartesian
c
      term = sqrt(0.75d0)
      do i = 1, n
         rpole(1,i) = mp(i)
         rpole(2,i) = dpx(i)
         rpole(3,i) = dpy(i)
         rpole(4,i) = dpz(i)
         rpole(5,i) = -0.5d0*q20(i) + term*q22c(i)
         rpole(6,i) = term*q22s(i)
         rpole(7,i) = term*q21c(i)
         rpole(8,i) = rpole(6,i)
         rpole(9,i) = -0.5d0*q20(i) - term*q22c(i)
         rpole(10,i) = term*q21s(i)
         rpole(11,i) = rpole(7,i)
         rpole(12,i) = rpole(10,i)
         rpole(13,i) = q20(i)
      end do
c
c     check for GDMA coordinate values in atomic units
c
      use_bohr = .false.
      rewind (unit=idma)
      do while (.true.)
         read (idma,80,err=90,end=90)  record
   80    format (a240)
         if (record(1:27) .eq. 'Positions and radii in bohr') then
            use_bohr = .true.
            goto 90
         end if
      end do
   90 continue
c
c     convert coordinates from Bohrs to Angstroms if needed
c
      if (use_bohr) then
         do i = 1, n
            x(i) = x(i) * bohr
            y(i) = y(i) * bohr
            z(i) = z(i) * bohr
         end do
      end if
c
c     find atomic numbers in verbose GDMA output if available
c
      done = .false.
      rewind (unit=idma)
      do while (.true.)
         read (idma,100,err=120,end=120)  record
  100    format (a240)
         if (record(1:16) .eq. 'Nuclear charges:') then
            k = min(n,20)
            read (record(17:240),*,err=120,end=120)  (atomic(i),i=1,k)
            do while (k .ne. n)
               j = k + 1
               k = min(n,k+20)
               read (idma,110,err=120,end=120)  record
  110          format (a240)
               read (record,*,err=120,end=120)  (atomic(i),i=j,k)
            end do
            done = .true.
         end if
      end do
  120 continue
      close (unit=idma)
c
c     attempt to get atomic numbers from GDMA atom names
c
      if (.not. done) then
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
               read (atmnam,*,err=130,end=130)  atomic(i)
  130          continue
            end if
         end do
      end if
c
c     print the global frame Cartesian atomic multipoles
c
      write (iout,140)
  140 format (/,' Global Frame Cartesian Multipole Moments :')
      do i = 1, n
         write (iout,150)  i,name(i),atomic(i)
  150    format (/,' Site:',i8,9x,'Name:',3x,a3,7x,'Atomic Number:',i8)
         write (iout,160)  x(i),y(i),z(i)
  160    format (/,' Coordinates:',5x,3f15.6)
         write (iout,170)  rpole(1,i)
  170    format (/,' Charge:',10x,f15.5)
         write (iout,180)  rpole(2,i),rpole(3,i),rpole(4,i)
  180    format (' Dipole:',10x,3f15.5)
         write (iout,190)  rpole(5,i)
  190    format (' Quadrupole:',6x,f15.5)
         write (iout,200)  rpole(8,i),rpole(9,i)
  200    format (18x,2f15.5)
         write (iout,210)  rpole(11,i),rpole(12,i),rpole(13,i)
  210    format (18x,3f15.5)
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
c     close the GDMA multipole analysis output file
c
      if (.not. opened)  close (unit=idma)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine match1  --  match first value from GDMA output  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "match1" finds and stores the first multipole component found
c     on a line of output from Stone's GDMA program
c
c
      subroutine match1 (i,record)
      use dma
      implicit none
      integer i
      character*240 record
c
c
c     store first multipole components on a line of GDMA output
c
      if (record(6:8) .eq. 'Q0 ') then
         read (record(13:23),*)  mp(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q00 ') then
         read (record(26:36),*)  mp(i)
      else if (record(20:23) .eq. 'Q10 ') then
         read (record(26:36),*)  dpz(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q11c') then
         read (record(26:36),*)  dpx(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q11s') then
         read (record(26:36),*)  dpy(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q20 ') then
         read (record(26:36),*)  q20(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q21c') then
         read (record(26:36),*)  q21c(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q21s') then
         read (record(26:36),*)  q21s(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q22c') then
         read (record(26:36),*)  q22c(i)
         call match2 (i,record)
      else if (record(20:23) .eq. 'Q22s') then
         read (record(26:36),*)  q22s(i)
         call match2 (i,record)
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine match2  --  match second value from GDMA output  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "match2" finds and stores the second multipole component found
c     on a line of output from Stone's GDMA program
c
c
      subroutine match2 (i,record)
      use dma
      implicit none
      integer i
      character*240 record
c
c
c     store second multipole component on a line of GDMA output
c
      if (record(29:31) .eq. 'Q1 ') then
         read (record(36:46),*)  dpz(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q11c') then
         read (record(45:55),*)  dpx(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q11s') then
         read (record(45:55),*)  dpy(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q21c') then
         read (record(45:55),*)  q21c(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q21s') then
         read (record(45:55),*)  q21s(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q22c') then
         read (record(45:55),*)  q22c(i)
         call match3 (i,record)
      else if (record(39:42) .eq. 'Q22s') then
         read (record(45:55),*)  q22s(i)
         call match3 (i,record)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine match3  --  match third value from GDMA output  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "match3" finds and stores the third multipole component found
c     on a line of output from Stone's GDMA program
c
c
      subroutine match3 (i,record)
      use dma
      implicit none
      integer i
      character*240 record
c
c
c     store third multipole component on a line of GDMA output
c
      if (record(52:54) .eq. 'Q2 ') then
         read (record(59:69),*)  q20(i)
      else if (record(58:61) .eq. 'Q11s') then
         read (record(64:74),*)  dpy(i)
      else if (record(58:61) .eq. 'Q21s') then
         read (record(64:74),*)  q21s(i)
      else if (record(58:61) .eq. 'Q22c') then
         read (record(64:74),*)  q22c(i)
      else if (record(58:61) .eq. 'Q22s') then
         read (record(64:74),*)  q22s(i)
      end if
      return
      end

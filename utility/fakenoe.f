c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program fakenoe  --  generate artificial NOE data set  ##
c     ##                                                         ##
c     #############################################################
c
c
      program fakenoe
      use sizes
      use atoms
      use iounit
      use pdb
      implicit none
      integer i,j,k
      integer map(maxatm)
      real*8 dist2,xi,yi,zi
      real*8 random,fraction,trial
      real*8 noemax,noemax2
      real*8 noemin,weight
      character*4 ixname,kxname
c
c
c     read coordinate sets from a Tinker file and a PDB file
c
      call initial
      call getxyz
      call getpdb
c
c     set defaults for NOE maximum, minimum and degeneracy factor
c
      noemin = 1.9d0
      noemax = 4.0d0
      weight = 1.0d0
c
c     get the maximum allowed NOE separation distance
c
      write (iout,10)
   10 format (/,' Enter Maximum Allowed NOE Distance in Angstroms',
     &           ' [4.0] :  ',$)
      read (input,20)  noemax
   20 format (f10.0)
      if (noemax .eq. 0.0d0)  noemax = 4.0d0
      noemax2 = noemax**2
c
c     get the random fraction of NOEs to be output
c
      fraction = 0.0d0
      write (iout,30)
   30 format (/,' Enter Percentage of Inter-Residue NOEs',
     &           ' to be Listed [100%] :  ',$)
      read (input,40)  fraction
   40 format (f10.0)
      fraction = fraction / 100.0d0
      if (fraction .eq. 0.0d0)  fraction = 1.0d0
c
c     map the PDB atoms onto the corresponding Tinker atoms
c
      do i = 1, npdb
         xi = xpdb(i)
         yi = ypdb(i)
         zi = zpdb(i)
         map(i) = 0
         do k = 1, n
            dist2 = (x(k)-xi)**2 + (y(k)-yi)**2 + (z(k)-zi)**2
            if (dist2 .lt. 0.01d0) then
               map(i) = k
               goto 50
            end if
         end do
   50    continue
      end do
c
c     get a list of the close interresidue hydrogen pairs
c
      do i = 1, npdb-1
         if (atmnam(i)(2:2) .eq. 'H') then
            xi = xpdb(i)
            yi = ypdb(i)
            zi = zpdb(i)
            do k = i+1, npdb
               if (atmnam(k)(2:2) .eq. 'H') then
                  dist2 = (xpdb(k)-xi)**2 + (ypdb(k)-yi)**2
     &                           + (zpdb(k)-zi)**2
                  if (dist2 .le. noemax2) then
                     if (resnum(i) .ne. resnum(k)) then
                        trial = random ()
                        if (trial .le. fraction) then
c
c     write the NOEs in the Tinker keyfile format
c
                           write (iout,60)  map(i),map(k),weight,noemin,
     &                                      noemax,resnum(i),atmnam(i),
     &                                      resnum(k),atmnam(k)
   60                      format ('restrain-distance',1x,2i7,
     &                                4x,f5.3,2f7.2,2x,'['i4,
     &                                1x,a4,1x,i4,1x,a4,']')
c
c     write the NOEs in the XPLOR distance restraint format
c
c                          call xplorit (atmnam(i),ixname)
c                          call xplorit (atmnam(k),kxname)
c                          write (iout,70)  resnum(i),ixname,
c    &                                      resnum(k),kxname,
c    &                                      (noemax+noemin)/2.0d0,
c    &                                      (noemax-noemin)/2.0d0,
c    &                                      (noemax-noemin)/2.0d0
c  70                      format ('assign (resid',i4,' and name ',
c    &                                a4,' ) (resid',i4,' and name ',
c    &                                a4,' )',3f6.2)
                        end if
c                    else
c                       write (iout,80)  map(i),map(k),resnum(i)
c  80                   format ('intraresidue',2x,3i7)
                     end if
                  end if
               end if
            end do
         end if
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine xplorit  --  convert PDB atom names to XPLOR  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine xplorit (pdb,xplor)
      implicit none
      integer i
      character*1 letter
      character*4 pdb,xplor
c
c
c     convert to the XPLOR atom naming format for hydrogens
c
      xplor = pdb
      if (xplor .eq. ' H  ')  xplor = ' HN '
      if (xplor .eq. '1H  ')  xplor = ' HT1'
      if (xplor .eq. '2H  ')  xplor = ' HT2'
      if (xplor .eq. '3H  ')  xplor = ' HT3'
      if (xplor(1:1).ge.'1' .and. xplor(1:1).le.'3') then
         letter = xplor(1:1)
         do i = 1, 3
            xplor(i:i) = xplor(i+1:i+1)
         end do
         xplor(4:4) = ' '
         do i = 1, 4
            if (xplor(i:i) .eq. ' ') then
               xplor(i:i) = letter
               goto 10
            end if
         end do
   10    continue
      end if
      return
      end

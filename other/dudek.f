c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program dudek  --  convert Dudek to Tinker multipoles  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "dudek" transforms multipoles in Mike Dudek's format
c     into the multipole format used by the Tinker package
c
c
      program dudek
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'iounit.i'
      include 'mpole.i'
      integer i,j,k,m,last,next
      integer iloc,iaxe,imdq
      integer freeunit,trimtext
      integer z1t,z1b,z2t,z2b,xt,xb
      integer zz1t(maxatm),zz1b(maxatm)
      integer zz2t(maxatm),zz2b(maxatm)
      integer xxt(maxatm),xxb(maxatm)
      real*8 a(3,3),a0(3,3),a1(3,3)
      character*1 letter
      character*3 z1tip,z1base,z2tip,z2base,xtip,xbase
      character*6 axename
      character*60 locfile,axefile,mdqfile
      character*80 record
c
c
c     get the user specified input structure
c
      write (iout,10)
   10 format (/,' Enter the Dudek Multipole File Name :  ',$)
      read (input,20)  filename
   20 format (a60)
c
c     remove any extension from the filename
c
      leng = trimtext (filename)
      last = leng
      do i = 1, leng
         letter = filename(i:i)
         if (letter .eq. '/')  last = leng
         if (letter .eq. ']')  last = leng
         if (letter .eq. ':')  last = leng
         if (letter .eq. '.')  last = i - 1
      end do
      leng = min(leng,last)
c
c     get Dudek's local axes and spherical harmonic multipoles
c
      call getpole
c
c     open the file containing the Dudek's local axes definitions
c
      locfile = filename(1:leng)
      call suffix (locfile,'loc')
      call version (locfile,'old')
      iloc = freeunit ()
      open (unit=iloc,file=locfile,status='old')
      rewind (unit=iloc)
c
c     read in the atoms defining Dudek's local axes
c
      read (iloc,30)
   30 format (/)
      do i = 1, n
         read (iloc,40)
   40    format ()
         read (iloc,50)  record
   50    format (a80)
         next = 1
         call getword (record,axename,next)
         call getword (record,z1tip,next)
         call getword (record,z1base,next)
         call getword (record,z2tip,next)
         call getword (record,z2base,next)
         do k = 1, n
            if (z1tip .eq. name(k))  zz1t(i) = k
            if (z1base .eq. name(k))  zz1b(i) = k
            if (z2tip .eq. name(k))  zz2t(i) = k
            if (z2base .eq. name(k))  zz2b(i) = k
         end do
         read (iloc,60)  record
   60    format (a80)
         next = 1
         call getword (record,axename,next)
         call getword (record,xtip,next)
         call getword (record,xbase,next)
         do k = 1, n
            if (xtip .eq. name(k))  xxt(i) = k
            if (xbase .eq. name(k))  xxb(i) = k
         end do
      end do
      close (iloc)
c
c     open the file containing the new local axes definitions
c
      axefile = filename(1:leng)
      call suffix (axefile,'axe')
      call version (axefile,'old')
      iaxe = freeunit ()
      open (unit=iaxe,file=axefile,status='old')
      rewind (unit=iaxe)
c
c     read in the atoms defining the new local axes
c
      read (iaxe,70)
   70 format (/)
      do i = 1, n
         read (iaxe,80)  record
   80    format (a80)
         read (record,*)  j,zaxis(i),xaxis(i)
      end do
      close (iaxe)
c
c     open the output file for the final rotated multipoles
c
      mdqfile = filename(1:leng)
      call suffix (mdqfile,'mdq')
      call version (mdqfile,'new')
      imdq = freeunit ()
      open (unit=imdq,file=mdqfile,status='new')
      rewind (unit=imdq)
c
c     get rotation matrix to convert old local frame to global
c
      do i = 1, n
         z1t = zz1t(i)
         z1b = zz1b(i)
         z2t = zz2t(i)
         z2b = zz2b(i)
         xt = xxt(i)
         xb = xxb(i)
         call rotmat (i,a0,z1t,z1b,z2t,z2b,xt,xb)
c
c     get rotation matrix to convert new local frame to global
c
         z1t = zaxis(i)
         z1b = i
         z2t = zaxis(i)
         z2b = i
         xt = xaxis(i)
         xb = i
         call rotmat (i,a1,z1t,z1b,z2t,z2b,xt,xb)
c
c     rotation to convert old to new local frame is the product
c
         do j = 1, 3
            do k = 1, 3
               a(j,k) = 0.0d0
               do m = 1, 3
                  a(j,k) = a(j,k) + a1(m,j)*a0(m,k)
               end do
            end do
         end do
c
c     rotate atomic multipoles from the old to new local frame
c
         call rotsite (i,a)
c
c     write out the Cartesian multipoles in the new frame
c
         write (imdq,90)  i,zaxis(i),xaxis(i),rpole(1,i)
   90    format ('multipole',5x,3i5,3x,f12.5)
         write (imdq,100)  (rpole(j,i),j=2,4)
  100    format (32x,3f12.5)
         write (imdq,110)  rpole(5,i)
  110    format (32x,f12.5)
         write (imdq,120)  (rpole(j,i),j=8,9)
  120    format (32x,2f12.5)
         write (imdq,130)  (rpole(j,i), j=11,13)
  130    format (32x,3f12.5)
      end do
      close (imdq)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine getpole  --  get atomic multipole parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "getpole" reads in the local coordinate axis definitions
c     and atomic multipole moments over spherical harmonics
c
c
      subroutine getpole
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'mpole.i'
c     include 'units.i'
      integer i,k,idma,next,freeunit
      real*8 q20,q21c,q21s,q22c,q22s
      character*60 dmafile
      character*80 record,string
c
c
c     open the file containing the multipole moments
c
      dmafile = filename(1:leng)
      call suffix (dmafile,'dma')
      call version (dmafile,'old')
      idma = freeunit ()
      open (unit=idma,file=dmafile,status='old')
      rewind (unit=idma)
c
c     read in the atom names and coordinates
c
      read (idma,10)
   10 format ()
      read (idma,20)  record
   20 format (a80)
      read (record,*)  n
      read (idma,30)
   30 format ()
      do i = 1, n
         read (idma,40)  record
   40    format (a80)
         next = 1
         call getword (record,name(i),next)
         string = record(next:80)
         read (string,*)  x(i),y(i),z(i)
      end do
c
c     read in the distributed multipole moments
c
      read (idma,50)
   50 format ()
      do i = 1, n
         read (idma,60)
   60    format ()
         read (idma,70)  record
   70    format (a80)
         read (record,*)  pole(1,i)
         read (idma,80)  record
   80    format (a80)
         read (record,*)  pole(4,i),pole(2,i),pole(3,i)
         read (idma,90)  record
   90    format (a80)
         read (record,*)  q20,q21c,q21s
         read (idma,100)  record
  100    format (a80)
         read (record,*)  q22c,q22s
         pole(13,i) = q20 / 3.0d0
         pole(7,i) = q21c * (0.8660254038d0/3.0d0)
         pole(10,i) = q21s * (0.8660254038d0/3.0d0)
         pole(5,i) = (q22c*0.8660254038d0 - q20*0.5d0) / 3.0d0
         pole(6,i) = q22s * (0.8660254038d0/3.0d0)
         pole(8,i) = pole(6,i)
         pole(9,i) = -(pole(5,i) + pole(13,i))
         pole(11,i) = pole(7,i)
         pole(12,i) = pole(10,i)
         do k = 1, 25
            read (idma,110)
  110       format ()
         end do
      end do
      close (idma)
c
c     convert the dipole and quadrupole moments to Angstroms
c
c     do i = 1, n
c        do k = 2, 4
c           pole(k,i) = pole(k,i) * bohr
c        end do
c        do k = 5, 13
c           pole(k,i) = pole(k,i) * bohr**2
c        end do
c     end do
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine rotmat  ##
c     ##                     ##
c     #########################
c
c
      subroutine rotmat (i,a,z1t,z1b,z2t,z2b,xt,xb)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,z1t,z1b,z2t,z2b,xt,xb
      real*8 dx,dy,dz
      real*8 r,dotxz,a(3,3)
      real*8 z1i,z1j,z1k
      real*8 z2i,z2j,z2k
      real*8 xi,xj,xk
c
c
c     find rotation matrix that converts local to global coordinates
c
      dx = x(z1t) - x(z1b)
      dy = y(z1t) - y(z1b)
      dz = z(z1t) - z(z1b)
      r = sqrt(dx*dx + dy*dy + dz*dz)
      z1i = dx / r
      z1j = dy / r
      z1k = dz / r
c
c
      dx = x(z2t) - x(z2b)
      dy = y(z2t) - y(z2b)
      dz = z(z2t) - z(z2b)
      r = sqrt(dx*dx + dy*dy + dz*dz)
      z2i = dx / r
      z2j = dy / r
      z2k = dz / r
c
c
      dx = z1i + z2i
      dy = z1j + z2j
      dz = z1k + z2k
      r = sqrt(dx*dx + dy*dy + dz*dz)
      a(1,3) = dx / r
      a(2,3) = dy / r
      a(3,3) = dz / r
c
c
      dx = x(xt) - x(xb)
      dy = y(xt) - y(xb)
      dz = z(xt) - z(xb)
      dotxz = dx*a(1,3) + dy*a(2,3) + dz*a(3,3)
      xi = dx - dotxz*a(1,3)
      xj = dy - dotxz*a(2,3)
      xk = dz - dotxz*a(3,3)
      r = sqrt(xi*xi + xj*xj + xk*xk)
      a(1,1) = xi / r
      a(2,1) = xj / r
      a(3,1) = xk / r
c
c
      a(1,2) = a(3,1)*a(2,3) - a(2,1)*a(3,3)
      a(2,2) = a(1,1)*a(3,3) - a(3,1)*a(1,3)
      a(3,2) = a(2,1)*a(1,3) - a(1,1)*a(2,3)
      return
      end

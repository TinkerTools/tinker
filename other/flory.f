c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program flory  --  get Flory polymer statistical averages  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "flory" reads in a set of structures and then computes the
c     average values for RMS superposition of Cartesian and dihedral
c     angle coordinates, end-to-end distance and radius of gyration
c
c
      program flory
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'zcoord.i'
      integer maxset
      parameter (maxset=1200)
      integer i,j,k
      integer ia,ib,ic,id
      integer ixyz,idat
      integer inn,inp
      integer ipn,ipp
      integer eps,lext
      integer next,start
      integer nset,npair
      integer nres,nrot
      integer trimtext,freeunit
      integer last,set(maxatm)
      real*8 rho,tau
      real*8 delta,sigma
      real*8 dist,rg,diff
      real8* rms,average
      real*8 stdev,pairs
      real*8 sum,sum2
      real*8 bond,size
      real*8 dihedral
      real*8 phi,psi
      real*8 qnn,qnp,qpn,qpp
      real*8 x1(maxatm)
      real*8 y1(maxatm)
      real*8 z1(maxatm)
      real*8 x2(maxatm)
      real*8 y2(maxatm)
      real*8 z2(maxatm)
      real*8 xs(maxgeo,maxset)
      real*8 ys(maxgeo,maxset)
      real*8 zs(maxgeo,maxset)
      real*8 rot(maxgeo,maxset)
      logical exist,query,oldverb
      logical alkane,peptide
      character*1 answer
      character*7 ext
      character*60 xyzfile,datafile
      character*80 record,string
c
c
c     get the name to use for the input coordinate files
c
      call initial
      call nextarg (filename,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Enter Cartesian Coordinate File Name :  ',$)
         read (input,20)  filename
   20    format (a60)
      end if
      call basefile (filename)
c
c     set the range of files to be used for chain statistics
c
      start = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=30,end=30)  start
         query = .false.
      end if
   30 continue
      if (query) then
         start = 1
         write (iout,40)
   40    format (/,' Number of the First File to Use [<CR>=1] :  ',$)
         read (input,50)  record
   50    format (a80)
         read (record,*,err=60,end=60)  start
   60    continue
      end if
c
c     molecule type is either polypeptide or polyethylene
c
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,70)
   70    format (/,' PolyPeptide or PolyEthylene (P or E) [P] :  ',$)
         read (input,80)  record
   80    format (a80)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'E') then
         alkane = .true.
         peptide = .false.
         bond = 0.54d0
         size = 1.53d0
      else
         peptide = .true.
         alkane = .false.
         bond = 1.34d0
         size = 3.80d0
      end if
c
c     get the alpha carbon coordinates and phi-psi angles
c
      nset = 0
      nres = 0
      exist = .true.
      dowhile (exist)
         lext = 3
         call numeral (nset+start,ext,lext)
         xyzfile = filename(1:leng)//'.'//ext(1:lext)
         inquire (file=xyzfile,exist=exist)
         if (exist) then
            ixyz = freeunit ()
            open (unit=ixyz,file=xyzfile,status='old')
            call readxyz (ixyz)
            close (unit=ixyz)
            nset = nset + 1
            if (nset .eq. 1) then
               if (peptide) then
                  do i = 1, n
                     if (name(i)(1:1) .eq. 'N') then
                        nres = nres + 1
                        set(nres) = i + 1
                     end if
                  end do
                  nrot = 2 * (nres-2)
               else if (alkane) then
                  do i = 1, n
                     if (name(i)(1:1) .eq. 'C') then
                        nres = nres + 1
                        set(nres) = i
                     end if
                     set(i) = i
                  end do
                  nrot = nres - 3
               end if
            end if
            do i = 1, nres
               xs(i,nset) = x(set(i))
               ys(i,nset) = y(set(i))
               zs(i,nset) = z(set(i))
            end do
            call makeint (0)
            if (peptide) then
               do i = 2, nres-1
                  ia = set(i-1) + 1
                  ib = set(i) - 1
                  ic = set(i)
                  id = set(i) + 1
                  rot(2*i-3,nset) = dihedral(ia,ib,ic,id)
                  ia = set(i) - 1
                  ib = set(i)
                  ic = set(i) + 1
                  id = set(i+1) - 1
                  rot(2*i-2,nset) = dihedral(ia,ib,ic,id)
               end do
            else if (alkane) then
               do i = 1, n-3
                  rot(i,nset) = ztors(i+3)
               end do
            end if
         end if
      end do
c
c     print number of structures and residues per structure
c
      write (iout,90)  nset,nres
   90 format (/,' Number of Structures :',i6,10x,
     &           'Residues per Structure :',i6)
c
c     set number of unique pairs of structures
c
      if (nset .lt. 3) then
         write (iout,100)
  100    format (/,' FLORY  --  Not Enough Structures for Statistics')
         call fatal
      end if
      npair = nset * (nset-1) / 2
      pairs = dble(npair)
      n = nres
c
c     compute the average end-to-end distance
c
      sum = 0.0d0
      sum2 = 0.0d0
      do i = 1, nset
         dist = (xs(1,i)-xs(nres,i))**2 + (ys(1,i)-ys(nres,i))**2
     &                      + (zs(1,i)-zs(nres,i))**2
         dist = sqrt(dist)
         sum = sum + dist
         sum2 = sum2 + dist**2
      end do
      average = sum / dble(nset)
      stdev = sqrt((dble(nset)*sum2-sum**2)/dble(nset*(nset-1)))
      eps = nint(100.0d0*stdev/average)
      write (iout,110)  average,stdev,eps
  110 format (/,' Average End-to-End Distance :',4x,f14.4,'  +/-',
     &           f8.4,2x,'('i3,'%)')
c
c     compute the average value of the radius of gyration
c
      sum = 0.0d0
      sum2 = 0.0d0
      do i = 1, nset
         do j = 1, nres
            x(j) = xs(j,i)
            y(j) = ys(j,i)
            z(j) = zs(j,i)
         end do
         call gyrate (rg)
         sum = sum + rg
         sum2 = sum2 + rg**2
      end do
      average = sum / dble(nset)
      stdev = sqrt((dble(nset)*sum2-sum**2)/dble(nset*(nset-1)))
      eps = nint(100.0d0*stdev/average)
      write (iout,120)  average,stdev,eps
  120 format (' Average Radius of Gyration :',5x,f14.4,'  +/-',
     &           f8.4,2x,'('i3,'%)')
c
c     compute the average root mean square superposition
c
      oldverb = verbose
      verbose = .false.
      sum = 0.0d0
      sum2 = 0.0d0
      do i = 1, nset-1
         do k = 1, nres
            x1(k) = xs(k,i)
            y1(k) = ys(k,i)
            z1(k) = zs(k,i)
         end do
         do j = i+1, nset
            do k = 1, nres
               x2(k) = xs(k,j)
               y2(k) = ys(k,j)
               z2(k) = zs(k,j)
            end do
            call impose (nres,x1,y1,z1,nres,x2,y2,z2,rms)
            rho = rms / (bond * sqrt(dble(nres-1)))
            sum = sum + rms
            sum2 = sum2 + rms**2
         end do
      end do
      verbose = oldverb
      average = sum / pairs
      stdev = sqrt((pairs*sum2-sum**2)/(pairs*(pairs-1.0d0)))
      eps = nint(100.0d0*stdev/average)
      write (iout,130)  average,stdev,eps
  130 format (' Average Pairwise RMSD :',10x,f14.4,'  +/-',
     &           f8.4,2x,'('i3,'%)')
      sum = sum  / (bond * sqrt(dble(nres-1)))
      sum2 = sum2 / (bond * sqrt(dble(nres-1)))**2
      average = sum / pairs
      stdev = sqrt((pairs*sum2-sum**2)/(pairs*(pairs-1.0d0)))
      eps = nint(100.0d0*stdev/average)
      write (iout,140)  average,stdev,eps
  140 format (/,' The Value of RHO :',15x,f14.4,'  +/-',
     &           f8.4,2x,'('i3,'%)')
c
c     compute the rms torsional angle variation
c
      sum = 0.0d0
      sum2 = 0.0d0
      do i = 1, nset-1
         do j = i+1, nset
            rms = 0.0d0
            do k = 1, nrot
               diff = abs(rot(k,j)-rot(k,i))
               if (diff .gt. 180.0d0)  diff = 360.0d0 - diff
               rms = rms + diff**2
            end do
            rms = sqrt(rms/dble(nrot))
            delta = rms / 103.9d0
            sum = sum + delta
            sum2 = sum2 + delta**2
         end do
      end do
      average = sum / pairs
      stdev = sqrt((pairs*sum2-sum**2)/(pairs*(pairs-1.0d0)))
      eps = nint(100.0d0*stdev/average)
      write (iout,150)  average,stdev,eps
  150 format (' The Value of DELTA :',13x,f14.4,'  +/-',
     &           f8.4,2x,'('i3,'%)')
c
c     compute the average squared end-to-end distance
c
      sum = 0.0d0
      sum2 = 0.0d0
      do i = 1, nset
         dist = (xs(1,i)-xs(nres,i))**2 + (ys(1,i)-ys(nres,i))**2
     &                      + (zs(1,i)-zs(nres,i))**2
         sigma = dist / (size**2 * dble(nres-1))
         sum = sum + sigma
         sum2 = sum2 + sigma**2
      end do
      average = sum / dble(nset)
      stdev = sqrt((dble(nset)*sum2-sum**2)/dble(nset*(nset-1)))
      eps = nint(100.0d0*stdev/average)
      write (iout,160)  average,stdev,eps
  160 format (' The Value of SIGMA^2 :',11x,f14.4,'  +/-',
     &           f8.4,2x,'('i3,'%)')
c
c     compute the average squared radius of gyration
c
      sum = 0.0d0
      sum2 = 0.0d0
      do i = 1, nset
         do j = 1, nres
            x(j) = xs(j,i)
            y(j) = ys(j,i)
            z(j) = zs(j,i)
         end do
         call gyrate (rg)
         tau = (6.0d0 * rg**2) / (size**2 * dble(nres-1))
         sum = sum + tau
         sum2 = sum2 + tau**2
      end do
      average = sum / dble(nset)
      stdev = sqrt((dble(nset)*sum2-sum**2)/dble(nset*(nset-1)))
      eps = nint(100.0d0*stdev/average)
      write (iout,170)  average,stdev,eps
  170 format (' The Value of TAU^2 :',13x,f14.4,'  +/-',
     &           f8.4,2x,'('i3,'%)')
c
c     compute peptide phi-psi quadrant occupancy and write angles
c
      if (peptide) then
         idat = freeunit ()
         datafile = filename(1:leng)//'.dat'
         call version (datafile,'new')
         open (unit=idat,file=datafile,status='new')
         inn = 0
         inp = 0
         ipn = 0
         ipp = 0
         do i = 1, nset
            do j = 2, nres-1
               phi = rot(2*j-3,i)
               psi = rot(2*j-2,i)
               if (phi.lt.0.0d0 .and. psi.lt.0.0d0)  inn = inn + 1
               if (phi.lt.0.0d0 .and. psi.gt.0.0d0)  inp = inp + 1
               if (phi.gt.0.0d0 .and. psi.lt.0.0d0)  ipn = ipn + 1
               if (phi.gt.0.0d0 .and. psi.gt.0.0d0)  ipp = ipp + 1
               write (idat,180)  i,j,phi,psi
  180          format (2i6,2f8.1)
            end do
         end do
         close (unit=idat)
         sum = dble(inn + inp + ipn + ipp)
         qnn = 100.0d0 * dble(inn) / sum
         qnp = 100.0d0 * dble(inp) / sum
         qpn = 100.0d0 * dble(ipn) / sum
         qpp = 100.0d0 * dble(ipp) / sum
         write (iout,190)  qnn,qnp,qpn,qpp
  190    format (/,' Phi-Psi Quadrant Occupancy :',8x,4(f7.2,'%'))
      end if
c
c     perform any final tasks before program exit
c
      call final
      end

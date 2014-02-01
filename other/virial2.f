c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2001 by Pengyu Ren and Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program virial2  --  calculate second virial coefficient  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "virial2" computes the second virial coefficient; use trajectory
c     frame as the first input file, and dimer structure as the second
c     input file; current version is specific for water
c
c
      program virial2
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'inter.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'units.i'
      integer maxbin
      parameter (maxbin=1000)
      integer i,j,k,m,ii
      integer freeunit,trimtext
      integer rint,nummol,molatm
      integer hist(maxbin)
      real*8 width,weigh
      real*8 conver,factor
      real*8 temp,kt,molmas
      real*8 bt2,btcf,btct
      real*8 xcm,ycm,zcm
      real*8 xicm,yicm,zicm
      real*8 xjcm,yjcm,zjcm
      real*8 rdist,rmax,rbin
      real*8 xdiff,ydiff,zdiff
      real*8 energy,hbar,bolt
      real*8 bterm,eterm1
      real*8 fterm1,tterm1
      real*8 trqx,trqy,trqz
      real*8 xx,xy,xz,yy,yz,zz
      real*8 xterm,yterm,zterm
      real*8 dist2,fsum,tsum
      real*8 moment(3),rvec(3)
      real*8 fcmi(3),fcmj(3)
      real*8 rcmx(3),rcmy(3),rcmz(3)
      real*8 vec(3,3),tensor(3,3)
      real*8 a(3,3),fatm(3,6)
      real*8 eterm(maxbin)
      real*8 fterm(maxbin)
      real*8 tterm(maxbin)
      real*8, allocatable :: xatm(:)
      real*8, allocatable :: yatm(:)
      real*8, allocatable :: zatm(:)
      real*8, allocatable :: pmi(:,:)
      logical exist,include
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     get the temperature value and width of the first bin
c
      temp = -1.0d0
      do while (temp .lt. 0.0d0)
         write (iout,10)
   10    format (/,' Enter the Desired Temperature in Degrees',
     &              ' K [298] :  ',$)
         read (input,20,err=30)  temp
   20    format (f20.0)
         if (temp .le. 0.0d0)  temp = 298.0d0
   30    continue
      end do
      width = -1.0d0
      do while (width .lt. 0.0d0)
         write (iout,40)
   40    format (/,' Enter the Initial Distance Bin Width in Ang',
     &              ' [1.0] :  ',$)
         read (input,50,err=60)  width
   50    format (f20.0)
         if (width .le. 0.0d0)  width = 1.0d0
   60    continue
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xatm(n))
      allocate (yatm(n))
      allocate (zatm(n))
      allocate (pmi(3,n))
c
c     initialize quantities needed for the calculation
c
      rmax = 20.0d0
      kt = gasconst * temp
      hbar = planck / (2.0d0*pi)
      conver = avogadro * 1.0d-27
      do j = 1, 3
         do i = 1, 6
            fatm(j,i) = 0.0d0
         end do
         fcmi(j) = 0.0d0
         fcmj(j) = 0.0d0
      end do
      do i = 1, maxbin
         hist(i) = 0
         eterm(i) = 0.0d0
         fterm(i) = 0.0d0
         tterm(i) = 0.0d0
      end do
      do i = 1, n
         do j = 1, 3
            pmi(j,i) = 0.0d0
         end do
      end do
      do i = 1, n
         xatm(i) = x(i)
         yatm(i) = y(i)
         zatm(i) = z(i)
      end do
c
c     initialize quantities needed for the calculation
c
      do j = 1, 3
         do i = 1, 6
            fatm(j,i) = 0.0d0
         end do
         fcmi(j) = 0.0d0
         fcmj(j) = 0.0d0
      end do
      do i = 1, maxbin
         hist(i) = 0
         eterm(i) = 0.0d0
         fterm(i) = 0.0d0
         tterm(i) = 0.0d0
      end do
      do i = 1, nmol
         do j = 1, 3
            pmi(j,i) = 0.0d0
         end do
      end do
      molatm = n / nmol
      molmas = 0.0d0
      do i = 1, molatm
         molmas = molmas + mass(i)
      end do
      factor = hbar**2 * avogadro**2 * 1.0d23
     &             / (24.0d0*(kt)**3*molmas) / 4182.0d0
c
c     compute and then diagonalize the inertial tensor
c
      do j = 1, nmol
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         ii = 1
         do i = 3*j-2, j*3
            weigh = mass(ii)
            xcm = xcm + xatm(i)*weigh
            ycm = ycm + yatm(i)*weigh
            zcm = zcm + zatm(i)*weigh
            ii = ii + 1
         end do
         xcm = xcm / molmas
         ycm = ycm / molmas
         zcm = zcm / molmas
         xx = 0.0d0
         xy = 0.0d0
         xz = 0.0d0
         yy = 0.0d0
         yz = 0.0d0
         zz = 0.0d0
         ii = 1
         do i = 3*j-2, 3*j
            weigh = mass(ii)
            xterm = xatm(i) - xcm
            yterm = yatm(i) - ycm
            zterm = zatm(i) - zcm
            xx = xx + xterm*xterm*weigh
            xy = xy + xterm*yterm*weigh
            xz = xz + xterm*zterm*weigh
            yy = yy + yterm*yterm*weigh
            yz = yz + yterm*zterm*weigh
            zz = zz + zterm*zterm*weigh
            ii = ii + 1
         end do
         tensor(1,1) = yy + zz
         tensor(2,1) = -xy
         tensor(3,1) = -xz
         tensor(1,2) = -xy
         tensor(2,2) = xx + zz
         tensor(3,2) = -yz
         tensor(1,3) = -xz
         tensor(2,3) = -yz
         tensor(3,3) = xx + yy
         call jacobi (3,tensor,moment,vec)
         do i = 1, 3
            pmi(i,j) = moment(i)
         end do
      end do
c
c     read a dimer structure to setup the BT(2) calculation
c
      nummol = nmol
      call getxyz
      call mechanic
c
c     compute the BT(2) integrand over distance-based bins
c
      write (iout,70)  nummol
   70 format (/,5x,'Number of Molecules :  ',i8)
      write (iout,80)  temp
   80 format (5x,'Temperature (K) :  ',f12.2)
      write (iout,90)
   90 format (/,7x,'R',6x,'Eterm',6x,'Fterm',6x,'Tterm',
     &          6x,'BT(2)',7x,'BTCf',7x,'BTCt')
      bt2 = 0.0d0
      btcf = 0.0d0
      btct = 0.0d0
      rdist = width
      rbin = width
      do rint = 1, maxbin
         do i = 1, nummol-1
            do j = i+1, nummol
               xicm = 0.0d0
               yicm = 0.0d0
               zicm = 0.0d0
               do m = 1, 3
                  k = (i-1)*3 + m
                  x(m) = xatm(k)
                  y(m) = yatm(k)
                  z(m) = zatm(k)
                  xicm = xicm + x(m)*mass(m)/molmas
                  yicm = yicm + y(m)*mass(m)/molmas
                  zicm = zicm + z(m)*mass(m)/molmas
               end do
               xjcm = 0.0d0
               yjcm = 0.0d0
               zjcm = 0.0d0
               do m = 4, 6
                  k = (j-2)*3 + m
                  x(m) = xatm(k)
                  y(m) = yatm(k)
                  z(m) = zatm(k)
                  xjcm = xjcm + x(m)*mass(m)/molmas
                  yjcm = yjcm + y(m)*mass(m)/molmas
                  zjcm = zjcm + z(m)*mass(m)/molmas
               end do
               xdiff = rdist + xicm - xjcm
               ydiff = yicm - yjcm
               zdiff = zicm - zjcm
c              call ranvec (rvec)
c              xdiff = rdist*rvec(1) + xicm - xjcm
c              ydiff = rdist*rvec(2) + yicm - yjcm
c              zdiff = rdist*rvec(3) + zicm - zjcm
               do m = 4, 6
                  x(m) = x(m) + xdiff
                  y(m) = y(m) + ydiff
                  z(m) = z(m) + zdiff
               end do
               include = .true.
               do k = 1, 3
                  do m = 4, 6
                     dist2 = (x(k)-x(m))**2 + (y(k)-y(m))**2
     &                          + (z(k)-z(m))**2
                     if (dist2 .lt. 0.25d0)  include = .false.
                  end do
               end do
               if (include) then
                  call gradient (energy,fatm)
                  hist(rint) = hist(rint) + 1
                  bolt = exp(-einter/kt)
                  eterm(rint) = eterm(rint) + bolt
                  trqx = 0.0d0
                  trqy = 0.0d0
                  trqz = 0.0d0
                  do ii = 1, 3
                     fcmi(ii) = fatm(ii,1) + fatm(ii,2) + fatm(ii,3)
c                    fcmj(ii) = fatm(ii,4) + fatm(ii,5) + fatm(ii,6)
                     rcmx(ii) = x(ii) - xicm
                     rcmy(ii) = y(ii) - yicm
                     rcmz(ii) = z(ii) - zicm
                     trqx = trqx+rcmy(ii)*fatm(3,ii)-rcmz(ii)*fatm(2,ii)
                     trqy = trqy+rcmz(ii)*fatm(1,ii)-rcmx(ii)*fatm(3,ii)
                     trqz = trqz+rcmx(ii)*fatm(2,ii)-rcmy(ii)*fatm(1,ii)
                  end do
                  fsum = bolt * (fcmi(1)**2+fcmi(2)**2+fcmi(3)**2)
                  tsum = bolt * (trqx**2/pmi(1,i)+trqy**2/pmi(2,i)
     &                                  +trqz**2/pmi(3,i))
                  fterm(rint) = fterm(rint) + fsum
                  tterm(rint) = tterm(rint) + tsum
               end if
            end do
         end do
         if (hist(rint) .gt. 0) then
            eterm1 = eterm(rint)/dble(hist(rint)) - 1.0d0
            fterm1 = fterm(rint)/dble(hist(rint)) * factor
            tterm1 = tterm(rint)/dble(hist(rint)) * factor * molmas
         else
            eterm(i) = 0.0d0
         end if
         bterm = 2.0d0 * pi * conver * rbin * rdist * rdist
         bt2 = bt2 - bterm*eterm1
         btcf = btcf + 2.0d0*bterm*fterm1
         btct = btct + 2.0d0*bterm*tterm1
         write (iout,100)  rdist,eterm1,fterm1,tterm1,
     &                     bt2+btcf+btct,btcf,btct
  100    format (f8.3,6f11.5)
         rbin = width / (((rdist-2.8d0)**3+10.0d0)
     &                    * exp(-0.15d0*(rdist-2.8d0)**2)+3.0d0)
         rdist = rdist + rbin
         if (rdist .ge. rmax)  goto 110
      end do
  110 continue
c
c     perform deallocation of some local arrays
c
      deallocate (xatm)
      deallocate (yatm)
      deallocate (zatm)
      deallocate (pmi)
c
c     perform any final tasks before program exit
c
      call final
      end

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
c     "virial2" computes the classical second virial coefficient with
c     quantum corrections; provide an MD trajectory frame as the first
c     input file, and any dimer structure as the second input file
c
c     note this version is specific for water, but caould be easily
c     generalized to other molecules
c
c
      program virial2
      use sizes
      use atomid
      use atoms
      use boxes
      use inter
      use iounit
      use math
      use molcul
      use units
      implicit none
      integer maxbin
      parameter (maxbin=1000)
      integer i,j,k,m,ii
      integer freeunit,trimtext
      integer rint,nummol,molatm
      integer hist(maxbin)
      real*8 kelvin,kt
      real*8 weigh,molmas
      real*8 conver,factor
      real*8 bt2,btcf,btct
      real*8 xcm,ycm,zcm
      real*8 xicm,yicm,zicm
      real*8 xjcm,yjcm,zjcm
      real*8 rdist,rbin
      real*8 rinit,rmax
      real*8 xdiff,ydiff,zdiff
      real*8 e,hbar,bolt,bterm
      real*8 etot,ftot,ttot
      real*8 trqx,trqy,trqz
      real*8 xx,xy,xz,yy,yz,zz
      real*8 xterm,yterm,zterm
      real*8 dist2,fsum,tsum
      real*8 moment(3),rvec(3)
      real*8 fcmi(3),fcmj(3)
      real*8 rcmx(3),rcmy(3),rcmz(3)
      real*8 vec(3,3),tensor(3,3)
      real*8 a(3,3),g(3,6)
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
c     get temperature value to be used for the second virial
c
      kelvin = -1.0d0
      do while (kelvin .lt. 0.0d0)
         write (iout,10)
   10    format (/,' Enter the Desired Temperature in Degrees',
     &              ' K [298] :  ',$)
         read (input,20,err=30)  kelvin
   20    format (f20.0)
         if (kelvin .le. 0.0d0)  kelvin = 298.0d0
   30    continue
      end do
c
c     set needed constants and unit conversion factors
c
      nummol = nmol
      molatm = n / nmol
      molmas = molmass(1)
      rinit = 1.0d0
      rmax = 20.0d0
      kt = gasconst * kelvin
      hbar = planck / (2.0d0*pi)
      conver = avogadro * 1.0d-24
      factor = hbar**2 * avogadro**2 * 1.0d23
     &            / (24.0d0*(kt)**3) / 4182.0d0
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
      do i = 1, n
         xatm(i) = x(i)
         yatm(i) = y(i)
         zatm(i) = z(i)
      end do
      do i = 1, n
         do j = 1, 3
            pmi(j,i) = 0.0d0
         end do
      end do
      do j = 1, 3
         do i = 1, 6
            g(j,i) = 0.0d0
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
c
c     compute and then diagonalize the inertial tensor
c
      do i = 1, nmol
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         do j = imol(1,i), imol(2,i)
            k = kmol(j)
            weigh = mass(k)
            xcm = xcm + xatm(k)*weigh
            ycm = ycm + yatm(k)*weigh
            zcm = zcm + zatm(k)*weigh
         end do
         xcm = xcm / molmass(i)
         ycm = ycm / molmass(i)
         zcm = zcm / molmass(i)
         xx = 0.0d0
         xy = 0.0d0
         xz = 0.0d0
         yy = 0.0d0
         yz = 0.0d0
         zz = 0.0d0
         ii = 1
         do j = imol(1,i), imol(2,i)
            k = kmol(j)
            weigh = mass(k)
            xterm = xatm(k) - xcm
            yterm = yatm(k) - ycm
            zterm = zatm(k) - zcm
            xx = xx + xterm*xterm*weigh
            xy = xy + xterm*yterm*weigh
            xz = xz + xterm*zterm*weigh
            yy = yy + yterm*yterm*weigh
            yz = yz + yterm*zterm*weigh
            zz = zz + zterm*zterm*weigh
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
         do j = 1, 3
            pmi(j,i) = moment(j)
         end do
      end do
c
c     reset box size to disable periodic boundary conditions
c
      xbox = 0.0d0
      ybox = 0.0d0
      zbox = 0.0d0
c
c     read any dimer structure to reset mechanics calculation
c
      call getxyz
      call mechanic
c
c     print initial information prior to numerical integration
c
      write (iout,40)  nummol
   40 format (/,' Number of Molecules :  ',i12)
      write (iout,50)  nummol*(nummol-1)/2
   50 format (' Dimer Orientations :  ',i13)
      write (iout,60)  kelvin
   60 format (' Temperature (K) :  ',f16.2)
      write (iout,70)
   70 format (/,' Cummulative Classical and Quantum Corrected',
     &           ' B(T) Values :')
      write (iout,80)
   80 format (/,7x,'R',5x,'Bterm',5x,'Fterm',5x,'Tterm',4x,
     &           'Bcl(T)',6x,'dBtr',5x,'dBrot',6x,'B(T)',/)
c
c     compute the B(T) integrands over distance-based bins
c
      bt2 = 0.0d0
      btcf = 0.0d0
      btct = 0.0d0
      rdist = rinit
      rbin = rinit
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
               call ranvec (rvec)
               xdiff = rdist*rvec(1) + xicm - xjcm
               ydiff = rdist*rvec(2) + yicm - yjcm
               zdiff = rdist*rvec(3) + zicm - zjcm
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
                     if (dist2 .lt. 1.0d0)  include = .false.
                  end do
               end do
               if (include) then
                  call gradient (e,g)
                  hist(rint) = hist(rint) + 1
                  bolt = exp(-einter/kt)
                  eterm(rint) = eterm(rint) + bolt
                  trqx = 0.0d0
                  trqy = 0.0d0
                  trqz = 0.0d0
                  do ii = 1, 3
                     fcmi(ii) = -g(ii,1) - g(ii,2) - g(ii,3)
                     fcmj(ii) = -g(ii,4) - g(ii,5) - g(ii,6)
                     rcmx(ii) = x(ii) - xicm
                     rcmy(ii) = y(ii) - yicm
                     rcmz(ii) = z(ii) - zicm
                     trqx = trqx - rcmy(ii)*g(3,ii) + rcmz(ii)*g(2,ii)
                     trqy = trqy - rcmz(ii)*g(1,ii) + rcmx(ii)*g(3,ii)
                     trqz = trqz - rcmx(ii)*g(2,ii) + rcmy(ii)*g(1,ii)
                  end do
                  fsum = bolt * (fcmi(1)**2+fcmi(2)**2+fcmi(3)**2)
                  tsum = bolt * (trqx**2/pmi(1,i)+trqy**2/pmi(2,i)
     &                                   +trqz**2/pmi(3,i))
                  fterm(rint) = fterm(rint) + fsum
                  tterm(rint) = tterm(rint) + tsum
               end if
            end do
         end do
         if (hist(rint) .gt. 0) then
            etot = eterm(rint)/dble(hist(rint)) - 1.0d0
            ftot = fterm(rint)/dble(hist(rint)) * factor / molmas
            ttot = tterm(rint)/dble(hist(rint)) * factor
         else
            etot = -1.0d0
            ftot = 0.0d0
            ttot = 0.0d0
         end if
         bterm = 4.0d0 * pi * conver * rbin * rdist**2
         bt2 = bt2 - 0.5d0*bterm*etot
         btcf = btcf + bterm*ftot
         btct = btct + bterm*ttot
         write (iout,90)  rdist,etot,ftot,ttot,
     &                     bt2,btcf,btct,bt2+btcf+btct
   90    format (f8.3,3f10.4,4f10.1)
         rbin = rinit / (((rdist-2.8d0)**3+10.0d0)
     &                    * exp(-0.15d0*(rdist-2.8d0)**2)+3.0d0)
         rdist = rdist + rbin
         if (rdist .ge. rmax)  goto 100
      end do
  100 continue
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

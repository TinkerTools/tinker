cccccccccccccccccccccccccccccccccccccccccccccccc
c                                              c
c   Copute second virial coefficient           c
c   starting bin width:  0.75~1                c
c   Need water trajectories frames (at least 3)c
c   need dimer xyz and key file at second coord prompt
c   don't input box file name untill first prompt  c
cccccccccccccccccccccccccccccccccccccccccccccccc
      program virial2 
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'boxes.i'
      include 'couple.i'
      include 'deriv.i'
      include 'energi.i'
      include 'files.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'mpole.i'
      include 'polar.i'
      include 'shunt.i'
      include 'units.i'
      include 'usage.i'
      include 'atmtyp.i'
      include 'iounit.i'
      include 'math.i'

      integer maxbin, maxatoms
      parameter(maxbin = 2000)
      parameter(maxatoms = 500000)
      integer hist(maxbin),rint,binno
      integer trimtext,leng2
      integer freeunit
      integer start,stop,step
      integer iarc,ixyz,last
      integer file_num,run
      integer lext
      integer i,j,m,k,molatm,totmol,ii,kk
      real*8 rdist2,rdist,width,energy,conver,conver2
      real*8 tempa,kt,molmas,bt2, btcf,btct, h
      real*8 xatm(maxatoms),yatm(maxatoms),zatm(maxatoms)
      real*8 eu(maxbin),edim,fatm(3,6),fcmi(3),fcmj(3)
      real*8 fterm(maxbin),tterm(maxbin)
      real*8 fterm1,tterm1,eterm1
      real*8 xicm,xjcm,yicm,yjcm,zicm,zjcm,tt,ttt
      real*8 xdif,ydif,zdif,bolt
      real*8 trqx, trqy,trqz,rcmx(3),rcmy(3),rcmz(3)

      real*8 weigh,total,dot
      real*8 xcm,ycm,zcm
      real*8 xx,xy,xz,yy,yz,zz
      real*8 xterm,yterm,zterm
      real*8 pmi(3,maxatoms/3),moment(3),vec(3,3)
      real*8 work1(3),work2(3)
      real*8 tensor(3,3),a(3,3)

      character*1 letter
      character*80 record,coordfile,indfile,output,tempas
      character*10 ext
      logical exist
      call initial
      call getxyz
      call mechanic    
 

      leng = trimtext (filename)
      iarc = freeunit ()
      last = leng
      do i = 1, leng
         letter = filename(i:i)
         if (letter .eq. '/')  last = leng
         if (letter .eq. ']')  last = leng
         if (letter .eq. ':')  last = leng
         if (letter .eq. '.')  last = i - 1
      end do
      leng = min(leng,last)
      start = 0
      stop = 0
      step = 0
      write (iout,70)
   70 format (/,' Numbers of First & Last File and Step',
     &              ' Increment :  ',$)
      read (input,80)  record
   80 format (a80)
      read (record,*,err=90,end=90)  start,stop,step
   90 continue
      if (stop .eq. 0)  stop = start
      if (step .eq. 0)  step = 1
  110 continue
      write (iout,116)
  116 format(/,' Temperature(K) and Distance bin width: ', $)
      read (input,80)  record
      read (record,*,err=90,end=90)  tempa, width
c
c     cycle over the user specified coordinate files
c
      kt = gasconst * tempa
      conver = avogadro * 1.0d-27
      file_num = start
      run = 0
      h = 6.6262d-34/(2.0d0*pi)
      do i = 1,maxatoms 
         xatm(i) =0.0d0
         yatm(i) =0.0d0
         zatm(i) =0.0d0
      end do
      dowhile (file_num .le. stop)
c         tt = mod(file_num, 100)
c         if (tt  .eq. 0) print *, "Processing ",file_num
         ixyz = freeunit ()
         lext = 3
         call numeral (file_num,ext,lext)
         coordfile = filename(1:leng)//'.'//ext(1:lext)
         call version (coordfile,'old')
         inquire (file=coordfile,exist=exist)
         if (.not.exist .and. file_num.lt.100) then
            lext = 2
            call numeral (file_num,ext,lext)
            coordfile = filename(1:leng)//'.'//ext(1:lext)
            call version (coordfile,'old')
            inquire (file=coordfile,exist=exist)
         end if
         if (.not.exist .and. file_num.lt.10) then
            lext = 1
            call numeral (file_num,ext,lext)
            coordfile = filename(1:leng)//'.'//ext(1:lext)
            call version (coordfile,'old')
            inquire (file=coordfile,exist=exist)
         end if
c    read in the coordinates and induced dipole from each frame       
         if (exist) then
            leng2 = trimtext (coordfile)
            open (unit=ixyz,file=coordfile,status='old')
            call readxyz (ixyz)
            close (unit=ixyz)
         else 
            go to 111
         end if
         do i = 1, n
            j=run*n+i
            xatm(j) = x(i)
            yatm(j) = y(i)
            zatm(j) = z(i)
         end do
         run=run + 1
         file_num = file_num + step
 111     continue
      end do
      print *, " "
      print *, " Finishing reading archive" 
      print *, " "
c
c  computer principle moemnets of inertia for each mol
c
c      binno = int(9.0/width+6.0/width/10.0)
      molmas = 0.0d0 
      edim = 0.0d0
      do j = 1, 3
         do i = 1, 6 
            fatm(j,i) = 0.0d0
         end do
         fcmi(j) = 0.0d0
         fcmj(j) = 0.0d0
      end do
      do i = 1, maxbin 
        eu(i) = 0.0d0
        hist(i) = 0
        fterm(i) = 0.0d0
        tterm(i) = 0.0d0
      end do
      do i = 1, n
        x(i) = 0.0d0
        y(i) = 0.0d0
        z(i) = 0.0d0
      end do
      molatm = n/nmol
      totmol = run * nmol
      do i = 1, molatm
        molmas = molmas + mass(i)
      end do
      conver2= h**2*avogadro**2*1.0d23/(24.d0*(kt)**3*molmas)/4182.0
c      print *, 'Conver2 ', conver2
      do i = 1, maxatoms/3 
         do j = 1, 3
           pmi(j,i) = 0.0d0
        end do
      end do
      do m = 1, totmol
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         ii = 1
         do i = 3*m-2, m*3
            weigh = mass(ii)
            xcm = xcm + xatm(i)*weigh
            ycm = ycm + yatm(i)*weigh
            zcm = zcm + zatm(i)*weigh
            ii = ii + 1
         end do
         xcm = xcm / molmas
         ycm = ycm / molmas
         zcm = zcm / molmas
c
c     compute and then diagonalize the inertial tensor
c
         xx = 0.0d0
         xy = 0.0d0
         xz = 0.0d0
         yy = 0.0d0
         yz = 0.0d0
         zz = 0.0d0
         ii = 1
         do i = 3*m-2, 3*m
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
         call jacobi (3,3,tensor,moment,vec,work1,work2)
         do kk = 1,3
            pmi(kk,m) = moment(kk)
         end do
      end do
      do i = 1, totmol 
         do j = 1, 3
        end do
      end do

c
c  compute BT(2) integrand
c
c
c  Swith to dimer
c
      call initial
      call getxyz
      call mechanic    
      write(tempas,'(f8.2)') tempa
      output=tempas(1:5)//".out"
      ixyz = freeunit ()
      inquire (file=output,exist=exist)
      if (exist) then
         open (unit=ixyz,file=output,status='old')
         dowhile (.true.)
            read (ixyz,120,err=130,end=130)
  120       format ()
         end do
  130    continue
      else
         open (unit=ixyz,file=output,status='new')
      end if
      write (ixyz,16) totmol
      write (ixyz,18) start,stop,step,tempa
      write (ixyz,17) 
  16  format (/,"#Total number of molecules: ",i10,$)
  17  format (/,"    r      Eterm      Fterm      Tterm      BT(2)"
     &         ,"      BTCf       BTCt ")   
  18  format (/,"     Start Stop  Step  Temp",/,3i6,f8.2)

      bt2 = 0.0d0
      btcf= 0.0d0
      btct= 0.0d0
      rdist = width
      ttt = width
      do rint = 1, maxbin 
         print *, "Processing r = ", rdist
         do i = 1, totmol
c            write(*,*) "moli ",i
            do j = i+1, totmol
               xicm = 0.0d0
               yicm = 0.0d0
               zicm = 0.0d0
               do m = 1,3
                  k = (i-1)*3 + m
                  x(m)= xatm(k)
                  y(m)= yatm(k)
                  z(m)= zatm(k)
                  xicm = xicm + x(m)*mass(m)/molmas
                  yicm = yicm + y(m)*mass(m)/molmas
                  zicm = zicm + z(m)*mass(m)/molmas
c               print * , m,"mass", mass(m)
               end do
               xjcm = 0.0d0
               yjcm = 0.0d0
               zjcm = 0.0d0
               do m = 4, 6 
                  k = (j-2)*3 + m
                  x(m)= xatm(k)
                  y(m)= yatm(k)
                  z(m)= zatm(k)
                  xjcm = xjcm + x(m)*mass(m)/molmas
                  yjcm = yjcm + y(m)*mass(m)/molmas
                  zjcm = zjcm + z(m)*mass(m)/molmas
c                  print * , m,"mass", mass(m)
               end do
               xdif = rdist + xicm - xjcm
               ydif = yicm - yjcm
               zdif = zicm - zjcm
               do m = 4, 6
                  x(m) = x(m) + xdif
                  y(m) = y(m) + ydif
                  z(m) = z(m) + zdif
               end do
c            do ii = 1, 6
c               print *, x(ii),' ',y(ii), '  ',z(ii)
c            end do 
               call gradient(edim, fatm)
c               write(*,*) "calling gradient"
               if ( einter.gt.-10.0d0) then
                  hist(rint) = hist(rint) + 1
                  bolt = exp(-einter/kt)
                  eu(rint) = eu(rint) + bolt 
                  trqx= 0.0d0
                  trqy= 0.0d0
                  trqz= 0.0d0
                  do ii = 1,3
                     fcmi(ii) = fatm(ii,1) + fatm(ii,2) + fatm(ii,3)
c                     fcmj(ii) = fatm(ii,4) + fatm(ii,5) + fatm(ii,6)
                     rcmx(ii) = x(ii) - xicm
                     rcmy(ii) = y(ii) - yicm
                     rcmz(ii) = z(ii) - zicm
                     trqx = trqx+rcmy(ii)*fatm(3,ii)-rcmz(ii)*fatm(2,ii)
                     trqy = trqy+rcmz(ii)*fatm(1,ii)-rcmx(ii)*fatm(3,ii)
                     trqz = trqz+rcmx(ii)*fatm(2,ii)-rcmy(ii)*fatm(1,ii)
                  end do
                  fterm(rint) = fterm(rint) + bolt * (fcmi(1)**2
     &                     + fcmi(2)**2 + fcmi(3)**2)
                  tterm(rint)= tterm(rint) + bolt*(trqx**2/pmi(1,i)
     &                    + trqy**2/pmi(2,i) + trqz**2/pmi(3,i))
               end if
            end do
         end do
         if (hist(rint) .gt. 0 ) then
            eterm1 = eu(rint)/hist(rint) - 1
            fterm1 = fterm(rint)/hist(rint)*conver2
            tterm1 = tterm(rint)/hist(rint)*conver2*molmas
         else
            eu(i) = 0.0d0
         end if
         bt2 = bt2 - 2*pi*ttt*(rdist)**2*eterm1*conver
         btcf = btcf + 4*pi*ttt*(rdist)**2*fterm1*conver
         btct = btct + 4*pi*ttt*(rdist)**2*tterm1*conver
         write (ixyz,14) rdist, eterm1, 
     &      fterm1,tterm1, bt2+btcf+btct, btcf,btct
 14      format(f8.3,6f11.5)
         ttt = width/(((rdist-2.8)**3+10.0)
     &          *exp(-0.15*(rdist-2.8)**2)+3.0)
         rdist = rdist + ttt 
         if (rdist .gt. 16.0 ) goto 19
      end do
 19   continue
      close (unit=ixyz)
      write(*,*) 'JOB IS FINISHED !!!'
      end

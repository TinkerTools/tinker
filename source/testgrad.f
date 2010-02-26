c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program testgrad  --  derivative test; Cartesian version  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "testgrad" computes and compares the analytical and numerical
c     gradient vectors of the potential energy function with respect
c     to Cartesian coordinates
c
c
      program testgrad
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'deriv.i'
      include 'energi.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'solute.i'
      include 'usage.i'
      integer i,j,next
      real*8 etot,f,f0,eps,eps0,old,energy
      real*8 eb0,ea0,eba0,eub0,eaa0,eopb0
      real*8 eopd0,eid0,eit0,et0,ept0,ebt0
      real*8 ett0,ev0,ec0,ecd0,ed0,em0,ep0
      real*8 er0,es0,elf0,eg0,ex0
      real*8 ndeb(3,maxatm),ndea(3,maxatm)
      real*8 ndeba(3,maxatm),ndeub(3,maxatm)
      real*8 ndeaa(3,maxatm),ndeopb(3,maxatm)
      real*8 ndeopd(3,maxatm),ndeid(3,maxatm)
      real*8 ndeit(3,maxatm),ndet(3,maxatm)
      real*8 ndept(3,maxatm),ndebt(3,maxatm)
      real*8 ndett(3,maxatm),ndev(3,maxatm)
      real*8 ndec(3,maxatm),ndecd(3,maxatm)
      real*8 nded(3,maxatm),ndem(3,maxatm)
      real*8 ndep(3,maxatm),nder(3,maxatm)
      real*8 ndes(3,maxatm),ndelf(3,maxatm)
      real*8 ndeg(3,maxatm),ndex(3,maxatm)
      real*8 detot(3,maxatm),ndetot(3,maxatm)
      real*8 denorm(maxatm),ndenorm(maxatm)
      real*8 totnorm,ntotnorm,rms,nrms
      logical exist,query
      logical doanalyt,donumer,dofull
      character*1 answer
      character*1 axis(3)
      character*120 record
      character*120 string
      data axis  / 'X','Y','Z' /
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     decide whether to do an analytical gradient calculation
c
      doanalyt = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Compute the Analytical Gradient Vector [Y] :  ',$)
         read (input,20)  record
   20    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  doanalyt = .false.
c
c     decide whether to do a numerical gradient calculation
c
      donumer = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,30)
   30    format (/,' Compute the Numerical Gradient Vector [Y] :   ',$)
         read (input,40)  record
   40    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  donumer = .false.
c
c     get the stepsize for numerical gradient calculation
c
      if (donumer) then
         eps = -1.0d0
         eps0 = 1.0d-5
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=50,end=50)  eps
            query = .false.
         end if
   50    continue
         if (query) then
            write (iout,60)  eps0
   60       format (/,' Enter a Numerical Stepsize [',d7.1,
     &                 ' Ang] :  ',$)
            read (input,70,err=50)  eps
   70       format (f20.0)
         end if
         if (eps .le. 0.0d0)  eps = eps0
      end if
c
c     decide whether to output results by gradient component
c
      dofull = .true.
      if (n .gt. 100) then
         dofull = .false.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,80)
   80       format (/,' Output Breakdown by Gradient Component',
     &                 ' [N] :  ',$)
            read (input,90)  record
   90       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  dofull = .true.
      end if
c
c     compute the analytical gradient components
c
      if (doanalyt) then
         call gradient (etot,detot)
      end if
c
c     print the total potential energy of the system
c
      if (doanalyt) then
         if (digits .ge. 8) then
            write (iout,100)  etot
  100       format (/,' Total Potential Energy :',8x,f20.8,' Kcal/mole')
         else if (digits .ge. 6) then
            write (iout,110)  etot
  110       format (/,' Total Potential Energy :',8x,f18.6,' Kcal/mole')
         else
            write (iout,120)  etot
  120       format (/,' Total Potential Energy :',8x,f16.4,' Kcal/mole')
         end if
c
c     print the energy breakdown over individual components
c
         write (iout,130)
  130    format (/,' Potential Energy Breakdown by Individual',
     &              ' Components :')
         if (digits .ge. 8) then
            write (iout,140)
  140       format (/,'  Energy',7x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
     &              /,'  Terms',8x,'EAA',13x,'EOPB',12x,'EOPD',
     &                 12x,'EID',
     &              /,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',
     &              /,15x,'ETT',13x,'EV',14x,'EC',14x,'ECD',
     &              /,15x,'ED',14x,'EM',14x,'EP',14x,'ER',
     &              /,15x,'ES',14x,'ELF',13x,'EG',14x,'EX')
            write (iout,150)  eb,ea,eba,eub,eaa,eopb,eopd,eid,
     &                        eit,et,ept,ebt,ett,ev,ec,ecd,ed,
     &                        em,ep,er,es,elf,eg,ex
  150       format (/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,
     &                 /,6x,4f16.8,/,6x,4f16.8)
         else if (digits .ge. 6) then
            write (iout,160)
  160       format (/,'  Energy',6x,'EB',12x,'EA',12x,'EBA',11x,'EUB',
     &                 11x,'EAA',
     &              /,'  Terms',7x,'EOPB',10x,'EOPD',10x,'EID',
     &                 11x,'EIT',11x,'ET',
     &              /,14x,'EPT',11x,'EBT',11x,'ETT',11x,'EV',12x,'EC',
     &              /,14x,'ECD',11x,'ED',12x,'EM',12x,'EP',12x,'ER',
     &              /,14x,'ES',12x,'ELF',11x,'EG',12x,'EX')
            write (iout,170)  eb,ea,eba,eub,eaa,eopb,eopd,eid,
     &                        eit,et,ept,ebt,ett,ev,ec,ecd,ed,
     &                        em,ep,er,es,elf,eg,ex
  170       format (/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6,
     &                 /,6x,4f14.6)
         else
            write (iout,180)
  180       format (/,'  Energy',6x,'EB',10x,'EA',10x,'EBA',9x,'EUB',
     &                 9x,'EAA',9x,'EOPB',
     &              /,'  Terms',7x,'EOPD',8x,'EID',9x,'EIT',9x,'ET',
     &                 10x,'EPT',9x,'EBT',
     &              /,14x,'ETT',9x,'EV',10x,'EC',10x,'ECD',9x,'ED',
     &                 10x,'EM',
     &              /,14x,'EP',10x,'ER',10x,'ES',10x,'ELF',9x,'EG',
     &                 10x,'EX')
            write (iout,190)  eb,ea,eba,eub,eaa,eopb,eopd,eid,
     &                        eit,et,ept,ebt,ett,ev,ec,ecd,ed,
     &                        em,ep,er,es,elf,eg,ex
  190       format (/,6x,6f12.4,/,6x,6f12.4,/,6x,6f12.4,/,6x,6f12.4)
         end if
      end if
c
c     print a header for the gradients of individual potentials
c
      if (dofull) then
         write (iout,200)
  200    format (/,' Cartesian Gradient Breakdown by Individual',
     &                 ' Components :')
         if (digits .ge. 8) then
            write (iout,210)
  210       format (/,2x,'Atom',9x,'d EB',12x,'d EA',12x,'d EBA',
     &                 11x,'d EUB',
     &              /,2x,'Axis',9x,'d EAA',11x,'d EOPB',10x,'d EOPD',
     &                 10x,'d EID',
     &              /,2x,'Type',9x,'d EIT',11x,'d ET',12x,'d EPT',
     &                 11x,'d EBT',
     &              /,15x,'d ETT',11x,'d EV',12x,'d EC',12x,'d ECD',
     &              /,15x,'d ED',12x,'d EM',12x,'d EP',12x,'d ER',
     &              /,15x,'d ES',12x,'d ELF',11x,'d EG',12x,'d EX')
         else if (digits .ge. 6) then
            write (iout,220)
  220       format (/,2x,'Atom',8x,'d EB',10x,'d EA',10x,'d EBA',
     &                 9x,'d EUB',9x,'d EAA',
     &              /,2x,'Axis',8x,'d EOPB',8x,'d EOPD',8x,'d EID',
     &                 9x,'d EIT',9x,'d ET',
     &              /,2x,'Type',8x,'d EPT',9x,'d EBT',9x,'d ETT',
     &                 9x,'d EV',10x,'d EC',
     &              /,14x,'d ECD',9x,'d ED',10x,'d EM',10x,'d EP',
     &                 10x,'d ER',
     &              /,14x,'d ES',10x,'d ELF',9x,'d EG',10x,'d EX')
         else
            write (iout,230)
  230       format (/,2x,'Atom',6x,'d EB',8x,'d EA',8x,'d EBA',
     &                 7x,'d EUB',7x,'d EAA',7x,'d EOPB',
     &              /,2x,'Axis',6x,'d EOPD',6x,'d EID',7x,'d EIT',
     &                 7x,'d ET',8x,'d EPT',7x,'d EBT',
     &              /,2x,'Type',6x,'d ETT',7x,'d EV',8x,'d EC',
     &                 8x,'d ECD',7x,'d ED',8x,'d EM',
     &              /,12x,'d EP',8x,'d ER',8x,'d ES',8x,'d ELF',
     &                 7x,'d EG',8x,'d EX')
         end if
      end if
c
c     get the Cartesian component two-sided numerical gradients
c
      do i = 1, n
         if (donumer .and. use(i)) then
            do j = 1, 3
               if (j .eq. 1) then
                  old = x(i)
                  x(i) = x(i) - 0.5d0*eps
               else if (j .eq. 2) then
                  old = y(i)
                  y(i) = y(i) - 0.5d0*eps
               else if (j .eq. 3) then
                  old = z(i)
                  z(i) = z(i) - 0.5d0*eps
               end if
               f0 = energy ()
               eb0 = eb
               ea0 = ea
               eba0 = eba
               eub0 = eub
               eaa0 = eaa
               eopb0 = eopb
               eopd0 = eopd
               eid0 = eid
               eit0 = eit
               et0 = et
               ept0 = ept
               ebt0 = ebt
               ett0 = ett
               ev0 = ev
               ec0 = ec
               ecd0 = ecd
               ed0 = ed
               em0 = em
               ep0 = ep
               er0 = er
               es0 = es
               elf0 = elf
               eg0 = eg
               ex0 = ex
               if (j .eq. 1) then
                  x(i) = x(i) + eps
               else if (j .eq. 2) then
                  y(i) = y(i) + eps
               else if (j .eq. 3) then
                  z(i) = z(i) + eps
               end if
               f = energy ()
               if (j .eq. 1) then
                  x(i) = old
               else if (j .eq. 2) then
                  y(i) = old
               else if (j .eq. 3) then
                  z(i) = old
               end if
               ndeb(j,i) = (eb - eb0) / eps
               ndea(j,i) = (ea - ea0) / eps
               ndeba(j,i) = (eba - eba0) / eps
               ndeub(j,i) = (eub - eub0) / eps
               ndeaa(j,i) = (eaa - eaa0) / eps
               ndeopb(j,i) = (eopb - eopb0) / eps
               ndeopd(j,i) = (eopd - eopd0) / eps
               ndeid(j,i) = (eid - eid0) / eps
               ndeit(j,i) = (eit - eit0) / eps
               ndet(j,i) = (et - et0) / eps
               ndept(j,i) = (ept - ept0) / eps
               ndebt(j,i) = (ebt - ebt0) / eps
               ndett(j,i) = (ett - ett0) / eps
               ndev(j,i) = (ev - ev0) / eps
               ndec(j,i) = (ec - ec0) / eps
               ndecd(j,i) = (ecd - ecd0) / eps
               nded(j,i) = (ed - ed0) / eps
               ndem(j,i) = (em - em0) / eps
               ndep(j,i) = (ep - ep0) / eps
               nder(j,i) = (er - er0) / eps
               ndes(j,i) = (es - es0) / eps
               ndelf(j,i) = (elf - elf0) / eps
               ndeg(j,i) = (eg - eg0) / eps
               ndex(j,i) = (ex - ex0) / eps
               ndetot(j,i) = (f - f0) / eps
            end do
         end if
c
c     print analytical gradients of each energy term for each atom
c
         if (dofull .and. use(i)) then
            do j = 1, 3
               if (doanalyt) then
                  if (digits .ge. 8) then
                     write (iout,240)  i,deb(j,i),dea(j,i),deba(j,i),
     &                                 deub(j,i),axis(j),deaa(j,i),
     &                                 deopb(j,i),deopd(j,i),deid(j,i),
     &                                 deit(j,i),det(j,i),dept(j,i),
     &                                 debt(j,i),dett(j,i),dev(j,i),
     &                                 dec(j,i),decd(j,i),ded(j,i),
     &                                 dem(j,i),dep(j,i),der(j,i),
     &                                 des(j,i),delf(j,i),deg(j,i),
     &                                 dex(j,i)
  240                format (/,i6,4f16.8,/,5x,a1,4f16.8,
     &                          /,' Anlyt',4f16.8,/,6x,4f16.8,
     &                          /,6x,4f16.8,/,6x,4f16.8)
                  else if (digits .ge. 6) then
                     write (iout,250)  i,deb(j,i),dea(j,i),deba(j,i),
     &                                 deub(j,i),deaa(j,i),axis(j),
     &                                 deopb(j,i),deopd(j,i),deid(j,i),
     &                                 deit(j,i),det(j,i),dept(j,i),
     &                                 debt(j,i),dett(j,i),dev(j,i),
     &                                 dec(j,i),decd(j,i),ded(j,i),
     &                                 dem(j,i),dep(j,i),der(j,i),
     &                                 des(j,i),delf(j,i),deg(j,i),
     &                                 dex(j,i)
  250                format (/,i6,5f14.6,/,5x,a1,5f14.6,/,' Anlyt',
     &                          5f14.6,/,6x,5f14.6,/,6x,4f14.6)
                  else
                     write (iout,260)  i,deb(j,i),dea(j,i),deba(j,i),
     &                                 deub(j,i),deaa(j,i),deopb(j,i),
     &                                 axis(j),deopd(j,i),deid(j,i),
     &                                 deit(j,i),det(j,i),dept(j,i),
     &                                 debt(j,i),dett(j,i),dev(j,i),
     &                                 dec(j,i),decd(j,i),ded(j,i),
     &                                 dem(j,i),dep(j,i),der(j,i),
     &                                 des(j,i),delf(j,i),deg(j,i),
     &                                 dex(j,i)
  260                format (/,i6,6f12.4,/,5x,a1,6f12.4,/,' Anlyt',
     &                          6f12.4,/,6x,6f12.4)
                  end if
               end if
c
c     print numerical gradients of each energy term for each atom
c
               if (donumer) then
                  if (digits .ge. 8) then
                     write (iout,270)  i,ndeb(j,i),ndea(j,i),ndeba(j,i),
     &                                 ndeub(j,i),axis(j),ndeaa(j,i),
     &                                 ndeopb(j,i),ndeopd(j,i),
     &                                 ndeid(j,i),ndeit(j,i),ndet(j,i),
     &                                 ndept(j,i),ndebt(j,i),ndett(j,i),
     &                                 ndev(j,i),ndec(j,i),ndecd(j,i),
     &                                 nded(j,i),ndem(j,i),ndep(j,i),
     &                                 nder(j,i),ndes(j,i),ndelf(j,i),
     &                                 ndeg(j,i),ndex(j,i)
  270                format (/,i6,4f16.8,/,5x,a1,4f16.8,
     &                          /,' Numer',4f16.8,/,6x,4f16.8,
     &                          /,6x,4f16.8,/,6x,4f16.8)
                  else if (digits .ge. 6) then
                     write (iout,280)  i,ndeb(j,i),ndea(j,i),ndeba(j,i),
     &                                 ndeub(j,i),ndeaa(j,i),axis(j),
     &                                 ndeopb(j,i),ndeopd(j,i),
     &                                 ndeid(j,i),ndeit(j,i),ndet(j,i),
     &                                 ndept(j,i),ndebt(j,i),ndett(j,i),
     &                                 ndev(j,i),ndec(j,i),ndecd(j,i),
     &                                 nded(j,i),ndem(j,i),ndep(j,i),
     &                                 nder(j,i),ndes(j,i),ndelf(j,i),
     &                                 ndeg(j,i),ndex(j,i)
  280                format (/,i6,5f14.6,/,5x,a1,5f14.6,/,' Numer',
     &                          5f14.6,/,6x,5f14.6,/,6x,4f14.6)
                  else
                     write (iout,290)  i,ndeb(j,i),ndea(j,i),
     &                                 ndeba(j,i),ndeub(j,i),ndeaa(j,i),
     &                                 ndeopb(j,i),axis(j),ndeopd(j,i),
     &                                 ndeid(j,i),ndeit(j,i),ndet(j,i),
     &                                 ndept(j,i),ndebt(j,i),ndett(j,i),
     &                                 ndev(j,i),ndec(j,i),ndecd(j,i),
     &                                 nded(j,i),ndem(j,i),ndep(j,i),
     &                                 nder(j,i),ndes(j,i),ndelf(j,i),
     &                                 ndeg(j,i),ndex(j,i)
  290                format (/,i6,6f12.4,/,5x,a1,6f12.4,/,' Numer',
     &                          6f12.4,/,6x,6f12.4)
                  end if
               end if
            end do
         end if
      end do
c
c     print the total gradient components for each atom
c
      write (iout,300)
  300 format (/,' Cartesian Gradient Breakdown over Individual Atoms :')
      if (digits .ge. 8) then
         write (iout,310)
  310    format (/,2x,'Type',4x,'Atom',10x,'dE/dX',11x,'dE/dY',
     &              11x,'dE/dZ',11x,'Norm',/)
      else if (digits .ge. 6) then
         write (iout,320)
  320    format (/,2x,'Type',6x,'Atom',11x,'dE/dX',9x,'dE/dY',
     &              9x,'dE/dZ',11x,'Norm',/)
      else
         write (iout,330)
  330    format (/,2x,'Type',6x,'Atom',14x,'dE/dX',7x,'dE/dY',
     &              7x,'dE/dZ',10x,'Norm',/)
      end if
      totnorm = 0.0d0
      ntotnorm = 0.0d0
      do i = 1, n
         if (doanalyt .and. use(i)) then
            denorm(i) = detot(1,i)**2 + detot(2,i)**2 +
     &                       detot(3,i)**2
            totnorm = totnorm + denorm(i)
            denorm(i) = sqrt(denorm(i))
            if (digits .ge. 8) then
               write (iout,340)  i,(detot(j,i),j=1,3),denorm(i)
  340          format (' Anlyt',i8,1x,3f16.8,f16.8)
            else if (digits .ge. 6) then
               write (iout,350)  i,(detot(j,i),j=1,3),denorm(i)
  350          format (' Anlyt',2x,i8,3x,3f14.6,2x,f14.6)
            else
               write (iout,360)  i,(detot(j,i),j=1,3),denorm(i)
  360          format (' Anlyt',2x,i8,7x,3f12.4,2x,f12.4)
            end if
         end if
         if (donumer .and. use(i)) then
            ndenorm(i) = ndetot(1,i)**2 + ndetot(2,i)**2 +
     &                        ndetot(3,i)**2
            ntotnorm = ntotnorm + ndenorm(i)
            ndenorm(i) = sqrt(ndenorm(i))
            if (digits .ge. 8) then
               write (iout,370)  i,(ndetot(j,i),j=1,3),ndenorm(i)
  370          format (' Numer',i8,1x,3f16.8,f16.8)
            else if (digits .ge. 6) then
               write (iout,380)  i,(ndetot(j,i),j=1,3),ndenorm(i)
  380          format (' Numer',2x,i8,3x,3f14.6,2x,f14.6)
            else
               write (iout,390)  i,(ndetot(j,i),j=1,3),ndenorm(i)
  390          format (' Numer',2x,i8,7x,3f12.4,2x,f12.4)
            end if
         end if
      end do
c
c     print the total norm for the analytical gradient
c
      write (iout,400)
  400 format (/,' Total Gradient Norm and RMS Gradient per Atom :')
      if (doanalyt) then
         totnorm = sqrt(totnorm)
         if (digits .ge. 8) then
            write (iout,410)  totnorm
  410       format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f20.8)
         else if (digits .ge. 6) then
            write (iout,420)  totnorm
  420       format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f18.6)
         else
            write (iout,430)  totnorm
  430       format (/,' Anlyt',6x,'Total Gradient Norm Value',6x,f16.4)
         end if
      end if
c
c     print the total norm for the numerical gradient
c
      if (donumer) then
         ntotnorm = sqrt(ntotnorm)
         if (digits .ge. 8) then
            write (iout,440)  ntotnorm
  440       format (' Numer',6x,'Total Gradient Norm Value',6x,f20.8)
         else if (digits .ge. 6) then
            write (iout,450)  ntotnorm
  450       format (' Numer',6x,'Total Gradient Norm Value',6x,f18.6)
         else
            write (iout,460)  ntotnorm
  460       format (' Numer',6x,'Total Gradient Norm Value',6x,f16.4)
         end if
      end if
c
c     print the rms per atom norm for the analytical gradient
c
      if (doanalyt) then
         rms = totnorm / sqrt(dble(nuse))
         if (digits .ge. 8) then
            write (iout,470)  rms
  470       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
     &                 4x,f20.8)
         else if (digits .ge. 6) then
            write (iout,480)  rms
  480       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
     &                 4x,f18.6)
         else
            write (iout,490)  rms
  490       format (/,' Anlyt',6x,'RMS Gradient over All Atoms',
     &                 4x,f16.4)
         end if
      end if
c
c     print the rms per atom norm for the numerical gradient
c
      if (donumer) then
         nrms = ntotnorm / sqrt(dble(nuse))
         if (digits .ge. 8) then
            write (iout,500)  nrms
  500       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f20.8)
         else if (digits .ge. 6) then
            write (iout,510)  nrms
  510       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f18.6)
         else
            write (iout,520)  nrms
  520       format (' Numer',6x,'RMS Gradient over All Atoms',4x,f16.4)
         end if
      end if
c
c     perform any final tasks before program exit
c
      call final
      end

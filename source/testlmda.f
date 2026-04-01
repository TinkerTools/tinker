c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung & Jay W. Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ####################################################
c     ##                                                ##
c     ##  program testlmda  --  lambda derivative test  ##
c     ##                                                ##
c     ####################################################
c
c
c     "testlmda" computes and compares the analytical and numerical
c     derivatives of the potential energy function with respect
c     to the lambda parameter
c
c
      program testlmda
      use atoms
      use deriv
      use dlmda
      use energi
      use files
      use inform
      use iounit
      use mutant
      use ost
      use usage
      use virial
      implicit none
      integer i,j,ixyz
      integer next,frame
      integer freeunit
      integer nask
      real*8 eval,energy
      real*8 eps,eps0
      real*8 lmda,lmda0
      real*8 adelmda,adevlmda
      real*8 ademlmda,adeplmda
      real*8 adelmda2,adevlmda2
      real*8 ademlmda2,adeplmda2
      real*8 ndelmda,ndevlmda
      real*8 ndemlmda,ndeplmda
      real*8 ndelmda2,ndevlmda2
      real*8 ndemlmda2,ndeplmda2
      real*8 esum2,esum1,esum0
      real*8 em2,em1,em0
      real*8 ep2,ep1,ep0
      real*8 ev2,ev1,ev0
      real*8 oldvdl,oldeml
      real*8 oldepl
      real*8 denorm,ndenorm
      real*8 totnorm,ntotnorm,rms,nrms
      real*8 adldvir(3,3)
      real*8 ndldvir(3,3)
      real*8 ndldepvir(3,3)
      real*8 ndldemvir(3,3)
      real*8 ndldevvir(3,3)
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: adldesum(:,:)
      real*8, allocatable :: adldev(:,:)
      real*8, allocatable :: adldem(:,:)
      real*8, allocatable :: adldep(:,:)
      real*8, allocatable :: ndldesum(:,:)
      real*8, allocatable :: ndldev(:,:)
      real*8, allocatable :: ndldem(:,:)
      real*8, allocatable :: ndldep(:,:)
      logical exist,query
      logical doanalyt,donumer
      character*1 answer
      character*240 xyzfile
      character*240 record
      character*240 string
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     set flag to use lambda derivative
c
      use_dlmda = .true.
      use_pol4i = .true.
      use_pol4f = .true.
c
c     decide whether to do an analytical derivative calculation
c
      doanalyt = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Compute the Analytical Derivative [Y] :  ',$)
         read (input,20)  record
   20    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  doanalyt = .false.
c
c     decide whether to do a numerical derivative calculation
c
      donumer = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,30)
   30    format (/,' Compute the Numerical Derivative [Y] :   ',$)
         read (input,40)  record
   40    format (a240)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  donumer = .false.
c
c     get the stepsize for numerical derivative calculation
c
      if (donumer) then
         eps = -1.0d0
         eps0 = 0.00001d0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=50,end=50)  eps
            query = .false.
         end if
   50    continue
         if (query) then
            write (iout,60)  eps0
   60       format (/,' Enter Finite Difference Stepsize [',d8.1,
     &                 ' Lambda] :  ',$)
            read (input,70,err=50)  eps
   70       format (f20.0)
         end if
         if (eps .le. 0.0d0)  eps = eps0
      end if
c
c     reopen the coordinates file and read the first structure
c
      frame = 0
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
      if (doanalyt) then
         allocate (adldesum(3,n))
         allocate (adldev(3,n))
         allocate (adldem(3,n))
         allocate (adldep(3,n))
      end if
      if (donumer) then
         allocate (ndldesum(3,n))
         allocate (ndldev(3,n))
         allocate (ndldem(3,n))
         allocate (ndldep(3,n))
      end if
c
c     perform analysis for each successive coordinate structure
c
      do while (.not. abort)
         frame = frame + 1
         if (frame .gt. 1) then
            write (iout,80)  frame
   80       format (/,' Analysis for Archive Structure :',8x,i8)
         end if
c
c     compute the analytical lambda derivatives
c
         if (doanalyt) then
            call altelec
            call gradient (eval,derivs)
            adelmda = delmda
            adevlmda = devlmda
            ademlmda = demlmda
            adeplmda = deplmda
            adelmda2 = delmda2
            adevlmda2 = devlmda2
            ademlmda2 = demlmda2
            adeplmda2 = deplmda2
            do i = 1, n
               do j = 1, 3
                  adldesum(j,i) = dldesum(j,i)
                  adldev(j,i) = dldev(j,i)
                  adldem(j,i) = dldem(j,i)
                  adldep(j,i) = dldep(j,i)
               end do
            end do
            do i = 1, 3
               do j = 1, 3
                  adldvir(j,i) = dldvir(j,i)
               end do
            end do
         end if
c
c     compute the numerical lambda derivatives
c
         if (donumer) then
            oldvdl = vlambda
            oldeml = elambda
            oldepl = plambda
            vlambda = oldvdl + eps
            elambda = oldeml + eps
            plambda = oldepl + eps
            call altelec
            esum2 = energy ()
            ev2 = ev
            em2 = em
            ep2 = ep
            vlambda = oldvdl - eps
            elambda = oldeml - eps
            plambda = oldepl - eps
            call altelec
            esum0 = energy ()
            ev0 = ev
            em0 = em
            ep0 = ep
            vlambda = oldvdl
            elambda = oldeml
            plambda = oldepl
            call altelec
            esum1 = energy ()
            ev1 = ev
            em1 = em
            ep1 = ep
            ndelmda = (esum2 - esum0) / (2.0d0 * eps)
            ndevlmda = (ev2 - ev0) / (2.0d0 * eps)
            ndemlmda = (em2 - em0) / (2.0d0 * eps)
            ndeplmda = (ep2 - ep0) / (2.0d0 * eps)
            ndelmda2 = (esum2 - 2.0d0 * esum1 + esum0) / (eps*eps)
            ndevlmda2 = (ev2 - 2.0d0 * ev1 + ev0) / (eps*eps)
            ndemlmda2 = (em2 - 2.0d0 * em1 + em0) / (eps*eps)
            ndeplmda2 = (ep2 - 2.0d0 * ep1 + ep0) / (eps*eps)
            vlambda = oldvdl + eps
            elambda = oldeml + eps
            plambda = oldepl + eps
            call altelec
            call gradient (eval,derivs)
            do i = 1, n
               do j = 1, 3
                  ndldesum(j,i) = derivs(j,i)
                  ndldev(j,i) = dev(j,i)
                  ndldem(j,i) = dem(j,i)
                  ndldep(j,i) = dep(j,i)
               end do
            end do
            do i = 1, 3
               do j = 1, 3
                  ndldvir(j,i) = vir(j,i)
                  ndldepvir(j,i) = epvir(j,i)
                  ndldemvir(j,i) = emvir(j,i)
                  ndldevvir(j,i) = evvir(j,i)
               end do
            end do
            vlambda = oldvdl - eps
            elambda = oldeml - eps
            plambda = oldepl - eps
            call altelec
            call gradient (eval,derivs)
            do i = 1, n
               do j = 1, 3
                  ndldesum(j,i) = (ndldesum(j,i)-derivs(j,i))
     &                                                 / (2.0d0 * eps)
                  ndldev(j,i) = (ndldev(j,i)-dev(j,i)) / (2.0d0 * eps)
                  ndldem(j,i) = (ndldem(j,i)-dem(j,i)) / (2.0d0 * eps)
                  ndldep(j,i) = (ndldep(j,i)-dep(j,i)) / (2.0d0 * eps)
               end do
            end do
            do i = 1, 3
               do j = 1, 3
                  ndldvir(j,i) = (ndldvir(j,i) - vir(j,i))
     &                         / (2.0d0 * eps)
                  ndldepvir(j,i) = (ndldepvir(j,i)-epvir(j,i))
     &                           / (2.0d0 * eps)
                  ndldemvir(j,i) = (ndldemvir(j,i)-emvir(j,i))
     &                           / (2.0d0 * eps)
                  ndldevvir(j,i) = (ndldevvir(j,i)-evvir(j,i))
     &                           / (2.0d0 * eps)
               end do
            end do
c
c     apply chain rule if using global lambda in ost
c
            if (use_ost) then
               ndeplmda2 = ndeplmda2 * dplambda*dplambda
     &                           + ndeplmda * d2plambda
               ndeplmda = ndeplmda * dplambda
               ndevlmda2 = ndevlmda2 * dvlambda*dvlambda
     &                           + ndevlmda * d2vlambda
               ndevlmda = ndevlmda * dvlambda
               ndemlmda2 = ndemlmda2 * delambda*delambda
     &                           + ndemlmda * d2elambda
               ndemlmda = ndemlmda * delambda
               ndelmda = ndeplmda + ndevlmda + ndemlmda
               ndelmda2 = ndeplmda2 + ndevlmda2 + ndemlmda2
               do i = 1, n
                  do j = 1, 3
                     ndldep(j,i) = ndldep(j,i) * dplambda
                     ndldem(j,i) = ndldem(j,i) * delambda
                     ndldev(j,i) = ndldev(j,i) * dvlambda
                     ndldesum(j,i) = ndldep(j,i)+ndldem(j,i)+ndldev(j,i)
                  end do
               end do
               do i = 1, 3
                  do j = 1, 3
                     ndldepvir(j,i) = ndldepvir(j,i) * dplambda
                     ndldemvir(j,i) = ndldemvir(j,i) * delambda
                     ndldevvir(j,i) = ndldevvir(j,i) * dvlambda
                     ndldvir(j,i) = ndldepvir(j,i) + ndldemvir(j,i)
     &                            + ndldevvir(j,i)
                  end do
               end do
            end if
            vlambda = oldvdl
            elambda = oldeml
            plambda = oldepl
            call altelec
         end if
c
c     print the analytical lambda derivatives
c
         if (doanalyt) then
            write (iout,90)  'dE/dL', 'dEV/dL', 'dEM/dL', 'dEP/dL',
     &                         adelmda, adevlmda, ademlmda, adeplmda
  90       format (/,' Analytical Lambda Derivatives :', 4a14, /, 32x,
     &                 4f14.6)
         end if
c
c     print the numerical lambda derivatives
c
         if (donumer) then
            write (iout,100)  'dE/dL', 'dEV/dL', 'dEM/dL', 'dEP/dL',
     &                         ndelmda, ndevlmda, ndemlmda, ndeplmda
  100       format (/,' Numerical Lambda Derivatives : ', 4a14, /, 32x,
     &                 4f14.6)
         end if
c
c     print the analytical second order lambda derivatives
c
         if (doanalyt) then
            write (iout,110)  'd2E/dL2', 'd2EV/dL2', 'd2EM/dL2', 
     &                         'd2EP/dL2',
     &                         adelmda2, adevlmda2, ademlmda2, adeplmda2
  110       format (/,' Analytical 2nd Lambda Derivatives :',
     &                 4a14, /, 36x, 4f14.6)
         end if
c
c     print the numerical second order lambda derivatives
c
         if (donumer) then
            write (iout,120)  'd2E/dL2', 'd2EV/dL2', 'd2EM/dL2', 
     &                         'd2EP/dL2',
     &                         ndelmda2, ndevlmda2, ndemlmda2, ndeplmda2
  120       format (/,' Numerical 2nd Lambda Derivatives : ',
     &                 4a14, /, 36x, 4f14.6)
         end if
c
c     print the total force lambda derivatives for each atom
c
         if (doanalyt .or. donumer) then
            write (iout,130)
  130       format (/,' Lambda Gradient Breakdown over Individual',
     &                 ' Atoms :')
            write (iout,140)
  140       format (/,2x,'Type',6x,'Atom',11x,'dFx/dL',8x,'dFy/dL',
     &                 8x,'dFz/dL',10x,'Norm',/)
         end if
         totnorm = 0.0d0
         ntotnorm = 0.0d0
         do i = 1, n
            if (doanalyt .and. use(i)) then
               denorm = adldesum(1,i)**2 + adldesum(2,i)**2
     &                        + adldesum(3,i)**2
               totnorm = totnorm + denorm
               denorm = sqrt(denorm)
               write (iout,150)  i,(adldesum(j,i),j=1,3),denorm
  150          format (' Anlyt',2x,i8,3x,3f14.6,2x,f14.6)
            end if
            if (donumer .and. use(i)) then
               ndenorm = ndldesum(1,i)**2 + ndldesum(2,i)**2
     &                         + ndldesum(3,i)**2
               ntotnorm = ntotnorm + ndenorm
               ndenorm = sqrt(ndenorm)
               write (iout,160)  i,(ndldesum(j,i),j=1,3),ndenorm
  160          format (' Numer',2x,i8,3x,3f14.6,2x,f14.6)
            end if
         end do
c
c     print the total norm for the analytical lambda derivative
c
         if (doanalyt .or. donumer) then
            write (iout,170)
  170       format (/,' Total Gradient Norm and RMS Gradient',
     &                 ' per Atom :',/)
         end if
         if (doanalyt) then
            totnorm = sqrt(totnorm)
            write (iout,180)  totnorm
  180       format (' Anlyt',6x,'Total Gradient Norm Value',
     &                 6x,f18.6)
         end if
c
c     print the total norm for the numerical lambda derivative
c
         if (donumer) then
            ntotnorm = sqrt(ntotnorm)
            write (iout,190)  ntotnorm
  190       format (' Numer',6x,'Total Gradient Norm Value',
     &                 6x,f18.6)
         end if
c
c     print the rms per atom norm for the analytical lambda derivative
c
         if (doanalyt .or. donumer) then
            write (iout,200)
  200       format ()
         end if
         if (doanalyt) then
            rms = totnorm / sqrt(dble(nuse))
            write (iout,210)  rms
  210       format (' Anlyt',6x,'RMS Gradient over All Atoms',
     &                 4x,f18.6)
         end if
c
c     print the rms per atom norm for the numerical lambda derivative
c
         if (donumer) then
            nrms = ntotnorm / sqrt(dble(nuse))
            write (iout,220)  nrms
  220       format (' Numer',6x,'RMS Gradient over All Atoms',
     &                 4x,f18.6)
         end if
c
c     print the components of the analytical internal virial
c
      write (iout,230)  (adldvir(1,i),adldvir(2,i),adldvir(3,i),i=1,3)
  230 format (/,' Analytical dV/dL :',8x,3f13.3,
     &           /,27x,3f13.3,/,27x,3f13.3)
c
c     print the components of the numerical internal virial
c
      write (iout,240)  (ndldvir(1,i),ndldvir(2,i),ndldvir(3,i),i=1,3)
  240 format (/,' Numerical dV/dL :',9x,3f13.3,
     &           /,27x,3f13.3,/,27x,3f13.3)
c
c     attempt to read next structure from the coordinate file
c
         call readxyz (ixyz)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      if (doanalyt) then
         deallocate (adldesum)
         deallocate (adldev)
         deallocate (adldem)
         deallocate (adldep)
      end if
      if (donumer) then
         deallocate (ndldesum)
         deallocate (ndldev)
         deallocate (ndldem)
         deallocate (ndldep)
      end if
c
c     perform any final tasks before program exit
c
      call final
      end

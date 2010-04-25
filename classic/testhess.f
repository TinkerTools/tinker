c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program testhess  --  Hessian matrix test; cart. version  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "testhess" computes and compares the analytical and numerical
c     Hessian matrices of the potential energy function with respect
c     to Cartesian coordinates
c
c
      program testhess
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'hescut.i'
      include 'inform.i'
      include 'iounit.i'
      include 'usage.i'
      integer maxnum
      parameter (maxnum=300)
      integer i,j,k,m,ii,jj,ihes
      integer index,next,freeunit
      integer hinit(3,maxatm)
      integer hstop(3,maxatm)
      integer hindex(maxhess)
      real*8 energy,e,old,eps,eps0
      real*8 diff,delta,sum
      real*8 g(3,maxatm),g0(3,maxatm)
      real*8 h(maxhess),hdiag(3,maxatm)
      real*8 nhess(3,maxnum,3,maxnum)
      logical doanalyt,donumer
      logical dograd,dofull
      logical exist,query
      logical identical
      character*1 answer
      character*1 axis(3)
      character*120 hessfile
      character*120 record
      character*120 string
      external energy
      data axis  / 'X','Y','Z' /
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     set difference threshhold via the energy precision
c
      delta = 1.0d-4
      if (digits .ge. 6)  delta = 1.0d-6
      if (digits .ge. 8)  delta = 1.0d-8
c
c     decide whether to do an analytical Hessian calculation
c
      doanalyt = .true.
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Compute Analytical Hessian Matrix [Y] :  ',$)
         read (input,20)  record
   20    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  doanalyt = .false.
c
c     decide whether to do a numerical Hessian calculation
c
      donumer = .false.
      if (n .le. maxnum) then
         donumer = .true.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,30)
   30       format (/,' Compute Numerical Hessian Matrix [Y] :   ',$)
            read (input,40)  record
   40       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'N')  donumer = .false.
      end if
c
c     get numerical Hessian from either gradient or energy
c
      if (donumer) then
         dograd = .true.
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,50)
   50       format (/,' Numerical Hessian from Gradient',
     &                 ' or Function [G] :  ',$)
            read (input,60)  record
   60       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'F')  dograd = .false.
c
c     get the stepsize for numerical Hessian calculation
c
         eps = -1.0d0
         eps0 = 1.0d-3
         if (dograd)  eps0 = 1.0d-5
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=70,end=70)  eps
            query = .false.
         end if
   70    continue
         if (query) then
            write (iout,80)  eps0
   80       format (/,' Enter a Numerical Stepsize [',d7.1,
     &                 ' Ang] :  ',$)
            read (input,90,err=70)  eps
   90       format (f20.0)
         end if
         if (eps .le. 0.0d0)  eps = eps0
      end if
c
c     decide whether to output results by Hessian component
c
      dofull = .false.
      if (n.le.20 .and. donumer) then
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,100)
  100       format (/,' List Individual Hessian Components [N] :   ',$)
            read (input,110)  record
  110       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  dofull = .true.
      end if
c
c     get the analytical Hessian matrix elements
c
      identical = .true.
      if (doanalyt) then
         hesscut = 0.0d0
         call hessian (h,hinit,hstop,hindex,hdiag)
         call hessian (h,hinit,hstop,hindex,hdiag)
         call hessian (h,hinit,hstop,hindex,hdiag)
      end if
c
c     get the two-sided numerical Hessian matrix elements
c
      do i = 1, n
         if (donumer .and. use(i)) then
            old = x(i)
            x(i) = x(i) - 0.5d0*eps
            if (dograd) then
               call gradient (e,g)
            else
               call numgrad (energy,g,eps)
            end if
            do k = 1, n
               do j = 1, 3
                  g0(j,k) = g(j,k)
               end do
            end do
            x(i) = x(i) + eps
            if (dograd) then
               call gradient (e,g)
            else
               call numgrad (energy,g,eps)
            end if
            x(i) = old
            do k = 1, n
               do j = 1, 3
                  nhess(j,k,1,i) = (g(j,k) - g0(j,k)) / eps
               end do
            end do
            old = y(i)
            y(i) = y(i) - 0.5d0*eps
            if (dograd) then
               call gradient (e,g)
            else
               call numgrad (energy,g,eps)
            end if
            do k = 1, n
               do j = 1, 3
                  g0(j,k) = g(j,k)
               end do
            end do
            y(i) = y(i) + eps
            if (dograd) then
               call gradient (e,g)
            else
               call numgrad (energy,g,eps)
            end if
            y(i) = old
            do k = 1, n
               do j = 1, 3
                  nhess(j,k,2,i) = (g(j,k) - g0(j,k)) / eps
               end do
            end do
            old = z(i)
            z(i) = z(i) - 0.5d0*eps
            if (dograd) then
               call gradient (e,g)
            else
               call numgrad (energy,g,eps)
            end if
            do k = 1, n
               do j = 1, 3
                  g0(j,k) = g(j,k)
               end do
            end do
            z(i) = z(i) + eps
            if (dograd) then
               call gradient (e,g)
            else
               call numgrad (energy,g,eps)
            end if
            z(i) = old
            do k = 1, n
               do j = 1, 3
                  nhess(j,k,3,i) = (g(j,k) - g0(j,k)) / eps
               end do
            end do
         end if
c
c     compare the analytical and numerical diagonal elements
c
         if (doanalyt .and. donumer) then
            do j = 1, 3
               diff = abs(hdiag(j,i)-nhess(j,i,j,i))
               if (diff .gt. delta) then
                  if (identical) then
                     identical = .false.
                     write (iout,120)
  120                format (/,' Comparison of Analytical and',
     &                         ' Numerical Hessian Elements :',
     &                       //,3x,'1st Atom',4x,'2nd Atom',
     &                         10x,'Analytical',8x,'Numerical',
     &                         6x,'Difference',/)
                  end if
                  if (digits .ge. 8) then
                     write (iout,130)  i,axis(j),i,axis(j),
     &                                 hdiag(j,i),nhess(j,i,j,i),
     &                                 hdiag(j,i)-nhess(j,i,j,i)
  130                format (1x,i6,' (',a1,') ',1x,i6,' (', a1,') ',
     &                          2x,2f17.8,f16.8)
                  else if (digits .ge. 6) then
                     write (iout,140)  i,axis(j),i,axis(j),
     &                                 hdiag(j,i),nhess(j,i,j,i),
     &                                 hdiag(j,i)-nhess(j,i,j,i)
  140                format (1x,i6,' (',a1,') ',1x,i6,' (', a1,') ',
     &                          2x,2f17.6,f16.6)
                  else
                     write (iout,150)  i,axis(j),i,axis(j),
     &                                 hdiag(j,i),nhess(j,i,j,i),
     &                                 hdiag(j,i)-nhess(j,i,j,i)
  150                format (1x,i6,' (',a1,') ',1x,i6,' (', a1,') ',
     &                          2x,2f17.4,f16.4)
                  end if
               end if
c
c     compare the analytical and numerical off-diagonal elements
c
               do k = hinit(j,i), hstop(j,i)
                  index = hindex(k)
                  jj = mod(index,3)
                  if (jj .eq. 0)  jj = 3
                  ii = (index+2) / 3
                  diff = abs(h(k)-nhess(jj,ii,j,i))
                  if (diff .gt. delta) then
                     if (identical) then
                        identical = .false.
                        write (iout,160)
  160                   format (/,' Comparison of Analytical and',
     &                            ' Numerical Hessian Elements :',
     &                          //,3x,'1st Atom',4x,'2nd Atom',
     &                            10x,'Analytical',8x,'Numerical',
     &                            6x,'Difference',/)
                     end if
                     if (digits .ge. 8) then
                        write (iout,170)  i,axis(j),ii,axis(jj),
     &                                    h(k),nhess(jj,ii,j,i),
     &                                    h(k)-nhess(jj,ii,j,i)
  170                   format (1x,i6,' (',a1,') ',1x,i6,' (', a1,') ',
     &                             2x,2f17.8,f16.8)
                     else if (digits .ge. 6) then
                        write (iout,180)  i,axis(j),ii,axis(jj),
     &                                    h(k),nhess(jj,ii,j,i),
     &                                    h(k)-nhess(jj,ii,j,i)
  180                   format (1x,i6,' (',a1,') ',1x,i6,' (', a1,') ',
     &                             2x,2f17.6,f16.6)
                     else
                        write (iout,190)  i,axis(j),ii,axis(jj),
     &                                    h(k),nhess(jj,ii,j,i),
     &                                    h(k)-nhess(jj,ii,j,i)
  190                   format (1x,i6,' (',a1,') ',1x,i6,' (', a1,') ',
     &                             2x,2f17.4,f16.4)
                     end if
                  end if
               end do
            end do
         end if
      end do
c
c     success if the analytical and numerical elements are the same
c
      if (doanalyt .and. donumer) then
         if (identical) then
            write (iout,200)
  200       format (/,' Analytical and Numerical Hessian Elements',
     &                 ' are Identical')
         end if
      end if
c
c     write out the diagonal Hessian elements for each atom
c
      if (doanalyt) then
         if (digits .ge. 8) then
            write (iout,210)
  210       format (/,' Diagonal Hessian Elements for Each Atom :',
     &              //,6x,'Atom',21x,'X',19x,'Y',19x,'Z',/)
         else if (digits .ge. 6) then
            write (iout,220)
  220       format (/,' Diagonal Hessian Elements for Each Atom :',
     &              //,6x,'Atom',19x,'X',17x,'Y',17x,'Z',/)
         else
            write (iout,230)
  230       format (/,' Diagonal Hessian Elements for Each Atom :',
     &              //,6x,'Atom',17x,'X',15x,'Y',15x,'Z',/)
         end if
         do i = 1, n
            if (digits .ge. 8) then
               write (iout,240)  i,(hdiag(j,i),j=1,3)
  240          format (i10,5x,3f20.8)
            else if (digits .ge. 6) then
               write (iout,250)  i,(hdiag(j,i),j=1,3)
  250          format (i10,5x,3f18.6)
            else
               write (iout,260)  i,(hdiag(j,i),j=1,3)
  260          format (i10,5x,3f16.4)
            end if
         end do
      end if
c
c     write out the Hessian trace as sum of diagonal elements
c
      if (doanalyt) then
         sum = 0.0d0
         do i = 1, n
            do j = 1, 3
               sum = sum + hdiag(j,i)
            end do
         end do
         if (digits .ge. 8) then
            write (iout,270)  sum
  270       format (/,' Sum of Diagonal Hessian Elements :',6x,f20.8)
         else if (digits .ge. 6) then
            write (iout,280)  sum
  280       format (/,' Sum of Diagonal Hessian Elements :',6x,f18.6)
         else
            write (iout,290)  sum
  290       format (/,' Sum of Diagonal Hessian Elements :',6x,f16.4)
         end if
      end if
c
c     write out the full matrix of numerical Hessian elements
c
      if (dofull .and. donumer) then
         do i = 1, n
            do k = 1, n
               write (iout,300)  i,k
  300          format (/,' 3x3 Hessian Block for Atoms :',3x,2i8,/)
               do j = 1, 3
                  if (digits .ge. 8) then
                     write (iout,310)  (nhess(m,i,j,k),m=1,3)
  310                format (' Numer',5x,3f20.8)
                  else if (digits .ge. 6) then
                     write (iout,320)  (nhess(m,i,j,k),m=1,3)
  320                format (' Numer',5x,3f18.6)
                  else
                     write (iout,330)  (nhess(m,i,j,k),m=1,3)
  330                format (' Numer',5x,3f16.4)
                  end if
               end do
            end do
         end do
      end if
c
c     write out the full matrix of analytical Hessian elements
c
      if (doanalyt .and. .not.donumer) then
         ihes = freeunit ()
         hessfile = filename(1:leng)//'.hes'
         call version (hessfile,'new')
         open (unit=ihes,file=hessfile,status='new')
         write (iout,340)  hessfile
  340    format (/,' Hessian Matrix written to File :  ',a40)
         write (ihes,350)
  350    format (/,' Diagonal Hessian Elements  (3 per Atom)',/)
         if (digits .ge. 8) then
            write (ihes,360)  ((hdiag(j,i),j=1,3),i=1,n)
  360       format (4f16.8)
         else if (digits .ge. 6) then
            write (ihes,370)  ((hdiag(j,i),j=1,3),i=1,n)
  370       format (5f14.6)
         else
            write (ihes,380)  ((hdiag(j,i),j=1,3),i=1,n)
  380       format (6f12.4)
         end if
         do i = 1, n
            do j = 1, 3
               if (hinit(j,i) .le. hstop(j,i)) then
                  write (ihes,390)  i,axis(j)
  390             format (/,' Off-diagonal Hessian Elements for Atom',
     &                       i6,1x,a1,/)
                  if (digits .ge. 8) then
                     write (ihes,400)  (h(k),k=hinit(j,i),hstop(j,i))
  400                format (4f16.8)
                  else if (digits .ge. 6) then
                     write (ihes,410)  (h(k),k=hinit(j,i),hstop(j,i))
  410                format (5f14.6)
                  else
                     write (ihes,420)  (h(k),k=hinit(j,i),hstop(j,i))
  420                format (6f12.4)
                  end if
               end if
            end do
         end do
         close (unit=ihes)
      end if
c
c     perform any final tasks before program exit
c
      call final
      end

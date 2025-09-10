c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mdrest  --  stop system translation & rotation  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mdrest" finds and removes any translational or rotational
c     kinetic energy of the center of mass of the overall system,
c     of rigid bodies or of user-defined atom groups
c
c
      subroutine mdrest (istep)
      use mdstuf
      implicit none
      integer istep
c
c
c     check steps between center of mass motion removal
c
      if (.not. dorest)  return
      if (mod(istep,irest) .ne. 0)  return
c
c     compute linear velocity of the system center of mass
c
      if (integrate .eq. 'RIGIDBODY') then
         call rgdrest
      else
         call xyzrest
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine xyzrest  --  remove Cartesian system inertia  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "xyzrest" removes any translational or rotational inertia
c     during molecular dynamics for the overall system or atom groups
c
c
      subroutine xyzrest
      use atomid
      use atoms
      use bound
      use group
      use inform
      use iounit
      use moldyn
      use units
      implicit none
      integer i,j,k,m
      real*8 weigh,eps
      real*8 xx,yy,zz,xy,xz,yz
      real*8 xdel,ydel,zdel
      real*8 mang(3),tensor(3,3)
      real*8, allocatable :: totmass(:)
      real*8, allocatable :: etrans(:)
      real*8, allocatable :: erot(:)
      real*8, allocatable :: xtot(:)
      real*8, allocatable :: ytot(:)
      real*8, allocatable :: ztot(:)
      real*8, allocatable :: vtot(:,:)
      real*8, allocatable :: vang(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (totmass(0:ngrp))
      allocate (etrans(0:ngrp))
      allocate (vtot(3,0:ngrp))
      if (.not.use_bounds .or. ngrp.ne.0) then
         allocate (erot(0:ngrp))
         allocate (xtot(0:ngrp))
         allocate (ytot(0:ngrp))
         allocate (ztot(0:ngrp))
         allocate (vang(3,0:ngrp))
      end if
c
c     zero out total mass and linear velocity of each group
c
      do i = 0, ngrp
         totmass(i) = 0.0d0
         do j = 1, 3
            vtot(j,i) = 0.0d0
         end do
      end do
c
c     compute linear velocity of each group center of mass
c
      do i = 0, ngrp
         do k = igrp(1,i), igrp(2,i)
            m = kgrp(k)
            weigh = mass(m)
            totmass = totmass + weigh
            do j = 1, 3
               vtot(j,i) = vtot(j,i) + v(j,m)*weigh
            end do
         end do
      end do
c
c     compute translational kinetic energy of each group
c
      do i = 0, ngrp
         etrans(i) = 0.0d0
         do j = 1, 3
            vtot(j,i) = vtot(j,i) / totmass(i)
            etrans(i) = etrans(i) + vtot(j,i)**2
         end do
         etrans(i) = 0.5d0 * etrans(i) * totmass(i) / ekcal
      end do
c
c     find the center of mass coordinates of each atom group
c
      if (.not.use_bounds .or. ngrp.ne.0) then
         do i = 0, ngrp
            xtot(i) = 0.0d0
            ytot(i) = 0.0d0
            ztot(i) = 0.0d0
            do k = igrp(1,i), igrp(2,i)
               m = kgrp(k)
               weigh = mass(m)
               xtot(i) = xtot(i) + x(m)*weigh
               ytot(i) = ytot(i) + y(m)*weigh
               ztot(i) = ztot(i) + z(m)*weigh
            end do
            xtot(i) = xtot(i) / totmass(i)
            ytot(i) = ytot(i) / totmass(i)
            ztot(i) = ztot(i) / totmass(i)
c
c     compute the angular momentum of each atom group
c
            do j = 1, 3
               mang(j) = 0.0d0
            end do
            do k = igrp(1,i), igrp(2,i)
               m = kgrp(k)
               weigh = mass(m)
               mang(1) = mang(1) + (y(m)*v(3,m)-z(m)*v(2,m))*weigh
               mang(2) = mang(2) + (z(m)*v(1,m)-x(m)*v(3,m))*weigh
               mang(3) = mang(3) + (x(m)*v(2,m)-y(m)*v(1,m))*weigh
            end do
            mang(1) = mang(1) - (ytot(i)*vtot(3,i)-ztot(i)*vtot(2,i))
     &                             *totmass(i)
            mang(2) = mang(2) - (ztot(i)*vtot(1,i)-xtot(i)*vtot(3,i))
     &                             *totmass(i)
            mang(3) = mang(3) - (xtot(i)*vtot(2,i)-ytot(i)*vtot(1,i))
     &                             *totmass(i)
c
c     calculate the moment of inertia tensor
c
            xx = 0.0d0
            xy = 0.0d0
            xz = 0.0d0
            yy = 0.0d0
            yz = 0.0d0
            zz = 0.0d0
            do k = igrp(1,i), igrp(2,i)
               m = kgrp(k)
               weigh = mass(m)
               xdel = x(m) - xtot(i)
               ydel = y(m) - ytot(i)
               zdel = z(m) - ztot(i)
               xx = xx + xdel*xdel*weigh
               xy = xy + xdel*ydel*weigh
               xz = xz + xdel*zdel*weigh
               yy = yy + ydel*ydel*weigh
               yz = yz + ydel*zdel*weigh
               zz = zz + zdel*zdel*weigh
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
c
c     fix to avoid singularity for one- or two-body groups
c
            if (igrp(2,i)-igrp(1,i) .le. 2) then
               eps = 0.000001d0
               tensor(1,1) = tensor(1,1) + eps
               tensor(2,2) = tensor(2,2) + eps
               tensor(3,3) = tensor(3,3) + eps
            end if
c
c     diagonalize the moment of inertia tensor
c
            call invert (3,tensor)
c
c     compute angular velocity and rotational kinetic energy
c
            erot(i) = 0.0d0
            do k = 1, 3
               vang(k,i) = 0.0d0
               do j = 1, 3
                  vang(k,i) = vang(k,i) + tensor(k,j)*mang(j)
               end do
               erot(i) = erot(i) + vang(k,i)*mang(k)
            end do
            erot(i) = 0.5d0 * erot(i) / ekcal
         end do
      end if
c
c     eliminate any translation of each atom group
c
      do i = 0, ngrp
         do k = igrp(1,i), igrp(2,i)
            m = kgrp(k)
            do j = 1, 3
               v(j,m) = v(j,m) - vtot(j,i)
            end do
         end do
      end do
c
c     print the translational velocity of each atom group
c
      if (debug) then
         write (iout,10)
   10    format ()
         if (ngrp .eq. 0) then
            write (iout,20)  (vtot(i,0),i=1,3),etrans(0)
   20       format (' System Linear Velocity :  ',3d12.2,
     &              /,' Translational Kinetic Energy :',10x,f12.4,
     &                 ' Kcal/mole')
         else
            do i = 0, ngrp
               write (iout,30)  i,(vtot(j,i),j=1,3),etrans(i)
   30          format (' Group',i4,' Linear Velocity :  ',3d12.2,
     &                 /,' Translational Kinetic Energy :',10x,f12.4,
     &                    ' Kcal/mole')
            end do
         end if
      end if
c
c     eliminate any rotation about each group center of mass
c
      if (.not.use_bounds .or. ngrp.ne.0) then
         do i = 0, ngrp
            do k = igrp(1,i), igrp(2,i)
               m = kgrp(k)
               xdel = x(m) - xtot(i)
               ydel = y(m) - ytot(i)
               zdel = z(m) - ztot(i)
               v(1,m) = v(1,m) - vang(2,i)*zdel + vang(3,i)*ydel
               v(2,m) = v(2,m) - vang(3,i)*xdel + vang(1,i)*zdel
               v(3,m) = v(3,m) - vang(1,i)*ydel + vang(2,i)*xdel
            end do
         end do
c
c     print the angular velocity of each atom group
c
         if (debug) then
            if (ngrp .eq. 0) then
               write (iout,40)  (vang(j,0),j=1,3),erot(0)
   40          format (' System Angular Velocity : ',3d12.2,
     &                 /,' Rotational Kinetic Energy :',13x,f12.4,
     &                    ' Kcal/mole')
            else
               do i = 0, ngrp
                  write (iout,50)  i,(vang(j,i),j=1,3),erot(i)
   50             format (' Group',i4,' Angular Velocity : ',3d12.2,
     &                    /,' Rotational Kinetic Energy :',13x,f12.4,
     &                       ' Kcal/mole')
               end do
            end if
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (totmass)
      deallocate (etrans)
      deallocate (vtot)
      if (.not.use_bounds .or. ngrp.ne.0) then
         deallocate (erot)
         deallocate (xtot)
         deallocate (ytot)
         deallocate (ztot)
         deallocate (vang)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine rgdrest  --  remove rigidbody system inertia  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "rgdrest" removes any translational or rotational inertia
c     during molecular dynamics over rigid body coordinates
c
c
      subroutine rgdrest
      use atomid
      use atoms
      use bound
      use group
      use inform
      use iounit
      use rgddyn
      use units
      implicit none
      integer i,j,k
      real*8 etrans,erot
      real*8 weigh,totmass,eps
      real*8 xx,yy,zz,xy,xz,yz
      real*8 xtot,ytot,ztot
      real*8 xdel,ydel,zdel
      real*8 mang(3),vang(3)
      real*8 vtot(3),tensor(3,3)
      real*8, allocatable :: xcm(:)
      real*8, allocatable :: ycm(:)
      real*8, allocatable :: zcm(:)
c
c
c     zero out the total mass and overall linear velocity
c
      totmass = 0.0d0
      do j = 1, 3
         vtot(j) = 0.0d0
      end do
c
c     compute linear velocity of the system center of mass
c
      do i = 1, ngrp
         weigh = grpmass(i)
         totmass = totmass + weigh
         do j = 1, 3
            vtot(j) = vtot(j) + vcm(j,i)*weigh
         end do
      end do
c
c     compute translational kinetic energy of overall system
c
      etrans = 0.0d0
      do j = 1, 3
         vtot(j) = vtot(j) / totmass
         etrans = etrans + vtot(j)**2
      end do
      etrans = 0.5d0 * etrans * totmass / ekcal
c
c     perform dynamic allocation of some local arrays
c
      if (.not. use_bounds) then
         allocate (xcm(ngrp))
         allocate (ycm(ngrp))
         allocate (zcm(ngrp))
      end if
c
c     find the center of mass coordinates of the overall system
c
      if (.not. use_bounds) then
         xtot = 0.0d0
         ytot = 0.0d0
         ztot = 0.0d0
         do i = 1, ngrp
            xcm(i) = 0.0d0
            ycm(i) = 0.0d0
            zcm(i) = 0.0d0
            do j = igrp(1,i), igrp(2,i)
               k = kgrp(j)
               weigh = mass(k)
               xcm(i) = xcm(i) + x(k)*weigh
               ycm(i) = ycm(i) + y(k)*weigh
               zcm(i) = zcm(i) + z(k)*weigh
            end do
            xtot = xtot + xcm(i)
            ytot = ytot + ycm(i)
            ztot = ztot + zcm(i)
            weigh = max(1.0d0,grpmass(i))
            xcm(i) = xcm(i) / weigh
            ycm(i) = ycm(i) / weigh
            zcm(i) = zcm(i) / weigh
         end do
         xtot = xtot / totmass
         ytot = ytot / totmass
         ztot = ztot / totmass
c
c     compute the angular momentum of the overall system
c
         do j = 1, 3
            mang(j) = 0.0d0
         end do
         do i = 1, ngrp
            weigh = grpmass(i)
            mang(1) = mang(1) + (ycm(i)*vcm(3,i)
     &                          -zcm(i)*vcm(2,i))*weigh
            mang(2) = mang(2) + (zcm(i)*vcm(1,i)
     &                          -xcm(i)*vcm(3,i))*weigh
            mang(3) = mang(3) + (xcm(i)*vcm(2,i)
     &                          -ycm(i)*vcm(1,i))*weigh
         end do
         mang(1) = mang(1) - (ytot*vtot(3)-ztot*vtot(2))*totmass
         mang(2) = mang(2) - (ztot*vtot(1)-xtot*vtot(3))*totmass
         mang(3) = mang(3) - (xtot*vtot(2)-ytot*vtot(1))*totmass
c
c     calculate the moment of inertia tensor
c
         xx = 0.0d0
         xy = 0.0d0
         xz = 0.0d0
         yy = 0.0d0
         yz = 0.0d0
         zz = 0.0d0
         do i = 1, ngrp
            weigh = grpmass(i)
            xdel = xcm(i) - xtot
            ydel = ycm(i) - ytot
            zdel = zcm(i) - ztot
            xx = xx + xdel*xdel*weigh
            xy = xy + xdel*ydel*weigh
            xz = xz + xdel*zdel*weigh
            yy = yy + ydel*ydel*weigh
            yz = yz + ydel*zdel*weigh
            zz = zz + zdel*zdel*weigh
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
c
c     fix to avoid singularity for one- or two-body systems
c
         if (ngrp .le. 2) then
            eps = 0.000001d0
            tensor(1,1) = tensor(1,1) + eps
            tensor(2,2) = tensor(2,2) + eps
            tensor(3,3) = tensor(3,3) + eps
         end if
c
c     diagonalize the moment of inertia tensor
c
         call invert (3,tensor)
c
c     compute angular velocity and rotational kinetic energy
c
         erot = 0.0d0
         do i = 1, 3
            vang(i) = 0.0d0
            do j = 1, 3
               vang(i) = vang(i) + tensor(i,j)*mang(j)
            end do
            erot = erot + vang(i)*mang(i)
         end do
         erot = 0.5d0 * erot / ekcal
      end if
c
c     eliminate any translation of the overall system
c
      do i = 1, ngrp
         do j = 1, 3
            vcm(j,i) = vcm(j,i) - vtot(j)
         end do
      end do
c
c     print the translational velocity of the overall system
c
      if (debug) then
         write (iout,10)  (vtot(j),j=1,3),etrans
   10    format (' System Linear Velocity :  ',3d12.2,
     &           /,' Translational Kinetic Energy :',10x,f12.4,
     &              ' Kcal/mole')
      end if
c
c     eliminate any rotation about the system center of mass
c
      if (.not. use_bounds) then
         do i = 1, ngrp
            xdel = xcm(i) - xtot
            ydel = ycm(i) - ytot
            zdel = zcm(i) - ztot
            vcm(1,i) = vcm(1,i) - vang(2)*zdel + vang(3)*ydel
            vcm(2,i) = vcm(2,i) - vang(3)*xdel + vang(1)*zdel
            vcm(3,i) = vcm(3,i) - vang(1)*ydel + vang(2)*xdel
         end do
c
c     print the angular velocity of the overall system
c
         if (debug) then
            write (iout,20)  (vang(j),j=1,3),erot
   20       format (' System Angular Velocity : ',3d12.2,
     &              /,' Rotational Kinetic Energy :',13x,f12.4,
     &                 ' Kcal/mole')
         end if
      end if
c
c     perform deallocation of some local arrays
c
      if (.not. use_bounds) then
         deallocate (xcm)
         deallocate (ycm)
         deallocate (zcm)
      end if
      return
      end

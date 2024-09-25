c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2000  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine moments  --  total electric multipole moments  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "moments" computes the total electric charge, dipole components
c     and quadrupole components as a sum over the partial charges,
c     bond dipoles and atomic multipole moments over active atoms or
c     the full system
c
c     literature reference:
c
c     C. Gray and K. E. Gubbins, "Theory of Molecular Fluids, Volume 1:
c     Fundamentals", Oxford University Press, (1984)  [factor of 3/2 in
c     conversion of traced to traceless quadrupoles; pages 50-51]
c
c
      subroutine moments (mode)
      use atomid
      use atoms
      use bound
      use charge
      use dipole
      use limits
      use moment
      use mpole
      use polar
      use potent
      use rigid
      use solpot
      use units
      use usage
      implicit none
      integer i,j,k
      real*8 xc,yc,zc
      real*8 xi,yi,zi,ri
      real*8 weigh,trace
      real*8 xmid,ymid,zmid
      real*8 xbnd,ybnd,zbnd
      real*8, allocatable :: xcm(:)
      real*8, allocatable :: ycm(:)
      real*8, allocatable :: zcm(:)
      real*8 a(3,3),b(3,3)
      logical, allocatable :: temp(:)
      character*6 mode
c
c
c     zero out total charge, dipole and quadrupole components
c
      netchg = 0.0d0
      netdpl = 0.0d0
      netqpl(1) = 0.0d0
      netqpl(2) = 0.0d0
      netqpl(3) = 0.0d0
      xdpl = 0.0d0
      ydpl = 0.0d0
      zdpl = 0.0d0
      xxqpl = 0.0d0
      xyqpl = 0.0d0
      xzqpl = 0.0d0
      yxqpl = 0.0d0
      yyqpl = 0.0d0
      yzqpl = 0.0d0
      zxqpl = 0.0d0
      zyqpl = 0.0d0
      zzqpl = 0.0d0
c
c     perform dynamic allocation of some local arrays
c
      if (mode .eq. 'FULL')  allocate (temp(n))
c
c     store active atom list, and make all atoms active
c
      if (mode.eq.'FULL' .and. nuse.ne.n) then
         do i = 1, n
            temp(i) = use(i)
            use(i) = .true.
         end do
      end if
c
c     maintain periodic boundaries and neighbor lists
c
      if (use_bounds .and. .not.use_rigid)  call bounds
      if (use_clist .or. use_mlist)  call nblist
c
c     perform dynamic allocation of some local arrays
c
      allocate (xcm(n))
      allocate (ycm(n))
      allocate (zcm(n))
c
c     find the center of mass of the set of active atoms
c
      weigh = 0.0d0
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      do i = 1, n
         if (use(i)) then
            weigh = weigh + mass(i)
            xmid = xmid + x(i)*mass(i)
            ymid = ymid + y(i)*mass(i)
            zmid = zmid + z(i)*mass(i)
         end if
      end do
      if (weigh .ne. 0.0d0) then
         xmid = xmid / weigh
         ymid = ymid / weigh
         zmid = zmid / weigh
      end if
      do i = 1, n
         xcm(i) = x(i) - xmid
         ycm(i) = y(i) - ymid
         zcm(i) = z(i) - zmid
      end do
c
c     alter partial charges and monopoles via charge flux
c
      if (use_chgflx)  call alterchg
c
c     set the multipole moment components due to partial charges
c
      do i = 1, nion
         k = iion(i)
         if (use(k) .and. momuse(k)) then
            netchg = netchg + pchg(k)
            xdpl = xdpl + xcm(k)*pchg(k)
            ydpl = ydpl + ycm(k)*pchg(k)
            zdpl = zdpl + zcm(k)*pchg(k)
            xxqpl = xxqpl + xcm(k)*xcm(k)*pchg(k)
            xyqpl = xyqpl + xcm(k)*ycm(k)*pchg(k)
            xzqpl = xzqpl + xcm(k)*zcm(k)*pchg(k)
            yxqpl = yxqpl + ycm(k)*xcm(k)*pchg(k)
            yyqpl = yyqpl + ycm(k)*ycm(k)*pchg(k)
            yzqpl = yzqpl + ycm(k)*zcm(k)*pchg(k)
            zxqpl = zxqpl + zcm(k)*xcm(k)*pchg(k)
            zyqpl = zyqpl + zcm(k)*ycm(k)*pchg(k)
            zzqpl = zzqpl + zcm(k)*zcm(k)*pchg(k)
         end if
      end do
c
c     set the multipole moment components due to bond dipoles
c
      do i = 1, ndipole
         j = idpl(1,i)
         k = idpl(2,i)
         if ((use(j).and.momuse(j)) .or. (use(k).and.momuse(k))) then
            xi = x(j) - x(k)
            yi = y(j) - y(k)
            zi = z(j) - z(k)
            ri = sqrt(xi*xi + yi*yi + zi*zi)
            xbnd = bdpl(i) * (xi/ri) / debye
            ybnd = bdpl(i) * (yi/ri) / debye
            zbnd = bdpl(i) * (zi/ri) / debye
            xc = x(j) - xi*sdpl(i)
            yc = y(j) - yi*sdpl(i)
            zc = z(j) - zi*sdpl(i)
            xdpl = xdpl + xbnd
            ydpl = ydpl + ybnd
            zdpl = zdpl + zbnd
            xxqpl = xxqpl + 2.0d0*xc*xbnd
            xyqpl = xyqpl + xc*ybnd + yc*xbnd
            xzqpl = xzqpl + xc*zbnd + zc*xbnd
            yxqpl = yxqpl + yc*xbnd + xc*ybnd
            yyqpl = yyqpl + 2.0d0*yc*ybnd
            yzqpl = yzqpl + yc*zbnd + zc*ybnd
            zxqpl = zxqpl + zc*xbnd + xc*zbnd
            zyqpl = zyqpl + zc*ybnd + yc*zbnd
            zzqpl = zzqpl + 2.0d0*zc*zbnd
         end if
      end do
c
c     find atomic multipoles and induced dipoles in global frame
c
      if (use_born)  call born
      call chkpole
      call rotpole ('MPOLE')
      call induce
      if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
         do i = 1, npole
            k = ipole(i)
            rpole(2,k) = rpole(2,k) + uinds(1,k)
            rpole(3,k) = rpole(3,k) + uinds(2,k)
            rpole(4,k) = rpole(4,k) + uinds(3,k)
         end do
      else
         do i = 1, npole
            k = ipole(i)
            rpole(2,k) = rpole(2,k) + uind(1,k)
            rpole(3,k) = rpole(3,k) + uind(2,k)
            rpole(4,k) = rpole(4,k) + uind(3,k)
         end do
      end if
c
c     set the moment components due to atomic monopoles and dipoles
c
      do i = 1, npole
         k = ipole(i)
         if (use(k) .and. momuse(k)) then
            netchg = netchg + rpole(1,k)
            xdpl = xdpl + xcm(k)*rpole(1,k) + rpole(2,k)
            ydpl = ydpl + ycm(k)*rpole(1,k) + rpole(3,k)
            zdpl = zdpl + zcm(k)*rpole(1,k) + rpole(4,k)
            xxqpl = xxqpl + xcm(k)*xcm(k)*rpole(1,k)
     &                 + 2.0d0*xcm(k)*rpole(2,k)
            xyqpl = xyqpl + xcm(k)*ycm(k)*rpole(1,k)
     &                 + xcm(k)*rpole(3,k) + ycm(k)*rpole(2,k)
            xzqpl = xzqpl + xcm(k)*zcm(k)*rpole(1,k)
     &                 + xcm(k)*rpole(4,k) + zcm(k)*rpole(2,k)
            yxqpl = yxqpl + ycm(k)*xcm(k)*rpole(1,k)
     &                 + ycm(k)*rpole(2,k) + xcm(k)*rpole(3,k)
            yyqpl = yyqpl + ycm(k)*ycm(k)*rpole(1,k)
     &                 + 2.0d0*ycm(k)*rpole(3,k)
            yzqpl = yzqpl + ycm(k)*zcm(k)*rpole(1,k)
     &                 + ycm(k)*rpole(4,k) + zcm(k)*rpole(3,k)
            zxqpl = zxqpl + zcm(k)*xcm(k)*rpole(1,k)
     &                 + zcm(k)*rpole(2,k) + xcm(k)*rpole(4,k)
            zyqpl = zyqpl + zcm(k)*ycm(k)*rpole(1,k)
     &                 + zcm(k)*rpole(3,k) + ycm(k)*rpole(4,k)
            zzqpl = zzqpl + zcm(k)*zcm(k)*rpole(1,k)
     &                 + 2.0d0*zcm(k)*rpole(4,k)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xcm)
      deallocate (ycm)
      deallocate (zcm)
c
c     convert the quadrupole from traced to traceless form
c
      trace = (xxqpl + yyqpl + zzqpl) / 3.0d0
      xxqpl = 1.5d0 * (xxqpl-trace)
      xyqpl = 1.5d0 * xyqpl
      xzqpl = 1.5d0 * xzqpl
      yxqpl = 1.5d0 * yxqpl
      yyqpl = 1.5d0 * (yyqpl-trace)
      yzqpl = 1.5d0 * yzqpl
      zxqpl = 1.5d0 * zxqpl
      zyqpl = 1.5d0 * zyqpl
      zzqpl = 1.5d0 * (zzqpl-trace)
c
c     add the traceless atomic quadrupoles to total quadrupole
c
      do i = 1, npole
         k = ipole(i)
         if (use(k) .and. momuse(k)) then
            xxqpl = xxqpl + 3.0d0*rpole(5,k)
            xyqpl = xyqpl + 3.0d0*rpole(6,k)
            xzqpl = xzqpl + 3.0d0*rpole(7,k)
            yxqpl = yxqpl + 3.0d0*rpole(8,k)
            yyqpl = yyqpl + 3.0d0*rpole(9,k)
            yzqpl = yzqpl + 3.0d0*rpole(10,k)
            zxqpl = zxqpl + 3.0d0*rpole(11,k)
            zyqpl = zyqpl + 3.0d0*rpole(12,k)
            zzqpl = zzqpl + 3.0d0*rpole(13,k)
         end if
      end do
c
c     revert to the original set of active atoms
c
      if (mode.eq.'FULL' .and. nuse.ne.n) then
         do i = 1, n
            use(i) = temp(i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      if (mode .eq. 'FULL')  deallocate (temp)
c
c     convert dipole to Debye and quadrupole to Buckingham
c
      xdpl = xdpl * debye
      ydpl = ydpl * debye
      zdpl = zdpl * debye
      xxqpl = xxqpl * debye
      xyqpl = xyqpl * debye
      xzqpl = xzqpl * debye
      yxqpl = yxqpl * debye
      yyqpl = yyqpl * debye
      yzqpl = yzqpl * debye
      zxqpl = zxqpl * debye
      zyqpl = zyqpl * debye
      zzqpl = zzqpl * debye
c
c     get dipole magnitude and diagonalize quadrupole tensor
c
      netdpl = sqrt(xdpl*xdpl + ydpl*ydpl + zdpl*zdpl)
      a(1,1) = xxqpl
      a(1,2) = xyqpl
      a(1,3) = xzqpl
      a(2,1) = yxqpl
      a(2,2) = yyqpl
      a(2,3) = yzqpl
      a(3,1) = zxqpl
      a(3,2) = zyqpl
      a(3,3) = zzqpl
      call jacobi (3,a,netqpl,b)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine dmoments  --  total electric dipole moments  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "dmoments" computes the total dipole moments over all atoms and
c     the unique atom type dipole moments;
c     called in mdsave, it is assumed bound is called
c
c
      subroutine dmoments (xustc,yustc,zustc,xuind,yuind,zuind)
      use atomid
      use atoms
      use moment
      use mpole
      use polar
      use potent
      use uatom
      use units
      implicit none
      integer i,j
      integer ut
      real*8 xmid,ymid,zmid,weigh
      real*8 xu,yu,zu
      real*8 xustc,yustc,zustc
      real*8 xuind,yuind,zuind
c
c
c     find the center of mass of the set of active atoms
c
      weigh = 0.0d0
      xmid = 0.0d0
      ymid = 0.0d0
      zmid = 0.0d0
      do i = 1, n
         weigh = weigh + mass(i)
         xmid = xmid + x(i)*mass(i)
         ymid = ymid + y(i)*mass(i)
         zmid = zmid + z(i)*mass(i)
      end do
      if (weigh .ne. 0.0d0) then
         xmid = xmid / weigh
         ymid = ymid / weigh
         zmid = zmid / weigh
      end if
c
c     zero out dipole moments
c
      xustc = 0.0d0
      yustc = 0.0d0
      zustc = 0.0d0
      xuind = 0.0d0
      yuind = 0.0d0
      zuind = 0.0d0
      do i = 1, nunique
         do j = 1, 3
            utv1(j,i) = 0.0d0
            utv2(j,i) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(n,x,y,z,rpole,uind,use_polar,momuse,xmid,ymid,zmid,
!$OMP& xustc,yustc,zustc,xuind,yuind,zuind,utv1,utv2,type,utypeinv)
!$OMP DO reduction(+:xustc,yustc,zustc,utv1) schedule(guided)
c
c     compute the static dipole moment
c
      do i = 1, n
         if (momuse(i)) then
            xu = (x(i)-xmid)*rpole(1,i) + rpole(2,i)
            yu = (y(i)-ymid)*rpole(1,i) + rpole(3,i)
            zu = (z(i)-zmid)*rpole(1,i) + rpole(4,i)
            xustc = xustc + xu
            yustc = yustc + yu
            zustc = zustc + zu
            ut = utypeinv(type(i))
            utv1(1,ut) = utv1(1,ut) + xu
            utv1(2,ut) = utv1(2,ut) + yu
            utv1(3,ut) = utv1(3,ut) + zu
         end if
      end do
!$OMP END DO
c
c     compute the induced dipole moment
c
      if (use_polar) then
!$OMP DO reduction(+:xuind,yuind,zuind,utv2) schedule(guided)
         do i = 1, n
            if (momuse(i)) then
               xu = uind(1,i)
               yu = uind(2,i)
               zu = uind(3,i)
               xuind = xuind + xu
               yuind = yuind + yu
               zuind = zuind + zu
               ut = utypeinv(type(i))
               utv2(1,ut) = utv2(1,ut) + xu
               utv2(2,ut) = utv2(2,ut) + yu
               utv2(3,ut) = utv2(3,ut) + zu
            end if
         end do
!$OMP END DO
      end if
!$OMP END PARALLEL
c
c     convert dipole to Debye
c
      xustc = xustc * debye
      yustc = yustc * debye
      zustc = zustc * debye
      xuind = xuind * debye
      yuind = yuind * debye
      zuind = zuind * debye
      do i = 1, nunique
         do j = 1, 3
            utv1(j,i) = utv1(j,i) * debye
            utv2(j,i) = utv2(j,i) * debye
         end do
      end do
      return
      end

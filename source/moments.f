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
c     "moments" computes the total electric charge, dipole and
c     quadrupole moments for the active atoms as a sum over the
c     partial charges, bond dipoles and atomic multipole moments
c
c
      subroutine moments
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
      real*8 weigh,qave
      real*8 xc,yc,zc
      real*8 xi,yi,zi,ri
      real*8 xmid,ymid,zmid
      real*8 xbnd,ybnd,zbnd
      real*8, allocatable :: xcm(:)
      real*8, allocatable :: ycm(:)
      real*8, allocatable :: zcm(:)
      real*8 a(3,3),b(3,3)
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
c     set the multipole moment components due to partial charges
c
      do i = 1, nion
         k = iion(i)
         if (use(k)) then
            netchg = netchg + pchg(i)
            xdpl = xdpl + xcm(k)*pchg(i)
            ydpl = ydpl + ycm(k)*pchg(i)
            zdpl = zdpl + zcm(k)*pchg(i)
            xxqpl = xxqpl + xcm(k)*xcm(k)*pchg(i)
            xyqpl = xyqpl + xcm(k)*ycm(k)*pchg(i)
            xzqpl = xzqpl + xcm(k)*zcm(k)*pchg(i)
            yxqpl = yxqpl + ycm(k)*xcm(k)*pchg(i)
            yyqpl = yyqpl + ycm(k)*ycm(k)*pchg(i)
            yzqpl = yzqpl + ycm(k)*zcm(k)*pchg(i)
            zxqpl = zxqpl + zcm(k)*xcm(k)*pchg(i)
            zyqpl = zyqpl + zcm(k)*ycm(k)*pchg(i)
            zzqpl = zzqpl + zcm(k)*zcm(k)*pchg(i)
         end if
      end do
c
c     set the multipole moment components due to bond dipoles
c
      do i = 1, ndipole
         j = idpl(1,i)
         k = idpl(2,i)
         if (use(j) .or. use(k)) then
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
      call rotpole
      call induce
      if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
         do i = 1, npole
            rpole(2,i) = rpole(2,i) + uinds(1,i)
            rpole(3,i) = rpole(3,i) + uinds(2,i)
            rpole(4,i) = rpole(4,i) + uinds(3,i)
         end do
      else
         do i = 1, npole
            rpole(2,i) = rpole(2,i) + uind(1,i)
            rpole(3,i) = rpole(3,i) + uind(2,i)
            rpole(4,i) = rpole(4,i) + uind(3,i)
         end do
      end if
c
c     set the multipole moment components due to atomic multipoles
c
      do i = 1, npole
         k = ipole(i)
         if (use(k)) then
            netchg = netchg + rpole(1,i)
            xdpl = xdpl + xcm(k)*rpole(1,i) + rpole(2,i)
            ydpl = ydpl + ycm(k)*rpole(1,i) + rpole(3,i)
            zdpl = zdpl + zcm(k)*rpole(1,i) + rpole(4,i)
            xxqpl = xxqpl + xcm(k)*xcm(k)*rpole(1,i)
     &                 + 2.0d0*xcm(k)*rpole(2,i)
            xyqpl = xyqpl + xcm(k)*ycm(k)*rpole(1,i)
     &                 + xcm(k)*rpole(3,i) + ycm(k)*rpole(2,i)
            xzqpl = xzqpl + xcm(k)*zcm(k)*rpole(1,i)
     &                 + xcm(k)*rpole(4,i) + zcm(k)*rpole(2,i)
            yxqpl = yxqpl + ycm(k)*xcm(k)*rpole(1,i)
     &                 + ycm(k)*rpole(2,i) + xcm(k)*rpole(3,i)
            yyqpl = yyqpl + ycm(k)*ycm(k)*rpole(1,i)
     &                 + 2.0d0*ycm(k)*rpole(3,i)
            yzqpl = yzqpl + ycm(k)*zcm(k)*rpole(1,i)
     &                 + ycm(k)*rpole(4,i) + zcm(k)*rpole(3,i)
            zxqpl = zxqpl + zcm(k)*xcm(k)*rpole(1,i)
     &                 + zcm(k)*rpole(2,i) + xcm(k)*rpole(4,i)
            zyqpl = zyqpl + zcm(k)*ycm(k)*rpole(1,i)
     &                 + zcm(k)*rpole(3,i) + ycm(k)*rpole(4,i)
            zzqpl = zzqpl + zcm(k)*zcm(k)*rpole(1,i)
     &                 + 2.0d0*zcm(k)*rpole(4,i)
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
      qave = (xxqpl + yyqpl + zzqpl) / 3.0d0
      xxqpl = 1.5d0 * (xxqpl-qave)
      xyqpl = 1.5d0 * xyqpl
      xzqpl = 1.5d0 * xzqpl
      yxqpl = 1.5d0 * yxqpl
      yyqpl = 1.5d0 * (yyqpl-qave)
      yzqpl = 1.5d0 * yzqpl
      zxqpl = 1.5d0 * zxqpl
      zyqpl = 1.5d0 * zyqpl
      zzqpl = 1.5d0 * (zzqpl-qave)
c
c     add the traceless atomic quadrupoles to total quadrupole
c
      do i = 1, npole
         k = ipole(i)
         if (use(k)) then
            xxqpl = xxqpl + 3.0d0*rpole(5,i)
            xyqpl = xyqpl + 3.0d0*rpole(6,i)
            xzqpl = xzqpl + 3.0d0*rpole(7,i)
            yxqpl = yxqpl + 3.0d0*rpole(8,i)
            yyqpl = yyqpl + 3.0d0*rpole(9,i)
            yzqpl = yzqpl + 3.0d0*rpole(10,i)
            zxqpl = zxqpl + 3.0d0*rpole(11,i)
            zyqpl = zyqpl + 3.0d0*rpole(12,i)
            zzqpl = zzqpl + 3.0d0*rpole(13,i)
         end if
      end do
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
c     #################################################################
c     ##                                                             ##
c     ##  subroutine momfull  --  multipole moments for full system  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "momfull" computes the electric moments for the full system
c     as a sum over the partial charges, bond dipoles and atomic
c     multipole moments
c
c
      subroutine momfull
      use atoms
      use usage
      implicit none
      integer i
      logical, allocatable :: temp(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (temp(n))
c
c     store active atom list, and make all atoms active
c
      if (nuse .ne. n) then
         do i = 1, n
            temp(i) = use(i)
            use(i) = .true.
         end do
      end if
c
c     compute the electric multipole moments for the system
c
      call moments
c
c     revert to the original set of active atoms
c
      if (nuse .ne. n) then
         do i = 1, n
            use(i) = temp(i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (temp)
      return
      end

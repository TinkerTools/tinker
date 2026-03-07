c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2003 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ktortor  --  tors-tors parameter assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ktortor" assigns torsion-torsion parameters to adjacent
c     torsion pairs and processes any new or changed values
c
c
      subroutine ktortor
      use atomid
      use atoms
      use bitor
      use inform
      use iounit
      use keys
      use ktrtor
      use potent
      use tortor
      implicit none
      integer i,j,k,m
      integer ia,ib,ic,id,ie
      integer ita,itb,itc,itd,ite
      integer size,next,ntt
      integer nx,ny,nxy
      integer tkey(maxtgrd2)
      real*8 eps
      real*8 tx(maxtgrd2)
      real*8 ty(maxtgrd2)
      real*8 tf(maxtgrd2)
      real*8 tind(maxtgrd2)
      real*8 bs(0:maxtgrd)
      real*8 cs(0:maxtgrd)
      real*8 ds(0:maxtgrd)
      real*8 tmp1(0:maxtgrd)
      real*8 tmp2(0:maxtgrd)
      real*8 tmp3(0:maxtgrd)
      real*8 tmp4(0:maxtgrd)
      real*8 tmp5(0:maxtgrd)
      real*8 tmp6(0:maxtgrd)
      real*8 tmp7(0:maxtgrd)
      logical header,cyclic
      character*3 ttag
      character*4 pa,pb,pc,pd,pe
      character*20 blank,pt
      character*20 pt1,pt2
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing torsion-torsion parameters
c
      blank = '                    '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'TORTORS ') then
            ia = 0
            ib = 0
            ic = 0
            id = 0
            ie = 0
            nx = 0
            ny = 0
            nxy = 0
            ttag = '   '
            do j = 1, maxtgrd2
               tx(j) = 0.0d0
               ty(j) = 0.0d0
               tf(j) = 0.0d0
            end do
            string = record(next:240)
            read (string,*,err=40,end=40)  ia,ib,ic,id,ie,nx,ny
            nxy = nx * ny
            call getword (record,ttag,next)
            m = i
            j = 0
            dowhile (j .lt. nxy)
               m = m + 1
               record = keyline(m)
               read (record,*,err=10,end=10)  tx(j+1),ty(j+1),tf(j+1),
     &                                        tx(j+2),ty(j+2),tf(j+2),
     &                                        tx(j+3),ty(j+3),tf(j+3)
               j = j + 3
               goto 30
   10          continue
               read (record,*,err=20,end=20)  tx(j+1),ty(j+1),tf(j+1),
     &                                        tx(j+2),ty(j+2),tf(j+2)
               j = j + 2
               goto 30
   20          continue
               read (record,*,err=40,end=40)  tx(j+1),ty(j+1),tf(j+1)
               j = j + 1
   30          continue
            end do
   40       continue
            if (.not. silent) then
               if (header) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Torsion-Torsion Parameters :',
     &                    //,5x,'Atom Classes',12x,'Grid-1',
     &                       6x,'Grid-2',6x,'Tier',/)
               end if
               write (iout,60)  ia,ib,ic,id,ie,nx,ny,ttag
   60          format (1x,5i4,5x,i8,4x,i8,8x,a3)
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            call numeral (id,pd,size)
            call numeral (ie,pe,size)
            pt = pa//pb//pc//pd//pe
            do j = 1, maxntt
               if (ktt(j).eq.blank .or. ktt(j).eq.pt) then
                  ktt(j) = pt
                  ttier(j) = ttag
                  do k = 1, nxy
                     tind(k) = 360.0d0*ty(k) + tx(k)
                     tkey(k) = k
                  end do
                  call sort2 (nxy,tind,tkey)
                  do k = 1, nxy
                     tbf(k,j) = tf(tkey(k))
                  end do
                  nx = nxy
                  call sort9 (nx,tx)
                  tnx(j) = nx
                  do k = 1, nx
                     ttx(k,j) = tx(k)
                  end do
                  ny = nxy
                  call sort9 (ny,ty)
                  tny(j) = ny
                  do k = 1, ny
                     tty(k,j) = ty(k)
                  end do
                  goto 80
               end if
            end do
            write (iout,70)
   70       format (/,' KTORTOR  --  Too many Torsion-Torsion',
     &                 ' Parameters')
            abort = .true.
   80       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      ntt = maxntt
      do i = maxntt, 1, -1
         if (ktt(i) .eq. blank)  ntt = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(itt))  deallocate (itt)
      allocate (itt(3,nbitor))
c
c     check whether each torsion-torsion parameter is periodic;
c     assumes the "tbf" array is sorted with both indices in
c     increasing order and the first index changing most rapidly
c
      do i = 1, ntt
         cyclic = .true.
         eps = 0.000001d0
         nx = tnx(i) - 1
         ny = tny(i) - 1
         if (abs(abs(ttx(1,i)-ttx(tnx(i),i))-360.0d0) .gt. eps)
     &      cyclic = .false.
         if (abs(abs(tty(1,i)-tty(tny(i),i))-360.0d0) .gt. eps)
     &      cyclic = .false.
         if (cyclic) then
            do j = 1, tny(i)
               k = (j-1)*tnx(i) + 1
               if (abs(tbf(k,i)-tbf(k+nx,i)) .gt. eps) then
                  write (iout,90)  tbf(k,i),tbf(k+nx,i)
   90             format (/,' KTORTOR  --  Warning, Unequal Tor-Tor',
     &                        ' Values',3x,2f12.5)
                  abort = .true.
               end if
            end do
            k = ny * tnx(i)
            do j = 1, tnx(i)
               if (abs(tbf(j,i)-tbf(j+k,i)) .gt. eps) then
                  write (iout,100)  tbf(j,i),tbf(j+k,i)
  100             format (/,' KTORTOR  --  Warning, Unequal Tor-Tor',
     &                        ' Values',3x,2f12.5)
                  abort = .true.
               end if
            end do
         end if
c
c     spline fit the derivatives about the first torsion
c
         do j = 1, tnx(i)
            tmp1(j-1) = ttx(j,i)
         end do
         m = 0
         do j = 1, tny(i)
            do k = 1, tnx(i)
               tmp2(k-1) = tbf(m+k,i)
            end do
            if (cyclic) then
               call cspline (nx,tmp1,tmp2,bs,cs,ds,tmp3,
     &                         tmp4,tmp5,tmp6,tmp7)
            else
               call nspline (nx,tmp1,tmp2,bs,cs,tmp3,
     &                         tmp4,tmp5,tmp6,tmp7)
            end if
            do k = 1, tnx(i)
               tbx(m+k,i) = bs(k-1)
            end do
            m = m + tnx(i)
         end do
c
c     spline fit the derivatives about the second torsion
c
         do j = 1, tny(i)
            tmp1(j-1) = tty(j,i)
         end do
         m = 1
         do j = 1, tnx(i)
            do k = 1, tny(i)
               tmp2(k-1) = tbf(m+(k-1)*tnx(i),i)
            end do
            if (cyclic) then
               call cspline (ny,tmp1,tmp2,bs,cs,ds,tmp3,
     &                         tmp4,tmp5,tmp6,tmp7)
            else
               call nspline (ny,tmp1,tmp2,bs,cs,tmp3,
     &                         tmp4,tmp5,tmp6,tmp7)
            end if
            do k = 1, tny(i)
               tby(m+(k-1)*tnx(i),i) = bs(k-1)
            end do
            m = m + 1
         end do
c
c     spline fit the cross derivatives about both torsions
c
         m = 1
         do j = 1, tnx(i)
            do k = 1, tny(i)
               tmp2(k-1) = tbx(m+(k-1)*tnx(i),i)
            end do
            if (cyclic) then
               call cspline (ny,tmp1,tmp2,bs,cs,ds,tmp3,
     &                          tmp4,tmp5,tmp6,tmp7)
            else
               call nspline (ny,tmp1,tmp2,bs,cs,tmp3,
     &                         tmp4,tmp5,tmp6,tmp7)
            end if
            do k = 1, tny(i)
               tbxy(m+(k-1)*tnx(i),i) = bs(k-1)
            end do
            m = m + 1
         end do
      end do
c
c     assign torsion-torsion parameters for each bitorsion
c
      ntortor = 0
      do i = 1, nbitor
         ia = ibitor(1,i)
         ib = ibitor(2,i)
         ic = ibitor(3,i)
         id = ibitor(4,i)
         ie = ibitor(5,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         itd = class(id)
         ite = class(ie)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         call numeral (itd,pd,size)
         call numeral (ite,pe,size)
         pt1 = pa//pb//pc//pd//pe
         pt2 = pe//pd//pc//pb//pa
c
c     assign tier-specific parameters for this torsion-torsion
c
         do j = 1, ntt
            if (ttier(j) .eq. tier(ic)) then
               if (ktt(j) .eq. pt1) then
                  ntortor = ntortor + 1
                  itt(1,ntortor) = i
                  itt(2,ntortor) = j
                  itt(3,ntortor) = 1
                  goto 110
               else if (ktt(j) .eq. pt2) then
                  ntortor = ntortor + 1
                  itt(1,ntortor) = i
                  itt(2,ntortor) = j
                  itt(3,ntortor) = -1
                  goto 110
               end if
            end if
         end do
c
c     assign nonspecific parameters for this torsion-torsion
c
         do j = 1, ntt
            if (ttier(j) .eq. '   ') then
               if (ktt(j) .eq. pt1) then
                  ntortor = ntortor + 1
                  itt(1,ntortor) = i
                  itt(2,ntortor) = j
                  itt(3,ntortor) = 1
                  goto 110
               else if (ktt(j) .eq. pt2) then
                  ntortor = ntortor + 1
                  itt(1,ntortor) = i
                  itt(2,ntortor) = j
                  itt(3,ntortor) = -1
                  goto 110
               end if
            end if
         end do
  110    continue
      end do
c
c     turn off the torsion-torsion potential if it is not used
c
      if (ntortor .eq. 0)  use_tortor = .false.
      return
      end

c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prtarc  --  output of a Tinker archive file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtarc" writes out a set of Cartesian coordinates for
c     all active atoms in the Tinker XYZ archive format
c
c
      subroutine prtarc (iarc)
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use titles
      use usage
      implicit none
      integer i,j,k,iarc
      integer size,crdsiz
      real*8 crdmin,crdmax
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
      character*240 arcfile
c
c
c     open output unit if not already done
c
      inquire (unit=iarc,opened=opened)
      if (.not. opened) then
         arcfile = filename(1:leng)//'.arc'
         call version (arcfile,'new')
         open (unit=iarc,file=arcfile,status='new')
      end if
c
c     check for large systems needing extended formatting
c
      atmc = 'i6'
      if (n .ge. 100000)  atmc = 'i7'
      if (n .ge. 1000000)  atmc = 'i8'
      crdmin = 0.0d0
      crdmax = 0.0d0
      do i = 1, n
         crdmin = min(crdmin,x(i),y(i),z(i))
         crdmax = max(crdmax,x(i),y(i),z(i))
      end do
      crdsiz = 6
      if (crdmin .le. -1000.0d0)  crdsiz = 7
      if (crdmax .ge. 10000.0d0)  crdsiz = 7
      if (crdmin .le. -10000.0d0)  crdsiz = 8
      if (crdmax .ge. 100000.0d0)  crdsiz = 8
      crdsiz = crdsiz + max(6,digits)
      size = 0
      call numeral (crdsiz,crdc,size)
      if (digits .le. 6) then
         digc = '6 '
      else if (digits .le. 8) then
         digc = '8'
      else
         digc = '10'
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         fstr = '('//atmc//')'
         write (iarc,fstr(1:4))  nuse
      else
         fstr = '('//atmc//',2x,a)'
         write (iarc,fstr(1:9))  nuse,title(1:ltitle)
      end if
c
c     write out the periodic cell lengths and angles
c
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (iarc,fstr)  xbox,ybox,zbox,alpha,beta,gamma
      end if
c
c     write out the coordinate line for each atom
c
      fstr = '('//atmc//',2x,a3,3f'//crdc//
     &          '.'//digc//',i6,8'//atmc//')'
      do i = 1, n
         if (use(i)) then
            k = n12(i)
            if (k .eq. 0) then
               write (iarc,fstr)  iuse(i),name(i),x(i),y(i),z(i),
     &                            type(i)
            else
               write (iarc,fstr)  iuse(i),name(i),x(i),y(i),z(i),
     &                            type(i),(iuse(i12(j,i)),j=1,k)
            end if
         end if
      end do
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=iarc)
      return
      end

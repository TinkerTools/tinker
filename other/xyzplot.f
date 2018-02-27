c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program xyzplot  --  plot "Rubenstein" model of structure  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "xyzplot" makes a plot of a Tinker Cartesian coordinates file
c     on a display terminal, Imagen LaserPrinter or HP 7475 plotter
c
c     this version is written to link against the ancient Unified
c     Graphics System library originally available from SLAC; the
c     code should be easily modified to use other graphics libraries
c     with simple line drawing capability
c
c
      program xyzplot
      implicit none
      include 'sizes.i'
      integer maxseg
      parameter (maxseg=1000)
      integer i,j,ixyz,ndraw
      integer segm,incr,trim
      integer nstruc,maxstruc
      integer last,next,freeunit
      real*8 xrot,yrot,zrot,rotval,scale
      real*8 bline,radian,xsize,ysize
      logical flag
      character*1 answer,axis
      character*2 label(maxatm)
      character*5 device
      character*20 timbuf
      character*80 record
      common /plot/ segm(maxseg)
c
c
c     get the user specified input structure
c
      call initial
      call getxyz
c
c     get the device type and set up its screen area
c
      call opendev (device)
c
c     initialize "radian" and "ndraw"; if using HP7550
c     then give the user a chance to increase "ndraw"
c
      radian = 180.0d0 / acos(-1.0d0)
      nstruc = 1
      ndraw = 1
      if (device .eq. 'HPP  ') then
         write (*,40)
   40    format (/,' Number of Overwrites [1] ?  ',$)
         read (*,50)  ndraw
   50    format (i)
         if (ndraw .eq. 0)  ndraw = 1
      end if
c
c     draw the initial structure
c
      write (*,60)
   60 format (/,' Draw Current Structure [Y] ?  ',$)
      read (*,70)  record
   70 format (a80)
      next = 1
      call gettext (record,answer,next)
      call upcase (answer)
      if (answer .ne. 'N') then
         call scale (xsize,ysize)
         call setug ('FULL ',xsize,ysize)
         call caption (ndraw,ysize)
         call molcule (nstruc,ndraw)
         if (device .eq. 'HPP  ')  call reset (device)
      end if
c     nstruc = 2
c
c     option to rotate to make desired plot
c
      answer = 'Y'
      dowhile (answer .eq. 'Y')
         write (*,80)
   80    format (/,' Rotate to New View [N] ?  ',$)
         read (*,90)  record
   90    format (a80)
         next = 1
         call gettext (record,answer,next)
         call upcase (answer)
         if (answer .eq. 'Y') then
            write (*,100)
  100       format (/,' Rotation Angles [X,Y,Z] ?  ',$)
            read (*,110)  xrot,yrot,zrot
  110       format (3f10.0)
            if (xrot .ne. 0.0)  call rotate ('X',xrot/radian)
            if (yrot .ne. 0.0)  call rotate ('Y',yrot/radian)
            if (zrot .ne. 0.0)  call rotate ('Z',zrot/radian)
            call scale (xsize,ysize)
            call setug ('FULL ',xsize,ysize)
            call caption (ndraw,ysize)
            call molcule (nstruc,ndraw)
            if (device .eq. 'HPP  ')  call reset (device)
            write (*,120)
  120       format (/,' Undo Last Rotation [N] ?  ',$)
            read (*,130)  record
  130       format (a80)
            next = 1
            call gettext (record,answer,next)
            call upcase (answer)
            if (answer .eq. 'Y') then
               if (zrot .ne. 0.0)  call rotate ('Z',-zrot/radian)
               if (yrot .ne. 0.0)  call rotate ('Y',-yrot/radian)
               if (xrot .ne. 0.0)  call rotate ('X',-xrot/radian)
            end if
         end if
      end do
c
c     option to make a full-sized stereo plot
c
      answer = 'Y'
      dowhile (answer .eq. 'Y')
         write (*,140)
  140    format (/,' Stereo View [N] ?  ',$)
         read (*,150)  record
  150    format (a80)
         next = 1
         call gettext (record,answer,next)
         call upcase (answer)
         if (answer .eq. 'Y') then
            write (*,160)
  160       format (' Rotation Angle [5.0] ?  ',$)
            read (*,170)  rotval
  170       format (f10.0)
            if (rotval .eq. 0.0)  rotval = 5.0
            write (*,180)
  180       format (/,' Rotation Axis, X, Y or Z [Y] ?  ',$)
            read (*,190)  record
  190       format (a80)
            next = 1
            call gettext (record,axis,next)
            call upcase (axis)
            if (axis.ne.'X' .and. axis.ne.'Z')  axis = 'Y'
            call rotate (axis,rotval*3.14159/180.0)
            call setug ('FULL ',xsize,ysize)
            call caption (ndraw,ysize)
            call molcule (nstruc,ndraw)
            call rotate (axis,-rotval*3.14159/180.0)
            if (device .eq. 'HPP  ')  call reset (device)
         end if
      end do
c
c     option to draw a single stereoscopic structure pair
c
      write (*,200)
  200 format (/,' Stereoscopic Pair [N] ?  ',$)
      read (*,210)  record
  210 format (a80)
      next = 1
      call gettext (record,answer,next)
      call upcase (answer)
      if (answer .eq. 'Y') then
         rotval = 5.0
         call scale (xsize,ysize)
         call setug ('FULL ',xsize,ysize)
         call caption (ndraw,ysize)
         call setug ('LEFT ',xsize,ysize)
         call molcule (nstruc,ndraw)
         call rotate ('Y',rotval*3.14159/180.0)
         call setug ('RIGHT',xsize,ysize)
         call molcule (nstruc,ndraw)
         call rotate ('Y',-rotval*3.14159/180.0)
         if (device .eq. 'HPP  ')  call reset (device)
      end if
c
c     close the graphics device
c
      call closdev (device)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################
c     ##                          ##
c     ##     subroutine scale     ##
c     ##                          ##
c     ##############################
c
c
      subroutine scale (xsize,ysize)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i
      real*8 xmin,ymin,zmin
      real*8 xmax,ymax,zmax
      real*8 xsize,ysize,zsize
      real*8 xcenter,ycenter,zcenter
      real*8 xdev,ydev,xydev,factor
      common /device/ xdev,ydev
c
c
c     find the max and min coordinate values; translate
c     the structure so (0,0,0) is at the center
c
      xmin = x(1)
      xmax = x(1)
      ymin = y(1)
      ymax = y(1)
      zmin = z(1)
      zmax = z(1)
      do i = 2, n
         if (x(i) .lt. xmin)  xmin = x(i)
         if (x(i) .gt. xmax)  xmax = x(i)
         if (y(i) .lt. ymin)  ymin = y(i)
         if (y(i) .gt. ymax)  ymax = y(i)
         if (z(i) .lt. zmin)  zmin = z(i)
         if (z(i) .gt. zmax)  zmax = z(i)
      end do
      xsize = xmax - xmin
      ysize = ymax - ymin
      zsize = zmax - zmin
      xcenter = 0.5 * (xmin+xmax)
      ycenter = 0.5 * (ymin+ymax)
      zcenter = 0.5 * (zmin+zmax)
      do i = 1, n
         x(i) = x(i) - xcenter
         y(i) = y(i) - ycenter
         z(i) = z(i) - zcenter
      end do
c
c     now do scaling of the structure size according
c     to the plotting area of the output device
c
      xdev = 11.0
      ydev = 8.5
      xydev = xdev / ydev
      if (xsize/ysize .ge. xydev) then
         ysize = xsize / xydev
      else
         xsize = ysize * xydev
      end if
      write (*,10)
   10 format (/,' Scale Factor [0.95] ?  ',$)
      read (*,20)  factor
   20 format (f10.0)
      if (factor .eq. 0.0)  factor = 0.95
      xsize = xsize / factor
      ysize = ysize / factor
      return
      end
c
c
c     ##############################
c     ##                          ##
c     ##     subroutine setug     ##
c     ##                          ##
c     ##############################
c
c
      subroutine setug (portion,xsize,ysize)
      implicit none
      integer maxseg
      parameter (maxseg=1000)
      integer segm
      real*8 xsize,ysize,xdev,ydev
      real*8 area(2,2),window(2,2)
      character*5 portion
      common /device/ xdev,ydev
      common /plot/ segm(maxseg)
c
c
c     call the "unified graphics system" routines
c
      if (portion .eq. 'FULL ') then
         call ugpict ('CLEAR',0)
         call ugdspc ('PUT',xdev,ydev,1.0)
      end if
      if (portion .eq. 'FULL ') then
         area(1,1) = 0.00
         area(1,2) = xdev
         area(2,1) = 0.00
         area(2,2) = ydev
      else
         if (portion .eq. 'LEFT ') then
            area(1,1) = 0.00 * xdev
            area(1,2) = 0.50 * xdev
            area(2,1) = 0.25 * ydev
            area(2,2) = 0.75 * ydev
         else
            area(1,1) = 0.50 * xdev
            area(1,2) = 1.00 * xdev
            area(2,1) = 0.25 * ydev
            area(2,2) = 0.75 * ydev
         end if
      end if
      window(1,1) = -0.5 * xsize
      window(1,2) =  0.5 * xsize
      window(2,1) = -0.5 * ysize
      window(2,2) =  0.5 * ysize
      call ugwdow ('PUT',area,window)
      call uginit ('CLEAR',segm,1000)
      return
      end
c
c
c     ################################
c     ##                            ##
c     ##     subroutine caption     ##
c     ##                            ##
c     ################################
c
c
      subroutine caption (ndraw,ysize)
      implicit none
      include 'sizes.i'
      include 'title.i'
      integer maxseg
      parameter (maxseg=1000)
      integer i,segm,ndraw,trim,textfile
      real*8 ysize,bline
      common /plot/ segm(maxseg)
c
c
      if (ltitle .ne. 0) then
         textfile = 1
         bline = -0.485 * ysize
         do i = 1, ndraw
            call ugtext ('DSIZE=0.020,CENTER,VBRIGHT,BLACK',
     &                      0.0,bline,title(1:ltitle),segm)
            call ugwrit (' ',textfile,segm)
            call uginit ('CLEAR',segm,1000)
         end do
      end if
      return
      end
c
c
c     ################################
c     ##                            ##
c     ##     subroutine molcule     ##
c     ##                            ##
c     ################################
c
c
      subroutine molcule (nstruc,ndraw)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer maxseg
      parameter (maxseg=1000)
      integer i,j,nstruc,ndraw,segm
      common /plot/ segm(maxseg)
c
c
      do i = 1, ndraw
         do j = 1, n
            call atom (j)
            call ugwrit (' ',nstruc,segm)
            call uginit ('CLEAR',segm,1000)
         end do
      end do
      return
      end
c
c
c
c     ###############################
c     ##                           ##
c     ##     subroutine rotate     ##
c     ##                           ##
c     ###############################
c
c
      subroutine rotate (axis,angle)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,j,iaxis
      real*8 angle,rot(3,3)
      character*1 axis
c
c
c     calculate rotation matrix; "iaxis" = 1,2,3 corresponds
c     to rotation about the x,y,z axes respectively
c
      if (axis .eq. 'X')  iaxis = 1
      if (axis .eq. 'Y')  iaxis = 2
      if (axis .eq. 'Z')  iaxis = 3
      do i = 1, 3
         do j = 1, 3
            if (i.ne.iaxis .and. j.ne.iaxis) then
               if (i .eq. j) then
                  rot(i,j) = cos(angle)
               else if (i .gt. j) then
                  rot(i,j) = sin(angle)
               else
                  rot(i,j) = -sin(angle)
               end if
            else
               rot(i,j) = 0.0
            end if
         end do
      end do
      rot(iaxis,iaxis) = 1.0
c
c     multiply each coordinate triplet by the rotation matrix
c
      do i = 1, n
         call matmul (i,rot)
      end do
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##     subroutine matmul     ##
c     ##                           ##
c     ###############################
c
c
      subroutine matmul (num,array)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,num
      real*8 temp(3),array(3,3)
c
c
      do i = 1, 3
         temp(i) = array(i,1)*x(num) + array(i,2)*y(num)
     &                + array(i,3)*z(num)
      end do
      x(num) = temp(1)
      y(num) = temp(2)
      z(num) = temp(3)
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##     subroutine atom     ##
c     ##                         ##
c     #############################
c
c
      subroutine atom (i)
      implicit none
      include 'sizes.i'
      include 'couple.i'
      integer i,j
c
c
      do j = 1, n12(i)
         if (i .lt. i12(j,i)) then
            call bond (i,i12(j,i))
         end if
      end do
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##     subroutine bond     ##
c     ##                         ##
c     #############################
c
c
      subroutine bond (i,j)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      integer maxseg
      parameter (maxseg=1000)
      integer i,j,k,segm
      real*8 depth
      character*3 name,element(2)
      character*30 option(2),temp
      common /plot/ segm(maxseg)
c
c
c     set options string for "ugline" by atom type
c
      element(1) = name(i)
      element(2) = name(j)
      option(1) = 'SOLID,VBRIGHT,BLACK'
      option(2) = 'SOLID,VBRIGHT,BLACK'
      do k = 1, 2
c        if (element(k) .eq. 'H  ')  element(k) = element(3-k)
c        if (element(k) .eq. 'C  ')  option(k) = 'SOLID,VBRIGHT,BLACK'
c        if (element(k) .eq. 'N  ')  option(k) = 'SOLID,VBRIGHT,BLUE'
c        if (element(k) .eq. 'N+ ')  option(k) = 'SOLID,VBRIGHT,BLUE'
c        if (element(k) .eq. 'O  ')  option(k) = 'SOLID,VBRIGHT,RED'
c        if (element(k) .eq. 'O- ')  option(k) = 'SOLID,VBRIGHT,RED'
c        if (element(k) .eq. 'H  ')  option(k) = 'SOLID,VBRIGHT,GREEN'
      end do
c
c     find forward or backward depth cueing for bond
c
      depth = z(j) - z(i)
c
c     draw bond as either straight line, dash or wedge
c
      if (abs(depth) .lt. 0.2) then
         call plain (i,j,option)
      else if (depth .gt. 0.0) then
         call wedge (i,j,option)
      else if (n12(j) .le. 1) then
         call dash (i,j,option)
      else
         temp = option(1)
         option(1) = option(2)
         option(2) = temp
         call wedge (j,i,option)
      end if
      return
      end
c
c
c     ########################
c     ##                    ##
c     ##  subroutine plain  ##
c     ##                    ##
c     ########################
c
c
      subroutine plain (i,j,option)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer maxseg
      parameter (maxseg=1000)
      integer i,j,segm
      real*8 xhalf,yhalf
      character*30 option(2)
      common /plot/ segm(maxseg)
c
c
      xhalf = 0.5 * (x(i)+x(j))
      yhalf = 0.5 * (y(i)+y(j))
      call ugline (option(1),x(i),y(i),0,segm)
      call ugline (option(1),xhalf,yhalf,1,segm)
      call ugline (option(2),x(j),y(j),1,segm)
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##     subroutine dash     ##
c     ##                         ##
c     #############################
c
c
      subroutine dash (i,j,option)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer maxseg
      parameter (maxseg=1000)
      integer i,j,k,segm
      real*8 r,dx,dy,dz,dnew
      real*8 dx1,dx2,dy1,dy2
      character*30 option(2)
      common /plot/ segm(maxseg)
c
c
      dx = x(j) - x(i)
      dy = y(j) - y(i)
      dz = z(j) - z(i)
      dnew = dz * 0.1/sqrt(dx*dx+dy*dy+dz*dz)
      dx1 = dx + dy*dnew
      dx2 = dx - dy*dnew
      dy1 = dy - dx*dnew
      dy2 = dy + dx*dnew
      do k = 1, 5
         r = float(k) / 10.0
         call ugline (option(1),x(i)+dx1*r,y(i)+dy1*r,0,segm)
         call ugline (option(1),x(i)+dx2*r,y(i)+dy2*r,1,segm)
      end do
      do k = 6, 10
         r = float(k) / 10.0
         call ugline (option(2),x(i)+dx1*r,y(i)+dy1*r,0,segm)
         call ugline (option(2),x(i)+dx2*r,y(i)+dy2*r,1,segm)
      end do
      return
      end
c
c
c     ##############################
c     ##                          ##
c     ##     subroutine wedge     ##
c     ##                          ##
c     ##############################
c
c
      subroutine wedge (i,j,option)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer maxseg
      parameter (maxseg=1000)
      integer i,j,segm
      real*8 dx,dy,dz,dnew
      real*8 x1,x2,y1,y2
      real*8 x1half,x2half,y1half,y2half
      character*30 option(2)
      common /plot/ segm(maxseg)
c
c
      dx = x(j) - x(i)
      dy = y(j) - y(i)
      dz = z(j) - z(i)
      dnew = dz * 0.1/sqrt(dx*dx+dy*dy+dz*dz)
      x1 = x(j) + dy*dnew
      x2 = x(j) - dy*dnew
      y1 = y(j) - dx*dnew
      y2 = y(j) + dx*dnew
      x1half = 0.5 * (x(i)+x1)
      x2half = 0.5 * (x(i)+x2)
      y1half = 0.5 * (y(i)+y1)
      y2half = 0.5 * (y(i)+y2)
      call ugline (option(1),x(i),y(i),0,segm)
      call ugline (option(1),x1half,y1half,1,segm)
      call ugline (option(2),x1,y1,1,segm)
      call ugline (option(2),x2,y2,1,segm)
      call ugline (option(2),x2half,y2half,1,segm)
      call ugline (option(1),x(i),y(i),1,segm)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine opendev  --  setup of UGSYS graphics device  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "opendev" opens a graphics device, clears the screen,
c     and defines the coordinates; coordinates are 11. x 8.5,
c     ie, like a sheet of standard sheet of paper
c
c
      subroutine opendev (device)
      implicit none
      integer maxseg
      parameter (maxseg=1000)
      integer segm,next
      character*5 device
      character*80 record
      common /plot/ segm(maxseg)
c
c
      write (*,10)
   10 format (/,' Graphics Device (NDS,HDS,VT340,HPP,IMA',
     &           ' or <CR>=EXIT) :  ',$)
      read (*,20)  record
   20 format (a80)
      next = 1
      call gettext (record,device,next)
      call upcase (device)
      if (device .eq. 'HPP') then
         call ugopen ('HP7550A',99)
      else if (device .eq. 'NDS') then
         call ugopen ('TEK4010',99)
      else if (device .eq. 'HDS') then
         call ugopen ('TEK4010',99)
      else if (device .eq. 'VT340') then
         call ugopen ('TEK4010',99)
      else if (device .eq. 'IMA') then
         call ugopen ('IMA3320',99)
      else
         stop
      end if
c
c     write any special terminal dependent control strings
c
      if (device .eq. 'VT340') then
         write (5,30)  char(27),char(91),char(63),    !! ESC[?38h
     &                 char(51),char(56),char(104)
   30    format (' ',6a1)
      end if
c
c     call UGSYS routines involved in graphics setup
c
      call ugpict ('CLEAR',0)              ! clear screen
      call uginit ('CLEAR',segm,1000)      ! clear grf
      call ugdspc ('PUT',11.0,8.5,1.0)     ! define drawing space
      call ugfont ('DUPLEX')
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine closdev  --  graceful close of UGSYS device  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine closdev (device)
      implicit none
      integer maxseg
      parameter (maxseg=1000)
      integer segm
      character*1 pause
      character*5 device
      common /plot/ segm(maxseg)
c
c
      call ugwrit (' ',0,segm)
      if (device .eq. 'NDS') then
         call ugclos ('NOCLEAR')
         read (5,10)  pause
   10    format (a1)
         write (5,20)  char(2)         !! STX
   20    format ('+',a1)
      else if (device.eq.'HDS' .or. device.eq.'VT340') then
         call ugclos ('NOCLEAR')
         read (5,30)  pause
   30    format (a1)
         write (5,40)  char(27),char(91),char(63),    !! ESC[?38l
     &                 char(51),char(56),char(108)
   40    format ('+',6a1)
      else
         call ugclos ('ALL')
      end if
      return
      end
c
c
c     ##############################
c     ##                          ##
c     ##     subroutine reset     ##
c     ##                          ##
c     ##############################
c
c
c     "reset" will close and then reopen the current graphics
c     device; at present this is used only for HP7550A plotter
c
c
      subroutine reset (device)
      implicit none
      integer maxseg
      parameter (maxseg=1000)
      integer segm
      character*5 device
      common /plot/ segm(maxseg)
c
c
c     first, close the current graphics device
c
      call closdev (device)
c
c     now, reopen the same graphics device
c
      if (device .eq. 'HPP') then
         call ugopen ('HP7550A',99)
      else if (device .eq. 'NDS') then
         call ugopen ('TEK4010',99)
      else if (device .eq. 'HDS') then
         call ugopen ('TEK4010',99)
      else if (device .eq. 'VT340') then
         call ugopen ('TEK4010',99)
      else if (device .eq. 'IMA') then
         call ugopen ('IMA3320',99)
      end if
c
c     call UGSYS routines involved in graphics setup
c
      call ugpict ('CLEAR',0)              ! clear screen
      call uginit ('CLEAR',segm,1000)      ! clear grf
      call ugdspc ('PUT',11.0,8.5,1.0)     ! define drawing space
      call ugfont ('DUPLEX')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ugxerr  --  error handler for UGSYS routines  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ugxerr" is the error processor subroutine for UGSYS
c
c     level    indicates the level number of the error received
c     sname    name of the subroutine which encountered the error
c     index    index number for the error; see UGSYS manual
c
c
c     subroutine ugxerr (level,sname,index)
c     implicit none
c     integer maxseg
c     parameter (maxseg=1000)
c     integer level,index,segm
c     character*8 sname
c     common /plot/ segm(maxseg)
c
c
c     if (index .eq. 11) then
c        call ugwrit (' ',0,segm)
c        call uginit ('CONTINUE',segm,500)
c        level = 0
c     end if
c     return
c     end

c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  program movie  --  display movie of molecular motion  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "movie" displays in rapid succession a series of molecular
c     coordinate sets on a Silicon Graphics workstation; requires
c     the TINKER package and the SGI GL graphics library
c
c     note that the parameter "maxatm" in subroutines "draw" and
c     "view" below must have the same value as in the "sizes.i"
c     dimensioning file for the TINKER package
c
c     command to link an executable:
c
c     f77 -o movie movie.f -lfgl -lgl libtinker.a
c
c
      program movie
      implicit none
      include 'fget.h'
c
c
c     initialize the TINKER related stuff
c
      call initial
c
c     read in the molecule frame coordinates
c
      call coord
c
c     initialize the screen for graphics
c
      call graphic
c
c     draw the molecular movie
c
      call draw
c
c     return the monitor to mono mode
c
      call setmon (HZ60)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ########################
c     ##                    ##
c     ##  subroutine coord  ##
c     ##                    ##
c     ########################
c
c
      subroutine coord
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'output.i'
      integer maxfrm
      parameter (maxfrm=200)
      integer i,j,k,iarc,label
      integer leng,trimtext,freeunit
      integer nframe,next,nline
      integer hnline,hline,line
      integer thousand,hundred,tens,ones
      integer bbcount
      real*4 xyz,xnew,ynew,znew
      logical exist
      character*1 digit(0:9)
      character*4 text,string
      character*60 arcfile
      character*80 record,substring
      common /bkbone/ bbcount
      common /coords/ nframe,xyz(3,maxatm,maxfrm)
      common /lines / nline,hnline,hline(2,maxbnd),line(2,maxbnd)
      common /text  / text(maxfrm)
      data digit  / '0','1','2','3','4','5','6','7','8','9' /
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (arcfile,exist)
      if (exist) then
         call basefile (arcfile)
         call suffix (arcfile,'arc')
         call version (arcfile,'old')
         inquire (file=arcfile,exist=exist)
      end if
c
c     ask for the user specified input archive filename
c
      dowhile (.not. exist)
         write (iout,10)
   10    format (/,' Enter the Name of the Archive File :  ',$)
         read (input,20)  arcfile
   20    format (a60)
         call basefile (arcfile)
         call suffix (arcfile,'arc')
         call version (arcfile,'old')
      end do
c
c     get the first molecule file
c
      leng = trimtext (arcfile)
      iarc = freeunit ()
      open (unit=iarc,file=arcfile(1:leng))
      call readxyz (iarc)
      rewind (unit=iarc)
c
c     generate array of connected atoms
c
      hnline = 0
      do i = 1, n
         do j = 1, n12(i)
            k = i12(j,i)
            if (i .lt. k ) then
               hnline = hnline + 1
               hline(1,hnline) = i
               hline(2,hnline) = k
            end if
         end do
      end do
c
c     generate array of connected atoms, without the hydrogens
c
      nline = 0
      do i = 1, n
         if (name(i)(1:1) .ne. 'H') then
            do j = 1, n12(i)
               k = i12(j,i)
               if (i.lt.k .and. name(k)(1:1) .ne. 'H') then
                  nline = nline + 1
                  line(1,nline) = i
                  line(2,nline) = k
               end if
            end do
         end if
      end do
c
c     determine the number of atoms in the amino acid backbone
c     (assumes an old TINKER format, which is no longer used)
c
      bbcount = 0
      do i = 1, n/3
         if (name((3*i)-2) .eq. 'N  ' .or.
     &       name((3*i)-2) .eq. 'N+ ' .and.
     &       name((3*i)-1) .eq. 'C  ' .and.
     &       name(3*i)     .eq. 'C  ') then
            bbcount = bbcount + 3
         else
            goto 30
         end if
      end do
   30 continue
c
c     set up the array, xyz, to store x y & z coordinates;
c     two loops over the atoms are included, the first is
c     completely general, the second is faster but depends
c     on use of the record format written by "prtxyz"
c
      nframe = 0
      dowhile (.true.)
         read (iarc,40,err=90,end=90)  record
   40    format (a80)
         nframe = nframe + 1
         do k = 1, n
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           read (iarc,50)  record                   c
c  50       format (a80)                             c
c           read (record,*)  label                   c
c           call getword (record,name,next)          c
c           substring = record(next:80)              c
c           read (substring,*)  xnew,ynew,znew       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
            read (iarc,60)  xnew,ynew,znew
   60       format (11x,3f12.6)
            xyz(1,k,nframe) = xnew
            xyz(2,k,nframe) = ynew
            xyz(3,k,nframe) = znew
         end do
         thousand = nframe / 1000
         hundred = (nframe - 1000*thousand) / 100
         tens = (nframe - 1000*thousand - 100*hundred) / 10
         ones = nframe - 1000*thousand - 100*hundred - 10*tens
         string(1:1) = digit(thousand)
         string(2:2) = digit(hundred)
         string(3:3) = digit(tens)
         string(4:4) = digit(ones)
         do k = 1, 4
            if (string(k:k) .ne. '0') goto 70
            string(k:k) = ' '
         end do
   70    continue
         text(nframe) = string
         if (mod(nframe,10) .eq. 0) then
            write (iout,80)  nframe
   80       format (' Processed',i5,' Coordinate Sets')
         end if
      end do
   90 continue
      if (mod(nframe,10) .ne. 0) then
         write (iout,100)  nframe
  100    format (' Processed',i5,' Coordinate Sets')
      end if
      close (unit=iarc)
      return
      end
c
c
c     ##########################
c     ##                      ##
c     ##  subroutine graphic  ##
c     ##                      ##
c     ##########################
c
c
      subroutine graphic
      implicit none
      include 'gl/fgl.h'
      include 'fget.h'
      integer xmaxscreen,ymaxscreen
      parameter (xmaxscreen=1279)
      parameter (ymaxscreen=1023)
      integer index
c
c
      call prefpo (0,xmaxscreen,0,ymaxscreen)
      index = winope ("MOVIE",5)
      call double
      call RGBmod
      call gconfi
      call viewpo (0,xmaxscreen,0,ymaxscreen)
      call frontb (.true.)
      call RGBcol (0,0,0)
      call clear
      call frontb (.false.)
      call swapin (1)
      call setmon (STRREC)
      call font (0)
      call depthc (.TRUE.)
      return
      end
c
c
c     #######################
c     ##                   ##
c     ##  subroutine draw  ##
c     ##                   ##
c     #######################
c
c
      subroutine draw
      implicit none
      integer maxatm,maxfrm
      parameter (maxatm=3000)
      parameter (maxfrm=200)
      include 'gl/fgl.h'
      include 'fdevice.h'
      include 'atoms.i'
      integer j,k,frame,nframe,value1
      integer dialpos0,dialpos1,dialpos2,dialpos3
      integer dialpos4,dialpos5,dialpos6,dialpos7
      integer color,colorvalue
      integer mainmenu,colormenu,sidemenu,bkbonemenu
      real*4 xyz,dist,mid,xrot,yrot,zrot
      real*4 xave,yave,zave,xmid,ymid,zmid
      real*4 xtiny,ytiny,ztiny,xhuge,yhuge,zhuge
      common /colors/ color,colorvalue(3,7)
      common /coords/ nframe,xyz(3,maxatm,maxfrm)
      common /dial  / dialpos0,dialpos1,dialpos2,dialpos3,
     &                dialpos4,dialpos5,dialpos6,dialpos7
      common /menus / mainmenu,colormenu,sidemenu,bkbonemenu
      common /rotate/ xrot,yrot,zrot
      common /scale / dist,mid
      common /trans / xave,yave,zave,xmid,ymid,zmid
c
c
c     set initial values for the dial positions
c
      call setval (DIAL0,320,0,1279)
      call setval (DIAL1,1,1,30)
      call setval (DIAL2,640,0,1279)
      call setval (DIAL3,640,0,1279)
      call setval (DIAL4,640,0,1279)
      call setval (DIAL5,640,0,1279)
      call setval (DIAL6,640,0,1279)
      call setval (DIAL7,640,0,1279)
      dialpos0 = 320
      dialpos1 = 1
      dialpos2 = 640
      dialpos3 = 640
      dialpos4 = 640
      dialpos5 = 640
      dialpos6 = 640
      dialpos7 = 640
c
c     determine dimensions for the molecule box
c
      xtiny = xyz(1,1,1)
      xhuge = xtiny
      ytiny = xyz(2,1,1)
      yhuge = ytiny
      ztiny = xyz(3,1,1)
      zhuge = ztiny
      do k = 1, nframe
         do j = 1, n
            if (xyz(1,j,k) .lt. xtiny)  xtiny = xyz(1,j,k)
            if (xyz(1,j,k) .gt. xhuge)  xhuge = xyz(1,j,k)
            if (xyz(2,j,k) .lt. ytiny)  ytiny = xyz(2,j,k)
            if (xyz(2,j,k) .gt. yhuge)  yhuge = xyz(2,j,k)
            if (xyz(3,j,k) .lt. ztiny)  ztiny = xyz(3,j,k)
            if (xyz(3,j,k) .gt. zhuge)  zhuge = xyz(3,j,k)
         end do
      end do
      mid = 0.5 * max(xhuge-xtiny,yhuge-ytiny,zhuge-ztiny)
      xmid = 0.5 * (xtiny+xhuge)
      ymid = 0.5 * (ytiny+yhuge)
      zmid = 0.5 * (ztiny+zhuge)
c
c     set initial scale factor, translation and rotation
c
      dist = mid
      xave = xmid
      yave = ymid
      zave = zmid
      xrot = 0.0
      yrot = 0.0
      zrot = 0.0
c
c     define the color submenu
c
      colormenu = newpup()
      call addtop (colormenu,'RED %x1',7,0)
      call addtop (colormenu,'GREEN %x2',9,0)
      call addtop (colormenu,'YELLOW %x3',10,0)
      call addtop (colormenu,'BLUE %x4',8,0)
      call addtop (colormenu,'MAGENTA %x5',11,0)
      call addtop (colormenu,'CYAN %x6',8,0)
      call addtop (colormenu,'WHITE %x7',9,0)
c
c     define the side(chain) submenu
c
      sidemenu = newpup()
      call addtop (sidemenu,'RESIDUE %x14',12,0)
c
c     define backbone submenu
c
      bkbonemenu = newpup()
      call addtop (bkbonemenu,'ONLY %x10',8,0)
      call addtop (bkbonemenu,'COLOR %x11',10,0)
c
c     define the main menu
c
      mainmenu = newpup()
      call addtop (mainmenu,'COLORS %m',9,colormenu)
      call addtop (mainmenu,'DEPTH %x8',9,0)
      call addtop (mainmenu,'SIDE %m',7,sidemenu)
      call addtop (mainmenu,'BKBONE %x10',11,bkbonemenu)
      call addtop (mainmenu,'H-ATOMS %x12',12,0)
      call addtop (mainmenu,'QUIT %x13',9,0)
c
c     initialize color scales
c
      colorvalue(1,1) = 255
      colorvalue(2,1) = 0
      colorvalue(3,1) = 0
      colorvalue(1,2) = 0
      colorvalue(2,2) = 255
      colorvalue(3,2) = 0
      colorvalue(1,3) = 255
      colorvalue(2,3) = 255
      colorvalue(3,3) = 0
      colorvalue(1,4) = 0
      colorvalue(2,4) = 0
      colorvalue(3,4) = 255
      colorvalue(1,5) = 255
      colorvalue(2,5) = 0
      colorvalue(3,5) = 255
      colorvalue(1,6) = 0
      colorvalue(2,6) = 255
      colorvalue(3,6) = 255
      colorvalue(1,7) = 255
      colorvalue(2,7) = 255
      colorvalue(3,7) = 255
c
c     draw the molecule, looping continuously over frames
c
      dowhile (.true.)
         do frame = 1, nframe
            value1 = getval(DIAL1)
            do k = 1, value1
               call dials
               call view ('RIGHT',frame)
               call view ('LEFT',frame)
            end do
         end do
         do frame = nframe-1, 2, -1
            value1 = getval(DIAL1)
            do k = 1, value1
               call dials
               call view ('RIGHT',frame)
               call view ('LEFT',frame)
            end do
         end do
      end do
      return
      end
c
c
c     ########################
c     ##                    ##
c     ##  subroutine dials  ##
c     ##                    ##
c     ########################
c
c
      subroutine dials
      implicit none
      include 'gl/fgl.h'
      include 'fdevice.h'
      integer orig(0:7)
      integer dialpos0,dialpos1,dialpos2,dialpos3
      integer dialpos4,dialpos5,dialpos6,dialpos7
      real*4 xrot,yrot,zrot,dist,mid
      real*4 xave,yave,zave,xmid,ymid,zmid
      common /dial  / dialpos0,dialpos1,dialpos2,dialpos3,
     &                dialpos4,dialpos5,dialpos6,dialpos7
      common /rotate/ xrot,yrot,zrot
      common /scale / dist,mid
      common /trans / xave,yave,zave,xmid,ymid,zmid
      save orig
      data orig / 320, 1, 640, 640, 640, 640, 640, 640 /
c
c
c     check for any change in the dial positions
c
      if (getval(DIAL0) .ne. dialpos0) then
         dist = mid * orig(0)/getval(DIAL0)
         dialpos0 = getval(DIAL0)
      end if
      if (getval(DIAL2) .ne. dialpos2) then
         zave = zmid + (getval(DIAL2) - orig(2))/10.0
         dialpos2 = getval(DIAL2)
      end if
      if (getval(DIAL4) .ne. dialpos4) then
         yave = ymid + (getval(DIAL4) - orig(4))/10.0
         dialpos4 = getval(DIAL4)
      end if
      if (getval(DIAL6) .ne. dialpos6) then
         xave = xmid + (getval(DIAL6) - orig(6))/10.0
         dialpos6 = getval(DIAL6)
      end if
      if (getval(DIAL3) .ne. dialpos3) then
         zrot = orig(3) - getval(DIAL3)
         dialpos3 = getval(DIAL3)
      end if
      if (getval(DIAL5) .ne. dialpos5) then
         yrot = orig(5) - getval(DIAL5)
         dialpos5 = getval(DIAL5)
      end if
      if (getval(DIAL7) .ne. dialpos7) then
         xrot = orig(7) - getval(DIAL7)
         dialpos7 = getval(DIAL7)
      end if
      return
      end
c
c
c     #######################
c     ##                   ##
c     ##  subroutine view  ##
c     ##                   ##
c     #######################
c
c
      subroutine view (eye,frame)
      implicit none
      include 'gl/fgl.h'
      include 'fdevice.h'
      integer maxatm,maxbnd,maxfrm
      parameter (maxatm=3000)
      parameter (maxbnd=6*maxatm/5)
      parameter (maxfrm=200)
      integer j,k,frame,nframe,nline,hnline,hline,line
      integer mainmenu,colormenu,sidemenu,bkbonemenu
      integer color,tempcolor,colorvalue,chaincolor
      integer bbcount
      real*4 xyz,dist,mid,xrot,yrot,zrot
      real*4 xave,yave,zave,xmid,ymid,zmid
      real*4 bottom,top,left,right,near,far
      character*4 text
      character*5 eye
      common /bkbone/ bbcount
      common /colors/ color,colorvalue(3,7)
      common /coords/ nframe,xyz(3,maxatm,maxfrm)
      common /lines / nline,hnline,hline(2,maxbnd),line(2,maxbnd)
      common /menus / mainmenu,colormenu,sidemenu,bkbonemenu
      common /rotate/ xrot,yrot,zrot
      common /scale / dist,mid
      common /text  / text(maxfrm)
      common /trans / xave,yave,zave,xmid,ymid,zmid
c
c     perspective variables; should not be changed
c
      real*4 fv,fl(2),fr(2),ft(2)
      save fv,fl,fr,ft
      data fv / 1.071797 /
      data fl / 1.185464, 1.052131 /
      data fr / 1.052131, 1.185464 /
      data ft / -0.1, 0.1 /
c
c     viewport sizes; should not be changed
c
      integer viewleft,viewright
      integer viewbot(2),viewtop(2),viewnum(2)
      save viewleft,viewright
      save viewbot,viewtop,viewnum
      data viewleft  / 128 /
      data viewright / 1151 /
      data viewbot   / 0, 532 /
      data viewtop   / 491, 1023 /
      data viewnum   / 466, 998 /
c
c     default color for molecule is YELLOW
c
      data color / 3 /
c
c     monitor setting - mono
c
      integer*4 HZ60
      data HZ60 / 1 /
c
c     depth cue default is OFF
c
      logical DEPTHON
      data DEPTHON /.FALSE./
      save DEPTHON
c
c     hydrogen atom display default is OFF
c
      logical HATOMON
      data HATOMON /.TRUE./
      save HATOMON
c
c     backbone atom display default is OFF
c
      logical BBONE
      data BBONE /.FALSE./
      save BBONE
c
c     set parameters for right or left eye
c
      if (eye .eq. 'RIGHT') then
         k = 1
      else
         k = 2
      end if
c
c     set the viewport and clear the screen
c
      call viewpo (viewleft,viewright,viewbot(k),viewtop(k))
      call RGBcol (0,0,0)
      call clear
c
c     set the window for the molecule
c
      left = xave - fl(k)*dist
      right = xave + fr(k)*dist
      bottom = yave - fv*dist
      top = yave + fv*dist
      near = 4.0 * dist
      far = 8.0 * dist
      call window (left,right,bottom,top,near,far)
      call lsetde (int(near),int(far)+1)
c
c     rotate and translate the molecule
c
      call transl (xmid+ft(k)*dist,ymid,zmid+zave-6.0*dist)
      call rot (xrot,'x')
      call rot (yrot,'y')
      call rot (zrot,'z')
      call transl (-xmid,-ymid,-zmid)
c
c     get user input of molecule color
c
      if ( getbut(LEFTMO) ) then
         tempcolor = dopup (mainmenu)
         if (tempcolor.ge.1 .and. tempcolor.le.7) then
            color = tempcolor
         else if (tempcolor .eq. 8) then
            if (DEPTHON) then
               DEPTHON = .FALSE.
            else
               DEPTHON = .TRUE.
            end if
         else if (tempcolor .eq. 12) then
            if (HATOMON) then
               HATOMON = .FALSE.
            else
               HATOMON = .TRUE.
            end if
         else if (tempcolor .eq. 9) then
            chaincolor = dopup (colormenu)
         else if (tempcolor .eq. 10) then
            if (BBONE) then
               BBONE = .FALSE.
            else if (bbcount .ne. 0) then
               BBONE = .TRUE.
            end if
c         else if (tempcolor .eq. 11) then
c            bbcolor = dopup(colormenu)
         else if (tempcolor .eq. 14) then
            call sidecolor
         else if ( tempcolor .eq. 13 ) then
            call setmon (HZ60)
            STOP
         end if
      end if
c
c     draw the molecule as a collection of bond lines
c
      if (DEPTHON) then
           call lRGBra (0,0,0,colorvalue(1,color),
     &                  colorvalue(2,color),colorvalue(3,color),
     &                  int(near),int(far)+1)
      else
           call lRGBra (colorvalue(1,color),colorvalue(2,color),
     &                  colorvalue(3,color),colorvalue(1,color),
     &                  colorvalue(2,color),colorvalue(3,color),
     &                  int(near),int(far)+1)
      end if
      if (BBONE) then
        do j = 1, bbcount
           call bgnlin
           call v3f (xyz(1,j,frame))
           call v3f (xyz(1,j+1,frame))
           call endlin
        end do
      else
        if (hatomon) then
           do j = 1, hnline
              call bgnlin
              call v3f (xyz(1,hline(1,j),frame))
              call v3f (xyz(1,hline(2,j),frame))
              call endlin
           end do
        else
           do j = 1, nline
              call bgnlin
              call v3f (xyz(1,line(1,j),frame))
              call v3f (xyz(1,line(2,j),frame))
              call endlin
           end do
        end if
      end if
c
c     draw the frame number in upper left corner
c
      call viewpo ( 10,110,viewnum(k),viewtop(k) )
      call RGBcol (0,0,0)
      call clear
      call RGBcol (  0,255,255 )
      call ortho2 (-1.0,1.0,-1.0,1.0)
      call cmov2 (-0.3,-0.2)
      call charst (text(frame),4)
c
c     swap the front and back buffers
c
      call swapbu
      return
      end
c
c
c     ############################
c     ##                        ##
c     ##  subroutine sidecolor  ##
c     ##                        ##
c     ############################
c
c
c     "sidecolor" retrieves the number of an amino acid
c     sidechain to be colored from the user
c
c
      subroutine sidecolor
      implicit none
      include 'gl/fgl.h'
      include 'fget.h'
      include 'iounit.i'
      integer*4 resnumber,color
      integer mainmenu,colormenu,sidemenu,bkbonemenu
      common /menus / mainmenu,colormenu,sidemenu,bkbonemenu
c
c
      call setmon (HZ60)
      write (iout,10)
   10 format (/,'What residue would you like to color?  ',$)
      read (input,20) resnumber
   20 format (i4)
c     color = dopup(colormenu)
      call setmon (STRREC)
      return
      end

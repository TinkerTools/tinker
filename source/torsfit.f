c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2010 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program torsfit  --  fit torsional force field parameters  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "torsfit" refines torsional force field parameters based on
c     a quantum mechanical potential surface and analytical gradient
c
c
      program torsfit 
      implicit none
      include 'sizes.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      integer i,nvar,next
      integer length
      integer torbnd(10)
      logical exist,query
      character*20 keyword
      character*120 record
      character*120 string
      character*120 xyzfile
c
c
c     initialization of the various modes of operation
c
      call initial
c
c     read the Cartesian coordinates and connectivity info
c
      call getxyz
      xyzfile = filename
      length = leng
c
c     read structure and vibrational data from Gaussian output
c
      call getkey
      call mechanic
c
c     choose the first torsion based on the center bond atoms
c
      do i = 1, 10
         torbnd(i) = 0
      end do
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  torbnd(1),torbnd(2)
         query = .false.
      end if
   10 continue
      if (query) then
         dowhile (torbnd(1).eq.0 .or. torbnd(2).eq.0) 
            write (iout,20)
   20       format (/,' Two Center Bond Atoms of the Torsion : ',$)
            read (input,*,err=30,end=30)  torbnd(1),torbnd(2) 
         end do
   30    continue
      end if
c
c     choose the second torsion based on the center bond atoms
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then 
         read (string,*,err=40,end=40)  torbnd(3),torbnd(4)
         query = .false.
      end if
   40 continue
      if (query) then 
         write(iout,50) 
   50    format (/,' Two Center Bond Atoms for the 2nd Torsion ',
     &              '[Optional, <CR>=Exit] : ',$)
         read (input,60,err=70,end=70)  record 
   60    format (a120)
         read (record,*,err=70,end=70)  torbnd(3),torbnd(4)
   70    continue
      end if
c
c     fit the torsional parameters based on potential surface
c
      call fittorsion (torbnd)
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine fittorsion  --  fit torsion parameters           ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "fittorsion" refines torsion parameters based
c     on a quantum mechanical optimized energy surface
c
c
      subroutine fittorsion(torbnd)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmtyp.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'ktorsn.i'
      include 'kgeoms.i'
      include 'math.i'
      include 'output.i'
      include 'potent.i'
      include 'qmstuf.i'
      include 'scales.i'
      include 'tors.i'
      include 'usage.i'
      integer maxfittor,maxconf
      parameter (maxfittor = 12)
      parameter (maxconf = 500)
      integer torbnd(10)
      integer torcrs(4,maxfittor)
      integer ctorid(maxfittor)
      integer ftorid(maxfittor)
      integer i,j,k,ii,jj,kk
      integer ia,ib,ic,id
      integer ita,itb,itc,itd
      integer itmpa,itmpb,otfix
      integer ntorfit,ntorcrs,nconf
      integer size
      real*8  tmpa,tmpb
      real*8  eqm(maxconf),emm(maxconf)
      real*8  eqmmin,emmmin
      real*8  erqm(maxconf),ermm(maxconf)
      real*8  delte(maxconf)
      real*8  ftv(maxconf,maxfittor),tv
      real*8  ctv(maxconf,9*maxfittor)
      real*8  rftv(maxconf,maxfittor)
      real*8  fwt(maxconf)
      real*8  vxx(6*maxfittor),vcon,vxxl(6*maxfittor)
      real*8  zrms,rms,torf(maxconf)
      real*8  tord(6*maxfittor,6*maxfittor)
      real*8  coeff(maxconf,6*maxfittor)
      real*8  minimum,grdmin
      integer istep,maxstep,nvxx,ivxx
      integer tflg(maxfittor),cflg(9*maxfittor)
      logical vflg (6,maxfittor),done
      logical confvisited(maxconf)
      character*8 pkt(maxfittor)
      character*4 pa,pb,pc,pd
      character*16 kft(maxfittor)
      character*16 kct(9*maxfittor)
      character*120 record
      character*120 oldfilename
      integer oldleng
      character*120 oldkeyline(maxkey)
      integer oldnkey
      real*8 mata(6*maxfittor,6*maxfittor),vectb(6*maxfittor)
      real*8 avedl
      integer nvar
      integer refconf(maxconf)
      real*8 xx(maxvar)
      external minimiz2
      external optsave
      real*8 energy, geometry
      integer ikey,freeunit
      integer trimtext
      character*120 keyfile
c
c
c     set initial values
c
      ntorfit = 0
      ntorcrs = 0   
      otfix = ntfix
      istep = 0   
      tv = 0.0d0
      vcon = 0.5d0
      do i = 1, maxfittor
         ftorid (i) = 0
         tflg (i) = 0
         do j = 1, 6
            vflg(j,i) = .false.
         end do
      end do
      do i = 1, 6*maxfittor
         vxx(i) = 0.0d0
         vxxl(i)= 0.1d0
         avedl = 0.0d0
         do j = 1, 6*maxfittor
           tord (i,j) = 0.0d0
         end do
      end do
      do i = 1, 9*maxfittor
         cflg(i) = 0
      end do
      do i = 1, maxconf
        fwt (i) = 1.0d0
        torf(i) = 0.0d0
        confvisited(i) = .false.
      end do
      do i = 1, maxconf
         do j = 1, 6*maxfittor
           coeff(i,j) = 0.0d0
         end do
         refconf = 0
      end do    
      grdmin = 0.01   
      if (torbnd(1) .gt. torbnd(2) ) then
         itmpa = torbnd(1)
         torbnd(1) = torbnd(2)
         torbnd(2) = itmpa
      end if
c
c     store the file system information
c
      oldfilename = filename
      oldleng = leng
      oldnkey = nkey
      do i = 1, nkey
         oldkeyline(i) = keyline(i)
      end do
c
c     check all the torsions cross the two center bond atoms
c                        
      write (*,10)
   10 format (/,' Torsions Crossing the Center Bond :',/) 
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         itmpa = ib
         itmpb = ic
         if (itmpa .gt. itmpb) then
            j = itmpa
            itmpa = itmpb
            itmpb = j
         end if
         if ((torbnd(1) .eq. itmpa .and. torbnd(2) .eq. itmpb)
     &  .or.(torbnd(3) .eq. itmpa .and. torbnd(4) .eq. itmpb)) 
     &   then
            ntorcrs = ntorcrs + 1
            torcrs(1,ntorcrs) = ia
            torcrs(2,ntorcrs) = ib
            torcrs(3,ntorcrs) = ic
            torcrs(4,ntorcrs) = id
            ctorid(ntorcrs) = i
            write (*,20)  ntorcrs,ia,ib,ic,id
   20       format (1x,'Torsion No. ',i5,' : ',3(i5,' ---'),i5)
         end if
       end do   
c
c      choose the specific torsion for fitting
c        
       write (*,30)
   30  format (/,' Choose the Torsions for Fitting from Above List ') 
       read (*,40,err=50,end=50)  record 
   40  format (a120)
   50  continue
       read (record,*,err=60,end=60)  (ftorid(i),i=1,ntorcrs)      
   60  continue
c
c      count the torsions to be fitted
c                                     
       do i = 1, ntorcrs
          if (ftorid(i) .gt. 0)  ntorfit = ntorfit + 1
       end do
c
c      get the number of conformations for fitting
c       
       write (*,70)
   70  format (/,' Enter Total Conformation Number: ',$)
       read (*,*)  nconf
c
c      read the QM coordinates and conformations energies
c     
       do i = 1, nconf
          call readgau
          write (*,80)  i 
   80     format (/ ,' Finished Reading Conformation',i4)
          do j = 1, n
             x(j) = gx(j)
             y(j) = gy(j)
             z(j) = gz(j)
          end do
          call makeref (i) 
          eqm(i) = egau
       end do
c
c      calculate the relative QM conformational energies
c 
       eqmmin = eqm(1)
       do i = 2, nconf
          if (eqm(i) .lt. eqmmin) eqmmin = eqm(i) 
       end do

       do i = 1, nconf
          erqm(i) = eqm(i) - eqmmin
          write (*,83)  i,erqm(i)
   83     format(' Relative Conformational Energy (QM) ',
     &              i5,1x,f8.5,' kcal/mole')
       end do
c
c      get fitting torsion type (atom classes)
c
       do i = 1, ntorfit
          j = ftorid(i)
          k = ctorid(j)
          ia = itors(1,k)
          ib = itors(2,k)
          ic = itors(3,k)
          id = itors(4,k)
          ita = class(ia)
          itb = class(ib)
          itc = class(ic)
          itd = class(id)
          size = 4
          call numeral (ita,pa,size)
          call numeral (itb,pb,size)
          call numeral (itc,pc,size)
          call numeral (itd,pd,size)
          if (itb .le. itc) then
             kft(i) = pa//pb//pc//pd
          else 
             kft(i) = pd//pc//pb//pa
          end if
       end do
c
c      get all the cross torsion types
c
       do i = 1, ntorcrs
          k = ctorid(i)
          ia = itors(1,k)
          ib = itors(2,k)
          ic = itors(3,k)
          id = itors(4,k)
          ita = class(ia)
          itb = class(ib)
          itc = class(ic)
          itd = class(id)
          size = 4
          call numeral (ita,pa,size)
          call numeral (itb,pb,size)
          call numeral (itc,pc,size)
          call numeral (itd,pd,size)
          if (itb .le. itc) then
             kct(i) = pa//pb//pc//pd
          else
             kct(i) = pd//pc//pb//pa
          end if
       end do
c
c      initialize the torsion and geometry restrain parameters 
c
       nvxx = 0
       do i = 1, ntorfit
          j = ftorid(i)
          k = ctorid(j)
          done = .false. 
          tflg(i) = 0
          ia = itors(1,k)
          ib = itors(2,k)
          ic = itors(3,k)
          id = itors(4,k)
          ita = class(ia)
          itb = class(ib)
          itc = class(ic)
          itd = class(id)
          write(*,13) ita,itb,itc,itd,
     &    tors1(1,k),tors2(1,k),tors3(1,k),
     &    tors4(1,k),tors5(1,k),tors6(1,k)
   13     format(1x,'torsion ',4i4,6f8.3)
          do ii = 1, i-1
             jj = ftorid(ii)
             kk = ctorid(jj)
             if (kft(i) .eq. kft(ii)) then
                done = .true.
                tflg(i) = ii
                goto 85
             end if                
          end do
          do ii = 1, ntorcrs
             if (kct(ii) .eq. kft(i)
     &          .and. ii .ne. j) cflg(ii) = j
          end do
          if (abs(tors1(1,k)) .gt. 0.0d0) then
             nvxx = nvxx +1
             vxx(nvxx) = tors1(1,k)
             vflg(1,i) = .true.
          end if
          tors1(1,k) = 0.0d0
          tors1(2,k) = 0.0d0
          if (abs(tors2(1,k)) .gt. 0.0d0) then
             nvxx = nvxx +1
             vxx(nvxx) = tors2(1,k)
             vflg(2,i) = .true.
          end if
          tors2(1,k) = 0.0d0
          tors2(2,k) = 180.0d0
          if (abs(tors3(1,k)) .gt. 0.0d0) then
             nvxx = nvxx +1
             vxx(nvxx) = tors3(1,k)
             vflg(3,i) = .true.
          end if
          tors3(1,k) = 0.0d0
          tors3(2,k) = 0.0d0
          if (abs(tors4(1,k)) .gt. 0.0d0) then
             nvxx = nvxx +1
             vxx(nvxx) = tors4(1,k)
             vflg(4,i) = .true.
          end if
          tors4(1,k) = 0.0d0
          tors4(2,k) = 180.0d0
          if (abs(tors5(1,k)) .gt. 0.0d0) then
             nvxx = nvxx +1
             vxx(nvxx) = tors5(1,k)
             vflg(5,i) = .true.
          end if
          tors5(1,k) = 0.0d0
          tors5(2,k) = 0.0d0
          if (abs(tors6(1,k)) .gt. 0.0d0) then
             nvxx = nvxx +1
             vxx(nvxx) = tors6(1,k)
             vflg(6,i) = .true.
          end if
          tors6(1,k) = 0.0d0
          tors6(2,k) = 180.0d0
          ntfix = ntfix+1
          itfix(1,ntfix) = ia
          itfix(2,ntfix) = ib
          itfix(3,ntfix) = ic
          itfix(4,ntfix) = id
          tfix(1,ntfix) = 5.0d0
          write(*,*) 'Fixed Torsion ',ia,ib,ic,id
   85     continue 
      end do
c
c     print torsion flags (check duplicated torsion types)
c
      do i = 1, ntorfit
         write (*,14)  i,tflg(i)
  14     format (' Fitting Torsion No. ',i4,' Flag ',i4)
         do j = 1, 6
            write (*,24)  i,j,vflg(j,i)
  24        format (' Variable ',2i4,' Variable Flag ',l4)
         end do
      end do
c
c     print torsion flags for all the torsions cross the bond
c   
      write (*,*)  'All the Torsions Cross the Bond'
      do i = 1, ntorcrs
         k = ctorid(i) 
         if (cflg(i) .gt. 0) then
            tors1(1,k) = 0.0d0
            tors2(1,k) = 0.0d0
            tors3(1,k) = 0.0d0
            tors4(1,k) = 0.0d0
            tors5(1,k) = 0.0d0
            tors6(1,k) = 0.0d0
         end if
         write (*,14)  i,cflg(i)
      end do
c
c     add one constant variable
c
      nvxx = nvxx + 1
      vxx(nvxx) = vcon
c
c     get inital energy difference
c       
      do i = 1, nconf
         call getref(i)
         kk = 0
         do j = 1, ntorfit
            k = ftorid(j)
            ia = torcrs(1,k)
            ib = torcrs(2,k)
            ic = torcrs(3,k)
            id = torcrs(4,k)
            ftv(i,j) = geometry(ia,ib,ic,id)
            write (*,34)  i,j,ftv(i,j)
   34       format (' Fitting Torsion Value',2i5,f12.4)
            if (tflg(j) .eq. 0) then
               kk = kk+1
               tfix(2,otfix+kk) = ftv(i,j)
               tfix(3,otfix+kk) = tfix(2,otfix+kk)
            end if
         end do
         do k = 1, ntorcrs
            ia = torcrs(1,k)
            ib = torcrs(2,k)
            ic = torcrs(3,k)
            id = torcrs(4,k)
            ctv(i,k) = geometry(ia,ib,ic,id)
         end do 
c
c      scale the coordinates of each active atom
c
         nvar = 0
         do j = 1, n
            if (use(j)) then
               nvar = nvar + 1
               scale(nvar) = 12.0d0
               xx(nvar) = x(j) * scale(nvar)
               nvar = nvar + 1
               scale(nvar) = 12.0d0
               xx(nvar) = y(j) * scale(nvar)
               nvar = nvar + 1
               scale(nvar) = 12.0d0
               xx(nvar) = z(j) * scale(nvar)
            end if
         end do
c
c      make the call to the optimization routine
c
         write (*,90) i
   90    format (/,' Minimizing Structure No. ',i3)
         coordtype = 'CARTESIAN'
         use_geom = .true.
         grdmin = 0.01d0
         iwrite = 0
         iprint = 0
         call lbfgs (nvar,xx,minimum,grdmin,minimiz2,optsave)
c
c     unscale the final coordinates for active atoms
c
         nvar = 0
         do j = 1, n
            if (use(j)) then
               nvar = nvar + 1
               x(j) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               y(j) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               z(j) = xx(nvar) / scale(nvar)
            end if
         end do
         emm(i) = energy()
      end do
c
c     calculate relative value for each torsional angle
c
      do i = 1, nconf
         do j = 1, ntorfit
            rftv(i,j) = ftv(i,j)-ftv(i,1)
         end do
      end do
c
c     calculate the relative MM energies
c   
      emmmin = emm(1)
      do i = 2, nconf
         if (emm(i) .lt. emmmin) emmmin = emm(i)
      end do 
c
c     calcualte the energy difference and RMS
c
      rms = 0.0d0
      zrms = 0.0d0
      do i = 1, nconf
         ermm (i) = emm(i) - emmmin 
         delte (i) = erqm (i) - ermm(i) 
         rms = rms + delte(i)*delte(i)
         write(*,93) i,ermm(i)
   93    format(1x,'Relative Conformational Energy (MM) ',
     &          i5,1x,f12.5,' kcal/mole')
      end do               
      rms = sqrt(rms/dble(nconf))
      zrms = rms
      write(*,*) "Energy Diff RMS ",rms
c
c     calculate the weights 
c
c     do i = 1,nconf
c        do j = 1,nconf
c           if (.not. confvisited(j)) then
c              tmpa = ftv(j,1)
c              itmpa = j
c              confvisited(j) = .true.
c              goto 16
c           end if
c        end do
c   16   continue
c        do j = 1,nconf
c           if ((ftv(j,1) .lt. tmpa) .and. (.not. confvisited(j))) then
c              confvisited(itmpa) = .false.
c              itmpa = j
c              tmpa = ftv(j,1)
c           end if
c        end do
c        refconf(itmpa) = i
c        confvisited(itmpa) = .true.
c        write(*,*) itmpa,' <===> ',refconf(itmpa)
c     end do
          
       if (nconf .gt. 1 .and. torbnd(3) .eq. 0) then
          if((ftv(nconf,1)+180.0d0) .lt. 1.0d0 
     &       .and. ftv(nconf-1,1) .gt. 0.0d0)
     &       ftv(nconf,1) = 180.0d0
          tmpa = erqm(2) - erqm(1)
          tmpb = (ftv(2,1) - ftv(1,1))/radian
          fwt(1) = 1.0d0/sqrt(1+(tmpa/tmpb)**2)
          if (nconf .gt. 2) then
             do i = 2, nconf-1
                tmpa = erqm(i+1) - erqm(i-1)
                tmpb = (ftv(i+1,1) - ftv(i-1,1))/radian
                fwt(i) = 1.0d0/sqrt(1+(tmpa/tmpb)**2)
             end do
          end if
          tmpa = erqm(nconf) - erqm(nconf-1)
          tmpb = (ftv(nconf,1) - ftv(nconf-1,1))/radian
          fwt(nconf) = 1.0d0/sqrt(1+(tmpa/tmpb)**2)
       end if
       do i = 1, nconf
          write(*,73) i,fwt(i)
   73     format(1x, 'Conformation ',i3,' Weight ',f8.4)
       end do
c
c      set initial values for torsions to be fitted
c
       ivxx = 0
       do i = 1, ntorfit
          j = ftorid(i)
          k = ctorid(j)
          do ii = 1, 6
             if(vflg(ii,i) .and. tflg(i) .eq. 0) then
                ivxx = ivxx +1
                if (ii .eq. 1) then
                   tors1(1,k) = vxx(ivxx)
                else if (ii .eq. 2) then
                   tors2(1,k) = vxx(ivxx)
                else if (ii .eq. 3) then
                   tors3(1,k) = vxx(ivxx)
                else if (ii .eq. 4) then
                   tors4(1,k) = vxx(ivxx)
                else if (ii .eq. 5) then
                   tors5(1,k) = vxx(ivxx)
                else if (ii .eq. 6) then
                   tors6(1,k) = vxx(ivxx)
                end if
                do jj = 1, ntorfit
                   kk = ctorid(ftorid(jj))
                   if(tflg(jj) .eq. i) then
                      if (ii .eq. 1) then
                         tors1(1,kk) = vxx(ivxx)
                      else if (ii .eq. 2) then
                         tors2(1,kk) = vxx(ivxx)
                      else if (ii .eq. 3) then
                         tors3(1,kk) = vxx(ivxx)
                      else if (ii .eq. 4) then
                         tors4(1,kk) = vxx(ivxx)
                      else if (ii .eq. 5) then
                         tors5(1,kk) = vxx(ivxx)
                      else if (ii .eq. 6) then
                         tors6(1,kk) = vxx(ivxx)
                      end if
                   end if
             end do
             end if
          end do
          ivxx = ivxx +1
          vcon = vxx(ivxx)
       end do
c
c      fitting the torsion parameters
c 
       maxstep = 1 
       avedl = 0.5d0 
       dowhile (avedl .gt. 0.1 .and. istep .lt. maxstep) 
          do i = 1, nconf
             ivxx = 0
             torf(i) = 0.0d0
             do j = 1, ntorfit
                jj = ftorid(j)
                kk = ctorid(jj) 
                ia = itors(1,kk)
                ib = itors(2,kk)
                ic = itors(3,kk)
                id = itors(4,kk)
                ita = class(ia)
                itb = class(ib)
                itc = class(ic)
                itd = class(id)
                tv = ftv(i,j)/radian
                tmpa = 
     &          tors1(1,kk)*(1+cos(tv))
     &          +tors2(1,kk)*(1-cos(2*tv))
     &          +tors3(1,kk)*(1+cos(3*tv))
     &          +tors4(1,kk)*(1-cos(4*tv))
     &          +tors5(1,kk)*(1+cos(5*tv))
     &          +tors6(1,kk)*(1-cos(6*tv))
                torf(i) = torf(i) + 0.5*tmpa             
                do ii = 1, 6
                   if(vflg(ii,j) .and. tflg(j) .eq. 0) then
                      ivxx = ivxx +1
                      coeff(i,ivxx) = 0.5*(1+(-1)**(ii+1)
     &                *cos(dble(ii)*tv))
                      do k = 1, ntorcrs
                         if(cflg(k) .gt. 0 
     &                     .and. cflg(k) .eq. jj) then
                            coeff(i,ivxx) = coeff(i,ivxx) 
     &                      +0.5*(1+(-1)**(ii+1)     
     &                      *cos(dble(ii)*ctv(i,k)/radian))
                         end if
                      end do
c                      write(*,33) i,ivxx,coeff(i,ivxx)
c   33                 format(1x,'Derivative ',2i4,f8.4)
                   end if
                end do
             end do
             torf(i) = torf(i) + vcon - delte(i)
             ivxx = ivxx +1
             coeff(i,ivxx) = 1.0d0
c             write(*,*) 'Energy Difference ',i,torf(i)
          end do
c
c         set Maxtrix elements for Matrix A
c
          do i = 1, nvxx
             do j = 1, nvxx
                tord(i,j) = 0.0d0
                do k = 1,nconf
                   tord(i,j) = tord(i,j)+coeff(k,i)*coeff(k,j)*fwt(k)
                end do
             end do
          end do
c         
c         Matrix A elements   
c  
          write(*,*) 'Total Variable Number ',nvxx
          write(*,*) 'Matrix A elements '
          do i = 1,nvxx
             do j = 1, nvxx
                 mata(i,j) = tord(i,j)
             end do
             write(*,23) (mata (i,j),j=1,nvxx)
   23        format(1x,5f12.4)
          end do            
c        
c         multiply vector: Yi*Coeff*Weight 
c
          do i = 1, nvxx
             torf(i) = 0.0d0
             do j = 1, nconf
                torf(i) = torf(i) + delte(j)*fwt(j)*coeff(j,i)
             end do
          end do
          do i = 1,nvxx
             mata(i,nvxx+1) = torf(i)
          end do
c 
c         solve the linear matrix equiation with Gauss-Jordan elimination 
c
          call gaussjordan (nvxx,mata)
c
c         get new torsion force constants
c
          do i = 1,nvxx
             vxx(i) = mata(i,nvxx+1)
          end do
          ivxx = 0
          do i = 1, ntorfit
             j = ftorid(i)
             k = ctorid(j)
             do ii = 1, 6
                if(vflg(ii,i) .and. tflg(i) .eq. 0) then
                   ivxx = ivxx +1
                   if (ii .eq. 1) then
                      tors1(1,k) = vxx(ivxx) 
                   else if (ii .eq. 2) then
                      tors2(1,k) = vxx(ivxx)
                   else if (ii .eq. 3) then
                      tors3(1,k) = vxx(ivxx)
                   else if (ii .eq. 4) then
                      tors4(1,k) = vxx(ivxx)  
                   else if (ii .eq. 5) then
                      tors5(1,k) = vxx(ivxx)
                   else if (ii .eq. 6) then 
                      tors6(1,k) = vxx(ivxx)
                   end if
                   do jj = 1, ntorcrs
                      kk = ctorid(jj)
                      if(cflg(j) .gt. 0 .and. 
     &                   cflg(jj) .eq. j) then
                         if (ii .eq. 1) then
                            tors1(1,kk) = vxx(ivxx)  
                         else if (ii .eq. 2) then
                            tors2(1,kk) = vxx(ivxx)
                         else if (ii .eq. 3) then
                            tors3(1,kk) = vxx(ivxx)
                         else if (ii .eq. 4) then
                            tors4(1,kk) = vxx(ivxx)  
                         else if (ii .eq. 5) then
                            tors5(1,kk) = vxx(ivxx)
                         else if (ii .eq. 6) then 
                            tors6(1,kk) = vxx(ivxx)
                         end if
                      end if
                   end do
                end if
             end do
             ivxx = ivxx +1
             vcon = vxx(ivxx)
          end do
          istep = istep +1
       end do 
c
c      validate the fitted results
c
       do i = 1, nconf
          call getref(i)
          kk = 0
          do j = 1, ntorfit
             k = ftorid(j)
             ia = torcrs(1,k)
             ib = torcrs(2,k)
             ic = torcrs(3,k)
             id = torcrs(4,k)
             ftv(i,j) = geometry(ia,ib,ic,id)
             if (tflg(j) .eq. 0) then
                kk = kk +1
                tfix(2,otfix+kk) = ftv(i,j)
                tfix(3,otfix+kk) = tfix(2,otfix+kk)
             end if
          end do
c
c      scale the coordinates of each active atom
c
         nvar = 0
         do j = 1, n
            if (use(j)) then
               nvar = nvar + 1
               xx(nvar) = x(j) * scale(nvar)
               nvar = nvar + 1
               xx(nvar) = y(j) * scale(nvar)
               nvar = nvar + 1
               xx(nvar) = z(j) * scale(nvar)
            end if
         end do
c
c      make the call to the optimization routine
c
         write (*,100)  i
  100    format (' Minimizing Structure',i4,' with New Parameters')
         coordtype = 'CARTESIAN'
         call lbfgs (nvar,xx,minimum,grdmin,minimiz2,optsave)
c
c     unscale the final coordinates for active atoms
c
         nvar = 0
         do j = 1, n
            if (use(j)) then
               nvar = nvar + 1
               x(j) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               y(j) = xx(nvar) / scale(nvar)
               nvar = nvar + 1
               z(j) = xx(nvar) / scale(nvar)
            end if
         end do
         emm(i) = energy ()
      end do
c
c     calculate the relative MM energies
c   
      emmmin = emm(1)
      do i = 2, nconf
         if (emm(i) .lt. emmmin) emmmin = emm(i)
      end do 
c
c     calcualte the energy difference and RMS
c
       rms = 0.0d0
       do i = 1, nconf
          ermm (i) = emm(i) - emmmin 
          delte (i) = erqm (i) - ermm(i) 
          rms = rms + delte(i)*delte(i)
          write (*,93)  i,ermm(i) 
       end do               
       rms = sqrt(rms/dble(nconf))
       write (*,*)  'Energy RMS With Fitting Parmeters = ',rms
       if (rms .gt. zrms ) then
          write (*,94)  zrms
   94     format (/,' Annihilating the torsions is better',
     &            /,' The Final RMS is :',f12.6,' kcal/mol',/)
       end if
c
c      output the fitted parameters
c
       filename = oldfilename
       leng = oldleng
       nkey = oldnkey
       do i = 1, nkey
          keyline(i) = oldkeyline(i)
       end do
c
c     output some definitions and parameters to a keyfile
c
      ikey = freeunit ()
      keyfile = filename(1:leng)//'.key'
      call version (keyfile,'new')
      open (unit=ikey,file=keyfile,status='new')
c
c     copy the contents of any previously existing keyfile
c
      do i = 1, nkey
         record = keyline(i)
         size = trimtext (record)
         write (ikey,120)  record(1:size)
  120    format (a)
      end do
c
c     list the valence parameters
c
      write (ikey,130)
  130 format (/,'#',/,'# Results of Valence Parameter Fitting',
     &              /,'#',/)

       write(iout,*)
       write(iout,*) 'Fitting Torsional Parameters:'
       do i = 1, ntorfit
          if (tflg(i) .eq. 0) then
             j = ftorid(i)
             k = ctorid(j)
             ia = itors(1,k)
             ib = itors(2,k)
             ic = itors(3,k)
             id = itors(4,k)
             ita = class(ia)
             itb = class(ib)
             itc = class(ic)
             itd = class(id)
             if (rms .gt. zrms) then
                tors1(1,k) = 0.0d0
                tors2(1,k) = 0.0d0
                tors3(1,k) = 0.0d0
             end if
             write(*,140) ita,itb,itc,itd,tors1(1,k),
     &       tors2(1,k),tors3(1,k)
             write(ikey,140) ita,itb,itc,itd,tors1(1,k),
     &       tors2(1,k),tors3(1,k)
  140        format(1x,'torsion ',4i4,f8.3,' 0.0 1 ',
     &       f8.3,' 180.0 2 ',f8.3,' 0.0 3')
          end if
       end do
       return
       end 
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine gaussjordan  --  Gauss-Jordan elimination  ##
c     ##                                                        ##
c     ############################################################
c
c     "gaussjordan" solves a system of linear equations by using
c     the method of Gaussian elimination with partial pivoting
c
c
      subroutine gaussjordan (n,a)
      implicit none
      integer maxn
      parameter (maxn=72)
      integer i,j,k,l,n
      real*8 t,a(maxn,maxn),av
c
c
      do k = 1, n-1
         av = 0.0d0
         do i = k, n
            if (abs(a(i,k)) .gt. abs(av)) then
               av = a(i,k)
               l=i
            end if
         end do
         if (abs(av) .lt. 1.0d-8) then
            write (*,*)  'Singular Coefficient Matrix'
            stop
         end if
         if (l .ne. k) then
            do j = k, n+1
               t = a(k,j)
               a(k,j) = a(l,j)
               a(l,j) = t
            end do
         end if
         av = 1.0d0 / av
         do j = k+1, n+1
            a(k,j) = a(k,j)*av
            do i = k+1,n
               a(i,j) = a(i,j) - a(i,k)*a(k,j)
            end do
         end do
      end do
      a(n,n+1) = a(n,n+1) / a(n,n)
      do k = 1, n-1
         i = n-k
         av = 0.0d0
         do j = i+1, n
            av = av + a(i,j)*a(j,n+1)
         end do
         a(i,n+1) = a(i,n+1) - av
      end do
      return
      end 
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  function minimiz2  --  energy and gradient for minimize  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "minimiz2" is a service routine that computes the energy and
c     gradient for a low storage BFGS optimization in Cartesian
c     coordinate space
c
c
      function minimiz2 (xx,g)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'scales.i'
      include 'usage.i'
      integer i,nvar
      real*8 minimiz2,e
      real*8 energy,eps
      real*8 xx(maxvar)
      real*8 g(maxvar)
      real*8 derivs(3,maxatm)
      logical analytic
c
c
c     use either analytical or numerical gradients
c
      analytic = .true.
      eps = 0.00001d0
c
c     translate optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            x(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar) / scale(nvar)
         end if
      end do
c
c     compute and store the energy and gradient
c
      if (analytic) then
         call gradient (e,derivs)
      else
         e = energy ()
         call numgrad (energy,derivs,eps)
      end if
      minimiz2 = e
c
c     store Cartesian gradient as optimization gradient
c
      nvar = 0
      do i = 1, n
         if (use(i)) then
            nvar = nvar + 1
            g(nvar) = derivs(1,i) / scale(nvar)
            nvar = nvar + 1
            g(nvar) = derivs(2,i) / scale(nvar)
            nvar = nvar + 1
            g(nvar) = derivs(3,i) / scale(nvar)
         end if
      end do
      return
      end

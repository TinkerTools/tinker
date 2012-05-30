c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 1991 by Shawn Huston & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program path  --  conformational interconversion pathway  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "path" locates a series of structures equally spaced along
c     a conformational pathway connecting the input reactant and
c     product structures; a series of constrained optimizations
c     orthogonal to the path is done via Lagrangian multipliers
c
c     literature reference:
c
c     R. Czerminski and R. Elber, "Reaction Path Study of
c     Conformational Transitions in Flexible Systems: Applications
c     to Peptides", Journal of Chemical Physics, 92, 5580-5601 (1990)
c
c
      program path
      implicit none
      include 'sizes.i'
      include 'align.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'linmin.i'
      include 'minima.i'
      include 'output.i'
      include 'paths.i'
      integer i,j,k,nvar
      integer ix,iy,iz
      integer ipath,npath
      real*8 rmsvalue,project
      real*8 epot,etot
      real*8 epot0,epot1
      real*8 sum,rplen
      real*8 path1
      real*8 grdmin,potnrg
      real*8 g(maxvar)
      real*8 p(maxvar)
      real*8 ge(maxvar)
      real*8 temp(maxvar,7)
      real*8 xtmp(maxatm)
      real*8 ytmp(maxatm)
      real*8 ztmp(maxatm)
      logical exist
      character*120 string
      external path1
      external optsave
c
c
c     initialize some constants and variables
c
      call initial
c
c     get and store the initial structure coordinates
c
      call getxyz
      do i = 1, n
         p0(3*i-2) = x(i)
         p0(3*i-1) = y(i)
         p0(3*i) = z(i)
         xtmp(i) = x(i)
         ytmp(i) = y(i)
         ztmp(i) = z(i)
      end do
c
c     get the coordinates for the final structure
c
      call getxyz
      call mechanic
c
c     set default values for some control variables
c
      nvar = 3 * n
      cyclesave = .true.
      stpmax = 1.0d0
      iwrite = 0
      if (verbose) then
         iprint = 1
      else
         iprint = 0
      end if
c
c     get the number of path points to be generated
c
      npath = -1
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  npath
   10 continue
      if (npath .le. 0) then
         write (iout,20)
   20    format (/,' Enter Number of Path Points to Generate [9] :  ',$)
         read (input,30)  npath
   30    format (i10)
      end if
      if (npath .le. 0)  npath = 9
c
c     get the termination criterion as RMS gradient along path
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  grdmin
   40 continue
      if (grdmin .le. 0.0d0) then
         write (iout,50)
   50    format (/,' Enter RMS Gradient per Atom Criterion',
     &              ' [0.1] :  ',$)
         read (input,60)  grdmin
   60    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.1d0
c
c     superimpose the reactant and product structures
c
      nfit = n
      do i = 1, n
         ifit(1,i) = i
         ifit(2,i) = i
         wfit(i) = mass(i)
      end do
      call impose (n,xtmp,ytmp,ztmp,n,x,y,z,rmsvalue)
      write (iout,70)  rmsvalue
   70 format (/,' RMS Fit for Reactant and Product :',f12.4)
c
c     store the coordinates for the superimposed product
c
      do i = 1, n
         p1(3*i-2) = x(i)
         p1(3*i-1) = y(i)
         p1(3*i) = z(i)
      end do
c
c     write out the starting potential energy values
c
      epot0 = potnrg (p0,g)
      epot1 = potnrg (p1,g)
      write (iout,80)  epot0,epot1
   80 format (/,' Reactant Potential Energy :',f12.4,
     &        /,' Product Potential Energy : ',f12.4)
c
c     construct step vector for getting
c     optimization-initial coordinates
c
      rplen = npath + 1
      pnorm = 0.0d0
      do i = 1, nvar
         pvect(i) = p1(i) - p0(i)
         pstep(i) = pvect(i) / rplen
         pnorm = pnorm + pvect(i)**2
      end do
      pnorm = sqrt(pnorm)
c
c     set the gradient of constraints array
c
      do i = 1, n
         ix = 3*(i-1) + 1
         iy = ix + 1
         iz = iy + 1
         gc(ix,1) = pvect(ix)
         gc(iy,1) = pvect(iy)
         gc(iz,1) = pvect(iz)
         gc(ix,2) = mass(i)
         gc(iy,2) = 0.0d0
         gc(iz,2) = 0.0d0
         gc(ix,3) = 0.0d0
         gc(iy,3) = mass(i)
         gc(iz,3) = 0.0d0
         gc(ix,4) = 0.0d0
         gc(iy,4) = 0.0d0
         gc(iz,4) = mass(i)
         gc(ix,5) = 0.0d0
         gc(iy,5) = mass(i) * p0(iz)
         gc(iz,5) = -mass(i) * p0(iy)
         gc(ix,6) = -mass(i) * p0(iz)
         gc(iy,6) = 0.0d0
         gc(iz,6) = mass(i) * p0(ix)
         gc(ix,7) = mass(i) * p0(iy)
         gc(iy,7) = -mass(i) * p0(ix)
         gc(iz,7) = 0.0d0
      end do
c
c     copy to temporary storage and orthogonalize
c
      do i = 1, 7
         do j = 1, nvar
            temp(j,i) = gc(j,i)
         end do
      end do
      call orthog (nvar,maxvar,7,gc)
c
c     set the A matrix to transform sigma into C space
c
      do i = 1, 7
         do k = 1, 7
            sum = 0.0d0
            do j = 1, nvar
               sum = sum + temp(j,i)*gc(j,k)
            end do
            acoeff(i,k) = sum
         end do
      end do
c
c     perform the matrix inversion to get A matrix
c     which transforms C into sigma space
c
      call invert (7,acoeff)
c
c     set the current path point to be the reactant
c
      do i = 1, nvar
         p(i) = p0(i)
      end do
c
c     loop over structures along path to be optimized
c
      do ipath = 1, npath
         write (iout,90)  ipath
   90    format (/,' Path Point :',i12)
c
c     get r(zeta), set initial path point and energy
c
         do i = 1, nvar
            pzet(i) = p0(i) + ipath*pstep(i)
            p(i) = p(i) + pstep(i)
         end do
         epot = potnrg (p,ge)
         write (iout,100)  epot
  100    format (' Initial Point :',12x,f12.4)
c
c     call optimizer to get constrained minimum
c
         call lbfgs (nvar,p,etot,grdmin,path1,optsave)
c        call ocvm (nvar,p,etot,grdmin,path1,optsave)
c
c     print energy and constraint value at the minimum
c
         epot = potnrg (p,ge)
         write (iout,110)  epot
  110    format (' Optimized Point :',10x,f12.4)
         write (iout,120)  etot-epot
  120    format (' Target-Energy Difference :',d13.3)
c
c     write coordinates of the current path point
c
         call optsave (ipath,epot,p)
c
c     find projection of the gradient along path direction
c
         project = 0.0d0
         do i = 1, nvar
            project = project + ge(i)*pvect(i)/pnorm
         end do
         write (iout,130)  project
  130    format (' Gradient along Path :',6x,f12.4)
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function path1  --  value and gradient of target function  ##
c     ##                                                             ##
c     #################################################################
c
c
      function path1 (p,gt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'paths.i'
      integer i,j,nvar
      integer ix,iy,iz
      real*8 xx,yy,zz
      real*8 path1,cterm
      real*8 potnrg
      real*8 gamma(7)
      real*8 cnst(7)
      real*8 sigma(7)
      real*8 p(maxvar)
      real*8 gt(maxvar)
      real*8 ge(maxvar)
c
c
c     get the value of the potential energy
c
      nvar = 3 * n
      path1 = potnrg (p,ge)
c
c     construct the Lagrangian multipliers
c
      do i = 1, 7
         gamma(i) = 0.0d0
         do j = 1, nvar
            gamma(i) = gamma(i) - ge(j)*gc(j,i)
         end do
      end do
c
c     set the path value, translation and rotation constraints
c
      do i = 1, 7
         cnst(i) = 0.0d0
      end do
      do i = 1, n
         ix = 3*(i-1) + 1
         iy = ix + 1
         iz = iy + 1
         xx = p(ix) - pzet(ix)
         yy = p(iy) - pzet(iy)
         zz = p(iz) - pzet(iz)
         cnst(1) = cnst(1) + xx*pvect(ix) + yy*pvect(iy) + zz*pvect(iz)
         cnst(2) = cnst(2) + mass(i) * (p(ix)-p0(ix))
         cnst(3) = cnst(3) + mass(i) * (p(iy)-p0(iy))
         cnst(4) = cnst(4) + mass(i) * (p(iz)-p0(iz))
         cnst(5) = cnst(5) + mass(i) * (p(iy)*p0(iz)-p(iz)*p0(iy))
         cnst(6) = cnst(6) + mass(i) * (p(iz)*p0(ix)-p(ix)*p0(iz))
         cnst(7) = cnst(7) + mass(i) * (p(ix)*p0(iy)-p(iy)*p0(ix))
      end do
c
c     construct the orthonormal "sigma" constraints
c
      do i = 1, 7
         sigma(i) = 0.0d0
         do j = 1, 7
            sigma(i) = sigma(i) + acoeff(i,j)*cnst(j)
         end do
      end do
c
c     find the target function value
c
      cterm = 0.0d0
      do i = 1, 7
         cterm = cterm + gamma(i)*sigma(i)
      end do
      path1 = path1 + cterm
c
c     construct the gradient of the target function
c
      do i = 1, nvar
         gt(i) = ge(i)
         do j = 1, 7
            gt(i) = gt(i) + gamma(j)*gc(i,j)
         end do
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  function potnrg  --  potential energy and gradient  ##
c     ##                                                      ##
c     ##########################################################
c
c
      function potnrg (xx,g)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'paths.i'
      integer i
      real*8 energy
      real*8 potnrg
      real*8 xx(maxvar)
      real*8 g(maxvar)
      real*8 deriv(3,maxatm)
c
c
c     copy position vector into atomic coordinates
c
      do i = 1, n
         x(i) = xx(3*i-2)
         y(i) = xx(3*i-1)
         z(i) = xx(3*i)
      end do
c
c     compute potential energy and Cartesian derivatives
c
      call gradient (energy,deriv)
c
c     set the energy value and gradient vector
c
      potnrg = energy
      do i = 1, n
         g(3*i-2) = deriv(1,i)
         g(3*i-1) = deriv(2,i)
         g(3*i) = deriv(3,i)
      end do
      return
      end

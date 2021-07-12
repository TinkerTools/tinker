c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2001 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program polarize  --  compute the molecular polarizability  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "polarize" computes the molecular polarizability by applying
c     an external field along each axis followed by diagonalization
c     of the resulting polarizability tensor
c
c
      program polarize
      use atoms
      use inform
      use iounit
      use molcul
      use mpole
      use polar
      use polpot
      use potent
      implicit none
      integer i
      real*8 addu,malpha
      real*8 external
      real*8 exfield(3)
      real*8 umol(3)
      real*8 umol0(3)
      real*8 dalpha(3)
      real*8 alpha(3,3)
      real*8 valpha(3,3)
      character*40 fstr
c
c
c     get the coordinates and required force field parameters
c
      call initial
      call getxyz
      call attach
      call field
      call cutoffs
      call katom
      call molecule
      call kmpole
      call kpolar
      call kchgtrn
      call mutate
c
c     sum atomic polarizabilities to get additive molecular value
c
      if (.not. use_polar) then
         write (iout,10)
   10    format (/,' POLARIZE  --  Dipole Polarizability',
     &              ' is Not in Use')
         call fatal
      end if
      addu = 0.0d0
      do i = 1, npole
         addu = polarity(i) + addu
      end do
      fstr = ' Additive Total Polarizability :    '
      if (nmol .eq. 1)  fstr = ' Additive Molecular Polarizability :'
      if (digits .ge. 8) then
         write (iout,20)  fstr(1:36),addu
   20    format (/,a36,f20.8)
      else if (digits .ge. 6) then
         write (iout,30)  fstr(1:36),addu
   30    format (/,a36,f18.6)
      else
         write (iout,40)  fstr(1:36),addu
   40    format (/,a36,f16.4)
      end if
c
c     find induced dipoles in absence of an external field
c
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      call moluind (exfield,umol0)
c
c     compute each column of the polarizability tensor
c
      external = 0.01d0
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(1) = external
      call moluind (exfield,umol)
      alpha(1,1) = (umol(1)-umol0(1)) / exfield(1)
      alpha(2,1) = (umol(2)-umol0(2)) / exfield(1)
      alpha(3,1) = (umol(3)-umol0(3)) / exfield(1)
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(2) = external
      call moluind (exfield,umol)
      alpha(1,2) = (umol(1)-umol0(1)) / exfield(2)
      alpha(2,2) = (umol(2)-umol0(2)) / exfield(2)
      alpha(3,2) = (umol(3)-umol0(3)) / exfield(2)
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(3) = external
      call moluind (exfield,umol)
      alpha(1,3) = (umol(1)-umol0(1)) / exfield(3)
      alpha(2,3) = (umol(2)-umol0(2)) / exfield(3)
      alpha(3,3) = (umol(3)-umol0(3)) / exfield(3)
c
c     print out the full polarizability tensor
c
      fstr = ' Total Polarizability Tensor :    '
      if (nmol .eq. 1)  fstr = ' Molecular Polarizability Tensor :'
      write (iout,50)  fstr(1:34)
   50 format (/,a34,/)
      if (digits .ge. 8) then
         write (iout,60)  alpha(1,1),alpha(1,2),alpha(1,3),
     &                    alpha(2,1),alpha(2,2),alpha(2,3),
     &                    alpha(3,1),alpha(3,2),alpha(3,3)
   60    format (13x,3f17.8,/,13x,3f17.8,/,13x,3f17.8)
      else if (digits .ge. 6) then
         write (iout,70)  alpha(1,1),alpha(1,2),alpha(1,3),
     &                    alpha(2,1),alpha(2,2),alpha(2,3),
     &                    alpha(3,1),alpha(3,2),alpha(3,3)
   70    format (13x,3f15.6,/,13x,3f15.6,/,13x,3f15.6)
      else
         write (iout,80)  alpha(1,1),alpha(1,2),alpha(1,3),
     &                    alpha(2,1),alpha(2,2),alpha(2,3),
     &                    alpha(3,1),alpha(3,2),alpha(3,3)
   80    format (13x,3f13.4,/,13x,3f13.4,/,13x,3f13.4)
      end if
c
c     diagonalize the tensor and get molecular polarizability
c
      call jacobi (3,alpha,dalpha,valpha)
      fstr = ' Polarizability Tensor Eigenvalues :'
      write (iout,90)  fstr(1:36)
   90 format (/,a36,/)
      if (digits .ge. 8) then
         write (iout,100)  dalpha(1),dalpha(2),dalpha(3)
  100    format (13x,3f17.8)
      else if (digits .ge. 6) then
         write (iout,110)  dalpha(1),dalpha(2),dalpha(3)
  110    format (13x,3f15.6)
      else
         write (iout,120)  dalpha(1),dalpha(2),dalpha(3)
  120    format (13x,3f13.4)
      end if
      malpha = (dalpha(1)+dalpha(2)+dalpha(3)) / 3.0d0
      fstr = ' Interactive Total Polarizability :    '
      if (nmol .eq. 1)  fstr = ' Interactive Molecular Polarizability :'
      if (digits .ge. 8) then
         write (iout,130)  fstr(1:39),malpha
  130    format (/,a39,f17.8)
      else if (digits .ge. 6) then
         write (iout,140)  fstr(1:39),malpha
  140    format (/,a39,f15.6)
      else
         write (iout,150)  fstr(1:39),malpha
  150    format (/,a39,f13.4)
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine moluind  --  molecular induced dipole in field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "moluind" computes the molecular induced dipole components
c     in the presence of an external electric field
c
c
      subroutine moluind (exfield,umol)
      use atoms
      use inform
      use iounit
      use mpole
      use polar
      use polopt
      use polpcg
      use polpot
      use units
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 eps,epsold
      real*8 polmin
      real*8 a,b,sum,term
      real*8 norm,exmax
      real*8 umol(3)
      real*8 exfield(3)
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: usum(:,:)
      logical header,done
      logical dodfield
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (poli(npole))
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (rsd(3,npole))
      allocate (zrsd(3,npole))
      allocate (conj(3,npole))
      allocate (vec(3,npole))
c
c     check for chiral multipoles and rotate to global frame
c
      call chkpole
      call rotpole
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
c
c     get the electrostatic field due to permanent multipoles
c
      dodfield = .true.
      if (dodfield)  call dfield0a (field,fieldp)
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
         end do
      end do
c
c     increment induced dipoles to account for external field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = udir(j,i) + polarity(i)*exfield(j)
            uind(j,i) = udir(j,i)
         end do
      end do
c
c     get induced dipoles via the OPT extrapolation method
c
      if (poltyp .eq. 'OPT') then
         do i = 1, npole
            if (douind(ipole(i))) then
               do j = 1, 3
                  uopt(0,j,i) = udir(j,i)
               end do
            end if
         end do
         do k = 1, optorder
            optlevel = k - 1
            call ufield0a (field,fieldp)
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     uopt(k,j,i) = polarity(i) * field(j,i)
                     uind(j,i) = uopt(k,j,i)
                  end do
               end if
            end do
         end do
         allocate (usum(3,n))
         do i = 1, npole
            if (douind(ipole(i))) then
               do j = 1, 3
                  uind(j,i) = 0.0d0
                  usum(j,i) = 0.0d0
                  do k = 0, optorder
                     usum(j,i) = usum(j,i) + uopt(k,j,i)
                     uind(j,i) = uind(j,i) + copt(k)*usum(j,i)
                  end do
               end do
            end if
         end do
         deallocate (usum)
      end if
c
c     compute mutual induced dipole moments via CG algorithm
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
         call ufield0a (field,fieldp)
         do i = 1, npole
            poli(i) = max(polmin,polarity(i))
            do j = 1, 3
               rsd(j,i) = field(j,i)
               zrsd(j,i) = rsd(j,i) * poli(i)
               conj(j,i) = zrsd(j,i)
            end do
         end do
c
c     iterate the mutual induced dipoles and check convergence
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               do j = 1, 3
                  vec(j,i) = uind(j,i)
                  uind(j,i) = conj(j,i)
               end do
            end do
            call ufield0a (field,fieldp)
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
               end do
            end do
            a = 0.0d0
            sum = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
               end do
            end do
            b = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  zrsd(j,i) = rsd(j,i) * poli(i)
                  b = b + rsd(j,i)*zrsd(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            eps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  eps = eps + rsd(j,i)*rsd(j,i)
               end do
            end do
            eps = debye * sqrt(eps/dble(npolar))
            epsold = eps
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debye)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. politer)  done = .true.
c
c     apply a "peek" iteration to the mutual induced dipoles
c
            if (done) then
               do i = 1, npole
                  if (douind(ipole(i))) then
                     term = pcgpeek * poli(i)
                     do j = 1, 3
                        uind(j,i) = uind(j,i) + term*rsd(j,i)
                     end do
                  end if
               end do
            end if
         end do
c
c     print a warning if induced dipoles failed to converge
c
         if (iter.ge.maxiter .or. eps.gt.epsold) then
            write (iout,30)
   30       format (/,' MOLUIND  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call fatal
         end if
      end if
c
c     sum up the total molecular induced dipole components
c
      do j = 1, 3
         umol(j) = 0.0d0
      end do
      do i = 1, npole
         umol(1) = umol(1) + uind(1,i)
         umol(2) = umol(2) + uind(2,i)
         umol(3) = umol(3) + uind(3,i)
      end do
c
c     print out a list of the final induced dipole moments
c
      if (verbose) then
         exmax = max(exfield(1),exfield(2),exfield(3))
         if (dodfield .or. exmax.ne.0.0d0) then
            write (iout,40)  (exfield(j),j=1,3)
   40       format (/,' Applied External Field :',//,13x,3f13.4)
            header = .true.
            do i = 1, npole
               if (polarity(i) .ne. 0.0d0) then
                  if (header) then
                     header = .false.
                     write (iout,50)
   50                format (/,' Induced Dipole Moments (Debye) :')
                     write (iout,60)
   60                format (/,4x,'Atom',15x,'X',12x,'Y',12x,'Z',
     &                          11x,'Total',/)
                  end if
                  norm = sqrt(uind(1,i)**2+uind(2,i)**2+uind(3,i)**2)
                  write (iout,70)  ipole(i),(debye*uind(j,i),j=1,3),
     &                             debye*norm
   70             format (i8,5x,3f13.4,1x,f13.4)
               end if
            end do
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (poli)
      deallocate (field)
      deallocate (fieldp)
      deallocate (rsd)
      deallocate (zrsd)
      deallocate (conj)
      deallocate (vec)
      return
      end

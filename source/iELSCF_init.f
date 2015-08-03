c
c
c     ######################################################
c     ##  COPYRIGHT (C)  2014  by  Alex Albaugh (THG Lab) ##
c     ##              No Rights Reserved                  ##
c     ######################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dinit  --  initialize Drude polarization      ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dinit" initializes the proper settings for Drude polarization
c     along with Drude data structures and initial 
c
c
      subroutine iELSCF_init(dt)
      use atoms
      use polar
      use mdstuf
      use keys
      use bound
      use iounit
      use atomid
      use usage
      use units
      use moldyn
      use mpole
      use files
      use iELSCF
      use uprior
      implicit none
      integer i,j
      integer next
      real*8 dt
      real*8 ekt_aux
      real*8 qterm_aux
      real*8 speed,maxwell2,maxwell
      real*8 vec(3)
      character*20 keyword
      character*120 record
      character*120 string
      real*8 eksum,ekin,temp
      integer idynaux
      logical exist
      character*120 dynfileaux
      integer freeunit

c
c     Allocate arrays.
c 
      allocate(a_aux(3,n))
      allocate(v_aux(3,n))
      allocate(uind_aux(3,n))
      allocate(uinp_aux(3,n))
      
      do i = 1, n
         do j = 1, 3
            uind_aux(j,i) = 0.0d0
            uinp_aux(j,i) = 0.0d0
         end do
      end do
      
      first = .false.
      

c
c     set default parameters for the Drude dynamics
c 
      if(integrate .ne. 'VERLET') then
         integrate = 'VERLET'
         print*, "iEL-SCF only available with Verlet integration."
      end if
      if(use_pred) then
         use_pred = .false.
         print*, "PCG predictor not available with iEL-SCF."
      end if
      auxstat = 'NONE'!BERENDSEN RESCALE NOSE-HOOVER
      omega = dsqrt(2.0d0)
c
c     set default values for Drude temperature control
c 
      aux_kelvin = 100000.0d0!e**2 * Ang**2 / ps**2
      aux_tautemp = 0.1d0!ELtautemp in ps
      do i = 1, maxnose
         vnh_aux(i) = 0.0d0
         qnh_aux(i) = 0.0d0
         gnh_aux(i) = 0.0d0
         pnh_aux(i) = 0.0d0
         
         pnh(i) = 0.0d0
      end do
c
c     check for keywords containing any altered Drude parameters
c 
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:13) .eq. 'AUX-TAU-TEMP ') then
            read (string,*,err=10,end=10)  aux_tautemp
         else if (keyword(1:9) .eq. 'AUX-TEMP ') then
            read (string,*,err=10,end=10)  aux_kelvin
         else if (keyword(1:6) .eq. 'OMEGA ') then
            read (string,*,err=10,end=10) omega
         else if (keyword(1:8) .eq. 'AUXSTAT ') then
            read (string,*,err=10,end=10) auxstat
         end if
   10    continue
      end do
      
      omega = omega/dt
      print*,omega,omega*dt
      
      if( (auxstat .ne. 'NOSE-HOOVER')  .and.
     &    (auxstat .ne. 'NOSE-HOOVER1') .and.
     &    (auxstat .ne. 'RESCALE')      .and.
     &    (auxstat .ne. 'BERENDSEN')    .and.
     &    (auxstat .ne. 'NONE') )       then
         print*, "Auxiliary thermostat not available. Reset to 'NONE'."
         auxstat = 'NONE'
      end if

c
c     set masses for Nose-Hoover thermostat and barostat
c 
      auxDoF = 3 * n
      ekt_aux = aux_kelvin!e**2 * Ang**2 / ps**2
      qterm_aux = ekt_aux * aux_tautemp * aux_tautemp !e**2 * Ang**2
      do j = 1, maxnose
         if (qnh_aux(j) .eq. 0.0d0) qnh_aux(j) = qterm_aux
      end do
      qnh_aux(1) = dble(auxDoF) * qnh_aux(1)!e**2 * Ang**2
c
c     set velocities and accelerations for auxiliary dipoles
c 
      dynfileaux = filename(1:leng)//'.auxdyn'
      call version (dynfileaux,'old')
      inquire (file=dynfileaux,exist=exist)
      if (exist) then
         print*,"RESTARTING AUXILIARY DIPOLES FROM SAVED FILES."
         idynaux = freeunit ()
         open (unit=idynaux,file=dynfileaux,status='old')
         rewind(unit=idynaux)
         call readdynaux (idynaux)
         close (unit=idynaux)
         do i = 1, n
            do j = 1, 3
               uinp(j,i) = uind(j,i)
               uinp_aux(j,i) = uind_aux(j,i)!Water approx.
            end do
         end do
      else
         first = .true.
         call induce
         first = .false.
         
         do i = 1, n
            do j = 1, 3
               uind_aux(j,i) = uind(j,i)
               uinp_aux(j,i) = uinp(j,i)
            end do
         end do
         
         if(.not. isothermal) then
            print*,"NVE velocities initialized at 298.0 K."
         end if
         
         do i = 1, n
            if(.not. isothermal) then
               speed = maxwell (mass(i),298.0d0)
               call ranvec (vec)
               do j = 1, 3
                  v(j,i) = speed*vec(j)
               end do
            end if
            
            if(auxstat .ne. 'NONE') then
               speed = maxwell2 (aux_kelvin)!ELmass in ps**2/Ang**3, factors convert to g/mol/electron**2
               call ranvec (vec)
               do j = 1, 3
                  v_aux(j,i) = speed*vec(j)
                  a_aux(j,i) = 0.0d0
               end do
            else
               do j = 1, 3
                  v_aux(j,i) = 0.0d0
                  a_aux(j,i) = 0.0d0
               end do
            end if
         end do
         call auxkinetic(.true.)
      end if
      
      return
      end

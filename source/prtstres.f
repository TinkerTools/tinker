c
c
c     ###################################################
c     ##                                               ##
c     ##                                               ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine prtstres  --  output of stress tensor          ##
c     ##                                                            ##
c     ################################################################
c
c
c     "prtstres" writes stress tensor components to an output file
c
c
      subroutine prtstres (istep,stress)
      use bound
      use boxes
      use files
      use inform
      use iounit
      use output
      use strvar
      use titles
      implicit none
      integer istep
      integer istr
      integer modstre
      real*8 stress(3,3)
      logical exist
      character*240 strfile
c
c
c     only necessary if user requests stress tensor output
c
      if (.not. stresav)  return
c
c     check number of steps between trajectory file saves
c
      modstre = mod(istep,istress)
      if (modstre .ne. 0)  return
c
c     open and format the output file 
c
      istr = freeunit ()
      strfile = filename(1:leng)
      call suffix (strfile,'str','old')
      inquire (file=strfile,exist=exist)
      if (exist) then
         call openend (istr,strfile)
      else
         open (unit=istr,file=strfile,status='new')
         write (istr, 10)
   10    format ('MD STEP',4x,
     &           'STRESS(1,1)',4x,'STRESS(2,1)',4x,'STRESS(3,1)',4x,
     &           'STRESS(1,2)',4x,'STRESS(2,2)',4x,'STRESS(3,2)',4x,
     &           'STRESS(1,3)',4x,'STRESS(2,3)',4x,'STRESS(3,3)')
      end if
c
c     write the value of the stress tensor components
c
      write (istr, 20) istep,stress
   20 format (i7,9f15.4)
c
c     close the stress tensor output file
c      
      close (unit=istr)
      return
      end
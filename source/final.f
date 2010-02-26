c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine final  --  final actions before program exit  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "final" performs any final program actions, prints a status
c     message, and then pauses if necessary to avoid closing the
c     execution window
c
c
      subroutine final
      implicit none
      include 'sizes.i'
      include 'inform.i'
      include 'iounit.i'
      include 'socket.i'
      include 'solute.i'
c
c
c     close any open socket used for external communication
c
      if (use_socket) then
         call sktkill
      end if
c
c     free memory used by the APBS Poisson-Boltzmann solver
c
      if (solvtyp .eq. 'PB') then
         call apbsfinal
      end if
c
c     print a final status message before exiting TINKER
c
      if (debug) then
         write (iout,10)
   10    format (/,' TINKER is Exiting following Normal Termination',
     &              ' of the Program',/)
      end if
c
c     may need a pause to avoid closing the execution window
c
      if (holdup) then
         read (input,20)
   20    format ()
      end if
      return
      end

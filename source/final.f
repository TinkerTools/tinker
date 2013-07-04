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
      include 'chunks.i'
      include 'disgeo.i'
      include 'inform.i'
      include 'iounit.i'
      include 'neigh.i'
      include 'paths.i'
      include 'pme.i'
      include 'socket.i'
      include 'solute.i'
      include 'uprior.i'
      include 'usage.i'
      include 'usolve.i'
      include 'vibs.i'
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
c     perform deallocation of associated pointers arrays
c
      if (associated(pmetable))  deallocate (pmetable)
      if (associated(bnd))  deallocate (bnd)
      if (associated(georad))  deallocate (georad)
      if (associated(xvold))  deallocate (xvold)
      if (associated(yvold))  deallocate (yvold)
      if (associated(zvold))  deallocate (zvold)
      if (associated(xcold))  deallocate (xcold)
      if (associated(ycold))  deallocate (ycold)
      if (associated(zcold))  deallocate (zcold)
      if (associated(xmold))  deallocate (xmold)
      if (associated(ymold))  deallocate (ymold)
      if (associated(zmold))  deallocate (zmold)
      if (associated(nvlst))  deallocate (nvlst)
      if (associated(vlst))  deallocate (vlst)
      if (associated(nelst))  deallocate (nelst)
      if (associated(elst))  deallocate (elst)
      if (associated(nulst))  deallocate (nulst)
      if (associated(ulst))  deallocate (ulst)
      if (associated(pc0))  deallocate (pc0)
      if (associated(pc1))  deallocate (pc1)
      if (associated(pvect))  deallocate (pvect)
      if (associated(pstep))  deallocate (pstep)
      if (associated(pzet))  deallocate (pzet)
      if (associated(gc))  deallocate (gc)
      if (associated(thetai1))  deallocate (thetai1)
      if (associated(thetai2))  deallocate (thetai2)
      if (associated(thetai3))  deallocate (thetai3)
      if (associated(qgrid))  deallocate (qgrid)
      if (associated(qfac))  deallocate (qfac)
      if (associated(udalt))  deallocate (udalt)
      if (associated(upalt))  deallocate (upalt)
      if (associated(usalt))  deallocate (usalt)
      if (associated(upsalt))  deallocate (upsalt)
      if (associated(iuse))  deallocate (iuse)
      if (associated(use))  deallocate (use)
      if (associated(mindex))  deallocate (mindex)
      if (associated(minv))  deallocate (minv)
      if (associated(phi))  deallocate (phi)
      if (associated(phik))  deallocate (phik)
      if (associated(pwork))  deallocate (pwork)
c
c     print a final status message before exiting TINKER
c
      if (debug) then
         write (iout,10)
   10    format (/,' TINKER is Exiting following Normal Termination',
     &              ' of the Program',/)
      end if
c
c     gracious exit from MPI
c
c      call cleanup_nlist_builder
c      call mpi_cleanup
c
c     may need a pause to avoid closing the execution window
c
      if (holdup) then
         read (input,20)
   20    format ()
      end if
      return
      end

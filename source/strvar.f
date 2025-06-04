c
c
c     ###################################################
c     ##                                               ##
c     ##                                               ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module strvar -- stress tensor output control parameters   ##
c     ##                                                             ##
c     #################################################################
c
c
c     stresav   logical flag to save stress tensor components
c     stresfrq  time (in femtoseconds) between stress tensor prints 
c     istress   steps between stress tensor component prints   
c
c
      module strvar
      implicit none
      integer istress
      real*8 stresfrq
      logical stresav
      save
      end
      

      module constantsmodule
        implicit none
        integer max_order
        real*8 kfc,min_gvol
        real*8 volmina,volminb
        parameter(max_order = 16)
        parameter (kfc = 2.2269859253d0)
        parameter (min_gvol = 0.0d0)
        parameter (volmina = 0.01d0)
        parameter (volminb = 0.1d0)
      end module constantsmodule

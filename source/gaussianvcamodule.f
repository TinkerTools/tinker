      module gaussianvcamodule
         implicit none
         type :: GaussianVca
            real*8 :: v ! Gaussian volume
            real*8 :: a ! Gaussian exponent
            real*8, dimension(3) :: c ! center
         end type GaussianVca
      end module gaussianvcamodule

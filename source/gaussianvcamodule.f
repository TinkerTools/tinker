      module gaussianvcamodule
        use constantsmodule
        use math
        type :: GaussianVca
          real*8 :: v ! Gaussian volume
          real*8 :: a ! Gaussian exponent
          real*8, dimension(3) :: c ! center
        end type GaussianVca

        contains

        ! overlap between two Gaussians represented by a (V,c,a) triplet
        !   V: volume of Gaussian
        !   c: position of Gaussian
        !   a: exponential coefficient
        !
        !   g(x) = V (a/pi)^(3/2) exp(-a(x-c)^2)
        !
        !   this version is based on V=V(V1,V2,r1,r2,alpha)
        !   alpha = (a1 + a2)/(a1 a2)
        !
        !   dVdr is (1/r)*(dV12/dr)
        !   dVdV is dV12/dV1 
        !   dVdalpha is dV12/dalpha
        !   d2Vdalphadr is (1/r)*d^2V12/dalpha dr
        !   d2VdVdr is (1/r) d^2V12/dV1 dr
        function ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp) 
     &      result(sgvol)
          class(GaussianVca), intent(in) :: g1, g2
          class(GaussianVca), intent(out) :: g12
          real*8, intent(out) :: dVdr, dVdV, sfp
          real*8 :: d2, deltai, gvol, p12, a12
          real*8 :: s, sp, df, dgvol, dgvolv, ef, dgvola2, dgvola1
          real*8 :: dgalpha, dgalpha2, dgvolvdr, sgvol
          real*8, dimension(3) :: c1, c2, dist

          c1 = g1%c
          c2 = g2%c
          dist = c2 - c1
          d2 = dot_product(dist, dist)

          a12 = g1%a + g2%a
          deltai = 1.0d0 / a12
          df = (g1%a) * (g2%a) * deltai

          ef = exp(-df * d2)
          gvol = ((g1%v * g2%v) / (PI / df)**1.5d0) * ef
          dgvol = -2.0d0 * df * gvol
          dgvolv = 0.0d0
          if (g1%v > 0.0d0) then
            dgvolv = gvol / g1%v
          end if

          g12%c = ((c1 * g1%a) + (c2 * g2%a)) * deltai
          g12%a = a12
          g12%v = gvol

          s = pol_switchfunc(gvol, volmina, volminb, sp)
          sfp = sp * gvol + s
          dVdr = dgvol
          dVdV = dgvolv

          sgvol = s * gvol
        end function ogauss_alpha

      end module gaussianvcamodule

      ! Overlap volume switching function + 1st derivative 
      function pol_switchfunc(gvol, volmina, volminb, sp) result(s)
        real*8, intent(in) :: gvol, volmina, volminb
        real*8, intent(out) :: sp
        real*8 :: swf, swfp, swd, swu, swu2, swu3, s

        swf = 0.0d0
        swfp = 1.0d0
        if (gvol > volminb) then
          swf = 1.0d0
          swfp = 0.0d0
        else if (gvol < volmina) then
          swf = 0.0d0
          swfp = 0.0d0
        end if

        swd = 1.0d0 / (volminb - volmina)
        swu = (gvol - volmina) * swd
        swu2 = swu * swu
        swu3 = swu * swu2
        s = swf + swfp * swu3 * (10.0d0 - 15.0d0 * swu + 6.0d0 * swu2)
        sp = swfp * swd * 30.0d0 * swu2 * (1.0d0 - 2.0d0 * swu + swu2)

        ! turn off switching function
        !sp = 0.0d0
        !s = 1.0d0
      end function pol_switchfunc

      function dot_product(a, b) result(dot)
        real*8, dimension(:), intent(in) :: a, b
        real*8 :: dot
        integer :: i
      
        dot = 0.0d0
        do i = 1, size(a)
          dot = dot + a(i) * b(i)
        end do
      end function dot_product

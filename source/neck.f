c
c
c     ############################################################
c     ##  COPYRIGHT (C)  2023  by Rae Corrigan & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine neck  --  neck contribution to effective radii  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "neck" calculates the descreening contribution of any neck 
c     formed between two atoms
c
c     literature reference:
c
c     B. Aguilar, R. Shadrach, and A. V. Onufriev, "Reducing the 
c     Secondary Structure Bias in the Generalized Born Model via 
c     R6 Effective Radii", Journal of Chemical Theory and 
c     Computation, 6, 3613-3630 (2010)
c
c     variables and parameters:
c
c     r          separation distance between two atoms
c     intstarti  start of descreening integral for atom i
c     desck      descreening radius of atom k
c     mixsn      mixed Sneck scale factor for atoms i and k
c
c
      subroutine neck (r,intstarti,desck,mixsn,neckval)
      use math
      use solute
      implicit none
      integer i
      real*8 r,intstarti
      real*8 desck,mixsn
      real*8 neckval
      real*8 usea,useb
      real*8 rhow,pi43
      real*8 rminb,radminr
      real*8 rminb4,radminr4
c
c
c     assign and initialize probe radius and constants
c
      rhow = 1.4d0
      usea = 0.0d0
      useb = 0.0d0
c
c     if atoms too far separated then no neck is formed
c
      if (r .gt. intstarti+desck+2.0d0*rhow) then
          neckval = 0.0d0
c
c     if atoms form a neck then calculate neck contribution
c
      else
         call neckcon (intstarti,desck,usea,useb)
         pi43 = 4.0d0 * third * pi 
         rminb = r - useb
         rminb4 = rminb**4
         radminr = intstarti + desck + 2.0d0*rhow - r
         radminr4 = radminr**4
         neckval = pi43 * mixsn * usea * rminb4 * radminr4
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine neckder  --  get neck descreening derivatives  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "neckder" returns the derivative of the neck descreening
c     integral for the Born chain rule term  
c
c
      subroutine neckder (r,intstarti,desck,mixsn,neckderi)
      use math
      use solute
      implicit none
      integer i
      real*8 r,intstarti
      real*8 desck,mixsn
      real*8 neckderi
      real*8 usea,useb
      real*8 rhow,pi43
      real*8 rminb,radminr
      real*8 rminb3,radminr3
      real*8 rminb4,radminr4

c
c     assign and initialize probe radius and constants
c
      rhow = 1.4d0
      usea = 0.0d0
      useb = 0.0d0
c
c     if atoms too far separated then no neck is formed
c
      if (r .gt. intstarti+desck+2.0d0*rhow) then
         neckderi = 0.0d0
c
c     if atoms form a neck then calculate neck contribution
c
      else
         call neckcon (intstarti,desck,usea,useb)
         pi43 = 4.0d0 * third * pi 
         rminb = r - useb
         rminb3 = rminb**3
         rminb4 = rminb3 * rminb
         radminr = intstarti + desck + 2.0d0*rhow - r
         radminr3 = radminr**3
         radminr4 = radminr3 * radminr
         neckderi = 4.0d0 * pi43 * (mixsn*usea*rminb3*radminr4
     &                            - mixsn*usea*rminb4*radminr3)
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine getbounds  --  get the radii array indices  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "getbounds" returns the indices in the radii array table    
c     to use for interpolation of Aij and Bij values  
c
c     variables and parameters:
c
c     rho     input radius
c     below   radsize index for radius next smaller than rho
c     above   radsize index for radius next larger than rho
c
c
      subroutine getbounds (rho,below,above)
      integer below,above
      integer numpoints
      real*8 rho
      real*8 calcindex
      real*8 minrad,maxrad
      real*8 space
c
c
      minrad = 0.80d0
      maxrad = 3.00d0
      space = 0.05d0
      numpoints = 45
      calcindex = 0.0d0
      below = 0
      above = 0
      calcindex = (rho-minrad) / space
      below = floor(calcindex) + 1
      above = below + 1
      if (above .ge. numpoints) then
         below = numpoints
         above = numpoints - 1
      else if (below .lt. 0) then
         below = 1
         above = 2
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine interp2d  --  interpolation of Aij/Bij values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "interp2d" returns the Aij and Bij values interpolated from
c     benchmark tables derived from Poisson-Boltzmann calculations
c
c     variables and parameters:
c
c     x1      radsize radius immediately smaller than descreened atom
c     x2      radsize radius immediately larger than descreened atom
c     y1      radsize radius immediately smaller than descreening atom
c     y2      radsize radius immediately larger than descreening atom
c     x       descreened atom radius + descreening offset, if used
c     y       descreening atom radius
c     fx1y1   constant for interacting atoms with radii x1 and y1
c     fx2y1   constant for interacting atoms with radii x2 and y1
c     fx1y2   constant for interacting atoms with radii x1 and y2
c     fx2y2   constant for interacting atoms with radii x2 and y2
c     val     returned interpolated constant value
c
c
      subroutine interp2d (x1,x2,y1,y2,x,y,fx1y1,fx2y1,fx1y2,fx2y2,val)
      real*8 x,x1,x2
      real*8 y,y1,y2
      real*8 fx1y1,fx1y2
      real*8 fx2y1,fx2y2
      real*8 fxy1,fxy2
      real*8 val
c
c
c     perform 2D interpolation of neck correction constant
c
      fxy1 = (x2-x)/(x2-x1)*fx1y1 + (x-x1)/(x2-x1)*fx2y1
      fxy2 = (x2-x)/(x2-x1)*fx1y2 + (x-x1)/(x2-x1)*fx2y2
      val = (y2-y)/(y2-y1)*fxy1 + (y-y1)/(y2-y1)*fxy2
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine neckcon  --  generate the neck constants  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "neckcon" returns the Aij and Bij values for a particular
c     pair of input radii
c
c     variables and parameters:
c
c     rhdsd    rho descreened - radius of descreened atom
c     rhdsg    rho descreening - radius of descreening atom
c     aijloc   returned value of Aij constant
c     bijloc   returned value of Bij constant
c
c
      subroutine neckcon (rhdsd,rhdsg,aijloc,bijloc)
      use solute
      integer lowi,highi
      integer lowj,highj
      real*8 rhdsd
      real*8 rhdsg
      real*8 aijloc
      real*8 bijloc
      real*8 rli,rhi
      real*8 rlj,rhj
      real*8 lla,hla
      real*8 lha,hha
      real*8 llb,hlb
      real*8 lhb,hhb
c
c
c     initialize some bounds and values
c
      lowi = 0
      lowj = 0
      highi = 0
      highj = 0
      aijloc = 0.0d0
      bijloc = 0.0d0
c
c     find Aij and Bij values via the rad array
c
      call getbounds (rhdsd,lowi,highi)
      call getbounds (rhdsg,lowj,highj)
      rli = radsize(lowi)
      rhi = radsize(highi)
      rlj = radsize(lowj)
      rhj = radsize(highj)
      lla = aij(lowi,lowj)
      hla = aij(highi,lowj)
      lha = aij(lowi,highj)
      hha = aij(highi,highj)
      llb = bij(lowi,lowj)
      hlb = bij(highi,lowj)
      lhb = bij(lowi,highj)
      hhb = bij(highi,highj)
      call interp2d (rli,rhi,rlj,rhj,rhdsd,rhdsg,lla,hla,lha,hha,aijloc)
      call interp2d (rli,rhi,rlj,rhj,rhdsd,rhdsg,llb,hlb,lhb,hhb,bijloc)
      if (aijloc .lt. 0.0d0) then
         aijloc = 0.0d0
      end if
      return
      end

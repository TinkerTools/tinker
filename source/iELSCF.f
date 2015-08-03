      module iELSCF
      use bath
      implicit none
      logical use_iELSCF
      logical first
      integer auxscale
      integer auxDoF
      real*8 omega
      real*8 aux_tautemp
      real*8 aux_kelvin
      real*8 aux_eksum
      real*8 aux_temp
      real*8 aux_ekin(3,3)
      real*8 vnh_aux(maxnose)
      real*8 qnh_aux(maxnose)
      real*8 gnh_aux(maxnose)
      real*8 pnh_aux(maxnose)
      real*8 pnh(maxnose)
      real*8, allocatable :: a_aux(:,:)
      real*8, allocatable :: v_aux(:,:)
      real*8, allocatable :: uind_aux(:,:)
      real*8, allocatable :: uinp_aux(:,:)
      character*120, auxstat
      save
      end

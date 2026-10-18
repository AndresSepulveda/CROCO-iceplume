! $Id: iceplume_const.h $
!
!======================================================================
! ICEPLUME for CROCO -- physical constants.
! Values from ROMS-ICEPLUME mod_iceplume.F (C. Wang). Only included
! in the plume physics routines of iceplume.F, which do NOT include
! scalars.h (that file also defines g and pi).
!======================================================================
!
      real alpha, tIce, sIce, rhoRef, rhoAir, g, cW, cI, L, tGade,
     &     lambda1, lambda2, lambda3, GamT, GamS, Cd, RiB, gRedBkg,
     &     CdBkg, velBkg, wIni, detSigma, avgAlpha, pi
      parameter (pi       = 3.14159265358979323846)
      parameter (alpha    = 0.1)        ! entrainment rate
      parameter (tIce     = -10.0)      ! ice temperature [degC]
      parameter (sIce     = 0.0)        ! ice salinity
      parameter (rhoRef   = 1020.0)     ! reference density
      parameter (rhoAir   = 1.225)      ! air density
      parameter (g        = 9.81)       ! gravity [m/s2]
      parameter (cW       = 3974.0)     ! heat capacity water
      parameter (cI       = 2000.0)     ! heat capacity ice
      parameter (L        = 335000.0)   ! latent heat of melting
      parameter (tGade    = -(L - cI*tIce)/cW)
      parameter (lambda1  = -0.0573)    ! freezing point slope
      parameter (lambda2  = 0.0832)     ! freezing point offset
      parameter (lambda3  = 0.000761)   ! freezing point depth slope
      parameter (GamT     = 0.022)      ! thermal transfer coeff.
      parameter (GamS     = 0.00062)    ! salt transfer coeff.
      parameter (Cd       = 0.065)      ! ice-plume drag
      parameter (RiB      = 1.0)        ! critical Richardson number
      parameter (gRedBkg  = 0.01)       ! background reduced gravity
      parameter (CdBkg    = 0.065)      ! background drag
      parameter (velBkg   = 0.03)       ! background velocity [m/s]
      parameter (wIni     = 0.1)        ! initial plume velocity
      parameter (detSigma = 0.5)        ! detrainment Gaussian sigma
      parameter (avgAlpha = 1.0)        ! exponential smoothing weight

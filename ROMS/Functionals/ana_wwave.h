      SUBROUTINE ana_wwave (ng, tile, model)
!
!! git $Id$
!! svn $Id: ana_wwave.h 1049 2020-11-30 04:34:51Z arango $
!!======================================================================
!! Copyright (c) 2002-2020 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets wind induced wave amplitude, direction and     !
!  period to be used in the bottom boundary layer formulation.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
#include "tile.h"
!
      CALL ana_wwave_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#if defined BBL_MODEL || defined NEARSHORE_MELLOR
     &                     FORCES(ng) % Dwave,                          &
#endif
#ifdef WAVES_HEIGHT
     &                     FORCES(ng) % Hwave,                          &
#endif
#ifdef WAVES_LENGTH
     &                     FORCES(ng) % Lwave,                          &
#endif
#ifdef WAVES_TOP_PERIOD
     &                     FORCES(ng) % Pwave_top,                      &
#endif
#ifdef WAVES_BOT_PERIOD
     &                     FORCES(ng) % Pwave_bot,                      &
#endif
#ifdef WAVES_UB
     &                     FORCES(ng) % Ub_swan,                        &
#endif
#ifdef TKE_WAVEDISS
     &                     FORCES(ng) % wave_dissip,                    &
#endif
     &                     GRID(ng) % angler,                           &
     &                     GRID(ng) % h)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(37)=MyFile
      END IF
!
      RETURN
      END SUBROUTINE ana_wwave
!
!***********************************************************************
      SUBROUTINE ana_wwave_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
#if defined BBL_MODEL || defined NEARSHORE_MELLOR
     &                           Dwave,                                 &
#endif
#ifdef WAVES_HEIGHT
     &                           Hwave,                                 &
#endif
#ifdef WAVES_LENGTH
     &                           Lwave,                                 &
#endif
#ifdef WAVES_TOP_PERIOD
     &                           Pwave_top,                             &
#endif
#ifdef WAVES_BOT_PERIOD
     &                           Pwave_bot,                             &
#endif
#ifdef WAVES_UB
     &                           Ub_swan,                               &
#endif
#ifdef TKE_WAVEDISS
     &                           wave_dissip,                           &
#endif
     &                           angler, h)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
# if defined BBL_MODEL || defined NEARSHORE_MELLOR
      real(r8), intent(inout) :: Dwave(LBi:,LBj:)
# endif
# ifdef WAVES_HEIGHT
      real(r8), intent(inout) :: Hwave(LBi:,LBj:)
# endif
# ifdef WAVES_LENGTH
      real(r8), intent(inout) :: Lwave(LBi:,LBj:)
# endif
# ifdef WAVES_TOP_PERIOD
      real(r8), intent(inout) :: Pwave_top(LBi:,LBj:)
# endif
# ifdef WAVES_BOT_PERIOD
      real(r8), intent(inout) :: Pwave_bot(LBi:,LBj:)
# endif
# ifdef WAVES_UB
      real(r8), intent(inout) :: Ub_swan(LBi:,LBj:)
# endif
# ifdef TKE_WAVEDISS
      real(r8), intent(inout) :: wave_dissip(LBi:,LBj:)
# endif

#else

      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
# if defined BBL_MODEL || defined NEARSHORE_MELLOR
      real(r8), intent(inout) :: Dwave(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_HEIGHT
      real(r8), intent(inout) :: Hwave(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_LENGTH
      real(r8), intent(inout) :: Lwave(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_TOP_PERIOD
      real(r8), intent(inout) :: Pwave_top(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_BOT_PERIOD
      real(r8), intent(inout) :: Pwave_bot(LBi:UBi,LBj:UBj)
# endif
# ifdef WAVES_UB
      real(r8), intent(inout) :: Ub_swan(LBi:UBi,LBj:UBj)
# endif
# ifdef TKE_WAVEDISS
      real(r8), intent(inout) :: wave_dissip(LBi:UBi,LBj:UBj)
# endif
#endif
!
!  Local variable declarations.
!
      integer :: i, j
!
      real(r8) :: cff, wdir
#if defined LAKE_SIGNELL
      real(r8) :: cff1, mxst, ramp_u, ramp_time, ramp_d
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set wind induced wave amplitude (m), direction (radians) and
!  period (s) at RHO-points.
!-----------------------------------------------------------------------
!
#if defined BL_TEST
      wdir=210.0_r8*deg2rad
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Hwave(i,j)=0.5_r8
          Dwave(i,j)=wdir
          Pwave_bot(i,j)=8.0_r8
        END DO
      END DO
#elif defined LAKE_SIGNELL
      mxst=0.25_r8         ! Wave amplitude (1/2 wave height) (meters)
      ramp_u=15.0_r8       ! start ramp UP at RAMP_UP (hours)
      ramp_time=10.0_r8    ! ramp from 0 to 1 over RAMP_TIME (hours)
      ramp_d=50.0_r8       ! start ramp DOWN at RAMP_DOWN (hours)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Dwave(i,j)=270.0_r8*deg2rad
          Pwave_bot(i,j)=5.0_r8    ! wave period (seconds)
           cff1=MIN((0.5_r8*(TANH((time(ng)/3600.0_r8-ramp_u)/          &
     &                            (ramp_time/5.0_r8))+1.0_r8)),         &
     &              (1.0_r8-(0.5_r8*(TANH((time(ng)/3600.0_r8-ramp_d)/  &
     &                                    (ramp_time/5.0_r8))+1.0_r8))))
          Hwave(i,j)=MAX((cff1*mxst),0.01_r8)
        END DO
      END DO
#elif defined NJ_BIGHT
!!    wdir=210.0_r8*deg2rad
      wdir=150.0_r8*deg2rad
      IF ((tdays(ng)-dstart).lt.1.5_r8) THEN
        cff=TANH(0.5_r8*(tdays(ng)-dstart))
        cff=1.0_r8
      ELSE
        cff=1.0_r8
      END IF
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Hwave(i,j)=0.12_r8
          Dwave(i,j)=wdir-angler(i,j)
          Pwave_bot(i,j)=10.0_r8
        END DO
      END DO
#elif defined SED_TOY
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          Hwave(i,j)=2.0_r8
          Dwave(i,j)=90.0_r8*deg2rad
          Pwave_bot(i,j)=8.0_r8
          Lwave(i,j)=20.0_r8
        END DO
      END DO
#else
      ana_wwave: no values provided for Hwave, Dwave, Pwave, Lwave.
#endif
!
!  Exchange boundary data.
!
#if defined WAVES_DIR
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Dwave)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Dwave)
# endif
#endif

#ifdef WAVES_HEIGHT
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Hwave)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Hwave)
# endif
#endif

#ifdef WAVES_LENGTH
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Lwave)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Lwave)
# endif
#endif

#ifdef WAVES_TOP_PERIOD
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Pwave_top)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Pwave_top)
# endif
#endif

#ifdef WAVES_BOT_PERIOD
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Pwave_bot)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Pwave_bot)
# endif
#endif

#ifdef WAVES_UB
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Ub_swan)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Ub_swan)
# endif
#endif

#ifdef TKE_WAVEDISS
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          wave_dissip)
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    wave_dissip)
# endif
#endif
!
      RETURN
      END SUBROUTINE ana_wwave_tile

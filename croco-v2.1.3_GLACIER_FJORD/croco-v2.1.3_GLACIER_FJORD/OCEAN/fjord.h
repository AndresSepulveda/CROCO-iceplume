!======================================================================
!  FJORD test case parameters (included in ana_grid.F, ana_initial.F
!  and analytical.F). Options are set in cppdefs.h:
!
!  FJORD_SILL      : sill at the mouth (undef: no sill)
!  FJORD_TIDES     : M2+O1 tide at the mouth (undef: no tide)
!  FJORD_EXCHANGE  : two-layer exchange at the mouth, zero net
!                    transport: fj_uin toward the head below fj_zex,
!                    compensating outflow above
!  FJORD_RIVER     : river at the head (undef: no river; ignored with
!                    FJORD_GLACIER)
!  FJORD_GLACIER   : marine-terminating glacier at the head
!                    (subglacial plume, ICEPLUME) instead of a river
!  FJORD_WIND_UP   : along-fjord wind toward the head (up-fjord)
!  FJORD_WIND_DOWN : along-fjord wind toward the mouth (down-fjord)
!  FJORD_GREENLAND : Greenland glacial fjord (Sermilik-like); default
!                    is a northern Patagonian fjord (Puyuhuapi-like)
!----------------------------------------------------------------------
!  fj_Hout   depth seaward of the sill (open boundary side)   [m]
!  fj_Hbas   depth of the inner fjord basin                    [m]
!  fj_Hsill  depth of the sill crest; also the reference depth of
!            the deep-water dye tpas02 (with or without sill)  [m]
!  fj_Hhead  depth at the fjord head                           [m]
!  fj_xsill  along-fjord position of the sill crest            [m]
!  fj_Lsill  e-folding half-width of the Gaussian sill         [m]
!  fj_Ltr    width of the outer-to-basin depth transition      [m]
!  fj_xhead  position of the head slope                        [m]
!  fj_Lhead  width of the head slope                           [m]
!  fj_f0     Coriolis parameter                              [1/s]
!  fj_aM2, fj_aO1   tidal amplitudes                          [m]
!  fj_pM2, fj_pO1   tidal phases                            [rad]
!  fj_TM2, fj_TO1   tidal periods                             [s]
!  fj_Tramp  tidal ramp time scale                             [s]
!  fj_tau0   along-fjord wind stress magnitude              [N/m2]
!             (FJORD_WIND_UP: toward the head, +x;
!              FJORD_WIND_DOWN: toward the mouth, -x)
!  fj_uin    exchange inflow speed below fj_zex (toward head) [m/s]
!  fj_zex    depth of the exchange interface at the mouth       [m]
!  fj_dex    thickness of the interface (tanh transition)       [m]
!  fj_wt0    wind onset time (model time)                      [s]
!  fj_wramp  wind ramp time scale                              [s]
!======================================================================
      real fj_Hout, fj_Hbas, fj_Hsill, fj_Hhead, fj_f0
      real fj_xsill, fj_Lsill, fj_Ltr, fj_xhead, fj_Lhead
      real fj_aM2, fj_pM2, fj_TM2, fj_aO1, fj_pO1, fj_TO1
      real fj_Tramp
#ifdef FJORD_GREENLAND
      parameter (fj_Hout=500., fj_Hbas=600., fj_Hsill=200.,
     &           fj_f0=1.328e-4)
#else
      parameter (fj_Hout=200., fj_Hbas=250., fj_Hsill=50.,
     &           fj_f0=-1.024e-4)
#endif
#ifdef FJORD_GLACIER
!  Vertical glacier front at the head: no shoaling (h=fj_Hbas)
      parameter (fj_Hhead=fj_Hbas)
#elif defined FJORD_GREENLAND
      parameter (fj_Hhead=100.)
#else
      parameter (fj_Hhead=30.)
#endif
      parameter (fj_xsill=12.e3, fj_Lsill=2.5e3, fj_Ltr=3.e3,
     &           fj_xhead=54.e3, fj_Lhead=3.e3)
      parameter (fj_aM2=0.82, fj_pM2=0., fj_TM2=12.4206012*3600.,
     &           fj_aO1=0.13, fj_pO1=0., fj_TO1=25.8193417*3600.,
     &           fj_Tramp=86400.)
      real fj_uin, fj_zex, fj_dex
      parameter (fj_uin=0.10, fj_zex=-0.5*fj_Hout, fj_dex=10.)
      real fj_tau0, fj_wt0, fj_wramp
      parameter (fj_tau0=0.1, fj_wt0=0., fj_wramp=86400.)

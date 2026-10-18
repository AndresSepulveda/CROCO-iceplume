! $Id: iceplume.h $
!
!======================================================================
! ICEPLUME for CROCO -- shared variables (COMMON blocks)
!
! Adapted from ROMS-ICEPLUME (mod_iceplume.F) by Chuning Wang,
! itself derived from IcePlume for MITgcm by Tom Cowton.
!
! Dimensions: Msrc (param.h) sources, N levels, NT tracers.
! Profiles on W points are (Msrc,0:N); on RHO points (Msrc,N).
!======================================================================
!
#ifdef ICEPLUME
!
! Inputs (read from croco.in, keyword "iceplume:")
!
! SG_typ   plume type: 0 none, 1 melt only, 2 half cone,
!          3 finite line, 4 sheet, 5 detached cone, 6 half ellipse
! SG_dep   grounding line depth [m, <0]; >=0 means bottom
! SG_len   length of subglacial channel/crack [m] (types 3, 4)
! SG_bar   subglacial discharge [m3/s]
! SG_trc   discharge tracer concentrations (T,S,passive)
! SG_I/Jrange  averaging window (global indices) for DET_AVERAGE
!
      integer SG_typ(Msrc)
      real SG_dep(Msrc), SG_len(Msrc), SG_bar(Msrc)
      real SG_trc(Msrc,NT)
      common /ipl_input_i/ SG_typ
      common /ipl_input_r/ SG_dep, SG_len, SG_bar, SG_trc
# ifdef ICEPLUME_DET_AVERAGE
      integer SG_Irange(Msrc,2), SG_Jrange(Msrc,2)
      common /ipl_input_avg/ SG_Irange, SG_Jrange
# endif
!
! Geometry: global indices of the wet plume cell (pl_Isrc,pl_Jsrc),
! plume direction pl_dir (+1/-1 across the face, 0 = mass source),
! cell lengths across (pl_dx) and along (pl_dy) the glacier front.
!
      integer pl_Isrc(Msrc), pl_Jsrc(Msrc)
      real pl_dir(Msrc), pl_dx(Msrc), pl_dy(Msrc)
      logical ipl_initialized
      common /ipl_geom_i/ pl_Isrc, pl_Jsrc
      common /ipl_geom_r/ pl_dir, pl_dx, pl_dy
      common /ipl_geom_l/ ipl_initialized
!
! Plume profiles on W levels
!
      real pl_zW(Msrc,0:N), pl_f(Msrc,0:N), pl_w(Msrc,0:N),
     &     pl_t(Msrc,0:N), pl_s(Msrc,0:N), pl_a(Msrc,0:N),
     &     pl_mInt(Msrc,0:N), pl_rho(Msrc,0:N),
     &     pl_lm(Msrc,0:N), pl_lc(Msrc,0:N)
      common /ipl_prof_w/ pl_zW, pl_f, pl_w, pl_t, pl_s, pl_a,
     &     pl_mInt, pl_rho, pl_lm, pl_lc
!
! Ambient profiles on RHO levels (tAm in situ, tpAm potential,
! rhoAm full in-situ density [kg/m3], level N+1 = air)
!
      real pl_zR(Msrc,N), pl_dz(Msrc,N), pl_sAm(Msrc,N),
     &     pl_tAm(Msrc,N), pl_tpAm(Msrc,N), pl_rhoAm(Msrc,N+1),
     &     pl_vAm(Msrc,N), pl_wAm(Msrc,N)
      common /ipl_prof_amb/ pl_zR, pl_dz, pl_sAm, pl_tAm, pl_tpAm,
     &     pl_rhoAm, pl_vAm, pl_wAm
!
! Fluxes per level [m3/s]: ent (<=0), det (>=0), plume melt m,
! background melt mB, background melt temperature tB; net transport
!
      integer pl_detI(Msrc,N)
      real pl_ent(Msrc,N), pl_det(Msrc,N), pl_detFrac(Msrc,N),
     &     pl_m(Msrc,N), pl_mB(Msrc,N), pl_tB(Msrc,N), pl_trs(Msrc)
      common /ipl_flux_i/ pl_detI
      common /ipl_flux_r/ pl_ent, pl_det, pl_detFrac, pl_m, pl_mB,
     &     pl_tB, pl_trs
# ifdef ICEPLUME_DET_NEUTRAL
      real pl_detF(Msrc,N), pl_detE(Msrc,N), pl_detTrc(Msrc,N,NT)
      common /ipl_flux_neu/ pl_detF, pl_detE, pl_detTrc
# endif
!
! Tracers: plume (trc), ambient (trcAm), background melt (trcB)
!
      real pl_trc(Msrc,NT), pl_trcAm(Msrc,N,NT), pl_trcB(Msrc,N,NT),
     &     pl_trcCum(Msrc,NT), pl_trcIni(Msrc,NT)
      common /ipl_trc/ pl_trc, pl_trcAm, pl_trcB, pl_trcCum,
     &     pl_trcIni
!
! Averaging of ambient density (ICEPLUME_DET_AVERAGE)
!
# ifdef ICEPLUME_DET_AVERAGE
      real pl_avgArea(Msrc,N), pl_rhoAm0(Msrc,N), pl_rhoSum(Msrc,N)
      common /ipl_avg/ pl_avgArea, pl_rhoAm0, pl_rhoSum
# endif
!
! Diagnostics: overshoot and neutral-buoyancy depths [m],
! integrated plume melt and background melt [m3/s]
!
      real pl_zOS(Msrc), pl_zNB(Msrc), pl_meltP(Msrc), pl_meltB(Msrc)
      common /ipl_diag/ pl_zOS, pl_zNB, pl_meltP, pl_meltB
!
#endif /* ICEPLUME */

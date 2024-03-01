#########################################################################
# This is a typical input file for OASIS3-MCT.
# Keywords used in previous versions of OASIS3 
# but now obsolete are marked "Not used"
# Don't hesitate to ask precisions or make suggestions (oasishelp@cerfacs.fr). 
#
# Any line beginning with # is ignored. Blank lines are not allowed.
#
#########################################################################
#
# NFIELDS: total number of fields being exchanged
 $NFIELDS
 3
#########################################################################
# NBMODEL: number of models and their names (6 characters) 
 $NBMODEL
 2  wrfexe   wwatch
###########################################################################
# RUNTIME: total simulated time for the actual run in seconds (<I8)
 $RUNTIME
 <runtime>
###########################################################################
# NLOGPRT: debug and time statistics informations printed in log file 
#          First number: 0 one log file for master, and one for other procs
#                        1 one log file for master, and one for other errors
#                        2 one file per proc with normal diagnostics
#                        5 as 2 + initial debug info
#                        10 as 5 + routine calling tree
#                        12 as 10 + some routine calling notes
#                        15 as 12 + even more debug diagnostics
#                        20 as 15 + some extra runtime analysis
#                        30 full debug information
#          Second number: time statistics
#          		 0 nothing calculated
#          		 1 one file for proc 0 and min/max of other procs
#          		 2 as 1 + one file per proc
#          		 3 as 2 + proc 0 writes all procs results in its file
 $NLOGPRT
 1 1
###########################################################################
# Beginning of fields exchange definition
 $STRINGS
#
# For each exchanged field:
#
# line 1: field in sending model, field in target model, unused, coupling 
#         period, number of transformation, restart file, field status
# line 2: nb of pts for sending model grid (without halo) first dim, and second dim,
#         for target grid first dim, and second dim, sending model grid name, target 
#         model grid name, lag = time step of sending model
# line 3: sending model grid periodical (P) or regional (R), and nb of overlapping 
#         points, target model grid periodical (P) or regional (R), and number of
#         overlapping points
# line 4: list of transformations performed
# line 5: parameters for each transformation
#
# See the correspondances between variables in models and in OASIS:
# Note: for CROCO and WRF  nesting capability is usable in coupled 
#       mode. For  CROCO, the domain  of  the  sent  variables  is 
#       defined  by the last number of coupled field name; for the
#       received  variable  the  CROCO  domain  is  defined by the
#       associated CPLMASK ( CROCO_VAR_CPLMASK* ),  and the domain
#       of  the  coupled model by the extension _EXT. For WRF, WRF
#       domain  is  defined  by  the  number  after WRF_d, and the
#       domain of the coupled model (CROCO for example)  is  given
#       by EXT_d in coupled field name 
#
# |--------------------------------------------------------------|
# | Possibly sent fields by CROCO:                 CROCO | OASIS |
# |--------------------------------------------------------------|
# |     t(:,:,N,nnew,itemp)  |    CROCO_SST                      |
# |                   zeta   |    CROCO_SSH                      |
# |     u,v (at rho points)  |    CROCO_UOCE CROCO_VOCE          |
# |     u,v (at rho points)                                      |
# | in E/N dir(rotated grids)|    CROCO_EOCE CROCO_NOCE          |
# |--------------------------------------------------------------|
# | Possibly received fields by CROCO:            CROCO | OASIS  |
# |--------------------------------------------------------------|
# |                  srflx   |    CROCO_SRFL                     |
# |       stflx(:,:,isalt)   |    CROCO_EVPR                     |
# |      stflx(,:,:,itemp)   |    CROCO_STFL                     |
# |                  sustr   |    CROCO_UTAU  CROCO_ETAU         |
# |                  svstr   |    CROCO_VTAU  CROCO_NTAU         |
# |                  smstr   |    CROCO_TAUM                     |
# |                patm2d    |    CROCO_PSFC                     |
# |                  whrm    |    CROCO_HS                       |
# |                  wfrq    |    CROCO_T0M1                     |
# |             wdrx,wdre    |    CROCO_DIR                      |
# |                   foc    |    CROCO_FOC                      |
# |                   wlm    |    CROCO_LM                       |
# |                   ubr    |    CROCO_UBRX  CROCO_UBRY         |
# |               ust_ext    |    CROCO_USSX  CROCO_USSY         |
# |             twox,twoy    |    CROCO_UTWO  CROCO_VTWO         |
# |             etwo,ntwo    |    CROCO_ETWO  CROCO_NTWO         |
# |             tawx,tawy    |    CROCO_UTAW  CROCO_VTAW         |
# |             etaw,ntaw    |    CROCO_ETAW  CROCO_NTAW         |
# |--------------------------------------------------------------|
# | Possibly sent fields by WW3:                    WW3 | OASIS  |
# |--------------------------------------------------------------|
# |            not defined   |    WW3_ODRY                       |
# |                   T0M1   |    WW3_T0M1                       |
# |                     HS   |    WW3__OHS WW3__AHS              |
# |                    THM   |    WW3__DIR                       |
# |                    TAW   |    WW3_TAWX WW3_TAWY              |
# |                    TWO   |    WW3_TWOX WW3_TWOY              |
# |                    BHD   |    WW3__BHD                       |
# |                    UBR   |    WW3__UBR                       |
# |                    FOC   |    WW3__FOC                       |
# |                     LM   |    WW3___LM                       |
# |                    CUR   |    WW3_WSSU WW3_WSSV              |
# |                   ACHA   |    WW3_ACHA                       |
# |                     FP   |    WW3___FP                       |
# |--------------------------------------------------------------|
# | Possibly received fields by WW3:                WW3 | OASIS  |
# |--------------------------------------------------------------|
# |            not defined   |    WW3_OWDH WW3_OWDU WW3_OWDV     |
# |                    SSH   |    WW3__SSH                       |
# |                    CUR   |    WW3_OSSU WW3_OSSV              |
# |                    WND   |    WW3__U10 WW3__V10              |
# |--------------------------------------------------------------|
# | Possibly sent fields by WRF:                    WRF | OASIS  |
# |--------------------------------------------------------------|
# |  QFX-(RAINCV+RAINNCV+                                        |
# |   SNOWNCV+HAILNCV+                                           |
# |   GRAUPELNCV)/DT     |    WRF_d01_EXT_d01_EVAP-PRECIP        |
# |  (RAINCV+RAINNCV)/DT |    WRF_d01_EXT_d01_LIQUID_PRECIP      |
# |  (SNOWNCV+HAILNCV+                                           |
# |       GRAUPELNCV)/DT |    WRF_d01_EXT_d01_SOLID_PRECIP       |
# |                QFX   |    WRF_d01_EXT_d01_TOTAL_EVAP         |
# |                GSW   |    WRF_d01_EXT_d01_SURF_NET_SOLAR     |
# |   GLW-STBOLT*EMISS                                           |
# |     *SST**4          |    WRF_d01_EXT_d01_SURF_NET_LONGWAVE  |
# |                -LH   |    WRF_d01_EXT_d01_SURF_LATENT        |
# |               -HFX   |    WRF_d01_EXT_d01_SURF_SENSIBLE      |
# |   GLW-STBOLT*EMISS                                           |
# |     *SST**4-LH-HFX   |    WRF_d01_EXT_d01_SURF_NET_NON-SOLAR |
# | taut * u_uo / wspd   |    WRF_d01_EXT_d01_TAUX               |
# | taut * u_uo / wspd   |    WRF_d01_EXT_d01_TAUY               |
# |               taut   |    WRF_d01_EXT_d01_TAUMOD             |
# | cosa*taui - sina*tauj|    WRF_d01_EXT_d01_TAUE               |
# | cosa*tauj - sina*taui|    WRF_d01_EXT_d01_TAUN               |
# |               u_uo   |    WRF_d01_EXT_d01_WND_U_01           |
# |               v_vo   |    WRF_d01_EXT_d01_WND_V_01           |
# | cosa*u_uo - sina*v_vo|    WRF_d01_EXT_d01_WND_E_01           |
# | cosa*v_vo - sina*u_uo|    WRF_d01_EXT_d01_WND_N_01           |
# |               psfc   |    WRF_d01_EXT_d01_PSFC               |
# |--------------------------------------------------------------|
# | Possibly received fields by WRF:                WRF | OASIS  |
# |--------------------------------------------------------------|
# |                    SST   |    WRF_d01_EXT_d01_SST            |
# |                   UOCE   |    WRF_d01_EXT_d01_UOCE           |
# |                   VOCE   |    WRF_d01_EXT_d01_VOCE           |
# |                   UOCE   |    WRF_d01_EXT_d01_EOCE           |
# |                   VOCE   |    WRF_d01_EXT_d01_NOCE           |
# |               CHA_COEF   |    WRF_d01_EXT_d01_CHA_COEF       |
# |--------------------------------------------------------------|
#
#                     ------------------------------------
#                        WRF (wrfexe) ==> WW3 (wwatch)
#                     ------------------------------------
#~~~~~~~~~~~
# U_01 : wind U component at first level (m/s)
#~~~~~~~~~~~
WRF_d01_EXT_d01_WND_E_01 WW3__U10 1 <cpldt> 2 atm.nc  EXPORTED
<atmnx> <atmny> <wavnx> <wavny> atmt ww3t LAG=<atmdt>
R 0 R 0
LOCTRANS SCRIPR
AVERAGES
DISTWGT LR SCALAR LATLON 1 4
#
#~~~~~~~~~~~
# V_01 : wind V component at first level (m/s)
#~~~~~~~~~~~
WRF_d01_EXT_d01_WND_N_01 WW3__V10 1 <cpldt> 2 atm.nc  EXPORTED
<atmnx> <atmny> <wavnx> <wavny> atmt ww3t LAG=<atmdt>
R 0 R 0
LOCTRANS SCRIPR
AVERAGES
DISTWGT LR SCALAR LATLON 1 4
#
#                     ------------------------------------
#                        WW3 (wwatch) ==> WRF (wrfexe)
#                     ------------------------------------
#
#~~~~~~~~~~~
# CHA_COEF : Charnock coefficient
#~~~~~~~~~~~
WW3_ACHA WRF_d01_EXT_d01_CHA_COEF 1 <cpldt> 2 wav.nc  EXPORTED
<wavnx> <wavny> <atmnx> <atmny> ww3t atmt LAG=<wavdt>
R  0  R  0
LOCTRANS SCRIPR
AVERAGES
DISTWGT LR SCALAR LATLON 1 4
#
###########################################################################
$END

#!/bin/bash
#
# git $Id$
# svn $Id: job_i4dvar.sh 1019 2020-04-30 20:46:51Z arango $
#######################################################################
# Copyright (c) 2002-2020 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
# Incremental strong constraint I4D-Var job BASH script:              #
#                                                                     #
# This script NEEDS to be run before any run:                         #
#                                                                     #
#   (1) It copies a new clean nonlinear model initial conditions      #
#       file. The nonlinear model is initialized from the             #
#       background or reference state.                                #
#   (2) Specify initial conditions, boundary conditions, and surface  #
#       forcing error convariance input standard deviations files.    #
#   (3) Specify initial conditions, boundary conditions, and surface  #
#       forcing error convariance input/output normalization factors  #
#       files.                                                        #
#   (4) Copy a clean copy of the observations NetCDF file.            #
#   (5) Create 4D-Var input script "i4dvar.in" from template and      #
#       specify the error covariance standard deviation, error        #
#       covariance normalization factors, and observation files to    #
#       be used.                                                      #
#                                                                     #
#######################################################################

# Set path definition to one directory up in the tree.

 Dir=`dirname ${PWD}`

# Set string manipulations perl script.

 SUBSTITUTE=${ROMS_ROOT}/ROMS/Bin/substitute

# Copy nonlinear model initial conditions file.

 cp -p ${Dir}/Data/wc13_ini.nc wc13_ini.nc

# Set initial conditions, boundary conditions and surface forcing
# error covariance standard deviations files.

 STDnameI=${Dir}/Data/wc13_std_i.nc
 STDnameB=${Dir}/Data/wc13_std_b.nc
 STDnameF=${Dir}/Data/wc13_std_f.nc

# Set initial conditions, boundary conditions and surface forcing
# error covariance normalization factors files.

 NRMnameI=${Dir}/Data/wc13_nrm_i.nc
 NRMnameB=${Dir}/Data/wc13_nrm_b.nc
 NRMnameF=${Dir}/Data/wc13_nrm_f.nc

# Set observations file.

 OBSname=wc13_obs.nc

# Get a clean copy of the observation file.  This is really
# important since this file is modified.

 cp -p ${Dir}/Data/${OBSname} .

# Modify 4D-Var template input script and specify above files.

 I4DVAR=i4dvar.in
 if [ -f $I4DVAR ]; then
   /bin/rm $I4DVAR
 fi
 cp s4dvar.in $I4DVAR

 $SUBSTITUTE $I4DVAR roms_std_i.nc $STDnameI
 $SUBSTITUTE $I4DVAR roms_std_b.nc $STDnameB
 $SUBSTITUTE $I4DVAR roms_std_f.nc $STDnameF
 $SUBSTITUTE $I4DVAR roms_nrm_i.nc $NRMnameI
 $SUBSTITUTE $I4DVAR roms_nrm_b.nc $NRMnameB
 $SUBSTITUTE $I4DVAR roms_nrm_f.nc $NRMnameF
 $SUBSTITUTE $I4DVAR roms_obs.nc $OBSname
 $SUBSTITUTE $I4DVAR roms_hss.nc wc13_hss.nc
 $SUBSTITUTE $I4DVAR roms_lcz.nc wc13_lcz.nc
 $SUBSTITUTE $I4DVAR roms_mod.nc wc13_mod.nc
 $SUBSTITUTE $I4DVAR roms_err.nc wc13_err.nc

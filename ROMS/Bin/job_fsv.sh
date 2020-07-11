#!/bin/bash
#
# svn $Id$
#######################################################################
# Copyright (c) 2002-2020 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
#  Generalized Stability Theory: Forcing Singular Vectors             #
#                                                                     #
#  This script is used to set-up ROMS Forcing Singular Vectors        #
#  algorithm.                                                         #
#                                                                     #
#######################################################################

# Set ROOT of the directory to run application.  The following
# "dirname" command returns a path by removing any suffix from
# the last slash ('/').  It returns a path above current diretory.

Dir=`dirname ${PWD}`

# Set basic state trajectory, forward file:

#HISname=${Dir}/Forward/gyre3d_his_00.nc
 HISname=${Dir}/Forward/gyre3d_his_01.nc

FWDname=gyre3d_fwd.nc

if [ -f $FWDname ]; then
  /bin/rm $FWDname
fi
ln -s -v $HISname $FWDname

# Set tangent linear model initial conditions file: zero fields.

ITLname=gyre3d_itl.nc
if [ -f $ITLname ]; then
  /bin/rm $ITLname
fi
ln -s -v ${Dir}/Data/gyre3d_ini_zero.nc $ITLname

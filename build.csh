#!/bin/csh -f

source /discover/nobackup/bzhao/geos_mom6/test-land-bcs/GEOSgcm/install-datm16/bin/g5_modules



ifort -o surf_layer surf_layer.F90 icepack_kinds.F90 icepack_parameters.F90


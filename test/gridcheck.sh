#!/bin/bash

vpid="$OMPI_MCA_ns_nds_vpid"
if test "$1" = "gdb" ; then
  shift 1
  xterm -e gdb --eval-command=run --args ./gridcheck "$@"
else
  ./gridcheck "$@" &> gridcheck.$vpid.out
fi

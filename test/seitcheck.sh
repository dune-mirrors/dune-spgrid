#!/bin/bash

base=`dirname $0`

vpid="$OMPI_MCA_ns_nds_vpid"
size="$OMPI_MCA_ns_nds_num_procs"
if test "$1" = "gdb" ; then
  shift 1
  xterm -e gdb --eval-command=run --args ./seitcheck "$@"
else
  mkdir -p $base/mpiout
  ./seitcheck "$@" &> $base/mpiout/seitcheck-$size-$vpid
fi

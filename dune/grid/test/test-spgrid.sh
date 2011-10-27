#!/bin/bash

base=`dirname $0`

vpid="$OMPI_COMM_WORLD_RANK"
size="$OMPI_COMM_WORLD_SIZE"

if test -z "$size" ; then
  echo "Unable to determine number of processors."
  set
  exit 1;
fi

if test -z "$vpid" ; then
  echo "Unable to determine process rank."
  exit 1;
fi

if test "$1" = "gdb" ; then
  shift 1
  xterm -e gdb --eval-command=run --args ./test-spgrid "$@"
else
  mkdir -p $base/mpiout
  ./test-spgrid "$@" &> $base/mpiout/gridcheck-$size-$vpid
fi

#!/bin/bash

vpid="$OMPI_MCA_ns_nds_vpid"
./gridcheck "$@" &> gridcheck.$vpid.out

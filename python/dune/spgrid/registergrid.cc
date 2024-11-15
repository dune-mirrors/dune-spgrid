#include <config.h>
#define INCLUDE_REGALUGRID_INLINE
#include "registergrid.hh"
#ifndef DIM
#error "DIM need to be defined!"
#endif
template void registerSPGrid<DIM>(pybind11::module);

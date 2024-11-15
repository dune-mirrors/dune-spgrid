#include <config.h>

#include <dune/python/pybind11/pybind11.h>

#ifdef DUNE_ENABLE_PYTHONMODULE_PRECOMPILE
#include "registergrid.hh"
#endif

PYBIND11_MODULE( _spgrid, module )
{
#ifdef DUNE_ENABLE_PYTHONMODULE_PRECOMPILE
  // pre-compiled objects for SPGrid
  registerSPGrid< 1 > ( module );
  registerSPGrid< 2 > ( module );
  registerSPGrid< 3 > ( module );
  registerSPGrid< 4 > ( module );
#endif
}

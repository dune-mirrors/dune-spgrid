#ifndef DUNE_SPGRIDTYPE_HH
#define DUNE_SPGRIDTYPE_HH

#include <dune/grid/utility/griddim.hh>

#if defined SPGRID
  #if HAVE_GRIDTYPE
    #error "Ambiguous definition of GRIDTYPE."
  #endif

  #include <dune/grid/spgrid.hh>
  namespace Dune
  {
    namespace GridSelector
    {
      typedef SPGrid< double, dimgrid > GridType;
    }
  }
  #define HAVE_GRIDTYPE 1
#endif

#endif // #ifndef DUNE_SPGRIDTYPE_HH

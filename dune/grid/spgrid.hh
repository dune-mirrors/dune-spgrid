#ifndef DUNE_SPGRID_HH
#define DUNE_SPGRID_HH

#include <dune/grid/spgrid/grid.hh>
#include <dune/grid/spgrid/hierarchicsearch.hh>
#include <dune/grid/spgrid/persistentcontainer.hh>

#if HAVE_DUNE_FEM 
#include <dune/fem/quadrature/caching/twistutility.hh>

namespace Dune { 

  // Specialization for SPGrid
  //
  template< class ct, int dim, SPRefinementStrategy strategy , class Comm > 
  struct TwistUtility< SPGrid< ct, dim, strategy, Comm > > 
    : public TwistFreeTwistUtility< SPGrid< ct, dim, strategy, Comm > > 
  {
  };
}

#endif // #if HAVE_DUNE_FEM 

#endif // #ifndef DUNE_SPGRID_HH

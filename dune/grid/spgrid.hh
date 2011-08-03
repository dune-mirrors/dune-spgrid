#ifndef DUNE_SPGRID_HH
#define DUNE_SPGRID_HH

#include <dune/grid/spgrid/grid.hh>
#include <dune/grid/spgrid/hierarchicsearch.hh>
#include <dune/grid/spgrid/persistentcontainer.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class Grid >
  struct TwistFreeTwistUtility;

  template< class Grid >
  struct TwistUtility;



  // TwistUtility for SPGrid
  // -----------------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm > 
  struct TwistUtility< SPGrid< ct, dim, strategy, Comm > > 
  : public TwistFreeTwistUtility< SPGrid< ct, dim, strategy, Comm > >
  {};

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_HH

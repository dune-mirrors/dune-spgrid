#ifndef DUNE_SPGRID_HH
#define DUNE_SPGRID_HH

#include <dune/grid/spgrid/backuprestore.hh>
#include <dune/grid/spgrid/grid.hh>
#include <dune/grid/spgrid/hierarchicsearch.hh>
#include <dune/grid/spgrid/persistentcontainer.hh>
#include <dune/grid/spgrid/tree.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< class Grid >
    struct TwistFreeTwistUtility;

    template< class Grid >
    struct TwistUtility;



    // TwistUtility for SPGrid
    // -----------------------

    template< class ct, int dim, template< int > class Ref, class Comm >
    struct TwistUtility< SPGrid< ct, dim, Ref, Comm > >
      : public TwistFreeTwistUtility< SPGrid< ct, dim, Ref, Comm > >
    {};

  } // end namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_HH

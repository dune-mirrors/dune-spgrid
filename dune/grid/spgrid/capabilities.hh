#ifndef DUNE_SPGRID_CAPABILITIES_HH
#define DUNE_SPGRID_CAPABILITIES_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/grid/extensions/superentityiterator.hh>

#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ct, int dim, SPRefinementStrategy strategy >
  class SPGrid;



  // Capabilities
  // ------------

  namespace Capabilities
  {

    template< class ct, int dim, SPRefinementStrategy strategy, int codim >
    struct hasEntity< SPGrid< ct, dim, strategy >, codim >
    {
      static const bool v = true;
    };

    template< class ct, int dim, SPRefinementStrategy strategy >
    struct isParallel< SPGrid< ct, dim, strategy > >
    {
      static const bool v = false;
    };

    template< class ct, int dim, SPRefinementStrategy strategy >
    struct isLevelwiseConforming< SPGrid< ct, dim, strategy > >
    {
      static const bool v = true;
    };

    template< class ct, int dim, SPRefinementStrategy strategy >
    struct isLeafwiseConforming< SPGrid< ct, dim, strategy > >
    {
      static const bool v = true;
    };

    template< class ct, int dim, SPRefinementStrategy strategy >
    struct hasBackupRestoreFacilities< SPGrid< ct, dim, strategy > >
    {
      static const bool v = true;
    };

    template< class ct, int dim, SPRefinementStrategy strategy >
    struct IsUnstructured< SPGrid< ct, dim, strategy > >
    {
      static const bool v = false;
    };

    template< class ct, int dim, SPRefinementStrategy strategy >
    struct threadSafe< SPGrid< ct, dim, strategy > >
    {
      static const bool v = false;
    };

    template< class ct, int dim, SPRefinementStrategy strategy >
    struct viewThreadSafe< SPGrid< ct, dim, strategy > >
    {
      static const bool v = false;
    };



    // non-standard capabilities (see dune-fem)
    // ----------------------------------------

    template< class Grid >
    struct hasHierarchicIndexSet;

    template< class ct, int dim, SPRefinementStrategy strategy >
    struct hasHierarchicIndexSet< SPGrid< ct, dim, strategy > >
    {
      static const bool v = true;
    };


    template< class Grid >
    struct supportsCallbackAdaptation;

    template< class ct, int dim, SPRefinementStrategy strategy >
    struct supportsCallbackAdaptation< SPGrid< ct, dim, strategy > >
    {
      static const bool v = true;
    };

  }



  // Extensions
  // ----------

  namespace Extensions
  {

    template< class ct, int dim, SPRefinementStrategy strategy, int codim >
    struct SuperEntityIterator< SPGrid< ct, dim, strategy >, codim >
    {
      static const bool v = true;
    };

  }

}

#endif // #ifndef DUNE_SPGRID_CAPABILITIES_HH

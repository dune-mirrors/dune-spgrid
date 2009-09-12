#ifndef DUNE_SPGRID_CAPABILITIES_HH
#define DUNE_SPGRID_CAPABILITIES_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/grid/extensions/superentityiterator.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ct, int dim >
  class SPGrid;



  // Capabilities
  // ------------

  namespace Capabilities
  {

    template< class ct, int dim, int codim >
    struct hasEntity< SPGrid< ct, dim >, codim >
    {
      static const bool v = true;
    };

    template< class ct, int dim >
    struct isParallel< SPGrid< ct, dim > >
    {
      static const bool v = false;
    };

    template< class ct, int dim >
    struct isLevelwiseConforming< SPGrid< ct, dim > >
    {
      static const bool v = true;
    };

    template< class ct, int dim >
    struct isLeafwiseConforming< SPGrid< ct, dim > >
    {
      static const bool v = true;
    };

    template< class ct, int dim >
    struct hasBackupRestoreFacilities< SPGrid< ct, dim > >
    {
      static const bool v = true;
    };

    template< class ct, int dim >
    struct IsUnstructured< SPGrid< ct, dim > >
    {
      static const bool v = false;
    };

    template< class ct, int dim >
    struct threadSafe< SPGrid< ct, dim > >
    {
      static const bool v = false;
    };

    template< class ct, int dim >
    struct viewThreadSafe< SPGrid< ct, dim > >
    {
      static const bool v = false;
    };



    // non-standard capabilities (see dune-fem)
    // ----------------------------------------

    template< class Grid >
    struct hasHierarchicIndexSet;

    template< class ct, int dim >
    struct hasHierarchicIndexSet< SPGrid< ct, dim > >
    {
      static const bool v = true;
    };


    template< class Grid >
    struct supportsCallbackAdaptation;

    template< class ct, int dim >
    struct supportsCallbackAdaptation< SPGrid< ct, dim > >
    {
      static const bool v = true;
    };

  }



  // Extensions
  // ----------

  namespace Extensions
  {

    template< class ct, int dim, int codim >
    struct SuperEntityIterator< SPGrid< ct, dim >, codim >
    {
      static const bool v = true;
    };

  }

}

#endif // #ifndef DUNE_SPGRID_CAPABILITIES_HH

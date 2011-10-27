#ifndef DUNE_CARTESIANGRID_CAPABILITIES_HH
#define DUNE_CARTESIANGRID_CAPABILITIES_HH

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid >
  class CartesianGrid;



  // Capabilities
  // ------------

  namespace Capabilities
  {

    // Capabilities from dune-grid
    // ---------------------------

    template< class HostGrid >
    struct hasSingleGeometryType< CartesianGrid< HostGrid > >
    {
      /** \brief all elements in \ref Dune::CartesianGrid "CartesianGrid" have the same geometry type */
      static const bool v = true;
      /** \brief \ref Dune::CartesianGrid "CartesianGrid" has only cube elements */
      static const unsigned int topologyId = GenericGeometry::CubeTopology< HostGrid::dimension >::type::id;
    };

    /** \brief Is the grid Cartesian?
     *
     *  Cartesian grids satisfy the following properties:
     *  - all geometries are affine
     *  - The unit outer normal can be computed by the following code:
     *  \code
     *  FieldVector< ctype, dim > n( 0 );
     *  n[ face / 2 ] = ctype( 2*(face % 2) - 1 );
     *  \endcode
     *  .
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class HostGrid > 
    struct isCartesian< CartesianGrid< HostGrid > >
    {
      /** \brief \ref Dune::CartesianGrid "CartesianGrid" is a Cartesian grid */
      static const bool v = true;
    };


    template< class HostGrid, int codim >
    struct hasEntity< CartesianGrid< HostGrid >, codim >
    {
      static const bool v = hasEntity< HostGrid, codim >::v;
    };

    
    template< class HostGrid >
    struct isParallel< CartesianGrid< HostGrid > >
    {
      static const bool v = isParallel< HostGrid >::v;
    };


    template< class HostGrid, int codim >
    struct canCommunicate< CartesianGrid< HostGrid >, codim >
    {
      static const bool v = canCommunicate< HostGrid, codim >::v;
    };


    template< class HostGrid >
    struct hasBackupRestoreFacilities< CartesianGrid< HostGrid > >
    {
      static const bool v = hasBackupRestoreFacilities< HostGrid >::v;
    };

    template< class HostGrid >
    struct isLevelwiseConforming< CartesianGrid< HostGrid > >
    {
      static const bool v = isLevelwiseConforming< HostGrid >::v;
    };

    template< class HostGrid >
    struct isLeafwiseConforming< CartesianGrid< HostGrid > >
    {
      static const bool v = isLeafwiseConforming< HostGrid >::v;
    };

    template< class HostGrid >
    struct threadSafe< CartesianGrid< HostGrid > >
    {
      static const bool v = false;
    };

    template< class HostGrid >
    struct viewThreadSafe< CartesianGrid< HostGrid > >
    {
      static const bool v = false;
    };



    // non-standard capabilities
    // -------------------------

    template< class Grid >
    struct hasHierarchicIndexSet;

    template< class HostGrid >
    struct hasHierarchicIndexSet< CartesianGrid< HostGrid > >
    {
      static const bool v = false;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_CARTESIANGRID_CAPABILITIES_HH

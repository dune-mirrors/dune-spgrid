#ifndef DUNE_SPGRID_CAPABILITIES_HH
#define DUNE_SPGRID_CAPABILITIES_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>

#include <dune/grid/extensions/superentityiterator.hh>

#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  class SPGrid;



  // Capabilities
  // ------------

  /** \brief namespace containing all capability */
  namespace Capabilities
  {

    /** \brief Do elements of a grid always have the same geometry type?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct hasSingleGeometryType< SPGrid< ct, dim, strategy, Comm > >
    {
      /** \brief all elements in \ref Dune::SPGrid "SPGrid" have the same
       *         geometry type */
      static const bool v = true;
      /** \brief \ref Dune::SPGrid "SPGrid" has only cube elements */
      static const unsigned int topologyId = GenericGeometry::CubeTopology< dim >::type::id;
    };

    /** \brief Does a grid implement entities of a codimension?
     *
     *  \tparam  Grid   grid for which the information is desired
     *  \tparam  codim  codimension in question
     */
    template< class ct, int dim, SPRefinementStrategy strategy, class Comm, int codim >
    struct hasEntity< SPGrid< ct, dim, strategy, Comm >, codim >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" implements entities for all
       *         codimensions */
      static const bool v = ((codim >= 0) && (codim <= dim));
    };

#if HAVE_MPI
    /** \brief Does a grid support parallel programs?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, SPRefinementStrategy strategy >
    struct isParallel< SPGrid< ct, dim, strategy, MPI_Comm > >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" with MPI_Comm supports
       *         parallelism */
      static const bool v = true;
    };

    /** \brief Can a parallel grid communicate on a given codimension?
     *
     *  \tparam  Grid   grid for which the information is desired
     *  \tparam  codim  codimension in question
     *
     *  \note In order to communicate on a given codimension, the grid has to
     *        implement entities for that codimension.
     */
    template< class ct, int dim, SPRefinementStrategy strategy, int codim >
    struct canCommunicate< SPGrid< ct, dim, strategy, MPI_Comm >, codim >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" with MPI_Comm can communicate on
       *         all codimensions */
      static const bool v = ((codim >= 0) && (codim <= dim));
    };
#endif // #if HAVE_MPI

    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct isLevelwiseConforming< SPGrid< ct, dim, strategy, Comm > >
    {
      static const bool v = true;
    };

    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct isLeafwiseConforming< SPGrid< ct, dim, strategy, Comm > >
    {
      static const bool v = true;
    };

    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct hasBackupRestoreFacilities< SPGrid< ct, dim, strategy, Comm > >
    {
      static const bool v = true;
    };

    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct threadSafe< SPGrid< ct, dim, strategy, Comm > >
    {
      static const bool v = false;
    };

    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct viewThreadSafe< SPGrid< ct, dim, strategy, Comm > >
    {
      static const bool v = false;
    };



    // non-standard capabilities (see dune-fem)
    // ----------------------------------------

    template< class Grid >
    struct hasHierarchicIndexSet;

    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct hasHierarchicIndexSet< SPGrid< ct, dim, strategy, Comm > >
    {
      static const bool v = true;
    };


    template< class Grid >
    struct supportsCallbackAdaptation;

    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct supportsCallbackAdaptation< SPGrid< ct, dim, strategy, Comm > >
    {
      static const bool v = true;
    };

  }



  // Extensions
  // ----------

  namespace Extensions
  {

    template< class ct, int dim, SPRefinementStrategy strategy, class Comm, int codim >
    struct SuperEntityIterator< SPGrid< ct, dim, strategy, Comm >, codim >
    {
      static const bool v = true;
    };

  }

}

#endif // #ifndef DUNE_SPGRID_CAPABILITIES_HH

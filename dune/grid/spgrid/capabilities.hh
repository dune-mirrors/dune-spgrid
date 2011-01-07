#ifndef DUNE_SPGRID_CAPABILITIES_HH
#define DUNE_SPGRID_CAPABILITIES_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/grid/genericgeometry/topologytypes.hh>

#include <dune/grid/extensions/superentityiterator.hh>

#include <dune/grid/spgrid/refinement.hh>

/** \file
 *  \author Martin Nolte
 *  \brief  capabilities for \ref Dune::SPGrid "SPGrid"
 */

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  class SPGrid;



  // Capabilities
  // ------------

  /** \brief namespace containing all capabilities */
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

    /** \brief Are all levels of a grid always conform?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct isLevelwiseConforming< SPGrid< ct, dim, strategy, Comm > >
    {
      /** \brief All levels of a \ref Dune::SPGrid "SPGrid" are always conform */
      static const bool v = true;
    };

    /** \brief Is the leaf level of a grid always conform?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct isLeafwiseConforming< SPGrid< ct, dim, strategy, Comm > >
    {
      /** \brief The leaf level of a \ref Dune::SPGrid "SPGrid" are always conform */
      static const bool v = true;
    };

    /** \brief Does a grid provide backup and restore facilities?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct hasBackupRestoreFacilities< SPGrid< ct, dim, strategy, Comm > >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" provides backup and restore facilities */
      static const bool v = true;
    };

    /** \brief Is a grid implementation thread safe?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct threadSafe< SPGrid< ct, dim, strategy, Comm > >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" is not thread safe */
      static const bool v = false;
    };

    /** \brief Is a grid implementation thread safe while not being modified?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct viewThreadSafe< SPGrid< ct, dim, strategy, Comm > >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" is not thread safe */
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

    /** Does a grid implementation support callback adaptation?
     *
     *  \tparam  Grid  grid for which the information is desired
     *
     *  \note This is not a standard dune-grid capability.
     */
    template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
    struct supportsCallbackAdaptation< SPGrid< ct, dim, strategy, Comm > >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" supports callback adaptation */
      static const bool v = true;
    };

  }



  // Extensions
  // ----------

  namespace Extensions
  {

    /** \brief Does a grid support superentity iterators of a codimension?
     *
     *  \tparam  Grid   grid for which the information is desired
     *  \tparam  codim  codimension in question
     */
    template< class ct, int dim, SPRefinementStrategy strategy, class Comm, int codim >
    struct SuperEntityIterator< SPGrid< ct, dim, strategy, Comm >, codim >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" supports superentity iterators for all
       *         codimensions */
      static const bool v = ((codim >= 0) && (codim <= dim));
    };

  }

}

#endif // #ifndef DUNE_SPGRID_CAPABILITIES_HH

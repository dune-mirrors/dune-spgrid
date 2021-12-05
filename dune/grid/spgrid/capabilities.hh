#ifndef DUNE_SPGRID_CAPABILITIES_HH
#define DUNE_SPGRID_CAPABILITIES_HH

#if HAVE_MPI
#include <mpi.h>
#endif

#include <dune/geometry/type.hh>

#include <dune/grid/common/capabilities.hh>

#include <dune/grid/extensions/superentityiterator.hh>

#include <dune/grid/spgrid/declaration.hh>

/** \file
 *  \author Martin Nolte
 *  \brief  capabilities for \ref Dune::SPGrid "SPGrid"
 */

namespace Dune
{

  // Capabilities
  // ------------

  /** \brief namespace containing all capabilities */
  namespace Capabilities
  {

    /** \brief Do elements of a grid always have the same geometry type?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, template< int > class Ref, class Comm >
    struct hasSingleGeometryType< SPGrid< ct, dim, Ref, Comm > >
    {
      /** \brief all elements in \ref Dune::SPGrid "SPGrid" have the same
       *         geometry type */
      static const bool v = true;
      /** \brief \ref Dune::SPGrid "SPGrid" has only cube elements */
      static const unsigned int topologyId = GeometryTypes::cube(dim).id();
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
    template< class ct, int dim, template< int > class Ref, class Comm >
    struct isCartesian< SPGrid< ct, dim, Ref, Comm > >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" is a Cartesian grid */
      static const bool v = true;
    };

    /** \brief Does a grid implement entities of a codimension?
     *
     *  \tparam  Grid   grid for which the information is desired
     *  \tparam  codim  codimension in question
     */
    template< class ct, int dim, template< int > class Ref, class Comm, int codim >
    struct hasEntity< SPGrid< ct, dim, Ref, Comm >, codim >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" implements entities for all
       *         codimensions */
      static const bool v = ((codim >= 0) && (codim <= dim));
    };

    /** \brief Does a grid implement entity iterators of a codimension?
     *
     *  \tparam  Grid   grid for which the information is desired
     *  \tparam  codim  codimension in question
     */
    template< class ct, int dim, template< int > class Ref, class Comm, int codim >
    struct hasEntityIterator< SPGrid< ct, dim, Ref, Comm >, codim >
     : public hasEntity< SPGrid< ct, dim, Ref, Comm >, codim >
    {
    };

#if HAVE_MPI
    /** \brief Can a parallel grid communicate on a given codimension?
     *
     *  \tparam  Grid   grid for which the information is desired
     *  \tparam  codim  codimension in question
     *
     *  \note In order to communicate on a given codimension, the grid has to
     *        implement entities for that codimension.
     */
    template< class ct, int dim, template< int > class Ref, int codim >
    struct canCommunicate< SPGrid< ct, dim, Ref, MPI_Comm >, codim >
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
    template< class ct, int dim, template< int > class Ref, class Comm >
    struct isLevelwiseConforming< SPGrid< ct, dim, Ref, Comm > >
    {
      /** \brief All levels of a \ref Dune::SPGrid "SPGrid" are always conform */
      static const bool v = true;
    };

    /** \brief Is the leaf level of a grid always conform?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, template< int > class Ref, class Comm >
    struct isLeafwiseConforming< SPGrid< ct, dim, Ref, Comm > >
    {
      /** \brief The leaf level of a \ref Dune::SPGrid "SPGrid" are always conform */
      static const bool v = true;
    };

    /** \brief Does a grid provide backup and restore facilities?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, template< int > class Ref, class Comm >
    struct hasBackupRestoreFacilities< SPGrid< ct, dim, Ref, Comm > >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" provides backup and restore facilities */
      static const bool v = true;
    };

    /** \brief Is a grid implementation thread safe?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, template< int > class Ref, class Comm >
    struct threadSafe< SPGrid< ct, dim, Ref, Comm > >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" is not thread safe */
      static const bool v = false;
    };

    /** \brief Is a grid implementation thread safe while not being modified?
     *
     *  \tparam  Grid  grid for which the information is desired
     */
    template< class ct, int dim, template< int > class Ref, class Comm >
    struct viewThreadSafe< SPGrid< ct, dim, Ref, Comm > >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" is not thread safe */
      static const bool v = true;
    };



    // non-standard capabilities (see dune-fem)
    // ----------------------------------------

    template< class Grid >
    struct hasHierarchicIndexSet;

    template< class ct, int dim, template< int > class Ref, class Comm >
    struct hasHierarchicIndexSet< SPGrid< ct, dim, Ref, Comm > >
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
    template< class ct, int dim, template< int > class Ref, class Comm >
    struct supportsCallbackAdaptation< SPGrid< ct, dim, Ref, Comm > >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" supports callback adaptation */
      static const bool v = true;
    };

  } // namespace Capabilities



  // Extensions
  // ----------

  namespace Extensions
  {

    /** \brief Does a grid support superentity iterators of a codimension?
     *
     *  \tparam  Grid   grid for which the information is desired
     *  \tparam  codim  codimension in question
     */
    template< class ct, int dim, template< int > class Ref, class Comm, int codim >
    struct SuperEntityIterator< SPGrid< ct, dim, Ref, Comm >, codim >
    {
      /** \brief \ref Dune::SPGrid "SPGrid" supports superentity iterators for all
       *         codimensions */
      static const bool v = ((codim >= 0) && (codim <= dim));
    };

  } // namespace Extensions

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_CAPABILITIES_HH

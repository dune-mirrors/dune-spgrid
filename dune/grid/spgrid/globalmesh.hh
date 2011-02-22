#ifndef DUNE_SPGRID_GLOBALMESH_HH
#define DUNE_SPGRID_GLOBALMESH_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>

/** \file
 *  \author Martin Nolte
 *  \brief  description of a mesh for the infinite space
 */

namespace Dune
{

  /** \class SPGlobalMesh
   *  \brief description of a mesh for the infinite space
   *
   *  \tparam  ct   coordinate type
   *  \tparam  dim  dimension of the space
   */
  template< class ct, int dim >
  class SPGlobalMesh
  {
    typedef SPGlobalMesh< ct, dim > This;

  public:
    typedef ct ctype;

    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef SPMultiIndex< dimension > MultiIndex;

    /** \brief constructor
     *
     *  \param[in]  h  mesh width
     *
     *  \note The origin will be zero in this case.
     */
    explicit SPGlobalMesh ( const GlobalVector &h )
    : origin_( 0 )
    {
      for( int i = 0; i < dimension; ++i )
        h_[ i ] = std::abs( h[ i ] );
    }

    /** \brief constructor
     *
     *  \param[in]  h  mesh width
     *  \param[in]  x  any mesh point
     *
     *  \note The origin will in general not coincide with x.
     */
    SPGlobalMesh ( const GlobalVector &h, const GlobalVector &x )
    {
      for( int i = 0; i < dimension; ++i )
      {
        h_[ i ] = std::abs( h[ i ] );
        origin_[ i ] = x[ i ] - std::floor( x[ i ] / h_[ i ] ) * h_[ i ];
      }
    }

    /** \brief construct a global mesh from a split rectangle
     *
     *  \param[in]  a      corner of the rectangle
     *  \param[in]  b      corner of the rectangle opposite to a
     *  \param[in]  cells  number of cells in each direction
     */
    SPGlobalMesh ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells )
    {
      for( int i = 0; i < dimension; ++i )
      {
        h_[ i ] = std::abs( (b[ i ] - a[ i ]) / ctype( cells[ i ] ) );
        origin_[ i ] = a[ i ] - std::floor( a[ i ] / h_[ i ] ) * h_[ i ];
      }
    }

    /** \brief obtain the coordinates of a mesh point
     *
     *  \param[in]  id  id of the mesh point
     *
     *  \returns coordinates of the mesh point in space
     */
    GlobalVector coordinate ( const MultiIndex &id ) const
    {
      GlobalVector coordinate( origin_ );
      for( int i = 0; i < dimension; ++i )
        coordinate[ i ] += (id[ i ] / 2) * h_[ i ];
      return coordinate;
    }

    /** \brief obtain the id of the nearest mesh point
     *
     *  \param[in]  x  position of a point in space
     *
     *  \returns id of the nearest mesh point
     */
    MultiIndex id ( const GlobalVector &x )
    {
      MultiIndex id;
      for( int i = 0; i < dimension; ++i )
        id[ i ] = 2*std::round( x[ i ] / h_[ i ] );
      return id;
    }

    /** \brief determine whether this mesh equals another
     *
     *  \param[in]  other    global mesh to test for equality with this one
     *  \param[in]  epsilon  tolerance when comparing floating point numbers
     *
     *  \returns true if both meshes are equal up to the tolerance epsilon
     */
    bool equals ( const This &other, const ctype epsilon ) const
    {
      const ctype eps2 = epsilon*epsilon;
      return ((h_ - other.h_).two_norm2() <= eps2) && ((origin_ - other.origin_).two_norm2() <= eps2);
    }

    /** \brief obtain the origin of the mesh
     */
    const GlobalVector &origin () const { return origin_; }

    /** \brief obtain a refined version of this mesh
     *
     *  \param[in]  refinement  type of refinement to apply
     *
     *  \returns a refined version of the global mesh
     */
    template< SPRefinementStrategy strategy >
    This refine ( const SPRefinement< dimension, strategy > &refinement ) const
    {
      GlobalVector h;
      for( int i = 0; i < dimension; ++i )
        h[ i ] = h_[ i ] / ctype( refinement.factor( i ) );
      return This( h, origin_ );
    }

  private:
    GlobalVector h_;
    GlobalVector origin_;
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GLOBALMESH_HH

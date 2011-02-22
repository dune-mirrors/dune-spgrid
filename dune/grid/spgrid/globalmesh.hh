#ifndef DUNE_SPGRID_GLOBALMESH_HH
#define DUNE_SPGRID_GLOBALMESH_HH

#include <cmath>

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  template< class ct, int dim >
  class GlobalMesh
  {
    typedef GlobalMesh< ct, dim > This;

  public:
    typedef ct ctype;

    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef SPMultiIndex< dimension > MultiIndex;

    explicit GlobalMesh ( const GlobalVector &h )
    : origin_( 0 )
    {
      for( int i = 0; i < dimension; ++i )
        h_[ i ] = std::abs( h[ i ] );
    }

    GlobalMesh ( const GlobalVector &h, const GlobalVector &origin )
    {
      for( int i = 0; i < dimension; ++i )
      {
        h_[ i ] = std::abs( h[ i ] );
        origin_[ i ] = origin[ i ] - std::floor( origin[ i ] / h_[ i ] ) * h_[ i ];
      }
    }

    GlobalMesh ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells )
    {
      for( int i = 0; i < dimension; ++i )
      {
        h_[ i ] = std::abs( (b[ i ] - a[ i ]) / cells[ i ] );
        origin_[ i ] = a[ i ] - std::floor( a[ i ] / h_[ i ] ) * h_[ i ];
      }
    }

    GlobalVector coordinate ( const MultiIndex &id ) const
    {
      GlobalVector coordinate( origin_ );
      for( int i = 0; i < dimension; ++i )
        coordinate[ i ] += (id[ i ] / 2) * h_[ i ];
      return coordinate;
    }

    MultiIndex id ( const GlobalVector &x )
    {
      MultiIndex id;
      for( int i = 0; i < dimension; ++i )
        id[ i ] = 2*std::round( x[ i ] / h_[ i ] );
      return id;
    }

    bool equals ( const This &other, const ctype epsilon ) const
    {
      const ctype eps2 = epsilon*epsilon;
      return ((h_ - other.h_).two_norm2() <= eps2) && ((origin_ - other.origin_).two_norm2() <= eps2);
    }

    const GlobalVector &origin () const { return origin_; }

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

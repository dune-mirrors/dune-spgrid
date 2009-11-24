#ifndef DUNE_SPGRID_DOMAIN_HH
#define DUNE_SPGRID_DOMAIN_HH

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // SPDomain
  // --------

  template< class ct, int dim >
  class SPDomain
  {
    typedef SPDomain< ct, dim > This;

  public:
    typedef ct ctype;

    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef SPMultiIndex< dimension > MultiIndex;

    SPDomain ()
    : periodic_( 0 )
    {
      for( int i = 0; i < dimension; ++i )
      {
        origin_[ i ] = ctype( 0 );
        width_[ i ] = ctype( 1 );
        cells_[ i ] = 1;
      }
    }

    SPDomain ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells,
               const unsigned int periodic = 0 )
    : cells_( cells ),
      periodic_( periodic )
    {
      for( int i = 0; i < dimension; ++i )
      {
        origin_[ i ] = std::min( a[ i ], b[ i ] );
        width_[ i ] = std::max( a[ i ], b[ i ] ) - origin_[ i ];
      }
    }

    template< SPRefinementStrategy strategy >
    SPDomain ( const This &other, const SPRefinement< dimension, strategy > &refinement )
    : origin_( other.origin_ ),
      width_( other.width_ ),
      periodic_( other.periodic_ )
    {
      for( int i = 0; i < dimension; ++i )
        cells_[ i ] = refinement.factor( i ) * other.cells_[ i ];
    }

    const GlobalVector &origin () const
    {
      return origin_;
    }

    const GlobalVector &width () const
    {
      return width_;
    }

    const MultiIndex &cells () const
    {
      return cells_;
    }

    GlobalVector h () const
    {
      GlobalVector h;
      for( int i = 0; i < dimension; ++i )
        h[ i ] = width_[ i ] / ctype( cells_[ i ] );
      return h;
    }

    bool periodic ( const int i ) const
    {
      assert( (i >= 0) && (i < dimension) );
      return ((periodic_ & (1 << i)) != 0);
    }

    unsigned int periodic () const
    {
      return periodic_;
    }

  private:
    GlobalVector origin_, width_;
    MultiIndex cells_;
    unsigned int periodic_;
  };

}

#endif // #ifndef DUNE_SPGRID_DOMAIN_HH

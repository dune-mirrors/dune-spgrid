#ifndef DUNE_SPGRID_DOMAIN_HH
#define DUNE_SPGRID_DOMAIN_HH

#include <dune/common/fvector.hh>

namespace Dune
{

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
    {
      for( int i = 0; i < dimension; ++i )
      {
        origin_[ i ] = ctype( 0 );
        width_[ i ] = ctype( 1 );
        cells_[ i ] = 1;
      }
    }

    SPDomain ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells )
    : cells_( cells )
    {
      for( int i = 0; i < dimension; ++i )
      {
        origin_[ i ] = std::min( a[ i ], b[ i ] );
        width_[ i ] = std::max( a[ i ], b[ i ] ) - origin_[ i ];
      }
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

  private:
    GlobalVector origin_, width_;
    MultiIndex cells_;
  };

}

#endif // #ifndef DUNE_SPGRID_DOMAIN_HH

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

    SPDomain ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells,
               const unsigned int periodic = 0 )
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

    GlobalVector h () const
    {
      GlobalVector h;
      for( int i = 0; i < dimension; ++i )
        h[ i ] = width_[ i ] / ctype( cells_[ i ] );
      return h;
    }

  private:
    GlobalVector origin_, width_;
    MultiIndex offset_, cells_;
  };


#if 0
  template< class ct, int dim >
  void SPDomain< ct, dim >::decompose ( int rank, int size )
  {
    undecompose();

    assert( (rank >= 0) && (rank < size) );
    if( size > 1 )
    {
      int dir = 0;
      for( int i = 1; i < dimension; ++i )
        dir = (localCells_[ i ] > localCells_[ dir ] ? i : dir);

      const int lsize = size / 2;
      const int cells = localCells_[ dir ];
      const int lcells = lsize * cells / size;
      const ctype lwidth = ctype( lcells )*localWidth_[ dir ] / ctype( cells );
      if( rank < lsize )
      {
        localCells_[ dir ] = lcells;
        localWidth_[ dir ] = lwidth;
        decompose( rank, lsize );
      }
      else
      {
        localCells_[ dir ] -= lcells;
        localWidth_[ dir ] -= lwidth;
        localOrigin_[ dir ] += lwidth;
        decompose( rank - lsize, size - lsize ); 
      }
    }
  }
#endif

}

#endif // #ifndef DUNE_SPGRID_DOMAIN_HH

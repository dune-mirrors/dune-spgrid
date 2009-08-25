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
        offset_[ i ] = 0;

        neighbor_[ 2*i ] = 0;
        neighbor_[ 2*i+1 ] = 0;
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
        offset_[ i ] = 0;

        This *neighbor = ((periodic & (1 << i)) != 0 ? this : 0);
        neighbor_[ 2*i ] = neighbor;
        neighbor_[ 2*i+1 ] = neighbor;
      }
    }

#if 0
    void decompose ( std::vector< This > &decomposition )
    {
      const unsigned int size = decomposition.size();
      assert( size > 0 );
      decomposition[ 0 ] = *this;
      decompose( decomposition, 0, size );
    }
#endif

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
#if 0
    static void decompose ( std::vector< This > &decomposition,
                            const unsigned int offset, const unsigned int size );
#endif

    GlobalVector origin_, width_;
    MultiIndex offset_, cells_;

    This *neighbor[ 2*dimension ];
  };


#if 0
  template< class ct, int dim >
  inline void
  SPDomain< ct, dim >::decompose ( std::vector< This > &decomposition,
                                   const unsigned int offset, const unsigned int size )
  {
    assert( offset+size <= decomposition.size() );
    if( size > 1 )
    {
      const int leftOffset = offset;
      const int leftSize = size / 2;

      const int rightOffset = leftOffset + leftSize;
      const int rightSize = size - leftSize;

      This &left = decomposition[ leftOffset ];
      This &right = decomposition[ rightOffset ];

      int dir = 0;
      for( int i = 1; i < dimension; ++i )
        dir = argmax( left.cells_, i, dir );

      const int cells = left.cells_[ dir ];
      const int leftCells = (leftSize * cells) / size;
      const int rightCells = cells - leftCells;

      right = left;

      left.cells_[ dir ] = leftCells;
      left.width_[ dir ] *= ctype( leftCells ) / ctype( cells );

      right.cells_[ dir ] = rightCells;
      right.width_[ dir ] *= ctype( leftCells ) / ctype( cells );

      right.offset_[ dir ] += left.width_[ dir ];

      decompose( decomposition, leftOffset, leftSize );
      decompose( decomposition, rightOffset, rightSize );
    }
  }
#endif

}

#endif // #ifndef DUNE_SPGRID_DOMAIN_HH

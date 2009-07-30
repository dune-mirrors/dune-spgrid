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
        globalOrigin_[ i ] = ctype( 0 );
        globalWidth_[ i ] = ctype( 1 );
        globalCells_[ i ] = 1;
      }
      undecompose();
    }

    SPDomain ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells )
    : globalCells_( cells )
    {
      for( int i = 0; i < dimension; ++i )
      {
        globalOrigin_[ i ] = std::min( a[ i ], b[ i ] );
        globalWidth_[ i ] = std::max( a[ i ], b[ i ] ) - globalOrigin_[ i ];
      }
      undecompose();
    }

    void undecompose ()
    {
      localOrigin_ = globalOrigin_;
      localWidth_ = globalWidth_;
      localCells_ = globalCells_;
    }

    void decompose ( int rank, int size );

    const GlobalVector &origin () const
    {
      return localOrigin_;
    }

    const GlobalVector &width () const
    {
      return localWidth_;
    }

    const MultiIndex &cells () const
    {
      return localCells_;
    }

    GlobalVector h () const
    {
      GlobalVector h;
      for( int i = 0; i < dimension; ++i )
        h[ i ] = globalWidth_[ i ] / ctype( globalCells_[ i ] );
      return h;
    }

  private:
    GlobalVector globalOrigin_, globalWidth_;
    GlobalVector localOrigin_, localWidth_;
    MultiIndex globalCells_, localCells_;
  };


  
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

      const int cells = localCells_[ dir ];
      const int lcells = cells / 2;
      const ctype lwidth = ctype( lcells )*localWidth_[ dir ] / ctype( cells );
      const int srank = size / 2;
      if( rank < srank )
      {
        localCells_[ dir ] = lcells;
        localWidth_[ dir ] = lwidth;
        decompose( rank, srank );
      }
      else
      {
        localCells_[ dir ] -= lcells;
        localWidth_[ dir ] -= lwidth;
        localOrigin_[ dir ] += lwidth;
        decompose( rank - srank, size - srank ); 
      }
    }
  }

}

#endif // #ifndef DUNE_SPGRID_DOMAIN_HH

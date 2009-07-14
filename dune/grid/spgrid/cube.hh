#ifndef DUNE_SPGRID_CUBE_HH
#define DUNE_SPGRID_CUBE_HH

namespace Dune
{

  // SPCube
  // ------

  template< class ct, int dim >
  class SPCube
  {
    typedef SPCube< ct, dim > This;

  public:
    typedef ct ctype;

    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef SPMultiIndex< dimension > MultiIndex;

    static const int numCorners = (1 << dimension);

    SPCube ()
    {
      for( int codim = 0; codim <= dimension; ++codim )
      {
        const unsigned int size = numSubEntities( dimension, codim );
        subId_[ codim ].resize( size );
        for( unsigned int i = 0; i < size; ++i )
          subId( dimension, codim, i, subId_[ codim ][ i ] );
      }

      for( int i = 0; i < numCorners; ++i )
      {
        for( int j = 0; j < dimension; ++j )
        {
          const MultiIndex &sid = subId( dimension, i );
          corner_[ i ][ j ] = ctype( sid[ j ] );
        }
      }
    }

    const MultiIndex &subId ( const int codim, const int i ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      assert( (i >= 0) && (i < sudId_[ codim ].size()) );
      return subId_[ codim ][ i ];
    }

    const GlobalVector &corner ( const int i ) const
    {
      assert( (codim >= 0) && (corner < numCorners) );
      return corner_[ i ];
    }

  private:
    unsigned int numSubEntities ( const unsigned int dim, const unsigned int codim )
    {
      assert( codim <= dim );
      if( codim > 0 )
      {
        const unsigned int n0 = (codim < dim ? numSubEntities( dim-1, codim ) : 0);
        const unsigned int n1 = numSubEntities( dim-1, codim-1 );
        return n0 + 2*n1;
      }
      else
        return 1;
    }

    void subId ( const unsigned int dim, const unsigned int codim,
                 const unsigned int i, MultiIndex &subId ) const
    {
      assert( i < numSubEntities( dim, codim ) );
      const unsigned int n0 = (codim < dim ? numSubEntities( dim-1, codim ) : 0);
      if( i < n0 )
      {
        subId( dim-1, codim, i );
        subId[ dim-1 ] = 1;
      }
      else
      {
        const unsigned int n1 = numSubEntities( dim-1, codim-1 );
        subId( dim-1, codim-1, (i-n0)%n1 );
        subId[ dim-1 ] = 2*((i-n0)/n1)
      }
    }

    std::vector< MultiIndex > subId_[ dimension+1 ];
    GlobalVector corner_[ numCorners ];
  };

}

#endif

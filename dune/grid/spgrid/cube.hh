#ifndef DUNE_SPGRID_CUBE_HH
#define DUNE_SPGRID_CUBE_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  namespace SPCubeHelper
  {
    inline unsigned int
    numSubEntities ( const unsigned int dimension, const unsigned int codim )
    {
      assert( codim <= dimension );
      if( codim > 0 )
      {
        const unsigned int n0 = (codim < dimension ? numSubEntities( dimension-1, codim ) : 0);
        const unsigned int n1 = numSubEntities( dimension-1, codim-1 );
        return n0 + 2*n1;
      }
      else
        return 1;
    }
  }



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
    static const int numFaces = 2*dimension;

    SPCube ();

    const MultiIndex &subId ( const int codim, const int i ) const
    {
      assert( (i >= 0) && (i < count( codim )) );
      return subId_[ codim ][ i ];
    }

    int count ( const int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return subId_[ codim ].size();
    }

    const GlobalVector &corner ( const int i ) const
    {
      assert( (i >= 0) && (i < numCorners) );
      return corner_[ i ];
    }

    const GlobalVector &center () const
    {
      return center_;
    }

    const GlobalVector &normal ( const int i ) const
    {
      assert( (i >= 0) && (i < numFaces) );
      return normal_[ i ];
    }

  private:
    void subId ( const unsigned int dimension, const unsigned int codim,
                 const unsigned int i, MultiIndex &sId ) const
    {
      assert( i < SPCubeHelper::numSubEntities( dimension, codim ) );
      if( dimension == 0 )
        return;

      const unsigned int n0 = (codim < dimension ? SPCubeHelper::numSubEntities( dimension-1, codim ) : 0);
      if( i < n0 )
      {
        subId( dimension-1, codim, i, sId );
        sId[ dimension-1 ] = 0;
      }
      else
      {
        const unsigned int n1 = SPCubeHelper::numSubEntities( dimension-1, codim-1 );
        subId( dimension-1, codim-1, (i-n0)%n1, sId );
        sId[ dimension-1 ] = 2*((i-n0)/n1) - 1;
      }
    }

    std::vector< MultiIndex > subId_[ dimension+1 ];
    GlobalVector corner_[ numCorners ];
    GlobalVector center_;
    GlobalVector normal_[ numFaces ];
  };



  template< class ct, int dim >
  inline SPCube< ct, dim >::SPCube ()
  {
    for( int codim = 0; codim <= dimension; ++codim )
    {
      const unsigned int size = SPCubeHelper::numSubEntities( dimension, codim );
      subId_[ codim ].resize( size );
      for( unsigned int i = 0; i < size; ++i )
        subId( dimension, codim, i, subId_[ codim ][ i ] );
    }

    for( int i = 0; i < numCorners; ++i )
    {
      for( int j = 0; j < dimension; ++j )
      {
        const MultiIndex &sid = subId( dimension, i );
        corner_[ i ][ j ] = ctype( (sid[ j ]+1)/2 );
      }
    }

    for( int j = 0; j < dimension; ++j )
      center_[ j ] = ctype( 1 ) / ctype( 2 );

    for( int i = 0; i < numFaces; ++i )
    {
      for( int j = 0; j < dimension; ++j )
      {
        const MultiIndex &sid = subId( 1, i );
        normal_[ i ][ j ] = ctype( sid[ j ] );
      }
    }
  }



  // SPCube (for dim = 0)
  // --------------------

  template< class ct >
  class SPCube< ct, 0 >
  {
    typedef SPCube< ct, 0 > This;

  public:
    typedef ct ctype;

    static const int dimension = 0;

    typedef FieldVector< ctype, 0 > GlobalVector;

    static const int numCorners = 1;

    SPCube ()
    : corner_( ctype( 0 ) )
    {}

    int count ( const int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return 1;
    }

    const GlobalVector &corner ( const int i ) const
    {
      assert( (i >= 0) && (i < numCorners) );
      return corner_;
    }

  private:
    GlobalVector corner_;
  };

}

#endif // #ifndef DUNE_SPGRID_CUBE_HH

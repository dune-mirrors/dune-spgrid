#ifndef DUNE_SPGRID_REFERENCECUBE_HH
#define DUNE_SPGRID_REFERENCECUBE_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/normal.hh>

namespace Dune
{

  namespace SPReferenceCubeHelper
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



  // SPReferenceCube
  // ---------------

  template< class ct, int dim >
  class SPReferenceCube
  {
    typedef SPReferenceCube< ct, dim > This;

  public:
    typedef ct ctype;

    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef SPNormalVector< ctype, dimension > NormalVector;
    typedef SPMultiIndex< dimension > MultiIndex;

    static const int numCorners = (1 << dimension);
    static const int numFaces = 2*dimension;

    SPReferenceCube ();

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

    NormalVector normal ( const int i ) const
    {
      assert( (i >= 0) && (i < numFaces) );
      NormalVector normal( i / 2, 2*(i&1)-1 );
      assert( normal == normal_[ i ] );
      return normal;
    }

  private:
    void subId ( const unsigned int dimension, const unsigned int codim,
                 const unsigned int i, MultiIndex &sId ) const
    {
      using SPReferenceCubeHelper::numSubEntities;
      assert( i < numSubEntities( dimension, codim ) );
      if( dimension == 0 )
        return;

      const unsigned int n0 = (codim < dimension ? numSubEntities( dimension-1, codim ) : 0);
      if( i < n0 )
      {
        subId( dimension-1, codim, i, sId );
        sId[ dimension-1 ] = 0;
      }
      else
      {
        const unsigned int n1 = numSubEntities( dimension-1, codim-1 );
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
  inline SPReferenceCube< ct, dim >::SPReferenceCube ()
  {
    for( int codim = 0; codim <= dimension; ++codim )
    {
      const unsigned int size = SPReferenceCubeHelper::numSubEntities( dimension, codim );
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



  // SPReferenceCube (for dim = 0)
  // -----------------------------

  template< class ct >
  class SPReferenceCube< ct, 0 >
  {
    typedef SPReferenceCube< ct, 0 > This;

  public:
    typedef ct ctype;

    static const int dimension = 0;

    typedef FieldVector< ctype, 0 > GlobalVector;

    static const int numCorners = 1;

    SPReferenceCube ()
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

#endif // #ifndef DUNE_SPGRID_REFERENCECUBE_HH

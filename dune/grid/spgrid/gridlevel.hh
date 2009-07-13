#ifndef DUNE_SPGRID_GRIDLEVEL_HH
#define DUNE_SPGRID_GRIDLEVEL_HH

#include <cassert>
#include <vector>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/codimtable.hh>

#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/domain.hh>

namespace Dune
{

  template< class ct, int dim >
  class SPGridLevel
  {
    typedef SPGridLevel< ct, dim > This;

  public:
    typedef SPDomain< ct, dim > Domain;
    typedef SPGridLevel< ct, dim > GridLevel;

    typedef typename Domain::ctype ctype;
    static const int dimension = Domain::dimension;

    typedef typename Domain::GlobalVector GlobalVector;
    typedef unsigned int MultiIndex[ dimension ];

  public:
    template< int codim >
    struct GeometryCache;
    
  private:
    typedef GenericGeometry::CodimTable< GeometryCache, dimension > GeometryCacheTable;

  public:
    SPGridLevel ( const Domain &domain, const MultiIndex &n )
    : father_( 0 ),
      domain_( &domain ),
      level_( 0 )
    {
      const GlobalVector &width  = domain.width();
      for( int i = 0; i < dimension; ++i )
      {
        cells_[ i ] = n[ i ];
        h_[ i ] = width[ i ] / (ctype)n[ i ];
      }
      GenericGeometry::ForLoop< GeometryCache, 0, dimension >::apply( h_, geometryCache_ );
    }

    SPGridLevel ( const GridLevel &father )
    : father_( &father ),
      domain_( father.domain_ ),
      level_( father.level_+1 )
    {
      for( int i = 0; i < dimension; ++i )
      {
        cells_[ i ] = 2*father.cells_[ i ];
        h_[ i ] = 0.5*father.h_[ i ];
      }
      for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
        geometryCache_[ dir ] = father.geometryCache_[ dir ];
    }

    ~SPGridLevel ()
    {
      if( father_ == 0 )
      {
        for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
          delete geometryCache_[ dir ];
      }
    }

    const Domain &domain () const
    {
      return *domain_;
    }

    const GridLevel &father () const
    {
      assert( father_ != 0 );
      return *father_;
    }

    const GlobalVector &h () const
    {
      return h_;
    }

    unsigned int level () const
    {
      return level_;
    }

    const MultiIndex &cells () const
    {
      return cells_;
    }

    template< int codim >
    const GeometryCache< codim > &geometryCache ( const unsigned int dir ) const
    {
      assert( bitcount( dir ) == dimension - codim );
      return (const GeometryCache< codim > &)( *geometryCache_[ dir ] );
    }

  private:
    const GridLevel *father_;
    const Domain *domain_;
    unsigned int level_;
    MultiIndex cells_;
    GlobalVector h_;
    void *geometryCache_[ 1 << dimension ];
  };



  template< class ct, int dim >
  template< int codim >
  class SPGridLevel< ct, dim >::GeometryCache
  {
    typedef GeometryCache< codim > This;

    friend class SPGridLevel< ct, dim >;

  public:
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;

    typedef FieldVector< ctype, mydimension > LocalVector;
    typedef FieldMatrix< ctype, dimension, mydimension > Jacobian;
    typedef FieldMatrix< ctype, mydimension, dimension > JacobianTransposed;

  private:
    static void
    apply ( const GlobalVector &h, void *(&geometryCache)[ 1 << dimension ] )
    {
      for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
      {
        const int mydim = bitcount( dir );
        if( mydim == dimension - codim )
          geometryCache[ dir ] = new This( h, dir );
      }
    }

    GeometryCache ( const GlobalVector &h, const unsigned int dir )
    : volume_( 1 ),
      jacobianTransposed_( 0 ),
      jacobianInverseTransposed_( 0 )
    {
      int k = 0;
      for( int j = 0; j < dimension; ++j )
      {
        if( ((dir >> j) & 1) == 0 )
          continue;
        volume_ *= h[ j ];
        jacobianTransposed_[ j ][ k ] = h[ j ];
        jacobianInverseTransposed_[ k ][ j ] = ctype( 1 ) / h[ j ];
        ++k;
      }
    }

  public:
    const ctype &volume () const
    {
      return volume_;
    }

    const JacobianTransposed &jacobianTransposed () const
    {
      return jacobianTransposed_;
    }

    const Jacobian &jacobianInverseTransposed () const
    {
      return jacobianInverseTransposed_;
    }

  private:
    ctype volume_;
    JacobianTransposed jacobianTransposed_;
    Jacobian jacobianInverseTransposed_;
  };

}

#endif // #ifndef DUNE_SPGRID_GRIDLEVEL_HH

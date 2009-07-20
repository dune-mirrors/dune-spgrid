#ifndef DUNE_SPGRID_GRIDLEVEL_HH
#define DUNE_SPGRID_GRIDLEVEL_HH

#include <cassert>
#include <vector>

#include <dune/common/smartpointer.hh>

#include <dune/grid/genericgeometry/misc.hh>

#include <dune/grid/spgrid/cube.hh>
#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/domain.hh>
#include <dune/grid/spgrid/geometrycache.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int mydim, int cdim, class Grid >
  class SPLocalGeometry;



  // SPGridLevel
  // -----------

  template< class Grid >
  class SPGridLevel
  {
    typedef SPGridLevel< Grid > This;

  public:
    typedef SPGridLevel< Grid > GridLevel;

    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::Cube Cube;
    typedef typename Traits::Domain Domain;

    typedef typename Cube::ctype ctype;
    static const int dimension = Cube::dimension;

    typedef typename Cube::GlobalVector GlobalVector;
    typedef typename Cube::MultiIndex MultiIndex;

    static const unsigned int numDirections = Cube::numCorners;

    template< int codim >
    struct Codim
    {
      typedef SPCube< ctype, dimension-codim > Cube;
      typedef SPGeometryCache< ctype, dimension, codim > GeometryCache;
    };

    typedef typename Traits::template Codim< 0 >::LocalGeometry LocalGeometry;

  private:
    typedef SPLocalGeometry< dimension, dimension, Grid > LocalGeometryImpl;

    template< int codim >
    struct BuildGeometryCache;

  public:
    SPGridLevel ( const Grid &grid, const MultiIndex &n );

    SPGridLevel ( GridLevel &father, const unsigned int refDir );

  private:
    SPGridLevel ( const This &other );

  public:
    ~SPGridLevel ()
    {
      for( unsigned int dir = 0; dir < numDirections; ++dir )
      {
        delete geometryInFather_[ dir ];
        delete geometryCache_[ dir ];
      }
    }

    const Grid &grid () const
    {
      return *grid_;
    }

    const Cube &cube () const
    {
      return grid().cube();
    }

    template< int codim >
    const typename Codim< codim >::Cube &cube () const
    {
      return grid().template cube< codim >();
    }

    const Domain &domain () const
    {
      return grid().domain();
    }

    const GridLevel &fatherLevel () const
    {
      assert( father_ != 0 );
      return *father_;
    }

    const GridLevel &childLevel () const
    {
      assert( !isLeaf() );
      return *child_;
    }

    bool isLeaf () const
    {
      return (child_ == 0);
    }

    const GlobalVector &h () const
    {
      return h_;
    }

    unsigned int level () const
    {
      return level_;
    }

    unsigned int refinementDirection () const
    {
      return refDir_;
    }

    const MultiIndex &cells () const
    {
      return cells_;
    }

    MultiIndex fatherId ( const MultiIndex &id ) const
    {
      MultiIndex fatherId;
      for( int i = 0; i < dimension; ++i )
        fatherId[ i ] = ((refDir_ >> i) & 1 ? (id[ i ] >> 1) | 1 : id[ i ]);
      return fatherId;
    }

    const LocalGeometry &geometryInFather ( const MultiIndex &id ) const
    {
      const unsigned int childIndex = 0;
      for( int i = 0; i < dimension; ++i )
        childIndex |= ((refDir_ >> i) & (id[ i ] >> 1) & 1) << i;
      return geometryInFather_[ childIndex ];
    }

    template< int codim >
    const typename Codim< codim >::GeometryCache &
    geometryCache ( const unsigned int dir ) const
    {
      typedef typename Codim< codim >::GeometryCache GeometryCache;
      assert( bitCount( dir ) == dimension - codim );
      return *((const GeometryCache *)geometryCache_[ dir ]);
    }

    const GlobalVector &volumeNormal ( const int i ) const
    {
      assert( (i >= 0) && (i < Cube::numFaces) );
      return normal_[ i ];
    }

  private:
    void buildGeometry ();

    const Grid *grid_;
    GridLevel *father_, *child_;

    unsigned int level_;
    unsigned int refDir_;
    MultiIndex cells_;
    GlobalVector h_;

    void *geometryCache_[ numDirections ];
    LocalGeometry *geometryInFather_[ numDirections ];
    GlobalVector normal_[ Cube::numFaces ];
  };


  template< class Grid >
  inline SPGridLevel< Grid >
    ::SPGridLevel ( const Grid &grid, const MultiIndex &n )
  : grid_( &grid ),
    father_( 0 ),
    child_( 0 ),
    level_( 0 ),
    refDir_( 0 )
  {
    const GlobalVector &width  = domain().width();
    for( int i = 0; i < dimension; ++i )
    {
      cells_[ i ] = n[ i ];
      h_[ i ] = width[ i ] / (ctype)n[ i ];
    }
    for( unsigned int dir = 0; dir < numDirections; ++dir )
      geometryInFather_[ dir ] = 0;
    buildGeometry();
  }


  template< class Grid >
  inline SPGridLevel< Grid >
    ::SPGridLevel ( GridLevel &father, const unsigned int refDir )
  : grid_( father.grid_ ),
    father_( &father ),
    child_( 0 ),
    level_( father.level_ + 1 ),
    refDir_( refDir )
  {
    assert( father.child_ == 0 );
    father.child_ = *this;
    GlobalVector hInFather;
    for( int i = 0; i < dimension; ++i )
    {
      const unsigned int factor = 2*((refDir >> i) & 1);
      cells_[ i ] = factor * father.cells_[ i ];
      hInFather[ i ] = ctype( 1 ) / ctype( factor );
      h_[ i ] = father.h_[ i ] * hInFather[ i ];
    }

    const typename Codim< 0 >::GeometryCache cacheInFather( hInFather, numDirections-1 );
    for( unsigned int dir = 0; dir < numDirections; ++dir )
    {
      geometryInFather_[ dir ] = 0;
      if( (dir & refDir) != dir )
        continue;

      GlobalVector origin;
      for( int i = 0; i < dimension; ++i )
        origin[ i ] = ctype( (refDir >> i) & 1 ) / ctype( 2 );
      geometryInFather_[ dir ] = new LocalGeometry( LocalGeometryImpl( cube(), cacheInFather, origin ) );
    }

    buildGeometry();
  }


  template< class Grid >
  inline void SPGridLevel< Grid >::buildGeometry ()
  {
    GenericGeometry::ForLoop< BuildGeometryCache, 0, dimension >::apply( h_, geometryCache_ );
    
    const ctype volume = geometryCache< 0 >( numDirections-1 ).volume();
    for( int face = 0; face < Cube::numFaces; ++face )
    {
      normal_[ face ] = cube().normal( face );
      const ctype hn = std::abs( normal_[ face ] * h_ );
      normal_[ face ] *= volume / hn;
    }
  }



  // SPGridLevel::BuildGeometryCache
  // -------------------------------

  template< class Grid >
  template< int codim >
  struct SPGridLevel< Grid >::BuildGeometryCache
  {
    static void
    apply ( const GlobalVector &h, void *(&geometryCache)[ 1 << dimension ] )
    {
      typedef typename Codim< codim >::GeometryCache GeometryCache;
      for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
      {
        const int mydim = bitCount( dir );
        if( mydim == dimension - codim )
          geometryCache[ dir ] = new GeometryCache( h, dir );
      }
    }
  };

}

#endif // #ifndef DUNE_SPGRID_GRIDLEVEL_HH

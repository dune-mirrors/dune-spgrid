#ifndef DUNE_SPGRID_GRIDLEVEL_HH
#define DUNE_SPGRID_GRIDLEVEL_HH

#include <cassert>
#include <vector>

#include <dune/common/forloop.hh>

#include <dune/grid/spgrid/cube.hh>
#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/refinement.hh>
#include <dune/grid/spgrid/domain.hh>
#include <dune/grid/spgrid/mesh.hh>
#include <dune/grid/spgrid/partitionlist.hh>
#include <dune/grid/spgrid/decomposition.hh>
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
    typedef typename Traits::Refinement Refinement;

    typedef typename Cube::ctype ctype;
    static const int dimension = Cube::dimension;

    typedef typename Cube::GlobalVector GlobalVector;
    typedef typename Cube::MultiIndex MultiIndex;

    typedef SPDecomposition< dimension > Decomposition;
    typedef SPPartitionPool< dimension > PartitionPool;

    typedef typename Decomposition::Mesh Mesh;

    typedef typename PartitionPool::PartitionList PartitionList;

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
    template< int codim >
    struct DestroyGeometryCache;

  public:
    SPGridLevel ( const Grid &grid, const Decomposition &decomposition );

    SPGridLevel ( GridLevel &father, const Refinement &refinement );

  private:
    SPGridLevel ( const This &other );

  public:
    ~SPGridLevel ();

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
      return domain_;
    }

    const Mesh &globalMesh () const
    {
      return globalMesh_;
    }

    const Mesh &localMesh () const
    {
      return localMesh_;
    }

    template< PartitionIteratorType pitype >
    const PartitionList &partition () const
    {
      return partitionPool_.template get< pitype >();
    }

    template< int codim >
    PartitionType partitionType ( const MultiIndex &id ) const
    {
      return partitionPool_.template partitionType< codim >( id );
    }

    const GridLevel &fatherLevel () const
    {
      assert( !isMacro() );
      return *father_;
    }

    const GridLevel &childLevel () const
    {
      assert( !isLeaf() );
      return *child_;
    }

    bool isMacro () const
    {
      return (father_ == 0 );
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

    const Refinement &refinement () const
    {
      assert( !isMacro() );
      return refinement_;
    }

    MultiIndex macroId ( const MultiIndex &id ) const
    {
      MultiIndex macroId;
      for( int i = 0; i < dimension; ++i )
        macroId[ i ] = (((id[ i ] >> 1) / macroFactor_[ i ]) << 1) | (id[ i ] & 1);
      return macroId;
    }

    const MultiIndex &cells () const
    {
      return domain().cells();
    }

    size_t boundaryIndex ( const MultiIndex &id, const int face ) const
    {
      // note: boundaryIndex ignores the last bit of macroId,
      //       hence we can use this fast computation
      MultiIndex macroId;
      for( int i = 0; i < dimension; ++i )
        macroId[ i ] = id[ i ] / macroFactor_[ i ];
      return grid().boundaryIndex( macroId, face );
    }

    const LocalGeometry &geometryInFather ( const MultiIndex &id ) const
    {
      assert( !isMacro() && (geometryInFather_ != 0) );
      return *(geometryInFather_[ refinement().childIndex( id ) ]);
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

    const int size () const
    {
      const int size = 1;
      for( int i = 0; i < dimension; ++i )
        size *= cells()[ i ];
      return size;
    }

  private:
    void buildGeometry ();

    const Grid *grid_;
    GridLevel *father_, *child_;

    unsigned int level_;
    const Refinement refinement_;
    unsigned int macroFactor_[ dimension ];
    Domain domain_;
    Mesh globalMesh_;
    Mesh localMesh_;
    PartitionPool partitionPool_;
    GlobalVector h_;

    void *geometryCache_[ numDirections ];
    LocalGeometry **geometryInFather_;
    GlobalVector normal_[ Cube::numFaces ];
  };


  template< class Grid >
  inline SPGridLevel< Grid >
    ::SPGridLevel ( const Grid &grid, const Decomposition &decomposition )
  : grid_( &grid ),
    father_( 0 ),
    child_( 0 ),
    level_( 0 ),
    refinement_(),
    domain_( grid.domain() ),
    globalMesh_( decomposition.mesh() ),
    localMesh_( decomposition.subMesh( grid.comm().rank() ) ),
    partitionPool_( localMesh_, globalMesh_, MultiIndex::zero() ),
    h_( domain().h() ),
    geometryInFather_( 0 )
  {
    for( int i = 0; i < dimension; ++i )
      macroFactor_[ i ] = 1;
    buildGeometry();
  }


  template< class Grid >
  inline SPGridLevel< Grid >
    ::SPGridLevel ( GridLevel &father, const Refinement &refinement )
  : grid_( father.grid_ ),
    father_( &father ),
    child_( 0 ),
    level_( father.level_ + 1 ),
    refinement_( refinement ),
    domain_( father.domain(), refinement ),
    globalMesh_( father.globalMesh().refine( refinement ) ),
    localMesh_( father.localMesh().refine( refinement ) ),
    partitionPool_( localMesh_, globalMesh_, MultiIndex::zero() ),
    h_( domain().h() )
  {
    assert( father.child_ == 0 );
    father.child_ = this;

    for( int i = 0; i < dimension; ++i )
      macroFactor_[ i ] = father.macroFactor_[ i ] * refinement.factor( i );

    const unsigned int numChildren = refinement.numChildren();
    geometryInFather_ = new LocalGeometry *[ numChildren ];
    const GlobalVector hInFather = refinement.template hInFather< ctype >();
    const typename Codim< 0 >::GeometryCache cacheInFather( hInFather, numDirections-1 );
    for( unsigned int index = 0; index < numChildren; ++index )
    {
      const GlobalVector origin = refinement.template originInFather< ctype >( index );
      geometryInFather_[ index ] = new LocalGeometry( LocalGeometryImpl( cube(), cacheInFather, origin ) );
    }

    buildGeometry();
  }


  template< class Grid >
  inline SPGridLevel< Grid >::~SPGridLevel ()
  {
    if( child_ != 0 )
      delete child_;
    assert( child_ == 0 );
    if( father_ != 0 )
      father_->child_ = 0;

    if( geometryInFather_ != 0 )
    {
      unsigned int numChildren = refinement().numChildren();
      for( unsigned int index = 0; index < numChildren; ++index )
        delete geometryInFather_[ index ];
      delete geometryInFather_;
    }

    ForLoop< DestroyGeometryCache, 0, dimension >::apply( geometryCache_ );
  }


  template< class Grid >
  inline void SPGridLevel< Grid >::buildGeometry ()
  {
    ForLoop< BuildGeometryCache, 0, dimension >::apply( h_, geometryCache_ );
    
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

  template< class Grid >
  template< int codim >
  struct SPGridLevel< Grid >::DestroyGeometryCache
  {
    static void
    apply ( void *(&geometryCache)[ 1 << dimension ] )
    {
      typedef typename Codim< codim >::GeometryCache GeometryCache;
      for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
      {
        const int mydim = bitCount( dir );
        if( mydim == dimension - codim )
        {
          delete (GeometryCache *)geometryCache[ dir ];
          geometryCache[ dir ] = 0;
        }
      }
    }
  };

}

#endif // #ifndef DUNE_SPGRID_GRIDLEVEL_HH

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
#include <dune/grid/spgrid/partitionpool.hh>
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
    // prohibit copying
    SPGridLevel ( const This &other );

  public:
    ~SPGridLevel ();

    const Grid &grid () const;

    const Cube &cube () const;

    template< int codim >
    const typename Codim< codim >::Cube &cube () const;

    const Domain &domain () const;

    const Mesh &globalMesh () const;
    const Mesh &localMesh () const;

    template< PartitionIteratorType pitype >
    const PartitionList &partition () const;

    template< int codim >
    PartitionType
    partitionType ( const MultiIndex &id, const unsigned int partitionNumber ) const;

    const GridLevel &fatherLevel () const;
    const GridLevel &childLevel () const;
    bool isMacro () const;
    bool isLeaf () const;
    unsigned int level () const;

    const GlobalVector &h () const;
    const Refinement &refinement () const;
    MultiIndex macroId ( const MultiIndex &id ) const;

    size_t boundaryIndex ( const MultiIndex &id,
                           const unsigned int partitionNumber,
                           const int face ) const;

    const LocalGeometry &geometryInFather ( const MultiIndex &id ) const;

    template< int codim >
    const typename Codim< codim >::GeometryCache &
    geometryCache ( const unsigned int dir ) const;

    int size () const;

    const GlobalVector &volumeNormal ( const int i ) const;

  private:
    void buildGeometry ();

    static MultiIndex coarseMacroFactor ();

    MultiIndex overlap () const;

    const Grid *grid_;
    GridLevel *father_, *child_;

    unsigned int level_;
    const Refinement refinement_;
    MultiIndex macroFactor_;
    Domain domain_;
    Mesh globalMesh_;
    Mesh localMesh_;
    PartitionPool partitionPool_;

    GlobalVector h_;
    void *geometryCache_[ numDirections ];
    LocalGeometry **geometryInFather_;
    GlobalVector normal_[ Cube::numFaces ];
  };



  // Implementation of SPGridLevel
  // -----------------------------

  template< class Grid >
  inline SPGridLevel< Grid >
    ::SPGridLevel ( const Grid &grid, const Decomposition &decomposition )
  : grid_( &grid ),
    father_( 0 ),
    child_( 0 ),
    level_( 0 ),
    refinement_(),
    macroFactor_( coarseMacroFactor() ),
    domain_( grid.domain() ),
    globalMesh_( decomposition.mesh() ),
    localMesh_( decomposition.subMesh( grid.comm().rank() ) ),
    partitionPool_( localMesh_, globalMesh_, overlap(), domain_.periodic() ),
    h_( domain().h() ),
    geometryInFather_( 0 )
  {
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
    macroFactor_( father.macroFactor_ * refinement ),
    domain_( father.domain(), refinement ),
    globalMesh_( father.globalMesh().refine( refinement ) ),
    localMesh_( father.localMesh().refine( refinement ) ),
    partitionPool_( localMesh_, globalMesh_, overlap(), domain_.periodic() ),
    h_( domain().h() )
  {
    assert( father.child_ == 0 );
    father.child_ = this;

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
  inline const Grid &SPGridLevel< Grid >::grid () const
  {
    return *grid_;
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::Cube &
  SPGridLevel< Grid >::cube () const
  {
    return grid().cube();
  }


  template< class Grid >
  template< int codim >
  inline const typename SPGridLevel< Grid >::template Codim< codim >::Cube &
  SPGridLevel< Grid >::cube () const
  {
    return grid().template cube< codim >();
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::Domain &
  SPGridLevel< Grid >::domain () const
  {
    return domain_;
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::Mesh &
  SPGridLevel< Grid >::globalMesh () const
  {
    return globalMesh_;
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::Mesh &
  SPGridLevel< Grid >::localMesh () const
  {
    return localMesh_;
  }


  template< class Grid >
  template< PartitionIteratorType pitype >
  inline const typename SPGridLevel< Grid >::PartitionList &
  SPGridLevel< Grid >::partition () const
  {
    return partitionPool_.template get< pitype >();
  }


  template< class Grid >
  template< int codim >
  inline PartitionType SPGridLevel< Grid >
    ::partitionType ( const MultiIndex &id, const unsigned int partitionNumber ) const
  {
    return partitionPool_.template partitionType< codim >( id, partitionNumber );
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::GridLevel &
  SPGridLevel< Grid >::fatherLevel () const
  {
    assert( !isMacro() );
    return *father_;
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::GridLevel &
  SPGridLevel< Grid >::childLevel () const
  {
    assert( !isLeaf() );
    return *child_;
  }


  template< class Grid >
  inline bool SPGridLevel< Grid >::isMacro () const
  {
    return (father_ == 0 );
  }


  template< class Grid >
  inline bool SPGridLevel< Grid >::isLeaf () const
  {
    return (child_ == 0);
  }


  template< class Grid >
  inline unsigned int SPGridLevel< Grid >::level () const
  {
    return level_;
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::GlobalVector &
  SPGridLevel< Grid >::h () const
  {
    return h_;
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::Refinement &
  SPGridLevel< Grid >::refinement () const
  {
    assert( !isMacro() );
    return refinement_;
  }

  
  template< class Grid >
  inline typename SPGridLevel< Grid >::MultiIndex
  SPGridLevel< Grid >::macroId ( const MultiIndex &id ) const
  {
    MultiIndex macroId;
    for( int i = 0; i < dimension; ++i )
      macroId[ i ] = (((id[ i ] >> 1) / macroFactor_[ i ]) << 1) | (id[ i ] & 1);
    return macroId;
  }


  template< class Grid >
  inline size_t SPGridLevel< Grid >
    ::boundaryIndex ( const MultiIndex &id,
                      const unsigned int partitionNumber,
                      const int face ) const
  {
    // note: boundaryIndex ignores the last bit of macroId,
    //       hence we can use this fast computation
    MultiIndex macroId;
    for( int i = 0; i < dimension; ++i )
      macroId[ i ] = id[ i ] / macroFactor_[ i ];
    return grid().boundaryIndex( macroId, partitionNumber, face );
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::LocalGeometry &
  SPGridLevel< Grid >::geometryInFather ( const MultiIndex &id ) const
  {
    assert( !isMacro() && (geometryInFather_ != 0) );
    return *(geometryInFather_[ refinement().childIndex( id ) ]);
  }


  template< class Grid >
  template< int codim >
  inline const typename SPGridLevel< Grid >::template Codim< codim >::GeometryCache &
  SPGridLevel< Grid >::geometryCache ( const unsigned int dir ) const
  {
    typedef typename Codim< codim >::GeometryCache GeometryCache;
    assert( bitCount( dir ) == dimension - codim );
    return *((const GeometryCache *)geometryCache_[ dir ]);
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::GlobalVector &
  SPGridLevel< Grid >::volumeNormal ( const int i ) const
  {
    assert( (i >= 0) && (i < Cube::numFaces) );
    return normal_[ i ];
  }


  template< class Grid >
  inline int SPGridLevel< Grid >::size () const
  {
    return globalMesh().volume();
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


  template< class Grid >
  inline typename SPGridLevel< Grid >::MultiIndex
  SPGridLevel< Grid >::coarseMacroFactor ()
  {
    MultiIndex macroFactor;
    for( int i = 0; i < dimension; ++i )
      macroFactor[ i ] = 1;
    return macroFactor;
  }


  template< class Grid >
  inline typename SPGridLevel< Grid >::MultiIndex
  SPGridLevel< Grid >::overlap () const
  {
    MultiIndex overlap;
    for( int i = 0; i < dimension; ++i )
      overlap[ i ] = macroFactor_[ i ] * grid().overlap()[ i ];
    return overlap;
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

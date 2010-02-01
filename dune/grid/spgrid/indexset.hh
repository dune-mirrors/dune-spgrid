#ifndef DUNE_SPGRID_INDEXSET_HH
#define DUNE_SPGRID_INDEXSET_HH

#include <dune/common/array.hh>

#include <dune/grid/common/indexidset.hh>

#include <dune/grid/spgrid/entityinfo.hh>
#include <dune/grid/spgrid/gridlevel.hh>

namespace Dune
{

  template< class Grid >
  class SPIndexSet
  : public IndexSet< Grid, SPIndexSet< Grid >, unsigned int >
  {
    typedef SPIndexSet< Grid > This;
    typedef IndexSet< Grid, This, unsigned int > Base;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Base::IndexType IndexType;

    static const int dimension = Traits::Cube::dimension;

    template< int codim >
    struct Codim
    {
      typedef SPEntityInfo< Grid, codim > EntityInfo;
      typedef typename Traits::template Codim< codim >::Entity Entity;
    };

    typedef SPGridLevel< Grid > GridLevel;

  private:
    typedef typename GridLevel::MultiIndex MultiIndex;

  public:
    SPIndexSet ()
    : gridLevel_( 0 )
    {
      makeGeomTypes();
    }

    explicit SPIndexSet ( const GridLevel &gridLevel )
    : gridLevel_( 0 )
    {
      makeGeomTypes();
      update( gridLevel );
    }

    void update ( const GridLevel &gridLevel );

  private:
    IndexType index ( const MultiIndex &id,
                      const unsigned int number ) const;

    void makeGeomTypes ()
    {
      for( int codim = 0; codim <= dimension; ++codim )
        geomTypes_[ codim ].push_back( GeometryType( GeometryType::cube, dimension-codim ) );
    }

  public:
    template< class Entity >
    IndexType index ( const Entity &entity ) const
    {
      return index< Entity::codimension >( entity );
    }

    template< int codim >
    IndexType index ( const typename Codim< codim >::Entity &entity ) const;

    IndexType subIndex ( const typename Codim< 0 >::Entity &entity,
                         const int i, const unsigned int codim ) const;

    const std::vector< GeometryType > &geomTypes ( const int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return geomTypes_[ codim ];
    }

    IndexType size ( const GeometryType &type ) const
    {
      return (type.isCube() ? size( dimension - type.dim() ) : 0);
    }

    IndexType size ( const int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return size_[ codim ];
    }

    template< class Entity >
    bool contains ( const Entity &entity ) const
    {
      return contains< Entity::codimension >( entity );
    }

    template< int codim >
    bool contains ( const typename Codim< codim >::Entity &entity ) const
    {
      const typename Codim< codim >::EntityInfo &entityInfo
        = Grid::getRealImplementation( entity ).entityInfo();
      return (&entityInfo.gridLevel() == &gridLevel());
    }

    const GridLevel &gridLevel () const
    {
      return *gridLevel_;
    }

  private:
    const GridLevel *gridLevel_;
    MultiIndex origin_, cells_;
    IndexType offsets_[ 1 << dimension ];
    IndexType size_[ dimension+1 ];
    std::vector< GeometryType > geomTypes_[ dimension+1 ];
  };



  // Implementation of SPIndexSet
  // ----------------------------

  template< class Grid >
  void SPIndexSet< Grid >::update ( const GridLevel &gridLevel )
  {
    gridLevel_ = &gridLevel;
    //origin_ = gridLevel.allPartition().begin()->begin();
    //cells_ = gridLevel.allPartition().begin()->width();
    origin_ = gridLevel.globalMesh().begin();
    cells_ = gridLevel.globalMesh().end() - origin_;

    for( int codim = 0; codim <= dimension; ++codim )
      size_[ codim ] = 0;

    for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
    {
      IndexType factor = 1;
      unsigned int codim = dimension;
      for( int j = 0; j < dimension; ++j )
      {
        const unsigned int d = (dir >> j) & 1;
        factor *= cells_[ j ] + (1-d);
        codim -= d;
      }
      offsets_[ dir ] = size_[ codim ];
      size_[ codim ] += factor;
    }
  }


  template< class Grid >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::index ( const MultiIndex &id,
                              const unsigned int number ) const
  {
    IndexType index = 0;
    IndexType factor = 1;
    unsigned int dir = 0;
    for( int j = 0; j < dimension; ++j )
    {
      const unsigned int d = id[ j ] & 1;
      const IndexType idLocal = (id[ j ] >> 1) - origin_[ j ];
      const IndexType width = cells_[ j ] + (1-d);
      assert( (idLocal >= 0) && (idLocal < width) );
      index += idLocal * factor;
      factor *= width;
      dir |= (d << j);
    }
    return offsets_[ dir ] + index;
  }


  template< class Grid >
  template< int codim >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::index ( const typename Codim< codim >::Entity &entity ) const
  {
    assert( contains( entity ) );
    const typename Codim< codim >::EntityInfo &entityInfo
      = Grid::getRealImplementation( entity ).entityInfo();
    return index( entityInfo.id(), entityInfo.partitionNumber() );
  }

  template< class Grid >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::subIndex ( const typename Codim< 0 >::Entity &entity,
                                 const int i, const unsigned int codim ) const
  {
    assert( contains( entity ) );
    const typename Codim< 0 >::EntityInfo &entityInfo
      = Grid::getRealImplementation( entity ).entityInfo();
    MultiIndex sid = entityInfo.id();
    sid += gridLevel().cube().subId( codim, i );
    // for the ghost approach, the partition number has to be corrected
    return index( sid, entityInfo.partitionNumber() );
  }



  // SPHierarchyIndexSet
  // -------------------

  template< class Grid >
  class SPHierarchyIndexSet
  : public IndexSet< Grid, SPHierarchyIndexSet< Grid >, unsigned int >
  {
    typedef SPHierarchyIndexSet< Grid > This;
    typedef IndexSet< Grid, This, unsigned int > Base;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Base::IndexType IndexType;

    static const int dimension = Traits::Cube::dimension;

    template< int codim >
    struct Codim
    {
      typedef SPEntityInfo< Grid, codim > EntityInfo;
      typedef typename Traits::template Codim< codim >::Entity Entity;
    };

    typedef SPGridLevel< Grid > GridLevel;

  private:
    typedef SPIndexSet< Grid > LevelIndexSet;

    typedef array< IndexType, dimension+1 > CodimIndexArray;

  public:
    explicit SPHierarchyIndexSet ( const Grid &grid )
    : grid_( &grid )
    {
      for( int codim = 0; codim <= dimension; ++codim )
        size_[ codim ] = 0;
    }

    void update ()
    {
      for( int codim = 0; codim <= dimension; ++codim )
        size_[ codim ] = 0;

      const int maxLevel = grid().maxLevel();
      levelIndexSets_.resize( maxLevel+1 );
      offsets_.resize( maxLevel+1 );
      for( int level = 0; level <= maxLevel; ++level )
      {
        const LevelIndexSet &levelIndexSet = grid().levelIndexSet( level );
        levelIndexSets_[ level ] = &levelIndexSet;
        for( int codim = 0; codim <= dimension; ++codim )
        {
          offsets_[ level ][ codim ] = size_[ codim ];
          size_[ codim ] += levelIndexSet.size( codim );
        }
      }
    }

    template< class Entity >
    IndexType index ( const Entity &entity ) const
    {
      return index< Entity::codimension >( entity );
    }

    template< int codim >
    IndexType index ( const typename Codim< codim >::Entity &entity ) const
    {
      const int level = entity.level();
      const IndexType offset = offsets_[ level ][ codim ];
      return offset + levelIndexSet( level ).index( entity );
    }

    template< int codim >
    IndexType DUNE_DEPRECATED
    subIndex ( const typename Codim< 0 >::Entity &entity, const int i ) const
    {
      DUNE_THROW( NotImplemented, "SPHierarchyIndexSet does not implement the old subIndex method." );
    }

    IndexType subIndex ( const typename Codim< 0 >::Entity &entity,
                         const int i, const unsigned int codim ) const
    {
      const int level = entity.level();
      const IndexType offset = offsets_[ level ][ codim ];
      return offset + levelIndexSet( level ).subIndex( entity, i, codim );
    }

    const std::vector< GeometryType > &geomTypes ( const int codim ) const
    {
      return levelIndexSet( 0 ).geomTypes( codim );
    }

    IndexType size ( const GeometryType &type ) const
    {
      return (type.isCube() ? size( dimension - type.dim() ) : 0);
    }

    IndexType size ( const int codim ) const
    {
      assert( (codim >= 0) && (codim <= dimension) );
      return size_[ codim ];
    }

    template< class Entity >
    bool contains ( const Entity &entity ) const
    {
      return contains< Entity::codimension >( entity );
    }

    template< int codim >
    bool contains ( const typename Codim< codim >::Entity &entity ) const
    {
      return true;
    }

    const Grid &grid () const
    {
      return *grid_;
    }

    const LevelIndexSet &levelIndexSet ( const int level ) const
    {
      assert( (level >= 0) && (level <= (int)levelIndexSets_.size()) );
      assert( (int)levelIndexSets_.size() == grid().maxLevel()+1 );
      return *levelIndexSets_[ level ];
    }

  private:
    const Grid *grid_;
    std::vector< const LevelIndexSet * > levelIndexSets_;
    std::vector< CodimIndexArray > offsets_;
    CodimIndexArray size_;
  };

}

#endif // #ifndef DUNE_SPGRID_INDEXSET_HH

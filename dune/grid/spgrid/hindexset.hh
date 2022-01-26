#ifndef DUNE_SPGRID_HINDEXSET_HH
#define DUNE_SPGRID_HINDEXSET_HH

#include <array>
#include <vector>
#include <type_traits>

#include <dune/grid/spgrid/indexset.hh>

namespace Dune
{

  // SPHierarchyIndexSet
  // -------------------

  template< class Grid >
  class SPHierarchyIndexSet
    : public IndexSet< Grid, SPHierarchyIndexSet< Grid >, unsigned int, std::array< GeometryType, 1 > >
  {
    typedef SPHierarchyIndexSet< Grid > This;
    typedef IndexSet< Grid, This, unsigned int, std::array< GeometryType, 1 > > Base;

    typedef typename std::remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Base::IndexType IndexType;
    typedef typename Base::Types Types;

    static const int dimension = Traits::ReferenceCube::dimension;

    template< int codim >
    struct Codim
    {
      typedef __SPGrid::EntityInfo< Grid, codim > EntityInfo;
      typedef typename Traits::template Codim< codim >::Entity Entity;
    };

    typedef SPGridLevel< typename std::remove_const< Grid >::type > GridLevel;

  private:
    typedef SPIndexSet< Grid > LevelIndexSet;

    typedef std::array< IndexType, dimension+1 > CodimIndexArray;

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

    template< class Entity >
    IndexType subIndex ( const Entity &entity,
                         const int i, const unsigned int codim ) const
    {
      return subIndex< Entity::codimension >( entity, i, codim );
    }

    template< int cd >
    IndexType subIndex ( const typename Codim< cd >::Entity &entity,
                         const int i, const unsigned int codim ) const
    {
      const int level = entity.level();
      const IndexType offset = offsets_[ level ][ codim ];
      return offset + levelIndexSet( level ).subIndex( entity, i, codim );
    }

    Types types ( int codim ) const { return {{ GeometryTypes::cube( dimension - codim ) }}; }

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

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_HINDEXSET_HH

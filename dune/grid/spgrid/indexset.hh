#ifndef DUNE_SPGRID_INDEXSET_HH
#define DUNE_SPGRID_INDEXSET_HH

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

    static const int dimension = Traits::dimension;

    template< int codim >
    struct Codim
    {
      typedef SPEntityInfo< Grid, codim > EntityInfo;
      typedef typename Traits::template Codim< codim >::Entity Entity;
    };

    typedef SPGridLevel< Grid > GridLevel;

  private:
    typedef typename GridLevel::MultiIndex MultiIndex;

    explicit SPIndexSet ( const GridLevel &gridLevel )
    {
      for( unsigned int codim = 0; codim <= dimension; ++codim )
        geomTypes_.push_back( GeometryType( GeometryType::cube, dimension-codim ) );
    }

    void update ( const GridLevel &gridLevel );

    IndexType index ( const MultiIndex &id ) const;

  public:
    template< class Entity >
    IndexType index ( const Entity &entity ) const
    {
      return index< Entity::codimension >( entity );
    }

    template< int codim >
    IndexType index ( const typename Codim< codim >::Entity &entity ) const
    {
      assert( contains( entity ) );
      const typename Codim< codim >::EntityInfo &entityInfo
        = Grid::getRealImplementation( entity ).entityInfo();
      return index( entityInfo.id() );
    }

    template< int codim >
    IndexType DUNE_DEPRECATED
    subIndex ( const typename Codim< 0 >::Entity &entity, const int i ) const
    {
      DUNE_THROW( NotImplemented, "SPIndexSet does not implement the old subIndex method." );
    }

    IndexType subIndex ( const typename Codim< 0 >::Entity &entity,
                         const int i, const unsigned int codim ) const
    {
      assert( contains( entity ) );
      const typename Codim< 0 >::EntityInfo &entityInfo
        = Grid::getRealImplementation( entity ).entityInfo();
      MultiIndex sid = entityInfo.id();
      sid += gridLevel().cube().subId( codim, i );
      return index( sid );
    }

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
      size_[ codim ];
    }

    template< class Entity >
    bool contains ( const Entity &entity ) const
    {
      return contains< Entity::codimension >( entity );
    }

    template< int codim >
    bool contains ( const typename Codim< codim >::Entity &entity )
    {
      const typename Codim< codim >::EntityInfo &entityInfo
        = Grid::getRealImplementation( entity ).entityInfo();
      return (&entity.gridLevel() == &gridLevel());
    }

    const GridLevel &gridLevel () const
    {
      return *gridLevel_;
    }

  private:
    const GridLevel *gridLevel_;
    IndexType offsets_[ 1 << dimension ];
    IndexType size_[ dimension+1 ];
    std::vector< GeometryType > geomTypes_[ dimension+1 ];
  };




  template< class Grid >
  void SPIndexSet< Grid >::update ( const GridLevel &gridLevel )
  {
    gridLevel_ = &gridLevel;

    const MultiIndex &cells = gridLevel.cells();
    for( unsigned int codim = 0; codim <= dimension; ++codim )
      size_[ codim ] = 0;

    for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
    {
      IndexType factor = 1;
      unsigned int codim = dimension;
      for( int j = 0; j < dimension; ++j )
      {
        const unsigned int d = (dir >> j) & 1;
        factor *= cells + (1-d);
        codim -= d;
      }
      offsets_[ dir ] = size_[ codim ];
      size_[ codim ] += factor;
    }
  }


  template< class Grid >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::index ( const MultiIndex &id ) const
  {
    const MultiIndex &cells = gridLevel().cells();
    IndexType index = 0;
    IndexType factor = 1;
    unsigned int dir = 0;
    for( int j = 0; j < dimension; ++j )
    {
      const unsigned int d = id[ j ] & 1;
      index += (id[ j ] >> 1) * factor;
      factor *= cells[ j ] + (1-d);
      dir |= (d << j);
    }
    return offsets_[ dir ] + index;
  }

}

#endif // #ifndef DUNE_SPGRID_INDEXSET_HH

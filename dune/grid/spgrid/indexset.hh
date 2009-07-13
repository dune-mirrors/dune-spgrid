#ifndef DUNE_SPGRID_INDEXSET_HH
#define DUNE_SPGRID_INDEXSET_HH

#include <dune/grid/common/indexset.hh>

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
      typedef typename Traits::template Codim< codim >::Entity Entity;
    };

    typedef SPGridLevel< ct, dim > GridLevel;

  private:
    SPIndexSet ( const GridLevel &gridLevel );

  public:
    using Base::index;

    template< int codim >
    IndexType index ( const typename Codim< codim >::Entity &entity ) const;

    template< int codim >
    IndexType DUNE_DEPRECATED
    subIndex ( const typename Codim< 0 >::Entity &entity, const int i ) const
    {
      DUNE_THROW( NotImplemented, "SPIndexSet does not implement the old subIndex method." );
    }

    IndexType subIndex ( const typename Codim< 0 >::Entity &entity,
                         const int i, const unsigned int codim ) const
    {
      // ...
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
      return entity.level() == gridLevel().level();
    }

    const GridLevel &gridLevel () const
    {
      return *gridLevel_;
    }

  private:
    const GridLevel *gridLevel_;
    IndexType size_[ dimension+1 ];
    std::vector< IndexType > offsets_[ dimension+1 ];
    std::vector< GeometryType > geomTypes_[ dimension+1 ];
  };


  template< class Grid >
  SPIndexSet< Grid >::SPIndexSet ( const GridLevel &gridLevel )
  : gridLevel_( &gridLevel )
  {
    for( int codim = 0; codim <= dimension; ++codim )
    {
      const unsigned int numDirections = gridLevel.numDirections( codim );
      offsets_.resize( numDirections );
      IndexType size = 0;
      for( int direction = 0; direction < numDirections; ++direction )
      {
        offsets_[ direction ] = size;
        IndexType factor = 1;
        for( int i = 0; i < dimension; ++i )
          factor *= gridLevel().n( direction, i );
        size += factor;
      }
      size_[ codim ] = size;

      geomTypes_.push_back( GeometryType( GeometryType::cube, dimension-codim ) );
    }
  }


  template< class Grid >
  template< int codim >
  typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::index ( const typename Codim< codim >::Entity &entity ) const
  {
    assert( contains( entity ) );
    const typename Entity::EntityInfo &entityInfo
      = Grid::getRealImplementation( entity ).entityInfo();

    const MultiIndex &multiIndex = entityInfo.multiIndex();
    const MultiIndex &multiDirection = entityInfo.multiDirection();
    const MultiIndex &n = gridLevel().n();
    IndexType index = 0;
    IndexType factor = 1;
    for( int i = 0; i < dimension; ++i )
    {
      index += multiIndex[ i ] * factor;
      factor *= n[ i ] + multiDirection[ i ];
    }
    const IndexType offset = offsets_[ codim ][ entityInfo.direction() ];
    return offset + index;
  }

}

#endif // #ifndef DUNE_SPGRID_INDEXSET_HH

#ifndef DUNE_SPGRID_IDSET_HH
#define DUNE_SPGRID_IDSET_HH

#include <dune/grid/common/indexidset.hh>

#include <dune/grid/spgrid/entityinfo.hh>
#include <dune/grid/spgrid/gridlevel.hh>

namespace Dune
{

  // SPLocalIdSet
  // ------------

  template< class Grid >
  class SPLocalIdSet
  : public IdSet< Grid, SPLocalIdSet< Grid >, unsigned int >
  {
    typedef SPLocalIdSet< Grid > This;
    typedef IdSet< Grid, This, unsigned int > Base;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Base::IdType IdType;

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
    typedef typename GridLevel::Mesh Mesh;

    IdType id ( const GridLevel &gridLevel, const MultiIndex &id ) const;

  public:
    template< class Entity >
    IdType id ( const Entity &entity ) const
    {
      return id< Entity::codimension >( entity );
    }

    template< int codim >
    IdType id ( const typename Codim< codim >::Entity &entity ) const
    {
      const typename Codim< codim >::EntityInfo &entityInfo
        = Grid::getRealImplementation( entity ).entityInfo();
      return id( entityInfo.gridLevel(), entityInfo.id() );
    }

    IdType subId ( const typename Codim< 0 >::Entity &entity, const int i, const unsigned int codim ) const
    {
      const typename Codim< 0 >::EntityInfo &entityInfo
        = Grid::getRealImplementation( entity ).entityInfo();
      const GridLevel &gridLevel = entityInfo.gridLevel();
      MultiIndex sid = entityInfo.id();
      sid +=  gridLevel.cube().subId( codim, i );
      return id( gridLevel, sid );
    }
  };


  template< class Grid >
  typename SPLocalIdSet< Grid >::IdType
  inline SPLocalIdSet< Grid >
    ::id ( const GridLevel &gridLevel, const MultiIndex &id ) const
  {
    const Mesh &globalMesh = gridLevel.globalMesh();
    const unsigned int level = gridLevel.level();
    
    IdType index = 0;
    IdType factor = 1;
    for( int i = 0; i < dimension; ++i )
    {
      index += id[ i ] * factor;
      factor *= 2*globalMesh.width( i ) + 1;
    }
    return index | (level << 26);
  }



  // SPGlobalIdSet
  // -------------

  template< class Grid >
  class SPGlobalIdSet
  : public SPLocalIdSet< Grid >
  {};

}

#endif // #ifndef DUNE_SPGRID_IDSET_HH

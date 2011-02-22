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
  : public IdSet< Grid, SPLocalIdSet< Grid >, unsigned long >
  {
    typedef SPLocalIdSet< Grid > This;
    typedef IdSet< Grid, This, unsigned long > Base;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Base::IdType IdType;

    static const int dimension = Traits::ReferenceCube::dimension;

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

    static const int levelShift = 8*sizeof( IdType ) - 8;

    IdType computeId ( const GridLevel &gridLevel, const MultiIndex &id ) const;

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
      return computeId( entityInfo.gridLevel(), entityInfo.id() );
    }

    IdType subId ( const typename Codim< 0 >::Entity &entity, const int i, const unsigned int codim ) const
    {
      const typename Codim< 0 >::EntityInfo &entityInfo
        = Grid::getRealImplementation( entity ).entityInfo();
      const GridLevel &gridLevel = entityInfo.gridLevel();
      MultiIndex sid = entityInfo.id();
      sid += gridLevel.referenceCube().subId( codim, i );
      return computeId( gridLevel, sid );
    }
  };


  template< class Grid >
  typename SPLocalIdSet< Grid >::IdType
  inline SPLocalIdSet< Grid >
    ::computeId ( const GridLevel &gridLevel, const MultiIndex &id ) const
  {
    const unsigned int level = gridLevel.level();
    if( (level > 0) && gridLevel.refinement().isCopy( id ) )
    {
      const Grid &grid = gridLevel.grid();
      MultiIndex fatherId( id );
      gridLevel.refinement().father( fatherId );
      return computeId( grid.gridLevel( level-1 ), fatherId );
    }

    const Mesh &globalMesh = gridLevel.globalMesh();

    IdType index = 0;
    IdType factor = 1;
    for( int i = 0; i < dimension; ++i )
    {
      index += IdType( id[ i ] ) * factor;
      factor *= IdType( 2*globalMesh.width( i ) + 1 );
    }
    return index | (IdType( level ) << levelShift);
  }



  // SPGlobalIdSet
  // -------------

  template< class Grid >
  class SPGlobalIdSet
  : public SPLocalIdSet< Grid >
  {};

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_IDSET_HH

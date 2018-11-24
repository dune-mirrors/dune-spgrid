#ifndef DUNE_SPGRID_IDSET_HH
#define DUNE_SPGRID_IDSET_HH

#include <type_traits>

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

    typedef typename std::remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Base::IdType IdType;

    static const int dimension = Traits::ReferenceCube::dimension;

    template< int codim >
    struct Codim
    {
      typedef __SPGrid::EntityInfo< Grid, codim > EntityInfo;
      typedef typename Traits::template Codim< codim >::Entity Entity;
    };

    typedef SPGridLevel< typename std::remove_const< Grid >::type > GridLevel;

  private:
    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename GridLevel::Mesh Mesh;

    static const int levelShift = 8*sizeof( IdType ) - 8;

    IdType computeId ( const GridLevel &gridLevel, const MultiIndex &id ) const;

    template< int cd >
    IdType computeSubId ( const GridLevel &gridLevel, const MultiIndex &id,
                          int i, int codim, std::integral_constant< int, cd > ) const;
    IdType computeSubId ( const GridLevel &gridLevel, const MultiIndex &id,
                          int i, int codim, std::integral_constant< int, 0 > ) const;
    IdType computeSubId ( const GridLevel &gridLevel, const MultiIndex &id,
                          int i, int codim, std::integral_constant< int, dimension > ) const;

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
        = entity.impl().entityInfo();
      return computeId( entityInfo.gridLevel(), entityInfo.id() );
    }

    template< class Entity >
    IdType subId ( const Entity &entity, int i, unsigned int codim ) const
    {
      return subId< Entity::codimension >( entity, i, codim );
    }

    template< int cd >
    IdType subId ( const typename Codim< cd >::Entity &entity, int i, unsigned int codim ) const
    {
      const typename Codim< cd >::EntityInfo &entityInfo
        = entity.impl().entityInfo();
      const GridLevel &gridLevel = entityInfo.gridLevel();
      return computeSubId( gridLevel, entityInfo.id(), i, codim, std::integral_constant< int, cd >() );
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


  template< class Grid >
  template< int cd >
  typename SPLocalIdSet< Grid >::IdType
  inline SPLocalIdSet< Grid >
    ::computeSubId ( const GridLevel &gridLevel, const MultiIndex &id,
                     int i, int codim, std::integral_constant< int, cd > ) const
  {
    const int mydim = dimension - cd;
    const SPMultiIndex< mydim > refId = gridLevel.template referenceCube< cd >().subId( codim - cd, i );
    MultiIndex subId( id );
    for( int k = 0, l = 0; k < dimension; ++k )
    {
      if( (id[ k ] & 1) != 0 )
        subId[ k ] += refId[ l++ ];
    }
    return computeId( gridLevel, subId );
  }

  template< class Grid >
  typename SPLocalIdSet< Grid >::IdType
  inline SPLocalIdSet< Grid >
    ::computeSubId ( const GridLevel &gridLevel, const MultiIndex &id,
                     int i, int codim, std::integral_constant< int, 0 > ) const
  {
    return computeId( gridLevel, id + gridLevel.referenceCube().subId( codim, i ) );
  }

  template< class Grid >
  typename SPLocalIdSet< Grid >::IdType
  inline SPLocalIdSet< Grid >
    ::computeSubId ( const GridLevel &gridLevel, const MultiIndex &id,
                     int i, int codim, std::integral_constant< int, dimension > ) const
  {
    assert( (codim == dimension) && (i == 0) );
    return computeId( gridLevel, id );
  }



  // SPGlobalIdSet
  // -------------

  template< class Grid >
  class SPGlobalIdSet
  : public SPLocalIdSet< Grid >
  {};

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_IDSET_HH

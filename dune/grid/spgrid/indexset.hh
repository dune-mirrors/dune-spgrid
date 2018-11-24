#ifndef DUNE_SPGRID_INDEXSET_HH
#define DUNE_SPGRID_INDEXSET_HH

#include <array>
#include <type_traits>
#include <vector>

#include <dune/grid/common/indexidset.hh>

#include <dune/grid/spgrid/entityinfo.hh>
#include <dune/grid/spgrid/gridlevel.hh>

namespace Dune
{

  template< class Grid >
  class SPIndexSet
    : public IndexSet< Grid, SPIndexSet< Grid >, unsigned int, std::array< GeometryType, 1 > >
  {
    typedef SPIndexSet< Grid > This;
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
    typedef typename GridLevel::PartitionList PartitionList;

  private:
    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename PartitionList::Partition Partition;

  public:
    SPIndexSet () = default;
    explicit SPIndexSet ( const GridLevel &gridLevel ) { update( gridLevel ); }

    void update ( const GridLevel &gridLevel );

  private:
    IndexType index ( const MultiIndex &id, unsigned int number ) const;

    template< int cd >
    IndexType subIndex ( const MultiIndex &id, int i, int codim, unsigned int number, std::integral_constant< int, cd > ) const;
    IndexType subIndex ( const MultiIndex &id, int i, int codim, unsigned int number, std::integral_constant< int, 0 > ) const;
    IndexType subIndex ( const MultiIndex &id, int i, int codim, unsigned int number, std::integral_constant< int, dimension > ) const;

  public:
    template< class Entity >
    IndexType index ( const Entity &entity ) const;

    template< int codim >
    IndexType index ( const typename Codim< codim >::Entity &entity ) const;

    template< class Entity >
    IndexType subIndex ( const Entity &entity, int i, unsigned int codim ) const;

    template< int cd >
    IndexType subIndex ( const typename Codim< cd >::Entity &entity, int i, unsigned int codim ) const;

    Types types ( int codim ) const { return {{ GeometryTypes::cube( dimension - codim ) }}; }

    IndexType size ( const GeometryType &type ) const;
    IndexType size ( const int codim ) const;

    template< class Entity >
    bool contains ( const Entity &entity ) const;

    template< int codim >
    bool contains ( const typename Codim< codim >::Entity &entity ) const;

    const GridLevel &gridLevel () const { return *gridLevel_; }

    const PartitionList &partitions () const { assert( partitions_ ); return *partitions_; }

  private:
    const GridLevel *gridLevel_ = nullptr;
    const PartitionList *partitions_ = nullptr;
    std::vector< std::array< IndexType, 1 << dimension > > offsets_;
    IndexType size_[ dimension+1 ];
  };



  // Implementation of SPIndexSet
  // ----------------------------

  template< class Grid >
  void SPIndexSet< Grid >::update ( const GridLevel &gridLevel )
  {
    gridLevel_ = &gridLevel;
    partitions_ = &gridLevel.template partition< All_Partition >();

    for( int codim = 0; codim <= dimension; ++codim )
      size_[ codim ] = 0;

    offsets_.resize( partitions().maxNumber() - partitions().minNumber() + 1 );
    for( typename PartitionList::Iterator pit = partitions().begin(); pit; ++pit )
    {
      for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
      {
        IndexType factor = 1;
        unsigned int codim = dimension;
        for( int j = 0; j < dimension; ++j )
        {
          const unsigned int d = (dir >> j) & 1;
          const int w = pit->bound( 1, j, d ) - pit->bound( 0, j, d );
          assert( w % 2 == 0 );
          factor *= (w / 2 + 1);
          codim -= d;
        }
        offsets_[ pit->number() - partitions().minNumber() ][ dir ] = size_[ codim ];
        size_[ codim ] += factor;
      }
    }
  }


  template< class Grid >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::index ( const MultiIndex &id, unsigned int number ) const
  {
    const Partition &partition = partitions().partition( number );

    IndexType index = 0;
    IndexType factor = 1;
    unsigned int dir = 0;
    for( int j = 0; j < dimension; ++j )
    {
      const unsigned int d = id[ j ] & 1;
      dir |= (d << j);

      const int begin = partition.bound( 0, j, d );
      const int end = partition.bound( 1, j, d );

      const IndexType idLocal = (id[ j ] - begin) >> 1;
      const IndexType width = ((end - begin) >> 1) + 1;
      assert( (idLocal >= 0) && (idLocal < width) );
      index += idLocal * factor;

      factor *= width;
    }
    return offsets_[ number - partitions().minNumber() ][ dir ] + index;
  }


  template< class Grid >
  template< int cd >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >
    ::subIndex ( const MultiIndex &id, int i, int codim, unsigned int number, std::integral_constant< int, cd > ) const
  {
    const int mydim = dimension - cd;
    const SPMultiIndex< mydim > refId = gridLevel().template referenceCube< cd >().subId( codim - cd, i );
    MultiIndex subId( id );
    for( int k = 0, l = 0; k < dimension; ++k )
    {
      if( (id[ k ] & 1) != 0 )
        subId[ k ] += refId[ l++ ];
    }
    return index( subId, number );
  }

  template< class Grid >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >
    ::subIndex ( const MultiIndex &id, int i, int codim, unsigned int number, std::integral_constant< int, 0 > ) const
  {
    return index( id + gridLevel().referenceCube().subId( codim, i ), number );
  }

  template< class Grid >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >
    ::subIndex ( const MultiIndex &id, int i, int codim, unsigned int number, std::integral_constant< int, dimension > ) const
  {
    assert( (codim == dimension) && (i == 0) );
    return index( id, number );
  }


  template< class Grid >
  template< class Entity >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::index ( const Entity &entity ) const
  {
    return index< Entity::codimension >( entity );
  }


  template< class Grid >
  template< int codim >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::index ( const typename Codim< codim >::Entity &entity ) const
  {
    assert( contains( entity ) );
    const typename Codim< codim >::EntityInfo &entityInfo
      = entity.impl().entityInfo();
    return index( entityInfo.id(), entityInfo.partitionNumber() );
  }


  template< class Grid >
  template< class Entity >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::subIndex ( const Entity &entity, int i, unsigned int codim ) const
  {
    return subIndex< Entity::codimension >( entity, i, codim );
  }


  template< class Grid >
  template< int cd >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >
    ::subIndex ( const typename Codim< cd >::Entity &entity, int i, unsigned int codim ) const
  {
    assert( contains( entity ) );
    const typename Codim< cd >::EntityInfo &entityInfo
      = entity.impl().entityInfo();
    // for the ghost approach, the partition number has to be corrected
    return subIndex( entityInfo.id(), i, codim, entityInfo.partitionNumber(), std::integral_constant< int, cd >() );
  }


  template< class Grid >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::size ( const GeometryType &type ) const
  {
    return (type.isCube() ? size( dimension - type.dim() ) : 0);
  }


  template< class Grid >
  inline typename SPIndexSet< Grid >::IndexType
  SPIndexSet< Grid >::size ( const int codim ) const
  {
    assert( (codim >= 0) && (codim <= dimension) );
    return size_[ codim ];
  }


  template< class Grid >
  template< class Entity >
  inline bool SPIndexSet< Grid >::contains ( const Entity &entity ) const
  {
    return contains< Entity::codimension >( entity );
  }


  template< class Grid >
  template< int codim >
  inline bool SPIndexSet< Grid >
    ::contains ( const typename Codim< codim >::Entity &entity ) const
  {
    const typename Codim< codim >::EntityInfo &entityInfo
      = entity.impl().entityInfo();
    assert( partitions().contains( entityInfo.partitionNumber() ) );
    return (&entityInfo.gridLevel() == &gridLevel());
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_INDEXSET_HH

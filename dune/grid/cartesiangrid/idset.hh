#ifndef DUNE_CARTESIANGRID_IDSET_HH
#define DUNE_CARTESIANGRID_IDSET_HH

#include <dune/grid/common/indexidset.hh>

namespace Dune
{

  // CartesianGridIdSet
  // -----------
  
  template< class Grid, class HostIdSet >
  class CartesianGridIdSet
  : public IdSet< Grid, CartesianGridIdSet< Grid, HostIdSet >, typename HostIdSet::IdType >
  {
    typedef IdSet< Grid, CartesianGridIdSet< Grid, HostIdSet >, typename HostIdSet::IdType > Base;

    typedef typename remove_const< Grid >::type::Traits Traits;
    
  public:
    typedef typename HostIdSet::IdType IdType;

    using Base::subId;

    CartesianGridIdSet ( const HostIdSet &hostIdSet )
    : hostIdSet_( hostIdSet )
    {}

    //! id meethod for entity and specific codim 
    template< int codim >
    IdType id ( const typename Traits::template Codim< codim >::Entity &entity ) const
    {
      return id( Grid::template getHostEntity< codim >( entity ) );
    }

    //! id method for host entity (e.g. in ParallelGrid)
    template< int codim >
    IdType id ( const typename Traits::HostGridType::template Codim< codim >::Entity &entity ) const
    {
      return hostIdSet_.id( entity );
    }

    //! id method of all entities 
    template< class Entity >
    IdType id ( const Entity &entity ) const
    {
      return id< Entity::codimension >( entity );
    }

    IdType subId ( const typename Traits::template Codim< 0 >::Entity &entity, int i, unsigned int codim ) const
    {
      return hostIdSet_.subId( Grid::template getHostEntity< 0 >( entity ), i, codim );
    }

  private:
    CartesianGridIdSet ( const CartesianGridIdSet & );
    CartesianGridIdSet &operator= ( const CartesianGridIdSet & );

    const HostIdSet &hostIdSet_;
  };

}

#endif // #ifndef DUNE_CARTESIANGRID_IDSET_HH

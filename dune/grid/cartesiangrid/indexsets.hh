#ifndef DUNE_CARTESIANGRID_INDEXSETS_HH
#define DUNE_CARTESIANGRID_INDEXSETS_HH

#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/grid/cartesiangrid/capabilities.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class HostGrid >
  class CartesianGrid;



  // CartesianGridIndexSet
  // --------------

  template< class Grid, class HostIndexSet >
  class CartesianGridIndexSet
  : public IndexSet< Grid, CartesianGridIndexSet< Grid, HostIndexSet >, typename HostIndexSet::IndexType >
  {
    typedef CartesianGridIndexSet< Grid, HostIndexSet > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::HostGrid HostGrid;

  public:
    typedef IndexSet< Grid, This, typename HostIndexSet::IndexType > Base;

    static const int dimension = Grid::dimension;

    typedef typename Base::IndexType IndexType;

    CartesianGridIndexSet ( const HostIndexSet &hostIndexSet )
    : hostIndexSet_( &hostIndexSet )
    {}

    using Base::index;
    using Base::subIndex;

    template< int cc >
    IndexType index ( const typename Traits::template Codim< cc >::Entity &entity ) const
    {
      return hostIndexSet().index( Grid::template getHostEntity< cc >( entity ) );
    }

    template< int cc >
    IndexType subIndex ( const typename Traits::template Codim< cc >::Entity &entity, int i, unsigned int codim ) const
    {
      return hostIndexSet().subIndex( Grid::template getHostEntity< cc >( entity ), i, codim );
    }

    IndexType size ( GeometryType type ) const
    {
      return hostIndexSet().size( type );
    }
         
    int size ( int codim ) const
    {
      return hostIndexSet().size( codim );
    }

    template< class Entity >
    bool contains ( const Entity &entity ) const
    {
      static const int cc = Entity::codimension;
      return hostIndexSet().contains( Grid::template getHostEntity< cc >( entity ) );
    }

    const std::vector< GeometryType > &geomTypes ( int codim ) const
    {
      return hostIndexSet().geomTypes( codim );
    }
   
  private:
    const HostIndexSet &hostIndexSet () const
    {
      assert( hostIndexSet_ );
      return *hostIndexSet_;
    }

    const HostIndexSet *hostIndexSet_;
  };

}

#endif // #ifndef DUNE_CARTESIANGRID_INDEXSETS_HH

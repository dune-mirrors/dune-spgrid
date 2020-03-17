#ifndef DUNE_SPGRID_CHECKIDCOMMUNICATION_HH
#define DUNE_SPGRID_CHECKIDCOMMUNICATION_HH

#include <dune/common/hybridutilities.hh>

#include <dune/geometry/dimension.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

namespace Dune
{

  template< class VT >
  struct CheckIdCommunicationDataHandle;


  template< InterfaceType iftype, class VT >
  inline void checkIdCommunication ( const GridView< VT > &gridView )
  {
    std::cout << "Checking communication ids for " << iftype << "..." << std::endl;
    CheckIdCommunicationDataHandle< VT > handle( gridView );
    gridView.communicate( handle, iftype, Dune::ForwardCommunication );
  }


  template< class VT >
  inline void checkIdCommunication ( const GridView< VT > &gridView )
  {
    checkIdCommunication< InteriorBorder_InteriorBorder_Interface >( gridView );
    checkIdCommunication< InteriorBorder_All_Interface >( gridView );
    checkIdCommunication< Overlap_OverlapFront_Interface >( gridView );
    checkIdCommunication< Overlap_All_Interface >( gridView );
    checkIdCommunication< All_All_Interface >( gridView );
  }


  template< class VT >
  struct CheckIdCommunicationDataHandle
  : public CommDataHandleIF< CheckIdCommunicationDataHandle< VT >, typename GridView< VT >::Grid::GlobalIdSet::IdType >
  {
    typedef Dune::GridView< VT > GridView;

    static const int dimension = GridView::dimension;

    typedef typename GridView::Grid Grid;
    typedef typename Grid::GlobalIdSet GlobalIdSet;
    typedef typename GlobalIdSet::IdType IdType;

    explicit CheckIdCommunicationDataHandle ( const GridView &gridView )
      : rank_( gridView.grid().comm().rank() ),
        idSet_( gridView.grid().globalIdSet() )
    {
      Hybrid::forEach( std::make_integer_sequence< int, dimension+1 >(), [ this ] ( auto codim ) {
          contains_[ codim ] = Dune::Capabilities::canCommunicate< Grid, codim >::v;
        } );
    }

    bool contains ( int const dim, const int codim ) const
    {
      return ((codim >= 0) && (codim <= dimension) ? contains_[ codim ] : false);
    }

    bool fixedSize ( const int dim, const int codim ) const
    {
      return true;
    }

    template< class Entity >
    size_t size ( const Entity &entity ) const
    {
      return (contains_[ Entity::codimension ] ? 1 : 0);
    }

    template< class Buffer, class Entity >
    void gather ( Buffer &buffer, const Entity &entity ) const
    {
      buffer.write( idSet_.id( entity ) );
    }

    template< class Buffer, class Entity >
    void scatter ( Buffer &buffer, const Entity &entity, size_t n )
    {
      IdType id;
      buffer.read( id );

      if( id != idSet_.id( entity ) )
      {
        std::cerr << "[ " << rank_ << " ] Error: receive data from entity " << id
                  << " on entity " << idSet_.id( entity ) << "." << std::endl;
      }
    }

  private:
    const int rank_;
    const GlobalIdSet &idSet_;
    bool contains_[ dimension+1 ];
  };

}

#endif // #ifndef DUNE_SPGRID_CHECKIDCOMMUNICATION_HH

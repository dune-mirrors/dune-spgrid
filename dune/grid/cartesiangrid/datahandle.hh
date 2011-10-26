#ifndef DUNE_CARTESIANGRID_DATAHANDLE_HH
#define DUNE_CARTESIANGRID_DATAHANDLE_HH

#include <dune/common/typetraits.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid/cartesiangrid/capabilities.hh>

namespace Dune
{
  template< class Grid, class HostEntity >
  struct CartesianGridEntityProxy
  {
    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::ExtraData ExtraData;

    static const int dimension = HostEntity::dimension;
    static const int codimension = HostEntity::codimension;
    
    typedef Dune::Entity< codimension, dimension, const Grid, CartesianGridEntity > Entity;

  private:
    typedef CartesianGridEntity< codimension, dimension, const Grid > EntityImpl;

  public:
    CartesianGridEntityProxy ( ExtraData data, const HostEntity &hostEntity )
    : entity_( EntityImpl( data, hostEntity ) )
    {}

    const Entity &operator* () const { return entity_; }
    Entity &operator* () { return entity_; }

  private:
    Entity entity_;
  };


  // CartesianGridDataHandle
  // -----------------------

  template< class Grid, class WrappedHandle >
  class CartesianGridDataHandle
  : public CommDataHandleIF< CartesianGridDataHandle< Grid, WrappedHandle >, typename WrappedHandle::DataType >
  {
    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::ExtraData ExtraData;

  public:
    CartesianGridDataHandle ( ExtraData data, WrappedHandle &handle )
    : data_( data ),
      wrappedHandle_( handle )
    {}

    bool contains ( int dim, int codim ) const
    {
      return wrappedHandle_.contains( dim, codim );
    }

    bool fixedsize ( int dim, int codim ) const
    {
      return wrappedHandle_.fixedsize( dim, codim );
    }

    template< class HostEntity >
    size_t size ( const HostEntity &hostEntity ) const
    {
      CartesianGridEntityProxy< Grid, HostEntity > proxy( data_, hostEntity );
      return wrappedHandle_.size( *proxy );
    }

    template< class MessageBuffer, class HostEntity >
    void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
    {
      CartesianGridEntityProxy< Grid, HostEntity > proxy( data_, hostEntity );
      wrappedHandle_.gather( buffer, *proxy );
    }

    template< class MessageBuffer, class HostEntity >
    void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
    {
      CartesianGridEntityProxy< Grid, HostEntity > proxy( data_, hostEntity );
      wrappedHandle_.scatter( buffer, *proxy, size );
    }

  private:
    ExtraData data_;
    WrappedHandle &wrappedHandle_;
  };

  // CartesianGridWrappedDofManager
  // ----------------

  template< class Grid, class WrappedHandle >
  class CartesianGridWrappedDofManager
  {
    typedef typename remove_const< Grid >::type::Traits Traits;
    typedef typename Traits::ExtraData ExtraData;

    typedef typename Grid :: HostGrid HostGrid ;
    typedef typename Grid :: InStreamType   InStreamType ;
    typedef typename Grid :: OutStreamType  OutStreamType ;
    typedef typename HostGrid :: template Codim< 0 > :: Entity  HostElement;

  public:
    CartesianGridWrappedDofManager( ExtraData data, WrappedHandle &handle )
    : data_( data ),
      wrappedHandle_( handle )
    {}

    void inlineData ( OutStreamType &stream, const HostElement &hostElement ) const
    {
      const CartesianGridEntityProxy< Grid, HostElement > proxy( data_, hostElement );
      wrappedHandle_.inlineData( stream, *proxy );
    }

    void xtractData ( InStreamType &stream, const HostElement &hostElement, size_t newElements )
    {
      const CartesianGridEntityProxy< Grid, HostElement > proxy( data_, hostElement );
      wrappedHandle_.xtractData( stream, *proxy, newElements );
    }

    void compress ()
    {
      wrappedHandle_.compress();
    }

  private:
    ExtraData data_;
    WrappedHandle &wrappedHandle_;
  };

}

#endif // #ifndef DUNE_CARTESIANGRID_DATAHANDLE_HH

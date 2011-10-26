#ifndef DUNE_CARTESIANGRID_ADAPTCALLBACK_HH
#define DUNE_CARTESIANGRID_ADAPTCALLBACK_HH

#include <dune/grid/common/adaptcallback.hh>
#include <dune/grid/cartesiangrid/datahandle.hh>

namespace Dune
{

  // CartesianGridAdaptDataHandle
  // ------------------------

  template< class Grid, class WrappedHandle >
  class CartesianGridAdaptDataHandle 
    : public AdaptDataHandle< typename Grid::HostGrid, CartesianGridAdaptDataHandle< Grid, WrappedHandle > >
  {
    typedef CartesianGridAdaptDataHandle< Grid, WrappedHandle > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::ExtraData ExtraData;

  public:
    typedef typename Grid::template Codim< 0 >::Entity Entity;

  private:
    CartesianGridAdaptDataHandle ( const This & );
    This &operator= ( const This & );

  public:
    CartesianGridAdaptDataHandle( ExtraData data,  WrappedHandle &handle )
      : data_( data ),
        wrappedHandle_( handle )
    {}

    void preAdapt ( const unsigned int estimateAdditionalElements )
    {
      wrappedHandle_.preAdapt( estimateAdditionalElements );
    }

    void postAdapt ()
    {
      wrappedHandle_.postAdapt();
    }

    template <class HostEntity> 
    void preCoarsening ( const HostEntity &father ) const
    {
      CartesianGridEntityProxy< Grid, HostEntity > proxy( data_, father );
      wrappedHandle_.preCoarsening( *proxy );
    }

    template <class HostEntity> 
    void postRefinement ( const HostEntity &father ) const
    {
      CartesianGridEntityProxy< Grid, HostEntity > proxy( data_, father );
      wrappedHandle_.postRefinement( *proxy );
    }

  protected:
    ExtraData data_;
    WrappedHandle &wrappedHandle_;
  };

} // end namespace Dune 

#endif

#ifndef DUNE_SPGRID_IDSET_HH
#define DUNE_SPGRID_IDSET_HH

#include <dune/grid/common/indexidset.hh>

namespace Dune
{

  template< class Grid >
  class SPLocalIdSet
  : public IdSet< Grid, SPLocalIdSet< Grid >, unsigned int >
  {
    typedef SPLocalIdSet< Grid > This;
    typedef IdSet< Grid, This, unsigned int > Base;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Base::IdType IdType;

    template< int codim >
    struct Codim
    {
      typedef typename Traits::template Codim< codim >::Entity Entity;
    }

    using Base::id;

    template< int codim >
    IdType id ( const typename Codim< codim >::Entity &entity ) const
    {
      // ...
    }

    template< int codim >
    IdType DUNE_DEPRECATED
    subId ( const typename Codim< 0 >::Entity &entity, const int i ) const
    {
      DUNE_THROW( NotImplemented, "SPLocalIdSet does not implement the old subId method." );
    }

    IdType subId ( const typename Codim< 0 >::Entity &entity, const int i, const unsigned int codim ) const
    {
      // ...
    }
  };


  template< class Grid >
  class SPGlobalIdSet
  : public SPLocalIdSet< Grid >
  {};

}

#endif // #ifndef DUNE_SPGRID_IDSET_HH

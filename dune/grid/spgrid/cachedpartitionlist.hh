#ifndef DUNE_SPGRID_CACHEDPARTITIONLIST_HH
#define DUNE_SPGRID_CACHEDPARTITIONLIST_HH

#include <limits>

#include <dune/grid/spgrid/partitionlist.hh>

namespace Dune
{

  // SPCachedPartitionList
  // ---------------------

  template< int dim >
  class SPCachedPartitionList
  : public SPPartitionList< dim >
  {
    typedef SPCachedPartitionList< dim > This;
    typedef SPPartitionList< dim > Base;

  protected:
    using Base::Node;

  public:
    SPCachedPartitionList ()
    : first_( std::numeric_limits< unsigned int >::max() ),
      last_( std::numeric_limits< unsigned int >::min() ),
      cache_( 0 )
    {}

    SPCachedPartitionList ( const This &other )
    : Base( other )
      cache_( 0 )
    {
      updateCache();
    }

    void updateCache ();

  private:
    unsigned int first_, last_;
    const Node **cache_;
  };



  // Implementation of SPCachedPartitionList
  // ---------------------------------------

  template< int dim >
  inline void SPCachedPartitionList< dim >::updateCache ()
  {
    // find new cache size
    first_ = std::numeric_limits< unsigned int >::max();
    last_ = std::numeric_limits< unsigned int >::min();
    for( const Node *it = Base::head_; it; it = it->next() )
    {
      first_ = std::min( first_, it->partition().number() );
      last_ = std::max( last_, it->partition().number() );
    }

    // create new empty cache
    delete[] cache_;
    const unsigned int cacheSize = last_ - first_ + 1;
    cache_ = new const Node *[ cacheSize ];
    for( unsigned int i = 0; i < cacheSize; ++i )
      cache_[ i ] = 0;

    // fill cache
    for( const Node *it = Base::head_; it; it = it->next() )
    {
      const unsigned int number = it->partition().number();
      if( cache_[ number - first ] != 0 )
        DUNE_THROW( GridError, "Partition number " << number << " is not unique." );
      cache_[ number - first_ ] = it;
    }
  }

}

#endif // #ifndef DUNE_SPGRID_CACHEDPARTITIONLIST_HH

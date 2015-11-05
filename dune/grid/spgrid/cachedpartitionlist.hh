#ifndef DUNE_SPGRID_CACHEDPARTITIONLIST_HH
#define DUNE_SPGRID_CACHEDPARTITIONLIST_HH

#include <limits>

#include <dune/grid/common/exceptions.hh>

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
    typedef typename Base::Node Node;

  public:
    typedef typename Base::MultiIndex MultiIndex;
    typedef typename Base::Partition Partition;

    SPCachedPartitionList ()
    : first_( std::numeric_limits< unsigned int >::max() ),
      last_( std::numeric_limits< unsigned int >::min() ),
      cache_( nullptr )
    {}

    SPCachedPartitionList ( const This &other )
    : Base( other ),
      cache_( nullptr )
    {
      updateCache();
    }

    ~SPCachedPartitionList () { delete[] cache_; }

    const This &operator= ( const This &other )
    {
      *(Base *)this = other;
      updateCache();
      return *this;
    }

    bool contains ( const unsigned int number ) const;
    bool contains ( const MultiIndex &id, const unsigned int number ) const;
    const Partition &partition ( const unsigned int number ) const;

    unsigned int minNumber () const;
    unsigned int maxNumber () const;

    void updateCache ();

  private:
    unsigned int first_, last_;
    const Node **cache_;
  };



  // Implementation of SPCachedPartitionList
  // ---------------------------------------

  template< int dim >
  inline bool
  SPCachedPartitionList< dim >::contains ( const unsigned int number ) const
  {
    const bool inRange = (number >= minNumber()) && (number <= maxNumber());
    return inRange && (cache_[ number - minNumber() ] != 0);
  }


  template< int dim >
  inline bool
  SPCachedPartitionList< dim >
    ::contains ( const MultiIndex &id, const unsigned int number ) const
  {
    if( contains( number ) )
      return cache_[ number - minNumber() ]->partition().contains( id );
    else
      return false;
  }


  template< int dim >
  inline const typename SPCachedPartitionList< dim >::Partition &
  SPCachedPartitionList< dim >::partition ( const unsigned int number ) const
  {
    assert( contains( number ) );
    return cache_[ number - minNumber() ]->partition();
  }


  template< int dim >
  inline unsigned int SPCachedPartitionList< dim >::minNumber () const
  {
    return first_;
  }


  template< int dim >
  inline unsigned int SPCachedPartitionList< dim >::maxNumber () const
  {
    return last_;
  }


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
      if( cache_[ number - first_ ] != 0 )
        DUNE_THROW( GridError, "Partition number " << number << " is not unique." );
      cache_[ number - first_ ] = it;
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_CACHEDPARTITIONLIST_HH

#ifndef DUNE_SPGRID_ITERATOR_HH
#define DUNE_SPGRID_ITERATOR_HH

#include <algorithm>
#include <type_traits>

#include <dune/grid/spgrid/direction.hh>
#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/entity.hh>

namespace Dune
{

  // SPPartitionIterator
  // -------------------

  template< int codim, class Grid >
  class SPPartitionIterator
  {
    typedef SPPartitionIterator< codim, Grid > This;

  public:
    typedef typename std::remove_const< Grid >::type::Traits Traits;

    static const int dimension = Traits::ReferenceCube::dimension;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;

    typedef typename Traits::template Codim< codimension >::Entity Entity;

  private:
    typedef SPEntity< codimension, dimension, Grid > EntityImpl;

  public:
    typedef typename EntityImpl::EntityInfo EntityInfo;
    typedef typename EntityImpl::GridLevel GridLevel;

    typedef SPPartitionList< dimension > PartitionList;

    typedef typename EntityInfo::Direction Direction;

    static const unsigned int numDirections = GridLevel::numDirections;

    struct Begin {};
    struct End {};

  protected:
    typedef typename EntityInfo::MultiIndex MultiIndex;

    typedef SPDirectionIterator< dimension, codimension > DirectionIterator;

  public:
    SPPartitionIterator () = default;

    SPPartitionIterator ( const GridLevel &gridLevel, const PartitionList &partitionList,
                          const Begin &b, const unsigned int sweepDir = 0 );
    SPPartitionIterator ( const GridLevel &gridLevel, const PartitionList &partitionList,
                          const End &e, const unsigned int sweepDir = 0 );

    operator bool () const { return bool( partition_ ); }

    Entity operator* () const { return dereference(); }

    bool operator== ( const This &other ) const { return equals( other ); }
    bool operator!= ( const This &other ) const { return !equals( other ); }

    This &operator++ () { increment(); return *this; }

    Entity dereference () const { return EntityImpl( entityInfo() ); }

    bool equals ( const This &other ) const { return entityInfo().equals( other.entityInfo() ); }

    void increment ();

    const EntityInfo &entityInfo () const { return entityInfo_; }
    EntityInfo &entityInfo () { return entityInfo_; }

    const GridLevel &gridLevel () const { return entityInfo().gridLevel(); }

  private:
    int begin ( int i, Direction dir ) const;
    int end ( int i, Direction dir ) const;

    void init ();

  private:
    EntityInfo entityInfo_;
    typename PartitionList::Iterator partition_;
    unsigned int sweepDirection_;
  };



  // Implementation of SPPartitionIterator
  // -------------------------------------

  template< int codim, class Grid >
  inline SPPartitionIterator< codim, Grid >
    ::SPPartitionIterator ( const GridLevel &gridLevel, const PartitionList &partitionList,
                            const Begin &b, unsigned int sweepDir )
    : entityInfo_( gridLevel ),
      partition_( partitionList.begin() ),
      sweepDirection_( sweepDir )
  {
    assert( sweepDir < numDirections );
    init();
  }


  template< int codim, class Grid >
  inline SPPartitionIterator< codim, Grid >
    ::SPPartitionIterator ( const GridLevel &gridLevel, const PartitionList &partitionList,
                            const End &e, unsigned int sweepDir )
    : entityInfo_( gridLevel ),
      partition_( partitionList.end() ),
      sweepDirection_( sweepDir )
  {
    assert( sweepDir < numDirections );
    init();
  }


  template< int codim, class Grid >
  inline void SPPartitionIterator< codim, Grid >::increment ()
  {
    MultiIndex &id = entityInfo().id();
    for( int i = 0; i < dimension; ++i )
    {
      const unsigned int sweep = (sweepDirection_ >> i) & 1;
      id[ i ] += (2 - 4*sweep);
      if( id[ i ] != end( i, entityInfo().direction() ) )
        return entityInfo().update();
      id[ i ] = begin( i, entityInfo().direction() );
    }

    DirectionIterator dirIt( entityInfo().direction() );
    ++dirIt;
    for( ; dirIt && partition_->empty( *dirIt ); ++dirIt )
      continue;
    if( dirIt )
    {
      for( int i = 0; i < dimension; ++i )
        id[ i ] = begin( i, *dirIt );
      entityInfo().update();
    }
    else
    {
      ++partition_;
      init();
    }
  }


  template< int codim, class Grid >
  inline int SPPartitionIterator< codim, Grid >::begin ( int i, Direction dir ) const
  {
    const unsigned int s = (sweepDirection_ >> i) & 1;
    return partition_->bound( s, i, dir[ i ] );
  }


  template< int codim, class Grid >
  inline int SPPartitionIterator< codim, Grid >::end ( int i, Direction dir ) const
  {
    const unsigned int s = (sweepDirection_ >> i) & 1;
    const int bnd = partition_->bound( 1-s, i, dir[ i ] );
    return bnd + 2*(2*(1-s) - 1);
  }


  template< int codim, class Grid >
  inline void SPPartitionIterator< codim, Grid >::init ()
  {
    MultiIndex &id = entityInfo().id();
    if( partition_ )
    {
      DirectionIterator dirIt;
      for( ; dirIt && partition_->empty( *dirIt ); ++dirIt )
        continue;
      if( dirIt )
      {
        for( int i = 0; i < dimension; ++i )
          id[ i ] = begin( i, *dirIt );
        entityInfo().update( partition_->number() );
      }
      else
      {
        ++partition_;
        init();
      }
    }
    else
      std::fill( id.begin(), id.end(), std::numeric_limits< int >::max() );
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_ITERATOR_HH

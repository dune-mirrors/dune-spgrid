#ifndef DUNE_SPGRID_ITERATOR_HH
#define DUNE_SPGRID_ITERATOR_HH

#include <algorithm>

#include <dune/grid/spgrid/direction.hh>
#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/entitypointer.hh>

namespace Dune
{

  // SPPartitionIterator
  // -------------------

  template< int codim, class Grid >
  class SPPartitionIterator
    : public SPEntityPointer< codim, Grid >
  {
    typedef SPPartitionIterator< codim, Grid > This;
    typedef SPEntityPointer< codim, Grid > Base;

  public:
    typedef typename Base::EntityInfo EntityInfo;
    typedef typename Base::GridLevel GridLevel;

    static const int dimension = Base::dimension;
    static const int codimension = Base::codimension;
    static const int mydimension = Base::mydimension;

    typedef SPPartitionList< dimension > PartitionList;

    typedef typename EntityInfo::Direction Direction;

    static const unsigned int numDirections = GridLevel::numDirections;

    struct Begin {};
    struct End {};

  protected:
    typedef typename EntityInfo::MultiIndex MultiIndex;

    typedef SPDirectionIterator< dimension, codimension > DirectionIterator;

  public:
    using Base::entityInfo;
    using Base::gridLevel;

    SPPartitionIterator () = default;

    SPPartitionIterator ( const GridLevel &gridLevel, const PartitionList &partitionList,
                          const Begin &b, const unsigned int sweepDir = 0 );
    SPPartitionIterator ( const GridLevel &gridLevel, const PartitionList &partitionList,
                          const End &e, const unsigned int sweepDir = 0 );

    operator bool () const { return bool( partition_ ); }
    This &operator++ ();

    void increment ();

  private:
    int begin ( int i, Direction dir ) const;
    int end ( int i, Direction dir ) const;

    void init ();

  private:
    typename PartitionList::Iterator partition_;
    unsigned int sweepDirection_;
  };



  // Implementation of SPPartitionIterator
  // -------------------------------------

  template< int codim, class Grid >
  inline SPPartitionIterator< codim, Grid >
    ::SPPartitionIterator ( const GridLevel &gridLevel, const PartitionList &partitionList,
                            const Begin &b, const unsigned int sweepDir )
  : Base( gridLevel ),
    partition_( partitionList.begin() ),
    sweepDirection_( sweepDir )
  {
    assert( sweepDir < numDirections );
    init();
  }


  template< int codim, class Grid >
  inline SPPartitionIterator< codim, Grid >
    ::SPPartitionIterator ( const GridLevel &gridLevel, const PartitionList &partitionList,
                            const End &e, const unsigned int sweepDir )
  : Base( gridLevel ),
    partition_( partitionList.end() ),
    sweepDirection_( sweepDir )
  {
    assert( sweepDir < numDirections );
    init();
  }


  template< int codim, class Grid >
  inline typename SPPartitionIterator< codim, Grid >::This &
  SPPartitionIterator< codim, Grid >::operator++ ()
  {
    increment();
    return *this;
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

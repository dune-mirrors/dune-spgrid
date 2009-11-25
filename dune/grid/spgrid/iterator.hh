#ifndef DUNE_SPGRID_ITERATOR_HH
#define DUNE_SPGRID_ITERATOR_HH

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

    typedef typename GridLevel::Partition Partition;

    static const int dimension = Base::dimension;
    static const int codimension = Base::codimension;
    static const int mydimension = Base::mydimension;

    static const unsigned int numDirections = GridLevel::numDirections;

  protected:
    typedef typename EntityInfo::MultiIndex MultiIndex;

  public:
    using Base::gridLevel;

    SPPartitionIterator ( const GridLevel &gridLevel, const Partition *partition, const unsigned int sweepDir = 0 );

    void increment ();

  private:
    int begin ( const int i, const unsigned int dir ) const;
    int end ( const int i, const unsigned int dir ) const;

    void init ( const Partition *partition );

  protected:
    using Base::entity_;

  private:
    const Partition *partition_;
    unsigned int sweepDirection_;
  };



  // Implementation of SPPartitionIterator
  // -------------------------------------

  template< int codim, class Grid >
  inline SPPartitionIterator< codim, Grid >
    ::SPPartitionIterator ( const GridLevel &gridLevel, const Partition *partition, const unsigned int sweepDir )
  : Base( gridLevel ),
    sweepDirection_( sweepDir )
  {
    assert( sweepDir < numDirections );

    init( partition );
  }


  template< int codim, class Grid >
  inline void SPPartitionIterator< codim, Grid >::increment ()
  {
    assert( partition_ != 0 );

    EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();
    MultiIndex &id = entityInfo.id();
    for( int i = 0; i < dimension; ++i )
    {
      const unsigned int sweep = (sweepDirection_ >> i) & 1;
      id[ i ] += (2 - 4*sweep);
      if( id[ i ] != end( i, entityInfo.direction() ) )
        return entityInfo.update();
      id[ i ] = begin( i, entityInfo.direction() );
    }

    unsigned int dir = entityInfo.direction()+1;
    const unsigned int mydim = mydimension;
    for( ; (dir < numDirections) && (bitCount( dir ) != mydim); ++dir );
    if( dir < numDirections )
    {
      for( int i = 0; i < dimension; ++i )
        id[ i ] = begin( i, dir );
      entityInfo.update();
    }
    else
      init( partition_->next() );
  }


  template< int codim, class Grid >
  inline int
  SPPartitionIterator< codim, Grid >::begin ( const int i, const unsigned int dir ) const
  {
    assert( partition_ != 0 );
    //const MultiIndex &cells = gridLevel().cells();
    const MultiIndex &begin = partition_->begin();
    const MultiIndex &end = partition_->end();
    const unsigned int sweep = (sweepDirection_ >> i) & 1;
    const unsigned int d = (dir >> i) & 1;
    //return d + sweep*2*(cells[ i ]-d);
    return (1-sweep)*(2*begin[ i ] + d) + sweep*(2*end[ i ] - d);
  }


  template< int codim, class Grid >
  inline int
  SPPartitionIterator< codim, Grid >::end ( const int i, const unsigned int dir ) const
  {
    assert( partition_ != 0 );
    //const MultiIndex &cells = gridLevel().cells();
    const MultiIndex &begin = partition_->begin();
    const MultiIndex &end = partition_->end();
    const unsigned int sweep = (sweepDirection_ >> i) & 1;
    const unsigned int d = (dir >> i) & 1;
    //return (d-2) + (1-sweep)*2*(cells[ i ]-(d-2));
    return (1-sweep)*(2*end[ i ] - (d-2)) + sweep*(2*begin[ i ] + (d-2));
  }


  template< int codim, class Grid >
  inline void
  SPPartitionIterator< codim, Grid >::init ( const Partition *partition )
  {
    partition_ = partition;

    EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();
    MultiIndex &id = entityInfo.id();
    if( partition_ != 0 )
    {
      unsigned int dir = 0;
      const unsigned int mydim = mydimension;
      for( ; (dir < numDirections) && (bitCount( dir ) != mydim); ++dir );
      assert( dir < numDirections );

      for( int i = 0; i < dimension; ++i )
        id[ i ] = begin( i, dir );
    }
    else
      id = std::numeric_limits< MultiIndex >::max();

    entityInfo.update();
  }



  // SPIterator
  // ----------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class SPIterator
  : public SPPartitionIterator< codim, Grid >
  {
    typedef SPIterator< codim, pitype, Grid > This;
    typedef SPPartitionIterator< codim, Grid > Base;

  public:
    typedef typename Base::GridLevel GridLevel;

    struct Begin {};
    struct End {};

    SPIterator ( const GridLevel &gridLevel, const Begin &b, const unsigned int sweepDir = 0 )
    : Base( gridLevel, &gridLevel.allPartition(), sweepDir )
    {}

    SPIterator ( const GridLevel &gridLevel, const End &e, const unsigned int sweepDir = 0 )
    : Base( gridLevel, 0, sweepDir )
    {}
  };

}

#endif // #ifndef DUNE_SPGRID_ITERATOR_HH

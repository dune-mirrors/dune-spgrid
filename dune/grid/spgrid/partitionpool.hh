#ifndef DUNE_SPGRID_PARTITIONPOOL_HH
#define DUNE_SPGRID_PARTITIONPOOL_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/exceptions.hh>

#include <dune/grid/spgrid/cachedpartitionlist.hh>

namespace Dune
{

  // SPPartitionPool
  // ---------------

  template< int dim >
  class SPPartitionPool
  {
    typedef SPPartitionPool< dim > This;

  public:
    static const int dimension = dim;

    typedef SPCachedPartitionList< dimension > PartitionList;

    typedef typename PartitionList::Partition Partition;
    typedef typename PartitionList::MultiIndex MultiIndex;
    typedef typename PartitionList::Mesh Mesh;
    
    SPPartitionPool ( const Mesh &localMesh, const Mesh &globalMesh,
                      const MultiIndex &overlap, unsigned int periodic = 0 );

    template< PartitionIteratorType pitype >
    const PartitionList &get () const;

    template< int codim >
    PartitionType
    partitionType ( const MultiIndex &id, const unsigned int number ) const;

    const Mesh &globalMesh () const;

  private:
    Partition
    closedPartition ( const Mesh &localMesh, const unsigned int number ) const;
    Partition
    openPartition ( const Mesh &localMesh, const unsigned int number ) const;

    Mesh globalMesh_;
    unsigned int periodic_;

    PartitionList interiorList_;
    PartitionList interiorBorderList_;
    PartitionList overlapList_;
    PartitionList overlapFrontList_;
    PartitionList allList_;
    PartitionList ghostList_;
  };



  // Implementation of SPPartitionPool
  // ---------------------------------

  template< int dim >
  inline SPPartitionPool< dim >
    ::SPPartitionPool ( const Mesh &localMesh, const Mesh &globalMesh,
                        const MultiIndex &overlap, unsigned int periodic )
  : globalMesh_( globalMesh ),
    periodic_( periodic )
  {
    // generate Interior and InteriorBorder
    interiorList_ += openPartition( localMesh, 0 );
    interiorBorderList_ += closedPartition( localMesh, 0 );
    interiorList_.updateCache();
    interiorBorderList_.updateCache();

    // detect which directions have to be split in the overlap partition
    const MultiIndex globalWidth = globalMesh.width();
    Mesh overlapMesh = localMesh.grow( overlap );
    const MultiIndex overlapWidth = overlapMesh.width();
    int n = 0;
    int shift[ dimension ];
    int dir[ dimension ];
    for( int i = 0; i < dimension; ++i )
    {
      if( (periodic & (1 << i)) == 0 )
        continue;

      if( overlapWidth[ i ] >= globalWidth[ i ] )
      {
        MultiIndex begin = overlapMesh.begin();
        MultiIndex end = overlapMesh.end();
        begin[ i ] = globalMesh.begin()[ i ];
        end[ i ] = globalMesh.end()[ i ];
        overlapMesh = Mesh( begin, end );
        continue;
      }

      shift[ n ] = 0;
      if( overlapMesh.begin()[ i ] < globalMesh.begin()[ i ] )
        shift[ n ] += globalWidth[ i ];
      if( overlapMesh.end()[ i ] > globalMesh.end()[ i ] )
        shift[ n ] -= globalWidth[ i ];
      if( shift[ n ] != 0 )
        dir[ n++ ] = i;
    }

    // generate Overlap and OverlapFront
    const unsigned int size = 1 << n;
    for( unsigned int d = 0; d < size; ++d )
    {
      MultiIndex s = MultiIndex::zero();
      for( int i = 0; i < n; ++i )
        s[ dir[ i ] ] = ((d >> i)&1)*shift[ i ];
      Partition open = openPartition( globalMesh.intersect( overlapMesh + s ), d );
      Partition closed = closedPartition( globalMesh.intersect( overlapMesh + s ), d );
      for( int i = 0; i < n; ++i )
      {
        const int j = (shift[ i ] < 0) ^ ((d >> i)&1);
        open.neighbor( 2*i + j ) = d ^ (1 << i);
        closed.neighbor( 2*i + j ) = d ^ (1 << i);
      }
      overlapList_ += open;
      overlapFrontList_ += closed;
    }
    overlapList_.updateCache();
    overlapFrontList_.updateCache();

    // generate All
    allList_ = overlapFrontList_;
  }


  template< int dim >
  template< PartitionIteratorType pitype >
  inline const typename SPPartitionPool< dim >::PartitionList &
  SPPartitionPool< dim >::get () const
  {
    switch( pitype )
    {
    case Interior_Partition:
      return interiorList_;

    case InteriorBorder_Partition:
      return interiorBorderList_;

    case Overlap_Partition:
      return overlapList_;

    case OverlapFront_Partition:
      return overlapFrontList_;

    case All_Partition:
      return allList_;

    case Ghost_Partition:
      return ghostList_;

    default:
      DUNE_THROW( GridError, "No such PartitionIteratorType." );
    }
  }


  template< int dim >
  template< int codim >
  inline PartitionType
  SPPartitionPool< dim >
    ::partitionType ( const MultiIndex &id, const unsigned int number ) const
  {
    assert( allList_.contains( id, number ) );
    if( interiorBorderList_.contains( id, number ) )
      return ((codim == 0) || interiorList_.contains( id, number )) ? InteriorEntity : BorderEntity;
    else if( overlapFrontList_.contains( id, number ) )
      return ((codim == 0) || overlapList_.contains( id, number )) ? OverlapEntity : FrontEntity;
    else
      return GhostEntity;
  }


  template< int dim >
  inline const typename SPPartitionPool< dim >::Mesh &
  SPPartitionPool< dim >::globalMesh () const
  {
    return globalMesh_;
  }


  template< int dim >
  inline typename SPPartitionPool< dim >::Partition
  SPPartitionPool< dim >
    ::closedPartition ( const Mesh &localMesh, const unsigned int number ) const
  {
    const MultiIndex &lbegin = localMesh.begin();
    const MultiIndex &lend = localMesh.end();
    const MultiIndex &gbegin = globalMesh().begin();
    const MultiIndex &gend = globalMesh().end();

    // create partition
    Partition partition( 2*lbegin, 2*lend, globalMesh(), number );
    
    // deal with self-neighborship (periodicity)
    for( int i = 0; i < dimension; ++i )
    {
      if( (periodic_ & (1 << i)) == 0 )
        continue;
      if( (lbegin[ i ] == gbegin[ i ]) && (lend[ i ] == gend[ i ]) )
        partition.neighbor( 2*i ) = partition.neighbor( 2*i+1 ) = number;
    }

    return partition;
  }


  template< int dim >
  inline typename SPPartitionPool< dim >::Partition
  SPPartitionPool< dim >
    ::openPartition ( const Mesh &localMesh, const unsigned int number ) const
  {
    const MultiIndex &lbegin = localMesh.begin();
    const MultiIndex &lend = localMesh.end();
    const MultiIndex &gbegin = globalMesh().begin();
    const MultiIndex &gend = globalMesh().end();

    // create partition
    MultiIndex begin, end;
    for( int i = 0; i < dimension; ++i )
    {
      begin[ i ] = 2*lbegin[ i ] + int( lbegin[ i ] != gbegin[ i ] );
      end[ i ] = 2*lend[ i ] - int( lend[ i ] != gend[ i ] );
    }
    Partition partition( begin, end, globalMesh(), number );
    
    // deal with self-neighborship (periodicity)
    for( int i = 0; i < dimension; ++i )
    {
      if( (periodic_ & (1 << i)) == 0 )
        continue;
      if( (lbegin[ i ] == gbegin[ i ]) && (lend[ i ] == gend[ i ]) )
        partition.neighbor( 2*i ) = partition.neighbor( 2*i+1 ) = number;
    }

    return partition;
  }

}

#endif // #ifndef DUNE_SPGRID_PARTITIONPOOL_HH

#ifndef DUNE_SPGRID_PARTITIONPOOL_HH
#define DUNE_SPGRID_PARTITIONPOOL_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/exceptions.hh>

#include <dune/grid/spgrid/cachedpartitionlist.hh>
#include <dune/grid/spgrid/topology.hh>

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
    typedef SPTopology< dimension > Topology;

    typedef typename PartitionList::Partition Partition;
    typedef typename PartitionList::MultiIndex MultiIndex;
    typedef typename PartitionList::Mesh Mesh;
    
    SPPartitionPool ( const Mesh &localMesh, const Mesh &globalMesh,
                      const MultiIndex &overlap, const Topology &topology );

    template< PartitionIteratorType pitype >
    const PartitionList &get () const;

    template< int codim >
    PartitionType
    partitionType ( const MultiIndex &id, const unsigned int number ) const;

    const Mesh &globalMesh () const { return globalMesh_; }
    const MultiIndex &overlap () const { return overlap_; }
    const Topology &topology () const { return topology_; }

  private:
    Partition makePartition ( const Mesh &localMesh, const unsigned int number,
                              const unsigned int open ) const;

    Mesh globalMesh_;
    MultiIndex overlap_;
    Topology topology_;

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
                        const MultiIndex &overlap, const Topology &topology )
  : globalMesh_( globalMesh ),
    overlap_( overlap ),
    topology_( topology )
  {
    // generate Interior and InteriorBorder
    interiorList_ += makePartition( localMesh, 0, (1 << dimension) - 1 );
    interiorBorderList_ += makePartition( localMesh, 0, 0 );
    interiorList_.updateCache();
    interiorBorderList_.updateCache();

    // detect which directions have to be split in the overlap partition
    const MultiIndex globalWidth = globalMesh.width();
    Mesh overlapMesh = localMesh.grow( overlap );
    const MultiIndex overlapWidth = overlapMesh.width();
    int n = 0;
    int shift[ dimension ];
    int dir[ dimension ];
    unsigned int openOverlap = 0;

    // initialize shift
    for( int i = 0; i < dimension; ++i )
      shift[ i ] = 0;

    for( int i = 0; i < dimension; ++i )
    {
      openOverlap |= (overlap[ i ] > 0 ? (1 << i) : 0);
      if( !topology.hasNeighbor( 0, 2*i ) )
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
      Partition open = makePartition( globalMesh.intersect( overlapMesh + s ), d, openOverlap );
      Partition closed = makePartition( globalMesh.intersect( overlapMesh + s ), d, 0 );
      for( int i = 0; i < n; ++i )
      {
        const int j = (shift[ i ] < 0) ^ ((d >> i)&1);
        open.neighbor( 2*dir[ i ] + j ) = d ^ (1 << i);
        closed.neighbor( 2*dir[ i ] + j ) = d ^ (1 << i);
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
  inline typename SPPartitionPool< dim >::Partition
  SPPartitionPool< dim >
    ::makePartition ( const Mesh &localMesh, const unsigned int number,
                      const unsigned int open ) const
  {
    const MultiIndex &lbegin = localMesh.begin();
    const MultiIndex &lend = localMesh.end();
    const MultiIndex &gbegin = globalMesh().begin();
    const MultiIndex &gend = globalMesh().end();

    // create partition
    MultiIndex begin, end;
    for( int i = 0; i < dimension; ++i )
    {
      const int o = ((open >> i) & 1);
      begin[ i ] = 2*lbegin[ i ] + o*int( lbegin[ i ] != gbegin[ i ] );
      end[ i ] = 2*lend[ i ] - o*int( lend[ i ] != gend[ i ] );
    }
    Partition partition( begin, end, globalMesh(), number );
    
    // deal with self-neighborship (periodicity)
    for( int i = 0; i < dimension; ++i )
    {
      if( !topology().hasNeighbor( 0, 2*i ) )
        continue;
      if( (lbegin[ i ] == gbegin[ i ]) && (lend[ i ] == gend[ i ]) )
        partition.neighbor( 2*i ) = partition.neighbor( 2*i+1 ) = number;
    }

    return partition;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_PARTITIONPOOL_HH

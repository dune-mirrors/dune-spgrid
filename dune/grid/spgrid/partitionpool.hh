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
    typedef SPCachedPartitionList< dim > PartitionList;

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

  private:
    Partition
    closedPartition ( const Mesh &localMesh, const unsigned int number ) const;
    Partition
    openPartition ( const Mesh &localMesh, const unsigned int number ) const;

    Mesh globalMesh_;

    PartitionList interior_;
    PartitionList interiorBorder_;
    PartitionList overlap_;
    PartitionList overlapFront_;
    PartitionList all_;
    PartitionList ghost_;
  };



  // Implementation of SPPartitionPool
  // ---------------------------------

  template< int dim >
  inline SPPartitionPool< dim >
    ::SPPartitionPool ( const Mesh &localMesh, const Mesh &globalMesh,
                        const MultiIndex &overlap, unsigned int periodic )
  : globalMesh_( globalMesh )
  {
    interiorBorder_ += closedPartition( localMesh, 0 );
    interior_ += openPartition( localMesh, 0 );
    interiorBorder_.updateCache();
    interior_.updateCache();

    MultiIndex globalWidth = globalMesh.width();
    std::vector< Mesh > overlapMesh( 1, localMesh.grow( overlap ) );
    for( int i = 0; i < dim; ++i )
    {
      if( (periodic & (1 << i)) == 0 )
        continue;

      MultiIndex shift = MultiIndex::zero();
      if( overlapMesh[ 0 ].begin()[ i ] < globalMesh.begin()[ i ] )
        shift[ i ] += globalWidth[ i ];
      if( overlapMesh[ 0 ].end()[ i ] > globalMesh.end()[ i ] )
        shift[ i ] -= globalWidth[ i ];
      if( shift[ i ] == 0 )
        continue;

      const size_t size = overlapMesh.size();
      overlapMesh.reserve( 2*size );
      for( size_t i = 0; i < size; ++i )
        overlapMesh.push_back( overlapMesh[ i ] + shift );
    }

    const size_t size = overlapMesh.size();
    for( size_t i = 0; i < size; ++i )
    {
      overlapFront_ += closedPartition( globalMesh.intersect( overlapMesh[ i ] ), i );
      overlap_ += openPartition( globalMesh.intersect( overlapMesh[ i ] ), i );
    }
    overlapFront_.updateCache();
    overlap_.updateCache();

    all_ = overlapFront_;
  }


  template< int dim >
  template< PartitionIteratorType pitype >
  inline const typename SPPartitionPool< dim >::PartitionList &
  SPPartitionPool< dim >::get () const
  {
    switch( pitype )
    {
    case Interior_Partition:
      return interior_;

    case InteriorBorder_Partition:
      return interiorBorder_;

    case Overlap_Partition:
      return overlap_;

    case OverlapFront_Partition:
      return overlapFront_;

    case All_Partition:
      return all_;

    case Ghost_Partition:
      return ghost_;

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
    assert( all_.contains( id, number ) );
    if( interiorBorder_.contains( id, number ) )
      return ((codim == 0) || interior_.contains( id, number )) ? InteriorEntity : BorderEntity;
    else if( overlapFront_.contains( id, number ) )
      return ((codim == 0) || overlap_.contains( id, number )) ? OverlapEntity : FrontEntity;
    else
      return GhostEntity;
  }


  template< int dim >
  inline typename SPPartitionPool< dim >::Partition
  SPPartitionPool< dim >
    ::closedPartition ( const Mesh &localMesh, const unsigned int number ) const
  {
    const MultiIndex &lbegin = localMesh.begin();
    const MultiIndex &lend = localMesh.end();

    MultiIndex begin = 2*lbegin;
    MultiIndex end = 2*lend;
    return Partition( begin, end, number );
  }


  template< int dim >
  inline typename SPPartitionPool< dim >::Partition
  SPPartitionPool< dim >
    ::openPartition ( const Mesh &localMesh, const unsigned int number ) const
  {
    const MultiIndex &lbegin = localMesh.begin();
    const MultiIndex &lend = localMesh.end();
    const MultiIndex &gbegin = globalMesh_.begin();
    const MultiIndex &gend = globalMesh_.end();

    MultiIndex begin, end;
    for( int i = 0; i < dim; ++i )
    {
      begin[ i ] = 2*lbegin[ i ] + int( lbegin[ i ] != gbegin[ i ] );
      end[ i ] = 2*lend[ i ] - int( lend[ i ] != gend[ i ] );
    }
    return Partition( begin, end, number );
  }

}

#endif // #ifndef DUNE_SPGRID_PARTITIONPOOL_HH

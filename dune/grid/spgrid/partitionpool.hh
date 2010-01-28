#ifndef DUNE_SPGRID_PARTITIONPOOL_HH
#define DUNE_SPGRID_PARTITIONPOOL_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/exceptions.hh>

#include <dune/grid/spgrid/partitionlist.hh>

namespace Dune
{

  // SPPartitionPool
  // ---------------

  template< int dim >
  class SPPartitionPool
  {
    typedef SPPartitionPool< dim > This;

  public:
    typedef SPPartitionList< dim > PartitionList;

    typedef typename PartitionList::Partition Partition;
    typedef typename PartitionList::MultiIndex MultiIndex;
    typedef typename PartitionList::Mesh Mesh;
    
    SPPartitionPool ( const Mesh &localMesh, const Mesh &globalMesh,
                      const MultiIndex &overlap, unsigned int periodic = 0 );

    template< PartitionIteratorType pitype >
    const PartitionList &get () const;

    template< int codim >
    PartitionType partitionType ( const MultiIndex &id ) const;

  private:
    static Partition
    removeBorder ( const Mesh &localMesh, const Mesh &globalMesh );

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
  {
    interior_ += removeBorder( localMesh, globalMesh );
    interiorBorder_ += Partition( localMesh );

    Mesh overlapMesh = localMesh.grow( overlap );
    overlap_ += removeBorder( globalMesh.intersect( overlapMesh ), globalMesh );
    overlapFront_ += Partition( globalMesh.intersect( overlapMesh ) );
    
    Mesh allMesh = localMesh.grow( overlap ).grow( 1 );
    all_ += Partition( globalMesh.intersect( allMesh ) );
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
  SPPartitionPool< dim >::partitionType ( const MultiIndex &id ) const
  {
    assert( all_.contains( id ) );
    // both, interior_ and interiorBorder_ contain exactly one partition, so use this information
    if( interiorBorder_.begin()->contains( id ) )
      return ((codim == 0) || interior_.begin()->contains( id )) ? InteriorEntity : BorderEntity;
    else if( overlapFront_.contains( id ) )
      return ((codim == 0) || overlap_.contains( id )) ? OverlapEntity : FrontEntity;
    else
      return GhostEntity;
  }


  template< int dim >
  inline typename SPPartitionPool< dim >::Partition
  SPPartitionPool< dim >::removeBorder ( const Mesh &localMesh, const Mesh &globalMesh )
  {
    const MultiIndex &lbegin = localMesh.begin();
    const MultiIndex &lend = localMesh.end();
    const MultiIndex &gbegin = globalMesh.begin();
    const MultiIndex &gend = globalMesh.end();

    MultiIndex begin, end;
    for( int i = 0; i < dim; ++i )
    {
      begin[ i ] = 2*lbegin[ i ] + int( lbegin[ i ] != gbegin[ i ] );
      end[ i ] = 2*lend[ i ] - int( lend[ i ] != gend[ i ] );
    }
    return Partition( begin, end );
  }

}

#endif // #ifndef DUNE_SPGRID_PARTITIONPOOL_HH

#ifndef DUNE_SPGRID_LINKAGE_HH
#define DUNE_SPGRID_LINKAGE_HH

#include <dune/grid/spgrid/partitionlist.hh>
#include <dune/grid/spgrid/partitionpool.hh>

namespace Dune
{

  // SPLinkage
  // ---------

  template< int dim >
  class SPLinkage
  {
    typedef SPLinkage< dim > This;

  public:
    typedef SPPartitionPool< dim > PartitionPool;
    typedef SPPartitionList< dim > PartitionList;

    typedef typename PartitionPool::Mesh Mesh;
    typedef typename PartitionPool::MultiIndex MultiIndex;

    struct Interface;

    SPLinkage ( const int rank,
                const PartitionPool &localPool,
                const std::vector< Mesh > &decomposition );

    const Interface &interface ( const InterfaceType interface ) const;

  private:
    template< InterfaceType interface >
    bool build ( const int rank, const PartitionPool &localPool, const PartitionPool &remotePool );

   const PartitionList *intersect ( const PartitionList &local, const PartitionList &remote ) const;

    // note: We use the knowledge that interfaces are numbered 0, ..., 4.
    Interface interface_[ 5 ];
  };



  // SPLinkage::Interface
  // --------------------

  template< int dim >
  struct SPLinkage< dim >::Interface
  {
  };



  // SPCommunicationInterface
  // ------------------------

  template< InterfaceType >
  struct SPCommunicationInterface;

  template<>
  struct SPCommunicationInterface< InteriorBorder_InteriorBorder_Interface >
  {
    static const PartitionIteratorType sendPartition = InteriorBorder_Partition;
    static const PartitionIteratorType receivePartition = InteriorBorder_Partition;
  };

  template<>
  struct SPCommunicationInterface< InteriorBorder_All_Interface >
  {
    static const PartitionIteratorType sendPartition = InteriorBorder_Partition;
    static const PartitionIteratorType receivePartition = All_Partition;
  };

  template<>
  struct SPCommunicationInterface< Overlap_OverlapFront_Interface >
  {
    static const PartitionIteratorType sendPartition = Overlap_Partition;
    static const PartitionIteratorType receivePartition = OverlapFront_Partition;
  };

  template<>
  struct SPCommunicationInterface< Overlap_All_Interface >
  {
    static const PartitionIteratorType sendPartition = Overlap_Partition;
    static const PartitionIteratorType receivePartition = All_Partition;
  };

  template<>
  struct SPCommunicationInterface< All_All_Interface >
  {
    static const PartitionIteratorType sendPartition = All_Partition;
    static const PartitionIteratorType receivePartition = All_Partition;
  };



  // Implementation of SPLinkage
  // ---------------------------

  template< int dim >
  inline SPLinkage< dim >
    ::SPLinkage ( const int rank,
                  const PartitionPool &localPool,
                  const std::vector< Mesh > &decomposition )
  {
    const int size = decomposition.size();

    const Mesh &globalMesh = localPool.globalMesh();
    const MultiIndex &overlap = localPool.overlap();
    const unsigned int periodic = localPool.periodic();

    for( int p = 0; p < size; ++p )
    {
      if( p == rank )
        continue;
      PartitionPool remotePool( decomposition[ p ], globalMesh, overlap, periodic );
      if( build< All_All_Interface >( p, localPool, remotePool ) )
      {
        build< InteriorBorder_InteriorBorder_Interface >( p, localPool, remotePool );
        build< InteriorBorder_All_Interface >( p, localPool, remotePool );
        build< Overlap_OverlapFront_Interface >( p, localPool, remotePool );
        build< Overlap_All_Interface >( p, localPool, remotePool );
      }
    }
  }


  template< int dim >
  inline const typename SPLinkage< dim >::Interface &
  SPLinkage< dim >::interface ( const InterfaceType interface ) const
  {
    assert( (int( interface ) >= 0) && (int( interface ) < 5) );
    return interface_[ int( interface ) ];
  }


  template< int dim >
  template< InterfaceType interface >
  inline bool SPLinkage< dim >
    ::build ( const int rank, const PartitionPool &localPool, const PartitionPool &remotePool )
  {
    const PartitionType piSend = SPCommunicationInterface< interface >::sendPartition;
    const PartitionType piReceive = SPCommunicationInterface< interface >::receivePartition;

    // build intersection lists
    const PartitionList *sendList, *receiveList;
    sendList = intersect ( localPool.template get< piSend >(), remotePool.template get< piReceive >() );
    if( piSend != piReceive )
      receiveList = intersect( localPool.template get< piReceive >(), remotePool.template get< piSend >() );
    else
      receiveList = sendList;

    // if both lists are empty, no communication is necessary
    if( sendList->empty() && receiveList->empty() )
    {
      if( sendList != receiveList )
        delete receiveList;
      delete sendList;
      return false;
    }



    return true;
  }


  template< int dim >
  inline const typename SPLinkage< dim >::PartitionList *
  SPLinkage< dim >::intersect ( const PartitionList &local, const PartitionList &remote ) const
  {
    typedef SPBasicPartition< dim > Intersection;
    typedef typename PartitionList::Iterator Iterator;
    typedef typename PartitionList::Partition Partition;

    PartitionList *link = new PartitionList;
    for( Iterator lit = local.begin(); lit; ++lit )
    {
      const int number = lit->number();
      for( Iterator rit = remote.begin(); rit; ++rit )
      {
        Intersection intersection = lit->intersect( *rit );
        if( !intersection.empty() )
          link += Partition( intersection, number );
      }
    }
    return link;
  }

}

#endif // #ifndef DUNE_SPGRID_LINKAGE_HH

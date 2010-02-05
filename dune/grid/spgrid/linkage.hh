#ifndef DUNE_SPGRID_LINKAGE_HH
#define DUNE_SPGRID_LINKAGE_HH

#include <dune/grid/spgrid/partitionpool.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< InterfaceType >
  struct SPCommunicationInterface;



  // SPLinkage
  // ---------

  template< int dim >
  class SPLinkage
  {
    typedef SPLinkage< dim > This;

  public:
    typedef SPPartitionPool< dim > PartitionPool;

    typedef typename PartitionPool::Mesh Mesh;
    typedef typename PartitionPool::MultiIndex MultiIndex;

    struct Interface;

    SPLinkage ( const int rank,
                const PartitionPool &localPool,
                const std::vector< Mesh > &decomposition );

    const Interface &interface ( const InterfaceType interface ) const;

  private:
    template< InterfaceType interface >
    bool build ( const int rank, const PartitionPool &localPool, const PartitionPool &removePool );

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
    ::build ( const int rank, const PartitionPool &localPool, const PartitionPool &removePool )
  {
    return false;
  }

}

#endif // #ifndef DUNE_SPGRID_LINKAGE_HH

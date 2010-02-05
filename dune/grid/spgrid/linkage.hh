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

    class Interface;

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
  class SPLinkage< dim >::Interface
  {
    struct Node;

    typedef std::vector< Node > NodeContainer;

  public:
    typedef typename NodeContainer::const_iterator Iterator;

    Interface ();
    Interface ( const This &other );

    ~Interface ();

    Iterator begin () const;
    Iterator end () const;

    void add ( const int rank, const PartitionList *sendList, const PartitionList *receiveList );

  private:
    NodeContainer nodes_;
  };



  // SPLinkage::Interface::Node
  // --------------------------

  template< int dim >
  struct SPLinkage< dim >::Interface::Node
  {
    dune_static_assert( (int( ForwardCommunication ) == 0) && (int( BackwardCommunication ) == 1),
                        "enumeration CommunicationDirection has changed." );

    Node ( const int rank, const PartitionList *sendList, const PartitionList *receiveList );

    int rank () const;
    const PartitionList &sendList ( const CommunicationDirection direction = ForwardCommunication ) const;
    const PartitionList &receiveList ( const CommunicationDirection direction = ForwardCommunication ) const;

    void destroy ();

  private:
    const int rank_;
    const PartitionList *partitionList_[ 2 ];
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

    assert( (interface >= 0) && (interface < 5) );
    interface_[ interface ].add( rank, sendList, receiveList );
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



  // Implementation of SPLinkage::Interface
  // --------------------------------------

  template< int dim >
  inline SPLinkage< dim >::Interface::Interface ()
  {}


  template< int dim >
  inline SPLinkage< dim >::Interface::Interface ( const This &other )
  {
    nodes_.reserve( other.nodes_.size() );
    const Iterator end = other.end();
    for( Iterator it = other.begin(); it != end; ++it )
    {
      const PartitionList *sendList = new PartitionList( it->sendList() );
      const PartitionList *receiveList = new PartitionList( it->receiveList() );
      nodes_->push_back( it->rank(), sendList, receiveList );
    }
  }


  template< int dim >
  inline SPLinkage< dim >::Interface::~Interface ()
  {
    typedef typename NodeContainer::iterator Iterator;
    const Iterator end = nodes_.end();
    for( Iterator it = nodes_.begin(); it != end; ++it )
      it->destroy();
  }


  template< int dim >
  inline typename SPLinkage< dim >::Interface::Iterator
  SPLinkage< dim >::Interface::begin () const
  {
    return nodes_.begin();
  }


  template< int dim >
  inline typename SPLinkage< dim >::Interface::Iterator
  SPLinkage< dim >::Interface::end () const
  {
    return nodes_.end();
  }


  template< int dim >
  inline void SPLinkage< dim >::Interface
    ::add ( const int rank, const PartitionList *sendList, const PartitionList *receiveList )
  {
    nodes_.push_back( Node( rank, sendList, receiveList ) );
  }



  // Implementation of SPLinkage::Interface
  // --------------------------------------

  template< int dim >
  inline SPLinkage< dim >::Interface::Node
    ::Node ( const int rank, const PartitionList *sendList, const PartitionList *receiveList )
  : rank_( rank )
  {
    partitionList_[ 0 ] = sendList;
    partitionList_[ 1 ] = receiveList;
  }


  template< int dim >
  inline int SPLinkage< dim >::Interface::Node::rank () const
  {
    assert( rank_ >= 0 );
    return rank_;
  }


  template< int dim >
  inline const typename SPLinkage< dim >::PartitionList &
  SPLinkage< dim >::Interface::Node::sendList ( const CommunicationDirection direction ) const
  {
    assert( (int( direction ) >= 0) && (int( direction ) <= 1) );
    assert( partitionList_[ direction ] != 0 );
    return *partitionList_[ direction ];
  }

  template< int dim >
  inline const typename SPLinkage< dim >::PartitionList &
  SPLinkage< dim >::Interface::Node::receiveList ( const CommunicationDirection direction ) const
  {
    assert( (int( direction ) >= 0) && (int( direction ) <= 1) );
    assert( partitionList_[ 1 - direction ] != 0 );
    return *partitionList_[ 1 - direction ];
  }


  template< int dim >
  inline void SPLinkage< dim >::Interface::Node::destroy ()
  {
    rank_ = -1;
    if( partitionList_[ 0 ] != partitionList_[ 1 ] )
      delete partitionList_[ 1 ];
    delete partitionList_[ 0 ];
    partitionList_[ 0 ] = partitionList_[ 1 ] = 0;
  }

}

#endif // #ifndef DUNE_SPGRID_LINKAGE_HH

#ifndef DUNE_SPGRID_COMMUNICATION_HH
#define DUNE_SPGRID_COMMUNICATION_HH

#include <dune/common/forloop.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/common/parallel/mpitraits.hh>

#include <dune/grid/common/datahandleif.hh>

#include <dune/grid/spgrid/iterator.hh>
#include <dune/grid/spgrid/messagebuffer.hh>

namespace Dune
{

  // SPCommunicationTraits
  // ---------------------

  template< class Comm >
  struct SPCommunicationTraits
  {
    typedef Dune::CollectiveCommunication< Comm > CollectiveCommunication;

    template< class C >
    static CollectiveCommunication comm ( const C & )
    {
      return defaultComm();
    }

    static CollectiveCommunication defaultComm ()
    {
      return CollectiveCommunication();
    }
  };

#if HAVE_MPI
  template<>
  struct SPCommunicationTraits< MPI_Comm >
  {
    typedef Dune::CollectiveCommunication< MPI_Comm > CollectiveCommunication;

    static CollectiveCommunication comm ( const MPI_Comm &mpiComm )
    {
      return CollectiveCommunication( mpiComm );
    }

    static CollectiveCommunication defaultComm ()
    {
      return comm( MPI_COMM_WORLD );
    }
  };
#endif // #if HAVE_MPI



  // SPCommunication
  // ---------------

  template< class Grid, class DataHandle >
  class SPCommunication
  {
    template< int codim >
    struct Codim;

  public:
    static const unsigned int dimension = Grid::dimension;

    typedef SPGridLevel< Grid > GridLevel;
    typedef SPPartitionList< dimension > PartitionList;

    typedef typename DataHandle::DataType DataType;

    typedef typename GridLevel::CommInterface Interface;

  private:
    typedef SPPackedMessageWriteBuffer< typename Grid::CollectiveCommunication > WriteBuffer;
    typedef SPPackedMessageReadBuffer< typename Grid::CollectiveCommunication > ReadBuffer;

  public:
    SPCommunication ( const GridLevel &gridLevel, DataHandle &dataHandle,
                      InterfaceType iftype, CommunicationDirection dir );

    SPCommunication ( const SPCommunication & ) = delete;
    SPCommunication ( SPCommunication &&other );

    ~SPCommunication () { wait(); }

    bool pending () const { return bool( interface_ ); }

    void wait ();

  private:
    static int getTag ()
    {
      static unsigned char counter = 0;
      return int( counter++ ) + 1536;
    };

    const GridLevel &gridLevel_;
    DataHandle &dataHandle_;
    const Interface *interface_;
    CommunicationDirection dir_;
    int tag_;
    std::vector< WriteBuffer > writeBuffers_;
  };



  // SPCommunication::Codim
  // ----------------------

  template< class Grid, class DataHandle >
  template< int codim >
  struct SPCommunication< Grid, DataHandle >::Codim
  {
    typedef SPPartitionIterator< codim, const Grid > Iterator;

    static void
    apply ( const GridLevel &gridLevel, DataHandle &dataHandle,
            const PartitionList &partitionList, WriteBuffer &buffer );
    static void
    apply ( const GridLevel &gridLevel, DataHandle &dataHandle,
            const PartitionList &partitionList, ReadBuffer &buffer );
  };



  // Implementation of SPCommunication
  // ---------------------------------

  template< class Grid, class DataHandle >
  inline SPCommunication< Grid, DataHandle >
    ::SPCommunication ( const GridLevel &gridLevel, DataHandle &dataHandle,
                        InterfaceType iftype, CommunicationDirection dir )
    : gridLevel_( gridLevel ),
      dataHandle_( dataHandle ),
      interface_( &gridLevel.commInterface( iftype ) ),
      dir_( dir ),
      tag_( getTag() )
  {
    const typename Interface::Iterator end = interface_->end();
    for( typename Interface::Iterator it = interface_->begin(); it != end; ++it )
    {
      writeBuffers_.emplace_back( gridLevel.grid().comm() );
      ForLoop< Codim, 0, dimension >::apply( gridLevel_, dataHandle_, it->sendList( dir ), writeBuffers_.back() );
      writeBuffers_.back().send( it->rank(), tag_ );
    }
  }


  template< class Grid, class DataHandle >
  inline SPCommunication< Grid, DataHandle >::SPCommunication ( SPCommunication &&other )
    : gridLevel_( other.gridLevel_ ),
      dataHandle_( other.dataHandle_ ),
      interface_( other.interface_ ),
      dir_( other.dir_ ),
      tag_( other.tag_ ),
      writeBuffers_( std::move( other.writeBuffers_ ) )
  {
    other.interface_ = nullptr;
  }


  template< class Grid, class DataHandle >
  inline void SPCommunication< Grid, DataHandle >::wait ()
  {
    if( !pending() )
      return;

    ReadBuffer buffer( gridLevel_.grid().comm() );
    const typename Interface::Iterator end = interface_->end();
    for( typename Interface::Iterator it = interface_->begin(); it != end; ++it )
    {
      buffer.receive( it->rank(), tag_ );
      ForLoop< Codim, 0, dimension >::apply( gridLevel_, dataHandle_, it->receiveList( dir_ ), buffer );
    }

    for( typename std::vector< WriteBuffer >::iterator it = writeBuffers_.begin(); it != writeBuffers_.end(); ++it )
      it->wait();
    writeBuffers_.clear();
  }



  // Implementation of SPCommunication::Codim
  // ----------------------------------------

  template< class Grid, class DataHandle >
  template< int codim >
  inline void SPCommunication< Grid, DataHandle >::Codim< codim >
    ::apply ( const GridLevel &gridLevel, DataHandle &dataHandle,
              const PartitionList &partitionList, WriteBuffer &buffer )
  {
    if( dataHandle.contains( dimension, codim ) )
    {
      const Iterator end( gridLevel, partitionList, typename Iterator::End() );
      for( Iterator it( gridLevel, partitionList, typename Iterator::Begin() ); it != end; ++it )
      {
        buffer.write( static_cast< int >( dataHandle.size( *it ) ) );
        dataHandle.gather( buffer, *it ); 
      }
    }
  }

  template< class Grid, class DataHandle >
  template< int codim >
  inline void SPCommunication< Grid, DataHandle >::Codim< codim >
    ::apply ( const GridLevel &gridLevel, DataHandle &dataHandle,
              const PartitionList &partitionList, ReadBuffer &buffer )
  {
    if( dataHandle.contains( dimension, codim ) )
    {
      const Iterator end( gridLevel, partitionList, typename Iterator::End() );
      for( Iterator it( gridLevel, partitionList, typename Iterator::Begin() ); it != end; ++it )
      {
        int n;
        buffer.read( n );
        dataHandle.scatter( buffer, *it, n ); 
      }
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_COMMUNICATION_HH

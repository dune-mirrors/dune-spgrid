#ifndef DUNE_SPGRID_COMMUNICATION_HH
#define DUNE_SPGRID_COMMUNICATION_HH

#include <dune/common/forloop.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/common/parallel/mpitraits.hh>

#include <dune/grid/common/datahandleif.hh>

#include <dune/grid/spgrid/iterator.hh>

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
    template< class T >
    struct WriteBuffer;

    template< class T >
    struct ReadBuffer;

    template< int codim >
    struct Codim;

  public:
    static const unsigned int dimension = Grid::dimension;

    typedef SPGridLevel< Grid > GridLevel;
    typedef SPPartitionList< dimension > PartitionList;

    typedef typename DataHandle::DataType DataType;

    typedef typename GridLevel::CommInterface Interface;

  private:
    typedef std::pair< WriteBuffer< int >, WriteBuffer< DataType > > Packet;

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
      return 2*int( counter++ ) + 1536;
    };

    const GridLevel &gridLevel_;
    DataHandle &dataHandle_;
    const Interface *interface_;
    CommunicationDirection dir_;
    int tag_;
    std::vector< Packet * > packets_;
  };



  // SPCommunication::WriteBuffer
  // ----------------------------

  template< class Grid, class DataHandle >
  template< class T >
  struct SPCommunication< Grid, DataHandle >::WriteBuffer
  {
    void write ( const T &value );

    template< class C >
    void send ( int rank, int tag, const CollectiveCommunication< C > &comm );
    template< class C >
    void wait ( const CollectiveCommunication< C > &comm );

#if HAVE_MPI
    void send ( int rank, int tag, const CollectiveCommunication< MPI_Comm > &comm );
    void wait ( const CollectiveCommunication< MPI_Comm > &comm );
#endif // #if HAVE_MPI

  private:
    std::vector< T > buffer_;
#if HAVE_MPI
    MPI_Request request_;
#endif // #if HAVE_MPI
  };



  // SPCommunication::ReadBuffer
  // ---------------------------

  template< class Grid, class DataHandle >
  template< class T >
  struct SPCommunication< Grid, DataHandle >::ReadBuffer
  {
    ReadBuffer () : read_( buffer_.begin() ) {}

    template< class C >
    void receive ( int rank, int tag, const CollectiveCommunication< C > &comm )
    {}

#if HAVE_MPI
    void receive ( int rank, int tag, const CollectiveCommunication< MPI_Comm > &comm );
#endif // #if HAVE_MPI

    void read ( T &value );

  private:
    std::vector< T > buffer_;
    typename std::vector< T >::const_iterator read_;
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
            const PartitionList &partitionList,
            WriteBuffer< int > &sizes, WriteBuffer< DataType > &buffer );
    static void
    apply ( const GridLevel &gridLevel, DataHandle &dataHandle,
            const PartitionList &partitionList,
            ReadBuffer< int > &sizes, ReadBuffer< DataType > &data );
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
      Packet *packet = new Packet;
      ForLoop< Codim, 0, dimension >::apply( gridLevel_, dataHandle_, it->sendList( dir ), packet->first, packet->second );
      packet->first.send( it->rank(), tag_, gridLevel_.grid().comm() );
      packet->second.send( it->rank(), tag_+1, gridLevel_.grid().comm() );
      packets_.push_back( packet );
    }
  }


  template< class Grid, class DataHandle >
  inline SPCommunication< Grid, DataHandle >::SPCommunication ( SPCommunication &&other )
    : gridLevel_( other.gridLevel_ ),
      dataHandle_( other.dataHandle_ ),
      interface_( other.interface_ ),
      dir_( other.dir_ ),
      tag_( other.tag_ ),
      packets_( std::move( other.packets_ ) )
  {
    other.interface_ = nullptr;
  }


  template< class Grid, class DataHandle >
  inline void SPCommunication< Grid, DataHandle >::wait ()
  {
    if( !pending() )
      return;

    ReadBuffer< int > sizes;
    ReadBuffer< DataType > data;
    const typename Interface::Iterator end = interface_->end();
    for( typename Interface::Iterator it = interface_->begin(); it != end; ++it )
    {
      sizes.receive( it->rank(), tag_, gridLevel_.grid().comm() );
      data.receive( it->rank(), tag_+1, gridLevel_.grid().comm() );
      ForLoop< Codim, 0, dimension >::apply( gridLevel_, dataHandle_, it->receiveList( dir_ ), sizes, data );
    }

    for( typename std::vector< Packet * >::iterator it = packets_.begin(); it != packets_.end(); ++it )
    {
      (*it)->first.wait( gridLevel_.grid().comm() );
      (*it)->second.wait( gridLevel_.grid().comm() );
      delete *it;
    }
    packets_.clear();
  }



  // Implementation of SPCommunication::WriteBuffer
  // ----------------------------------------------

  template< class Grid, class DataHandle >
  template< class T >
  inline void
  SPCommunication< Grid, DataHandle >::WriteBuffer< T >::write ( const T &value )
  {
    buffer_.push_back( value );
  }


  template< class Grid, class DataHandle >
  template< class T >
  template< class C >
  inline void SPCommunication< Grid, DataHandle >::WriteBuffer< T >
    ::send ( int rank, int tag, const CollectiveCommunication< C > &comm )
  {}


  template< class Grid, class DataHandle >
  template< class T >
  template< class C >
  inline void SPCommunication< Grid, DataHandle >::WriteBuffer< T >
    ::wait ( const CollectiveCommunication< C > &comm )
  {}


#if HAVE_MPI
  template< class Grid, class DataHandle >
  template< class T >
  inline void SPCommunication< Grid, DataHandle >::WriteBuffer< T >
   ::send ( int rank, int tag, const CollectiveCommunication< MPI_Comm > &comm )
  {
    MPI_Datatype mpitype = MPITraits< T >::getType();
    MPI_Isend( &(buffer_[ 0 ]), buffer_.size(), mpitype, rank, tag, comm, &request_ );
  }


  template< class Grid, class DataHandle >
  template< class T >
  inline void SPCommunication< Grid, DataHandle >::WriteBuffer< T >
    ::wait ( const CollectiveCommunication< MPI_Comm > &comm )
  {
    MPI_Status status;
    MPI_Wait( &request_, &status );
  }
#endif // #if HAVE_MPI



  // Implementation of SPCommunication::ReadBuffer
  // ---------------------------------------------

#if HAVE_MPI
  template< class Grid, class DataHandle >
  template< class T >
  inline void SPCommunication< Grid, DataHandle >::ReadBuffer< T >
    ::receive ( int rank, int tag, const CollectiveCommunication< MPI_Comm > &comm )
  {
    MPI_Status status;
    MPI_Probe( rank, tag, comm, &status );

    MPI_Datatype mpitype = MPITraits< T >::getType();

    int count;
    MPI_Get_count( &status, mpitype, &count );
    buffer_.resize( count );

    MPI_Recv( &(buffer_[ 0 ]), buffer_.size(), mpitype, rank, tag, comm, &status );

    read_ = buffer_.begin();
  }
#endif // #if HAVE_MPI


  template< class Grid, class DataHandle >
  template< class T >
  inline void
  SPCommunication< Grid, DataHandle >::ReadBuffer< T >::read ( T &value )
  {
    if( read_ != buffer_.end() )
    {
      value = *read_;
      ++read_;
    }
    else
      DUNE_THROW( IOError, "Cannot read beyond the buffer's end." );
  }



  // Implementation of SPCommunication::Codim
  // ----------------------------------------

  template< class Grid, class DataHandle >
  template< int codim >
  inline void SPCommunication< Grid, DataHandle >::Codim< codim >
    ::apply ( const GridLevel &gridLevel, DataHandle &dataHandle,
              const PartitionList &partitionList,
              WriteBuffer< int > &sizes, WriteBuffer< DataType > &buffer )
  {
    if( dataHandle.contains( dimension, codim ) )
    {
      const Iterator end( gridLevel, partitionList, typename Iterator::End() );
      for( Iterator it( gridLevel, partitionList, typename Iterator::Begin() ); it != end; ++it )
      {
        sizes.write( dataHandle.size( *it ) );
        dataHandle.gather( buffer, *it ); 
      }
    }
  }

  template< class Grid, class DataHandle >
  template< int codim >
  inline void SPCommunication< Grid, DataHandle >::Codim< codim >
    ::apply ( const GridLevel &gridLevel, DataHandle &dataHandle,
              const PartitionList &partitionList,
              ReadBuffer< int > &sizes, ReadBuffer< DataType > &data )
  {
    if( dataHandle.contains( dimension, codim ) )
    {
      const Iterator end( gridLevel, partitionList, typename Iterator::End() );
      for( Iterator it( gridLevel, partitionList, typename Iterator::Begin() ); it != end; ++it )
      {
        int n;
        sizes.read( n );
        dataHandle.scatter( data, *it, n ); 
      }
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_COMMUNICATION_HH

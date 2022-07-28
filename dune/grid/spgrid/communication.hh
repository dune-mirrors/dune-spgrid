#ifndef DUNE_SPGRID_COMMUNICATION_HH
#define DUNE_SPGRID_COMMUNICATION_HH

#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/common/parallel/mpitraits.hh>
#include <dune/common/visibility.hh>

#include <dune/grid/common/exceptions.hh>
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
    typedef Dune::Communication< Comm > Communication;

    template< class C >
    static Communication comm ( const C & )
    {
      return defaultComm();
    }

    static Communication defaultComm ()
    {
      return Communication();
    }
  };

#if HAVE_MPI
  template<>
  struct SPCommunicationTraits< MPI_Comm >
  {
    typedef Dune::Communication< MPI_Comm > Communication;

    static Communication comm ( const MPI_Comm &mpiComm )
    {
      return Communication( mpiComm );
    }

    static Communication defaultComm ()
    {
      return comm( MPI_COMM_WORLD );
    }
  };
#endif // #if HAVE_MPI



  namespace __SPGrid
  {

    DUNE_EXPORT inline int getCommTag ()
    {
      static unsigned char counter = 0;
      return int( counter++ ) + 1536;
    }

  } // namespace __SPGrid



  // SPCommunication
  // ---------------

  template< class Grid, class DataHandle >
  struct SPCommunication
  {
    static const int dimension = Grid::dimension;

    typedef SPGridLevel< Grid > GridLevel;
    typedef SPPartitionList< dimension > PartitionList;

    typedef typename DataHandle::DataType DataType;

    typedef typename GridLevel::CommInterface Interface;

  private:
    typedef SPPackedMessageWriteBuffer< typename Grid::Communication > WriteBuffer;
    typedef SPPackedMessageReadBuffer< typename Grid::Communication > ReadBuffer;

  public:
    SPCommunication ( const GridLevel &gridLevel, DataHandle &dataHandle,
                      InterfaceType iftype, CommunicationDirection dir );

    SPCommunication ( const SPCommunication & ) = delete;
    SPCommunication ( SPCommunication &&other );

    ~SPCommunication () { wait(); }

    bool ready () const { return !( bool( interface_ )); }

    void wait ();

    [[deprecated]]
    bool pending () const { return !ready(); }
  private:
    const GridLevel &gridLevel_;
    DataHandle &dataHandle_;
    const Interface *interface_;
    CommunicationDirection dir_;
    int tag_;
    bool fixedSize_;
    std::vector< WriteBuffer > writeBuffers_;
    std::vector< ReadBuffer > readBuffers_;
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
      tag_( __SPGrid::getCommTag() ),
      fixedSize_( true )
  {
    for( int codim = 0; codim <= dimension; ++codim )
      fixedSize_ &= !dataHandle_.contains( dimension, codim ) || dataHandle_.fixedSize( dimension, codim );

    const std::size_t numLinks = interface_->size();
    readBuffers_.reserve( numLinks );

    if( fixedSize_ )
    {
      for( typename Interface::Iterator it = interface_->begin(); it != interface_->end(); ++it )
      {
        readBuffers_.emplace_back( gridLevel.grid().comm() );
        std::size_t size = 0;
        const PartitionList &partitionList = it->receiveList( dir );
        Hybrid::forEach( std::make_integer_sequence< int, dimension+1 >(), [ this, &partitionList, &size ] ( auto codim ) {
            typedef SPPartitionIterator< codim, const Grid > Iterator;

            if( !dataHandle_.contains( dimension, codim ) )
              return;

            const Iterator end( gridLevel_, partitionList, typename Iterator::End() );
            for( Iterator it( gridLevel_, partitionList, typename Iterator::Begin() ); it != end; ++it )
              size += dataHandle_.size( *it );
          } );
        size *= sizeof( DataType );
        readBuffers_.back().receive( it->rank(), tag_, size );
      }
    }

    writeBuffers_.reserve( numLinks );
    for( typename Interface::Iterator it = interface_->begin(); it != interface_->end(); ++it )
    {
      writeBuffers_.emplace_back( gridLevel.grid().comm() );
      const PartitionList &partitionList = it->sendList( dir );
      Hybrid::forEach( std::make_integer_sequence< int, dimension+1 >(), [ this, &partitionList ] ( auto codim ) {
          typedef SPPartitionIterator< codim, const Grid > Iterator;

          if( !dataHandle_.contains( dimension, codim ) )
            return;

          const bool fixedSize = dataHandle_.fixedSize( dimension, codim );
          const Iterator end( gridLevel_, partitionList, typename Iterator::End() );
          for( Iterator it( gridLevel_, partitionList, typename Iterator::Begin() ); it != end; ++it )
          {
            const auto &entity = *it;
            if( !fixedSize )
              writeBuffers_.back().write( static_cast< int >( dataHandle_.size( entity ) ) );
#ifndef NDEBUG
            const std::size_t posBeforeGather = writeBuffers_.back().position();
#endif // #ifndef NDEBUG
            dataHandle_.gather( writeBuffers_.back(), entity );
#ifndef NDEBUG
            const std::size_t posAfterGather = writeBuffers_.back().position();
            const std::size_t sizeInBytes = dataHandle_.size( entity ) * sizeof( DataType );
            if( posAfterGather - posBeforeGather != sizeInBytes )
              DUNE_THROW( GridError, "Number of bytes written (" << (posAfterGather - posBeforeGather) << ") does not coincide with reported size (" << sizeInBytes << ")" );
#endif // #ifndef NDEBUG
          }
        } );
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
      fixedSize_( other.fixedSize_ ),
      writeBuffers_( std::move( other.writeBuffers_ ) ),
      readBuffers_( std::move( other.readBuffers_ ) )
  {
    other.interface_ = nullptr;
  }


  template< class Grid, class DataHandle >
  inline void SPCommunication< Grid, DataHandle >::wait ()
  {
    if( ready() )
      return;

    const std::size_t numLinks = interface_->size();

    if( !fixedSize_ )
    {
      for( std::size_t i = 0; i < numLinks; ++i )
      {
        readBuffers_.emplace_back( gridLevel_.grid().comm() );
        readBuffers_.back().receive( tag_ );
      }
    }

    for( std::size_t i = 0; i < numLinks; ++i )
    {
      const typename std::vector< ReadBuffer >::iterator buffer = waitAny( readBuffers_ );
      for( typename Interface::Iterator it = interface_->begin(); it != interface_->end(); ++it )
      {
        if( it->rank() == buffer->rank() )
        {
          const PartitionList &partitionList = it->receiveList( dir_ );
          Hybrid::forEach( std::make_integer_sequence< int, dimension+1 >(), [ this, &partitionList, buffer ] ( auto codim ) {
              typedef SPPartitionIterator< codim, const Grid > Iterator;

              if( !dataHandle_.contains( dimension, codim ) )
                return;

              const bool fixedSize = dataHandle_.fixedSize( dimension, codim );
              const Iterator end( gridLevel_, partitionList, typename Iterator::End() );
              for( Iterator it( gridLevel_, partitionList, typename Iterator::Begin() ); it != end; ++it )
              {
                const auto &entity = *it;

                int size;
                if( !fixedSize )
                  buffer->read( size );
                else
                  size = dataHandle_.size( entity );
#ifndef NDEBUG
                const std::size_t posBeforeGather = buffer->position();
#endif // #ifndef NDEBUG
                dataHandle_.scatter( *buffer, entity, size );
#ifndef NDEBUG
                const std::size_t posAfterGather = buffer->position();
                const std::size_t sizeInBytes = static_cast< std::size_t >( size ) * sizeof( DataType );
                if( posAfterGather - posBeforeGather != sizeInBytes )
                  DUNE_THROW( GridError, "Number of bytes read (" << (posAfterGather - posBeforeGather) << ") does not coincide with reported size (" << sizeInBytes << ")" );
#endif // #ifndef NDEBUG
              }
            } );
          break;
        }
      }
    }
    readBuffers_.clear();

    for( typename std::vector< WriteBuffer >::iterator it = writeBuffers_.begin(); it != writeBuffers_.end(); ++it )
      it->wait();
    writeBuffers_.clear();

    interface_ = nullptr;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_COMMUNICATION_HH
